# _*_ coding : utf-8 _*_
# import Deepcircm6a.predict as predict
import multiprocessing
from itertools import count

import numpy as np
import torch
from Bio import SeqIO
from ...deep_learning.Deepcircm6a.model_one_hot_NCP_EIIP import *
from ...deep_learning.Deepcircm6a.seq_load_one_hot_NCP_EIIP import convert_seq_to_bicoding
from ...utils.logs import logger


# 读入fasta
def fasta_read(file):
    Seqs = SeqIO.parse(file, "fasta")
    return Seqs


# 依据模式，提取序列中A位点前后25bp
def predict_sequence_generate(sequence, mode):
    """
    提取A位点前后 25bp 序列用以符合 DeepCircm6a 的序列输入
    liner 模式    按照线性提取，所以序列前后 25bp 会提取不到
    circ 模式     按照环状提取，对于序列前后的 25bp，通过前后补齐提取

    sequence: 要处理的序列
    mode: 提取的模式
    return: {} 字典，键为A位点位置，值为提取出的序列
    """
    if len(sequence) < 51:
        return {}

    if mode == "linear":
        sequence_dict = {}
        for i in range(25, len(sequence)-25):
            if sequence[i] == "A":
                sequence_dict[i] = sequence_dict.get(i, sequence[i-25:i+26])

    elif mode == "circular":
        sequence_dict = {}
        length = len(sequence)
        sequence = sequence*3
        for i in range(length):
            if sequence[i] == "A":
                sequence_dict[i] = sequence_dict.get(i, sequence[length+i-25:length+i+26])
    else:
        raise ValueError(f"{mode}Invalid parameter, please enter 'linear' or 'circular'")
    return sequence_dict


def sequence_encode(data):
    encode_list = []
    for i in data:
        bicoding = convert_seq_to_bicoding(i)
        encode_list.append(bicoding)
    return encode_list


# 预测模块
def predict(model, x):
    model.eval()
    fx = model.forward(x)
    return fx


def predict_score(sequence_dict, model, drive):
    fa_header = list(sequence_dict.keys())

    if fa_header == []:
        return {}

    X_test = sequence_encode(sequence_dict.values())
    X_test = np.array(X_test)
    X_test = X_test.reshape(X_test.shape[0], int(X_test.shape[1] / 8), 8)
    X_test = torch.from_numpy(X_test).float()
    X_test = X_test.to(drive)


    batch_size = 256
    i = 0
    N = X_test.shape[0]
    y_pred_test = {}
    while i + batch_size < N:
        x_batch = X_test[i:i + batch_size]
        header_batch = fa_header[i:i + batch_size]

        fx = predict(model, x_batch)
        prob_data = F.log_softmax(fx, dim=1).cpu().data.numpy()
        for m in range(len(prob_data)):
            y_pred_test[header_batch[m]] = np.exp(prob_data)[m][1]
        i += batch_size

    x_batch = X_test[i:N]
    header_batch = fa_header[i:N]
    fx = predict(model, x_batch)
    prob_data = F.log_softmax(fx, dim=1).cpu().data.numpy()
    for m in range(len(prob_data)):
        y_pred_test[header_batch[m]] = np.exp(prob_data)[m][1]

    return y_pred_test


def m6a_score(recode, model_path, mode):
    torch.set_num_threads(1)
    HIDDEN_NUM = 128
    LAYER_NUM = 3
    FC_DROPOUT = 0.5
    RNN_DROPOUT = 0.5
    CELL = 'LSTM'
    # 初始化model
    if torch.cuda.is_available():
        drive = torch.device('cuda')
        checkpoint = torch.load(model_path, map_location=torch.device('cuda'))
    else:
        drive = torch.device('cpu')
        checkpoint = torch.load(model_path, map_location=torch.device('cpu'))
    model = CNN51_RNN(HIDDEN_NUM, LAYER_NUM, FC_DROPOUT, RNN_DROPOUT, CELL)
    model.load_state_dict(checkpoint['state_dict'])
    model.to(drive)

    result = []
    for Seq in recode:
        sequence_dict = predict_sequence_generate(Seq.seq, mode)
        y_pred_score = predict_score(sequence_dict, model, drive)
        sequence = list(Seq.seq)
        score = ["0"]*len(sequence)
        for k in y_pred_score.keys():
            score[k] = str(y_pred_score[k])
        y_pred_score[-1] = 0
        cnt = sum(1 for v in y_pred_score.values() if v > 0.5)
        result.append([Seq.id, cnt, sequence, score])
    return result


# @profile
def DeepCircm6a_score(file, mode, model_path, outfile, thread):

    logger.info("Loading m6A data")
    Seqs = fasta_read(file)
    recode = list(Seqs)
    total_number = len(recode)
    step = int(total_number / thread) + 1

    logger.info(f"Loading Deepcircm6A model and prediction results, Number of processes used: {thread}")
    pool = multiprocessing.Pool(thread)
    pool_result = []
    torch.set_num_threads(1)
    for i in range(0, total_number, step):
        result = pool.apply_async(m6a_score, (recode[i:i+step], model_path, mode))
        pool_result.append(result)
    pool.close()
    pool.join()

    logger.info("Write data")
    final_result = []
    for result in pool_result:
        final_result.extend(result.get())
    fi = open(outfile, "w")
    for line in final_result:
        fi.write(">" + line[0] + "\t" + str(line[1]) + "\n")
        fi.write("\t".join(line[2]) + "\n")
        fi.write("\t".join(line[3]) + "\n")
    fi.close()
    logger.info("End of Deepcircm6A operation")


