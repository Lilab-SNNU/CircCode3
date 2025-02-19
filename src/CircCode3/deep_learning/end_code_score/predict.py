# _*_ coding : utf-8 _*_

from ...deep_learning.end_code_score.model_one_hot import *
from ...utils.logs import logger

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
wordvec_len = 4
FC_DROPOUT = 0.5

def predict(model, x):
    model.eval() #evaluation mode do not use drop out
    fx = model.forward(x)
    return fx


def DLMSC(predict_fa, model_path, outfile):
    checkpoint = torch.load(model_path, map_location=torch.device('cpu'))
    logger.info("Loading DLMSC model")
    model = CNN51_RNN(FC_DROPOUT)
    model.load_state_dict(checkpoint['state_dict'])

    logger.info("Loading termination codon sequence data")
    X_test, fa_header = load_data_bicoding_with_header(predict_fa)
    X_test=np.array(X_test)
    X_test = X_test.reshape(X_test.shape[0], int(X_test.shape[1] / wordvec_len), wordvec_len)
    X_test = torch.from_numpy(X_test).float()


    batch_size = 256
    i = 0
    N = X_test.shape[0]

    logger.info("Start predicting")
    with open(outfile, 'w') as fw:
        while i + batch_size < N:
            x_batch = X_test[i:i + batch_size]
            header_batch = fa_header[i:i + batch_size]

            fx = predict(model, x_batch)
            prob_data = F.log_softmax(fx, dim=1).cpu().data.numpy()
            for m in range(len(prob_data)):
                fw.write(header_batch[m] + '\t' + str(np.exp(prob_data)[m][1]) + '\n')
            i += batch_size

        x_batch = X_test[i:N]
        header_batch = fa_header[i:N]
        fx = predict(model, x_batch)
        prob_data = F.log_softmax(fx, dim=1).cpu().data.numpy()
        for m in range(len(prob_data)):
            fw.write(header_batch[m] + '\t' + str(np.exp(prob_data)[m][1]) + '\n')
    logger.info("End of DLMSC operation!")











