# -*- coding = utf-8 -*-
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def seq_translate(name, seq):
    """
    获取junction序列的转录序列
    return: 蛋白的SeqRecord对象
    """
    trans_seq = seq.translate()
    return SeqRecord(id=name, seq=trans_seq, description="")


def seq_reverse_translate(name, seq):
    """
    获取junction序列的反向互补序列的转录序列
    return: 蛋白的SeqRecord对象
    """
    reverse_trans_seq = seq.reverse_complement().translate()
    return SeqRecord(id=name, seq=reverse_trans_seq, description="")


def file_witre(seqs_liat, outfile):
    """
    写入文件
    seqs_liat: 前面生成的六种读码序列列表
    outfile: 输出文件地址
    """
    with open(outfile, "w") as fi:
        for S in seqs_liat:
            fi.write(">"+S.id+"\n")
            fi.write(str(S.seq+"\n"))


def junctions(circrnafile, outfile):
    """
    用于生成MS文件前置需要的六种 junction 文件
    circrnafile: circRNA的fasta文件
    outfile: 要输出的 junction 文件地址
    """
    Seqs = SeqIO.parse(circrnafile, "fasta")
    junction_list = []
    junction_seq = []

    # 生成junction位点附近序列
    for S in Seqs:
        S.name = str(S.name).split(":")[0]          # 去除host基因
        length = len(S.seq)
        if length > 200:
            S.seq = S.seq[-100:] + S.seq[:100]
        elif length <= 200:
            S.seq = S.seq[-int(length/2):] + S.seq[:int(length/2)]
        junction_seq.append(S)

    # 生成六种读码框序列
    for S in junction_seq:
        S1 = seq_translate("1_"+S.name, S.seq[0:])
        S2 = seq_translate("2_"+S.name, S.seq[1:])
        S3 = seq_translate("3_"+S.name, S.seq[2:])
        S4 = seq_reverse_translate("r1_"+S.name, S.seq[:])
        S5 = seq_reverse_translate("r2_"+S.name, S.seq[:-1])
        S6 = seq_reverse_translate("r3_"+S.name, S.seq[:-2])

        junction_list.extend([S1, S2, S3, S4, S5, S6])

    file_witre(junction_list, outfile)
