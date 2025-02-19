# -*- coding:utf-8 -*-

import re
import pandas as pd
from ..utils.logs import logger
from ..utils.path_utils import *


def MAXQuant(path: PathManager) -> str:
    """
    MAXQuant 内容抓取，储存为文件
    :param path:
    :return:
    """
    outfile = os.path.join(path.get_path("tmp"), "peptides.txt")
    peptide_table = pd.read_table(path.get_path("peptide"))
    result = peptide_table[["Sequence", "Proteins"]]
    result.columns = ["sequence", "circRNA"]
    result.loc[:, "circRNA"] = result["circRNA"].str.replace(';', '/')
    s = result["circRNA"].str.contains("CON__").fillna(True)
    result = result.drop(index=result[s].index)
    result.to_csv(outfile, sep="\t", index=False)
    return outfile


def remove_modifications(sequence: str) -> str:
    """
    去除肽段序列中的修饰标记
    :param sequence: 含有修饰标记的肽段序列
    :return: 清理后的肽段序列
    """
    # 使用正则表达式去除修饰，这里简单地移除以 '+' 开头并跟随数字和 '.' 的部分
    import re
    return re.sub(r'\+\d+(\.\d+)?', '', sequence)

def MSGF(path: 'PathManager') -> str:
    """
    MSGF 内容抓取，储存为文件，并去除肽段序列中的修饰标记
    :param path: 路径管理对象
    :return: 输出文件路径
    """
    outfile = os.path.join(path.get_path("tmp"), "peptides.txt")
    peptide_table = pd.read_table(path.get_path("peptide"))
    
    result = peptide_table[["Peptide", "Protein"]]
    result.columns = ["sequence", "circRNA"]
    result.loc[:, "circRNA"] = result["circRNA"].str.replace(',', '/')

    result.loc[:, "sequence"] = result["sequence"].apply(remove_modifications)
    result.loc[:, "sequence"] = result['sequence'].str.split(".").apply(lambda x: x[1] if len(x) > 1 else x[0])    
    result.to_csv(outfile, sep="\t", index=False)
    
    return outfile


def format_adapt(path: PathManager) -> str:
    """
    文件类型判断并调用对应函数。没有则报错。
    :param path: PathManager
    :return:
    """
    with open(path.get_path("peptide"), "r", encoding="utf-8") as fi:
        title = fi.readline()
    if re.match(r"#SpecFile\tSpecID\tScanNum", title):
        logger.info("Determine the type of input result as: MSGF")
        outfile = MSGF(path)
    elif re.match(r"Sequence\tN-term cleavage window\tC-term cleavage window", title.lstrip('\ufeff')):
        logger.info("Determine the type of input result as: MAXQuant")
        outfile = MAXQuant(path)
    else:
        logger.error("Failed to identify the type of rustic software for the input file.Please enter the MAXQuant or MSGF results. or use --clear_peptide")
        exit(1)

    return outfile
