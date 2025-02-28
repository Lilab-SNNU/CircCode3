# _*_ coding : utf-8 _*_

import pandas as pd
from Bio import SeqIO
from ..utils.logs import logger
from ..utils.path_utils import *

def rem_peptide(file, outfile):

    peptide_table = pd.read_table(file)
    group = peptide_table.groupby("sequence")
    number = group["circRNA"].value_counts().to_frame()
    number.to_csv(outfile, sep="\t")


def merge_peptide(peptide_list):
    sum_num = str(sum([peptide_list[i] for i in range(1, len(peptide_list), 3)]))
    peptide_list = [[peptide_list[i], peptide_list[i + 2]] for i in range(0, len(peptide_list), 3)]
    peptide_list = sorted(peptide_list, key=lambda x: x[1])
    peptide = peptide_list[0][0]
    site = peptide_list[0][1]
    peptide_list.pop(0)
    while len(peptide_list) != 0:
        if site + len(peptide) >= peptide_list[0][1]+len(peptide_list[0][0]) and peptide_list[0][1] >= site:
            peptide_list.pop(0)
        elif peptide_list[0][1]+len(peptide_list[0][0])> site + len(peptide) > peptide_list[0][1]:
            a = site + len(peptide)-peptide_list[0][1]
            peptide = peptide + peptide_list[0][0][a:]
            peptide_list.pop(0)
        elif peptide_list[0][1]+len(peptide_list[0][0])> site > peptide_list[0][1]:
            a = site - peptide_list[0][1]
            peptide = peptide_list[0][0][:a] + peptide
            site = peptide_list[0][1]
            peptide_list.pop(0)

    return [str(sum_num), peptide, str(site*3)]


def mapping_to_junction(junction, remove_peptide, peptide_file):

    dict_junction = {}
    fi = SeqIO.parse(junction, "fasta")
    for Seq in fi:
        dict_junction[Seq.id] = Seq.seq
    table = {}
    with open(remove_peptide) as fi:
        for line in fi.readlines()[1:]:
            li = line.strip().split("\t")
            peptide = li[0]
            if peptide.startswith("r"):
                continue
            circRNA = li[1].split("/")
            num = int(li[2])
            for name in circRNA:
                circ_junction = dict_junction[name]
                start = circ_junction.find(peptide)
                if start+len(peptide) > len(circ_junction)/2 > start >=0:
                    l = table.get(name, [])
                    l.extend([peptide, num, int(start-len(circ_junction)/2)])
                    table[name] = l

    # 获取最长肽段
    for k in table.keys():
        peptide_list = table[k]
        peptide_site = merge_peptide(peptide_list)
        peptide_site.extend([peptide_list[i] for i in range(0, len(peptide_list), 3)])
        table[k] = peptide_site

    with open(peptide_file, "w") as fi:
        for k in table.keys():
            fi.write(k+"\t"+"\t".join(table[k])+"\n")


def filter_coding_circ_peptide(peptide_path: PathManager):
    peptidefile = peptide_path.get_path("peptide")
    remove_peptide = os.path.join(peptide_path.get_path("tmp"), "remove_peptide.txt")
    peptides_table = os.path.join(peptide_path.get_path("tmp"), "peptide_table.txt")
    junction = peptide_path.get_path("junction")
    logger.info("Removing repetitive peptide segments")
    rem_peptide(peptidefile, remove_peptide)
    logger.info("Match to BSJ site")
    mapping_to_junction(junction, remove_peptide, peptides_table)
    return peptides_table

