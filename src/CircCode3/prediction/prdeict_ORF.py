# -*- coding = utf-8 -*-
import orfipy_core
import pandas as pd
import importlib.resources
from typing import Callable, List, Any
from Bio import Seq, SeqIO
from Bio.SeqRecord import SeqRecord
from concurrent.futures import as_completed, ProcessPoolExecutor
from ..command.command import Command
from ..utils.path_utils import *
from ..deep_learning.Deepcircm6a.sequence_m6a_predict import DeepCircm6a_score
from ..deep_learning.end_code_score.predict import DLMSC


def orf_predict(name: str, sequence: str) -> list[list[str | Any]]:
    """
    Predict all ORFs of circRNA
    :param name: circRNA name
    :param sequence: circRNA Sequence
    :return: orf list
    """
    length = len(sequence)
    sequence = str(sequence)
    quadruple_sequence = sequence * 4
    orf_list = []
    possible_orf = orfipy_core.orfs(quadruple_sequence, minlen=60, starts=["ATG"], strand="f", partial3=True)
    for start, stop, strand, description in possible_orf:
        info = description.split()
        if start < length < stop:
            if info[-1] == "Stop:NA":
                stop = "NA"
            orf_list.append([name, start, stop])
    return orf_list


def read_coding_reads(table_file: str) -> dict[str, list[int]]:
    """
    Read the reads file to obtain the number of reads and the longest length after concatenation
    :param table_file: clear reads file
    :return: Dictionary with circRNA names as keys
    """
    reads_site = {}
    with open(table_file) as fi:
        reads = fi.readlines()
        for i in range(0, len(reads), 2):
            line_name = reads[i].strip().split("\t")
            name = line_name[0]
            num = int(line_name[1])
            line_reads = reads[i+1].strip().split("\t")
            start = int(line_reads[1])
            end = int(line_reads[2])
            for i in range(0, len(line_reads), 3):
                start = min(int(line_reads[i+1]), start)
                end = max(int(line_reads[i+2]), end)
            reads_site[name] = [start, end, num]
    return reads_site


def orf_info_expansion_Ribo(name: str, sequence: str, reads: dict[str, list[int]]) -> list[list]:
    """
    Combining orf information with reads information
    :param name: circRNA name
    :param sequence: circRNA sequence
    :param reads: Dictionary with circRNA names as keys
    :return: list of orf info
    """
    orf_list = orf_predict(name, sequence)
    all_orf_info = []
    for info in orf_list:
        seq = sequence*4
        start = info[1]
        stop = -1 if info[2] == "NA" else info[2]
        orf = seq[start:stop]
        corf = Seq.Seq(orf).translate()
        all_orf_info.append([name, reads[2], reads[0], reads[1], info[1], info[2], corf])
    return all_orf_info


def mult_thread_orf_fillter_Ribo(result_path: PathManager, threads: int) -> tuple[str, str]:
    """
    Using multiple processes for ORF prediction
    :param result_path: PathManager
    :param threads: threads number
    :return: CircRNA fasta pathway and ORF information with reads support
    """
    logger.info(f"Start find Coding circRNA orf")
    table_file = result_path.get_path("coding_frag")
    circRNA_file =  result_path.get_path("circRNA")
    true_circRNA = os.path.join(result_path.get_path("tmp"), "true_circRNA.fasta")
    reads_info = os.path.join(result_path.get_path("tmp"), "Ribo_frag_ORF_info.xlsx")

    dict_reads = read_coding_reads(table_file)
    Seqs = SeqIO.parse(circRNA_file, "fasta")
    orf_fillter_list = []
    seq_list = []
    tasks = []
    with ProcessPoolExecutor(max_workers=threads) as executor:
        for S in Seqs:
            name = S.id
            if name not in dict_reads:
                continue
            seq_list.append(S)
            tasks.append(executor.submit(orf_info_expansion_Ribo, name, str(S.seq), dict_reads[name]))

        for future in as_completed(tasks):
            orf_fillter_list.extend(future.result())

    SeqIO.write(seq_list, true_circRNA, "fasta")
    if len(orf_fillter_list) == 0:
        logger.warning(f"Found zero CircRNA specific ORFs, program has exited")
        exit(0)

    logger.info(f"Discovery of {len(orf_fillter_list)} ORFs from circRNAs with translated fragments")
    circRNA_info_pd = pd.DataFrame(orf_fillter_list, columns=["circRNA_name", "Reads_number", "reads_Start", "reads_Stop", "ORF_Start", "ORF_Stop", "cORF"])
    circRNA_info_pd.to_excel(reads_info, index=False)
    return reads_info, true_circRNA


def read_coding_peptide(table_file: str) -> dict[str, Any]:
    """
    Read the peptide file to obtain the number of peptide and the longest length after concatenation
    :param table_file: clear reads file
    :return: Dictionary with circRNA names as keys
    """
    peptide_site = {}
    with open(table_file) as fi:
        for line in fi.readlines():
            li = line.strip().split("\t")
            name = li[0].split("_", 1)[1]
            value = peptide_site.get(name, [])
            value.append([li[1], li[2], int(li[3])])
            peptide_site[name] = value
    return peptide_site


def orf_info_expansion_MS(name: str, sequence: str, peptides: list[Any]) -> list[list]:
    """
    Combining orf information with peptide information
    :param name: circRNA name
    :param sequence: circRNA sequence
    :param peptide: Dictionary with circRNA names as keys
    :return: list of orf info
    """
    orf_list = orf_predict(name, sequence)
    all_orf_info = []
    for info in orf_list:
        seq = sequence*4
        start = info[1]
        stop = -1 if info[2] == "NA" else info[2]
        orf = seq[start:stop]
        corf = Seq.Seq(orf).translate()
        for peptide in peptides:
            if corf.find(peptide[1]) != -1:
                all_orf_info.append([name, peptide[0], peptide[2], peptide[1], info[1], info[2], corf])
    return all_orf_info


def mult_thread_orf_fillter_MS(result_path:PathManager, threads: int) -> tuple[str, str]:
    """
    Using multiple processes for ORF prediction
    :param result_path: PathManager
    :param threads: threads number
    :return: CircRNA fasta pathway and ORF information with peptide support
    """
    logger.info(f"Start find Coding circRNA orf")
    table_file = result_path.get_path("coding_frag")
    circRNA_file =  result_path.get_path("circRNA")
    true_circRNA = os.path.join(result_path.get_path("tmp"), "true_circRNA.fasta")
    peptide_info = os.path.join(result_path.get_path("tmp"), "MS_frag_ORF_info.xlsx")

    dict_reads = read_coding_peptide(table_file)
    Seqs = SeqIO.parse(circRNA_file, "fasta")
    orf_fillter_list = []
    seq_list = []
    tasks = []
    with ProcessPoolExecutor(max_workers=threads) as executor:
        for S in Seqs:
            name = S.id
            if name not in dict_reads:
                continue
            seq_list.append(S)
            tasks.append(executor.submit(orf_info_expansion_MS, name, str(S.seq), dict_reads[name]))

        for future in as_completed(tasks):
            orf_fillter_list.extend(future.result())

    SeqIO.write(seq_list, true_circRNA, "fasta")
    if len(orf_fillter_list) == 0:
        logger.warning(f"Found zero CircRNA specific ORFs, program has exited")
        exit(0)

    logger.info(f"Discovery of {len(orf_fillter_list)} ORFs from circRNAs with translated fragments")
    circRNA_info_pd = pd.DataFrame(orf_fillter_list, columns=["circRNA_name", "peptide_number", "peptide_Start", "peptide_total_sequence", "ORF_Start", "ORF_Stop", "cORF"])
    circRNA_info_pd.to_excel(peptide_info, index=False)
    return [peptide_info, true_circRNA]


def upstream_seq_extract(sequence: str, site: tuple[int, int]) -> str:
    """
    Extract the sequence of ORF upstream 101
    :param sequence: sequence
    :param site: orf initiator codon and termination codon positions
    :return: sequence for ires
    """
    l = len(sequence)
    sequence = sequence*4
    end = int(site[0] + l*3)
    start = end - 101
    return sequence[start: end]


# 提取上游位点151bp,end为ORF起始、上游结束位点，用于m6a位点预测
def upstream_m6a_seq_extract(sequence: str, site: tuple[int, int]) -> str:
    """
    Extract the sequence of ORF upstream 151
    :param sequence: sequence
    :param site: orf initiator codon and termination codon positions
    :return: sequence for m6a
    """
    l = len(sequence)
    sequence = sequence*5
    end = int(site[0] + l*3)
    start = end - 126
    return sequence[start: end + 25]


# 提取ORF序列，site为位置元组
def orf_seq_extract(sequence: str, site: tuple[int, int]) -> str:
    """
    Extract the sequence of ORF
    :param sequence: sequence
    :param site: orf initiator codon and termination codon positions
    :return: sequence of orf
    """
    sequence = sequence*4
    return sequence[site[0]: site[1]]


# 提取终止密码子前后序列，site为ORF终止位置
def end_codon_extract(sequence: str, site: tuple[int, int]) -> str:
    """
    Extract the sequence of ORF upstream 103
    :param sequence: sequence
    :param site: orf initiator codon and termination codon positions
    :return: sequence for termination codon
    """
    if site[1] == -1:
        return None
    sequence = sequence*2
    return sequence[site[1]-51: site[1]+54]


def read_orf_info(path: str) -> list[list]:
    """
    Extract the orf info from orf_info file
    :param path: path of xxx_frag_ORF_info.xlsx
    :return: info list
    """
    data_table = pd.read_excel(path)
    info = data_table[["circRNA_name", "ORF_Start", "ORF_Stop"]]
    return info.values.tolist()


def read_sequence(path: str):
    """
    Read the fasta file and return a simple index dictionary and sequence list
    :param path: path of circRNA sequence
    :return: index dictionary and sequence list
    """
    Seqs = SeqIO.parse(path, "fasta")
    index = {}
    sequence = []
    for i, S in enumerate(Seqs):
        index[S.id] = i
        sequence.append(S.seq)
    return index, sequence


def sequence_operation(orf_info:list, sequence: list, index: dict, filename: str, func: Callable[[str, tuple[int]], str]):
    """
    Encapsulate sequence extraction function for multi-threaded operations
    :param orf_info: orf info
    :param sequence: circRNA sequence
    :param index: index dictionary
    :param filename: output file name
    :param func: sequence function
    :return:
    """
    target_sequence = []

    for circRNA in orf_info:
        name = circRNA[0]
        start = circRNA[1]
        end = -1 if circRNA[2] == "NA" else circRNA[2]
        seq = sequence[index[name]]
        mark_name = name+f"|({start},{end})"
        seq = func(seq, (start, end))
        if seq == None:
            continue
        seqrecode = SeqRecord(id=mark_name, seq=Seq.Seq(seq), description="")
        target_sequence.append(seqrecode)

    SeqIO.write(target_sequence, filename, format="fasta")


def sequence_extract(orf_info_path: PathManager, threads: int) -> tuple[str, str, str, str]:
    """
    Extract the sequences required for the three tools
    :param orf_info_path: PathManager
    :param threads: threads
    :return: path of file
    """
    circ_info_path = orf_info_path.get_path("circRNA info")
    true_circRNA = orf_info_path.get_path("Coding circRNA Sequence")

    IRES_sequence = os.path.join(orf_info_path.get_path("tmp"), "candidate_IRES_Sequence.fasta")
    m6A_sequence = os.path.join(orf_info_path.get_path("tmp"), "candidate_m6A_Sequence.fasta")
    termination_codon_sequence = os.path.join(orf_info_path.get_path("tmp"), "candidate_termination_codon_Sequence.fasta")
    orf_sequence = os.path.join(orf_info_path.get_path("tmp"), "orf_Sequence.fasta")

    orf_info =  read_orf_info(circ_info_path)
    index, sequence = read_sequence(true_circRNA)

    logger.info("Extract ORF information sequence.")
    with ProcessPoolExecutor(max_workers=threads) as executor:
        executor.submit(sequence_operation,
                        orf_info, sequence, index, IRES_sequence,
                        upstream_seq_extract)
        executor.submit(sequence_operation,
                        orf_info, sequence, index, m6A_sequence,
                        upstream_m6a_seq_extract)
        executor.submit(sequence_operation,
                        orf_info, sequence, index, orf_sequence,
                        orf_seq_extract)
        executor.submit(sequence_operation,
                        orf_info, sequence, index, termination_codon_sequence,
                        end_codon_extract)

    return IRES_sequence, m6A_sequence, orf_sequence, termination_codon_sequence


def get_resource_path(package, file):
    with importlib.resources.path(package, file) as p:
        return p

def invoke_IRES_finder(sequence_file, IRES_score_file):
    """
    Call IRES finder
    :param sequence_file: The sequence file to be predicted
    :param IRES_score_file: output result
    :return:
    """
    ires_path = get_resource_path("CircCode3.deep_learning.IRESfinder_final","IRESfinder.py")

    ires_finder = Command("python3", ires_path,
                          "-f", sequence_file,
                          "-o", IRES_score_file)
    ires_finder.command_running()


def invoke_DLMSC(sequence_file, end_score_file):
    """
    Call DLMSC
    :param sequence_file: The sequence file to be predicted
    :param end_score_file: output result
    :return:
    """
    model_path = get_resource_path("CircCode3.resources.stop_codon_model","checkpoint_ATCG.pth.tar")

    DLMSC(sequence_file, model_path, end_score_file)


def invoke_DeepCircm6A(sequence_file, m6A_score_file, threads):
    """
    Call DeepCircm6A
    :param sequence_file: The sequence file to be predicted
    :param m6A_score_file: output result
    :param threads: threads
    :return:
    """
    model_path = get_resource_path("CircCode3.resources.DeepCircm6A","checkpoint.pth.tar")
    DeepCircm6a_score(sequence_file, "linear", model_path, m6A_score_file, threads)


def estimate_circRNA_ORF(estimate_path: PathManager, threads):
    """
    Call three tools to evaluate ORF
    :param estimate_path: PathManager
    :param threads: number of threads
    :return:
    """
    IRES_sequence_path = estimate_path.get_path("IRES Sequence")
    m6A_sequence_path = estimate_path.get_path("m6A sequence")
    termination_codon_sequence_path = estimate_path.get_path("termination codon sequence")

    IRES_sequence_result = os.path.join(estimate_path.get_path("tmp"), "IRES_score_file.txt")
    m6A_sequence_result = os.path.join(estimate_path.get_path("tmp"), "m6a_sequence_file.txt")
    termination_codon_sequence_result = os.path.join(estimate_path.get_path("tmp"), "end_score_file.txt")

    invoke_IRES_finder(IRES_sequence_path, IRES_sequence_result)
    invoke_DeepCircm6A(m6A_sequence_path, m6A_sequence_result, threads)
    invoke_DLMSC(termination_codon_sequence_path, termination_codon_sequence_result)


    return [IRES_sequence_result, m6A_sequence_result, termination_codon_sequence_result]


def read_table(file: str, head=False) -> dict[str,dict[str]]:
    """
    Read the results of IRES and DLMSC
    :param file: file path
    :param head:
    :return: Dictionary with circRNA names as keys
    """
    with open(file, "r") as fi:
        text = fi.readlines()
        if head:
            context = [t.strip().split("\t") for t in text[1:]]
        else:
            context = [t.strip().split("\t") for t in text]

    context_dict = {}
    for l in context:
        li = l[0].rsplit("|", 1)
        circname, site = li[0], li[1]
        value = context_dict.get(circname, {})
        value[site] = l[-1]
        context_dict[circname] = value
    return context_dict

def m6a_sorce_read(file: str) -> dict[str,dict[str, int]]:
    """
    Read the results of DeepCircm6A
    :param file: file path
    :return:
    """
    with open(file, "r") as fi:
        text = fi.readlines()
        context = [t[1:].strip().split("\t") for t in text[::3]]
    context_dict = {}
    for l in context:
        li = l[0].rsplit("|", 1)
        circname, site = li[0], li[1]
        value = context_dict.get(circname, {})
        value[site] = l[-1]
        context_dict[circname] = value
    return context_dict


def estimate_info_collect(estimate_path: PathManager, type: str) -> str:
    """
    Read the ORF evaluation results and merge the table
    :param estimate_path:PathManager
    :param type: data type : ms or ribo
    :return: file path
    """
    IRES_score_path = estimate_path.get_path("IRES score")
    m6A_score_path = estimate_path.get_path("m6A score")
    termination_codon_score_path = estimate_path.get_path("termination codon score")
    circ_info_path = estimate_path.get_path("circRNA info")

    orf_info_path = os.path.join(estimate_path.get_path("tmp"), f"{type}_ORF_evaluate.xlsx")

    ires = read_table(IRES_score_path, True)
    end_codon = read_table(termination_codon_score_path)
    m6a = m6a_sorce_read(m6A_score_path)

    orf_info =  pd.read_excel(circ_info_path)
    orf_info_list = orf_info.values.tolist()
    final_table = []
    for i in range(len(orf_info_list)):
        circ_name = orf_info_list[i][0]
        site = f"({orf_info_list[i][4]},{orf_info_list[i][5]})"
        line = [circ_name, orf_info_list[i][4], orf_info_list[i][5],
                ires[circ_name].get(site, None), m6a[circ_name].get(site, None), end_codon[circ_name].get(site, None)]
        final_table.append(line)

    circRNA_info_pd = pd.DataFrame(final_table, columns=["circRNA_name", "ORF_Start", "ORF_Stop", "IRES_Score", "m6A_Count", "Termination_Codon_Score"])
    total_info = pd.merge(orf_info, circRNA_info_pd, on=["circRNA_name", "ORF_Start", "ORF_Stop"], how="outer")
    cORF = total_info.pop("cORF")
    total_info.insert(total_info.shape[1], "cORF", cORF)
    total_info.to_excel(orf_info_path, index=False)
    return orf_info_path