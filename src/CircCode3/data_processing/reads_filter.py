# -*- coding = utf-8 -*-
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from ..command.command import Command
from ..utils.logs import logger
from ..utils.path_utils import *


def clear_line_reads(ribo_paths: PathManager, threads: int) -> list[str]:
    """
    :param threads:
    :param paths:
    :return:
    """
    transcript_fasta = ribo_paths.get_path("transcript")
    transcript_result_folder = os.path.join(ribo_paths.get_path("tmp"), "transcript")
    create_folder(transcript_result_folder)
    ribo_paths.add_path("transcript_index", os.path.join(transcript_result_folder, "transcript_index"))

    ribosome_fasta = ribo_paths.get_path("ribosome")
    ribosome_result_folder = os.path.join(ribo_paths.get_path("tmp"), "ribosome")
    create_folder(ribosome_result_folder)
    ribo_paths.add_path("ribosome_index", os.path.join(ribosome_result_folder, "ribosome_index"))

    transcript_index_build_command = Command("bowtie2-build", "--threads", threads, transcript_fasta,
                                             ribo_paths["transcript_index"])
    transcript_index_build_command.command_running()

    ribosome_index_build_command = Command("bowtie2-build", "--threads", threads, ribosome_fasta,
                                           ribo_paths["ribosome_index"])
    ribosome_index_build_command.command_running()
    reads_path = ribo_paths.get_path("reads")

    if len(reads_path) == 1:
        ribo_paths.add_path("transcript_unmap_reads", [os.path.join(transcript_result_folder, "transcript_unmap.fastq")])
        transcript_mapping_command = Command("bowtie2", "-p", threads, "--un", ribo_paths["transcript_unmap_reads"][0],
                                             "--norc", "-x", ribo_paths["transcript_index"],
                                             "-U", reads_path[0])
        transcript_mapping_command.command_running()

        ribo_paths.add_path("ribosome_unmap_reads", [os.path.join(ribosome_result_folder, "ribosome_unmap.fastq")])

        ribosome_mapping_command = Command("bowtie2", "-p", threads, "--un", ribo_paths["ribosome_unmap_reads"][0],
                                           "--norc", "-x", ribo_paths["ribosome_index"],
                                           "-U", ribo_paths["transcript_unmap_reads"][0])
        ribosome_mapping_command.command_running()

    elif len(reads_path) == 2:
        ribo_paths.add_path("transcript_unmap_reads", os.path.join(transcript_result_folder, "transcript_unmap"))
        transcript_mapping_command = Command("bowtie2", "-p", threads, "--un-conc", ribo_paths["transcript_unmap_reads"],
                                             "--norc", "-x", ribo_paths["transcript_index"],
                                             "-1", reads_path[0], "-2", reads_path[1])
        transcript_mapping_command.command_running()

        unmap_path = [ribo_paths["transcript_unmap_reads"] + "_1.fastq", ribo_paths["transcript_unmap_reads"] + "_1.fastq"]

        ribo_paths.add_path("ribosome_unmap_reads", os.path.join(ribosome_result_folder, "ribosome_unmap"))

        ribosome_mapping_command = Command("bowtie2", "-p", threads, "--un-conc", ribo_paths["ribosome_unmap_reads"],
                                           "--norc", "-x", ribo_paths["ribosome_index"],
                                           "-1", unmap_path[0], "-2", unmap_path[1])
        ribosome_mapping_command.command_running()

        ribo_paths["ribosome_unmap_reads"] = [ribo_paths["ribosome_unmap_reads"] + "_1.fastq",
                                         ribo_paths["ribosome_unmap_reads"] + "_2.fastq"]

    return ribo_paths.get_path("ribosome_unmap_reads")


def make_junction_seq(paths: PathManager) -> str:
    length = 50
    logger.info("make BSJ sequence...")
    junction_result_folder = os.path.join(paths.get_path("tmp"), "junction")
    create_folder(junction_result_folder)
    junction_file = os.path.join(junction_result_folder, "CircRNA_Junction.fasta")
    junction_sequence_list = []
    for sequence in SeqIO.parse(paths["circRNA"], 'fasta'):
        bsj_sequence = sequence.seq[-length:]+sequence.seq[:length]
        bsj_seqRecode = SeqRecord(bsj_sequence,
                                  id = sequence.id,
                                  description = '')
        junction_sequence_list.append(bsj_seqRecode)
    logger.info("write BSJ sequence...")
    SeqIO.write(junction_sequence_list, junction_file, "fasta")
    return junction_file


def filter_circ_reads(BSJ_paths: PathManager, threads: int) -> list[str]:
    junction_fasta = BSJ_paths.get_path("BSJ")
    junction_result_folder = os.path.join(BSJ_paths.get_path("tmp"), "junction")
    BSJ_paths.add_path("BSJ_index", junction_fasta)
    reads_path = BSJ_paths.get_path("reads")
    BSJ_paths.add_path("final_clear_result", os.path.join(junction_result_folder, "final_clear_result.sam"))

    junction_index_build_command = Command("bowtie2-build", "--threads", threads, junction_fasta,
                                             BSJ_paths["BSJ_index"])
    junction_index_build_command.command_running()

    if len(reads_path) == 1:
        junction_mapping_command = Command("bowtie2", "-p", threads, '-x', BSJ_paths["BSJ_index"],
                                           '-U', reads_path[0],
                                           "--norc", "--end-to-end",
                                           "-S", BSJ_paths["final_clear_result"])
        junction_mapping_command.command_running()
    elif len(reads_path) == 2:

        junction_mapping_command = Command("bowtie2", "-p", threads, '-x', BSJ_paths["BSJ_index"],
                                           '-1', reads_path[0], '-2', reads_path[1],
                                           "--norc", "--end-to-end",
                                           "-S", BSJ_paths["final_clear_result"])
        junction_mapping_command.command_running()
    return BSJ_paths["final_clear_result"]


def bam2bed(circ_ribo_path: PathManager) -> str:
        
    final_clear_result = circ_ribo_path.get_path("final_clear_result")
    bamfile = os.path.join(circ_ribo_path['tmp'], "candidate_reads.bam")
    bedfile = os.path.join(circ_ribo_path['tmp'], "candidate_reads.bed")

    sam_to_bam_command = Command("samtools", "view", "-bS", final_clear_result, "|", "samtools", "sort", "-o", bamfile)
    sam_to_bam_command.command_running()

    bam_to_bed_command = Command("bedtools", "bamtobed", "-i", bamfile, "-cigar", ">", bedfile)
    bam_to_bed_command.command_running()

    return bedfile


def read_bed(file):
    with open(file, "r") as fi:
        for line in fi:
            info = line.strip().split("\t")
            yield info[0], info[3], info[5], info[6], int(info[1])+1, int(info[2])


def cigar_analysis(cigar) -> tuple[int, int]:
    cigar_match = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    full_length = 0
    effective_length = 0
    for length, op_type in cigar_match:
        full_length += int(length)
        if op_type in "M=X":
            effective_length = int(length)
    return full_length, effective_length


def is_circ_reads(start, stop, cigar, strand, full_length, effective_length) -> bool:
    check = [start < 50 < stop,
             25 <= stop - start <= 35,
             strand == "+",
             "N" not in cigar,
             effective_length/full_length > 0.9]
    return all(check)


def filter_coding_circ_reads(circ_ribo_path: PathManager):

    logger.info(f"Start filter coding reads")
    bedfile = circ_ribo_path.get_path("bedfile")
    coding_reads_file = os.path.join(circ_ribo_path.get_path('tmp'), "coding_reads.txt")

    bed_generator = read_bed(bedfile)
    coding_reads:dict[str, list] = {}

    for circId, readsId, strand, cigar, start, stop in bed_generator:
        
        full_length, effective_length = cigar_analysis(cigar)
        if is_circ_reads(start, stop, cigar, strand, full_length, effective_length):
            coding_reads[circId] = coding_reads.get(circId, []) + [readsId, str(start-50), str(stop-50)]

    number = 0
    with open(coding_reads_file, "w") as fi:
        for circ_id, reads_info in coding_reads.items():
            n = len(reads_info)/3
            if n >= 3:
                number += 1
                fi.write(circ_id + f"\t{int(n)}" + "\n" + "\t".join(reads_info) + "\n")
    if number == 0:
        logger.info(f"No CircRNA has fragments that support translation, the program has exited")
        exit(0)
    logger.info(f"Filter coding reads Finshed\nTotal find {number} circRNA")
    return coding_reads_file