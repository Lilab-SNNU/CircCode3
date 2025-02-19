# -*- coding = utf-8 -*-
import importlib
import pandas as pd
from ..data_processing import MS_format_adapter
from ..data_processing import  reads_filter
from ..data_processing import junction
from ..data_processing import peptide_match
from ..prediction import prdeict_ORF
from ..visualization import drawcirc
from ..utils.path_utils import *
from ..deep_learning.Deepcircm6a.sequence_m6a_predict import DeepCircm6a_score
from ..deep_learning.end_code_score.predict import DLMSC


def sub_command_Ribo_seq_process(args):
    threads = args.threads
    base_path = PathManager()
    base_path.add_path("circRNA", args.circRNA)
    base_path.add_path("output", project_name(args.output))
    base_path.add_path("result", os.path.join(base_path["output"], "Result"))
    base_path.add_path("transcript", args.transcript)
    base_path.add_path("ribosome", args.ribosome)
    base_path.add_path("reads", args.single if args.single else args.paired)
    create_folder(base_path["output"])
    Ribo_seq_process(args , base_path, threads)


def Ribo_seq_process(args , base_path:PathManager, threads):
    logger.info("\n" + "="*50 + "\nStart processing Ribo-seq data\n" + "="*50)
    base_path.add_path("tmp", os.path.join(base_path["output"], "RiboTemp"))
    create_folder(base_path["tmp"])

    logger.info("\n" + "="*30 + "\nClean up linear transcript reads from Ribo seq data")
    BSJ_paths = PathManager(base_path.paths)
    ribosome_unmap_reads = reads_filter.clear_line_reads(base_path, threads)
    BSJ_paths.add_path("reads", ribosome_unmap_reads)
    junction_file = reads_filter.make_junction_seq(base_path)
    BSJ_paths.add_path("BSJ", junction_file)
    final_clear_result = reads_filter.filter_circ_reads(BSJ_paths, threads)

    logger.info("\n" + "="*30 + "\nMapping clean data to CircRNA")
    circ_ribo_path = PathManager(base_path.paths)
    circ_ribo_path.add_path("final_clear_result", final_clear_result)
    bedfile = reads_filter.bam2bed(circ_ribo_path)
    circ_ribo_path.add_path("bedfile", bedfile)
    coding_reads_file = reads_filter.filter_coding_circ_reads(circ_ribo_path)

    logger.info("\n" + "="*30 + "\nSearching and evaluating ORF information of CircRNA")
    predict_path = PathManager(base_path.paths)
    predict_path.add_path("coding_frag", coding_reads_file)
    reads_info, true_circRNA = prdeict_ORF.mult_thread_orf_fillter_Ribo(predict_path, threads)
    predict_path.add_path("circRNA info", reads_info)
    predict_path.add_path("Coding circRNA Sequence", true_circRNA)
    IRES_sequence, m6A_sequence, orf_sequence, termination_codon_sequence = prdeict_ORF.sequence_extract(predict_path, threads)
    predict_path.add_path("IRES Sequence", IRES_sequence)
    predict_path.add_path("m6A sequence", m6A_sequence)
    predict_path.add_path("termination codon sequence", termination_codon_sequence)
    IRES_sequence_result, m6A_sequence_result, termination_codon_sequence_result = prdeict_ORF.estimate_circRNA_ORF(predict_path, threads)
    predict_path.add_path("IRES score", IRES_sequence_result)
    predict_path.add_path("m6A score", m6A_sequence_result)
    predict_path.add_path("termination codon score", termination_codon_sequence_result)
    orf_info_path = prdeict_ORF.estimate_info_collect(predict_path, "Ribo")

    create_folder(base_path["result"])

    if args.draw_visualization != "None":
        logger.info("\n" + "="*30 + "\nDraw a visual display")
        draw_path = PathManager(base_path.paths)
        draw_path.add_path("m6A score", m6A_sequence_result)
        draw_path.add_path("orf_info_path", orf_info_path)
        drawcirc.Ribo_draw_circ(draw_path, threads, True if args.draw_visualization == "Longest" else False )

    copyfile(orf_info_path, base_path["result"])
    copyfile(m6A_sequence_result, base_path["result"])
    copyfile(coding_reads_file, base_path["result"])

    if not args.retain_temp_file:
        del_all_file(base_path["tmp"])


def sub_command_mass_spectrum_process(args):
    threads = args.threads

    base_path = PathManager()
    base_path.add_path("circRNA", args.circRNA)
    base_path.add_path("output", project_name(args.output))
    base_path.add_path("result", os.path.join(base_path["output"], "Result"))
    base_path.add_path("junction", args.junction)
    base_path.add_path("peptide", args.clear_peptide if args.clear_peptide else args.peptide)
    create_folder(base_path["output"])
    mass_spectrum_process( args, base_path, threads)


def mass_spectrum_process(args, base_path:PathManager, threads):
    logger.info("\n" + "="*50 + "\nStart processing mass spectrometry data\n" + "="*50)
    base_path.add_path("tmp", os.path.join(base_path["output"], "MSTemp"))
    create_folder(base_path["tmp"])
    logger.info("\n" + "="*30 + "\nMerge duplicate peptide segments with statistical comparison information")

    base_path["peptide"] = MS_format_adapter.format_adapt(base_path)
    clear_pepide = peptide_match.filter_coding_circ_peptide(base_path)

    logger.info("\n" + "="*30 + "\nSearching and evaluating ORF information of CircRNA")
    predict_path = PathManager(base_path.paths)
    predict_path.add_path("coding_frag", clear_pepide)
    peptide_info, true_circRNA = prdeict_ORF.mult_thread_orf_fillter_MS(predict_path, threads)
    predict_path.add_path("circRNA info", peptide_info)
    predict_path.add_path("Coding circRNA Sequence", true_circRNA)
    IRES_sequence, m6A_sequence, orf_sequence, termination_codon_sequence = prdeict_ORF.sequence_extract(predict_path, threads)
    predict_path.add_path("IRES Sequence", IRES_sequence)
    predict_path.add_path("m6A sequence", m6A_sequence)
    predict_path.add_path("termination codon sequence", termination_codon_sequence)
    IRES_sequence_result, m6A_sequence_result, termination_codon_sequence_result = prdeict_ORF.estimate_circRNA_ORF(predict_path, threads)
    predict_path.add_path("IRES score", IRES_sequence_result)
    predict_path.add_path("m6A score", m6A_sequence_result)
    predict_path.add_path("termination codon score", termination_codon_sequence_result)
    orf_info_path = prdeict_ORF.estimate_info_collect(predict_path, "MS")

    create_folder(base_path["result"])

    if args.draw_visualization != "None":
        logger.info("\n" + "="*30 + "\nDraw a visual display")
        draw_path = PathManager(base_path.paths)
        draw_path.add_path("m6A score", m6A_sequence_result)
        draw_path.add_path("orf_info_path", orf_info_path)
        drawcirc.MS_draw_circ(draw_path, threads, True if args.draw_visualization == "Longest" else False )

    copyfile(orf_info_path, base_path["result"])
    copyfile(m6A_sequence_result, base_path["result"])
    copyfile(clear_pepide, base_path["result"])

    if not args.retain_temp_file:
        del_all_file(base_path["tmp"])


def sub_command_both_process(args):

    threads = args.threads
    base_path = PathManager()
    base_path.add_path("output", project_name(args.output))
    base_path.add_path("result", os.path.join(base_path["output"], "Result"))

    Ribo_path = PathManager(base_path.paths)
    Ribo_path.add_path("circRNA", args.circRNA)
    Ribo_path.add_path("transcript", args.transcript)
    Ribo_path.add_path("ribosome", args.ribosome)
    Ribo_path.add_path("reads", args.single if args.single else args.paired)

    MS_path = PathManager(base_path.paths)
    MS_path.add_path("circRNA", args.circRNA)
    MS_path.add_path("junction", args.junction)
    MS_path.add_path("peptide", args.clear_peptide if args.clear_peptide else MS_format_adapter.format_adapt(args.peptide))
    create_folder(base_path["output"])

    Ribo_seq_process(args, Ribo_path, threads)
    mass_spectrum_process(args, MS_path, threads)

    logger.info("\n" + "="*30 +  "\nMerge mass spectrometry data with Ribo seq data.")
    MS_orf_info_path = os.path.join(MS_path.get_path("result"), "MS_ORF_evaluate.xlsx")
    Ribo_orf_info_path = os.path.join(Ribo_path.get_path("result"), "Ribo_ORF_evaluate.xlsx")
    Summary_orf_info_path = os.path.join(Ribo_path.get_path("result"), "Summary_Table_MS_Ribo.xlsx")
    MS_orf_info = pd.read_excel(MS_orf_info_path, engine="openpyxl")
    Ribo_orf_info = pd.read_excel(Ribo_orf_info_path, engine="openpyxl")
    merge_MS_Ribo = pd.merge(MS_orf_info, Ribo_orf_info, on=["circRNA_name", "ORF_Start", "ORF_Stop", "IRES_Score", "m6A_Score", "Termination_Codon_Score", "cORF"], how="outer")
    merge_MS_Ribo = merge_MS_Ribo.drop(["reads_Start", "reads_Stop", "peptide_Start", "peptide_total_sequence"], axis=1)
    merge_MS_Ribo = merge_MS_Ribo.reindex(labels=["circRNA_name", "Reads_number", "peptide_number", "ORF_Start", "ORF_Stop", "IRES_Score", "m6A_Score", "Termination_Codon_Score", "cORF"], axis=1)
    merge_MS_Ribo.to_excel(Summary_orf_info_path, index=False)
    Ribo_m6A_path = Ribo_path.get_path("m6A score")
    MS_m6A_path = MS_path.get_path("m6A score")
    result_m6A_path = os.path.join(base_path.get_path("result"), "m6a_sequence_file.txt")
    mrege_m6A(MS_m6A_path, Ribo_m6A_path, result_m6A_path)
    logger.info("Program running completed.")


def mrege_m6A(m6A1, m6A2, result):
    m6A1_context = {}
    with open(m6A1, "r") as fi:
        text = fi.readlines()
        for i in range(0, len(text), 3):
            m6A1_context[text[i][1:].strip().split("\t")[0]] = text[i] + text[i+1] + text[i+2]
    m6A2_context = {}
    with open(m6A2, "r") as fi:
        text = fi.readlines()
        for i in range(0, len(text), 3):
            m6A2_context[text[i][1:].strip().split("\t")[0]] = text[i] + text[i+1] + text[i+2]
    m6A1_context.update(m6A2_context)
    with open(result, "a+") as fi:
        fi.writelines(m6A1_context.values())


def sub_command_junction_process(args):
    circRNA = args.circRNA
    outfile = args.output
    logger.info("Start generating junction region protein sequence files")
    junction.junctions(circRNA, outfile)
    logger.info("Generation completed")


def get_resource_path(package, file):
    with importlib.resources.path(package, file) as p:
        return p

def sub_command_DLMSC_process(args):
    if not args.model_path:
        model_path = get_resource_path("CircCode3.resources.stop_codon_model","checkpoint_ATCG.pth.tar")
    else:
        model_path = args.model_path
    predict_fa, outfile = args.predict_fa,  args.output
    DLMSC(predict_fa, model_path, outfile)


def sub_command_DeepCircm6A_process(args):
    if not args.model_path:
        model_path = get_resource_path("CircCode3.resources.DeepCircm6A","checkpoint.pth.tar")
    else:
        model_path = args.model_path

    file, mode, outfile, thread = args.predict_fa, args.mode, args.output, args.threads
    DeepCircm6a_score(file, mode, model_path, outfile, thread)


def sub_command_Draw_process(args):
    threads = args.threads
    draw_path = PathManager()
    draw_path.add_path("circRNA", args.circRNA)
    draw_path.add_path("result", args.output)
    draw_path.add_path("m6A score", args.m6A_file)
    draw_path.add_path("orf_info_path", args.orf_info)
    create_folder(draw_path["result"])
    if args.type == "MS":
        drawcirc.MS_draw_circ(draw_path, threads, True if args.draw_visualization == "Longest" else False )
    elif args.type == "Ribo":
        drawcirc.Ribo_draw_circ(draw_path, threads, True if args.draw_visualization == "Longest" else False )



def command_main(args):
    title = """
_________ .__               _________            .___     ________  
\_   ___ \|__|______   ____ \_   ___ \  ____   __| _/____ \_____  \ 
/    \  \/|  \_  __ \_/ ___\/    \  \/ /  _ \ / __ |/ __ \  _(__  < 
\     \___|  ||  | \/\  \___\     \___(  <_> ) /_/ \  ___/ /       \\
 \______  /__||__|    \___  >\______  /\____/\____ |\___  >______  /
        \/                \/        \/            \/    \/       \/ 
                        """
    tips = """There's nothing here. \nYou can use CircCode3 -h or choose a mode and use the -h parameter to get help."""
    print(title, tips)