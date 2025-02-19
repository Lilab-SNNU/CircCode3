# _*_ coding : utf-8 _*_
import math
import multiprocessing
from typing import Sequence
import pandas as pd
from Bio import SeqIO
from multiprocessing import Value
from PIL import Image, ImageDraw, ImageFont


from ..utils.logs import logger
from ..utils.path_utils import *


class DrawCirc:
    """
    Base class for drawing;
    The function consists of two parts:
    1. Function used to draw line segment diagrams;
    2. Functions for drawing text graphs
    3. Text function for adding legends
    """

    def __init__(self):
        self.size = (5000, 5000)
        self.myseq = Image.new('RGB', self.size, 'white')
        self.DRAW = ImageDraw.Draw(self.myseq)
        self.CENTER = (2500, 3000)                      # 弧形的中心
        self.dpi = (300, 300)
        self.arial100 = ImageFont.truetype('src/resources/LiberationSans-Regular.ttf', 100)

    def get_angle(self, bp: int, length: int) -> float:
        """
        bp position and Angle conversion
        :param bp: sequence position
        :param length: total length of the sequence
        :return: The angle corresponding to the float bp position
        """
        return bp * 360 / length

    def coord(self, angle: float, center: tuple[int, int], radius: float) -> tuple[int, int]:
        """
        Obtain x and y coordinates from angles
        :param angle: angle
        :param center: Center
        :param radius: radius
        :return: float, x, y coordinates corresponding to float
        """
        rad = math.radians(angle)
        x = int(center[0] + math.cos(rad) * radius)
        y = int(center[1] + math.sin(rad) * radius)
        return x, y

    def draw_circle(self, radius_i: float, radius_o: float, start_angle: float, stop_angle: float, color1: str, color2: str):
        """
        Draw an arc formed by two sectors covering each other. Therefore, when drawing, the outer circle should be drawn first, followed by the inner circle.

        :param radius_i: internal sector radius .
        :param radius_o: external sector radius.
        :param start_angle: start angle
        :param stop_angle: stop angle
        :param color1: The colors used to fill the outer layers
        :param color2: The colors used to fill the inner layers
        :return:
        """
        x1 = self.CENTER[0] - radius_o
        y1 = self.CENTER[1] - radius_o
        x2 = self.CENTER[0] + radius_o
        y2 = self.CENTER[1] + radius_o
        self.DRAW.pieslice((x1, y1, x2, y2), start_angle, stop_angle, fill=color1)

        x3 = self.CENTER[0] - radius_i
        y3 = self.CENTER[1] - radius_i
        x4 = self.CENTER[0] + radius_i
        y4 = self.CENTER[1] + radius_i
        self.DRAW.pieslice((x3, y3, x4, y4), start_angle, stop_angle, fill=color2)

    def draw_arrow_tip(self, start: float, direction: float, radius_mid: float, color: str):
        """
        Draw a triangular arrow

        :param start: starting position
        :param direction: length
        :param radius_mid: End radius
        :param color: Fill color
        :return:
        """
        p1 = self.coord(start + direction, self.CENTER, radius_mid)
        p2 = self.coord(start, self.CENTER, radius_mid - 50)
        p3 = self.coord(start, self.CENTER, radius_mid + 50)
        self.DRAW.polygon((p1, p2, p3), fill=color)

    def draw_helix(self, radius_i: float, radius_o: float, start_angle: float, stop_angle: float, color1: str, color2: str):
        """
        Draw a spiral arc
        The spiral arc is composed of multiple segments of arcs spliced together, with each inward translation of 0.2 * step pixels.
        The step parameter controls the length of each arc, reducing the amount of detail that can be drawn and increasing it can speed up the drawing process.

        :param radius_i: internal sector radius .
        :param radius_o: external sector radius.
        :param start_angle: start angle
        :param stop_angle: stop angle
        :param color1: The colors used to fill the outer layers
        :param color2: The colors used to fill the inner layers
        :return:
        """
        step = 3                        # 调整绘图速度,与精度
        for i in range(int(start_angle/step), int(stop_angle/step)):
            j = 0.2 * step * i
            x1 = self.CENTER[0] - radius_o + j
            y1 = self.CENTER[1] - radius_o + j
            x2 = self.CENTER[0] + radius_o - j
            y2 = self.CENTER[1] + radius_o - j

            x3 = self.CENTER[0] - radius_i + j
            y3 = self.CENTER[1] - radius_i + j
            x4 = self.CENTER[0] + radius_i - j
            y4 = self.CENTER[1] + radius_i - j
            self.DRAW.pieslice((x1, y1, x2, y2), step*i, step*(i+1), fill=color1)
            self.DRAW.pieslice((x3, y3, x4, y4), step*i, step*(i+1), fill=color2)

    def draw_legend(self, site: Sequence[float], text: str, color1: str, color2: str):
        """
        Used to draw legends
        :param site: Location
        :param text: indicates the legend text
        :param color1: indicates the color of the legend
        :param color2: indicates the text color
        :return:
        """
        self.DRAW.rectangle(site, fill=color1)
        self.DRAW.text((site[0] + 300, site[1]), text, fill=color2, font=self.arial100)

    def draw_seq_text(self, i: int, base: str, color: str, length: int, radius: float):
        """
         Used to specify arc position to draw text

        :param i: bp position
        :param base: base text
        :param color: indicates the text color
        :param length: indicates the total length of the sequence
        :param radius: indicates the radius of an arc
        :return:
        """
        angle = self.get_angle(i, length)-90
        p = self.coord(angle, self.CENTER, radius)
        self.DRAW.text(p, base, fill=color, font=self.arial100)

    def draw_circle_text(self, sequence: str, start: int, end: int, radius: float, length: int, color: str):
        '''
        Used to draw an arc of text
        :param sequence: A sequence of text that needs to be drawn
        :param start: indicates the start position bp
        :param end: indicates the end position bp
        :param radius: Radius
        :param length: indicates the total length of the sequence
        :param color: sequence text color
        :return:
        '''
        for i in range(start, end):
            self.draw_seq_text(i, sequence[i], color, length, radius)


class DrawCircrna(DrawCirc):
    """
    Used to draw circrnas represented by lines
    """
    def __init__(self, circ_name: str, start: int, end: int, length: int, ires_score: float, m6a_score: float, endcodon_score: float, reads_count: int):
        super().__init__()
        self.circ_name = circ_name
        self.start = start
        self.end = end
        self.ires_score = round(ires_score, 3)
        self.m6a_score = round(m6a_score, 3)
        self.endcodon_score = round(endcodon_score, 3)
        self.length = length
        self.color = ("#002c53", "#ff6600", "#00994e", "#e30039", "#ffc000", "#99cc00", "#D8E6E7")
        self.reads_count = reads_count

    def draw_ires_circ(self):
        """
        Used to draw ires sequences
        Includes two parts of arc and arrow, wherein arc has included ires full length
        :return:
        """
        start = self.get_angle(self.start-101, self.length) - 90
        end = self.get_angle(self.start, self.length) - 90
        self.draw_circle(1600, 1650, start, end, self.color[1], "white")
        self.draw_arrow_tip(end, 10, 1625, self.color[1])

    def draw_orf_circ(self):
        """
        For ORF rendering, in order to be compatible with ORF with different number of turns, spiral rendering is used.
        :return:
        """
        start = self.get_angle(self.start, self.length) - 90
        end = self.get_angle(self.end, self.length) - 90
        self.draw_helix(1250, 1300, start, end, self.color[4], "white")

    def draw_m6a_circ(self, m6a_site: Sequence[int]):
        """
        Loop through the m6a list and plot the m6a sites.
        :param m6a_site: list stores the m6a location
        :return:
        """
        for i in m6a_site:
            start = self.get_angle(i, self.length)-90
            end = self.get_angle(i+1, self.length)-90
            end = max(end, start+1)                     # 最小占1度
            self.draw_circle(1450, 1500, start, end, self.color[2], "white")

    def draw_seq(self):
        """
        Draw the entire sequence as a circle with one notch (2 degrees).
        :return:
        """
        self.draw_circle(1450, 1500, -89, 269, self.color[6], "white")  # 3

    def draw_reads(self, reads_site: tuple[int]):
        """
        Plot the total length of the matched reads or peptides
        :param reads_site: indicates the start and end positions of the total length
        :return:
        """
        start = self.get_angle( reads_site[0], self.length) - 90
        end = self.get_angle(reads_site[1], self.length) - 90
        self.draw_circle(1700, 1750, start, end, self.color[0], "white")

    def draw_legends(self):
        """
        Draw legend
        :return:
        """
        self.draw_legend((3200, 300, 3400, 400), "CircRNA Sequence", self.color[6], "black")
        self.draw_legend((3200, 500, 3400, 600), "ORF Sequence", self.color[4], "black")
        self.draw_legend((3200, 700, 3400, 800), f"IRES Score : {self.ires_score}", self.color[1], "black")
        self.draw_legend((3200, 900, 3400, 1000), f"m6A Site Count : {self.m6a_score}", self.color[2], "black")
        self.draw_legend((3200, 1100, 3400, 1200), "Translation Support  Sequence", self.color[0], "black")

    def draw_text(self):
        """
        Draw the text information that needs to be represented
        :return:
        """
        self.DRAW.text((300,300), f'circRNA_lenth = {self.length}', fill='black', font=self.arial100)
        self.DRAW.text((300,500), f'ORF_lenth = {self.end-self.start}', fill='black', font=self.arial100)
        self.DRAW.text((300,700), f'ORF start site : {self.start}', fill='black', font=self.arial100)
        self.DRAW.text((300,900), f'ORF stop site : {self.end%self.length}', fill='black', font=self.arial100)
        self.DRAW.text((300,1100), f'reads count : {self.reads_count}', fill='black', font=self.arial100)
        # 标注3‘ 和5’ 在序列绘制那，会影响其他圈绘制
        self.DRAW.text((2420, 1400), "3'", fill='black', font=self.arial100)
        self.DRAW.text((2540, 1400), "5'", fill='black', font=self.arial100)


class DrawTextCircrna(DrawCirc):
    """
    Used to draw CircRNA diagrams in text form
    """
    def __init__(self, circ_name: str, start: int, end: int, sequence: str, length: int):
        super().__init__()
        self.circ_name = circ_name
        self.start = start
        self.end = end
        self.length = length
        self.sequence = sequence
        self.color = ("#002c53", "#ff6600", "#00994e", "#e30039", "#ffc000", "#99cc00", "#D8E6E7")
        self.arial100 = ImageFont.truetype('src/resources/LiberationSans-Regular.ttf', int(10000 / length))

    def draw_sequence(self):
        """
        Draw the full sequence, leaving a gap so the first base is not shown.
        :return:
        """
        self.draw_circle_text(self.sequence, 1, self.length, 1700, self.length, self.color[6])
        self.draw_circle(1700-int(10000/self.length), 1700, -90, -90+(100/self.length), self.color[6], "white")


    def draw_orf_text(self, translate_sequence: str):
        """
        The ORF sequence is drawn and the translated sequence is drawn inside it, with the translated protein corresponding to the middle bit of the triple codon.
        :param translate_sequence: translated sequence
        :return:
        """
        amino_acid = iter(translate_sequence)
        for i in range(self.start, self.end):
            site = i%self.length
            length = 200*i/self.length
            self.draw_seq_text(site, self.sequence[site], self.color[4], self.length, 1650-length)
            if (i-self.start+2)%3 == 0:
                self.draw_seq_text(site, next(amino_acid), self.color[5], self.length, 1570-length)

    def draw_ires_text(self):
        """
        Draw IRES text
        :return:
        """
        self.draw_circle_text(self.sequence, self.start-101, self.start, 1800, self.length, self.color[3])

    def draw_reads_text(self, site: list[int]):
        """
        Draws the matching reads text
        :param site: List of start and end locations
        :return:
        """
        self.draw_circle_text(self.sequence, site[0], site[1], 1900, self.length, self.color[0])

    def draw_peptide_text(self, site: list[int], petide: str):
        """
        Draw the matching peptide text
        :param site: List of start and end locations
        :param petide: str peptide text, corresponding to the first digit of the triple codon
        :return:
        """
        amino_acid = iter(petide)
        for i in range(site[0], site[1], 3):
            self.draw_seq_text(i, next(amino_acid), self.color[0], self.length, 1900)

    def draw_m6a_text(self, m6a_site: list):
        """
        Loop through the m6a list and plot the m6a sites.
        :param m6a_site: list Stores the m6a location
        :return:
        """
        for i in m6a_site:
            self.draw_seq_text(i, self.sequence[i], self.color[2], self.length, 1700)

    def draw_legends(self):
        """
        Draw legend
        """
        self.arial100 = ImageFont.truetype('src/resources/LiberationSans-Regular.ttf', 100)
        self.draw_legend((200, 300, 300, 400), "Translation Support Sequence", self.color[0], "black")
        self.draw_legend((200, 500, 300, 600), "ORF Sequence", self.color[4], "black")
        self.draw_legend((200, 700, 300, 800), "IRES Sequence", self.color[3], "black")
        self.draw_legend((200, 900, 300, 1000), "m6A Modification Sites", self.color[2], "black")
        self.draw_legend((200, 1100, 300, 1200), "Translate Sequence", self.color[5], "black")
        self.draw_legend((200, 1300, 300, 1400), "CircRNA Sequence", self.color[6], "black")


def MS_liner_circ_picture(orf_info_path: str, dict_m6a: dict[str, list[int]], dict_sequence: dict[str, str], outfile: str, threads: int, filter_flag: bool):
    """
    Plot a linear representation of mass spectrum data

    :param orf_info_path: orf info path
    :param dict_m6a: Dictionary of circRNA names and m6A locations
    :param dict_sequence:Dictionary of circRNA names and sequence
    :param outfile: Output directory
    :param threads: Number of processes used
    :param filter_flag: Whether to draw only the longest orf
    :return:
    """
    orf_info = pd.read_excel(orf_info_path, engine="openpyxl")
    if filter_flag:
        orf_info = filter_max_length(orf_info)
    total_number = orf_info.shape[0]
    pool = multiprocessing.Pool(threads)
    for row in orf_info.itertuples():
        name = row.circRNA_name
        start = row.ORF_Start
        end = row.ORF_Stop
        ires_score = row.IRES_Score
        m6a_score = row.m6A_Count
        endcodon_score = row.Termination_Codon_Score
        reads_num = row.peptide_number
        length = len(dict_sequence[name])

        index_name = name + f"|({start},{end})"
        m6a_site = dict_m6a[index_name]
        reads_site = (row.peptide_Start, row.peptide_Start + len(row.peptide_total_sequence)*3)
        # 实例化绘图
        pool.apply_async(ms_draw_liner_picture,
                         (name, start, end, length, ires_score, m6a_score, endcodon_score, reads_num, reads_site,
                          m6a_site, outfile, total_number), error_callback=error_break)
    pool.close()
    pool.join()
    print("\n")




def ms_draw_liner_picture(name: str, start: int, end: int, length: int, ires_score: float, m6a_score: float, endcodon_score: float, reads_num: int, peptide_site: tuple[int], m6a_site: list[int], outresult: str, total_number: int):
    """
    The drawing is separate for multiprocess encapsulation
    :param name:
    :param start:
    :param end:
    :param length:
    :param ires_score:
    :param m6a_score:
    :param endcodon_score:
    :param reads_num:
    :param peptide_site:
    :param m6a_site:
    :param outresult:
    :param total_number:
    :return:
    """
    draw = DrawCircrna(name, start, end, length, ires_score, m6a_score, endcodon_score, reads_num)
    outname = name + f"-({start},{end})"
    name = name + f"|({start},{end})"
    draw.draw_reads(peptide_site)
    draw.draw_ires_circ()
    draw.draw_seq()
    draw.draw_m6a_circ(m6a_site)
    draw.draw_orf_circ()
    draw.draw_legends()
    draw.draw_text()
    draw.DRAW.text((500, 4700), name, fill='black', font=draw.arial100)
    draw.myseq.save('{}/{}.pdf'.format(outresult, outname), "PDF", dpi=draw.dpi)

    num.value += 1
    print("\r", "{:.2f}%".format(round(num.value*100/total_number, 2)), end="", flush=True)


def Ribo_liner_circ_picture(orf_info_path, dict_m6a, dict_sequence, outfile, threads, filter_flag):
    """
    Plot a linear representation of Ribo-seq data

    :param orf_info_path: orf info path
    :param dict_m6a: Dictionary of circRNA names and m6A locations
    :param dict_sequence:Dictionary of circRNA names and sequence
    :param outfile: Output directory
    :param threads: Number of processes used
    :param filter_flag: Whether to draw only the longest orf
    :return:
    """

    orf_info = pd.read_excel(orf_info_path, engine="openpyxl")
    if filter_flag:
        orf_info = filter_max_length(orf_info)
    total_number = orf_info.shape[0]

    pool = multiprocessing.Pool(threads)
    for row in orf_info.itertuples():
        name = row.circRNA_name
        start = row.ORF_Start
        end = row.ORF_Stop
        ires_score = row.IRES_Score
        m6a_score = row.m6A_Count
        endcodon_score = row.Termination_Codon_Score
        reads_num = row.Reads_number
        length = len(dict_sequence[name])

        index_name = name + f"|({start},{end})"
        m6a_site = dict_m6a[index_name]
        reads_site = (row.reads_Start, row.reads_Stop)
        # 实例化绘图
        pool.apply_async(ribo_draw_liner_picture,
                         (name, start, end, length, ires_score, m6a_score, endcodon_score, reads_num, reads_site,
                          m6a_site, outfile, total_number), error_callback=error_break)
    pool.close()
    pool.join()
    print("\n")



def ribo_draw_liner_picture(name: str, start: int, end: int, length: int, ires_score: float, m6a_score: float, endcodon_score: float, reads_num: int, reads_site: tuple[int], m6a_site: list[int], outresult: str, total_number: int):
    """
    The drawing is separate for multiprocess encapsulation
    :param name:
    :param start:
    :param end:
    :param length:
    :param ires_score:
    :param m6a_score:
    :param endcodon_score:
    :param reads_num:
    :param reads_site:
    :param m6a_site:
    :param outresult:
    :param total_number:
    :return:
    """
    draw = DrawCircrna(name, start, end, length, ires_score, m6a_score, endcodon_score, reads_num)
    outname = name + f"-({start},{end})"
    name = name + f"|({start},{end})"
    draw.draw_reads(reads_site)
    draw.draw_ires_circ()
    draw.draw_seq()
    draw.draw_m6a_circ(m6a_site)
    draw.draw_orf_circ()
    draw.draw_legends()
    draw.draw_text()
    draw.DRAW.text((500, 4700), name, fill='black', font=draw.arial100)

    draw.myseq.save('{}/{}.pdf'.format(outresult, outname), "PDF", dpi=draw.dpi)
    num.value += 1
    print("\r", "{:.2f}%".format(round(num.value*100/total_number, 2)), end="", flush=True)


def sequence_dict(sequence: str) -> dict[str, str]:
    """
    Gets a dictionary of sequences and sequence names
    :param sequence: indicates the original fasta file address
    :return: Dictionary of names and sequences
    """

    Seqs = SeqIO.parse(sequence, "fasta")
    dict_seq = {}
    for Seq in Seqs:
        dict_seq[Seq.id] = Seq.seq
    return dict_seq


def m6a_site_encode(m6a_mapping_file: str) -> dict[str, list[int]]:
    """
    Position m6a. The threshold is set to 0,5
    :param m6a_mapping_file: Deepcircm6a output location file path
    :return: Names with m6a modified list of locations
    """
    with open(m6a_mapping_file, "r") as fi:
        context = fi.readlines()
    encode_dict = {}
    for i in range(0, len(context), 3):
        splice_name = context[i].split("\t")[0][1::]
        name = splice_name
        site = eval(splice_name.rsplit("|", 1)[1])
        scores = context[i+2].strip().split("\t")
        iter_sequence = iter(scores)                # 拿迭代器返回列表
        encode = []
        try:
            for j in range(site[0]-126, site[0]+25):
                if float(next(iter_sequence)) > 0.5:
                    encode.append(j)
        except StopIteration:
            raise "There are items in the m6a file that are less than 150bp, please check the m6a file length."
        encode_dict[name] = encode

    return encode_dict


def MS_text_circ_picture(orf_info_path, dict_m6a, dict_sequence, outfile, threads, filter_flag):
    """
    Draw a Text image of the MS result
    :param orf_info_path: orf info path
    :param dict_m6a: Dictionary of circRNA names and m6A locations
    :param dict_sequence:Dictionary of circRNA names and sequence
    :param outfile: Output directory
    :param threads: Number of processes used
    :param filter_flag: Whether to draw only the longest orf
    :return:
    """
    pool = multiprocessing.Pool(threads)

    orf_info = pd.read_excel(orf_info_path, engine="openpyxl")
    if filter_flag:
        orf_info = filter_max_length(orf_info)
    total_number = orf_info.shape[0]

    for row in orf_info.itertuples():
        name = row.circRNA_name
        start = row.ORF_Start
        end = row.ORF_Stop
        corf = row.cORF
        sequence = dict_sequence[name]
        peptide = row.peptide_total_sequence
        length = len(dict_sequence[name])
        index_name = name + f"|({start},{end})"
        m6a_site = dict_m6a[index_name]
        reads_site = (row.peptide_Start, row.peptide_Start + len(peptide)*3)
        pool.apply_async(ms_draw_text_picture,
                         (name, start, end, sequence, length, corf, m6a_site, reads_site, peptide, outfile, total_number), error_callback=error_break)
    pool.close()
    pool.join()
    print("\n")


def ms_draw_text_picture(name, start, end, sequence, length, corf, m6a_site, peptide_site, peptide, outfile, total_number):
    draw = DrawTextCircrna(name, start, end, sequence, length)
    outname = name + f"-({start},{end})"
    draw.draw_sequence()
    draw.draw_orf_text(corf)
    draw.draw_ires_text()
    draw.draw_m6a_text(m6a_site)
    draw.draw_peptide_text(peptide_site, peptide)
    draw.draw_legends()
    draw.myseq.save('{}/{}.pdf'.format(outfile, outname), "PDF", dpi=draw.dpi)

    num.value += 1
    print("\r", "{:.2f}%".format(round(num.value*100/total_number, 2)), end="", flush=True)


def Ribo_text_circ_picture(orf_info_path, dict_m6a, dict_sequence, outfile, threads, filter_flag):
    """
    Draw a Text image of the Ribo-seq result
    :param orf_info_path: orf info path
    :param dict_m6a: Dictionary of circRNA names and m6A locations
    :param dict_sequence:Dictionary of circRNA names and sequence
    :param outfile: Output directory
    :param threads: Number of processes used
    :param filter_flag: Whether to draw only the longest orf
    :return:
    """

    pool = multiprocessing.Pool(threads)

    orf_info = pd.read_excel(orf_info_path, engine="openpyxl")
    if filter_flag:
        orf_info = filter_max_length(orf_info)
    total_number = orf_info.shape[0]

    for row in orf_info.itertuples():
        name = row.circRNA_name
        start = row.ORF_Start
        end = row.ORF_Stop
        corf = row.cORF
        sequence = dict_sequence[name]
        length = len(dict_sequence[name])

        index_name = name + f"|({start},{end})"
        m6a_site = dict_m6a[index_name]
        reads_site = (row.reads_Start, row.reads_Stop)
        # 实例化绘图
        pool.apply_async(ribo_draw_text_picture,
                         (name, start, end, sequence, length, corf, m6a_site, reads_site, outfile, total_number))
    pool.close()
    pool.join()
    print("\n")


def ribo_draw_text_picture(name, start, end, sequence, length, corf, m6a_site, reads_site, outfile, total_number):
    draw = DrawTextCircrna(name, start, end, sequence, length)
    outname = name + f"-({start},{end})"
    draw.draw_sequence()
    draw.draw_orf_text(corf)
    draw.draw_ires_text()
    draw.draw_m6a_text(m6a_site)
    draw.draw_reads_text(reads_site)
    draw.draw_legends()
    draw.myseq.save('{}/{}.pdf'.format(outfile, outname), "PDF", dpi=draw.dpi)
    num.value += 1
    print("\r", "{:.2f}%".format(round(num.value*100/total_number, 2)), end="", flush=True)


def MS_draw_circ(draw_path:PathManager, threads, filter_flag=False):
    """
    Mass spectrum data structure diagram and text diagram drawing
    :param draw_path: A PathManager object that includes basic parameters
    :param threads: The number of threads used for drawing
    :param filter_flag: Whether to draw only the longest orf
    :return:
    """
    orf_info = draw_path.get_path("orf_info_path")
    m6A_sequence_result = draw_path.get_path("m6A score")
    circRNA = draw_path.get_path("circRNA")

    MS_CIRC_graph = os.path.join(draw_path.get_path("result"), "MS_CIRC_graph")
    create_folder(MS_CIRC_graph)
    MS_Text_graph = os.path.join(draw_path.get_path("result"), "MS_TEXT_graph")
    create_folder(MS_Text_graph)

    dict_m6a = m6a_site_encode(m6A_sequence_result)
    dict_sequence = sequence_dict(circRNA)

    global num
    logger.info("Draw the linear graph")
    num = Value("i", 0)
    MS_liner_circ_picture(orf_info, dict_m6a, dict_sequence, MS_CIRC_graph, threads, filter_flag)
    logger.info("Draw the text graph")
    num = Value("i", 0)
    MS_text_circ_picture(orf_info, dict_m6a, dict_sequence, MS_Text_graph, threads, filter_flag)
    logger.info("Visual image drawing completed")


def Ribo_draw_circ(draw_path:PathManager, threads, filter_flag=False):
    """
    Ribo-seq data structure diagram and text diagram drawing
    :param draw_path: A PathManager object that includes basic parameters
    :param threads: The number of threads used for drawing
    :param filter_flag: Whether to draw only the longest orf
    :return:
    """
    orf_info_path = draw_path.get_path("orf_info_path")
    m6A_sequence_result = draw_path.get_path("m6A score")
    circRNA = draw_path.get_path("circRNA")

    Ribo_CIRC_graph = os.path.join(draw_path.get_path("result"), "Ribo_CIRC_graph")
    create_folder(Ribo_CIRC_graph)
    Ribo_Text_graph = os.path.join(draw_path.get_path("result"), "Ribo_TEXT_graph")
    create_folder(Ribo_Text_graph)

    dict_m6a = m6a_site_encode(m6A_sequence_result)                   # m6a位置
    dict_sequence = sequence_dict(circRNA)

    global num

    logger.info("Draw the linear graph")
    num = Value("i", 0)
    Ribo_liner_circ_picture(orf_info_path, dict_m6a, dict_sequence, Ribo_CIRC_graph, threads, filter_flag)
    num = Value("i", 0)
    logger.info("Draw the text graph")
    Ribo_text_circ_picture(orf_info_path, dict_m6a, dict_sequence, Ribo_Text_graph, threads, filter_flag)
    logger.info("Visual image drawing completed")


def error_break(err):
    logger.error(err)


def filter_max_length(data:pd.DataFrame):
    """
    Filter the longest ORF of each CircRNA
    :param data:
    :return:
    """
    data["length"] = data["ORF_Start"] - data["ORF_Stop"]
    result = data.loc[data.groupby("circRNA_name")["length"].idxmax()]
    result = result.drop("length", axis=1)
    return result