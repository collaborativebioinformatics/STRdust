from argparse import ArgumentParser
import sys
import os
import shutil
from concurrent.futures import ProcessPoolExecutor
import pysam
import re
import subprocess
import logging
import pandas as pd
from .version import __version__

from spoa import poa
from itertools import groupby, repeat

import STRdust.utils as utils


class Insertion(object):
    def __init__(self, chrom, start, length, haplotype, seq):
        self.chrom = chrom
        self.start = start
        self.end = start + length
        self.length = length
        self.haplotype = haplotype
        self.seq = seq
        self.count = 1

    def __lt__(self, other):
        """
        Class method for sorting our inserts based on haplotype, chromosome and reference coordinate
        return a list with inserts from lowest haplotype, then chromosome number to highest,
        and then sorted based on reference coordinate
        """
        return (self.haplotype < other.haplotype) \
            or (self.chrom < other.chrom) \
            or (self.haplotype == other.haplotype and self.chrom == other.chrom and self.start < other.start)

    def is_overlapping(self, other, distance=15):
        """
        Class method used for comparing 2 inserts (with some wobble distance)
        Assuming sorted input with self < other
        """
        condition = [self.haplotype == other.haplotype,
                     self.chrom == other.chrom,
                     self.start + distance > other.start]
        return all(condition)


def main():
    args = get_args()

    args.out_dir = utils.create_output_directory(args.out_dir)
    utils._enable_logging(args.out_dir, args.debug, overwrite=True)

    try:
        utils._check_bam_files(args.bam)
    except utils.InputException as err:
        logging.error(f"Problem with input files: {err}")

    args.ins_dir, args.vcf_dir = utils.setup_temp_dirs(args.out_dir)

    vcf_final_file = os.path.join(args.out_dir, "strdust-list.vcf")
    if args.region:
        temporary_files = [run(args, args.region)]
    else:
        chromosomes = [c for c in pysam.AlignmentFile(args.bam, "rb").references if '_' not in c]
        with ProcessPoolExecutor(max_workers=args.threads) as executor:
            temporary_files = [f for f in executor.map(run, repeat(args), chromosomes)]

    concatenate_output(temporary_files, vcf_final_file)

    if not args.save_temp:
        logging.info("Cleaning up output directory.")
        shutil.rmtree(args.ins_dir)
        shutil.rmtree(args.vcf_dir)
    logging.info("Enjoy your annotation.")


def run(args, region):
    logging.info(f"-- Start processing: {region} --")
    insertions = extract_insertions(args.bam, region, minlen=15,
                                    mapq=10, merge_distance=args.distance, flank_distance=50)
    insertions = merge_overlapping_insertions(sorted(insertions), merge_distance=args.distance)

    region_string = region.replace(':', '_').replace('-', '_')
    ins_chr_file = os.path.join(args.ins_dir, f"ins_{region_string}.fa")
    write_ins_file(insertions, ins_chr_file)
    if not args.save_temp:
        os.remove(ins_chr_file)

    mreps_dict = parse_mreps_result(run_mreps(ins_chr_file, args.mreps_res))
    vcf_temporary_file = os.path.join(args.vcf_dir, f"strdust-{region_string}.tsv")
    if mreps_dict:
        vcfy(mreps_dict, vcf_temporary_file)
    return vcf_temporary_file


def extract_insertions(bamf, region, minlen, mapq, merge_distance, flank_distance):
    """
    Extract insertions and softclips from a bam file based on parsing CIGAR strings

    Depending on the CIGAR operation, the 'cursor' is moved forward in either the read coordinates,
    the reference coordinates or both ('consumes query' vs 'consumes reference').
    See also https://samtools.github.io/hts-specs/SAMv1.pdf, page8 [on 20211012], CIGAR

    BAM_CMATCH  0
    BAM_CINS    1
    BAM_CDEL    2
    BAM_CREF_SKIP   3
    BAM_CSOFT_CLIP  4
    BAM_CHARD_CLIP  5
    BAM_CPAD    6
    BAM_CEQUAL  7
    BAM_CDIFF   8
    BAM_CBACK   9
    """

    logging.info(f"{region}: Start extraction of insertions and softclips")

    insertions = []
    bam = pysam.AlignmentFile(bamf, "rb")
    for read in bam.fetch(region=region, multiple_iterators=True):
        insertions_per_read = []
        read_position = 0
        reference_position = read.reference_start + 1
        if read.mapping_quality > mapq:
            for operation, length in read.cigartuples:
                if operation in [0, 7, 8]:
                    read_position += length
                    reference_position += length
                elif operation in [2, 3]:
                    reference_position += length
                elif operation in [1, 4]:
                    if length >= minlen:
                        insertions_per_read.append(
                            Insertion(chrom=read.reference_name,
                                      start=reference_position,
                                      length=length,
                                      haplotype=get_haplotype(read),
                                      seq=read.query_sequence[read_position - flank_distance:read_position + length + flank_distance])
                        )
                    read_position += length

        if len(insertions_per_read) > 1:
            insertions_per_read = horizontal_merge(insertions_per_read, merge_distance)
        insertions.extend(insertions_per_read)

    logging.info(f"{region}: End with extraction of insertions and softclips.")
    return insertions


def get_haplotype(read):
    """Return the haplotype to which the read is assigned
    Or 'un' for reads that are unphased"""
    return str(read.get_tag('HP')) if read.has_tag('HP') else 'un'


def horizontal_merge(insertions, merge_distance):
    """Merge insertions occuring in the same read if they are within merge_distance"""
    insertions.sort()
    while True:
        distances = [insertions[n].start-insertions[n-1].start for n in range(1, len(insertions))]
        distance_below_cutoff = [d < merge_distance for d in distances]
        if any(distance_below_cutoff):
            new_ins = []
            skip = False
            for i, m in enumerate(distance_below_cutoff):
                if skip:
                    skip = False
                    continue
                if m:
                    new_ins.append(Insertion(
                        chrom=insertions[i].chrom,
                        start=(insertions[i].start + insertions[i+1].start) / 2,
                        length=insertions[i].length + insertions[i+1].length,
                        haplotype=insertions[i].haplotype,
                        seq=insertions[i].seq + insertions[i+1].seq,
                    ))
                    skip = True  # skip next insertion because we merged that one in the current
                else:
                    new_ins.append(insertions[i])
            insertions = new_ins
        else:
            return insertions


def merge_overlapping_insertions(insertions, merge_distance):
    logging.info("Start with merging overlapping insertions")

    merged_insertions = []
    to_merge = []

    for i in range(len(insertions)):
        to_merge.append(insertions[i])

        if (i == len(insertions) - 1) or (not insertions[i].is_overlapping(insertions[i + 1], distance=merge_distance)):
            logging.debug(f"Merging {i} {len(to_merge)}")
            cons_ins = create_consensus(to_merge)

            if cons_ins is not None:
                merged_insertions.append(cons_ins)

            to_merge = []

    logging.info("End with merging overlapping insertions")
    return merged_insertions


def create_consensus(insertions_to_merge, max_ins_length=7500):
    logging.debug(f"Start merging insertions {len(insertions_to_merge)}")

    if len(insertions_to_merge) == 1:
        return insertions_to_merge[0]
    else:
        count = len(insertions_to_merge)

        length_above_cutoff = [len(i.seq) > max_ins_length for i in insertions_to_merge]
        if any(length_above_cutoff):
            return None

        # logging.info(f"Seq lengths: {[len(i.seq) for i in insertions_to_merge]}")

        consensus_seq = assemble([i.seq for i in insertions_to_merge])

        merged = Insertion(
            chrom=insertions_to_merge[0].chrom,
            start=sum([i.start for i in insertions_to_merge]) / count,
            length=len(consensus_seq),
            haplotype=insertions_to_merge[0].haplotype,
            seq=consensus_seq)
        merged.count = count

        return merged


# REVIEW PARAMETERS
def assemble(seqs):
    """
    Create a consensus of the inserted fragments

    algorithm: 0 - local (Smith-Waterman) 1 - global (Needleman-Wunsch) 2 - semi-global
    m: match score
    n: score for mismatching bases
    g: gap opening penalty
    e: gap extension penalty
    q: gap opening penalty of the second affine function
    c: gap extension penalty of the second affine function
    """
    consensus, _ = poa(seqs, algorithm=1, m=2, n=-4, g=-4, e=-2, q=-24, c=-1)
    return consensus


def write_ins_file(insertions, file_name):
    """
    Output an fasta file contains all insertions for mreps
    """
    with open(file_name, "w") as ins_file:
        for i in insertions:
            ins_file.writelines(f">{i.chrom}_{int(i.start)}_{int(i.end)}\n{i.seq}\n")


def run_mreps(file_name, mreps_res):
    mreps_result = subprocess.run(
        ["mreps", "-fasta", "-res", str(mreps_res), file_name], capture_output=True)
    return(mreps_result.stdout.decode("utf-8"))


def parse_mreps_result(mreps_output_str):
    """
    Input: the output string from mreps
    Output: a dictionary that contains
    key: chrN_start_end (the location of the insertion in chromosome)
    value:[[start, end, string], [start, end, string], ...]
        start is the start position of the repeat
        end is the end position of the repeat
        string is the repeat string
    """
    mreps_split_str = ' ---------------------------------------------------------------------------------------------'

    result_dict = {}
    mreps_output_str = re.split("Processing sequence", mreps_output_str)[1:]
    for output_str in mreps_output_str:
        if "RESULTS: There are no repeats in the processed sequence" in output_str \
                or 'Processed sequence is too short' in output_str:
            continue
        else:
            output_list = output_str.split('\n')
            ins_loc = output_list[0]
            temp = []
            all_repeat_info = [list(g) for k, g in groupby(
                output_list, key=lambda x: x != mreps_split_str) if k][1]
            for info in all_repeat_info:
                info_list = info.split("\t")
                loc_list = re.findall(r'\d+', info_list[0])
                seq_list = info_list[-1].split()
                seq = max(seq_list, key=seq_list.count)
                # seq = info_list[-1].split()[0]
                ins_info = loc_list + [seq]
                temp.append(ins_info)
            result_dict.update({ins_loc: temp})
    return result_dict


def vcfy(mrep_dict, oufvcf):
    """
    Input: a dictionary that contains
    key: chrN_start_end (the location of the insertion in chromosome)
    value:[[start, end, string], [start, end, string], ...]
        start is the start position of the repeat based on the ins.fa
        end is the end position of the repeat based on the ins.fa
        string is the repeat string
    """
    strdust_vcf = open(oufvcf, "w")
    logging.info("Writing results to %s" % oufvcf)
    strdust_vcf.write("chrom\tstart\tend\trepeat_seq\tsize\n")

    for dustspec in mrep_dict.keys():
        try:
            [chrom, start_ins, end_ins] = dustspec.split("'")[1].split("_")
        except ValueError:
            sys.exit(dustspec)
        start_ins = int(start_ins)
        end_ins = int(end_ins)
        # mreps can find more than on repeated seq
        if len(mrep_dict[dustspec]) > 1:
            for eachstr in mrep_dict[dustspec]:
                [start_mrep, end_mrep, seq] = eachstr
                start_mrep = int(start_mrep)
                end_mrep = int(end_mrep)
                # skip homopolymers
                if len(seq) > 1:
                    strdust_vcf.write("%s\t%s\t%s\t%s\t%s\n" % (chrom, str(
                        start_mrep+start_ins), str(end_mrep+start_ins), seq, str(end_mrep-start_mrep)))
        else:
            [[start_mrep, end_mrep, seq]] = mrep_dict[dustspec]
            start_mrep = int(start_mrep)
            end_mrep = int(end_mrep)
            # skip homopolymers
            if len(seq) > 1:
                strdust_vcf.write("%s\t%s\t%s\t%s\t%s\n" % (chrom, str(
                    start_mrep+start_ins), str(end_mrep+start_ins), seq, str(end_mrep-start_mrep)))

    strdust_vcf.close()


def concatenate_output(temporary_files, output_file):
    """
    Concatentate files in temporary_files,and sort by chromosome and start position

    The run function returns just the file name, and the file may not have been created
    when no variants were called. Therefore, checking if file exists before reading.
    """
    pd.concat([pd.read_csv(f, sep="\t") for f in temporary_files if os.path.isfile(f)],
              ignore_index=True) \
        .sort_values(by=['chrom', 'start']) \
        .to_csv(output_file, sep="\t", index=False)


def get_args():
    parser = ArgumentParser("Genotype STRs from long reads")
    parser.add_argument("bam", help="phased bam file")
    parser.add_argument("-o", "--out_dir", help="output directory",
                        type=str, default=os.getcwd())
    parser.add_argument("-d", "--distance",
                        help="distance across which two events should be merged",
                        type=int,
                        default=50)
    parser.add_argument("-r", "--mreps_res",
                        help="tolerent error rate in mreps repeat finding",
                        type=int,
                        default=1)
    parser.add_argument("-t", "--threads",
                        help="number of threads to use",
                        type=int,
                        default=8)
    parser.add_argument("--save_temp", action="store_true",
                        dest="save_temp", default=False,
                        help="enable saving temporary files in output directory")
    parser.add_argument("--debug", action="store_true",
                        dest="debug", default=False,
                        help="enable debug output")
    parser.add_argument("--region", help="run on a specific interval only")
    parser.add_argument("-v", "--version",
                        help="Print version and exit.",
                        action="version",
                        version=f'STRdust {__version__}')

    return parser.parse_args()


if __name__ == '__main__':
    main()
