from argparse import ArgumentParser

from concurrent.futures import ProcessPoolExecutor
import pysam, re, subprocess, logging
# from Bio import SeqIO

from spoa import poa
from itertools import groupby


logger = logging.getLogger()


def _enable_logging(log_file, debug, overwrite):
    """
    Turns on logging, sets debug levels and assigns a log file
    """
    log_formatter = logging.Formatter("[%(asctime)s] %(name)s: %(levelname)s: "
                                      "%(message)s", "%Y-%m-%d %H:%M:%S")
    console_formatter = logging.Formatter("[%(asctime)s] %(levelname)s: "
                                          "%(message)s", "%Y-%m-%d %H:%M:%S")

    console_log = logging.StreamHandler()
    console_log.setFormatter(console_formatter)
    if not debug:
        console_log.setLevel(logging.INFO)

    if overwrite:
        open(log_file, "w").close()

    file_handler = logging.FileHandler(log_file, mode="a")
    file_handler.setFormatter(log_formatter)

    logger.setLevel(logging.DEBUG)
    logger.addHandler(console_log)
    logger.addHandler(file_handler)


class Insertion(object):
    def __init__(self, chrom, start, length, haplotype, seq):
        self.chrom = chrom
        self.start = start
        self.end = start + length
        self.length = length
        self.haplotype = haplotype
        self.seq = seq
        self.merged = False
        self.count = 1

    def __lt__(self, other):
        """
        Class method for sorting our inserts based on haplotype, chromosome and reference coordinate
        return a list with inserts from lowest haplotype, then chromosome number to highest,
        and then sorted based on reference coordinate
        """
        return (self.haplotype < other.haplotype) \
            or (self.chrom < other.chrom) \
            or (self.chrom == other.chrom and self.start < other.start)

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

    log_file = "STRdust.log"

    _enable_logging(log_file, args.debug, overwrite=False)

    dust = {}
    for chrom in pysam.AlignmentFile(args.bam, "rb").references:
        insertions = extract_insertions(args.bam, chrom, minlen=15,
                                        mapq=10, merge_distance=args.distance)
        insertions = merge_overlapping_insertions(sorted(insertions), merge_distance=args.distance)
        get_merged_ins_file(insertions, args.ins_file)
        mreps_dict = parse_mreps_result(run_mreps(args.ins_file, args.mreps_res))
        for key in mreps_dict.keys():
            dust[key] = mreps_dict[key]
        # Group those insertions that are at approximately the same location and the same haplotype
        # Create a consensus out of those by simple counting or local assembly
        # Assess if an insertion is repetitive (mreps?) and extract the unit motif

    vcfy(dust, args.outfile)


def extract_insertions(bamf, chrom, minlen, mapq, merge_distance):
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

    logging.info("Start extraction of insertions and softclips")

    insertions = []
    bam = pysam.AlignmentFile(bamf, "rb")
    for read in bam.fetch(contig=chrom, multiple_iterators=True):
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
                                      seq=read.query_sequence[read_position:read_position + length])
                        )
                    read_position += length

        if len(insertions_per_read) > 1:
            insertions_per_read = horizontal_merge(insertions_per_read, merge_distance)
        insertions.extend(insertions_per_read)

    logging.info("End with extraction of insertions and softclips.")
    return insertions


def get_haplotype(read):
    """Return the haplotype to which the read is assigned
    Or 'un' for reads that are unphased"""
    return str(read.get_tag('HP')) if read.has_tag('HP') else 'un'


# PLEASE REVIEW FUNCTION BELOW
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
            # logging.info(f"{i} {len(to_merge)}")
            cons_ins = create_consensus(to_merge)

            if cons_ins is not None:
                merged_insertions.append(cons_ins)

            to_merge = []

    logging.info("End with merging overlapping insertions")
    return merged_insertions


def create_consensus(insertions_to_merge, max_ins_length=7500):
    logging.info(f"Start merging insertions {len(insertions_to_merge)}")

    if len(insertions_to_merge) == 1:
        return insertions_to_merge[0]
    else:
        count = len(insertions_to_merge)
        # logging.info("Start assembling things.")

        length_above_cutoff = [len(i.seq) > max_ins_length for i in insertions_to_merge]
        if any(length_above_cutoff):
            return None

        # logging.info(f"Seq lengths: {[len(i.seq) for i in insertions_to_merge]}")

        consensus_seq = assemble([i.seq for i in insertions_to_merge])
        # logging.info("End assembling things")
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
    consensus, _ = poa(seqs, algorithm=1, m=2, n=-4, g=-4, e=-2, q=-24, c=-1)
    return consensus

def get_merged_ins_file(insertions, file_name):
    """
    Output an fasta file contains all insertions for mreps
    """
    with open(file_name, "w") as ins_file:
        for i in insertions:
            ins_file.writelines(f">{i.chrom}_{int(i.start)}_{int(i.end)}\n{i.seq}\n")

def run_mreps(file_name, mreps_res):
    mreps_result = subprocess.run(["mreps", "-fasta", "-res", str(mreps_res), file_name], capture_output=True)
    return(mreps_result.stdout.decode("utf-8"))

def parse_mreps_result(mreps_output_str):
    """
    Input: the output string from mreps
    Output: a dictionary that contains
    key: chr22_start_end (the location of the insertion in chromosome)
    value:[[start, end, string], [start, end, string], ...]
        start is the start position of the repeat
        end is the end position of the repeat
        string is the repeat string
    """
    mreps_split_str = ' ---------------------------------------------------------------------------------------------'

    result_dict = {}
    mreps_output_str = re.split("Processing sequence", mreps_output_str)[1:]
    for output_str in  mreps_output_str:
        if "RESULTS: There are no repeats in the processed sequence" in output_str:
            continue
        else:
            output_list = output_str.split('\n')
            ins_loc = output_list[0]
            temp = []
            all_repeat_info = [list(g) for k, g in groupby(output_list, key=lambda x: x != mreps_split_str) if k][1]
            for info in all_repeat_info:
                info_list = info.split("\t")
                loc_list = re.findall(r'\d+', info_list[0])
                seq = info_list[-1].split()[0]
                ins_info = loc_list + [seq]
                temp.append(ins_info)
            result_dict.update({ins_loc:temp})
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
    strdust_vcf.write("#chrom\tstart\tend\trepeat_seq\tsize\n")

    for dustspec in mrep_dict.keys():
        [chrom, start_ins, end_ins] = dustspec.split("'")[1].split("_")
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
                    strdust_vcf.write("%s\t%s\t%s\t%s\t%s\n" % (chrom, str(start_mrep+start_ins), str(end_mrep+start_ins), seq, str(end_mrep-start_mrep)))
        else:
            [[start_mrep, end_mrep, seq]] = mrep_dict[dustspec]
            start_mrep = int(start_mrep)
            end_mrep = int(end_mrep)
            # skip homopolymers
            if len(seq) > 1:
                strdust_vcf.write("%s\t%s\t%s\t%s\t%s\n" % (chrom, str(start_mrep+start_ins), str(end_mrep+start_ins), seq, str(end_mrep-start_mrep)))

    strdust_vcf.close()


def get_args():
    parser = ArgumentParser("Genotype STRs from long reads")
    parser.add_argument("bam", help="phased bam file")
    parser.add_argument("-d", "--distance",
                        help="distance across which two events should be merged",
                        type=int,
                        default=50)
    parser.add_argument("-f", "--ins_file",
                        help="location of merged insertion consensus fasta file",
                        type=str,
                        default="./ins.fa")
    parser.add_argument("-r", "--mreps_res",
                        help="tolerent error rate in mreps repeat finding",
                        type=int,
                        default=1)
    parser.add_argument("-o", "--outfile",
                        help="name of the vcf output",
                        type=str,
                        default="strdust-list.vcf")
    parser.add_argument("--debug", action="store_true",
                        dest="debug", default=False,
                        help="enable debug output")
    return parser.parse_args()


if __name__ == '__main__':
    main()
