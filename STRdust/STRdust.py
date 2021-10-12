from argparse import ArgumentParser
import pysam


class Insertion(object):
    def __init__(self, chrom, start, length, haplotype, seq):
        self.chrom = chrom
        self.start = start
        self.end = start + length
        self.haplotype = haplotype
        self.seq = seq

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
    for chrom in pysam.AlignmentFile(args.bam, "rb").references:
        insertions = extract_insertions(args.bam, chrom, minlen=15,
                                        mapq=10, merge_distance=args.distance)
        insertions = merge_overlapping_insertions(sorted(insertions), merge_distance=args.distance)
        # Group those insertions that are at approximately the same location and the same haplotype
        # Create a consensus out of those by simple counting or local assembly
        # Assess if an insertion is repetitive (mreps?) and extract the unit motif


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
    return insertions


def get_haplotype(read):
    """Return the haplotype to which the read is assigned
    Or 'un' for reads that are unphased"""
    return read.get_tag('HP') if read.has_tag('HP') else 'un'


# PLEASE REVIEW FUNCTION BELOW
def horizontal_merge(insertions, merge_distance):
    """Merge insertions occuring in the same read if they are within merge_distance"""
    insertions.sort()
    while True:
        distances = [insertions[n]-insertions[n-1] for n in range(1, len(insertions))]
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


def get_args():
    parser = ArgumentParser("Genotype STRs from long reads")
    parser.add_argument("bam", help="phased bam file")
    parser.add_argument("-d", "--distance",
                        help="distance across which two events should be merged",
                        type=int,
                        default=50)
    return parser.parse_args()


if __name__ == '__main__':
    main()
