from argparse import ArgumentParser
import pysam


class Insertion(object):
    def __init__(self, chrom, start, length, seq, type):
        self.chrom = chrom
        self.start = start
        self.end = start + length
        self.seq = seq
        self.type = type

    def __lt__(self, other):
        """
        Class method for sorting our inserts based on chromosome and reference coordinate
        return a list with inserts from lowest chromosome number to highest,
        and then sorted based on reference coordinate
        """
        return self.chrom < other.chrom or (self.chrom == other.chrom and self.start < other.start)

    def __eq__(self, other, distance=15):
        """
        Class method used for comparing 2 inserts (with some wobble distance)
        """
        condition = [self.chrom == other.chrom, self.ref_coords[0] + distance > other.ref_coords[0]]
        return all(condition)


def main():
    args = get_args()
    insertions = extract_insertions(args.bam, minlen=15, mapq=10)
    # Group those insertions that are at approximately the same location and the same haplotype
    # Create a consensus out of those by simple counting or local assembly
    # Assess if an insertion is repetitive (mreps?) and extract the unit motif


def extract_insertions(bamf, minlen, mapq):
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
    for read in bam.fetch():
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
                                      seq=read.query_sequence[read_position:read_position + length],
                                      type="SOFTCLIP" if operation == 4 else "INS")
                        )
                    read_position += length

        if len(insertions_per_read) != 0:
            insertions_per_read.sort()
            # TODO Merging of inserts that are close together on the same read
            insertions.extend(insertions_per_read)
    return insertions


def get_haplotype(read):
    """Return the haplotype to which the read is assigned
    Or 'un' for reads that are unphased"""
    return read.get_tag('HP') if read.has_tag('HP') else 'un'


def get_args():
    parser = ArgumentParser("Genotype STRs from long reads")
    parser.add_argument("bam", help="phased bam file")
    return parser.parse_args()


if __name__ == '__main__':
    main()
