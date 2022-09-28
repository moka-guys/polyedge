import pysam
import sys
import re
import statistics
from collections import defaultdict, Counter
from functools import reduce

CHROM = "2"
POS = 47641559
POLY = (POS, 47641586)
ANCHOR_LENGTH = 2


def main(bamfile):
    # list of poly lenght counts by first base allele
    counts = defaultdict(list)  # lengths of repeat given allele (one per read)
    quals = defaultdict(list)  # quality of first base
    purity = defaultdict(list)  # most common base count in poly w/o first base
    # open BAM file
    reader = pysam.Samfile(bamfile, 'rb')
    # fetch reads intersecting the poly repeat
    roi_reads = reader.fetch(CHROM, *POLY)
    for read in roi_reads:
        # define start and end of segment to be matched
        roi_start = POLY[0]-ANCHOR_LENGTH
        roi_end = POLY[1]+ANCHOR_LENGTH
        # select that span poly with defined ANCHOR sequence length
        if read.reference_start < roi_start and roi_end < read.reference_end:
            # find sequence segment with poly and anchor sequence
            # -> first matching position in the query->reference mapping
            seq_start = next(x for x in read.get_aligned_pairs()
                             if x[1] == roi_start)
            seq_end = next(x for x in read.get_aligned_pairs()
                           if x[1] == roi_end)
            # if anchored
            if seq_start[0] and seq_end[0]:
                seq_segment = read.query_sequence[seq_start[0]:seq_end[0]]
                first_base_pos = seq_start[0]+ANCHOR_LENGTH
                first_base_qual = read.query_qualities[first_base_pos]
                # extract poly sequence (excluding anchor bases)
                # assuming no indel in anchor sequence
                regex = r'.{' + re.escape(str(ANCHOR_LENGTH)) + \
                    r'}(.+).{' + re.escape(str(ANCHOR_LENGTH)) + r'}'
                seq_match = re.match(regex, seq_segment)
                if seq_match:
                    # partition poly by first base allele
                    allele = seq_match.group(1)[0]
                    counts[allele].append(len(seq_match.group(1)))
                    # count polyN purity
                    remainder_poly = Counter(seq_match.group(1)[1:])
                    poly_most_abundant = remainder_poly.most_common()[0][1]
                    purity[allele].append(poly_most_abundant)
                    # quality of first base of ROI (allele in question)
                    quals[allele].append(first_base_qual)
    # calculate statistics
    total_count = reduce(lambda a, c: a + len(c), counts.values(), 0)
    for allele in sorted(counts.keys()):
        data = (allele,
                len(counts[allele]),
                int(statistics.mean(quals[allele])),
                len(counts[allele])/total_count,
                statistics.mean(counts[allele]),
                statistics.stdev(counts[allele]) if len(counts[allele]) > 1 else 0,
                statistics.mode(counts[allele]),
                sum(purity[allele])/(sum(counts[allele])-len(counts[allele]))
                )
        print('{:>} {:4n} {:2n} {:.2f} {:2.2f} {:.2f} {:2n} {:.2f}'
              .format(*data))


if __name__ == "__main__":
    if len(sys.argv) == 3:
        ANCHOR_LENGTH = int(sys.argv[2])
    main(sys.argv[1])
