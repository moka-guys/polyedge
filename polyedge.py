import pysam
import sys
import re
import statistics
from collections import defaultdict
from functools import reduce

CHROM = "2"
POS = 47641559
POLY = (POS, 47641586)
ANCHOR_LENGTH = 3


def main(bamfile):
    counts = defaultdict(list)
    # open BAM file
    reader = pysam.Samfile(bamfile, 'rb')
    roi_reads = reader.fetch(CHROM, *POLY)
    for read in roi_reads:
        roi_start = POLY[0]-ANCHOR_LENGTH
        roi_end = POLY[1]+ANCHOR_LENGTH
        if read.reference_start < roi_start and roi_end < read.reference_end:
            seq_start = next(x for x in read.get_aligned_pairs()
                             if x[1] == roi_start)
            seq_end = next(x for x in read.get_aligned_pairs()
                           if x[1] == roi_end)
            segment = read.query_sequence[seq_start[0]:seq_end[0]]
            # extract poly sequence
            regex = r'[ATGC]{' + re.escape(str(ANCHOR_LENGTH)) + \
                r'}(.+)[ATGC]{' + re.escape(str(ANCHOR_LENGTH)) + r'}'
            m = re.match(regex, segment)
            if m:
                # partition poly by first base allele
                counts[m.group(1)[0]].append(len(m.group(1)))

    total_count = reduce(lambda a, c: a + len(c), counts.values(), 0)
    for allele in sorted(counts.keys()):
        data = (allele, len(counts[allele]), len(counts[allele])/total_count,
                statistics.mean(counts[allele]),
                statistics.stdev(counts[allele]),
                statistics.mode(counts[allele]))
        print('{:>} {:4n} {:.2f} {:.2f} {:.2f} {:2n}'.format(*data))


if __name__ == "__main__":
    if len(sys.argv) == 3:
        ANCHOR_LENGTH = int(sys.argv[2])
    main(sys.argv[1])
