import pysam
import sys
import re
import statistics
from collections import defaultdict
from functools import reduce

CHROM = "2"
POS = 47641559
POLY = (POS, 47641586)
ANCHOR_LENGTH = 2


def main(bamfile):
    # list of poly lenght counts by first base allele
    counts, quals = defaultdict(list), defaultdict(list)
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
            seq_start = next(x for x in read.get_aligned_pairs()
                             if x[1] == roi_start)
            seq_end = next(x for x in read.get_aligned_pairs()
                           if x[1] == roi_end)
            # if anchored
            if seq_start[0] and seq_end[0]:
                seq_segment = read.query_sequence[seq_start[0]:seq_end[0]]
                first_base_pos = seq_start[0]+ANCHOR_LENGTH
                first_base_qual = read.query_qualities[first_base_pos]
                # extract poly sequence assuming no indel in anchor sequence)
                regex = r'.{' + re.escape(str(ANCHOR_LENGTH)) + \
                    r'}(.+).{' + re.escape(str(ANCHOR_LENGTH)) + r'}'
                seq_match = re.match(regex, seq_segment)
                if seq_match:
                    # partition poly by first base allele
                    allele = seq_match.group(1)[0]
                    counts[allele].append(len(seq_match.group(1)))
                    quals[allele].append(first_base_qual)
    # calculate statistics
    total_count = reduce(lambda a, c: a + len(c), counts.values(), 0)
    for allele in sorted(counts.keys()):
        data = (allele,
                len(counts[allele]),
                int(statistics.mean(quals[allele])),
                len(counts[allele])/total_count,
                statistics.mean(counts[allele]),
                statistics.stdev(counts[allele]),
                statistics.mode(counts[allele]))
        print('{:>} {:4n} {:2n} {:.2f} {:2.2f} {:.2f} {:2n}'.format(*data))


if __name__ == "__main__":
    if len(sys.argv) == 3:
        ANCHOR_LENGTH = int(sys.argv[2])
    main(sys.argv[1])
