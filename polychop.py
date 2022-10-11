import pysam

CHROM = "2"
POLY = (47641559, 47641586)

'''Chop reads that span a predefined sequence (eg. polyN repeat)
preserving either up-stream sequence up to last aligned base of the polyN stretch'''

def main():
    # open BAM file
    reader = pysam.AlignmentFile('-', 'rb')
    writer = pysam.AlignmentFile('-', "wb", template=reader)
    # fetch reads intersecting the poly repeat
    for read in reader:
        # select that span poly with defined ANCHOR sequence length
        try:
            left = read.reference_start < POLY[0]
            rite = POLY[1] < read.reference_end
            assert read.reference_name == CHROM and left and rite
        except (AssertionError, TypeError):
            writer.write(read)
            continue
        except:
            raise

        # find cut position (last aligned base within poly)
        cut_pos = next(x for x in read.get_aligned_pairs()[::-1]
                        if x[0] and x[1] and x[1] <= POLY[1])

        # cut read
        if cut_pos[0]:
            # trim sequence and qualities
            quals = read.query_qualities  # pysam invalidates query_qualities when changing query_sequence
            read.query_sequence = read.query_sequence[:cut_pos[0]]
            read.query_qualities = quals[:cut_pos[0]]

            assert len(read.query_sequence) == len(read.query_qualities)

            # correct cigar string
            newcigar = []
            for c in read.cigartuples:
                # current query position
                qlen = sum([ x[1] for x in newcigar if x[0] in (0,1,4,7,8) ])
                if qlen < cut_pos[0]:
                    if c[0] in [0,1,4,7,8]:  # M, I, S, =, X operations (query consuming)
                        newtuple = (c[0], min(c[1], cut_pos[0]-qlen))
                        newcigar += [newtuple]
                    else:
                        newcigar += [c]
                else:
                    break
            read.cigartuples = newcigar

        # write read
        writer.write(read)

if __name__ == "__main__":
    main()
