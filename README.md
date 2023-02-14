# PolyEdge

Finds variants at the edge of a poly that is varying length.

## Usage

The script takes the following required non-optional command line arguments:
```
  -B BAM, --bam BAM     Bam file to analyse
  -I BAI, --bai BAI     Bam index file
  -S POLY_START, --poly_start POLY_START
                        Start position of poly stretch
  -E POLY_END, --poly_end POLY_END
                        End position of poly stretch
  -A ANCHOR_LENGTH, --anchor_length ANCHOR_LENGTH
                        Length of anchor sequence
  -C CHROM, --chrom CHROM
                        Chromosome of interest
```

 run with an indexed BAM file as first argument, like:
For example:
```
python polyedge.py your_bam_file.bam
```

NB: An optional second argument can be used to modify the default anchor length of 2 (number of aligned bases required on either side of the repeat). It is _not_ recommended to change this.

## Output

The result is a table written to STDOUT, like this:
```
A 204 33 0.38 23.52 2.62 25 1.00
T 328 20 0.62 25.51 2.34 27 0.99
|  |  |   |    |     |    |  +- Average purity of polyN repeat
|  |  |   |    |     |    +- Mode of poly length
|  |  |   |    |     +- Standard deviation of poly length
|  |  |   |    +- Mean poly length
|  |  |   +- Fraction of reads
|  |  +- Mean quality of first base
|  +- Read Count
+- First base of poly repeat allele
```
