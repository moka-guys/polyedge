# PolyEdge

Finds variants at the edge of a poly that is is varying length.

It is currently hard-coded to find the varying base at the 5'-end of intron 5 in MSH2.

## Usage 
Simply run with an indexed BAM file as first argument, like:
```
python polyedge.py your_bam_file.bam
```

The result is a table written to STDOUT, like this:
```
A  33 204 0.38 23.52 2.62 25
T  20 328 0.62 25.51 2.34 27
|  |  |   |    |     |    +- Mode of poly length
|  |  |   |    |     +- Standard deviation of poly length
|  |  |   |    +- Mean poly length
|  |  |   +- Fraction of reads
|  |  +- Read Count
|  +- Mean quality if first base
+- First base of poly repeat allele
```
