# PolyEdge

Finds variants at the 5' end of a poly stretch of varying length in a BAM file

This provides a solution to the alignment issues caused by poly repeat stretches, whereby the aligner spreads the variant across adjacent bases and dilutes it to a point where the standard variant caller cannot detect it.

## Inputs

The script takes the following required non-optional command line arguments:

```bash
  -B BAM,         --bam BAM                 Bam file to analyse
  -I BAI,         --bai BAI                 Bam index file
  -G GENE,        --gene GENE               Gene of interest
  -S POLY_START,  --poly_start POLY_START   Start position of poly stretch (0-based)
  -E POLY_END,    --poly_end POLY_END       End position of poly stretch (0-based)
  -C CHROM,       --chrom CHROM             Chromosome of interest
```

The following arguments are optional:

```bash
  -A ANCHOR_LENGTH, --anchor_length ANCHOR_LENGTH Length of anchor sequence
```

NB: Anchor_length (number of aligned bases required on either side of the repeat) of 2 is recommended, and is set as default. It is _not_ recommended to change this.

## Usage

The script can be run as follows:

```python
python polyedge.py -B your_bam_file.bam -I your_bam_file_index.bai -G MSH2 -S 47641559 -E 47641586 -C 2
```

## Output

The script outputs several files which describe the alleles seen at the position of interest (this position is the `poly_start` position, and whilst the script input is a 0-based coordinate, it is displayed in the output files as a 1-based coordinate).

* CSV file containing raw data - `$SAMPLENAME.refined.$GENE_polyedge.csv`
* HTML report - `$SAMPLENAME.refined.$GENE_polyedge.html`
* PDF report - `$SAMPLENAME.refined.$GENE_polyedge.pdf`

The HTML and PDF reports are identical.

The outputs contain tables with the calculated metrics required to interpret the variant. The HTML and PDF reports contain a list of interpretation thresholds. These thresholds are conservative and therefore can be applied to any case.

Columns containing metrics that have specified interpretation thresholds (Read count, mean quality of first base, fraction of reads, average purity of polyN repeat) are highlighted in green or red dependent upon whether they meet or fail to meet the specified threshold.

## Docker image

The docker image can be can be built, tagged and saved as a .tar.gz file using the Makefile:

```bash
sudo make build
```

The docker image can be run as follows:

```python

sudo docker run -v $PATH_TO_INPUTS:$PATH_TO_INPUTS -v $PATH_TO_OUTPUTS:/outputs seglh-polyedge:v1.1.0 -B $PATH_TO_INPUTS/NGS506B_96_286171_DH_M_VCP2R211Via_Pan4130_S68_R1_001.refined.bam -I $PATH_TO_INPUTS/NGS506B_96_286171_DH_M_VCP2R211Via_Pan4130_S68_R1_001.refined.bam.bai -G MSH2 -S 47641559 -E 47641586 -C 2


python polyedge.py -B your_bam_file.bam -I your_bam_file_index.bai -G MSH2 -S 47641559 -E 47641586 -C 2
```

N.B. $PATH_TO_INPUTS must be the same inside and outside the docker image
$PATH_TO_OUTPUTS should be wherever you want the outputs stored on your local system


The current and all previous verisons of the tool are stored as dockerised versions in 001_ToolsReferenceData project as .tar.gz files.

### This app was produced by the Synnovis Genome Informatics Team
