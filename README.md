# PolyEdge

Finds variants at the edge of a poly that is varying length.

## Inputs

The script takes the following required non-optional command line arguments:

```
  -B BAM,         --bam BAM                 Bam file to analyse
  -I BAI,         --bai BAI                 Bam index file
  -G GENE,        --gene GENE               Gene of interest
  -S POLY_START,  --poly_start POLY_START   Start position of poly stretch
  -E POLY_END,    --poly_end POLY_END       End position of poly stretch
  -C CHROM,       --chrom CHROM             Chromosome of interest
```

The following arguments are optional:

```
  -A ANCHOR_LENGTH, --anchor_length ANCHOR_LENGTH Length of anchor sequence
```

NB: Anchor_length (number of aligned bases required on either side of the repeat) of 2 is recommended, and is set as default. It is _not_ recommended to change this.

## Usage

The script can be run as follows:

```python
python polyedge.py -B your_bam_file.bam -I your_bam_file_index.bai -G MSH2 -S 47641559 -E 47641586 -C 2
```

## Output

The script outputs several files:
* CSV file containing raw data
* HTML report
* PDF report

The HTML and PDF reports are identical.

The outputs contain tables with the calculated metrics required to interpret the variant. The HTML and PDF reports contain a list of interpretation thresholds. Columns containing metrics that have specified interpretation thresholds (Read count, mean quality of first base, fraction of reads, average purity of polyN repeat) are highlighted in green or red dependent upon whether they meet or fail to meet the specified threshold.

## Docker image

A docker image can be built, tagged and saved as a .tar.gz file using the following commands:

```
sudo docker build -f Dockerfile .
sudo docker tag $IMAGE_ID seglh-polyedge:$VERSION
sudo docker save seglh-polyedge:$VERSION | gzip > seglh-polyedge:$VERSION.tar.gz
```

The current and all previous verisons of the tool are stored as dockerised versions in 001_ToolsReferenceData project as .tar.gz files.

### This app was produced by the Synnovis Genome Informatics Team
