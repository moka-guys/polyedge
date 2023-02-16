""" Config for polyedge.py so as not to clog up script
"""
import os

# General settings
APP_NAME = "moka-guys/polyedge"
TITLE = "POLYEDGE REPORT"

TABLE_HEADERS = [
    "Gene", "Chrom", "Position", "First base of poly repeat allele", "Read count",
    "Mean quality of first base", "Fraction of reads", "Mean poly length",
    "Standard deviation of poly length", "Mode of poly length", "Average purity of polyN repeat"
    ]

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))  # Project root
TEMPLATE_DIR = os.path.join(ROOT_DIR, 'templates')
LOGOPATH = os.path.join(ROOT_DIR, "images/logo.png")

# HTML/PDF-specific settings
HTML_TEMPLATE = os.path.join(TEMPLATE_DIR, 'template_report.html')
HTML_ROW = "<tr>{}{}{}{}{}{}{}{}{}{}{}</tr>"
CELL = '<td>{} </td>'
CELL_PASS = '<td style="background-color: #90EE90">{} </td>'
CELL_FAIL = '<td style="background-color: #e77e7e">{} </td>'
LIST = "<ul><li>{} {}</li><li>{} {}</li><li>{} {}</li><li>{} {}</li><li>{} {}</li><li>{} {}</li><li>{} {}</li></ul>"

THRESHOLDS = {
    "read_count": 100,
    "mean_quality": 20,
    "read_fraction": 0.2,
    "poly_purity": 0.95,
    }

INTERPRETATION_THRESHOLDS = [
    "Interpretation Thresholds:", "Read count >= 100X", "Mean BQ >= 20",
    "Het if both fractions of reads >=0.2",
    "Purity of polyN should be >=0.95 (max 5% variation within polyN repeat)"
    ]

PARAMETERS_STR = "The app was run with the following parameters:"

BAM_PARAM = "--bam"
BAI_PARAM = "--bai"
GENE_PARAM = "--gene"
CHROM_PARAM = "--chrom"
POLY_S_PARAM = "--poly_start"
POLY_E_PARAM = "--poly_end"
ANCHOR_LEN_PARAM = "--anchor_length"

# FOOTER_CONTENT = "BAM file: {}\n"\
#                  "BAI file: {}\n"\
#                  "Input parameters: --gene {}, --chrom {}, "\
#                  "--poly_start {}, --poly_end {}, --anchor_length {}"
