""" Config for polyedge.py so as not to clog up script
"""
import os

# General settings
APP_NAME = "moka-guys/polyedge"
TITLE = "POLYEDGE REPORT"
ROOT_DIR = os.path.dirname(os.path.abspath(__file__))  # Project root
TEMPLATE_DIR = os.path.join(ROOT_DIR, "templates")
LOGOPATH = os.path.join(ROOT_DIR, "images/logo.png")

TABLE_HEADERS = [
    "Gene",
    "Chrom",
    "Position",
    "First base of poly repeat allele",
    "Read count",
    "Mean quality of first base",
    "Fraction of reads",
    "Mean poly length",
    "Standard deviation of poly length",
    "Mode of poly length",
    "Average purity of polyN repeat",
]

# HTML/PDF-specific settings
HTML_TEMPLATE = os.path.join(TEMPLATE_DIR, "template_report.html")
HTML_TBL_HEADER = (
    "<tr><th>{}</th><th>{}</th><th>{}</th><th>{}</th><th>{}</th><th>{}</th>"
    "<th>{}</th><th>{}</th><th>{}</th><th>{}</th><th>{}</th></tr>"
)

HTML_TBL_ROW = "<tr>{}{}{}{}{}{}{}{}{}{}{}</tr>"
HTML_TBL_CELL = "<td>{} </td>"
HTML_TBL_CELL_PASS = '<td style="background-color: #90EE90">{} </td>'
HTML_TBL_CELL_FAIL = '<td style="background-color: #e77e7e">{} </td>'
HTML_LIST = (
    "<ul><li>{} {}</li><li>{} {}</li><li>{} {}</li><li>{} {}</li>"
    "<li>{} {}</li><li>{} {}</li><li>{} {}</li></ul>"
)

THRESHOLDS = {
    "read_count": 100,
    "mean_quality": 20,
    "read_fraction": 0.2,
    "poly_purity": 0.95,
}

INTERP_THRESHS = {
    "thresh_title": "Interpretation Thresholds",
    "thresh_str": "The above metrics are highlighted in green if they meet the below "
    "thresholds and in red if they do not.",
    "read_count_thresh": "Read count per allele >= 100X",
    "mean_bq_thresh": "Mean BQ >= 20",
    "fraction_reads_thresh": "Het if both fractions of reads >=0.2",
    "polyn_pur_thresh": "Purity of polyN should be >=0.95 (max 5% variation within polyN repeat)",
}

PARAMS_STR = "The app was run with the following parameters:"

PARAMS = {
    "bam": "--bam",
    "bai": "--bai",
    "gene": "--gene",
    "chrom": "--chrom",
    "poly_start": "--poly_start",
    "poly_end": "--poly_end",
    "anchor_length": "--anchor_length",
}
