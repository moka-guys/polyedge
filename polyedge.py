""" Finds variants at the edge of a poly that is varying length
"""
import subprocess
import re
import os
import argparse
import statistics
from collections import defaultdict, Counter
from functools import reduce
import csv
import datetime
import jinja2
import pdfkit
import pysam
import config


class PolyEdge:
    """
    Class for processing bam files to determine whether a variant is present at the position before
    a poly stretch, and write to file

        generate_output()
            Call methods required to generate all output files for the sample provided

        roi_reads(self):
            Open bam file and fetch reads from input bam (using pysam) intersecting the poly repeat
            :return roi_reads (obj):    Object containing reads

        calculate_metrics()
            Calculate metrics per allele.
            :return allele_list (list): List of dictionaries, one per allele,
                                        containing calculated metrics
        create_csvfile()
            Write metrics to csv. Copy template csv to output dest, write data to end of file
            :param allele_list (list):  List of dictionaries, one per allele,
                                        containing calculated metrics
        construct_table_html()
            Construct html with formatting dependent on whether metrics meet defined thresholds
            :param allele_list (list):  List of dictionaries, one per allele,
                                        containing calculated metrics
            :return table_html (str):   String containing the table html

        create_htmlfile()
            Write metrics to html
            :param table_html (str):    String containing the table html

        create_pdffile()
            Write metrics to pdf
    """

    def __init__(self, bamfile, bam_index, gene, chrom, poly, anchor_length):
        """
        Constructor for PolyEdge class

            :param bamfile (str):       Bam file to analyse
            :param gene (str) :         Gene of interest
            :param chrom (str):         Chromosome of interest
            :param poly (tuple):        Contains start and end position of poly stretch
            :param anchor_length (int): Length of anchor sequence
            :return data (list):        List of tuples per allele, with each tuple containing the
                                        output stats as defined in the readme
        """
        self.bamfile = bamfile
        self.gene = gene
        self.chrom = chrom
        self.chrchrom = f"chr{self.chrom}"
        self.anchor_length = anchor_length
        self.poly = poly
        self.variant_pos = poly[0] + 1
        # Define start and end of segment to be matched
        self.roi_start = self.poly[0] - self.anchor_length
        self.roi_end = self.poly[1] + self.anchor_length
        # Variables required for creating output files
        self.bamfile_prefix = os.path.splitext(bamfile)[0].split("/")[-1]
        self.baifile_prefix = os.path.splitext(bam_index)[0].split("/")[-1]
        self.sample_name = self.bamfile_prefix.split(".")[0]
        self.output_string = f"{self.bamfile_prefix}.{self.gene}_polyedge"
        self.timestamp = datetime.datetime.now().strftime("%d-%B-%Y %H:%M")
        self.csv_path = os.path.join(os.getcwd(), f"{self.output_string}.csv")
        self.html_path = os.path.join(os.getcwd(), f"{self.output_string}.html")
        self.pdf_path = os.path.join(os.getcwd(), f"{self.output_string}.pdf")
        self.footer_html = config.HTML_LIST.format(
            config.PARAMS["bam"],
            self.bamfile_prefix,
            config.PARAMS["bai"],
            self.baifile_prefix,
            config.PARAMS["gene"],
            self.gene,
            config.PARAMS["chrom"],
            self.chrom,
            config.PARAMS["poly_start"],
            self.poly[0],
            config.PARAMS["poly_end"],
            self.poly[1],
            config.PARAMS["anchor_length"],
            self.anchor_length,
        )

    def generate_output(self):
        """
        Call methods required to generate all output files for the sample provided
        """
        allele_list = self.calculate_metrics()
        self.create_csvfile(allele_list)
        table_html = self.construct_table_html(allele_list)
        self.create_htmlfile(table_html)
        self.create_pdffile()

    def calculate_metrics(self):
        """
        Calculate metrics per allele: first base of poly repeat allele, read count, mean quality
        of first base, fraction of reads, mean poly length, standard deviation of poly length,
        mode of poly length, average purity of polyN repeat

            :return allele_list (list): List of dictionaries, one per allele,
                                        containing calculated metrics
        """
        # Poly length counts by first base allele
        counts = defaultdict(
            list
        )  # Store lengths of repeat for given allele (one per read)
        quals = defaultdict(list)  # Store quality of first base
        purity = defaultdict(
            list
        )  # Store most common base count in poly w/o first base

        for read in self.roi_reads():
            # Select reads that span poly with defined ANCHOR sequence length
            if (
                read.reference_start < self.roi_start
                and self.roi_end < read.reference_end
            ):
                # Find sequence segment with poly and anchor sequence
                # -> first matching position in the query->reference mapping
                seq_start = next(
                    x for x in read.get_aligned_pairs() if x[1] == self.roi_start
                )
                seq_end = next(
                    x for x in read.get_aligned_pairs() if x[1] == self.roi_end
                )
                if seq_start[0] and seq_end[0]:  # If anchored
                    seq_segment = read.query_sequence[seq_start[0] : seq_end[0]]
                    first_base_pos = seq_start[0] + self.anchor_length
                    first_base_qual = read.query_qualities[first_base_pos]
                    # Extract poly sequence (excluding anchor bases)
                    # Assuming no indel in anchor sequence
                    regex = (
                        r".{"
                        + re.escape(str(self.anchor_length))
                        + r"}(.+).{"
                        + re.escape(str(self.anchor_length))
                        + r"}"
                    )
                    seq_match = re.match(regex, seq_segment)
                    if seq_match:
                        # Partition poly by first base allele
                        allele = seq_match.group(1)[0]
                        counts[allele].append(len(seq_match.group(1)))
                        # Count polyN purity
                        remainder_poly = Counter(seq_match.group(1)[1:])
                        poly_most_abundant = remainder_poly.most_common()[0][1]
                        purity[allele].append(poly_most_abundant)
                        # Quality of first base of ROI (allele in question)
                        quals[allele].append(first_base_qual)
        # Calculate statistics
        total_count = reduce(lambda a, c: a + len(c), counts.values(), 0)
        allele_list = []
        for allele in sorted(counts.keys()):
            data = {
                "gene": self.gene,
                "chrom": self.chrom,
                "pos": self.variant_pos,
                "first_base": allele,
                "read_count": len(counts[allele]),
                "mean_quality": int(statistics.mean(quals[allele])),
                "read_fraction": round(len(counts[allele]) / total_count, 2),
                "mean_polylen": round(statistics.mean(counts[allele]), 2),
                "stdev_polylen": round(statistics.stdev(counts[allele]), 2)
                if len(counts[allele]) > 1
                else 0,
                "mode_polylen": statistics.mode(counts[allele]),
                "poly_purity": round(
                    sum(purity[allele]) / (sum(counts[allele]) - len(counts[allele]))
                ),
            }
            allele_list.append(data)
        return allele_list

    def roi_reads(self):
        """
        Open bam file and fetch reads from input bam (using pysam) intersecting the poly repeat

            :return roi_reads (obj):    Object containing reads
        """
        reader = pysam.Samfile(self.bamfile, "rb")  # Open BAM file
        try:
            roi_reads = reader.fetch(str(self.chrom), *self.poly)
        except ValueError as exception:
            print(
                f"Using CHROM {str(self.chrchrom)} raised an exception ({exception}). "
                f"Trying: CHROM '{self.chrchrom}'"
            )
            roi_reads = reader.fetch(self.chrchrom, *self.poly)
        return roi_reads

    def create_csvfile(self, allele_list):
        """
        Write metrics to csv file

            :param allele_list (list):  List of dictionaries, one per allele,
                                        containing calculated metrics
        """
        with open(self.csv_path, "wt", encoding="utf-8") as file:
            writer = csv.writer(file, delimiter=",")
            writer.writerow(["Created by", f"{config.APP_NAME}"])
            writer.writerow(["Version", git_tag()])
            writer.writerow(["Date", self.timestamp])
            writer.writerow(["Sample", self.sample_name])
            writer.writerow([])
            writer.writerow([config.PARAMS_STR])
            writer.writerow([config.PARAMS["bam"], self.bamfile_prefix])
            writer.writerow([config.PARAMS["bai"], self.baifile_prefix])
            writer.writerow([config.PARAMS["gene"], self.gene])
            writer.writerow([config.PARAMS["chrom"], self.chrom])
            writer.writerow([config.PARAMS["poly_start"], self.poly[0]])
            writer.writerow([config.PARAMS["poly_end"], self.poly[1]])
            writer.writerow([config.PARAMS["anchor_length"], self.anchor_length])
            writer.writerow([])
            writer.writerow(config.TABLE_HEADERS)

            for row in allele_list:
                writer.writerow(row.values())
        file.close()

    def construct_table_html(self, allele_list):
        """
        Construct html with formatting dependent on whether metrics meet defined thresholds.
        (Formatting is defined in the config file as variables).

            :param allele_list (list):  List of dictionaries, one per allele,
                                        containing calculated metrics
            :return table_html (str):   String containing the table html
        """
        # Generate the table body, with a row per allele
        for allele_dict in allele_list:
            if allele_dict["read_count"] >= config.THRESHOLDS["read_count"]:
                read_count = config.HTML_TBL_CELL_PASS.format(allele_dict["read_count"])
            else:
                read_count = config.HTML_TBL_CELL_FAIL.format(allele_dict["read_count"])

            if allele_dict["mean_quality"] >= config.THRESHOLDS["mean_quality"]:
                mean_quality = config.HTML_TBL_CELL_PASS.format(
                    allele_dict["mean_quality"]
                )
            else:
                mean_quality = config.HTML_TBL_CELL_FAIL.format(
                    allele_dict["mean_quality"]
                )

            if allele_dict["read_fraction"] >= config.THRESHOLDS["read_fraction"]:
                read_fraction = config.HTML_TBL_CELL_PASS.format(
                    allele_dict["read_fraction"]
                )
            else:
                read_fraction = config.HTML_TBL_CELL_FAIL.format(
                    allele_dict["read_fraction"]
                )

            if allele_dict["poly_purity"] >= config.THRESHOLDS["poly_purity"]:
                poly_purity = config.HTML_TBL_CELL_PASS.format(
                    allele_dict["poly_purity"]
                )
            else:
                poly_purity = config.HTML_TBL_CELL_FAIL.format(
                    allele_dict["poly_purity"]
                )

            table_body = config.HTML_TBL_ROW.format(
                config.HTML_TBL_CELL.format(allele_dict["gene"]),
                config.HTML_TBL_CELL.format(allele_dict["chrom"]),
                config.HTML_TBL_CELL.format(allele_dict["pos"]),
                config.HTML_TBL_CELL.format(allele_dict["first_base"]),
                read_count,
                mean_quality,
                read_fraction,
                config.HTML_TBL_CELL.format(allele_dict["mean_polylen"]),
                config.HTML_TBL_CELL.format(allele_dict["stdev_polylen"]),
                config.HTML_TBL_CELL.format(allele_dict["mode_polylen"]),
                poly_purity,
            )
        table_header = config.HTML_TBL_HEADER.format(
            *config.TABLE_HEADERS
        )  # Generate table headers
        table_html = f"{table_header}{table_body}"

        return table_html

    def create_htmlfile(self, table_html):
        """
        Write metrics to html. Load template using jinja2, then write to new file, filling the
        place holders in the template with the specified placeholder values

            :param table_html (str): String containing the table html
        """
        template = jinja2.Environment(
            loader=jinja2.FileSystemLoader(config.TEMPLATE_DIR), autoescape=True
        ).get_template("template_report.html")

        html_placeholders = {
            "timestamp": self.timestamp,
            "app_version": git_tag(),
            "gene": self.gene,
            "sample_name": self.sample_name,
            "table": table_html,
            "logo": config.LOGOPATH,
            "parameters_str": config.PARAMS_STR,
            "footer_html": self.footer_html,
        } | config.INTERP_THRESHS

        with open(self.html_path, "w", encoding="utf-8") as html_file:
            html_file.write(template.render(html_placeholders))
        html_file.close()

    def create_pdffile(self):
        """
        Write metrics to pdf, specifying pdfkit options to turn off standard out and to allow
        pdfkit access to logo image
        """
        pdfkit.from_file(
            self.html_path,
            self.pdf_path,
            options={
                "enable-local-file-access": None,
                "quiet": "",
                "encoding": "UTF-8",
            },
        )


def arg_parse():
    """
    Parse arguments supplied by the command line. Create argument parser, define command line
    arguments, then parse supplied command line arguments using the created argument parser

        :return (Namespace object): parsed command line attributes
    """
    parser = argparse.ArgumentParser(
        description="Determine whether a known variant is present at "
        "the base preceding a poly stretch in a BAM file"
    )
    parser.add_argument(
        "-A",
        config.PARAMS["anchor_length"],
        type=int,
        nargs="?",
        default=2,
        const=2,
        help="Length of anchor sequence",
        required=False,
    )  # Optional arg
    requirednamed = parser.add_argument_group("Required named arguments")
    requirednamed.add_argument(
        "-B",
        config.PARAMS["bam"],
        type=validate_path,
        help="Bam file to analyse",
        required=True,
    )
    requirednamed.add_argument(
        "-I",
        config.PARAMS["bai"],
        type=validate_path,
        help="Bam index file",
        required=True,
    )
    requirednamed.add_argument(
        "-G", config.PARAMS["gene"], type=str, help="Gene of interest", required=True
    )
    requirednamed.add_argument(
        "-S",
        config.PARAMS["poly_start"],
        type=int,
        help="Start position of poly stretch",
        required=True,
    )
    requirednamed.add_argument(
        "-E",
        config.PARAMS["poly_end"],
        type=int,
        help="End position of poly stretch",
        required=True,
    )
    requirednamed.add_argument(
        "-C",
        config.PARAMS["chrom"],
        type=int,
        help="Chromosome of interest",
        required=True,
    )
    return vars(parser.parse_args())


def validate_path(path):
    """
    Check path exists

        :param path (str):  Input path to validate
    """
    if os.path.exists(path):
        return path
    else:
        raise argparse.ArgumentTypeError(f"{path} is not a valid path")


def git_tag():
    """Obtain git tag from current commit

    :return stdout (str):  String containing stdout, with newline characters removed
    """
    filepath = os.path.dirname(os.path.realpath(__file__))
    cmd = f"git -C {filepath} describe --tags"

    proc = subprocess.Popen(
        [cmd], stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True
    )
    out, _ = proc.communicate()
    return out.rstrip().decode("utf-8")


if __name__ == "__main__":
    args = arg_parse()
    bamfile = args["bam"]
    bam_index = args["bai"]
    gene = args["gene"]
    chrom = args["chrom"]
    poly = (args["poly_start"], args["poly_end"])
    anchor_length = args["anchor_length"]
    polyedge = PolyEdge(bamfile, bam_index, gene, chrom, poly, anchor_length)
    polyedge.generate_output()
