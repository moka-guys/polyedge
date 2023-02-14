import pysam
import argparse
import subprocess
import re
import os
import statistics
from collections import defaultdict, Counter
from functools import reduce
import csv
import shutil
import jinja2
import datetime
# from xhtml2pdf import pisa
import pdfkit

class PolyEdge():
    """
    Class for processing bam files to determine whether a variant is present at the position before
    a poly stretch, and write to file

    generate_output()

    calculation()

    create_pdf()

    create_csv()

    """
    def __init__(self, bamfile, gene, chrom, poly, anchor_length):
        """
        :param bamfile (str):       Bam file to analyse
        :param gene (str) :         Gene of interest
        :param chrom (str):         Chromosome of interest
        :param position (int):      Position of interest
        :param anchor_length (int): Length of anchor sequence
        :return data (list):        List of tuples per allele, with each tuple containing the
                                    output stats as defined in the readme
        """
        self.root_dir = os.path.dirname(os.path.abspath(__file__)) # This is your Project Root
        self.bamfile = bamfile
        self.gene = gene
        self.chrom = chrom
        self.chrCHROM = f"chr{self.chrom}"
        self.anchor_length = anchor_length
        self.poly = poly
        # Define start and end of segment to be matched
        self.roi_start = self.poly[0]-self.anchor_length
        self.roi_end = self.poly[1]+self.anchor_length
        # Variables required for creating output files
        self.bamfile_prefix = os.path.splitext(bamfile)[0].split("/")[-1]
        self.output_string = f"{self.bamfile_prefix}.{self.gene}_polyedge"
        self.template_dir = os.path.join(self.root_dir, 'templates')
        self.timestamp = datetime.datetime.now().strftime('%d-%B-%Y %H:%M')

        self.csv_template = os.path.join(self.template_dir, 'template_report.csv')
        self.csv_path = os.path.join(os.getcwd(), f"{self.output_string}.csv")
        self.csv_data_format = '{:>},{:4n},{:2n},{:.2f},{:2.2f},{:.2f},{:2n},{:.2f}'

        self.html_template = os.path.join(self.template_dir, 'template_report.html')
        self.html_path = os.path.join(os.getcwd(), f"{self.output_string}.html")
        self.html_data_format = "<tr><td>{:>}</td><td>{:4n}</td><td>{:2n}</td><td>{:.2f}</td>"\
                                "<td>{:2.2f}</td><td>{:.2f}</td><td>{:2n}</td><td>{:.2f}</td></tr>"
        self.pdf_path = os.path.join(os.getcwd(), f"{self.output_string}.pdf")


    def generate_output(self):
        """
        Call methods required to generate all output files for the sample provided
        """
        allele_list = self.calculate_metrics()
        print(self.bamfile_prefix)
        # self.create_csv(allele_list)
        self.create_html(allele_list)

    def calculate_metrics(self):
        """
        Calculate metrics: first base of poly repeat allele, read count, mean quality of first base,
        fraction of reads, mean poly length, standard deviation of poly length, mode of poly length,
        average purity of polyN repeat
        """
        # Poly length counts by first base allele
        counts = defaultdict(list)  # Lengths of repeat given allele (one per read)
        quals = defaultdict(list)  # Quality of first base
        purity = defaultdict(list)  # Most common base count in poly w/o first base
        reader = pysam.Samfile(self.bamfile, 'rb')  # Open BAM file
        # Fetch reads intersecting the poly repeat
        try:
            roi_reads = reader.fetch(str(self.chrom), *self.poly)
        except ValueError as e:
            print(f"Using CHROM {str(self.chrCHROM)} raised an exception ({e}). "
                  f"Trying: CHROM '{self.chrCHROM}'")
            roi_reads = reader.fetch(self.chrCHROM, *self.poly)
        for read in roi_reads:
            # Select that span poly with defined ANCHOR sequence length
            if read.reference_start < self.roi_start and self.roi_end < read.reference_end:
                # Find sequence segment with poly and anchor sequence
                # -> first matching position in the query->reference mapping
                seq_start = next(x for x in read.get_aligned_pairs()
                                 if x[1] == self.roi_start)
                seq_end = next(x for x in read.get_aligned_pairs()
                               if x[1] == self.roi_end)
                if seq_start[0] and seq_end[0]:  # If anchored
                    seq_segment = read.query_sequence[seq_start[0]:seq_end[0]]
                    first_base_pos = seq_start[0]+self.anchor_length
                    first_base_qual = read.query_qualities[first_base_pos]
                    # Extract poly sequence (excluding anchor bases)
                    # Assuming no indel in anchor sequence
                    regex = r'.{' + re.escape(str(self.anchor_length)) + \
                        r'}(.+).{' + re.escape(str(self.anchor_length)) + r'}'
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
            data = (allele,
                    len(counts[allele]),
                    int(statistics.mean(quals[allele])),
                    len(counts[allele])/total_count,
                    statistics.mean(counts[allele]),
                    statistics.stdev(counts[allele]) if len(counts[allele]) > 1 else 0,
                    statistics.mode(counts[allele]),
                    sum(purity[allele])/(sum(counts[allele])-len(counts[allele]))
                    )
            allele_list.append(data)
        return allele_list

    def create_csv(self, allele_list):
        """
        Write metrics to csv. Copy template csv to output dest, write data to end of file
        """
        shutil.copyfile(self.csv_template, self.csv_path)
        with open(self.csv_path, 'a', newline='', encoding='utf-8') as file:
            writer = csv.writer(file, delimiter=',')
            for row in allele_list:
                line_to_write = self.csv_data_format.format(*row).split(",")
                writer.writerow(line_to_write)
        file.close()

    def create_html(self, allele_list):
        """
        Write metrics to html
        """
        table_lines = ""
        for row in allele_list:
            table_lines += self.html_data_format.format(*row)

        logopath = os.path.join(self.root_dir, "images/logo.png")

        jinja_env = jinja2.Environment(loader=jinja2.FileSystemLoader(self.template_dir),
                                       autoescape=True)
        template = jinja_env.get_template('template_report.html')
        place_holder_values = {"timestamp": self.timestamp,
                               "app_version": git_tag(),
                               "gene": self.gene,
                               "sample_name": self.bamfile_prefix.split(".")[0],
                               "results": table_lines,
                               "logo": logopath}

        with open(self.html_path, "w", encoding="utf-8") as html_file:
            html_file.write(template.render(place_holder_values))

        # Specify pdfkit options to turn off standard out and also allow access to the images
        pdfkit.from_file(self.html_path, self.pdf_path,
                         options={'enable-local-file-access': None, "quiet": '', 'encoding': "UTF-8"})

def arg_parse():
    """
    Parses arguments supplied by the command line. Creates argument parser, defines command line
    arguments, then parses supplied command line arguments using the created argument parser.

        :return (Namespace object): parsed command line attributes
    """
    parser = argparse.ArgumentParser(description='Determine whether a known variant is present at '
                                                 'the base preceding a poly stretch in a BAM file')
    requiredNamed = parser.add_argument_group('Required named arguments')
    requiredNamed.add_argument('-B', '--bam', type=dir_path, help='Bam file to analyse',
                               required=True)
    requiredNamed.add_argument('-I', '--bai', type=dir_path, help='Bam index file', required=True)
    requiredNamed.add_argument('-G', '--gene', type=str, help='Gene of interest', required=True)
    requiredNamed.add_argument('-S', '--poly_start', type=int,
                               help='Start position of poly stretch', required=True)
    requiredNamed.add_argument('-E', '--poly_end', type=int,
                               help='End position of poly stretch', required=True)
    requiredNamed.add_argument('-A', '--anchor_length', type=int,
                               help='Length of anchor sequence', required=True)
    requiredNamed.add_argument('-C', '--chrom', type=int,
                               help="Chromosome of interest", required=True)
    return vars(parser.parse_args())


def dir_path(path):
    """
    Checks the command line argument provided is an existing directory
    """
    if os.path.exists(path):
        return path
    else:
        raise argparse.ArgumentTypeError(f"readable_dir:{path} is not a valid path")

def git_tag():
    """Obtain the git tag of the current commit"""
    filepath = os.path.dirname(os.path.realpath(__file__))
    cmd = f"git -C {filepath} describe --tags"

    proc = subprocess.Popen(
        [cmd], stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True
    )
    out, _ = proc.communicate()
    #  Return standard out, removing any new line characters
    return out.rstrip().decode("utf-8")


if __name__ == "__main__":
    git_tag()
    print(git_tag())



if __name__ == "__main__":
    args = arg_parse()
    bamfile = args['bam']
    # Whilst not explicitly specified within the script
    # this is required to fetch the reads (reader.fetch())
    bam_index = args['bai']
    gene = args['gene']
    chrom = args['chrom']
    poly = (args['poly_start'], args['poly_end'])
    anchor_length = args['anchor_length']

    polyedge = PolyEdge(bamfile, gene, chrom, poly, anchor_length)
    polyedge.generate_output()

    # GENE = "MSH2"
    # CHROM = "2"
    # POS = 47641559
    # POLY = (POS, 47641586)
    # ANCHOR_LENGTH = 2
