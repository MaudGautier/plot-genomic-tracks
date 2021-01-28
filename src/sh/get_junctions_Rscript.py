#!/usr/bin/env python

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                                                              #
# get_junctions_Rscript.py                                                     #
#                                                                              #
# By: Maud Gautier <https://github.com/MaudGautier>, 2021                      #
# Original code by https://github.com/guigolab/ggsashimi                       #
# Subsetted and modified by Maud Gautier                                       #
# Credit for the original code to guigolab                                     #
#                                                                              #
# Broad permissions are granted to use, modify, and distribute this software   #
# as specified in the MIT License included in this distribution's LICENSE file.#
#                                                                              #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Creates an R script containing a dataframe listing all junctions that are
# supported by more than a given number of reads (`-M` argument).
# The dataframe contains the coordinates of the arcs that will be used to draw
# the junctions and the number of reads supporting the junction.
# The dataframe has the following pieces of information:
# - x, y: coordinates of one side of the junction arc
# - xend, yend: coordinates of the other side of the junction arc
# - count: number of reads supporting the junction


## Usage
#  python3 src/sh/get_junctions_Rscript.py \
	#  -b <BAM_FILE> \
	#  -c <CHROM:START-STOP> \
	#  -s <STRAND:MATE1_SENSE|MATE2_SENSE> \
	#  -M <MIN_JUNCTIONS> \
	#  -r <OUTPUT_R_SCRIPT>


# Import modules
from argparse import ArgumentParser
import subprocess as sp
import sys, re, copy, os, codecs
from collections import OrderedDict

def define_options():
	# Argument parsing
	parser = ArgumentParser(description='Create sashimi plot for a given genomic region')
	parser.add_argument("-b", "--bam", type=str, required=True,
			help="""
			Individual bam file or file with a list of bam files.
			In the case of a list of files the format is tsv:
			1col: id for bam file,
			2col: path of bam file,
			3+col: additional columns
			""")
	parser.add_argument("-c", "--coordinates", type=str, required=True,
			help="Genomic region. Format: chr:start-end. Remember that bam coordinates are 0-based")
	parser.add_argument("-r", "--rfile", type=str, dest="rfile", default="R_script",
			help="Name for R script [default=%(default)s]")
	parser.add_argument("-M", "--min-coverage", type=int, default=1, dest="min_coverage",
			help="Minimum number of reads supporting a junction to be drawn [default=1]")
	parser.add_argument("-s", "--strand", default="NONE", type=str,
			help="Strand specificity: <NONE> <SENSE> <ANTISENSE> <MATE1_SENSE> <MATE2_SENSE> [default=%(default)s]")
	parser.add_argument("-L", "--labels", type=int, dest="labels", default=1,
			help="Index of column with labels (1-based) [default=%(default)s]")
	return parser



def parse_coordinates(c):
	c = c.replace(",", "")
	chr = c.split(":")[0]
	start, end = c.split(":")[1].split("-")
	# Convert to 0-based
	start, end = int(start) - 1, int(end)
	return chr, start, end



def count_operator(CIGAR_op, CIGAR_len, pos, start, end, a, junctions):

	# Match
	if CIGAR_op == "M":
		for i in range(pos, pos + CIGAR_len):
			if i < start or i >= end:
				continue
			ind = i - start
			a[ind] += 1

	# Insertion or Soft-clip
	if CIGAR_op == "I" or CIGAR_op == "S":
		return pos

	# Deletion
	if CIGAR_op == "D":
		pass

	# Junction
	if CIGAR_op == "N":
		don = pos
		acc = pos + CIGAR_len
		if don > start and acc < end:
			junctions[(don,acc)] = junctions.setdefault((don,acc), 0) + 1

	pos = pos + CIGAR_len

	return pos


def flip_read(s, samflag):
	if s == "NONE" or s == "SENSE":
		return 0
	if s == "ANTISENSE":
		return 1
	if s == "MATE1_SENSE":
		if int(samflag) & 64:
			return 0
		if int(samflag) & 128:
			return 1
	if s == "MATE2_SENSE":
		if int(samflag) & 64:
			return 1
		if int(samflag) & 128:
			return 0


def read_bam(f, c, s):

	_, start, end = parse_coordinates(c)

	# Initialize coverage array and junction dict
	a = {"+" : [0] * (end - start)}
	junctions = {"+": OrderedDict()}
	if s != "NONE":
		a["-"] = [0] * (end - start)
		junctions["-"] = OrderedDict()

	p = sp.Popen("samtools view %s %s " %(f, c), shell=True, stdout=sp.PIPE)
	for line in p.communicate()[0].decode('utf8').strip().split("\n"):

		if line == "":
			continue

		line_sp = line.strip().split("\t")
		samflag, read_start, CIGAR = line_sp[1], int(line_sp[3]), line_sp[5]

		# Ignore reads with more exotic CIGAR operators
		if any(map(lambda x: x in CIGAR, ["H", "P", "X", "="])):
			continue

		read_strand = ["+", "-"][flip_read(s, samflag) ^ bool(int(samflag) & 16)]
		if s == "NONE": read_strand = "+"

		CIGAR_lens = re.split("[MIDNS]", CIGAR)[:-1]
		CIGAR_ops = re.split("[0-9]+", CIGAR)[1:]

		pos = read_start

		for n, CIGAR_op in enumerate(CIGAR_ops):
			CIGAR_len = int(CIGAR_lens[n])
			pos = count_operator(CIGAR_op, CIGAR_len, pos, start, end, a[read_strand], junctions[read_strand])

	p.stdout.close()
	return a, junctions

def get_bam_path(index, path):
	if os.path.isabs(path):
		return path
	base_dir = os.path.dirname(index)
	return os.path.join(base_dir, path)

def read_bam_input(f, label):
	if f.endswith(".bam"):
		bn = f.strip().split("/")[-1].strip(".bam")
		yield bn, f, bn
		return
	with codecs.open(f, encoding='utf-8') as openf:
		for line in openf:
			line_sp = line.strip().split("\t")
			bam = get_bam_path(f, line_sp[1])
			label_text = line_sp[label-1] if label else None
			yield line_sp[0], bam, label_text


def prepare_for_R(a, junctions, c, m):

	_, start, _ = parse_coordinates(args.coordinates)

	# Convert the array index to genomic coordinates
	x = list(i+start for i in range(len(a)))
	y = a

	# Arrays for R
	dons, accs, yd, ya, counts = [], [], [], [], []

	# Prepare arrays for junctions (which will be the arcs)
	for (don, acc), n in junctions.items():

		# Do not add junctions with less than defined coverage
		if n < m:
			continue

		dons.append(don)
		accs.append(acc)
		counts.append(n)

		yd.append( a[ don - start -1 ])
		ya.append( a[ acc - start +1 ])

	return x, y, dons, accs, yd, ya, counts


def make_R_lists(id_list, d):
	s = ""
	# Iterate over ids to get bam signal and junctions
	for k in id_list:
		x, y, dons, accs, yd, ya, counts = d[k]
		s += """data.frame(x=c(%(dons)s), xend=c(%(accs)s), y=c(%(yd)s), yend=c(%(ya)s), count=c(%(counts)s))

""" %({
	'x' : ",".join(map(str, x)),
	'y' : ",".join(map(str, y)),
	'dons' : ",".join(map(str, dons)),
	'accs' : ",".join(map(str, accs)),
	'yd' : ",".join(map(str, yd)),
	'ya' : ",".join(map(str, ya)),
	'counts' : ",".join(map(str, counts))
	})
	return s


if __name__ == "__main__":

	strand_dict = {"plus": "+", "minus": "-"}

	parser = define_options()
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = parser.parse_args()

	bam_dict, id_list, label_dict = {"+":OrderedDict()}, [], OrderedDict()
	if args.strand != "NONE": bam_dict["-"] = OrderedDict()

	for id, bam, label_text in read_bam_input(args.bam, args.labels):
		if not os.path.isfile(bam):
			continue
		a, junctions = read_bam(bam, args.coordinates, args.strand)
		if a.keys() == ["+"] and all(map(lambda x: x==0, list(a.values()[0]))):
			print("WARN: Sample {} has no reads in the specified area.".format(id))
			continue
		id_list.append(id)
		label_dict[id] = label_text
		for strand in a:
			bam_dict[strand][id] = prepare_for_R(a[strand], junctions[strand], args.coordinates, args.min_coverage)
	# No bam files
	if not bam_dict["+"]:
		print("ERROR: No available bam files.")
		exit(1)

	# Iterate for plus and minus strand
	for strand in bam_dict:
		R_script = make_R_lists(id_list, bam_dict[strand])
		with open(args.rfile, 'w') as r:
			r.write(R_script)
	exit()


