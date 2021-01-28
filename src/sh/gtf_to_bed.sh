#!/usr/bin/env sh

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                                                              #
# gtf_to_bed.sh                                                                #
#                                                                              #
# By: Maud Gautier <https://github.com/MaudGautier>, 2021                      #
#                                                                              #
# Broad permissions are granted to use, modify, and distribute this software   #
# as specified in the MIT License included in this distribution's LICENSE file.#
#                                                                              #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Creates BED file from specified lines of a GTF file.
# The start-stop limits can be extended with parameter `extension`.


## Usage
# sh ./src/sh/gtf_to_bed.sh \
	# -i <INPUT_GTF> \
	# -o <OUTPUT_BED> \
	# -e <EXTENSION_IN_BP> \
	# -u <UNIT:exon|transcript|gene>


## Requirements:
# - GNU awk
# - GNU sort
# - groupBy (bedtools)


## Input
# The input file must be of this form:
# chr10 scallop transcript	  19012211	19206483	1000  - . gene_id "Ew_NG4"; transcript_id "Ew_NG4.5"; RPKM "0.7306"; cov "11.2000";
# chr10 scallop exon	      19012211  19018797	1000  -	. gene_id "Ew_NG4"; transcript_id "Ew_NG4.5"; exon "1";
# chr10 scallop exon		  19130037  19130094	1000  -	. gene_id "Ew_NG4"; transcript_id "Ew_NG4.5"; exon "2";
# chr10 scallop exon		  19149907  19149966	1000  -	. gene_id "Ew_NG4"; transcript_id "Ew_NG4.5"; exon "3";
# chr10 scallop exon		  19206021  19206483	1000  -	. gene_id "Ew_NG4"; transcript_id "Ew_NG4.5"; exon "4";

# In particular, the 9-th column must have gene_id and transcript_id labels as 
# first and second labels (on 'transcript' lines), and gene_id, transcript_id 
# and exon as first, second and third labels (on 'exon' lines).


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
####                              PARAMETERS                               ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

while [[ $# -gt 1 ]]
do
	key="$1"

	case $key in
		-i)
			input_gtf="$2"
			shift
			;;
		-o)
			output_bed="$2"
			shift
			;;
		-u|--unit)
			unit="$2"
			shift
			;;
		-e|--extension)
			extension="$2"
			shift
			;;
		*)
			# unknown option
			;;
	esac
	shift # past argument or value
done

echo GENOMIC UNIT    = "${unit}"
echo INPUT GTF       = "${input_gtf}"
echo OUTPUT BED      = "${output_bed}"
echo EXTENSION       = "${extension}"



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
####                          DEFAULT PARAMETERS                           ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# If no extension => it is 0
if [ -v $extension ] ; then
	extension=0
fi


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
####                          CREATE BED FROM GTF                          ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

if [ $unit = "transcript" ] ; then

	# If transcripts
	awk -v FS="\t" -v OFS="\t" -v EXT=$extension '
	$3=="transcript" {
	split($9, a, ";") 
	split(a[2], b, " ") 
	gsub("\"","",b[2])
	print $1, $4-EXT, $5+EXT, b[2]}' ${input_gtf} \
		| sort -k1,1 -k2,2n \
		> ${output_bed}

elif [ $unit = "gene" ] ; then

	# If genes
	awk -v FS="\t" -v OFS="\t" -v EXT=$extension '
	$3=="transcript" {
	split($9, a, ";")
	split(a[1], b, " ")
	gsub("\"","",b[2])
	print $1, $4, $5, b[2], $7}' ${input_gtf} \
		| sort -k1,1 -k2,2n | uniq \
		| sort -k4,4 - \
		| groupBy -g 1,4 -c 2,3,5 -o min,max,first \
		| awk -v OFS="\t" -v EXT=$extension '{print $1, $3-EXT, $4+EXT, $2}' \
		| sort -k1,1V -k2,2n \
		> ${output_bed}

elif [ $unit = "exon" ] ; then
	awk -v FS="\t" -v OFS="\t" -v EXT=$extension '
	$3=="exon" {
	split($9, a, ";")
	split(a[2], b, " ")
	split(a[3], c, " ")
	gsub("\"","",b[2])
	gsub("\"","",c[2])
	print $1, $4-EXT, $5+EXT, b[2]"-"c[2]
	}' ${input_gtf} \
	| sort -k1,1 -k2,2n \
	> ${output_bed}

fi

