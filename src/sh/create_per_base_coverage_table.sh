#!/usr/bin/env sh

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                                                              #
# create_per_base_coverage_table.sh                                            #
#                                                                              #
# By: Maud Gautier <https://github.com/MaudGautier>, 2021                      #
#                                                                              #
# Broad permissions are granted to use, modify, and distribute this software   #
# as specified in the MIT License included in this distribution's LICENSE file.#
#                                                                              #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Creates per-base coverage file (`-o`) on regions defined by input BED file 
# (`--bed|-b`) from input file that can be in either a BED or a BAM format
# (`-i`). The process is sped up if a chrom_sizes file is provided (`--gs|-g`).
# If the input file is in a BIGWIG format instead of BED or BAM and is to be
# processed as a GGAA track (i.e. make repeats visible), the option 
# `--rep|--repeat` should be turned on. 


## Usage
# sh ./src/sh/create_per_base_coverage_table.sh \
	# [--bigwig] \
	# -i <INPUT_BAM_OR_BED> \
	# --bed <BED_REGIONS> \
	# --gs <CHROM_SIZES> \
	# -o <OUTPUT_PREFIX>


## Requirements:
# - GNU awk
# - GNU sort
# - bedtools
# - bigWigToBedGraph


## Input
# The input chrom_sizes should be in this form:
# chr1  249250621
# chr2  243199373
# chr3  198022430
# chr4  191154276
# chr5  180915260


## Output
# The output file will contain 4 columns defining the depth (vertical coverage)
# at each position of the selected bed regions. It will be in this format:
# #Chromosome  Position  Gene_ID  Coverage
# chr1         73171083  Ew_NG1   7
# chr1         73171084  Ew_NG1   7
# chr1         73171085  Ew_NG1   7
# chr1         73171086  Ew_NG1   8


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
####                              PARAMETERS                               ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Default:
repeat="false"

while [[ $# -gt 1 ]]
do
	key="$1"

	case $key in
		-i)
			input_file="$2"
			shift
			;;
		-o)
			output_prefix="$2"
			shift
			;;
		-b|--bed)
			bed_regions="$2"
			shift
			;;
		-g|--gs)
			genome_sizes="$2"
			shift
			;;
		--rep|--repeat)
			repeat="true"
			;;
		*)
			# unknown option
			;;
	esac
	shift # past argument or value
done

echo GENOME SIZES    = "${genome_sizes}"
echo INPUT - BAM/BED = "${input_file}"
echo OUTPUT PREFIX   = "${output_prefix}"
echo BED REGIONS     = "${bed_regions}"
echo INPUT IS REPEAT = "${repeat}"



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
####                          CREATE BED FROM GTF                          ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# (Optional) If the input file is a bigwig file, tranform it to BED format
if [ $repeat = "true" ] ; then
	# Transform to bedgraph
	bigWigToBedGraph ${input_file} ${input_file}.bedgraph

	# Transform bedgraph to bed
	awk '{ \
			if ($1 ~ /^chr/) { \
		   print $1"\t"$2"\t"$3"\tid-"NR"\t"$4; \
			} \
		}' ${input_file}.bedgraph \
	> ${input_file}.bed

	# Sort
	sort -k1,1V -k2,2n ${input_file}.bed > ${input_file}.sorted.bed

	# Remove intermediary files
	rm -f ${input_file}.bedgraph
	rm -f ${input_file}.bed

	# Rename input
	temp_file=${input_file}.sorted.bed
	input_file=${input_file}.sorted.bed
fi


# Get per-base coverage
bedtools coverage -sorted -split -d \
	-g ${genome_sizes} \
	-a ${bed_regions} \
	-b ${input_file} \
	> ${output_prefix}.cov_per_base


# Get final file
echo -e "#Chromosome\tPosition\tGene_ID\tCoverage" \
	> ${output_prefix}.cov_per_base_final

awk -v OFS="\t" '{
	print $1, $2 + $5 - 1, $4, $6
   }' ${output_prefix}.cov_per_base \
   >> ${output_prefix}.cov_per_base_final


# (Optional) Get height proportional to number of repeats
if [ $repeat = "true" ] ; then

	mv ${output_prefix}.cov_per_base_final \
		${output_prefix}.cov_per_base_no_height

	awk -v OFS="\t" '
	NR==1 { print }
	NR > 1 {
		if ($4==1) { c++ }

		if ($4 == 1 && prev_cov == 0) {
			chrom = $1
			start = $2
			stop = $2
			ID = $3
			c = 1
		}

		if ($4 == 0 && prev_cov == 1 && ID == $3) {
			stop = prev_pos
			tot_cov=c/4
			for (i=start; i<=stop; i++) {
				print prev_chrom,i,prev_ID,tot_cov
			}
		}

		if ($4 == 0) { print }

		prev_chrom=$1
		prev_pos=$2
		prev_ID=$3
		prev_cov=$4
	}' ${output_prefix}.cov_per_base_no_height \
	> ${output_prefix}.cov_per_base_final

fi


## Delete intermediary file
rm -f ${output_prefix}.cov_per_base
if [ $repeat = "true" ] ; then
	rm -f ${temp_file}
	rm -f ${output_prefix}.cov_per_base_no_height
fi

