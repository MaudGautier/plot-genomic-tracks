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
# If the input file is in a BIGWIG format instead of BED or BAM, the option 
# `--bw|--bigwig` should be turned on. 


## Usage
# sh ./src/sh/create_per_base_coverage_table.sh \
	# -i <INPUT_BAM_OR_BED> \
	# --bed <BED_REGIONS> \
	# --gs <CHROM_SIZES> \
	# -o <OUTPUT_PREFIX> \
	# [--bigwig]


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
bigwig="false"

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
		--bigwig)
			bigwig="true"
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



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
####                          CREATE BED FROM GTF                          ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

## Add option bigwig




# Get per-base coverage
bedtools coverage -sorted -split -d \
	-g ${genome_sizes} \
	-a ${bed_regions} \
	-b ${input_file} \
	> ${output_prefix}.cov_per_base



## Add option bigwig



## Get final file
echo -e "#Chromosome\tPosition\tGene_ID\tCoverage" > ${output_prefix}.cov_per_base_final 
awk -v OFS="\t" '{
	print $1, $2 + $5 - 1, $4, $6
   }' ${output_prefix}.cov_per_base \
   >> ${output_prefix}.cov_per_base_final


## Delete intermediary file
rm -f ${output_prefix}.cov_per_base


