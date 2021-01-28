#!/usr/bin/env sh

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                                                              #
# config_prep_EwS.sh                                                           #
#                                                                              #
# By: Maud Gautier <https://github.com/MaudGautier>, 2021                      #
#                                                                              #
# Broad permissions are granted to use, modify, and distribute this software   #
# as specified in the MIT License included in this distribution's LICENSE file.#
#                                                                              #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Config plot for EwS neotranscripts.


# Generic parameters
output_folder=/data/kdi_prod/project_result/726/27.01/results/tests/
if [ ! -d $output_folder ] ; then mkdir $output_folder ; fi
input_folder=/data/kdi_prod/project_result/726/27.01/results/1_Data/
gtf_name=201201_Ewing_neos_NEW_IDs

# Step 1: Prepare BED file of regions to cover for each neotranscript
sh ./src/sh/gtf_to_bed.sh \
	-i ${input_folder}/${gtf_name}.gtf \
	-o ${output_folder}/genomic_regions.bed \
	-e 50000 \
	-u "gene"


# Sort chrom sizes in correct order
chrom_sizes_hg19="/data/annotations/pipelines/Human/hg19/genome/chrom_hg19.sizes"
chrom_sizes_sorted=${input_folder}/chrom_hg19.sizes.sorted
grep "chr[123456789]" ${chrom_sizes_hg19} | grep -v random | grep -v hap > ${chrom_sizes_sorted}
grep "chrM" ${chrom_sizes_hg19} >> ${chrom_sizes_sorted}
grep "chr[XY]" ${chrom_sizes_hg19} >> ${chrom_sizes_sorted}


# Step 2: Prepare per-base coverage files
name=EW_merged
sh ./src/sh/create_per_base_coverage_table.sh \
	-i ${input_folder}/${name}.bam \
	--bed ${output_folder}/genomic_regions.bed \
	--gs ${input_folder}/chrom_sizes_in_sorted_order.sizes \
	-o ${output_folder}/${name}.against.${gtf_name}



