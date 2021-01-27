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

# Step 1: Prepare BED file of regions to cover for each neotranscript
sh ./src/sh/gtf_to_bed.sh \
	-i $input_folder/201201_Ewing_neos_NEW_IDs.gtf \
	-o $output_folder/genomic_regions.bed \
	-e 50000 \
	-u "gene"

# Step 2: Prepare per-base coverage files


