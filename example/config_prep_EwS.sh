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



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                  Parameters                                  #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Generic parameters
output_folder=/data/kdi_prod/project_result/726/27.01/results/tests/
if [ ! -d $output_folder ] ; then mkdir $output_folder ; fi
input_folder=/data/kdi_prod/project_result/726/27.01/results/1_Data/
gtf_name=201201_Ewing_neos_NEW_IDs
chrom_sizes_hg19="/data/annotations/pipelines/Human/hg19/genome/chrom_hg19.sizes"
chrom_sizes_sorted=${input_folder}/chrom_hg19.sizes.sorted # Prepared below

# List of tracks to prepare for plotting -- density plots
list_files_merged=(
	EW_merged 
	FLI1_merged 
	H3K27ac_merged 
	H3K4me3_merged
)
list_files_d0_d7=(
	ASP14-d0-FLI1    ASP14-d7-FLI1 
	ASP14-d0-H3K27ac ASP14-d7-H3K27ac
	ASP14-d0-H3K4me3 ASP14-d7-H3K4me3
)
list_files_tum_CL=(
	B69T10
	B69T13
	G402T05
	G315T01
	G312T05
)





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                              Prepare track files                             #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# (Optional) Sort chrom sizes in correct order
if [ ! -f $chrom_sizes_sorted ] ; then
	grep "chr[123456789]" ${chrom_sizes_hg19} | grep -v random | grep -v hap > ${chrom_sizes_sorted}
	grep "chrM" ${chrom_sizes_hg19} >> ${chrom_sizes_sorted}
	grep "chr[XY]" ${chrom_sizes_hg19} >> ${chrom_sizes_sorted}
fi


# Step 1: Prepare BED file of regions to cover for each neotranscript
sh ./src/sh/gtf_to_bed.sh \
	-i ${input_folder}/${gtf_name}.gtf \
	-o ${output_folder}/genomic_regions.bed \
	-e 50000 \
	-u "gene"


# Step 2: Prepare per-base coverage files
##  2a -- Density tracks
for name in ${list_files_merged[@]} ${list_files_d0_d7[@]} ${list_files_tum_CL[@]} ; do
	echo "Processing ${name}..."
	sh ./src/sh/create_per_base_coverage_table.sh \
		-i ${input_folder}/${name}.bam \
		--bed ${output_folder}/genomic_regions.bed \
		--gs ${input_folder}/chrom_sizes_in_sorted_order.sizes \
		-o ${output_folder}/${name}.against.${gtf_name}
done

##  2b -- GGAA track (from BIGWIG file)  -  *requires `--bigwig` option*
sh ./src/sh/create_per_base_coverage_table.sh \
	--bigwig \
	-i ${input_folder}/hg19_GGAA_TTCC_20151221.bw \
	--bed ${output_folder}/genomic_regions.bed \
	--gs ${input_folder}/chrom_sizes_in_sorted_order.sizes \
	-o ${output_folder}/hg19_GGAA_TTCC_20151221.against.${gtf_name}


# Step 3: Prepare Rscript junction files for sashimi tracks
tsv_list_of_bams_MERGED=$SASHIMI/list_bams_MERGED.tsv
/data/kdi_prod/project_result/726/27.01/results/4_Sashimi/list_bams_MERGED.tsv
# echo -e "Pile-up\t$DATA/EW_merged.bam" > $tsv_list_of_bams_MERGED
while read chrom start stop ID ; do
	echo $ID

	# Get position
	position=$chrom":"$start"-"$stop
	
	# Get junctions on + strand
	python3 src/sh/get_junctions_Rscript.py \
		-b ${input_folder}/EW_merged.bam \
		-c $position \
		-s MATE1_SENSE \
		-M 1 \
		-r $output_folder/junctions_plus/${ID}.R

	# Get junctions on - strand
	python3 src/sh/get_junctions_Rscript.py \
		-b ${input_folder}/EW_merged.bam \
		-c $position \
		-s MATE2_SENSE \
		-M 1 \
		-r $output_folder/junctions_minus/${ID}.R

done < ${output_folder}/genomic_regions.bed



# Optional: Modify names in gene_id of GTF tracks

