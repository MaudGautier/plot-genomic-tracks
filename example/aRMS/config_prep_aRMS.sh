#!/usr/bin/env sh

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                                                              #
# config_prep_aRMS.sh                                                          #
#                                                                              #
# By: Maud Gautier <https://github.com/MaudGautier>, 2021                      #
#                                                                              #
# Broad permissions are granted to use, modify, and distribute this software   #
# as specified in the MIT License included in this distribution's LICENSE file.#
#                                                                              #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Config plot for aRMS neotranscripts.



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                  Parameters                                  #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Generic parameters
output_folder=/data/kdi_prod/project_result/726/27.02/results/tests/
if [ ! -d $output_folder ] ; then mkdir $output_folder ; fi
input_folder=/data/kdi_prod/project_result/726/27.02/results/1_Data/
chrom_sizes_hg19="/data/annotations/pipelines/Human/hg19/genome/chrom_hg19.sizes"
chrom_sizes_sorted=${input_folder}/chrom_hg19.sizes.sorted
hg19_gencode=/data/annotations/pipelines/Human/hg19/gtf/gencode.v19.annotation.sorted.gtf
hg19_gencode_final=${output_folder}/gencode.v19.annotation.gtf

# Parameters for bed regions
gtf_name=200624_aRMS_neos_NEW_IDs
gtf_regions=${input_folder}/${gtf_name}.gtf
bed_regions=${output_folder}/genomic_regions.bed
extension=50000

# Parameters for sashimi
folder_junctions_plus=$output_folder/junctions_plus/
folder_junctions_minus=$output_folder/junctions_minus/
if [ ! -d $folder_junctions_plus ] ; then mkdir $folder_junctions_plus ; fi
if [ ! -d $folder_junctions_minus ] ; then mkdir $folder_junctions_minus ; fi
min_junctions=1

# List of tracks to prepare for plotting -- density plots
list_files_density=(
	ChIP-PF
	H3K27ac
	H3K4me3
	G284T01
	aRMS_merged
)

# List of tracks to prepare for plotting -- sashimi plots
list_sashimis=(aRMS_merged)




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                              Prepare track files                             #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# (Optional) Sort chrom sizes in correct order
if [ ! -f $chrom_sizes_sorted ] ; then
	grep "chr[123456789]" ${chrom_sizes_hg19} | grep -v random | grep -v hap \
		> ${chrom_sizes_sorted}
	grep "chrM" ${chrom_sizes_hg19} >> ${chrom_sizes_sorted}
	grep "chr[XY]" ${chrom_sizes_hg19} >> ${chrom_sizes_sorted}
fi


# Step 1: Prepare BED file of regions to cover for each neotranscript
if [ ! -f ${bed_regions} ] ; then
	echo "Creating ${bed_regions}..."
	sh ./src/sh/gtf_to_bed.sh \
		-i ${gtf_regions} \
		-o ${bed_regions} \
		-e ${extension} \
		-u "gene"
fi


# Step 2: Prepare per-base coverage files
for name in ${list_files_density[@]} ; do
	echo "Processing ${name} track..."
	sh ./src/sh/create_per_base_coverage_table.sh \
		-i ${input_folder}/${name}.bam \
		--bed ${bed_regions} \
		--gs ${chrom_sizes_sorted} \
		-o ${output_folder}/${name}.against.${gtf_name}
done


# Step 3: Prepare Rscript junction files for sashimi tracks
for name in ${list_sashimis[@]} ; do
	while read chrom start stop ID ; do
		echo $ID

		# Get position
		position=$chrom":"$start"-"$stop
		
		# Get junctions on + strand
		python3 ./src/sh/get_junctions_Rscript.py \
			-b ${input_folder}/${name}.bam \
			-c $position \
			-s MATE1_SENSE \
			-M ${min_junctions} \
			-r ${folder_junctions_plus}/${ID}.R

		# Get junctions on - strand
		python3 ./src/sh/get_junctions_Rscript.py \
			-b ${input_folder}/${name}.bam \
			-c $position \
			-s MATE2_SENSE \
			-M ${min_junctions} \
			-r ${folder_junctions_minus}/${ID}.R

	done < ${bed_regions}
done


# (Optional): Modify names in gene_id of GTF tracks for hg19 transcripts
if [ ! -f ${hg19_gencode_final} ] ; then
	grep "^##" ${hg19_gencode} > ${hg19_gencode_final}
	grep -v "^##" ${hg19_gencode} | awk -v FS="\t" -v OFS="\t" '
	{
		split($9,a,";")
		split(a[5],b," ")
		split(a[1], c, " ") 
		printf $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tgene_id "b[2] 
		for (i=2; i<=length(a); i++) {
			printf ";"a[i]
		}
		printf("\n") 
	}' >> ${hg19_gencode_final}
fi

