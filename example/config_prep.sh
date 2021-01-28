#!/usr/bin/env sh

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                                                              #
# config_prep.sh                                                               #
#                                                                              #
# By: Maud Gautier <https://github.com/MaudGautier>, 2021                      #
#                                                                              #
# Broad permissions are granted to use, modify, and distribute this software   #
# as specified in the MIT License included in this distribution's LICENSE file.#
#                                                                              #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Example config plot



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                  Parameters                                  #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Generic parameters
output_folder=./test/
if [ ! -d $output_folder ] ; then mkdir $output_folder ; fi
input_folder=./test/input/
chrom_sizes_sorted=${input_folder}/chrom_hg19.sizes.sorted

# Parameters for bed regions
gtf_name=201201_Ewing_neos_NEW_IDs
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
	EW_merged 
	ASP14-d7-FLI1 
)

# List of tracks to prepare for plotting -- GGAA plots
list_repeats=(hg19_GGAA_TTCC_20151221)

# List of tracks to prepare for plotting -- sashimi plots
list_sashimis=(EW_merged)




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
##  2a -- Density tracks
for name in ${list_files_density[@]} ; do
	echo "Processing ${name} track..."
	sh ./src/sh/create_per_base_coverage_table.sh \
		-i ${input_folder}/${name}.bam \
		--bed ${bed_regions} \
		--gs ${chrom_sizes_sorted} \
		-o ${output_folder}/${name}.against.${gtf_name}
done

##  2b -- GGAA track (from BIGWIG file)  -  *requires `--bigwig` option*
for name in ${list_repeats[@]} ; do
	sh ./src/sh/create_per_base_coverage_table.sh \
		--repeat \
		-i ${input_folder}/${name}.bw \
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

