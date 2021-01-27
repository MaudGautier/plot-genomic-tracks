#!/usr/bin/env R

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                                                              #
# config_plot_cell_lines.R                                                     #
#                                                                              #
# By: Maud Gautier <https://github.com/MaudGautier>, 2021                      #
#                                                                              #
# Broad permissions are granted to use, modify, and distribute this software   #
# as specified in the MIT License included in this distribution's LICENSE file.#
#                                                                              #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# This is the configuration file used to plot Ewing cell lines



# Usage -------------------------------------------------------------------

# From the CLI:
# Rscript src/R/make_genomic_plot.R example/EwS/config_plot_cell_lines.R




# Configuration for EwS neotranscripts ------------------------------------


## General parameters
output_folder <- "/Users/maudgautier/Documents/github-savings/plot-genomic-tracks/plots/"
main_folder <- "/Users/maudgautier/Documents/github-savings/plot-genomic-tracks/"
gtf_genes_path <- '/Users/maudgautier/Documents/data/tmp_from_calcsub/figs_EwS_100kb/1_Data/201201_Ewing_neos_NEW_IDs.gtf'
tracks_folder <- "/Users/maudgautier/Documents/data/tmp_from_calcsub/figs_EwS_100kb/3d_IGV_New_IDs/Against_genes_of_201201_Ewing_neos_NEW_IDs/"

## Tracks
tracks_info <- list(
  "GGAA" = list(path = paste0(tracks_folder, "hg19_GGAA_TTCC_20151221.against.201201_Ewing_neos_NEW_IDs.cov_per_base_final"), 
                group = "GGAA", 
                color = "black", 
                height = 0.2, 
                type = "GGAA"),
  "EF_low" = list(path = paste0(tracks_folder, "ASP14-d7-FLI1.against.201201_Ewing_neos_NEW_IDs.cov_per_base_final"), 
                   group = "FLI1", 
                   color = "red", 
                   height = 0.3, 
                   type = "density"),
  "EF_high" = list(path = paste0(tracks_folder, "FLI1_merged.against.201201_Ewing_neos_NEW_IDs.cov_per_base_final"), 
                   group = "FLI1",
                   color = "red", 
                   height = 0.3, 
                   type = "density"),
  "H3K27ac" = list(path = paste0(tracks_folder, "H3K27ac_merged.against.201201_Ewing_neos_NEW_IDs.cov_per_base_final"), 
                   group = "H3K27", 
                   color = "green", 
                   height = 0.3, 
                   type = "density"),
  "H3K4me3" = list(path = paste0(tracks_folder, "H3K4me3_merged.against.201201_Ewing_neos_NEW_IDs.cov_per_base_final"), 
                   group = "H3K4", 
                   color = "blue", 
                   height = 0.3, 
                   type = "density"),
  "Cell line" = list(path = paste0(tracks_folder, "B69T10.against.201201_Ewing_neos_NEW_IDs.cov_per_base_final"), 
                   group = "cl", 
                   color = "grey", 
                   height = 0.3, 
                   type = "density"),
  "Prediction" = list(path = "/Users/maudgautier/Documents/data/tmp_from_calcsub/figs_EwS_100kb/1_Data/201201_Ewing_neos_NEW_IDs.gtf", 
                      group = "neotranscripts", 
                      height = 0.4, 
                      type = "transcripts"),
  "All genes" = list(path = "/Users/maudgautier/Documents/data/tmp_from_calcsub/figs_EwS/gencode.v19.annotation.sorted.gtf_MODIF", 
                      group = "neotranscripts", 
                      height = 0.4, 
                      type = "transcripts"),
  "Pile-up" = list(path = "/Users/maudgautier/Documents/data/tmp_from_calcsub/figs_EwS_100kb/3d_IGV_New_IDs/Against_genes_of_201201_Ewing_neos_NEW_IDs/EW_MERGED.against.201201_Ewing_neos_NEW_IDs.cov_per_base_final", 
                 group = "pileup", 
                 sashimi_plus = "/Users/maudgautier/Documents/data/tmp_from_calcsub/figs_EwS_100kb/4_junctions/junctions_100kb_Rscript_PLUS/", 
                 sashimi_minus = "/Users/maudgautier/Documents/data/tmp_from_calcsub/figs_EwS_100kb/4_junctions/junctions_100kb_Rscript_MINUS/", 
                 color = "grey", 
                 height = 0.6, 
                 type = "sashimi", 
                 default_min_junctions = 200,
                 adapted_min_junctions = list(
                   "Ew_NG6" = 100,
                   "Ew_NG8" = 500,
                   "Ew_NG11" = 400,
                   "Ew_NG20" = 3000
                 )
  )
)



## Genomic regions
genes <- sprintf("Ew_NG%s",c(13, 16, 20))
default_half_width <- 10000
half_widths <- list("Ew_NG20" = 50000)


