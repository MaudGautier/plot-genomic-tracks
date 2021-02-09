#!/usr/bin/env R

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                                                              #
# config_plot_PACBIO.R                                                         #
#                                                                              #
# By: Maud Gautier <https://github.com/MaudGautier>, 2021                      #
#                                                                              #
# Broad permissions are granted to use, modify, and distribute this software   #
# as specified in the MIT License included in this distribution's LICENSE file.#
#                                                                              #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# This is the configuration file used to plot Ewing tumors.



# Usage -------------------------------------------------------------------

# From the CLI:
# Rscript src/R/make_genomic_plot.R example/EwS/config_plot_PACBIO.R




# Configuration for EwS neotranscripts ------------------------------------


## General parameters
output_folder <- "/Users/maudgautier/Documents/github-savings/plot-genomic-tracks/plots/EwS_PACBIO/"
main_folder <- "/Users/maudgautier/Documents/github-savings/plot-genomic-tracks/"
gtf_genes_path <- '/Users/maudgautier/Documents/github-savings/plot-genomic-tracks/input/EwS/PACBIO.gtf'
tracks_folder <- "/Users/maudgautier/Documents/github-savings/plot-genomic-tracks/input/EwS/"
label_size <- 1.3

## Tracks
tracks_info <- list(
  "GGAA" = list(path = paste0(tracks_folder, "hg19_GGAA_TTCC_20151221.against.201201_Ewing_neos_NEW_IDs.cov_per_base_final"), 
                group = "GGAA", 
                color = "black", 
                height = 0.2, 
                type = "GGAA"),
  "EF d0" = list(path = paste0(tracks_folder, "ASP14-d0-FLI1.against.201201_Ewing_neos_NEW_IDs.cov_per_base_final"), 
                 group = "FLI1",
                 color = "red", 
                 height = 0.3, 
                 type = "density"),
  "EF d7" = list(path = paste0(tracks_folder, "ASP14-d7-FLI1.against.201201_Ewing_neos_NEW_IDs.cov_per_base_final"), 
                  group = "FLI1", 
                  color = "red", 
                  height = 0.3, 
                  type = "density"),
  "H3K27ac\nd0" = list(path = paste0(tracks_folder, "ASP14-d0-H3K27ac.against.201201_Ewing_neos_NEW_IDs.cov_per_base_final"), 
                      group = "H3K27", 
                      color = "green", 
                      height = 0.3, 
                      type = "density"),
  "H3K27ac\nd7" = list(path = paste0(tracks_folder, "ASP14-d7-H3K27ac.against.201201_Ewing_neos_NEW_IDs.cov_per_base_final"), 
                   group = "H3K27", 
                   color = "green", 
                   height = 0.3, 
                   type = "density"),
  "H3K4me3\nd0" = list(path = paste0(tracks_folder, "ASP14-d0-H3K4me3.against.201201_Ewing_neos_NEW_IDs.cov_per_base_final"), 
                      group = "H3K4", 
                      color = "blue", 
                      height = 0.3, 
                      type = "density"),
  "H3K4me3\nd7" = list(path = paste0(tracks_folder, "ASP14-d7-H3K4me3.against.201201_Ewing_neos_NEW_IDs.cov_per_base_final"), 
                   group = "H3K4", 
                   color = "blue", 
                   height = 0.3, 
                   type = "density"),
  "Prediction" = list(path = paste0(tracks_folder,"PACBIO.gtf"), 
                      group = "neotranscripts", 
                      height = 0.3, 
                      type = "transcripts"),
  "Cell line 1" = list(path = paste0(tracks_folder, "B69T10.against.201201_Ewing_neos_NEW_IDs.cov_per_base_final"), 
                   group = "tumor", 
                   color = "grey", 
                   height = 0.3, 
                   type = "density"),
  "Cell line 2" = list(path = paste0(tracks_folder, "B69T13.against.201201_Ewing_neos_NEW_IDs.cov_per_base_final"), 
                   group = "tumor", 
                   color = "grey", 
                   height = 0.3, 
                   type = "density"),
  "Tumor 1" = list(path = paste0(tracks_folder, "G312T05.against.201201_Ewing_neos_NEW_IDs.cov_per_base_final"), 
                 group = "tumor", 
                 color = "grey", 
                 height = 0.3, 
                 type = "density"),
  "Tumor 2" = list(path = paste0(tracks_folder, "G315T01.against.201201_Ewing_neos_NEW_IDs.cov_per_base_final"), 
                 group = "tumor", 
                 color = "grey", 
                 height = 0.3, 
                 type = "density"),
  "Tumor 3" = list(path = paste0(tracks_folder, "G402T05.against.201201_Ewing_neos_NEW_IDs.cov_per_base_final"), 
                   group = "tumor", 
                   color = "grey", 
                   height = 0.3, 
                   type = "density")
  # "All genes" = list(path = paste0(tracks_folder, "gencode.v19.annotation.gtf"), 
  #                    group = "neotranscripts", 
  #                    height = 0.4, 
  #                    type = "transcripts"),
  # "Pile-up" = list(path = paste0(tracks_folder, "EW_MERGED.against.201201_Ewing_neos_NEW_IDs.cov_per_base_final"), 
  #                  group = "pileup", 
  #                  sashimi_plus = paste0(tracks_folder, "junctions_plus/"),
  #                  sashimi_minus = paste0(tracks_folder, "junctions_minus/"),
  #                  color = "grey", 
  #                  height = 0.6, 
  #                  type = "sashimi", 
  #                  default_min_junctions = 200,
  #                  adapted_min_junctions = list()
  #)
)



## Genomic regions
genes <- sprintf("Ew_NG%s",c(1:4))
default_half_width <- 10000
half_widths <- list()

