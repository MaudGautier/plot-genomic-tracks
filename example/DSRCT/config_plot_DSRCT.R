#!/usr/bin/env R

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                                                              #
# config_plot_DSRCT.R                                                          #
#                                                                              #
# By: Maud Gautier <https://github.com/MaudGautier>, 2021                      #
#                                                                              #
# Broad permissions are granted to use, modify, and distribute this software   #
# as specified in the MIT License included in this distribution's LICENSE file.#
#                                                                              #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# This is the configuration file used to plot DSRCT tumors.



# Usage -------------------------------------------------------------------

# From the CLI:
# Rscript src/R/make_genomic_plot.R example/DSRCT/config_plot_DSRCT.R




# Configuration for EwS neotranscripts ------------------------------------


## General parameters
output_folder <- "/Users/maudgautier/Documents/github-savings/plot-genomic-tracks/plots/DSRCT/"
main_folder <- "/Users/maudgautier/Documents/github-savings/plot-genomic-tracks/"
gtf_genes_path <- "/Users/maudgautier/Documents/github-savings/plot-genomic-tracks/input/DSRCT/200624_DSRCT_neos_NEW_IDs.gtf"
tracks_folder <- "/Users/maudgautier/Documents/github-savings/plot-genomic-tracks/input/DSRCT/"

## Tracks
tracks_info <- list(
  "EWS_WT1" = list(path = paste0(tracks_folder, "ChIP-EWS-WT1.against.200624_DSRCT_neos_NEW_IDs.cov_per_base_final"), 
                group = "WT1", 
                color = "red", 
                height = 0.3, 
                type = "density"),
  "PolII" = list(path = paste0(tracks_folder, "ChIP-PolII.against.200624_DSRCT_neos_NEW_IDs.cov_per_base_final"), 
                group = "PolII", 
                color = "blue", 
                height = 0.3, 
                type = "density"),
  "Tumor" = list(path = paste0(tracks_folder, "G315T08.against.200624_DSRCT_neos_NEW_IDs.cov_per_base_final"), 
                group = "Tumor", 
                color = "grey", 
                height = 0.3, 
                type = "density"),
  "Prediction" = list(path = paste0(tracks_folder, "200624_DSRCT_neos_NEW_IDs.gtf"), 
                      group = "neotranscripts", 
                      height = 0.4, 
                      type = "transcripts"),
  "All genes" = list(path = paste0(tracks_folder,"gencode.v19.annotation.gtf"),
                      group = "neotranscripts", 
                      height = 0.4, 
                      type = "transcripts"),
  "Pile-up" = list(path = paste0(tracks_folder, "DSRCT_merged.against.200624_DSRCT_neos_NEW_IDs.cov_per_base_final"), 
                 group = "pileup", 
                 sashimi_plus = paste0(tracks_folder, "junctions_plus/"),
                 sashimi_minus = paste0(tracks_folder, "junctions_minus/"),
                 color = "grey", 
                 height = 0.6, 
                 type = "sashimi", 
                 default_min_junctions = 200,
                 adapted_min_junctions = list(
                   "DSRCT_NG5" = 5000,
                   "DSRCT_NG21" = 500,
                   "DSRCT_NG22" = 1000,
                   "DSRCT_NG28" = 500,
                   "DSRCT_NG29" = 500,
                   "DSRCT_NG36" = 2000,
                   "DSRCT_NG37" = 2000
                 )
  )
)



## Genomic regions
genes <- sprintf("DSRCT_NG%s",c(1:37))
default_half_width <- 10000
half_widths <- list("DSRCT_NG2" = 50000,
					"DSRCT_NG4" = 50000,
					"DSRCT_NG5" = 50000,
					"DSRCT_NG7" = 50000,
					"DSRCT_NG8" = 50000,
					"DSRCT_NG9" = 50000,
					"DSRCT_NG10" = 50000,
					"DSRCT_NG12" = 50000,
					"DSRCT_NG13" = 50000,
					"DSRCT_NG14" = 50000,
					"DSRCT_NG15" = 50000,
					"DSRCT_NG16" = 50000,
					"DSRCT_NG18" = 50000,
					"DSRCT_NG19" = 50000,
					"DSRCT_NG21" = 50000,
					"DSRCT_NG22" = 50000,
					"DSRCT_NG23" = 50000,
					"DSRCT_NG24" = 50000,
					"DSRCT_NG25" = 50000,
					"DSRCT_NG26" = 50000,
					"DSRCT_NG27" = 50000,
					"DSRCT_NG28" = 50000,
					"DSRCT_NG29" = 50000,
					"DSRCT_NG32" = 50000,
					"DSRCT_NG33" = 50000,
					"DSRCT_NG34" = 50000,
					"DSRCT_NG35" = 50000,
					"DSRCT_NG36" = 50000,
					"DSRCT_NG37" = 50000)

