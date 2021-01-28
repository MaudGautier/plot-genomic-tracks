#!/usr/bin/env R

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                                                              #
# config_plot.R                                                                #
#                                                                              #
# By: Maud Gautier <https://github.com/MaudGautier>, 2021                      #
#                                                                              #
# Broad permissions are granted to use, modify, and distribute this software   #
# as specified in the MIT License included in this distribution's LICENSE file.#
#                                                                              #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# This a configuration file to plot genomic tracks.



# Usage -------------------------------------------------------------------

# From the CLI:
# Rscript src/R/make_genomic_plot.R example/config_plot.R




# Configuration for EwS neotranscripts ------------------------------------


## General parameters
output_folder <- "./plots/"
main_folder <- "./"
gtf_genes_path <- './test/201201_Ewing_neos_NEW_IDs.gtf'


## Tracks
tracks_info <- list(
  "GGAA" = list(path = "./test/hg19_GGAA_TTCC_20151221.against.201201_Ewing_neos_NEW_IDs.cov_per_base_final", 
                group = "GGAA", 
                color = "black", 
                height = 0.2, 
                type = "GGAA"),
  "EF_low" = list(path = "./test/ASP14-d7-FLI1.against.201201_Ewing_neos_NEW_IDs.cov_per_base_final", 
                   group = "FLI1", 
                   color = "red", 
                   height = 0.3, 
                   type = "density"),
  "Prediction" = list(path = "./test/201201_Ewing_neos_NEW_IDs.gtf", 
                      group = "neotranscripts", 
                      height = 0.4, 
                      type = "transcripts"),
  "Pile-up" = list(path = "./test/EW_MERGED.against.201201_Ewing_neos_NEW_IDs.cov_per_base_final", 
                 group = "pileup", 
                 sashimi_plus = "./test/junctions_plus/", 
                 sashimi_minus = "./test/junctions_minus/", 
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
genes <- c("Ew_NG1", "Ew_NG2", "Ew_NG3", "Ew_NG4", "Ew_NG9")
default_half_width <- 10000
half_widths <- list("Ew_NG9" = 50000)



