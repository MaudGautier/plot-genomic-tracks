#!/usr/bin/env R

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                                                              #
# make_genomic_plot.R                                                                  #
#                                                                              #
# By: Maud Gautier <https://github.com/MaudGautier>, 2021                      #
#                                                                              #
# Broad permissions are granted to use, modify, and distribute this software   #
# as specified in the MIT License included in this distribution's LICENSE file.#
#                                                                              #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# This script allows to plot genomic tracks based on the parameters given in 
# input in a config file.


####### ANNOTATIONS POUR MOI #######

# Pour lancer depuis CLI
#Rscript /Users/maudgautier/Documents/github-savings/plot-genomic-tracks/src/R/make_genomic_plot.R



# CONFIG ------------------------------------------------------------------
## Config neotranscripts

## 1. Config : path des fichiers + noms associés + couleurs associées

# Name # File # Sashimi # height

## CONFIG
tracks_info <- list(
  "GGAA" = list(path = "/Users/maudgautier/Documents/data/tmp_from_calcsub/figs_EwS_100kb/3d_IGV_New_IDs/Against_genes_of_201201_Ewing_neos_NEW_IDs/hg19_GGAA_TTCC_20151221.against.201201_Ewing_neos_NEW_IDs.cov_per_base_final", 
                group = "GGAA", 
                color = "black", 
                height = 0.2, 
                type = "GGAA"),
  "EF_low" = list(path = "/Users/maudgautier/Documents/data/tmp_from_calcsub/figs_EwS_100kb/3d_IGV_New_IDs/Against_genes_of_201201_Ewing_neos_NEW_IDs/ASP14-d7-FLI1.against.201201_Ewing_neos_NEW_IDs.cov_per_base_final", 
                   group = "FLI1", 
                   color = "red", 
                   height = 0.3, 
                   type = "density"),
  "EF_high" = list(path = "/Users/maudgautier/Documents/data/tmp_from_calcsub/figs_EwS_100kb/3d_IGV_New_IDs/Against_genes_of_201201_Ewing_neos_NEW_IDs/FLI1_merged.against.201201_Ewing_neos_NEW_IDs.cov_per_base_final", 
                   group = "FLI1",
                   color = "red", 
                   height = 0.3, 
                   type = "density"),
  "H3K27ac" = list(path = "/Users/maudgautier/Documents/data/tmp_from_calcsub/figs_EwS_100kb/3d_IGV_New_IDs/Against_genes_of_201201_Ewing_neos_NEW_IDs/H3K27ac_merged.against.201201_Ewing_neos_NEW_IDs.cov_per_base_final", 
                   group = "H3K27", 
                   color = "green", 
                   height = 0.3, 
                   type = "density"),
  "H3K4me3" = list(path = "/Users/maudgautier/Documents/data/tmp_from_calcsub/figs_EwS_100kb/3d_IGV_New_IDs/Against_genes_of_201201_Ewing_neos_NEW_IDs/H3K4me3_merged.against.201201_Ewing_neos_NEW_IDs.cov_per_base_final", 
                   group = "H3K4", 
                   color = "blue", 
                   height = 0.3, 
                   type = "density"),
  "Tumor" = list(path = "/Users/maudgautier/Documents/data/tmp_from_calcsub/figs_EwS_100kb/3d_IGV_New_IDs/Against_genes_of_201201_Ewing_neos_NEW_IDs/G402T05.against.201201_Ewing_neos_NEW_IDs.cov_per_base_final", 
                   group = "tumor", 
                   color = "grey", 
                   height = 0.3, 
                   type = "density"),
  # "Cell line" = list(path = "/Users/maudgautier/Documents/data/tmp_from_calcsub/figs_EwS_100kb/3d_IGV_New_IDs/Against_genes_of_201201_Ewing_neos_NEW_IDs/B69T10.against.201201_Ewing_neos_NEW_IDs.cov_per_base_final", 
  #                group = "cl", 
  #                color = "grey", 
  #                height = 0.3, 
  #                type = "density"),
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
output_folder_neos <- "/Users/maudgautier/Documents/github-savings/plot-genomic-tracks/plots/"
main_folder <- "/Users/maudgautier/Documents/github-savings/plot-genomic-tracks/"
list_genes_neos_newIDs <- c("Ew_NG3")
list_longers <- c()
# idem : indiquer la taille de fenetre dans une liste des noms des genes
### Un autre config pour 13, 16 et 20 avec non pas tumor mais cell-line



## Mix
# gtf_file_neos <- '/Users/maudgautier/Documents/data/tmp_from_calcsub/figs_EwS_100kb/1_Data/201201_Ewing_neos_NEW_IDs.gtf'
# tab_gtf_neos <- read.table(gtf_file_neos, sep="\t")
# txdb_neos <- GenomicFeatures::makeTxDbFromGFF(file = gtf_file_neos, format="gtf")




## IN SCRIPT make_genomic_plot.R



# Prepare files -----------------------------------------------------------


# Deal with case where repo not given 
# (in this case, assumes that it is the working directory)
if (!exists(deparse(substitute(main_folder)))) {
  main_folder <- getwd()
}
source(paste0(main_folder, "/src/R/functions.R"))
source(config_file)

# Import libraries
import_libraries()

# Check that all files exist
for (track_id in 1:length(tracks_info)) {
  if (!file.exists(tracks_info[[track_id]]$path)) {
    stop(paste("Track", names(tracks_info)[[track_id]], 
               ": File", tracks_info[[track_id]]$path, "not found"))
    quit(save="ask")
  }
}

# Initialise lists and read tracks
cat("\nInitialising and reading tracks...\n")
list_tables <- list()
list_txdb <- list()
list_gtf <- list()
list_groups <- lapply(tracks_info, function(x) {return(x$group)})
heights <- as.numeric(lapply(tracks_info, function(x) {return(x$height)}))
for (track_id in 1:length(tracks_info)) {
  
  # Get track parameters
  track <- tracks_info[[track_id]]
  name <- names(tracks_info)[track_id]
  
  # Read track
  if (track$type %in% c("density", "GGAA", "sashimi")) {
    list_tables[[name]] <- read_track_table(track$path)
    in_same_group <- names(which(list_groups == track$group))
    tracks_info[[name]]$list_heights_group <- in_same_group
  } else if (track$type == "transcripts") {
    list_gtf[[name]] <- read.table(track$path, sep="\t")
    list_txdb[[name]] <- makeTxDbFromGFF(file = track$path, format="gtf")
  } 
}

# Verify that the dimensions of all tables are the same
if (!check_dim_equality(list_tables[which(names(list_tables) == "GGAA")])) {
  print("WARNING:")
  print("Dimensions of all files are not equal.")
  print("Make sure the input files are correctly created!")
}

# Read gtf genes
gtf_genes <- read.table(gtf_genes_path, sep="\t")
txdb_genes <- GenomicFeatures::makeTxDbFromGFF(file = gtf_genes_path, 
                                               format="gtf")





# Create plots ------------------------------------------------------------


for (gene in genes) {
  cat(paste0("\nCreating plot for ", gene, "...\n"))
  
  # Get genomic range
  gtf_lines <- gtf_genes[which(gtf_genes$V3 == "transcript" & 
                                 grepl(paste0(gene,";"), gtf_genes$V9)),]
  if (!is.null(half_widths[[gene]])) {
    half_width <- half_widths[[gene]]
  } else { half_width <- default_half_width }
  genom_range <- GRanges(gtf_lines$V1, IRanges(gtf_lines$V4 - half_width,
                                               gtf_lines$V5 + half_width))

  
  # Make plot for each track
  list_plots <- list()
  for (track_id in 1:length(tracks_info)) {
    # Parameters
    track <- tracks_info[[track_id]]
    name <- names(tracks_info)[track_id]
    
    # Find maximal height for the group 
    in_same_group <- tracks_info[[name]]$list_heights_group
    max_height_group <- get_max_height(list_tables[in_same_group], 
                                       gene, 50)
    
    # Add plot to list
    if (track$type == "transcripts") {
      list_plots[[name]] <- plot_transcripts(list_txdb[[name]], 
                                             genom_range, 
                                             "arrowrect")
      
    } else if (track$type == "density") {
      list_plots[[name]] <- plot_coverage(list_tables[[name]], 
                                          gene, 
                                          max_height_group, 
                                          track$color)
    
    } else if (track$type == "GGAA") {
      list_plots[[name]] <- plot_GGAA_repeats(list_tables[[name]], 
                                              genom_range, 
                                              gene,
                                              half_width,
                                              track$color)
      
    } else if (track$type == "sashimi") {
      
      # Select the file containing the junctions
      if (gtf_lines$V7[1] == "-") { 
        file_junctions <- paste0(track$sashimi_minus, '/', gene, '.R')
      } else if (gtf_lines$V7[1] == "+") { 
        file_junctions <- paste0(track$sashimi_plus, '/', gene, '.R')
      }
      
      # Select minimum number of junctions
      if (gene %in% names(track$adapted_min_junctions)) { 
        min_sashimi <- track$adapted_min_junctions[[gene]]
      } else { min_sashimi <- track$default_min_junctions }
      
      # Plot sashimi
      color_junction_arcs <- list()
      color_junction_arcs[[name]] <- "darkgrey"
      list_plots[[name]] <- plot_sashimi(list_tables[[name]], gene, 
                                         file_junctions, min_sashimi, 
                                         density_color = "black", 
                                         color_list = color_junction_arcs)
      
    } 
  }
  
  # Combine plots
  chromosome <- gtf_lines$V1[1]
  left <- min(genom_range@ranges@start)
  right <- max(genom_range@ranges@start + genom_range@ranges@width)
  png(paste0(output_folder, gene, ".png"), 
      width = 1800, height = 1050)
  print(ggbio::tracks(list_plots,
                      heights = heights,
                      label.text.cex = 1.5,
                      title = paste0(gene, " (", chromosome, ")"),
                      xlim = c(left, right)
  ) + ylab("") + scale_x_continuous(labels=comma) + 
    theme(axis.text.x = element_text(size = 15))) 
   dev.off()

}

