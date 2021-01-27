#!/usr/bin/env R

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                                                              #
# make_genomic_plot.R                                                          #
#                                                                              #
# By: Maud Gautier <https://github.com/MaudGautier>, 2021                      #
#                                                                              #
# Broad permissions are granted to use, modify, and distribute this software   #
# as specified in the MIT License included in this distribution's LICENSE file.#
#                                                                              #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# This script allows to plot genomic tracks based on the parameters given in 
# input in a config file.



# Parse arguments ---------------------------------------------------------


args <- commandArgs(trailingOnly = TRUE)
config_file <- args[1]




# Prepare files -----------------------------------------------------------


# Deal with case where repo not given 
# (in this case, assumes that it is the working directory)
if (!exists(deparse(substitute(main_folder)))) {
  main_folder <- getwd()
}
source(file.path(main_folder, "/src/R/functions.R"))
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

# Check that output folder exists
if (!dir.exists(output_folder)) {
  dir.create(output_folder, showWarnings = F)
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
        file_junctions <- file.path(track$sashimi_minus, paste0(gene, '.R'))
      } else if (gtf_lines$V7[1] == "+") { 
        file_junctions <- file.path(track$sashimi_plus, paste0(gene, '.R'))
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
  png(file.path(output_folder, paste0(gene, ".png")), 
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

