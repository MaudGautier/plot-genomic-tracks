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
                 sashimi_plus = "/Users/maudgautier/Documents/data/tmp_from_calcsub/figs_EwS_100kb/4_junctions/junctions_100kb_Rscript_MATE1_SENSE/", 
                 sashimi_minus = "/Users/maudgautier/Documents/data/tmp_from_calcsub/figs_EwS_100kb/4_junctions/junctions_100kb_Rscript_MATE2_SENSE/", 
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


# Deal with case where repo not given (in this case, assumes that it is the working directory)
if (!exists(deparse(substitute(main_folder)))) {
  main_folder <- getwd()
}
source(paste0(main_folder, "/src/R/functions.R"))


import_libraries()

# Print warning if files do not exist
for (track_id in 1:length(tracks_info)) {
  if (!file.exists(tracks_info[[track_id]]$path)) {
    stop(paste("Track", names(tracks_info)[[track_id]], ": File", tracks_info[[track_id]]$path, "not found"))
    quit(save="ask")
  }
}


# Initialise
list_tables <- list()
list_txdb <- list()
list_gtf <- list()
list_groups <- lapply(tracks_info, function(x) {return(x$group)})
heights <- as.numeric(lapply(tracks_info, function(x) {return(x$height)}))
for (track_id in 1:length(tracks_info)) {
  track <- tracks_info[[track_id]]
  name <- names(tracks_info)[track_id]
  
  # Read track
  if (track$type == "density" || track$type == "GGAA" || track$type == "sashimi") {
    list_tables[[name]] <- read_track_table(track$path)
    tracks_info[[name]]$list_heights_group <- names(which(list_groups == track$group))
  } else if (track$type == "transcripts") {
    list_gtf[[name]] <- read.table(track$path, sep="\t")
    list_txdb[[name]] <- GenomicFeatures::makeTxDbFromGFF(file = track$path, format="gtf")
  } 
}

# Print warning if not same size
if (!check_dim_equality(list_tables[which(names(list_tables) == "GGAA")])) {
  print("WARNING: Dimensions of all files are not equal. Make sure the input files re correctly created before carrying on!")
}





# Create plots ------------------------------------------------------------

for (selected_gene in list_genes_neos_newIDs) {
  cat(paste0("\nCreating plot for ", selected_gene, "...\n"))
  
  
  # REFAIRE CA SANS TAB_GTF_NEOS + avoir longueur des genes dans la liste des genes + min_sashimi dans liste des genes
  # Get genomic range
  selected_gtf_lines <- tab_gtf_neos[which(tab_gtf_neos$V3 == "transcript" & grepl(paste0(selected_gene,";"), tab_gtf_neos$V9)),]
  # Increase extension for certain genes
  if (selected_gene %in% list_longers) {
    extension <- 50000
  } else {
    extension <- 10000
  }
  genom_range <- GenomicRanges::GRanges(selected_gtf_lines$V1, 
                                        IRanges::IRanges(selected_gtf_lines$V4 - extension, 
                                                         selected_gtf_lines$V5 + extension))
  
  
  
  # Make plot for each track
  list_plots <- list()
  for (track_id in 1:length(tracks_info)) {
    # Parameters
    track <- tracks_info[[track_id]]
    name <- names(tracks_info)[track_id]
    
    # Find maximal height for the group 
    max_height_group <- get_max_height(list_tables[tracks_info[[name]]$list_heights_group], 
                                       selected_gene, 50)
    
    # Add plot to list
    if (track$type == "transcripts") {
      list_plots[[name]] <- plot_transcripts(list_txdb[[name]], genom_range, "arrowrect")
    } else if (track$type == "density") {
      list_plots[[name]] <- plot_coverage(list_tables[[name]], selected_gene, max_height_group, track$color)
    } else if (track$type == "GGAA") {
      list_plots[[name]] <- plot_GGAA_repeats(list_tables[[name]], genom_range, selected_gene, track$color)
    } else if (track$type == "sashimi") {
      density_list = list()
      junction_list = list()
      tab_dt <- as.data.table(list_tables[[name]][,c("Position", "Coverage")])
      colnames(tab_dt) = c("x", "y")
      density_list[[name]] = tab_dt
  
      # Source junctions file
      if (gtf_lines$V7[1] == "-") { 
        source(paste0(track$sashimi_minus, '/', selected_gene, '_pileup.R')) 
      } else if (gtf_lines$V7[1] == "+") { 
        source(paste0(track$sashimi_plus, '/', selected_gene, '_pileup.R')) 
      }
      
      # Limits
      if (selected_gene %in% names(track$adapted_min_junctions)) { 
        min_sashimi <- track$adapted_min_junctions[[selected_gene]]
      } else { min_sashimi <- track$default_min_junctions }
      
      junction_list_subset = list()
      junction_list_subset[[name]] <- junction_list[[name]][which(junction_list[[name]]$count > min_sashimi),]
      list_plots[[name]] <- plot_sashimi(list_tables[[name]], junction_list_subset, selected_gene, density_list)
    } 
  }
  
  # Combine plots
  chromosome <- selected_gtf_lines$V1[1]
  png(paste0(output_folder_neos, selected_gene, ".png"), 
      width = 1800, height = 1050)
  print(ggbio::tracks(list_plots,
                      heights = heights,
                      label.text.cex = 1.5,
                      title = paste0(selected_gene, " (", chromosome, ")"),
                      xlim = c(min(genom_range@ranges@start), max(genom_range@ranges@start + genom_range@ranges@width))
  ) + ylab("") + scale_x_continuous(labels=comma) + theme(axis.text.x = element_text(size = 15))) 
   dev.off()

}

