#!/usr/bin/env R

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                                                              #
# functions.R                                                                  #
#                                                                              #
# By: Maud Gautier <https://github.com/MaudGautier>, 2021                      #
#                                                                              #
# Broad permissions are granted to use, modify, and distribute this software   #
# as specified in the MIT License included in this distribution's LICENSE file.#
#                                                                              #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# This script contains the functions required to plot genomic tracks



# Import necessary libraries ----------------------------------------------
import_libraries <- function() {
  
  # Generic
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(ggbio))
  suppressPackageStartupMessages(library(scales))
  
  # For genes tracks
  suppressPackageStartupMessages(library(GenomicFeatures))
  
  # For sashimi tracks
  suppressPackageStartupMessages(library(grid))
  suppressPackageStartupMessages(library(gridExtra))
  suppressPackageStartupMessages(library(data.table))
  suppressPackageStartupMessages(library(gtable))
  suppressPackageStartupMessages(library(IRanges))
  suppressPackageStartupMessages(library(GenomicRanges))
}


# Plots read coverage along one selected gene -----------------------------
plot_coverage <- function(table_cov, 
                          selected_gene, 
                          max_height = NA, 
                          color = "grey", 
                          type = "bar") {
  
  # Select gene
  sub_table_cov <- table_cov[which(table_cov$Gene_ID == selected_gene), ]
  
  # Plot
  plot_cov <- ggplot(sub_table_cov, aes(x = Position, y = Coverage))
  
  # If bar
  if (type=="bar") {
    plot_cov <- plot_cov + 
      geom_bar(stat="identity", position="identity", fill = color)
  } else if (type == "line") {
    plot_cov <- plot_cov + 
      geom_line()
  }
  
  # Add labels
  plot_cov <- plot_cov +
    theme_classic() + 
    labs(x = '', y = "Coverage") +
    theme(text = element_text(size = 13))
  
  # Add max height
  if (!is.na(max_height)) { 
    plot_cov <- plot_cov +
      scale_y_continuous(limits = c(0, max_height))
  }
  
  # Return plot
  return(plot_cov)
  
}


# Plots all transcripts of one selected gene ------------------------------
plot_transcripts <- function(txdb, 
                             genom_range, 
                             geom_type = "rect") {
  
  # Plot transcripts
  trans_plot <- try(ggbio::autoplot(txdb, 
                                    which=genom_range, 
                                    names.expr = "gene_id", 
                                    range.geom = geom_type) + 
                      theme_classic())
  
  # In case of error (i.e. when no transcript in the genomic region),
  # plot a blank graphic
  if (class(trans_plot) == "try-error") {
    left <- min(genom_range@ranges@start)
    right <- max(genom_range@ranges@start + genom_range@ranges@width)
    
    trans_plot <- ggplot(data.frame(Position = seq(left, right), 
                                    Coverage = rep(0, right-left+1)), 
                         aes(x = Position, y = Coverage)) +
      theme_classic() + 
      theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(),
            axis.ticks.x = element_blank(), axis.text.x = element_blank())
  }
  
  return (trans_plot)
}


# Get maximum height of a given list of subtables -------------------------
get_max_height <- function(list_tables, 
                           selected_gene, 
                           default_max = 0) {
  
  # Get subsets for all tables
  new_max <- default_max
  for (tab in list_tables) {
    sub_tab <- tab[which(tab$Gene_ID == selected_gene), ]
    new_max <- max(sub_tab$Coverage, new_max)
  }
  
  # Return maximum
  return(new_max)
  
}


# Checks if dimensions are equal for all tables in a given list -----------
check_dim_equality <- function(list_tables) {
  for (t1 in list_tables) {
    for (t2 in list_tables) {
      if (dim(t1)[1] != dim(t2)[1] || dim(t1)[2] != dim(t2)[2]) {
        return(FALSE)
      }
    }
  }
  return(TRUE)
}


# Select tables with top coverage -----------------------------------------
# NB: This works *ONLY* if ONE table needs to be removed (not multiple)
select_tables <- function(list_tables_cov, 
                          selected_gene) {
  
  # Select gene
  list_means <- NULL
  for (tab in list_tables_cov) {
    sub_tab <- tab[which(tab$Gene_ID == selected_gene), ]
    list_means <- c(list_means, mean(sub_tab$Coverage))
  }
  # Only remove one
  # Need to loop if more than one is to be removed
  list_tables_cov[which.min(list_means)] <- NULL
  
  # Return list
  return(list_tables_cov)
  
}


# Defines line width as proportional to number of junctions ---------------
scale_lwd <- function(r) {
  lmin = 1
  lmax = 10
  return( r*(lmax-lmin)+lmin )
}


# Plots sashimi track for one gene ----------------------------------------
# NB: Modified from the original https://github.com/guigolab/ggsashimi
plot_sashimi <- function(table_cov, 
                         selected_gene, 
                         file_junctions,
                         min_sashimi,
                         density_color = "black", 
                         color_list = list("Pile-up"="darkgrey")) {
  
  # Fix problems with ggplot2 vs >3.0.0
  if(packageVersion('ggplot2') >= '3.0.0'){
    vs = 1
  } else {
    vs = 0
  }
  
  # Check that file_junctions exists
  if (!file.exists(file_junctions)) {
    stop(paste("Junction file", file_junctions, "not found. Please provide the correct path."))
    quit(save="ask")
  }
  
  tab_dt <- as.data.table(table_cov[,c("Position", "Coverage")])
  colnames(tab_dt) = c("x", "y")
  
  density_list = list()
  density_list[[name]] = tab_dt
  
  all_junction_list = list()
  all_junction_list[[name]] <- source(file_junctions)$value
  
  junction_list = list()
  junction_list[[name]] <- all_junction_list[[name]][which(all_junction_list[[name]]$count > min_sashimi),]
  
  
  # Initialise lists
  density_grobs = list()
  sub_table_cov <- table_cov[which(table_cov$Gene_ID == selected_gene), ]
  plot_cov <- ggplot(sub_table_cov, aes(x = Position, y = Coverage))
  maxheight <- max(sub_table_cov$Coverage)
  
  # Add all junctions
  for (bam_index in 1:length(density_list)) {
    
    # Get parameters
    id = names(density_list)[bam_index]
    d = data.table(density_list[[id]])
    junctions = data.table(junction_list[[id]])
    
    # Plot density (no junction yet)
    plot_cov <- plot_cov + geom_bar(stat="identity", 
                                    position="identity", 
                                    fill = density_color)
    
    # Aggregate junction counts
    row_i = c()
    if (nrow(junctions) >0 ) {
      junctions$jlabel = as.character(junctions$count)
      junctions = setNames(junctions[,.(max(y), 
                                        max(yend),
                                        round(mean(count)),
                                        paste(jlabel,collapse=",")
      ), keyby=.(x,xend)],
      names(junctions))
      
      # The number of rows (unique junctions per bam) has to be 
      # calculated after aggregation
      row_i = 1:nrow(junctions)
    }
    
    # Add each junction to plot
    for (i in row_i) {
      j_tot_counts = sum(junctions[['count']])
      j = as.numeric(junctions[i,1:5])
      
      # Find intron midpoint
      xmid = round(mean(j[1:2]), 1)
      ymid = max(j[3:4]) * 1.2
      
      # Thickness of the arch
      lwd = scale_lwd(j[5]/j_tot_counts)
      curve_par = gpar(lwd=lwd, col=color_list[[id]])
      
      # Arc grobs
      # Choose position of the arch (top or bottom)
      nss = i
      # Plot every other junction on the bottom (even)
      if (nss%%2 == 0) {  
        ymid = -0.3 * ymid
        
        # Draw the arcs
        # Left
        curve = xsplineGrob(x=c(0, 0, 1, 1), 
                            y=c(1, 0, 0, 0), 
                            shape=1, gp=curve_par)
        plot_cov = plot_cov + 
          annotation_custom(grob = curve, j[1], xmid, 0, ymid)
        
        # Right
        curve = xsplineGrob(x=c(1, 1, 0, 0), 
                            y=c(1, 0, 0, 0), 
                            shape=1, gp=curve_par)
        plot_cov = plot_cov + 
          annotation_custom(grob = curve, xmid, j[2], 0, ymid)
      }
      
      # Plot every other junction on the top (odd)
      if (nss%%2 != 0) {
        
        # Draw the arcs
        # Left
        curve = xsplineGrob(x=c(0, 0, 1, 1), 
                            y=c(0, 1, 1, 1), 
                            shape=1, gp=curve_par)
        plot_cov = plot_cov + 
          annotation_custom(grob = curve, j[1], xmid, j[3], ymid)
        
        # Right
        curve = xsplineGrob(x=c(1, 1, 0, 0), 
                            y=c(0, 1, 1, 1), 
                            shape=1, gp=curve_par)
        plot_cov = plot_cov + 
          annotation_custom(grob = curve, xmid, j[2], j[4], ymid)
      }
      
      # Add junction labels
      plot_cov = plot_cov + 
        annotate("label", x = xmid, y = ymid, 
                 label = as.character(junctions[i,6]),
                 vjust=0.5, hjust=0.5, 
                 label.padding=unit(0.01, "lines"),
                 label.size=NA, size=4
        )
      
    }
  }
  
  # Add theme and axis labels
  plot_cov = plot_cov +
    theme_classic() +
    labs(x = '', y = "Coverage") +
    theme(text = element_text(size = 13))
  
  # Return sashimi plot
  return(plot_cov)
  
}


# Create GGAA-repeats plot ------------------------------------------------
plot_GGAA_repeats <- function(tb, 
                              genom_range, 
                              selected_gene,
                              half_width,
                              color = "black") {
  
  # Plot only if >= 4 consecutive GGAA repeats
  if (max(tb[which(tb$Coverage >= 4 & 
                   tb$Gene_ID == selected_gene),]$Coverage) != -Inf) {
    
    # Fix: avoid invisible bars when the plotting window is very wide
    left <- min(genom_range@ranges@start)
    right <- max(genom_range@ranges@start + genom_range@ranges@width)
    if (right - left < 3*half_width) {
      plot_GGAA <- plot_coverage(tb[which(tb$Coverage >= 4 | 
                                            tb$Coverage == 0),], 
                                 selected_gene, NA, color, "bar")
    } else {
      plot_GGAA <- plot_coverage(tb[which(tb$Coverage >= 4 | 
                                            tb$Coverage == 0),], 
                                 selected_gene, NA, color, "line") 
    }
  } else {
    # Blank plot if fewer than 4 consecutive repeats
    plot_GGAA <- plot_coverage(tb, 
                               selected_gene, NA, "white", "bar") 
  }
  
  # Return GGAA plot
  return(plot_GGAA)
}


# Read track table --------------------------------------------------------
read_track_table <- function(track_file) {
  tab_track <- read.table(track_file, header=T, comment = "")
}


# -------------------------------------------------------------------------
