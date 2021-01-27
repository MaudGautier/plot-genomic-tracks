# plot-genomic-tracks

This repo contains scripts that will allow you to automatically create plots of genomic tracks at defined genomic regions in an IGV-like fashion.


## Description


## Usage


First, clone the repository
```
git clone https://github.com/MaudGautier/plot-genomic-tracks.git ./
cd ./plot-genomic-tracks
```

Second, prepare a configuration file in a similar fashion to that located in the [``example/`` directory](https://github.com/MaudGautier/plot-genomic-tracks/tree/main/example/EwS).

Once this is done, you can plot your tracks using the following command line:
```
Rscript ./src/R/make_genomic_plot.R your_config_plot.R
```




## Requirements

Hereunder is the list of applications that are necessary:

* R >= 4.0.3


In R, you *must* install the following libraries:

* ggplot2
* ggbio
* scales
* GenomicFeatures
* grid
* gridExtra
* data.table
* gtable
* IRanges
* GenomicRanges

These libraries can be installed as follows (from the R CLI):
```
# Base packages
install.packages("ggplot2")
install.packages("scales")
install.packages("grid")
install.packages("gridExtra")
install.packages("data.table")
install.packages("gtable")

# BiocManager
install.packages("BiocManager")

# Installation with BiocManager
BiocManager::install("ggbio")
BiocManager::install("GenomicFeatures")
BiocManager::install("IRanges")
BiocManager::install("GenomicRanges")
```



## TL;DR
### Plotting tracks

To plot the tracks, you must prepare the configuration file.
The configuration file ([see an example here](https://github.com/MaudGautier/plot-genomic-tracks/tree/main/example/config_plot.R)) must contain at least:

1. The list of tracks to plot and their characteristics
2. The list of genomic regions
3. The path to the output folder containing the plots







