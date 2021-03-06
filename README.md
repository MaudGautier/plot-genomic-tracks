# plot-genomic-tracks

This repo contains scripts that will allow you to automatically create plots of genomic tracks at defined genomic regions in an IGV-like fashion.



## Usage


**Preliminary step: clone the code**

Clone the repository with this command line:
```
git clone https://github.com/MaudGautier/plot-genomic-tracks.git ./
cd ./plot-genomic-tracks
```


**Step 1: Parse input files**

Parse the input files by preparing a configuration file in a similar fashion to the [``example/config_prep.sh`` file](https://github.com/MaudGautier/plot-genomic-tracks/tree/main/example/config_prep.sh) and executing it using this command line:
```
sh your_config_prep_file.sh
```



**Step 2: Plot tracks**

Prepare a configuration file in a similar fashion to the [``example/config_plot.R`` file](https://github.com/MaudGautier/plot-genomic-tracks/tree/main/example/config_plot.R).

Once this is done, you can plot your tracks using the following command line:
```
Rscript ./src/R/make_genomic_plot.R your_config_plot.R
```




## Requirements

Hereunder is the list of applications that are necessary and versions I used:

* R, version 4.0.3
* Python, version 3.7.6
* GNU bash, version 4.2.46
* GNU Awk, version 4.0.2
* GNU sort, version 8.22
* GNU grep, version 2.20
* bedtools, version 2.27.1
* bigWigToBedGraph, ucsc-bigwigtobedgraph version 357


For the plotting step, in R, you *must* install the following libraries:

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

### Description

To plot genomic tracks, two successive steps are required: first, the files containing the tracks to be plotted must be parsed adequately and second, the plots can be generated these parsed files.


### Step 1: Parse input files

To parse the input files, you must prepare the configuration file and execute it.
The configuration file ([see an example here](https://github.com/MaudGautier/plot-genomic-tracks/tree/main/example/config_prep.R)) proceeds as follows:

* If there is no existing BED file defining the span of the regions to plot, a GTF file can be provided and it will be parsed to create a correct BED file with a given extended spanning on both the left and right sides
* The BAM or BED files corresponding to all tracks to be plotted as density or sashimi tracks are parsed to create a per-base coverage file (containing the sequencing depth at each position between any designed BED region defined precedently)
* For each BAM file that is to be plotted as a sashimi track, an R script containing the junctions' arcs positions is created
* Tracks containing repeats (like GGAA repeats) provided in the BIGWIG format are parsed to create a per-base height file (where the height given at each position equals the number of repeats of the 4-sequence motif)



### Step 2: Plot tracks

To plot the tracks, you must prepare the configuration file.
The configuration file ([see an example here](https://github.com/MaudGautier/plot-genomic-tracks/tree/main/example/config_plot.R)) must contain at least:

1. The list of tracks to plot and their characteristics
2. The list of genomic regions
3. The path to the output folder containing the plots


#### List of tracks

Each track must be given a number of characteristics associated to it.
These include:

* (Required) `path`: the path to the file
* `group`: a string or numeric indicating the group to which the track belongs. This piece of information is used to find the scale on which certain tracks must be based (case where multiple tracks must be shown on the same scale to be compared to one another)
* `color`: the color of the track
* `height`: the relative height of the track when combined with all of the tracks
* (Required) `type`: the category of the track. This defines the actions that will be taken on that track to plot it. The categories can be: 
	* `'density'` to plot only the per-base coverage of each file
	* `'sashimi'` to add the junctions
	* `'transcripts'` to draw the transcripts


**Additional variables for `sashimi` types:**

* `sashimi_plus`: path to the R script file containing the junction reads located on the forward strand
* `sashimi_minus`: path to the R script file containing the junction reads located on the reverse strand
* `default_min_junctions`: the minimum number of junction reads required to draw the junction
* `adapted_min_junctions`: a list of minimum number of junctions for certain specific genes (only necessary when different from the `default_min_junctions` value)


**Format of input track files**

There are several types of input files that can be used as tracks. They all must be in one of the following formats:

* Per-base coverage plots (`type = 'density'`): 
```
#Chromosome  Position  Gene_ID  Coverage
chr1         73121083  Ew_NG1   3
chr1         73121084  Ew_NG1   3
chr1         73121085  Ew_NG1   3
chr1         73121086  Ew_NG1   3
```

* GTF files containing a list of transcripts to plot (`type = 'transcripts'`):
```
chr1   scallop  transcript  73171083   73221254   1000  +  .  gene_id  "Ew_NG1";   transcript_id  "Ew_NG1.1";   RPKM  "1.7032";  cov  "30.4365";
chr1   scallop  exon        73171083   73171273   1000  +  .  gene_id  "Ew_NG1";   transcript_id  "Ew_NG1.1";   exon  "1";
chr1   scallop  exon        73174742   73174919   1000  +  .  gene_id  "Ew_NG1";   transcript_id  "Ew_NG1.1";   exon  "2";
chr1   scallop  exon        73219351   73221254   1000  +  .  gene_id  "Ew_NG1";   transcript_id  "Ew_NG1.1";   exon  "3";
```

* R scripts containing a data.frame of junctions positions (`type = 'sashimi'`):
```
data.frame(x=c(17297203), xend=c(17310080), y=c(2), yend=c(3), count=c(2))
```




#### List of genomic regions

Each genomic region to be plotted must be listed in the `genes` variable by its ID.
In addition, a few pieces of information should be indicated:

* `default_half_width`: The width that must be spanned on the 5' and 3' sides of the transcripts
* `half_widths`: A list of adapted widths for certain specific genes (only necessary when different from the `default_half_width` value)





## TODO

- [ ] Allow to give BED file rather than GTF file to define regions (step 2 -- R plotting)
- [ ] Add a test folder with infos from 1 small gene
- [ ] Add an option for BIGWIG files
- [ ] Add possibility to have combined plots with horizontal label
- [ ] Add check that all names given to tracks_info are different from one another

