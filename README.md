# scMitoMut: a R package for mitochondiral mutation analysis in single cell omics data

## Installation

You can install the pacakge from source code:

```r
# install.packages("devtools")
devtools::install_github("TeamPerie/scMitoMut")
```

## Example

This is a basic example which shows you the main function of the package:

```r
library(scMitoMut)

# load the data
## Use the example data
f = system.file("extdata", "mini_dataset.tsv.gz", package = "scMitoMut")

## Load the data with parse_table function
f_h5 = parse_table(f, sep = "\t", h5_file = "./mut.h5")

## open the h5f file
x = open_h5_file(f_h5)

# run the model fit
# You can increase the cpu core to accelerate
# This step need some time, so the result will be kept in h5 file,
# you do not need to re-run this step, when you load the h5 file next time.
run_model_fit(x, mc.cores = 1)

# Filter the loci based on the model fit results
# The filter options will be keeped in the object by memory
# Next time you re-load the h5 file, the filter will be initiated as default
x = filter_loc(x, 
    min_cell = 10, 
    model = "bb", 
    p_threshold = 0.01, 
    p_adj_method = "fdr"
    )
x

# Set the cell annotation
f = system.file("extdata", "mini_dataset_cell_ann.csv", package = "scMitoMut")
cell_ann = read.csv(f, row.names=1)
# Prepare the color for cell annotation
colors = c(
    "Cancer Epi" = "#f28482",
    Blood = "#f6bd60")
ann_colors = list("SeuratCellTypes" = colors)

# plot the heatmap for binary mutation
plot_heatmap(x, type = "binary", cell_ann = cell_ann, ann_colors = ann_colors, percent_interp = 0.2)

# plot the heatmap for p-value
plot_heatmap(x, type = "p", cell_ann = cell_ann, ann_colors = ann_colors, percent_interp = 0.2)

# plot the heatmap for allele frequency
plot_heatmap(x, type = "af", cell_ann = cell_ann, ann_colors = ann_colors, percent_interp = 0.2)

# check af~coverage for one loci
plot_af_coverage(x, "chrM.200")
```



