PADOG: Pathway Analysis with Down-weighting of Overlapping Genes
================================================================

This package implements a general purpose gene set 
analysis method called PADOG that downplays the importance of
genes that apear often accross the sets of genes to be
analyzed. The package provides also a benchmark for gene set
analysis methods in terms of sensitivity and ranking using 24
public datasets from KEGGdzPathwaysGEO package.

Tarca AL (2012). PADOG: Pathway Analysis with Down-weighting of Overlapping Genes (PADOG). R package version 1.8.0.

Installation
------------

-   the latest released version from Bioconductor:

    ``` r
    source("http://bioconductor.org/biocLite.R")
    biocLite("PADOG")
    ```

-   the latest development version from github:

    ``` r
    devtools::install_github("pacifly/PADOG") 
    ```

You may also want to install the data package: `biocLite("KEGGdzPathwaysGEO")`

Documentation
-------------

To view documentation for the version of this package installed in your system, start R and enter: 
`browseVignettes("PADOG")`


