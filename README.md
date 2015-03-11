PADOG: Pathway Analysis with Down-weighting of Overlapping Genes
================================================================

This package implements a general purpose gene set 
analysis method called PADOG that downplays the importance of
genes that apear often accross the sets of genes to be
analyzed. The package provides also a benchmark for gene set
analysis methods in terms of sensitivity and ranking using 24
public datasets from KEGGdzPathwaysGEO package.

[Tarca, Adi L., et al. "Down-weighting overlapping genes improves gene set analysis." BMC bioinformatics 13.1 (2012): 136.](http://www.biomedcentral.com/1471-2105/13/136)

Installation
------------

-   the latest released version from [**Bioconductor**](http://www.bioconductor.org/packages/release/bioc/html/PADOG.html):

    ``` r
    source("http://bioconductor.org/biocLite.R")
    biocLite("PADOG")
    ```

-   the latest development version from github:

    ``` r
    ups = installed.packages()[,"Package"]
    biocp = c("KEGGdzPathwaysGEO", "KEGG.db", "limma", "AnnotationDbi", 
              "Biobase", "hgu133plus2.db", "hgu133a.db")
    todos = setdiff(biocp, ups)
    if (length(todos) > 0) {
        source("http://bioconductor.org/biocLite.R")
        biocLite(todos)
    }

    if (! "devtools" %in% ups) install.packages("devtools")
    devtools::install_github("pacifly/PADOG") 
    ```

Documentation
-------------

To view documentation for the version of this package installed in your system, start R and enter: 
`browseVignettes("PADOG")`


