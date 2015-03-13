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

Usage
------

-   To run PADOG gene set analysis:
    
    ``` r
    ### run padog on a colorectal cancer dataset of the 24 benchmark datasets
    set = "GSE9348"
    data(list=set, package="KEGGdzPathwaysGEO")
    x = get(set)
    # extract from the dataset the required info
    exp = experimentData(x)
    dataset = exp@name
    dat.m = exprs(x)
    ano = pData(x)
    design = notes(exp)$design
    annotation = paste(x@annotation,".db",sep="")
    targetGeneSets = notes(exp)$targetGeneSets

    # run without parallel
    myr = padog(esetm = dat.m, group = ano$Group, paired = design=="Paired",
        block = ano$Block, targetgs = targetGeneSets, annotation = annotation,
        gslist = "KEGG.db", organism = "hsa", verbose = TRUE, Nmin = 3, NI = 200,
        plots = FALSE, dseed = 1)

    head(myr$res, 20)

    # run with parallel
    myr2=padog(esetm = dat.m, group = ano$Group, paired = design=="Paired",
        block = ano$Block, targetgs = targetGeneSets, annotation = annotation,
        gslist = "KEGG.db", organism = "hsa", verbose = TRUE, Nmin = 3, NI = 200,
        plots = FALSE, dseed = 1, paral = TRUE, ncr = 2)

    # verify that the result is the same which is a built-in feature
    all.equal(myr, myr2)
    ```

-   To compare your gene set analysis method with those implemented in PADOG on benchmark datasets:
    
    Write a "constructor" function to return your gene set analysis method; here we construct
    some artificial methods that assign uniform (or skewed Beta distributed) p values to gene sets.

    ``` r
    setF = function(type=c("pos","uni","neg"), NI=1000, dseed=1, FDRmeth=c("BH","Permutation","holm")) {
        ### type -- uni for uniform p values, pos / neg for excessive false positive / negative 
        type = match.arg(type)
        FDRmeth = match.arg(FDRmeth)
        return(function(set, mygslist, minsize) {
        # define your new gene set analysis method that takes as input:
        # set -- the name of dataset file from the PADOGsetspackage
        # mygslist -- a list of the genesets
        # minsize -- minimum number of genes in a geneset to be considered for analysis 
            set.seed(dseed)
            data(list=set, package="KEGGdzPathwaysGEO", envir=environment())
            x = get(set)

            # extract from the dataset the required info
            exp = experimentData(x)
            dat.m = exprs(x)
            ano = pData(x)
            dataset = exp@name
            design = notes(exp)$design
            annotation = paste(x@annotation,".db",sep="")
            targetGeneSets = notes(exp)$targetGeneSets

            # get rid of duplicates probesets per ENTREZ ID
            aT1 = filteranot(esetm=dat.m, group=ano$Group, paired=(design == "Paired"),
                block=ano$Block, annotation=annotation)
                
            # remove from gene sets the genes absent in the current dataset
            mygslist = lapply(mygslist, function(z) z[z %in% (aT1$ENTREZID)])
            # drop gene sets with less than minsize genes in the current dataset
            mygslist = mygslist[sapply(mygslist, length) >= minsize]
            
            # this toy method simulate some artificial p values for gene sets
            pres = switch(type,
                        uni = runif(length(mygslist)),
                        pos = rbeta(length(mygslist), 1, 5),
                        neg = rbeta(length(mygslist), 5, 1)
            )
            res = data.frame(ID=names(mygslist), P=pres, stringsAsFactors=FALSE) 
                            
            # obtain null p value matrix (i.e. applying your gene set analysis 
            # method on permuted group labels for NI times) that will be used
            # in estimating fdr if FDRmeth = "Permutation"
            nullp = switch(type,
                        uni = matrix(runif(nrow(res) * NI), ncol = NI),
                        pos = matrix(rbeta(nrow(res) * NI, 1, 5), ncol = NI),
                        neg = matrix(rbeta(nrow(res) * NI, 5, 1), ncol = NI)
            )
            rownames(nullp) = res$ID
            
            # multiple test adjustment        
            res$FDR = switch(FDRmeth,
                            p.adjust(res$P, FDRmeth),
                            Permutation = PADOG:::getFDR(nullp, res$P)
            )
            
            # compute normalized ranks
            res$Rank = rank(res$P, na.last="keep") / sum(!is.na(res$P)) * 100
            # record method name and dataset name
            res$Method = type
            res$Dataset = dataset

            # return relevant gene sets result and null p values
            return(list(targ = res[res$ID %in% targetGeneSets,], pval = nullp))
        }
        )
    }

    mymnam = c("pos","uni","neg")
    mym = lapply(mymnam, setF, NI=200, dseed=1, FDRmeth="P")
    names(mym) = mymnam

    # compare your methods with existing ones
    out = compFDR(datasets=NULL, existingMethods=c("AbsmT", "GSA", "PADOG"),
        mymethods=mym, gslist="KEGG.db", Nmin=3, NI=200, Npsudo=20, FDRmeth="P", 
        plots=TRUE, verbose=TRUE, use.parallel=TRUE, dseed=1, pkgs=NULL)
    ```

    Here is an output figure from above running:


Documentation
-------------

To view documentation for the version of this package installed in your system, start R and enter: 
`browseVignettes("PADOG")`


