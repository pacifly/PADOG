% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/compFDR.R
\name{compFDR}
\alias{compFDR}
\title{Benchmark for gene set analysis methods using 24 datasets}
\usage{
compFDR(datasets = NULL, existingMethods = c("GSA", "PADOG"),
  mymethods = NULL, gs.names = NULL, gslist = "KEGG.db",
  organism = "hsa", Nmin = 3, NI = 1000, use.parallel = TRUE,
  ncr = NULL, pkgs = NULL, expVars = NULL, dseed = NULL, Npsudo = 20,
  FDRmeth = c("BH", "Permutation", "holm"), plots = FALSE,
  verbose = FALSE)
}
\arguments{
\item{datasets}{A character vector with valid names of datasets to use from the PADOGsets
package. If left NULL all datasets avalibale in PADOGsets will be used.}

\item{existingMethods}{A character vector with one or more of the predefined methods
c("GSA","PADOG"). The first is used as reference method.}

\item{mymethods}{A list whose elements are valid functions implementing gene set analysis
methods. See \code{Details} and \code{Examples} for what arguments the functions have to
take in and what kind of output they need to return.}

\item{gs.names}{A character vector giving additional information about each gene set. For
instance when gene seta are pathways, the full name of the pathway would be a meaningful
gene set name.}

\item{gslist}{Either the value "KEGG.db" or a list with the gene sets. If set to "KEGG.db",
then gene sets will be made of all KEGG pathways for human since all datasets available
in PADOG are for human.}

\item{organism}{A three letter string giving the name of the organism supported by the
"KEGG.db" package.}

\item{Nmin}{The minimum size of gene sets to be included in the analysis for all methods.}

\item{NI}{Number of iterations to determine the gene set score significance p-values in
PADOG and GSA methods.}

\item{use.parallel}{If set to TRUE and multiple CPU cores are available, the parallelization
will be distributed either over the \code{NI} iterations for "PADOG", or over
datasets/methods combination otherwise.}

\item{ncr}{The number of CPU cores used when \code{use.parallel} set to TRUE. Default is
to use all CPU cores detected.}

\item{pkgs}{Character vector of packages that the \code{existingMethods} and \code{mymethods}
depend on (\code{NULL} for "PADOG" and "GSA"). Consult the \code{.packages}
argument in \code{foreach} function from \code{foreach} package.}

\item{expVars}{Character vector of variables to export. Consult the \code{.export} argument
in \code{foreach} function from \code{foreach} package.}

\item{dseed}{Optional initial seed for random number generator (integer) used in \code{padog}.}

\item{Npsudo}{The number of permutations on phenotype (i.e. group label) to obtain null p
values and estimate fdr. Set to 0 if not interested in estimating fdr from permutation.}

\item{FDRmeth}{The method used to estimate fdr; "BH" for Benjamini & Hochberg, "Permutation"
for estimation based on observed and null p values. You can use "holm" for FWER as well.}

\item{plots}{If set to TRUE will plot the p values, ranks, the ranks differences w.r.t.
reference, null p values and fdr when applicable.}

\item{verbose}{If set to TRUE, for "PADOG" and "AbsmT" methods it will show the iterations
performed so far, and for others show the method and data set currently being executed
(when \code{use.parallel = TRUE}, this only works for command line R, not Rgui).}
}
\value{
A list of elements "summary","ranks","pvalues","qvalues" and "nullp" (if \code{Npsudo > 0}).
  "summary" is a data frame containing: \code{Method} is the name of the gene set analysis method;
    \code{p geomean} geometric mean of nominal p-values for the target gene sets (gene sets expected
    to be relevant); \code{p med} median of nominal p-values for the target genesets; \code{ \% p<0.05}
    is the fraction of target gene sets significant at 0.05 level (this is the sensitivity);
    \code{ \% q<0.05} is the fraction of target gene sets significant at 0.05 level after FDR correction;
    \code{rank mean} mean rank of the target gene sets; \code{rank med} median rank of the target gene sets;
    \code{p Wilcox.} p value from a Wilcoxon test paired at dataset level comparing the rank of target gene sets;
    \code{p LME} p value from a linear mixed effects (LME) model which unlike the Wilcoxon test above
      accounts for the fact that ranks for the same gene set may be correlated;
    \code{coef LME} Coefficient from the LME model giving the difference in ranks of the target gene
      sets between the current gene set analysis method and the reference method (the first method
      in the \code{existingMethods} argument.
  "ranks", "pvalues", "qvalues": a list of elements containing for each method the ranks (
    nominal p values, fdr q values) of the target gene sets in the corresponding data sets.
  "nullp": a list of elements containing for each method the null p values (a list of two
    elements "ap" for all gene sets and "tp" for target gene sets) pooled over all data sets.
}
\description{
This is a general purpose function to compare a given gene set analysis method in terms of
sensitivity, ranking and fdr against PADOG and GSA (if installed) using 24 public datasets.
}
\details{
The user defined gene set analysis method should return a list with elements "targ" containing
the analysis statistics for the target gene sets on a given data set, and "pval" containing
a matrix of null p values each column of which corresponds to the p values of all gene sets
applying the gene set analysis method on permutated phenotype (i.e. group label). The rownames
of the matrix are the gene set IDs; the number of column corresponds to argument \code{Npsudo}.
If \code{Npsudo = 0}, the "pval" element can be omitted and null p values / fdr will not be
compared between the methods. For "PADOG" and "AbsmT", two methods ("BH" and "Permutation") of
estimating fdr are available regardless of the value of \code{Npsudo} since null p values can be
accessed from the method algorithm itself. For "GSA", the "Permutation" method is only available
when \code{Npsudo > 0}.
}
\examples{
setF = function(type=c("pos","uni","neg"), NI=1000, dseed=1, FDRmeth=c("BH","Permutation","holm")) {
    ### use this caller function to return some artificial gene set analysis methods:
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
        mygslist = lapply(mygslist, function(z) z[z \%in\% (aT1$ENTREZID)])
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
        return(list(targ = res[res$ID \%in\% targetGeneSets,], pval = nullp))
      }
    )
}

mymnam = c("pos","uni","neg")
mym = lapply(mymnam, setF, NI=20, dseed=1, FDRmeth="P")
names(mym) = mymnam

# run the analysis on 3 datasets and compare the new methods with PADOG
# mysets = data(package="PADOGsets")$results[,"Item"]
mysets = c("GSE9348","GSE8671","GSE1297")
out = compFDR(datasets=mysets, existingMethods=c("AbsmT", "PADOG"),
    mymethods=mym, gslist="KEGG.db", Nmin=3, NI=20, Npsudo=0, FDRmeth="P", 
    plots=FALSE, verbose=FALSE, use.parallel=FALSE, dseed=1, pkgs=NULL)

}
\author{
Adi Laurentiu Tarca \email{atarca@med.wayne.edu},
   Zhonghui Xu \email{zhonghui.xu@gmail.com}
}
\references{
Adi L. Tarca, Sorin Draghici, Gaurav Bhatti, Roberto Romero. Down-weighting
   overlapping genes improves gene set analysis. BMC Bioinformatics, 2012.
}
\seealso{
See \code{\link{compPADOG}} for back-compatibility (i.e. no specificity/fdr comparison).
   Refer to latest development at \url{https://github.com/pacifly/PADOG}
}
\keyword{methods}
\keyword{nonparametric}

