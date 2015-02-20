\name{padog}
\alias{padog}
\title{Pathway Analysis with Down-weighting of Overlapping Genes (PADOG)}
\description{
This is a general purpose gene set analysis method that downplays the importance of genes that apear often accross the sets of genes analyzed. 
The package provides also a benchmark for gene set analysis in terms of sensitivity and ranking using 24 public datasets.
}
\usage{
padog(esetm=NULL,group=NULL,paired=FALSE,block=NULL,gslist="KEGG.db",organism="hsa",annotation=NULL,gs.names=NULL,NI=1000,plots=FALSE,targetgs=NULL,Nmin=3,verbose=TRUE,paral=FALSE,ncr=NULL)
}
\arguments{
  \item{esetm}{A matrix containing log transfomed and normalized gene expression data. Rows correspond to genes and columns to samples.
  }
  \item{group}{A character vector with the class labels of the samples. It can only contain "c" for control samples or "d" for disease samples.}
  \item{paired}{A logical value to indicate if the samples in the two groups are paired.}
  \item{block}{A character vector indicating the block ids of the samples classified by the group variable, if  \code{paired=TRUE}.  The paired samples must have 
  the same block value.}
  \item{gslist}{Either the value "KEGG.db" or a list with the gene sets. If set to "KEGG.db", then gene sets will be made of all KEGG pathways for the \code{organism} 
  specified. If a list is provided, instead, each element of the list should be a character vector with the identifiers for the genes. The identifiers can be probe(sets) 
  ids if the \code{annotation} argument is set to a valid annotation package, otherwise the gene identifiers must be of the same kind as the rownames of the matrix esetm.}
  \item{annotation}{A valid chip annotation package if the rownames of \code{esetm} are probe(set) ids and \code{gslist} contains ENTREZ identifiers or \code{gslist} is set to "KEGG.db".
  If the rownames are other gene identifies, then \code{annotation} has tyo be set to NULL, and the row names of \code{esetm} needs to be unique and be found among elements of \code{gslist}}
  \item{organism}{A three letter string giving the name of the organism supported by the "KEGG.db" package.}
    \item{gs.names}{Character vector with the names of the gene sets. If specified, must have the same length as gslist.}
  \item{NI}{Number of iterations to determine the gene set score significance p-values.}
  \item{plots}{If set to TRUE then the distribution of the PADOG scores with and without weighting the genes in raw and standardized form are shown using boxplots. 
  A pdf file will be created in the current directory having the name provided in the \code{targetgs} field. The scores for the \code{targetgs} gene set will be shown
   in red. }
  \item{targetgs}{The identifier of a traget gene set for which the scores will be highlighted in the plots produced if \code{plots=TRUE}  }
  \item{Nmin}{The minimum size of gene sets to be included in the analysis.}
  \item{verbose}{If set to TRUE, displays the number of iterations elapsed is displayed.}
  \item{paral}{If set to TRUE, the \code{NI} iterations will be executed in parallel if multiple CPU cores are available.}
  \item{ncr}{The number of CPU cores used when \code{paral} set to TRUE. Default is to use all CPU cores detected.}
}

\details{
See cited documents for more details.
}


\value{
 A data frame containing the ranked pathways and various statistics: \code{Name} is the name of the gene set;
 \code{ID} is the gene set identifier; \code{Size} is the number of genes in the geneset; \code{meanAbsT0} is the mean of absolute t-scores;
  \code{padog0} is the mean of weighted absolute t-scores;
 \code{PmeanAbsT} significance of the meanAbsT0;   \code{Ppadog} is the significance of the padog0 score;
}

\references{
Adi L. Tarca, Sorin Draghici, Gaurav Bhatti, Roberto Romero, Down-weighting overlapping genes improves gene set analysis, BMC Bioinformatics, 2012, submitted.  \cr

}

\author{Adi Laurentiu Tarca <atarca@med.wayne.edu> \and Zhonghui Xu <zhonghui.xu@gmail.com>}

\seealso{\code{\link{padog}}}

\examples{

#run padog on a colorectal cancer dataset of the 24 datasets benchmark GSE9348
#use NI=1000 for accurate results.
set="GSE9348"
data(list=set,package="KEGGdzPathwaysGEO")
x=get(set)
#Extract from the dataset the required info
exp=experimentData(x);
dataset= exp@name
dat.m=exprs(x)
ano=pData(x)
design= notes(exp)$design
annotation= paste(x@annotation,".db",sep="")
targetGeneSets= notes(exp)$targetGeneSets


myr=padog(
esetm=dat.m,
group=ano$Group,
paired=design=="Paired",
block=ano$Block,
targetgs=targetGeneSets,
annotation=annotation,
gslist="KEGG.db",
organism="hsa",
verbose=TRUE,
Nmin=3,
NI=25,
plots=FALSE)


myr2=padog(
esetm=dat.m,
group=ano$Group,
paired=design=="Paired",
block=ano$Block,
targetgs=targetGeneSets,
annotation=annotation,
gslist="KEGG.db",
organism="hsa",
verbose=TRUE,
Nmin=3,
NI=25,
plots=FALSE,
paral=TRUE)


myr[1:20,]

all.equal(myr, myr2)


}

\keyword{nonparametric}% at least one, from doc/KEYWORDS
\keyword{methods}% __ONLY ONE__ keyword per line

