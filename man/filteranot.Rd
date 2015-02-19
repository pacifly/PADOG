\name{filteranot}
\alias{filteranot}
\title{Remove duplicate probesets/probes from an gene expression matrix based on p-values from a moderated t-test, in order to apply a gene set analysis.}
\description{
    This function helps to deal with multiple probesets/probes per gene prior to geneset analysis.
}
\usage{
filteranot(esetm=NULL,group=NULL,paired=FALSE,block=NULL,annotation=NULL,include.details=FALSE)
}
\arguments{
  \item{esetm}{A matrix containing log transfomed and normalized gene expression data. Rows correspond to genes and columns to samples.  Rownames of esetm need to be valid 
  probeset or probe names.
  }
  \item{group}{A character vector with the class labels of the samples. It can only contain "c" for control samples or "d" for disease samples.}
  \item{paired}{A logical value to indicate if the samples in the two groups are paired.}
  \item{block}{A character vector indicating the block ids of the samples classified by the group variable, if  \code{paired=TRUE}.  The paired samples must have 
  the same block value.}
   \item{annotation}{A valid chip annotation package name (e.g. "hgu133plus2.db")}
   \item{include.details}{If set to true, will include all columns from limma's topTable for this dataset.}
  }

\details{
See cited documents for more details.
}


\value{
 A data frame containing the probeset IDs (and corresponding ENTREZ IDs) of the best probesets for each gene ;
}

\references{
Adi L. Tarca, Sorin Draghici, Gaurav Bhatti, Roberto Romero, Down-weighting overlapping genes improves gene set analysis, BMC Bioinformatics, 2012, submitted.  \cr

}

\author{Adi Laurentiu Tarca <atarca@med.wayne.edu>}

\seealso{\code{\link{padog}}}

\examples{

#run padog on a colorectal cancer dataset of the 24 datasets benchmark GSE9348
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

dim(dat.m)
#get rid of duplicates in the same way as is done for PADOG and assign probesets to ENTREZ IDS
#get rid of duplicates by choosing the probe(set) with lowest p-value; get ENTREZIDs for probes
aT1=filteranot(esetm=dat.m,group=ano$Group,paired=(design=="Paired"),block=ano$Block,annotation) 

#filtered expression matrix
filtexpr=dat.m[rownames(dat.m)\%in\%aT1$ID,]
dim(filtexpr)


}

\keyword{nonparametric}% at least one, from doc/KEYWORDS
\keyword{methods}% __ONLY ONE__ keyword per line
