\name{compPADOG}
\alias{compPADOG}
\title{Benchmark for gene set analysis methods using 24 datasets}
\description{
This is a general purpose function to compare a given gene set analysis method in terms of sensitivity and ranking agains PADOG and GSA (if installed) using 24 
public datasets.
}
\usage{
compPADOG(datasets=NULL,existingMethods=c("GSA","PADOG"),mymethods=NULL,gs.names=NULL,
          gslist="KEGG.db",organism="hsa",Nmin=3,NI=1000,use.parallel=TRUE,ncr=NULL,
          pkgs="GSA", expVars=NULL, dseed=NULL,plots=FALSE,verbose=FALSE)
}
\arguments{
  \item{datasets}{A character vector with valid names of datasets to use from the PADOGsets package. If left NULL all datasets avalibale in PADOGsets will be used.
  }
  \item{existingMethods}{A character vector with one or more of the predefined methods c("GSA","PADOG"). The first is used as reference method.}
  \item{mymethods}{A list whose elements are valid functions implementing gene set analysis methods. See the example to see what arguments the functions have to take in and 
  what kind of output they need to produce.}
  \item{gs.names}{A character vector giving additional information about each gene set. For instance when gene seta are pathways, the full name of the pathway would be a meaningful gene set name.}
  \item{gslist}{Either the value "KEGG.db" or a list with the gene sets. If set to "KEGG.db", then gene sets will be made of all KEGG pathways for human since all datasets available in 
  PADOG are for human.}
   \item{organism}{A three letter string giving the name of the organism supported by the "KEGG.db" package.}
  \item{Nmin}{The minimum size of gene sets to be included in the analysis for all methods.}
  \item{NI}{Number of iterations to determine the gene set score significance p-values in PADOG and GSA methods.}
  \item{use.parallel}{If set to TRUE and multiple CPU cores are available, the parallelization will be distributed either over the \code{NI} iterations for "PADOG", or over datasets/methods combination otherwise.}
  \item{ncr}{The number of CPU cores used when \code{use.parallel} set to TRUE. Default is to use all CPU cores detected.}
  \item{pkgs}{Character vector of packages that the \code{existingMethods} and \code{mymethods} depend on (e.g. \code{NULL} for "PADOG", "GSA" for "GSA"). Consult the \code{.packages} argument in \code{foreach} function from \code{foreach} package.}
  \item{expVars}{Character vector of variables to export. Consult the \code{.export} argument in \code{foreach} function from \code{foreach} package.}
  \item{dseed}{Optional initial seed for random number generator (integer) used in \code{padog}.}
  \item{plots}{If set to TRUE will plot the ranks of the target genesets and the ranks differences between a methods and the reference method.}
  \item{verbose}{This argument will be passed to PADOG and AbsmT methods. If set to TRUE
  it will show the interations performed so far.}
}

\details{
See cited documents for more details.
}


\value{
 A data frame containing the : \code{Method} is the name of the geneset analysis method;
 \code{p geomean} geometric mean of nominal p-values for the target genesets (genesets expected to be relevant);
  \code{p med} median of nominal p-values for the target genesets;
  \code{ \% p<0.05} is the fraction of target genesets significant at 0.05 level (this is the sensitivity);
  \code{ \% q<0.05} is the fraction of target genesets significant at 0.05 level after FDR correction;
  \code{rank mean} mean rank of the target genesets;
  \code{rank med} median rank of the target genesets;
  \code{p Wilcox.} p value from a Wilcoxon test paired at dataset level comparing the rank of target genesets ;
  \code{p LME} p value from a linear mixed effects (LME) model which unlike the Wilcoxon test above accounts for the fact that ranks for the same pathway may be correlated;
  \code{coef LME} Coefficient from the LME model giving the difference in ranks of the target genesets between the current geneset analysis Method and the reference method
  chose to be the first method in the \code{existingMethods} argument;
  
}

\references{
Adi L. Tarca, Sorin Draghici, Gaurav Bhatti, Roberto Romero, Down-weighting overlapping genes improves gene set analysis, BMC Bioinformatics, 2012, submitted.  \cr

}

\author{Adi Laurentiu Tarca <atarca@med.wayne.edu> and Zhonghui Xu <zhonghui.xu@gmail.com>}

\seealso{\code{\link{compPADOG}}}

\examples{

#compare a new geneset analysis method with PADOG and GSA

#define your new gene set analysis method that takes as input:
#set- the name of dataset file from the PADOGsetspackage
#mygslist - a list with the genesets
#minsize- minimum number of genes in a geneset to be considered for analysis 

randomF=function(set,mygslist,minsize){
set.seed(1)
#this loads the dataset in an ExpressionSet object called x
data(list=set,package="KEGGdzPathwaysGEO")
x=get(set)

#Extract from the dataset the required info to be passed to padog
exp=experimentData(x);
dat.m=exprs(x)
ano=pData(x)
dataset= exp@name
design= notes(exp)$design
annotation= paste(x@annotation,".db",sep="")
targetGeneSets= notes(exp)$targetGeneSets


#get rid of duplicates probesets per ENTREZ ID by keeping the probeset 
#with smallest p-value (computed using limma) 
aT1=filteranot(esetm=dat.m,group=ano$Group,paired=(design=="Paired"),
 block=ano$Block,annotation=annotation)
#create an output dataframe for this toy method with random gene set p-values
mygslistSize=unlist(lapply(mygslist,function(x){length(intersect(aT1$ENTREZID,x))}))
res=data.frame(ID=names(mygslist),P=runif(length(mygslist)),
 Size=mygslistSize,stringsAsFactors=FALSE)
res$FDR=p.adjust(res$P,"fdr")
#drop genesets with less than minsize genes in the current dataset 
res=res[res$Size>=minsize,]
#compute ranks
res$Rank=rank(res$P)/dim(res)[1]*100
#needed to compare ranks between methods; must be the same as given 
#in mymethods argument "list(myRand="
res$Method="myRand";
#needed because comparisons of ranks between methods is paired at dataset level
res$Dataset<-dataset;
#output only result for the targetGeneSets 
#which are gene sets expected to be relevant in this dataset
return(res[res$ID \%in\% targetGeneSets,])
}

#run the analysis on all 24 datasets and compare the new method "myRand" with 
#PADOG and GSA (if installed) (chosen as reference since is listed first in the existingMethods)
#out=compPADOG(datasets=NULL,existingMethods=c("GSA","PADOG"),
 #mymethods=list(myRand=randomF),gslist="KEGG.db",Nmin=3,NI=1000,
 #plots=FALSE,verbose=FALSE,use.parallel=TRUE,dseed=1,pkgs=c("GSA","PADOG"))

#compare myRand against PADOG on 3 datasets only
#mysets=data(package="PADOGsets")$results[,"Item"]
mysets=c("GSE9348","GSE8671","GSE1297")
out=compPADOG(datasets=mysets,existingMethods=c("PADOG"),
  mymethods=list(myRand=randomF),gslist="KEGG.db",Nmin=3,NI=20,
  plots=FALSE,verbose=FALSE,use.parallel=FALSE,dseed=1,pkgs=NULL)

}

\keyword{nonparametric}% at least one, from doc/KEYWORDS
\keyword{methods}% __ONLY ONE__ keyword per line
