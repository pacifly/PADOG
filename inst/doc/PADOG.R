### R code from vignette source 'PADOG.Rnw'

###################################################
### code chunk number 1: fig1
###################################################
library(KEGGdzPathwaysGEO)
library(PADOG)
library(Biobase)
set.seed(1)
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

#run padog function on KEGG pathways
#use NI=1000 for accurate results
myr=padog(
esetm=dat.m,
group=ano$Group,
paired=design=="Paired",
block=ano$Block,
targetgs=targetGeneSets,
annotation=annotation,
gslist="KEGG.db",
organism="hsa",
verbose=FALSE,
Nmin=3,
NI=50,
plots=TRUE)$targ

myr[1:15,-c(4,5)]


###################################################
### code chunk number 2: PADOG.Rnw:170-226
###################################################

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
return(res[res$ID %in% targetGeneSets,])
}

#run the analysis on all 24 datasets and compare the new method "myRand" with 
#PADOG and GSA (if installed) (chosen as reference since is listed first in the existingMethods)
#if the package parallel is installed datasets are analyzed in parallel.
#out=compPADOG(datasets=NULL,existingMethods=c("GSA","PADOG"),
 #mymethods=list(myRand=randomF),
 #gslist="KEGG.db",Nmin=3,NI=1000,plots=FALSE,verbose=FALSE,use.parallel=TRUE,dseed=1,pkgs=c("GSA","PADOG"))

#compare myRand against PADOG on 4 datasets only
#mysets=data(package="KEGGdzPathwaysGEO")$results[,"Item"]
mysets=c("GSE9348","GSE8671","GSE1297")
out=compPADOG(datasets=mysets,existingMethods=c("PADOG"),
 mymethods=list(myRand=randomF),
 gslist="KEGG.db",Nmin=3,NI=40,plots=TRUE,verbose=FALSE,use.parallel=FALSE,dseed=1,pkgs=NULL)

print(out)



