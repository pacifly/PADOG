

compPADOG=function(datasets=NULL,existingMethods=c("GSA","PADOG"),mymethods=NULL,gs.names=NULL,gslist="KEGG.db",organism="hsa",Nmin=3,NI=1000,use.parallel=TRUE,plots=FALSE,verbose=FALSE){

if(is.null(datasets)){
files=data(package="KEGGdzPathwaysGEO")$results[,"Item"]
}else{
 files=datasets
}


getdataaslist= function(x)
{
data(list=x,package="KEGGdzPathwaysGEO")
x=get(x)
exp=experimentData(x);dataset= exp@name;disease= notes(exp)$disease
dat.m= exprs(x);ano=pData(x);design= notes(exp)$design;annotation= paste(x@annotation,".db",sep="")
targetGeneSets= notes(exp)$targetGeneSets
list=list(dataset,disease,dat.m,ano,design,annotation,targetGeneSets)
names(list)= c("dataset","disease","dat.m","ano","design","annotation","targetGeneSets")
return(list)
}



absmtF=function(set,mygslist,minsize){

list=getdataaslist(set)

set.seed(1)
res=padog(
esetm=list$dat.m,
group=list$ano$Group,
paired=list$design=="Paired",
block=list$ano$Block,
annotation=list$annotation,
gslist=mygslist,
verbose=verbose,
Nmin=minsize,
NI=NI,
plots=FALSE)
res$Dataset<-list$dataset
res=res[order(res$PmeanAbsT,-res$meanAbsT0),]
res$Method<-"AbsmT"
res$Rank=(1:dim(res)[1])/dim(res)[1]*100
res$P=res$PmeanAbsT
res$FDR=p.adjust(res$P,"fdr")
res[res$ID%in%list$targetGeneSets,]
}

padogF=function(set,mygslist,minsize){
list=getdataaslist(set)

set.seed(1)
res=padog(
esetm=list$dat.m,
group=list$ano$Group,
paired=list$design=="Paired",
block=list$ano$Block,
annotation=list$annotation,
gslist=mygslist,
verbose=verbose,
Nmin=minsize,
NI=NI,
plots=FALSE)
res$Dataset<-list$dataset
res$Method<-"PADOG"
res$Rank=(1:dim(res)[1])/dim(res)[1]*100
res$P=res$Ppadog
res$FDR=p.adjust(res$P,"fdr")
rownames(res)<-NULL
res[res$ID%in%list$targetGeneSets,]
}


gsaF=function(set,mygslist,minsize){
require(GSA)
list=getdataaslist(set)

group=list$ano$Group
block=list$ano$Block
esetm=list$dat.m
if(!is.null(list$annotation)){
#get rid of duplicates in the same way as is done for PADOG and assign probesets to ENTREZ IDS
 #get rid of duplicates by choosing the probe(set) with lowest p-value; get ENTREZIDs for probes
aT1=filteranot(esetm=esetm,group=group,paired=(list$design=="Paired"),block,list$annotation) 
aT1<-aT1[aT1$ENTREZID%in%(unlist(mygslist)),]
#drop from esetm all duplicate genes and genes not in the genesets
esetm=esetm[rownames(esetm)%in%aT1$ID,]
rownames(esetm)<-aT1$ENTREZID[match(rownames(esetm),aT1$ID)]
}#
#Run GSA maxmean
nc=table(list$ano$Group)["c"]
nd=table(list$ano$Group)["d"]
if(list$design=="Not Paired"){yy=c(rep(1,nc),rep(2,nd))}else {
block=as.numeric(factor(list$ano$Block))
block[duplicated(block)]<-(-block[duplicated(block)])
yy=block
}

resgsa=GSA(x=esetm,y=yy, genesets=mygslist, genenames=rownames(esetm),
method="maxmean", resp.type=ifelse(list$design=="Not Paired","Two class unpaired","Two class paired"),
censoring.status=NULL,random.seed=1,  knn.neighbors=10,
s0=NULL, s0.perc=NULL,minsize=minsize,maxsize=1000,
restand=TRUE,restand.basis=c("catalog","data"),
 nperms=NI)

res=data.frame(ID=names(gslist),P=2*apply(cbind(resgsa$pvalues.lo,resgsa$pvalues.hi),1,min),
Dataset=list$dataset,stringsAsFactors=FALSE)
res$Method<-"GSA"
res=res[order(res$P),]
res$Rank=rank(res$P)/dim(res)[1]*100
res$FDR=p.adjust(res$P,"fdr")
rownames(res)<-NULL
res[res$ID%in%list$targetGeneSets,]
}


if(require(GSA)){
defGSmethods=list(GSA=gsaF,PADOG=padogF,AbsmT=absmtF)
}else{
defGSmethods=list(ABSMT=absmtF,PADOG=padogF)
}

GSmethods=defGSmethods[intersect(names(defGSmethods),existingMethods)]

GSmethods=c(GSmethods,mymethods)

refMethod=names(GSmethods)[1]



#check GS
if(length(gslist)==1 && gslist=="KEGG.db"){
require(KEGG.db)
pathsids=names(as.list(KEGGPATHID2EXTID))[grep(organism,names(as.list(KEGGPATHID2EXTID)))]
gslist=as.list(KEGGPATHID2EXTID)[pathsids]
names(gslist)=substr(names(gslist),4,nchar(names(gslist)))
gs.names=unlist(as.list(KEGGPATHID2NAME)[names(gslist)])
stopifnot(length(gslist)>=3)
}
stopifnot(class(gslist)=="list")
stopifnot(length(gslist)>=3)
if(!is.null(gs.names)){stopifnot(length(gslist)==length(gs.names))}







reslist=list()
for(i in 1:length(GSmethods)){
if(require(parallel)&use.parallel){
reslist[[names(GSmethods)[i]]]<- mclapply(files,GSmethods[[i]],mygslist=gslist,minsize=Nmin)
}else{
 reslist[[names(GSmethods)[i]]]<- lapply(files,GSmethods[[i]],mygslist=gslist,minsize=Nmin)
}
}



dfs=list();
for(i in 1:length(reslist))
{y=NULL
  for(j in 1:length(reslist[[i]])){
  tm=as.data.frame(reslist[[i]][[j]])
  rownames(tm)<-NULL
 y=rbind(y,tm[,c("ID","Rank","P","FDR","Dataset","Method")])}
 dfs[[names(reslist)[i]]]<-y
}


#identify common genesets/datasets combinations among the 3 methods.
nmsp=lapply(dfs,function(x){paste(x$Dataset,x$ID,sep="_")})

for( i in 2:length(nmsp)){
 stopifnot(all(nmsp[[i]]==nmsp[[i-1]]))
}

psList<-lapply(dfs,function(x){x$P})
fdrList<-lapply(dfs,function(x){x$FDR})
rankList<-lapply(dfs,function(x){x$Rank})
targetgsList<-lapply(dfs,function(x){x$ID})
dsList<-lapply(dfs,function(x){x$Dataset})



wi=function(x){
if(!all(x==rankList[[refMethod]])){
wilcox.test(x,rankList[[refMethod]],paired=TRUE,alternative="less")$"p.value"
}else{1}
}

require(nlme)
wioright=function(x){
if(!all(x==rankList[[refMethod]])){
dset=data.frame(Method=rep(c(1,0),each=length(rankList[[refMethod]])),
Y=c(x,rankList[[refMethod]]),Dataset=factor(dsList[[refMethod]]),
Path=factor(c(targetgsList[[refMethod]],targetgsList[[refMethod]])))
md=lme(Y~Method+Dataset,random = ~1|Path,data=dset)
re=summary(md)$tTable["Method",c(1,5)]
if(re[1]<0){c(re[1],re[2]/2)}else{c(re[1],1-re[2]/2)}
}else{c(0,1)}
}




repo=data.frame(matrix(unlist(lapply(rankList,wioright)),length(psList),2,byrow=TRUE))
names(repo)<-c("coef. LME","p LME")
repo$"p Wilcox."<-lapply(rankList,wi)

l05<-function(x){round(sum(x<0.05)/length(x)*100,2)}
geomean<-function(x){prod(x)^(1/length(x))}

repo$"% p.value<0.05"<-lapply(psList,l05)
repo$"% q.value<0.05"<-lapply(fdrList,l05)
repo$"p geomean"<-lapply(psList,geomean)
repo$"p med"<-lapply(psList,median)

repo$"rank mean"<-lapply(rankList,mean)
repo$"rank med"<-lapply(rankList,median)
repo$Method<-names(psList)

nmets=length(psList)

somecols=c("lightgrey","lightblue","orange","red","blue","grey")
set.seed(1)
if(nmets>6){somecols=c(somecols,sample(colors())[1:(nmets-6)])}

if(plots){

par(mfrow=c(1,3))
boxplot(psList,ylab=paste("p-value"),las = 3,col=somecols[1:nmets])
boxplot(rankList,ylab="Rank(%)",las = 3,col=somecols[1:nmets])

mff2=function(x){x-rankList[[refMethod]]}
newranks=lapply(rankList,mff2)
newranks[refMethod]<-NULL

if(length(newranks)==1){xlb=names(newranks)}else{xlb=NULL};
boxplot(newranks,ylab=paste("Rank(%)-Rank ",refMethod," (%)") ,las = 3,col=somecols[2:nmets],names=names(newranks),xlab=xlb)
abline(h=0)

}

out=repo[,c("Method","p geomean","p med","% p.value<0.05","% q.value<0.05","rank mean","rank med","p Wilcox.","p LME","coef. LME")]
list(summary=out,ranks=rankList,pvalues=psList,qvalues=fdrList)

}




