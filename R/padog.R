
padog<-function(esetm=NULL,group=NULL,paired=FALSE,block=NULL,gslist="KEGG.db",organism="hsa",annotation=NULL,gs.names=NULL,NI=1000,plots=FALSE,targetgs=NULL,Nmin=3,verbose=TRUE, paral = FALSE, ncr = NULL){

#load needed packages
    require(limma)

#initialize the gslist if using KEGG
    if (length(gslist) == 1 && gslist == "KEGG.db") {
        stopifnot(nchar(organism) == 3)
        require(KEGG.db)
        pw2id = as.list(KEGGPATHID2EXTID) 
        gslist = pw2id[grep(organism, names(pw2id))]
        names(gslist) = sub(paste("^",organism,sep=""), "", names(gslist))
        gs.names = unlist(as.list(KEGGPATHID2NAME)[names(gslist)])
        stopifnot(length(gslist) >= 3)
        rm(pw2id)
    }

#check arguments
stopifnot(class(esetm)=="matrix")
stopifnot(all(dim(esetm)>4))

stopifnot(class(group)%in%c("factor","character"))
stopifnot(length(group)==dim(esetm)[2])
stopifnot(all(group%in%c("c","d")))
stopifnot(all(table(group)>2))
if(paired){stopifnot(length(block)==length(group))
stopifnot(all(table(block)==2))
}

stopifnot(class(gslist)=="list")
stopifnot(length(gslist)>=3)
if(!is.null(gs.names)){stopifnot(length(gslist)==length(gs.names))}

stopifnot(class(NI)=="numeric")
stopifnot(NI>5)

if(plots){
stopifnot(targetgs%in%names(gslist))
}
if(!is.null(annotation)){
stopifnot(require(annotation,character.only=TRUE))
stopifnot( sum(rownames(esetm)%in%mappedkeys(get(paste(substr(annotation,1,nchar(annotation)-3),"ENTREZID",sep=""))))>4)
}else{
 stopifnot(sum(rownames(esetm)%in%as.character(unlist(gslist)))>10 & !any(duplicated(rownames(esetm))))
}



#substitute some names
Block=block




#


#compute gene frequencies accross genesets
gf=table(unlist(gslist))
if(!all(gf==1)){
if(quantile(gf,0.99)>mean(gf)+3*sd(gf)){
gf[gf>quantile(gf,0.99)]<-quantile(gf,0.99)
}
gff<-function(x){1+((max(x)-x)/(max(x)-min(x)))^0.5}
#compute weights
gf=gff(gf)
}else{

 fdfd=unique(unlist(gslist))
 gf=rep(1,length(fdfd))
 names(gf)<-fdfd
}


allGallP=unique(unlist(gslist))



if(!is.null(annotation)){
#get rid of duplicates in the esetm by choosing the probe(set) with lowest p-value; get ENTREZIDs for probes
aT1=filteranot(esetm,group,paired,block,annotation)
#drop genes not in any geneset
#drop from esetm all duplicate genes and genes not in the genesets
esetm=esetm[rownames(esetm)%in%aT1$ID,]
rownames(esetm)<-aT1$ENTREZID[match(rownames(esetm),aT1$ID)]
}
restg=setdiff(rownames(esetm),names(gf))
appendd=rep(1,length(restg))
names(appendd)<-restg
gf=c(gf,appendd)



 stopifnot(all(!duplicated(rownames(esetm))))
 stopifnot(sum(rownames(esetm)%in%allGallP)>10)

if(verbose){
cat(paste("Starting with ",length(gslist)," gene sets!",sep=""));
cat("\n");
}

#drop pathways with less than Nmin genes
gslist=gslist[unlist(lapply(gslist,function(x){length(intersect(rownames(esetm),x))>=Nmin}))]
gs.names=gs.names[names(gslist)]
stopifnot(length(gslist)>=3)

if(verbose){
cat(paste("Analyzing ",length(gslist)," gene sets with ",Nmin, " or more genes!",sep=""));
cat("\n");
}

##############################################
#compute scores for iterations 
MSabsT<-MSTop<-matrix(NA,length(gslist),NI+1)

for(ite in 1:(NI+1)){
 SabsT<-STop<-NULL


if(!paired){
G=factor(group)
if(ite>1){
G=factor(sample(group))
}
design <- model.matrix(~0+G)
colnames(design) <- levels(G)
fit <- lmFit(esetm, design)
cont.matrix <- makeContrasts(contrasts="d-c",levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
aT1<-topTable(fit2,coef=1, number=dim(esetm)[1])
}else{
G=group
block=factor(Block)
if(ite>1){
for(ss in 1:length(levels(block))){
 G[block==levels(block)[ss]]<-sample(G[block==levels(block)[ss]])
}
}
G=factor(G)
design <- model.matrix(~0+G+block)
colnames(design)<-substr(colnames(design),2,100)
fit <- lmFit(esetm, design)
cont.matrix <- makeContrasts(contrasts="d-c",levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
aT1<-topTable(fit2,coef=1, number=dim(esetm)[1])
}
aT1$ID=rownames(aT1)

de=abs(aT1$t)
names(de)<-aT1$ID


meande=mean(de)
sdde=sd(de)
meangfde=mean(de*gf[names(de)])
sdgfde=sd(de*gf[names(de)])

#for each pathway
psizes=NULL
 for(i in 1:length(names(gslist))){
  path<-names(gslist)[i]
  X<-(de[gslist[[path]]])
  X=na.omit(X)

 SabsT[i]<-mean(X)
 STop[i]<-mean((X*gf[names(X)]))
 psizes[i]=length(X)
 SabsT[i]<-(SabsT[i]-meande)/(sdde/sqrt(psizes[i]))
 STop[i]<-(STop[i]-meangfde)/(sdgfde/sqrt(psizes[i]))
 }
 
 

MSabsT[,ite]<-SabsT
MSTop[,ite]<-STop
if(ite==1){
 meanAbsT0=SabsT
 padog0=STop
}

if(verbose & ite/10==round(ite/10,0)){cat(paste(ite,"/",NI));cat("\n")}
}

plotIte=min(NI,21)

MSabsT_raw=MSabsT
MSTop_raw=MSTop

#standardize scores
for(k in 1:(NI+1)){
 MSabsT[,k]<-(MSabsT[,k]-mean(MSabsT[,k],na.rm=TRUE))/sd(MSabsT[,k],na.rm=TRUE)
 MSTop[,k]<-(MSTop[,k]-mean(MSTop[,k],na.rm=TRUE))/sd(MSTop[,k],na.rm=TRUE)
}

#compute p-values 
mff=function(x){if(!all(is.na(x))){sum(x[-1]>x[1],na.rm=TRUE)/sum(!is.na(x[-1]))}else{NA}}
PSabsT=apply(MSabsT,1,mff)
PSTop=apply(MSTop,1,mff)

PSabsT[PSabsT==0]<-1/NI/100
PSTop[PSTop==0]<-1/NI/100

#do plot the scores for the star pathway
if(plots){


par(mfrow=c(2,2))
boxplot(MSabsT_raw[,1:plotIte]~col(MSabsT_raw[,1:plotIte]),col=c("lightblue",rep("whitesmoke",NI)),names=c("0",1:(plotIte-1)),cex.axis=0.8,main="ABSmT scores after first standardization",cex.main=0.6)
points(1:plotIte,MSabsT_raw[names(gslist)==targetgs,1:plotIte],col="red",pch=19)
abline(h=MSabsT_raw[names(gslist)==targetgs,1],col="red")
boxplot(MSTop_raw[,1:plotIte]~col(MSTop_raw[,1:plotIte]),col=c("lightblue",rep("whitesmoke",NI)),names=c("0",1:(plotIte-1)),cex.axis=0.8,main="PADOG after first standardization",cex.main=0.6)
points(1:plotIte,MSTop_raw[names(gslist)==targetgs,1:plotIte],col="red",pch=19)
abline(h=MSTop_raw[names(gslist)==targetgs,1],col="red")

boxplot(MSabsT[,1:plotIte]~col(MSabsT[,1:plotIte]),col=c("lightblue",rep("whitesmoke",NI)),names=c("0",1:(plotIte-1)),cex.axis=0.8,main="ABSmT scores second standardization",cex.main=0.6)
points(1:plotIte,MSabsT[names(gslist)==targetgs,1:plotIte],col="red",pch=19)
abline(h=MSabsT[names(gslist)==targetgs,1],col="red")
boxplot(MSTop[,1:plotIte]~col(MSTop[,1:plotIte]),col=c("lightblue",rep("whitesmoke",NI)),names=c("0",1:(plotIte-1)),cex.axis=0.8,main="PADOG after second standardization",cex.main=0.6)
points(1:plotIte,MSTop[names(gslist)==targetgs,1:plotIte],col="red",pch=19)
abline(h=MSTop[names(gslist)==targetgs,1],col="red")

}

if(!is.null(gs.names)){myn=gs.names}else{myn=names(gslist)} 

SIZE=unlist(lapply(gslist,function(x){length(intersect(rownames(esetm),x))}))           
res=data.frame(Name=myn,ID=names(gslist),Size=SIZE,meanAbsT0,padog0,PmeanAbsT=PSabsT,Ppadog=PSTop,stringsAsFactors=FALSE)

res=res[order(res$Ppadog,-res$padog0),]
res

}


