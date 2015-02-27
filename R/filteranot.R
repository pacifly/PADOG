
#' @export
#'
#' @importFrom limma lmFit makeContrasts contrasts.fit eBayes topTable
#' @importFrom AnnotationDbi get as.list keys
#'
filteranot=function(esetm=NULL,group=NULL,paired=FALSE,block=NULL,annotation=NULL,include.details=FALSE){

Block=block
stopifnot(require(annotation,character.only=TRUE))

if(!paired){
G=factor(group)
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
G=factor(G)
design <- model.matrix(~0+G+block)
colnames(design)<-substr(colnames(design),2,100)
fit <- lmFit(esetm, design)
cont.matrix <- makeContrasts(contrasts="d-c",levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
aT1<-topTable(fit2,coef=1, number=dim(esetm)[1])
}
#annotate ids to ENTREZID
aT1$ID=rownames(aT1) 
anpack=paste(unlist(strsplit(annotation,split=".db")),"ENTREZID",sep="")
aT1<-aT1[aT1$ID%in%keys(get(anpack)),]
aT1$ENTREZID=unlist(as.list(get(anpack)[aT1$ID]))
aT1<-aT1[!is.na(aT1$ENTREZID),]
aT1<-aT1[!duplicated(aT1$ENTREZID),]

if(include.details){
 return(aT1)
}else{
return(aT1[,c("ID","ENTREZID")])
}

}
