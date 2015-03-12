setF = function(type=c("pos","uni","neg"), NI=1000, dseed=1, FDRmeth=c("BH","Permutation","holm")) {
    ### use this caller function to return some artificial gene set analysis methods:
    ### type -- uni for uniform p values, pos / neg for excessive false positive / negative respectively
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
mym = lapply(mymnam, setF, NI=20, dseed=1, FDRmeth="P")
names(mym) = mymnam

# run the analysis on 3 datasets and compare the new methods with PADOG
# mysets = data(package="PADOGsets")$results[,"Item"]
mysets = c("GSE9348","GSE8671","GSE1297")
out = compFDR(datasets=mysets, existingMethods=c("AbsmT", "PADOG"),
    mymethods=mym, gslist="KEGG.db", Nmin=3, NI=20, Npsudo=0, FDRmeth="P", 
    plots=FALSE, verbose=FALSE, use.parallel=FALSE, dseed=1, pkgs=NULL)

