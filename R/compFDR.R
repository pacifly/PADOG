
if(getRversion() >= "2.15.1")  utils::globalVariables(c("ini", "outi"))

#' @export
#'
#' @importFrom AnnotationDbi get as.list 
#' @importFrom Biobase experimentData notes exprs pData 
#' @importFrom GSA GSA
#' @importFrom foreach foreach %dopar% %:%
#'
compFDR = function(datasets = NULL, existingMethods = c("GSA", "PADOG"), mymethods = NULL, 
    gs.names = NULL, gslist = "KEGG.db", organism = "hsa", Nmin = 3, NI = 1000, use.parallel = TRUE, 
    ncr = NULL, pkgs = "GSA", expVars = NULL, dseed = NULL, Npsudo = 20, plots = FALSE, verbose = FALSE) {
    
    if (is.null(datasets)) {
        files = data(package = "KEGGdzPathwaysGEO")$results[, "Item"]
    } else {
        files = datasets
    }
    
    data(list = files, package = "KEGGdzPathwaysGEO", envir=environment())
    
    padogF = function(set, mygslist, minsize) {
        list = getdataaslist(set)
        
        out = padog(esetm = list$dat.m, group = list$ano$Group, paired = list$design == 
            "Paired", block = list$ano$Block, annotation = list$annotation, gslist = mygslist, 
            verbose = verbose, Nmin = minsize, NI = NI, plots = FALSE, paral = use.parallel, 
            ncr = ncr, dseed = dseed)

        res = out$res
        res$Dataset = list$dataset
        
        res$Method = "PADOG"
        res$Rank = (1:nrow(res))/nrow(res) * 100
        res$P = res$Ppadog
        # res$FDR = p.adjust(res$P, "fdr")
        res$FDR = res$FDRpadog
        
        rownames(res) = NULL
        pidx = res$ID %in% list$targetGeneSets
        list(targ = rbind(res[pidx, ], padog2absmt(res, list, estFDR=TRUE)), pval = out$pval)
    }
    
    
    gsaF = function(set, mygslist, minsize) {
        list = getdataaslist(set)
        
        group = list$ano$Group
        block = list$ano$Block
        esetm = list$dat.m
        if (!is.null(list$annotation)) 
            {
                # get rid of duplicates in the same way as is done for PADOG and assign probesets
                # to ENTREZ IDS get rid of duplicates by choosing the probe(set) with lowest
                # p-value; get ENTREZIDs for probes
                aT1 = filteranot(esetm = esetm, group = group, paired = (list$design == 
                  "Paired"), block, list$annotation)
                aT1 <- aT1[aT1$ENTREZID %in% (unlist(mygslist)), ]
                # drop from esetm all duplicate genes and genes not in the genesets
                esetm = esetm[rownames(esetm) %in% aT1$ID, ]
                rownames(esetm) <- aT1$ENTREZID[match(rownames(esetm), aT1$ID)]
            } 

        runGSA = function(group) { 
            # Run GSA maxmean
            tab = table(group)
            nc = tab["c"]
            nd = tab["d"]
            if (list$design == "Not Paired") {
                yy = c(rep(1, nc), rep(2, nd))
            } else {
                yy = as.numeric(factor(block))
                yy[duplicated(yy)] = (-yy[duplicated(yy)])
            }
            
            resgsa = GSA(x = esetm, y = yy, genesets = mygslist, genenames = rownames(esetm), 
                method = "maxmean", resp.type = ifelse(list$design == "Not Paired", "Two class unpaired", 
                    "Two class paired"), censoring.status = NULL, random.seed = 1, knn.neighbors = 10, 
                s0 = NULL, s0.perc = NULL, minsize = minsize, maxsize = 1000, restand = TRUE, 
                restand.basis = c("catalog", "data"), nperms = NI)
            2 * apply(cbind(resgsa$pvalues.lo, resgsa$pvalues.hi), 1, min)
        }
        
        B = factor(block)
        o = order(B)
        pgrps = lapply(1:(Npsudo + 1), function(x) {
            G = group
            if (x > 1) {
                G[o] = unlist(tapply(G, B, sample, simplify=FALSE))
            }
            G
        })

        pres = lapply(pgrps, runGSA)
        pres = do.call(cbind, pres)
        
        res = data.frame(ID = names(mygslist), P = pres[,1], Dataset = list$dataset, stringsAsFactors = FALSE)
        rownames(pres) = res$ID
        res$Method = "GSA"
        ord = order(res$P)
        res = res[ord, ]
        res$Rank = rank(res$P)/nrow(res) * 100
        res$FDR = p.adjust(res$P, "fdr")
        rownames(res) = NULL

        pidx = res$ID %in% list$targetGeneSets
        list(targ = res[pidx, ], pval = pres[,-1,drop=FALSE])
    }
    
    defGSmethods = list(GSA = gsaF, PADOG = padogF)
    GSmethods = c(as.list(existingMethods), mymethods)
    names(GSmethods) = c(existingMethods, names(mymethods))
    defMeth = intersect(names(defGSmethods), names(GSmethods))
    GSmethods[defMeth] = defGSmethods[defMeth]
    GSMok = GSmethods
    GSMok[names(GSMok) == "AbsmT"] = defGSmethods["PADOG"]
    names(GSMok)[names(GSMok) == "AbsmT"] = "PADOG"
    GSMok = GSMok[!duplicated(names(GSMok))]
    refMethod = names(GSmethods)[1]
    
    
    # check GS
    if (length(gslist) == 1 && gslist == "KEGG.db") {
        pw2id = as.list(KEGGPATHID2EXTID)
        gslist = pw2id[grep(organism, names(pw2id))]
        names(gslist) = sub(paste("^", organism, sep = ""), "", names(gslist))
        gs.names = unlist(as.list(KEGGPATHID2NAME)[names(gslist)])
        rm(pw2id)
    }
    stopifnot(class(gslist) == "list")
    stopifnot(length(gslist) >= 3)
    if (!is.null(gs.names)) {
        stopifnot(length(gslist) == length(gs.names))
    }
    
    
    aggFun = function(zdat) {
        tmp = do.call(rbind, zdat)[, c("ID", "Rank", "P", "FDR", "Dataset", "Method")]
        rownames(tmp) = NULL
        tmp
    }
    
    dfr = list()
    
    if ("PADOG" %in% names(GSMok)) {
        dfr[["PADOG"]] = aggFun(lapply(files, GSMok[["PADOG"]], mygslist = gslist, 
            minsize = Nmin))
        GSMok = GSMok[names(GSMok) != "PADOG"]
    }
    
    if (use.parallel && requireNamespace("doParallel", quietly = TRUE) && requireNamespace("parallel", 
        quietly = TRUE)) {
        ncores = parallel::detectCores()
        if (!is.null(ncr)) 
            ncores = min(c(ncores, ncr))
        clust = parallel::makeCluster(ncores)
        doParallel::registerDoParallel(clust)
        tryCatch({
            parRes <- foreach(outi = seq_along(GSMok), .combine = "c", .packages = pkgs, .export = expVars) %:% 
                      foreach(ini = seq_along(files),  .combine = "c", .packages = pkgs, .export = expVars) %dopar% {
                              lapply(files[ini], GSMok[[outi]], mygslist = gslist, minsize = Nmin)
            }
            parRes = aggFun(parRes)
            parRes = split(parRes, parRes$Method)
            dfr[names(GSMok)] = parRes[names(GSMok)]
            rm(parRes)
        }, finally = parallel::stopCluster(clust))
    } else {
        if (use.parallel) message("Execute in sequential! Packages 'doParallel' and 'parallel' 
                                  needed for parallelization!")
        dfr[names(GSMok)] = lapply(GSMok, function(m) aggFun(lapply(files, m, mygslist = gslist, 
            minsize = Nmin)))
    }
    
    
    shared = Reduce(merge, lapply(dfr, function(z) {
        z = z[complete.cases(z), ]
        z = z[, c("Dataset", "ID")]
        z = z[!duplicated(z), ]
        z
    }))
    
    dfs = list()
    dfs[names(GSmethods)] = lapply(names(GSmethods), function(m) {
        if (m == "AbsmT") {
            retn = dfr[["PADOG"]]
            retn = retn[retn$Method == "AbsmT", ]
        } else if (m == "PADOG") {
            retn = dfr[["PADOG"]]
            retn = retn[retn$Method == "PADOG", ]
        } else {
            retn = dfr[[m]]
        }
        stopifnot(!any(duplicated(retn[, c("Dataset", "ID")])))
        retn = merge(shared, retn, all.x = TRUE)
        retn[order(retn$Dataset, retn$ID), ]
    })
    rm(dfr)
    
    
    psList <- lapply(dfs, function(x) {
        x$P
    })
    fdrList <- lapply(dfs, function(x) {
        x$FDR
    })
    rankList <- lapply(dfs, function(x) {
        x$Rank
    })
    targetgsList <- lapply(dfs, function(x) {
        x$ID
    })
    dsList <- lapply(dfs, function(x) {
        x$Dataset
    })
    
    
    
    wi = function(x) {
        if (!all(x == rankList[[refMethod]])) {
            wilcox.test(x, rankList[[refMethod]], paired = TRUE, alternative = "less")$p.value
        } else {
            1
        }
    }
    
    wioright = function(x) {
        if (!all(x == rankList[[refMethod]])) {
            dset = data.frame(Method = rep(c(1, 0), each = length(rankList[[refMethod]])), 
                Y = c(x, rankList[[refMethod]]), Dataset = factor(dsList[[refMethod]]), 
                Path = factor(c(targetgsList[[refMethod]], targetgsList[[refMethod]])))
            md = lme(Y ~ Method + Dataset, random = ~1 | Path, data = dset)
            re = summary(md)$tTable["Method", c(1, 5)]
            if (re[1] < 0) {
                c(re[1], re[2]/2)
            } else {
                c(re[1], 1 - re[2]/2)
            }
        } else {
            c(0, 1)
        }
    }
    
    
    
    
    repo = data.frame(t(sapply(rankList, wioright)))
    names(repo) <- c("coef. LME", "p LME")
    repo$"p Wilcox." <- sapply(rankList, wi)
    
    l05 <- function(x) round(mean(x < 0.05) * 100, 2)
    geomean <- function(x) {
        x = ifelse(x == 0, 1e-16, x)
        exp(mean(log(x)))
    }
    
    repo$"% p.value<0.05" <- lapply(psList, l05)
    repo$"% q.value<0.05" <- lapply(fdrList, l05)
    repo$"p geomean" <- lapply(psList, geomean)
    repo$"p med" <- lapply(psList, median)
    
    repo$"rank mean" <- lapply(rankList, mean)
    repo$"rank med" <- lapply(rankList, median)
    repo$Method <- names(psList)
    
    nmets = length(psList)
    
    somecols = c("lightgrey", "lightblue", "orange", "red", "blue", "grey")
    set.seed(1)
    if (nmets > 6) {
        somecols = c(somecols, sample(setdiff(colors(), somecols))[1:(nmets - 6)])
    }
    
    if (plots) {
        
        par(mfrow = c(1, 3))
        boxplot(psList, ylab = paste("p-value"), las = 3, col = somecols[1:nmets])
        boxplot(rankList, ylab = "Rank(%)", las = 3, col = somecols[1:nmets])
        
        mff2 = function(x) {
            x - rankList[[refMethod]]
        }
        newranks = lapply(rankList, mff2)
        newranks[refMethod] <- NULL
        
        if (length(newranks) == 1) {
            xlb = names(newranks)
        } else {
            xlb = NULL
        }
        boxplot(newranks, ylab = paste("Rank(%)-Rank ", refMethod, " (%)"), las = 3, 
            col = somecols[2:nmets], names = names(newranks), xlab = xlb)
        abline(h = 0)
        
    }
    
    out = repo[, c("Method", "p geomean", "p med", "% p.value<0.05", "% q.value<0.05", 
        "rank mean", "rank med", "p Wilcox.", "p LME", "coef. LME")]
    list(summary = out, ranks = rankList, pvalues = psList, qvalues = fdrList)
    
}



 
