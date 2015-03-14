
if(getRversion() >= "2.15.1")  utils::globalVariables(c("ini", "outi"))


#' Benchmark for gene set analysis methods using 24 datasets
#'
#' This is a general purpose function to compare a given gene set analysis method in terms of 
#' sensitivity, ranking and fdr against PADOG and GSA (if installed) using 24 public datasets.
#'
#' The user defined gene set analysis method should return a list with elements "targ" containing 
#' the analysis statistics for the target gene sets on a given data set, and "pval" containing
#' a matrix of null p values each column of which corresponds to the p values of all gene sets
#' applying the gene set analysis method on permutated phenotype (i.e. group label). The rownames
#' of the matrix are the gene set IDs; the number of column corresponds to argument \code{Npsudo}.
#' If \code{Npsudo = 0}, the "pval" element can be omitted and null p values / fdr will not be
#' compared between the methods. For "PADOG" and "AbsmT", two methods ("BH" and "Permutation") of 
#' estimating fdr are available regardless of the value of \code{Npsudo} since null p values can be
#' accessed from the method algorithm itself. For "GSA", the "Permutation" method is only available 
#' when \code{Npsudo > 0}. 
#'
#' @param datasets A character vector with valid names of datasets to use from the PADOGsets 
#'   package. If left NULL all datasets avalibale in PADOGsets will be used.
#' @param existingMethods A character vector with one or more of the predefined methods 
#'   c("GSA","PADOG"). The first is used as reference method.
#' @param mymethods A list whose elements are valid functions implementing gene set analysis 
#'   methods. See \code{Details} and \code{Examples} for what arguments the functions have to 
#'   take in and what kind of output they need to return.
#' @param gs.names A character vector giving additional information about each gene set. For 
#'   instance when gene seta are pathways, the full name of the pathway would be a meaningful
#'   gene set name.
#' @param gslist Either the value "KEGG.db" or a list with the gene sets. If set to "KEGG.db",
#'   then gene sets will be made of all KEGG pathways for human since all datasets available 
#'   in PADOG are for human.
#' @param organism A three letter string giving the name of the organism supported by the 
#'   "KEGG.db" package.
#' @param Nmin The minimum size of gene sets to be included in the analysis for all methods.
#' @param NI Number of iterations to determine the gene set score significance p-values in 
#'   PADOG and GSA methods.
#' @param use.parallel If set to TRUE and multiple CPU cores are available, the parallelization 
#'   will be distributed either over the \code{NI} iterations for "PADOG", or over 
#'   datasets/methods combination otherwise.
#' @param ncr The number of CPU cores used when \code{use.parallel} set to TRUE. Default is 
#'   to use all CPU cores detected.
#' @param pkgs Character vector of packages that the \code{existingMethods} and \code{mymethods} 
#'   depend on (\code{NULL} for "PADOG" and "GSA"). Consult the \code{.packages} 
#'   argument in \code{foreach} function from \code{foreach} package.
#' @param expVars Character vector of variables to export. Consult the \code{.export} argument
#'   in \code{foreach} function from \code{foreach} package.
#' @param dseed Optional initial seed for random number generator (integer) used in \code{padog}.
#' @param Npsudo The number of permutations on phenotype (i.e. group label) to obtain null p
#'   values and estimate fdr. Set to 0 if not interested in estimating fdr from permutation.
#' @param FDRmeth The method used to estimate fdr; "BH" for Benjamini & Hochberg, "Permutation"
#'   for estimation based on observed and null p values. You can use "holm" for FWER as well.
#' @param plots If set to TRUE will plot the p values, ranks, the ranks differences w.r.t. 
#'   reference, null p values and fdr when applicable.
#' @param verbose If set to TRUE, for "PADOG" and "AbsmT" methods it will show the iterations 
#'   performed so far, and for others show the method and data set currently being executed
#'   (when \code{use.parallel = TRUE}, this only works for command line R, not Rgui).
#' @return A list of elements "summary","ranks","pvalues","qvalues" and "nullp" (if \code{Npsudo > 0}). 
#'   "summary" is a data frame containing: \code{Method} is the name of the gene set analysis method;
#'     \code{p geomean} geometric mean of nominal p-values for the target gene sets (gene sets expected 
#'     to be relevant); \code{p med} median of nominal p-values for the target genesets; \code{ \% p<0.05} 
#'     is the fraction of target gene sets significant at 0.05 level (this is the sensitivity); 
#'     \code{ \% q<0.05} is the fraction of target gene sets significant at 0.05 level after FDR correction;
#'     \code{rank mean} mean rank of the target gene sets; \code{rank med} median rank of the target gene sets;
#'     \code{p Wilcox.} p value from a Wilcoxon test paired at dataset level comparing the rank of target gene sets;
#'     \code{p LME} p value from a linear mixed effects (LME) model which unlike the Wilcoxon test above 
#'       accounts for the fact that ranks for the same gene set may be correlated;
#'     \code{coef LME} Coefficient from the LME model giving the difference in ranks of the target gene 
#'       sets between the current gene set analysis method and the reference method (the first method 
#'       in the \code{existingMethods} argument.
#'   "ranks", "pvalues", "qvalues": a list of elements containing for each method the ranks (
#'     nominal p values, fdr q values) of the target gene sets in the corresponding data sets.
#'   "nullp": a list of elements containing for each method the null p values (a list of two
#'     elements "ap" for all gene sets and "tp" for target gene sets) pooled over all data sets. 
#' @example /inst/examples/compFDR.R 
#'
#' @references Adi L. Tarca, Sorin Draghici, Gaurav Bhatti, Roberto Romero. Down-weighting 
#'    overlapping genes improves gene set analysis. BMC Bioinformatics, 2012.
#' @author
#'    Adi Laurentiu Tarca \email{atarca@@med.wayne.edu},
#'    Zhonghui Xu \email{zhonghui.xu@@gmail.com}
#'
#' @seealso See \code{\link{compPADOG}} for back-compatibility (i.e. no specificity/fdr comparison).
#'    Refer to latest development at \url{https://github.com/pacifly/PADOG}
#' @keywords nonparametric methods
#'
#' @export
#'
#' @importFrom AnnotationDbi get as.list 
#' @importFrom Biobase experimentData notes exprs pData 
#' @importFrom GSA GSA
#' @importFrom foreach foreach %dopar% %:%
#'
compFDR = function(datasets = NULL, existingMethods = c("GSA", "PADOG"), mymethods = NULL, 
    gs.names = NULL, gslist = "KEGG.db", organism = "hsa", Nmin = 3, NI = 1000, 
    use.parallel = TRUE, ncr = NULL, pkgs = NULL, expVars = NULL, dseed = NULL, 
    Npsudo = 20, FDRmeth = c("BH","Permutation","holm"), plots = FALSE, verbose = FALSE) {
   
    Npsudo = as.integer(Npsudo[1])
    Nmin = Nmin[1]
    NI = NI[1]
    stopifnot(Npsudo >= 0, Nmin > 0, NI > 1)
    FDRmeth = match.arg(FDRmeth)

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
        res$Rank = sapply(1:nrow(res), function(n) {
            p = res$Ppadog
            s = res$padog0
            ifelse(is.na(p[n]), NA, (sum(p < (p[n]), na.rm=TRUE) + 
                   sum(p == (p[n]) & s >= (s[n]), na.rm=TRUE)) / sum(! is.na(p), na.rm=TRUE) * 100 
            )
        }) 
        res$P = res$Ppadog
        res$FDR = switch(FDRmeth,
                         p.adjust(res$P, FDRmeth),
                         Permutation = res$FDRpadog
        )
        
        rownames(res) = NULL
        pidx = res$ID %in% list$targetGeneSets
        list(targ = rbind(res[pidx, ], padog2absmt(res, list, estFDR=FDRmeth)), pval = out$pval)
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

        KK = min(ncol(esetm), max(3, min(10, table(group))))

        runGSA = function(group) { 
            # Run GSA maxmean
            if (list$design == "Not Paired") {
                yy = as.numeric(factor(group))
            } else {
                yy = as.numeric(factor(block))
                yy = yy * as.numeric(as.character(factor(group, labels=c(-1, 1)) ))
            }
            
            resgsa = GSA(x = esetm, y = yy, genesets = mygslist, genenames = rownames(esetm), 
                method = "maxmean", resp.type = ifelse(list$design == "Not Paired", "Two class unpaired", 
                "Two class paired"), censoring.status = NULL, random.seed = dseed, knn.neighbors = KK, 
                s0 = NULL, s0.perc = NULL, minsize = minsize, maxsize = 1000, restand = TRUE, 
                restand.basis = c("catalog", "data"), nperms = NI)
            2 * apply(cbind(resgsa$pvalues.lo, resgsa$pvalues.hi), 1, min)
        }
        
        if (list$design == "Paired") {
            B = factor(block)
            o = order(B)
            pgrps = lapply(1:(Npsudo + 1), function(x) {
                G = group
                if (x > 1) {
                    G[o] = unlist(tapply(G, B, sample, simplify=FALSE))
                }
                G
            })
        } else {
            pgrps = c(list(group), replicate(Npsudo, sample(group), simplify=FALSE))
        }

        pres = lapply(pgrps, runGSA)
        pres = do.call(cbind, pres)
        
        res = data.frame(ID = names(mygslist), P = pres[,1], Dataset = list$dataset, stringsAsFactors = FALSE)

        if (Npsudo > 0) {
            res$FDRgsa = getFDR(pres[,-1,drop=FALSE], pres[,1]) 
        } else {
            res$FDRgsa = NA
        }

        rownames(pres) = res$ID
        res$Method = "GSA"
        ord = order(res$P)
        res = res[ord, ]
        res$Rank = sapply(res$P, function(p) ifelse(is.na(p), NA, mean(res$P <= p, na.rm=TRUE) * 100))
        res$FDR = switch(FDRmeth,
                         p.adjust(res$P, FDRmeth),
                         Permutation = res$FDRgsa
        )
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
        zdat = lapply(zdat, `[`, c("ID", "Rank", "P", "FDR", "Dataset", "Method"))
        tmp = do.call(rbind, zdat)
        rownames(tmp) = NULL
        tmp
    }
   
    aggP = function(zdat) {
        d1 = zdat[[1]]$pval
        ismat = is.matrix(d1)
        stopifnot(ismat || is.list(d1))
        atp = function(xdat) {
            list(ap = unlist(lapply(xdat, function(z) c(z$pval))),
                 tp = unlist(lapply(xdat, function(z) c(z$pval[z$targ$ID,]))))
        }
        if (ismat) {
            ret = atp(zdat)
        } else {
            ms = names(d1)
            ret = lapply(ms, function(m) {
                atp(lapply(zdat, function(z) list(targ = z$targ[z$targ$Method == m,], 
                                                  pval = z$pval[[m]])))
            })
            names(ret) = ms
        }
        ret
    }

    dfr = list()
    psim = list()

    if ("PADOG" %in% names(GSMok)) {
        tmpr = lapply(files, GSMok[["PADOG"]], mygslist = gslist, minsize = Nmin)
        dfr[["PADOG"]] = aggFun(lapply(tmpr, `[[`, "targ"))
        if (Npsudo > 0L) { psim[["PADOG"]] = aggP(tmpr) }
        GSMok = GSMok[names(GSMok) != "PADOG"]
    }
    
    if (use.parallel && requireNamespace("doParallel", quietly = TRUE) && requireNamespace("parallel", 
        quietly = TRUE)) {
        ncores = parallel::detectCores()
        if (!is.null(ncr)) 
            ncores = min(c(ncores, ncr))
        if (verbose) {
            clust = parallel::makeCluster(ncores, outfile="")
        } else {
            clust = parallel::makeCluster(ncores)
        }
        doParallel::registerDoParallel(clust)
        tryCatch({
            parRes <- foreach(outi = seq_along(GSMok), .combine = "c", .packages = pkgs, .export = expVars) %:% 
                      foreach(ini = seq_along(files),  .combine = "c", .packages = pkgs, .export = expVars) %dopar% {
                              inres = lapply(files[ini], GSMok[[outi]], mygslist = gslist, minsize = Nmin)
			      if (verbose) {
			          cat("Finish:", names(GSMok)[outi], " ------> ", files[ini], "\n")
			      }
                              inres
            }
            prt = aggFun(lapply(parRes, `[[`, "targ"))
            prt = split(prt, prt$Method)
            dfr[names(GSMok)] = prt[names(GSMok)]
            rm(prt)
            if (Npsudo > 0L) {
                ms = sapply(parRes, function(x) {
                    m = unique(x$targ$Method)
                    stopifnot(length(m) == 1)
                    m
                })
                parRes = lapply(names(GSMok), function(m) parRes[ms == m])
                psim[names(GSMok)] = lapply(parRes, aggP)
            }
            rm(parRes)
        }, finally = parallel::stopCluster(clust))
    } else {
        if (use.parallel) message("Execute in serial! Packages 'doParallel' and 'parallel' 
                                  needed for parallelization!")
        tmpr = lapply(GSMok, function(m) lapply(files, m, mygslist = gslist, minsize = Nmin))
        dfr[names(GSMok)] = lapply(tmpr, function(mres) aggFun(lapply(mres, `[[`, "targ")))
        if (Npsudo > 0L) {
            psim[names(GSMok)] = lapply(tmpr, function(mres) aggP(mres))
        }
    }
    
    
    shared = Reduce(merge, lapply(dfr, function(z) {
        z = z[complete.cases(z), ]
        z = z[, c("Dataset", "ID")]
        z = z[!duplicated(z), ]
        z
    }))
    
    dfs = list()
    dfs[names(GSmethods)] = lapply(names(GSmethods), function(m) {
        if (m %in% c("AbsmT", "PADOG")) {
            retn = dfr[["PADOG"]]
            retn = retn[retn$Method == m, ]
        } else {
            retn = dfr[[m]]
        }
        stopifnot(!any(duplicated(retn[, c("Dataset", "ID")])))
        retn = merge(shared, retn, all.x = TRUE)
        retn[order(retn$Dataset, retn$ID), ]
    })
    rm(dfr)
   
    if (Npsudo > 0L) {
        pperm = list()
        pperm[names(GSmethods)] = lapply(names(GSmethods), function(m) {
            if (m %in%  c("AbsmT", "PADOG")) {
                retn = psim[["PADOG"]][[m]]
            } else {
                retn = psim[[m]]
            }
            retn
        })
        rm(psim)
    }
    
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
            dset = data.frame(Method = gl(2, length(rankList[[refMethod]])), 
                Y = c(rankList[[refMethod]], x), Dataset = factor(rep(dsList[[refMethod]], 2)), 
                Path = factor(rep(targetgsList[[refMethod]], 2)))
            md = lme(Y ~ Method, random = ~1 | Path/Dataset, data = dset)
            re = summary(md)$tTable[2, c(1, 5)]
            if (re[1] < 0) {
                c(re[1], re[2]/2)
            } else {
                c(re[1], 1 - re[2]/2)
            }
        } else {
            c(0, 1)
        }
    }
    
    
    if (length(unique(targetgsList[[refMethod]])) == 1) {
        repo = data.frame(matrix(NA, nrow=length(rankList), ncol=2))
    } else {
        repo = data.frame(t(sapply(rankList, wioright)))
    }
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
    
    somecols = c("lightgrey", "lightblue", "orange", "red", "purple", "lightgreen")
    set.seed(1)
    if (nmets > 6) {
        somecols = c(somecols, sample(setdiff(colors(), somecols))[1:(nmets - 6)])
    }
    
    if (plots) {
        
        usrPar <- par(mfrow = c(1 + (Npsudo > 0L), 3))
        on.exit(par(usrPar))

        boxplot(psList, ylab = paste("p-value"), main="Target Gene Sets", las = 3, col = somecols[1:nmets])
        boxplot(rankList, ylab = "Rank(%)", main="Target Gene Sets", las = 3, col = somecols[1:nmets])
        
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
        boxplot(newranks, ylab = paste("Rank(%)-Rank ", refMethod, " (%)"), main="Target Gene Sets", las = 3, 
            col = somecols[2:nmets], names = names(newranks), xlab = xlb)
        abline(h = 0)
        


        if (Npsudo > 0L) {
            aps = lapply(pperm, `[[`, "ap")
            tps = lapply(pperm, `[[`, "tp")
            boxplot(aps, ylab = "null p-value", main="All Gene Sets", las = 3, col = somecols[1:nmets])
            boxplot(tps, ylab = "null p-value", main="Target Gene Sets", las = 3, col = somecols[1:nmets])
            boxplot(fdrList, ylab = "fdr", main="Target Gene Sets", las = 3, col = somecols[1:nmets])
        }


    }
    
    out = repo[, c("Method", "p geomean", "p med", "% p.value<0.05", "% q.value<0.05", 
        "rank mean", "rank med", "p Wilcox.", "p LME", "coef. LME")]
    if (Npsudo > 0L) {
        return(list(summary = out, ranks = rankList, pvalues = psList, qvalues = fdrList, nullp = pperm))
    } else {
        return(list(summary = out, ranks = rankList, pvalues = psList, qvalues = fdrList))
    }
    
}



 
