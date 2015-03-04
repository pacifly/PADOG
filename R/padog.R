
#' @import hgu133a.db hgu133plus2.db KEGG.db nlme
#' @import KEGGdzPathwaysGEO methods
NULL

#' @export
#'
#' @importFrom limma lmFit makeContrasts contrasts.fit eBayes topTable
#' @importFrom AnnotationDbi get as.list mappedkeys
#' @importFrom foreach foreach
#' @importFrom doRNG %dorng%
#'
padog <- function(esetm = NULL, group = NULL, paired = FALSE, block = NULL, gslist = "KEGG.db", 
    organism = "hsa", annotation = NULL, gs.names = NULL, NI = 1000, plots = FALSE, 
    targetgs = NULL, Nmin = 3, verbose = TRUE, paral = FALSE, dseed = NULL, ncr = NULL) {
    
    # initialize the gslist if using KEGG
    if (length(gslist) == 1 && gslist == "KEGG.db") {
        stopifnot(nchar(organism) == 3)
        pw2id = as.list(KEGGPATHID2EXTID)
        gslist = pw2id[grep(organism, names(pw2id))]
        names(gslist) = sub(paste("^", organism, sep = ""), "", names(gslist))
        gs.names = unlist(as.list(KEGGPATHID2NAME)[names(gslist)])
        stopifnot(length(gslist) >= 3)
        rm(pw2id)
    }
    
    # check arguments
    stopifnot(class(esetm) == "matrix")
    stopifnot(all(dim(esetm) > 4))
    
    stopifnot(class(group) %in% c("factor", "character"))
    stopifnot(length(group) == dim(esetm)[2])
    stopifnot(all(group %in% c("c", "d")))
    stopifnot(all(table(group) > 2))
    if (paired) {
        stopifnot(length(block) == length(group))
        stopifnot(all(table(block) == 2))
    }
    
    stopifnot(class(gslist) == "list")
    stopifnot(length(gslist) >= 3)
    if (!is.null(gs.names)) {
        stopifnot(length(gslist) == length(gs.names))
    }
    
    stopifnot(class(NI) == "numeric")
    stopifnot(NI > 5)
    
    if (plots) {
        stopifnot(targetgs %in% names(gslist))
    }
    if (!is.null(annotation)) {
        if (! annotation %in% c("hgu133a.db","hgu133plus2.db")) {
            stopifnot(require(annotation, character.only = TRUE))
        }
        stopifnot(sum(rownames(esetm) %in% mappedkeys(get(paste(substr(annotation, 
            1, nchar(annotation) - 3), "ENTREZID", sep = "")))) > 4)
    } else {
        stopifnot(sum(rownames(esetm) %in% as.character(unlist(gslist))) > 10 & !any(duplicated(rownames(esetm))))
    }
    
    # substitute some names
    Block = block
    
    # compute gene frequencies accross genesets
    gf = table(unlist(gslist))
    if (!all(gf == 1)) {
        if (quantile(gf, 0.99) > mean(gf) + 3 * sd(gf)) {
            gf[gf > quantile(gf, 0.99)] <- quantile(gf, 0.99)
        }
        gff <- function(x) {
            1 + ((max(x) - x)/(max(x) - min(x)))^0.5
        }
        # compute weights
        gf = gff(gf)
    } else {
        
        fdfd = unique(unlist(gslist))
        gf = rep(1, length(fdfd))
        names(gf) <- fdfd
    }
    
    
    allGallP = unique(unlist(gslist))
    
    
    if (!is.null(annotation)) {
        # get rid of duplicates in the esetm by choosing the probe(set) with lowest
        # p-value; get ENTREZIDs for probes
        aT1 = filteranot(esetm, group, paired, block, annotation)
        # drop genes not in any geneset drop from esetm all duplicate genes and genes not
        # in the genesets
        esetm = esetm[rownames(esetm) %in% aT1$ID, ]
        rownames(esetm) <- aT1$ENTREZID[match(rownames(esetm), aT1$ID)]
    }
    restg = setdiff(rownames(esetm), names(gf))
    appendd = rep(1, length(restg))
    names(appendd) <- restg
    gf = c(gf, appendd)
    
    
    
    stopifnot(all(!duplicated(rownames(esetm))))
    stopifnot(sum(rownames(esetm) %in% allGallP) > 10)
    
    if (verbose) {
        cat(paste("Starting with ", length(gslist), " gene sets!", sep = ""))
        cat("\n")
    }
    
    # drop pathways with less than Nmin genes
    gslist = gslist[unlist(lapply(gslist, function(x) {
        length(intersect(rownames(esetm), x)) >= Nmin
    }))]
    gs.names = gs.names[names(gslist)]
    stopifnot(length(gslist) >= 3)
    
    if (verbose) {
        cat(paste("Analyzing ", length(gslist), " gene sets with ", Nmin, " or more genes!", 
            sep = ""))
        cat("\n")
    }
    
    ############################################## 
    if (!is.null(dseed)) 
        set.seed(dseed)
    # compute scores for iterations
    G = factor(group)
    Glen = length(G)
    tab = table(G)
    idx = which.min(tab)
    minG = names(tab)[idx]
    minGSZ = tab[idx]
    bigG = rep(setdiff(levels(G), minG), length(G))
    block = factor(Block)
    # blockOrd = order(block)
    topSigNum = dim(esetm)[1]
    combFun = function(gi, countn = TRUE) {
        g = G[gi]
        tab = table(g)
        if (countn) {
            minsz = min(tab)
            ifelse(minsz > 10, -1, choose(length(g), minsz))
        } else {
            dup = which(g == minG)
            cms = combn(length(g), tab[minG])
            del = apply(cms, 2, setequal, dup)
            if (paired) {
                cms = cms[, order(del, decreasing = TRUE), drop = FALSE]
                cms[] = gi[c(cms)]
                cms
            } else {
                cms[, !del, drop = FALSE]
            }
        }
    }
    if (paired) {
        bct = tapply(seq_along(G), block, combFun, simplify = TRUE)
        nperm = ifelse(any(bct < 0), -1, prod(bct))
        if (nperm < 0 || nperm > NI) {
            ## combidx = replicate(NI, unlist(tapply(G, block, sample, simplify=FALSE))) #too
            ## slow
            btab = tapply(seq_along(G), block, `[`, simplify = FALSE)
            bSamp = function(gi) {
                g = G[gi]
                tab = table(g)
                bsz = length(g)
                minsz = tab[minG]
                cms = do.call(cbind, replicate(NI, sample.int(bsz, minsz), simplify = FALSE))
                cms[] = gi[c(cms)]
                cms
            }
            combidx = do.call(rbind, lapply(btab, bSamp))
        } else {
            bcomb = tapply(seq_along(G), block, combFun, countn = FALSE, simplify = FALSE)
            colb = expand.grid(lapply(bcomb, function(x) 1:ncol(x)))[-1, , drop = FALSE]
            combidx = mapply(function(x, y) x[, y, drop = FALSE], bcomb, colb, SIMPLIFY = FALSE)
            combidx = do.call(rbind, combidx)
        }
    } else {
        nperm = combFun(seq_along(G))
        if (nperm < 0 || nperm > NI) {
            combidx = do.call(cbind, replicate(NI, sample.int(Glen, minGSZ), simplify = FALSE))
        } else {
            combidx = combFun(seq_along(G), countn = FALSE)
        }
    }
    
    NI = ncol(combidx)
    if (verbose) {
        cat("# of permutations used:", NI, "\n")
    }
    
    deINgs = intersect(rownames(esetm), unlist(gslist))
    gslistINesetm = lapply(gslist, match, table = deINgs, nomatch = 0)
    MSabsT <- MSTop <- matrix(NA, length(gslistINesetm), NI + 1)
    gsScoreFun <- function(G, block) {
        # these two arguments are needed in parallel computing for the environment in
        # model.matrix
        force(G)
        force(block)
        if (ite > 1) {
            G = bigG
            G[combidx[, ite - 1]] = minG
            G = factor(G)
        }
        if (paired) {
            design <- model.matrix(~0 + G + block)
            colnames(design) <- substr(colnames(design), 2, 100)
        } else {
            design <- model.matrix(~0 + G)
            colnames(design) <- levels(G)
        }
        
        fit <- lmFit(esetm, design)
        cont.matrix <- makeContrasts(contrasts = "d-c", levels = design)
        fit2 <- contrasts.fit(fit, cont.matrix)
        fit2 <- eBayes(fit2)
        aT1 <- topTable(fit2, coef = 1, number = topSigNum)
        aT1$ID = rownames(aT1)
        de = abs(aT1$t)
        names(de) <- aT1$ID
        degf = scale(cbind(de, de * gf[names(de)]))
        rownames(degf) = names(de)
        degf = degf[deINgs, , drop = FALSE]
        sapply(gslistINesetm, function(z) {
            X = na.omit(degf[z, , drop = FALSE])
            colMeans(X, na.rm = TRUE) * sqrt(nrow(X))
        })
    }
    if (paral && requireNamespace("doParallel", quietly = TRUE) && requireNamespace("parallel", 
        quietly = TRUE)) {
        ncores = parallel::detectCores()
        if (!is.null(ncr)) 
            ncores = min(ncores, ncr)
        if (verbose) {
            clust = parallel::makeCluster(ncores, outfile="")
        } else {
            clust = parallel::makeCluster(ncores)
        }
        doParallel::registerDoParallel(clust)
        tryCatch({
            parRes = foreach(ite = 1:(NI + 1), .combine = "c", .packages = "limma") %dorng% 
                {
                  if (verbose && (ite %% 10 == 0)) {
                      cat(ite, "/", NI, "\n")
                  }
                  Sres <- gsScoreFun(G, block)
                  tmp <- list(t(Sres))
                  names(tmp) <- ite
                  tmp
                }
            parRes = do.call(cbind, parRes[order(as.numeric(names(parRes)))])
            evenCol = (1:ncol(parRes))%%2 == 0
            MSabsT[] = parRes[, !evenCol]
            MSTop[] = parRes[, evenCol]
            rm(parRes)
        }, finally = parallel::stopCluster(clust))
    } else {
        if (paral) message("Execute in sequential! Packages 'doParallel' and 'parallel' 
                           needed for parallelization!")
        for (ite in 1:(NI + 1)) {
            Sres <- gsScoreFun(G, block)
            MSabsT[, ite] <- Sres[1, ]
            MSTop[, ite] <- Sres[2, ]
            if (verbose && (ite %% 10 == 0)) {
                cat(ite, "/", NI, "\n")
            }
        }
    }
    ########################################################### 
    meanAbsT0 = MSabsT[, 1]
    padog0 = MSTop[, 1]
    plotIte = min(NI, 21)
    MSabsT_raw = MSabsT
    MSTop_raw = MSTop
    
    # standardize scores
    MSabsT = scale(MSabsT)
    MSTop = scale(MSTop)
    
    # compute p-values
    mff = function(x) {
        if (!all(is.na(x))) {
            mean(x[-1] > x[1], na.rm = TRUE)
        } else {
            NA
        }
    }
    PSabsT = apply(MSabsT, 1, mff)
    PSTop = apply(MSTop, 1, mff)
    PSabsT[PSabsT == 0] <- 1/NI/100
    PSTop[PSTop == 0] <- 1/NI/100
    
    # do plot the scores for the star pathway
    if (plots) {
        par(mfrow = c(2, 2))
        boxplot(MSabsT_raw[, 1:plotIte] ~ col(MSabsT_raw[, 1:plotIte]), col = c("lightblue", 
            rep("whitesmoke", NI)), names = c("0", 1:(plotIte - 1)), cex.axis = 0.8, 
            main = "ABSmT scores after first standardization", cex.main = 0.6)
        points(1:plotIte, MSabsT_raw[names(gslist) == targetgs, 1:plotIte], col = "red", 
            pch = 19)
        abline(h = MSabsT_raw[names(gslist) == targetgs, 1], col = "red")
        boxplot(MSTop_raw[, 1:plotIte] ~ col(MSTop_raw[, 1:plotIte]), col = c("lightblue", 
            rep("whitesmoke", NI)), names = c("0", 1:(plotIte - 1)), cex.axis = 0.8, 
            main = "PADOG after first standardization", cex.main = 0.6)
        points(1:plotIte, MSTop_raw[names(gslist) == targetgs, 1:plotIte], col = "red", 
            pch = 19)
        abline(h = MSTop_raw[names(gslist) == targetgs, 1], col = "red")
        boxplot(MSabsT[, 1:plotIte] ~ col(MSabsT[, 1:plotIte]), col = c("lightblue", 
            rep("whitesmoke", NI)), names = c("0", 1:(plotIte - 1)), cex.axis = 0.8, 
            main = "ABSmT scores second standardization", cex.main = 0.6)
        points(1:plotIte, MSabsT[names(gslist) == targetgs, 1:plotIte], col = "red", 
            pch = 19)
        abline(h = MSabsT[names(gslist) == targetgs, 1], col = "red")
        boxplot(MSTop[, 1:plotIte] ~ col(MSTop[, 1:plotIte]), col = c("lightblue", 
            rep("whitesmoke", NI)), names = c("0", 1:(plotIte - 1)), cex.axis = 0.8, 
            main = "PADOG after second standardization", cex.main = 0.6)
        points(1:plotIte, MSTop[names(gslist) == targetgs, 1:plotIte], col = "red", 
            pch = 19)
        abline(h = MSTop[names(gslist) == targetgs, 1], col = "red")
    }
    if (!is.null(gs.names)) {
        myn = gs.names
    } else {
        myn = names(gslist)
    }
    SIZE = unlist(lapply(gslist, function(x) {
        length(intersect(rownames(esetm), x))
    }))
    res = data.frame(Name = myn, ID = names(gslist), Size = SIZE, meanAbsT0, padog0, 
        PmeanAbsT = PSabsT, Ppadog = PSTop, stringsAsFactors = FALSE)
    res = res[order(res$Ppadog, -res$padog0), ]
    res
} 
