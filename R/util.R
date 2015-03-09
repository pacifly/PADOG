
    getdataaslist = function(x) {
        x = get(x, envir=parent.frame())
        exp = experimentData(x)
        dataset = exp@name
        disease = notes(exp)$disease
        dat.m = exprs(x)
        ano = pData(x)
        design = notes(exp)$design
        annotation = paste(x@annotation, ".db", sep = "")
        targetGeneSets = notes(exp)$targetGeneSets
        list = list(dataset, disease, dat.m, ano, design, annotation, targetGeneSets)
        names(list) = c("dataset", "disease", "dat.m", "ano", "design", "annotation", 
            "targetGeneSets")
        return(list)
    }
    
    padog2absmt = function(res, list, estFDR = FALSE) {
        res = res[complete.cases(res[,c("PmeanAbsT", "meanAbsT0")]),]
        res = res[order(res$PmeanAbsT, -res$meanAbsT0), ]
        res$Method = "AbsmT"
        res$Rank = (1:nrow(res))/nrow(res) * 100
        res$P = res$PmeanAbsT
        if (estFDR) {
            res$FDR = res$FDRmeanAbsT
        } else {
            res$FDR = p.adjust(res$P, "fdr")
        }
        pidx = res$ID %in% list$targetGeneSets
        res[pidx, ]
    }

    getFDR = function(p0, p1) {
        # p1: observed p value vector; p0: permutation p value matrix.
        fdr0 = sapply(p1, function(z) {
                   ifelse(is.na(z), na, max(sum(p0 <= z, na.rm=TRUE) * sum(!is.na(p1))
                   / sum(p1 <= z, na.rm=TRUE) / sum(!is.na(p0)), 1))
               })
        nna = !is.na(fdr0)
        fdr = fdr0[nna]
        ord = order(p1[nna], decreasing=TRUE)
        fdr[ord] = cummin(fdr[ord])
        fdr0[nna] = fdr
        fdr0
    }


