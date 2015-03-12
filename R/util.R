
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
    
    padog2absmt = function(res, list, estFDR = c("BH", "Permutation", "holm")) {
        estFDR = match.arg(estFDR)
        res = res[complete.cases(res[,c("PmeanAbsT", "meanAbsT0")]),]
        res = res[order(res$PmeanAbsT, -res$meanAbsT0), ]
        res$Method = "AbsmT"
        res$Rank = sapply(1:nrow(res), function(n) {
            p = res$PmeanAbsT
            s = res$meanAbsT0
            ifelse(is.na(p[n]), NA, (sum(p < (p[n]), na.rm=TRUE) + 
                   sum(p == (p[n]) & s >= (s[n]), na.rm=TRUE)) / sum(! is.na(p), na.rm=TRUE) * 100 
            )
        }) 
        res$P = res$PmeanAbsT
        res$FDR = switch(estFDR, 
                         p.adjust(res$P, estFDR),
                         Permutation = res$FDRmeanAbsT
        )
        pidx = res$ID %in% list$targetGeneSets
        res[pidx, ]
    }

    getFDR = function(p0, p1) {
        # p1: observed p value vector; p0: permutation p value matrix.
        fdr0 = sapply(p1, function(z) {
                   ifelse(is.na(z), NA, min(sum(p0 <= z, na.rm=TRUE) * sum(!is.na(p1))
                   / sum(p1 <= z, na.rm=TRUE) / sum(!is.na(p0)), 1))
               })
        nna = !is.na(fdr0)
        fdr = fdr0[nna]
        ord = order(p1[nna], decreasing=TRUE)
        fdr[ord] = cummin(fdr[ord])
        fdr0[nna] = fdr
        fdr0
    }


