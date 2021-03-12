
IH_CDF <- function(x, n) {
  if(n <= 20){
    X <-  floor(x)
    k <- seq(from = 0, to = X)
    # compute the cdf or p-value for n <= 100
    s <-  (-1)^k * choose(n, k)*( (x-k)^n)
    return(sum(s)/factorial(n))
  }else{
    # approximation for large n
    return(pnorm(x, n/2,sqrt(n/12), lower.tail = TRUE))
  }
}

# selectCelllines <- function(panCancerData, text, filter=NA) {
selectCelllines <- function(panCancerData, patho_annot) {
  if(patho_annot == "Pan_cancer") {
    return(panCancerData)
  }
  # CellLine_annot <-  read.csv("../external/InputData/ProjectDRIVE/TableS2.csv",stringsAsFactors = F)
  # CellLine_annot = CellLine_annot[CellLine_annot$PATHOLOGIST_ANNOTATION == patho_annot, , drop=FALSE]
  CellLine_annot = panCancerData$CellLine_annot
  CellLine_annot = CellLine_annot[CellLine_annot$PATHOLOGIST_ANNOTATION == patho_annot, , drop=FALSE]
  
  if(nrow(CellLine_annot) == 0 ) {
    return(NA)
  }
  
  
  indexes = which(colnames(panCancerData$viabilities) %in% paste0(CellLine_annot$CELLLINE, "_", CellLine_annot$PRIMARY_SITE))
  
  panCancerData$viabilities = panCancerData$viabilities[, indexes, drop=FALSE]
  panCancerData$D2_scores = panCancerData$D2_scores[, indexes, drop=FALSE]
  panCancerData$mutations_all = panCancerData$mutations_all[, indexes, drop=FALSE]
  panCancerData$silentMutations = panCancerData$silentMutations[, indexes, drop=FALSE]
  panCancerData$missenseMutations = panCancerData$missenseMutations[, indexes, drop=FALSE]
  panCancerData$missenseMutationsHotspot = panCancerData$missenseMutationsHotspot[, indexes, drop=FALSE]
  panCancerData$allTruncatingMutations = panCancerData$allTruncatingMutations[, indexes, drop=FALSE]
  panCancerData$lofMutations = panCancerData$lofMutations[, indexes, drop=FALSE]
  panCancerData$gofTruncatingMutations = panCancerData$gofTruncatingMutations[, indexes, drop=FALSE]
  panCancerData$copyNumbers = panCancerData$copyNumbers[, indexes, drop=FALSE]
  # panCancerData$exprData = panCancerData$exprData[, indexes]
  panCancerData$rankedViabilities = panCancerData$rankedViabilities[, indexes, drop=FALSE]
  panCancerData$CRISPR_data = panCancerData$CRISPR_data[, indexes, drop=FALSE]
  panCancerData$wildtypeMatrix = panCancerData$wildtypeMatrix[, indexes, drop=FALSE]
  
  panCancerData$mutationFreqMatrix$missenseMutNr = panCancerData$mutationFreqMatrix$missenseMutNr[, indexes, drop=FALSE]
  panCancerData$mutationFreqMatrix$deleteriousMutNr = panCancerData$mutationFreqMatrix$deleteriousMutNr[, indexes, drop=FALSE]
  panCancerData$mutationFreqMatrix$silentMutNr = panCancerData$mutationFreqMatrix$silentMutNr[, indexes, drop=FALSE]
  panCancerData$mutationFreqMatrix$nonsenseMutNr = panCancerData$mutationFreqMatrix$nonsenseMutNr[, indexes, drop=FALSE]
  panCancerData$mutationFreqMatrix$indelSplicMutNr = panCancerData$mutationFreqMatrix$indelSplicMutNr[, indexes, drop=FALSE]
  
  panCancerData$rankedCRISPR_data = panCancerData$rankedCRISPR_data[, indexes, drop=FALSE]
  panCancerData$rankedD2_scores = panCancerData$rankedD2_scores[, indexes, drop=FALSE]
  
  panCancerData
}

CDF_sum_weighted_uniform <- function(x, y, direction) {
  m = length(x)
  n = length(y)
  if((m+n) >= 20) {
    diff = mean(x) - mean(y)
    sd = sqrt( ((1/m ) + (1/n))/12)   
    
    pvalue = pnorm(diff, 0, sd = sd)
    if(direction == "high") {
      pvalue = 1 - pvalue
    }
  } else {
    # sd = sqrt((1/m + 1/n)/12)
    sd = 1
    stat_obs = (mean(x) - mean(y))/sd
    
    a_bar = (-1)^n/(m^m * n^n *sd^(n+m) )
    coef = (-1)^(m+n)/(a_bar * factorial(m+n))
    S_tot = 0
    
    if( n == 0 | m == 0 ) {
      stop("error")
    }
    
    for(k in 0:m) {
      for(p in 0:n) {
        S_pk =   (p/n - k/m ) / sd
        t = stat_obs  - S_pk
        
        theta_t = 0 
        if(t == 0)
          theta_t = 1/2
        else if(t > 0)
          theta_t = 1
        
        s_bar_pk = ((-1)^(p+k))
        
        S_tot = S_tot + choose(m, k) * choose(n, p)* s_bar_pk * (t)^(m+n) * theta_t
      }
    }
    pvalue = coef * S_tot
    if(direction == "high") {
      pvalue =  1- pvalue
    }
    
  }
  pvalue
}


computeRankViabilities <- function(viabilities) {
  apply(viabilities, 2, function(x) {
    if(all(is.na(x))) {
      r = rep(NA, length(x))
      names(r) = names(x)
      return(r)
    }
    n = length(x);
    normalize(rank(x, na.last = "keep"), method = "range", range = c(1/n, 1-1/n))})
}

getMutatedIndexes <- function(gene, apsic_type, geneticData, viabilities) {
  if(apsic_type == "missense") {
    mutatedIndexes = which(geneticData$missenseMutations[gene, ] == 1 & !is.na(viabilities[gene, ]) ) 
  } else if(apsic_type == "truncating_all") {
    mutatedIndexes = which(geneticData$allTruncatingMutations[gene, ] == 1 & !is.na(viabilities[gene, ]) ) 
  } else if(apsic_type == "truncating_lof") {
    mutatedIndexes = which(geneticData$lofMutations[gene, ] == 1 & !is.na(viabilities[gene, ]) ) 
  } else if(apsic_type == "truncating_gof") {
    mutatedIndexes = which(geneticData$gofTruncatingMutations[gene, ] == 1 & !is.na(viabilities[gene, ]) ) 
  } else if(apsic_type == "amplification") { 
    mutatedIndexes = which(geneticData$copyNumbers[gene, ] == 2 & !is.na(viabilities[gene, ])) # high amplification
  } else if(apsic_type == "non-genetic") {
    mutatedIndexes = NULL
  } else {
    stop("unknown apsic_type or non-genetic")
  }
  mutatedIndexes
}

getWildtypeIndexes <- function(gene, geneticData, viabilities) {
  which(geneticData$wildtypeMatrix[gene, ] == 1 & !is.na(viabilities[gene, ]) ) 
}


identifyDependencies <- function(cancerData,  dependencyType, excludeNotApplicable=TRUE ) {
  if(dependencyType == "tumor-suppressive-effectors") {
    type = "non-genetic"; direction = "high"
  }
  if(dependencyType == "tumor-promoting-effectors") {
    type = "non-genetic"; direction = "low"
  }
  if(dependencyType == "neomorphic-mutational-oncogenes") {
    type = "truncating_all"; direction = "low"
  }
  if(dependencyType == "mutational-oncogenes") {
    type = "missense"; direction = "low"
  }
  if(dependencyType == "amplified-oncogenes") {
    type = "amplification"; direction = "low"
  }
  
  geneticData = cancerData
  rViabilities = allRanks = computeRankViabilities(cancerData$viabilities)
  genes = rownames(cancerData$viabilities)
  results = data.frame(matrix(NA, 0, 4))
  colnames(results) = c("gene", "freq_wt", "freq_mut", "pvalue")
  
  for(gene in genes) {
    wtIndexes = getWildtypeIndexes(gene, cancerData, rViabilities) 
    
    if(type != "non-genetic") {
      mutatedIndexes = getMutatedIndexes(gene, type, geneticData, rViabilities)
      pvalue = NA
      
      if(length(mutatedIndexes) > 0 & length(wtIndexes) > 0) {
        ## compare mutated and wild-types
        V_mut = rViabilities[gene, mutatedIndexes]
        V_wt = rViabilities[gene, wtIndexes]
        pvalue=CDF_sum_weighted_uniform(V_mut, V_wt, direction)
      }
      
      # results = rbind(results, c(gene,  length(wtIndexes), length(mutatedIndexes), pvalue, NA))
      results[nrow(results)+1, ] = c(gene,  length(wtIndexes), length(mutatedIndexes), pvalue)
      
    } else if( type == "non-genetic" ) {
      p_wt = IH_CDF(sum(rViabilities[gene, wtIndexes]), length(wtIndexes) )
      
      if(direction == "high") {
        p_wt  = 1 - p_wt
      }
      
      freq_mut = sum(!is.na(rViabilities[gene, ])) - length(wtIndexes)
      p_wt = ifelse(length(wtIndexes) == 0, NA, p_wt)
      
      results[nrow(results)+1, ] = c(gene,  length(wtIndexes), freq_mut,  p_wt)
    }
  }
  row.names(results) = results$gene
  results$freq_wt = as.numeric(results$freq_wt)
  results$freq_mut = as.numeric(results$freq_mut)
  results$pvalue = as.numeric(results$pvalue)
  
  results_df = results[ , !(colnames(results) %in% "gene")]
  
  nrOfCelllines = sum(apply(rViabilities, 2, function(x){!all(is.na(x)) }))
  # nrOfCelllinesThr = max(2, ceiling(nrOfCelllines/100))
  nrOfCelllinesThr = 2
  if(type != "non-genetic") {
    indexes = which(results_df$freq_mut >= nrOfCelllinesThr & results_df$freq_wt >= nrOfCelllinesThr)  
    
    results_df[setdiff(1:nrow(results_df), indexes), 3] = NA
    
    if(excludeNotApplicable)
      results_df = results_df[indexes, ]
  } else {
    indexes = which(results_df$freq_wt >= nrOfCelllinesThr)
    results_df[setdiff(1:nrow(results_df), indexes), 3] = NA
    
    if(excludeNotApplicable)
      results_df = results_df[indexes, ]
  }
  
  nrTests = sum(!is.na(results_df$pvalue))
  pval_significance_thr = min(0.05, 1/nrTests)
  
  results_df$is_significant = ifelse(results_df$pvalue < pval_significance_thr, "yes", "no")
  colnames(results_df)[1:2] = c("#wt", "#mut")
  results_df = results_df[order(results_df$pvalue), ]
  # results_df$qvalue = p.adjust(results_df$pvalue, "fdr")
  
  results_df
}


############# plotting
getWaterfallSettings <- function(type = "mut", cols=NULL) {
  if(type == "mut") {
    if(is.null(cols)) {
      cols =  c("darkgray","#2c7bb6", "#fdae61", "#d7191c")
    }
    res =   list(colors = cols, legendText = c("WT", "Others", "Missense", "Non-missense"))
  } else if(type == "cna") {
    if(is.null(cols)) {
      cols = c("darkgray","darkblue", "lightblue",  "orange", "darkred")
    }
    res = list(colors = cols, 
               legendText = c("No change", "2 copy deletion", "1 copy deletion", "Amplification", "High-amplification"))
  }
  res
}

plotRandomRankBands <- function(n, handle, sig_alpha) {
  if(is.na(sig_alpha)) {
    sig_alpha = min(0.05, 1/n)
  }
  
  alpha = 1:n
  beta = n:1
  l = qbeta(sig_alpha/2, alpha, beta)
  h = qbeta(1-sig_alpha/2, alpha, beta)
  
  polygon(c(handle, rev(handle)), c( rev(l)-0.5,  rev(rev(h)-0.5) ), lty = 2 , lwd=2, col =rgb(1, 0, 0,0.35))
}


waterfallForGene <- function(panCancerData, gene, title, rank, legenedPos="bottomleft", 
                             cols=NULL, type="all", sig_alpha = NA, cex.axis=1.1) {
  indexes = which(!is.na(panCancerData$viabilities[gene, ]))
  if(length(indexes) > 0) {
    panCancerData$mutations_all = panCancerData$mutations_all[, indexes]
    panCancerData$silentMutations = panCancerData$silentMutations[, indexes]  
    panCancerData$missenseMutations = panCancerData$missenseMutations[, indexes]  
    panCancerData$truncatingMutations = panCancerData$allTruncatingMutations[, indexes]  
    panCancerData$viabilities = panCancerData$viabilities[, indexes]  
    
  } else {
    return(NULL)
  }
  
  
  
  mut_tmp = rep(1, ncol(panCancerData$mutations_all))
  mut_tmp[ which(panCancerData$silentMutations[gene, ] == 1) ] = 2
  mut_tmp[ which(panCancerData$missenseMutations[gene, ] == 1) ] = 3
  mut_tmp[ which(panCancerData$truncatingMutations[gene, ] == 1) ] = 4
  
  
  if(rank == TRUE){
    allRanks = apply(panCancerData$viabilities, 2, rank) / nrow(panCancerData$viabilities) - 0.5
    via_tmp = as.numeric(allRanks[gene,])
  } else {
    via_tmp = as.numeric(panCancerData$viabilities[gene,])
  }
  
  via = via_tmp[order(via_tmp, decreasing=TRUE)]
  mut = mut_tmp[order(via_tmp, decreasing=TRUE)]
  
  if(type == "only_wt") {
    indexes = which(mut==1)  # 1 means wt
    mut = mut[indexes]
    via = via[indexes]
    if(length(indexes) == 0) {
      # warnings("no observation!")
      return(NULL)
    }
    
  } else if(type == "only_mut") {
    indexes = which(mut!=1) # !=1 (i.e. 2, 3, 4) means mutated cases
    if(length(indexes) == 0) {
      # warning("no observation!")
      return(NULL)
    }
    mut = mut[indexes]
    via = via[indexes]
  } else if(type == "only_truncating") {
    indexes = which(mut==4) 
    if(length(indexes) == 0) {
      # warning("no observation!")
      return(NULL)
    } 
    mut = mut[indexes]
    via = via[indexes]
  } else if(type == "only_missense") {
    indexes = which(mut==3) 
    if(length(indexes) == 0) {
      # warning("no observation!")
      return(NULL)
    } 
    mut = mut[indexes]
    via = via[indexes]
  }
  
  wfSettings = getWaterfallSettings(cols=cols)
  col = wfSettings$colors[mut]
  
  legText = wfSettings$legendText[as.numeric(names(table(mut)> 0))]
  legColors = wfSettings$colors[as.numeric(names(table(mut)> 0))]
  
  if(rank == TRUE){
    ylim = c(-0.6, 0.6)
    ylab = "Ranked Gene Score"
    yaxt = "n"
  } else {
    ylim=c(-max(abs(via)), max(abs(via)))
    ylab = "Gene Score"
    yaxt = NULL
  }
  handle = barplot(via, col=col, border=col, space=0.5, ylim=ylim,
                   main =title, ylab=ylab, yaxt=yaxt,
                   cex.axis=1.2, cex.lab=1.4, legend.text=legText,
                   args.legend=list(title="", fill=legColors, border=NA, cex=0.6, x = legenedPos,  bty = "n",xpd=FALSE ))
  
  n =  length(via)
  if(rank == TRUE){
    # add rank labels    
    axis_y = seq(-0.5, 0.5, length.out=5) 
    axis(2, at=axis_y, labels=c("0", "0.25", "0.5", "0.75", "1"), cex.axis=cex.axis)
    
    plotRandomRankBands(n, handle, sig_alpha)
  }
  TRUE
}




waterfallForGene_CNA <- function(panCancerData, gene, title, rank, legenedPos="bottomleft", cols=NULL, type="all",
                                 sig_alpha = NA, cex.axis=1.1) {
  indexes = which(!is.na(panCancerData$viabilities[gene, ]))
  if(length(indexes) > 0) {
    panCancerData$copyNumbers = panCancerData$copyNumbers[, indexes]
    panCancerData$viabilities = panCancerData$viabilities[, indexes]  
  } else {
    return(NULL)
  }
  
  cna_tmp = rep(1, ncol(panCancerData$copyNumbers))
  cna_tmp[ panCancerData$copyNumbers[gene,] == -2 ] = 2
  cna_tmp[ panCancerData$copyNumbers[gene,] == -1 ] = 3
  cna_tmp[ panCancerData$copyNumbers[gene,] == 1 ] = 4
  cna_tmp[ panCancerData$copyNumbers[gene,] == 2 ] = 5
  
  
  if(rank == TRUE){
    allRanks = apply(panCancerData$viabilities, 2, rank) / nrow(panCancerData$viabilities) - 0.5
    via_tmp = as.numeric(allRanks[gene,])
  } else {
    via_tmp = as.numeric(panCancerData$viabilities[gene,])
  }
  
  
  via = via_tmp[order(via_tmp, decreasing=TRUE)]
  cna = cna_tmp[order(via_tmp, decreasing=TRUE)]
  
  
  if(type == "only_wt") {
    indexes = which(cna==1)  # 1 means wt
    via = via[indexes]
    cna = cna[indexes]
    if(length(indexes) == 0) {
      # warnings("no observation!")
      return(NULL)
    }
    
  } else if(type == "only_cna") {
    indexes = which(cna!=1) # !=1 (i.e. 2, 3, 4, 5) means with CNA
    if(length(indexes) == 0) {
      # warning("no observation!")
      return(NULL)
    }
    cna = cna[indexes]
    via = via[indexes]
  }
  
  
  wfSettings = getWaterfallSettings("cna",cols=cols)
  col = wfSettings$colors[cna]
  
  
  legText = wfSettings$legendText[as.numeric(names(table(cna)> 0))]
  legColors = wfSettings$colors[as.numeric(names(table(cna)> 0))]
  
  
  if(rank == TRUE){
    ylim = c(-0.6, 0.6)
    ylab = "Ranked Gene Score"
    yaxt = "n"
  } else {
    ylim=c(-max(abs(via)), max(abs(via)))
    ylab = "Gene Score"
    yaxt = NULL
  }
  
  handle = barplot(via, col=col, border=col, space=0.5, ylim=ylim,
                   main =title, ylab=ylab, yaxt=yaxt,
                   cex.axis=1.2, cex.lab=1.4, legend.text=legText,
                   args.legend=list(title="", fill=legColors, border=NA, cex=0.6, x = legenedPos,  bty = "n",xpd=FALSE ))
  
  n =  length(via)
  
  # compute bands according to the beta distribution  
  if(rank == TRUE){
    axis_y = seq(-0.5, 0.5, length.out=5) 
    axis(2, at=axis_y, labels=c("0", "0.25", "0.5", "0.75", "1"), cex.axis=cex.axis)
    
    plotRandomRankBands(n, handle, sig_alpha)
  }
  TRUE
}
