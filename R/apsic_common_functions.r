
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

getCancerTypes <- function() {
  base_folder = "../Data/"
  meta_data <- read.csv(paste0(base_folder, "MutationFiles/File_metadata.csv"), header = TRUE, stringsAsFactors = FALSE)
  meta_data$Primary_site
}




# selectCelllines <- function(panCancerData, text, filter=NA) {
selectCelllines <- function(panCancerData, patho_annot) {
  CellLine_annot <-  read.csv("../Data/ProjectDRIVE/TableS2.csv",stringsAsFactors = F)
  CellLine_annot = CellLine_annot[CellLine_annot$PATHOLOGIST_ANNOTATION == patho_annot, , drop=FALSE]
  
  if(nrow(CellLine_annot) == 0 ) {
    return(NA)
  }
  
  # if(is.na(filter)) {
  #   CellLine_annot = CellLine_annot[CellLine_annot$PRIMARY_SITE == text , ]
  #   
  # } else {
  #   CellLine_annot = CellLine_annot[CellLine_annot$PRIMARY_SITE == text & CellLine_annot$PATHOLOGIST_ANNOTATION == filter, ]
  # }
  
  
  
  indexes = which(colnames(panCancerData$viabilities) %in% paste0(CellLine_annot$CELLLINE, "_", CellLine_annot$PRIMARY_SITE))
  
  panCancerData$viabilities = panCancerData$viabilities[, indexes, drop=FALSE]
  panCancerData$mutations_all = panCancerData$mutations_all[, indexes, drop=FALSE]
  panCancerData$silentMutations = panCancerData$silentMutations[, indexes, drop=FALSE]
  panCancerData$missenseMutations = panCancerData$missenseMutations[, indexes, drop=FALSE]
  panCancerData$truncatingMutations = panCancerData$truncatingMutations[, indexes, drop=FALSE]
  
  panCancerData$copyNumbers = panCancerData$copyNumbers[, indexes, drop=FALSE]
  # panCancerData$exprData = panCancerData$exprData[, indexes]
  panCancerData
}


##### plotting tcga data
##### plotting tcga data
boxplot_gene <- function(tcga_data, gene) {
  gene_index = which(rownames(tcga_data$normal_counts)== gene) 
  
  tumor_cpm = cpm(tcga_data$tumor_counts)
  normal_cpm = cpm(tcga_data$normal_counts)
  
  N = normal_cpm[gene_index, ] + 1
  T = tumor_cpm[gene_index, ] + 1
  
  dat = data.frame(group = c(rep( "Normal", length(N)), rep( "Tumor", length(T)) ) , CPM = c(N, T) ) 
  
  # geom_boxplot proposes several arguments to custom appearance
  gPlot = ggplot(dat, aes(x=group, y=CPM)) + 
    geom_boxplot(
      
      # custom boxes
      color="blue",
      fill="blue",
      alpha=0.2,
      
      # Notch?
      # notch=TRUE,
      # notchwidth = 0.8,
      
      # custom outliers
      outlier.colour="red",
      outlier.fill="red",
      outlier.size=3
      
    ) + scale_y_continuous(trans='log10') + theme(axis.title.x=element_blank())
  
  print(gPlot)
}


mapPrimarySiteToTCGAproject <- function(primary_site) {
  data = read.csv("../Data/TCGA/TCGAProjectsDesc.csv", stringsAsFactors = F)
  rownames(data) = data$primary_site
  
  data[primary_site, ]$tcga_project
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

identifyDependencies <- function(cancerData,  dependencyType ) {
  if(dependencyType == "non-genetic-tsg") {
    type = "non-genetic"; direction = "high"
  }
  if(dependencyType == "non-genetic-oncogene") {
    type = "non-genetic"; direction = "low"
  }
  if(dependencyType == "mutation-tsg") {
    type = "deleterious"; direction = "high"
  }
  if(dependencyType == "mutation-oncogene") {
    type = "missense"; direction = "low"
  }
  if(dependencyType == "amplification-oncogene") {
    type = "amplification"; direction = "low"
  }
    
      
  allRanks = apply(cancerData$viabilities, 2, function(x) {n = length(x);
  normalize(rank(x), method = "range", range = c(1/n, 1-1/n))})
  genes = rownames(cancerData$viabilities)
  results = matrix("", 0, 4)
  colnames(results) = c("gene", "freq_wt", "freq_mut", "pvalue")
  
  for(gene in genes) {
    
    if(type != "non-genetic") {
      
      if(type == "missense") {
        mutatedIndexes = which(cancerData$missenseMutations[gene, ] == 1  ) 
      } else if(type == "deleterious") {
        mutatedIndexes = which(cancerData$truncatingMutations[gene, ] == 1  ) 
      } else if(type == "amplification") { 
        mutatedIndexes = which(cancerData$copyNumbers[gene, ] == 2  ) # high amplification
      } 
      wtIndexes = setdiff(1:ncol(cancerData$viabilities), mutatedIndexes)
      
      if(length(mutatedIndexes) == 0 | length(wtIndexes) == 0) 
      {
        pvalue = NA
      } else {
        ## compare mutated and wild-types
        V_mut = allRanks[gene, mutatedIndexes]
        V_wt = allRanks[gene, wtIndexes]
        pvalue=CDF_sum_weighted_uniform(V_mut, V_wt, direction)
      }
      
      results = rbind(results, c(gene,  length(wtIndexes), length(mutatedIndexes), pvalue))
    } else if( type == "non-genetic" ) {

      wtIndexes = which(cancerData$truncatingMutations[gene, ] == 0 & cancerData$missenseMutations[gene, ] == 0 &
                          ( abs(cancerData$copyNumbers[gene, ]) != 2 | is.na(cancerData$copyNumbers[gene, ]) ) ) # copy neutral
      
      p_wt = IH_CDF(sum(allRanks[gene, wtIndexes]), length(wtIndexes) )

      if(direction == "high") {
        p_wt  = 1 - p_wt
      }
      
      freq_mut = ncol(cancerData$viabilities) - length(wtIndexes)
      
      p_wt = ifelse(length(wtIndexes) == 0, NA, p_wt)
      
      results = rbind(results, c(gene,  length(wtIndexes), freq_mut,  p_wt))
    }
  }
  
  
  results_df = data.frame(results, stringsAsFactors = F)
  row.names(results_df) = results_df$gene
  results_df = results_df[ , !(names(results_df) %in% "gene")]
  results_df$freq_wt = as.numeric(results_df$freq_wt)
  results_df$freq_mut = as.numeric(results_df$freq_mut)
  results_df$pvalue = as.numeric(results_df$pvalue)
  
  results_df = results_df[order(as.numeric(results_df$pvalue)), ]
  
  nrOfCelllines = ncol(cancerData$viabilities)
  nrOfCelllinesThr = max(2, ceiling(nrOfCelllines/100))
  if(type != "non-genetic") {
    results_df = results_df[which(results_df$freq_mut >= nrOfCelllinesThr), ]
  } else {
    results_df = results_df[which(results_df$freq_wt >= nrOfCelllinesThr), ]
  }
  
  pval_significance_thr = min(0.05, 1/nrow(results_df))
  results_df$is_significant = ifelse(results_df$pvalue < pval_significance_thr, "yes", "no")
  
  colnames(results_df) = c("#wt", "#mut", "pvalue",	"is_significant")
  results_df
}


############# plotting
getWaterfallSettings <- function(type = "mut", cols=NULL) {
  if(type == "mut") {
    if(is.null(cols)) {
      cols =  c("darkgray","#2c7bb6", "#fdae61", "#d7191c")
    }
    res =   list(colors = cols, legendText = c("WT", "Others", "Missense", "Truncating"))
  } else if(type == "cna") {
    if(is.null(cols)) {
      cols = c("darkgray","darkblue", "lightblue",  "orange", "darkred")
    }
    res = list(colors = cols, 
               legendText = c("no change", "2 copy deletion", "1 copy deletion", "amplification", "high-amplification"))
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
      warnings("no observation!")
      return(NULL)
    }
    
  } else if(type == "only_mut") {
    indexes = which(mut!=1) # !=1 (i.e. 2, 3, 4) means mutated cases
    if(length(indexes) == 0) {
      warning("no observation!")
      return(NULL)
    }
    mut = mut[indexes]
    via = via[indexes]
  } else if(type == "only_truncating") {
    indexes = which(mut==4) 
    if(length(indexes) == 0) {
      warning("no observation!")
      return(NULL)
    } 
    mut = mut[indexes]
    via = via[indexes]
  } else if(type == "only_missense") {
    indexes = which(mut==3) 
    if(length(indexes) == 0) {
      warning("no observation!")
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
}




waterfallForGene_CNA <- function(panCancerData, gene, title, rank, legenedPos="bottomleft", cols=NULL, type="all",
                                 sig_alpha = NA, cex.axis=1.1) {
  cna_tmp = rep(1, ncol(panCancerData$copyNumbers))
  cna_tmp[ panCancerData$copyNumbers[gene,] < 1 ] = 2
  cna_tmp[ panCancerData$copyNumbers[gene,] == 1 ] = 3
  cna_tmp[ panCancerData$copyNumbers[gene,] == 3 ] = 4
  cna_tmp[ panCancerData$copyNumbers[gene,] > 3 ] = 5
  
  
  if(rank == TRUE){
    allRanks = apply(panCancerData$viabilities, 2, rank) / nrow(panCancerData$viabilities) - 0.5
    via_tmp = as.numeric(allRanks[gene,])
  } else {
    via_tmp = as.numeric(panCancerData$viabilities[gene,])
  }
  
  
  via = via_tmp[order(via_tmp, decreasing=TRUE)]
  cna = cna_tmp[order(via_tmp, decreasing=TRUE)]
  
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
}
