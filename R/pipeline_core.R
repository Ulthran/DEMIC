#' A function representing the pipeline of four steps
#' including GC bias correction, sample filtration, PCA and contig filtration
#' @param Y a matrix of coverages
#' @param i cutoff of filtering samples changes according to parameter i; i=1, cutoffRatio is 0.5; i=2, cutoffRatio is 1 as contig is clean
#' @return final filtered samples, matrix of sample, contig and corrected coverages,
#' filtered contigs with PC1 values, PC1 range, preliminary filtered samples
#'
#' @importFrom stats coef cor cor.test ks.test p.adjust aggregate
pipeline <- function(Y, i) {
  lmeModel <- lme4::lmer(log_cov ~ GC_content + (1 | sample:contig), data = Y, REML = FALSE)
  summeryMeanY <- aggregate(GC_content ~ (sample:contig), Y, FUN = "mean")
  summeryMeanY$s_c <- paste(summeryMeanY$sample, summeryMeanY$contig, sep = ":")

  lmeModelCoef <- coef(lmeModel)$`sample:contig`
  lmeModelCoef$s_c <- rownames(lmeModelCoef)

  summeryMeanYSort <- merge(lmeModelCoef, summeryMeanY, by = "s_c")
  summeryMeanYSort$correctY <- summeryMeanYSort$GC_content.x * mean(summeryMeanYSort$GC_content.y) + summeryMeanYSort$`(Intercept)` ###

  # remove samples with no coverage for most of contigs
  summeryMeanYSort2 <- summeryMeanYSort

  ### cutoff of filtering samples changes according to parameter i
  ### i=1, cutoffRatio is 0.5; i=2, cutoffRatio is 1 as contig is clean
  i <- 1
  Samples_filteredY <- filterSample(summeryMeanYSort2, 0, 1 / (3 - i))
  if (length(Samples_filteredY) < 2) {
    return("too few (<2) samples with reliable coverages for the set of contigs")
  }

  # summeryMeanYSortFilterWide is the relatively clean matrix with more confident samples and the corresponding contigs with available values
  summeryMeanYSortFilterWide <- reshapeFiltered(Samples_filteredY, summeryMeanYSort2)

  # do PCA for contigs
  pca <- contig_pca(summeryMeanYSortFilterWide)
  ksResult <- ks.test(pca$PC1, "punif", min(pca$PC1), max(pca$PC1))

  # all good contigs follow uniform distribution
  range <- select_by_ks_test(sort(pca$PC1))

  if (range[3] == TRUE) {
    contigPCAPC1Filtered <- subset(pca, PC1 >= range[1] & PC1 <= range[2])
  } else {
    return("cannot find a continuous set of contigs in uniform distribution")
  }
  # largerClusterContig contains contigs within the range consistent with uniform distribution
  largerClusterContig <- rownames(contigPCAPC1Filtered)

  summeryMeanYSort2 <- subset(summeryMeanYSort, contig %in% largerClusterContig)
  summeryMeanYSortWide <- reshape2::dcast(subset(summeryMeanYSort2, select = c("sample", "contig", "correctY")), contig ~ sample)
  summeryMeanYSortWide2 <- subset(summeryMeanYSortWide, select = -c(contig))
  rownames(summeryMeanYSortWide2) <- summeryMeanYSortWide$contig
  # summeryMeanYSortWide <- reshapeRmNA(summeryMeanYSort2)
  ### debug to skip cor.test when too few values
  summaryMeanCor <- data.frame(Mean = apply(summeryMeanYSortWide2, 2, function(x) mean(x, na.rm = TRUE)), NAprtg = apply(summeryMeanYSortWide2, 2, function(x) length(x[is.na(x)])), cor.p = apply(summeryMeanYSortWide2, 2, function(x) ifelse(length(x[!is.na(x)]) >= 5, cor.test(contigPCAPC1Filtered$PC1, x)$p.value, NA)), KS.p = apply(summeryMeanYSortWide2, 2, function(x) ifelse(length(x[!is.na(x)]) >= 5, ks(x), NA)))
  ### use the cutoff adjust.cor.p
  ### should have few missing contigs(cleaned), cutoff of percentage 0.05 here

  summaryMeanCor$adjust.cor.p <- p.adjust(summaryMeanCor$cor.p, "BH")
  summaryMeanCor$adjust.KS.p <- p.adjust(summaryMeanCor$KS.p, "BH")
  SamplesFilteredFinal <- rownames(summaryMeanCor[summaryMeanCor$NAprtg <= 0.05 * nrow(summeryMeanYSortWide2) & (summaryMeanCor$Mean > 0 | summaryMeanCor$adjust.cor.p < 0.05), ])

  summeryMeanYSortFilteredSampleContig <- subset(summeryMeanYSort2, sample %in% SamplesFilteredFinal, select = c(sample, contig, correctY))

  return(list(SamplesFilteredFinal, summeryMeanYSortFilteredSampleContig, contigPCAPC1Filtered, range, Samples_filteredY))
}

#' A function for iteration of pipeline until convergence
#' @param Z a matrix of coverages
#' @return final filtered samples, matrix of sample, contig and corrected coverages,
#' filtered contigs with PC1 values, PC1 range, preliminary filtered samples
iterate_pipelines <- function(Z) {
  pipeline1 <- pipeline(Z, 1)
  if (length(pipeline1) == 1) {
    return(pipeline1)
  }
  Samples_filtered2 <- pipeline1[[1]]
  summeryMeanYSortFilteredSampleContig2 <- pipeline1[[2]]
  contigPCAPC1Filtered2 <- pipeline1[[3]]
  range2 <- pipeline1[[4]]
  Samples_filteredY2 <- pipeline1[[5]]
  rm(pipeline1)

  ### until convergence
  while ((length(unique(Z$sample)) != length(Samples_filtered2)) | (length(unique(Z$contig)) != length(unique(contigPCAPC1Filtered2$contig)))) {
    Z <- subset(Z, sample %in% Samples_filtered2 & contig %in% contigPCAPC1Filtered2$contig)
    pipeline2 <- pipeline(Z, 1)
    if (length(pipeline2) == 1) {
      return(pipeline2)
    }

    Samples_filtered2 <- pipeline2[[1]]
    summeryMeanYSortFilteredSampleContig2 <- pipeline2[[2]]
    contigPCAPC1Filtered2 <- pipeline2[[3]]
    range2 <- pipeline2[[4]]
  }

  list(Samples_filtered2, summeryMeanYSortFilteredSampleContig2, contigPCAPC1Filtered2, range2, Samples_filteredY2)
}
