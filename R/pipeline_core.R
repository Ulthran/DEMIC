#' Run mixed linear model with random effect using lme4
#'
#' @param X input data frame
#' @return the model coefficients in a dataframe
#' \itemize{
#'   \item (Intercept): intercept
#'   \item GC_content: coefficient for GC content
#'   \item s_c: sample:contig
#' }
#'
#' @importFrom lme4 lmer
#' @importFrom stats coef
lme4_model <- function(X) {
  lmeModel <- lmer(log_cov ~ GC_content + (1 | sample:contig), data = X, REML = FALSE)

  lmeModelCoef <- coef(lmeModel)$`sample:contig`
  lmeModelCoef$s_c <- rownames(lmeModelCoef)

  lmeModelCoef
}

#' Corrects for the GC bias in coverage calculations
#'
#' @param X input data frame
#' @param lmeModelCoef output from lme4_model
#' @return a dataframe with corrections
#' \itemize{
#'   \item s_c: sample:contig
#'   \item (Intercept): intercept
#'   \item GC_content.x: LME's GC content
#'   \item sample: sample
#'   \item contig: contig
#'   \item GC_content.y: mean GC content
#'   \item correctY: corrected coverage
#' }
#'
#' @importFrom stats aggregate
correct_GC_bias <- function(X, lmeModelCoef) {
  s_c <- aggregate(GC_content ~ (sample:contig), X, FUN = "mean")
  s_c$s_c <- paste(s_c$sample, s_c$contig, sep = ":")
  s_c <- merge(lmeModelCoef, s_c, by = "s_c")
  s_c$correctY <- s_c$GC_content.x * mean(s_c$GC_content.y) + s_c$`(Intercept)`

  s_c
}

pipeline_hq <- function(Y) {
  # Assuming:
  #   A few, large contigs
  #   All from the same species (do filtering beforehand to map contigs back to ref genome if metagenomic)
  #
  # Correct for sequencing bias (LMM) (will this work with contigs the size of the genome? Hopefully yes, it should be per window. Wait no but then the actual correction is averaged over the whole contig...)
  # Infer relative distances of contigs from replication origin
  # Fit contigs to lin reg per sample and use eptr formula




  PC1 <- contig <- correctY <- NULL

  lmeModelCoef <- lme4_model(Y) # Get adjustments to cov per sample/contig (s_c)

  summeryMeanY <- aggregate(GC_content ~ (sample:contig), Y, FUN = "mean") # Avg GC_content per s_c
  summeryMeanY$s_c <- paste(summeryMeanY$sample, summeryMeanY$contig, sep = ":") # add s_c column

  summeryMeanYSort <- merge(lmeModelCoef, summeryMeanY, by = "s_c") # merge dataframes
  summeryMeanYSort$correctY <- summeryMeanYSort$GC_content.x * mean(summeryMeanYSort$GC_content.y) + summeryMeanYSort$`(Intercept)` # Adjustment (intercept) plus average times adjustment (slope)

  # remove samples with no coverage for most of contigs
  summeryMeanYSort2 <- summeryMeanYSort

  ### cutoff of filtering samples changes according to parameter i
  ### i=1, cutoffRatio is 0.5; i=2, cutoffRatio is 1 as contig is clean
  i <- 1
  Samples_filteredY <- filter_sample(summeryMeanYSort2, 0, 1 / (3 - i))
  #browser()
  if (length(Samples_filteredY) < 2) {
    return("too few (<2) samples with reliable coverages for the set of contigs")
  }

  # summeryMeanYSortFilterWide is the relatively clean matrix with more confident samples and the corresponding contigs with available values
  summeryMeanYSortFilterWide <- reshape_filtered(Samples_filteredY, summeryMeanYSort2)
  summeryMeanYSortFilterWide <- reshape_filtered(Y$sample, summeryMeanYSort2)
  browser

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
  summeryMeanYSortWide <- reshape2::dcast(subset(summeryMeanYSort2, select = c("sample", "contig", "correctY")), contig ~ sample, value.var = "correctY")
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

  return(list(samples = SamplesFilteredFinal, correct_ys = summeryMeanYSortFilteredSampleContig, pc1 = contigPCAPC1Filtered, pc1_range = range, samples_y = Samples_filteredY))
}

#' A function representing the pipeline of four steps:
#'   GC bias correction
#'   sample filtration
#'   PCA adjustment
#'   contig filtration
#'
#' @param Y a matrix of coverages
#' @param i cutoff of filtering samples changes according to parameter i; i=1, cutoffRatio is 0.5; i=2, cutoffRatio is 1 as contig is clean
#' @return a named list
#' \itemize{
#'   \item samples: final list of filtered samples
#'   \item correct_ys: dataframe with correct Y values per contig/sample
#'   \item pc1: PC1 results of PCA per contig
#'   \item pc1_range: range of PC1
#'   \item samples_y: samples filtered for reliable coverage
#' }
#'
#' @importFrom stats coef cor cor.test ks.test p.adjust aggregate
pipeline <- function(Y, i = 1) {
  PC1 <- contig <- correctY <- NULL

  # Correct for GC bias
  lmeModelCoef <- lme4_model(Y)
  corrected_covs <- correct_GC_bias(Y, lmeModelCoef)

  # Filter samples with too few contigs with sufficient coverage
  filtered_samples <- filter_samples(corrected_covs, 0, 1/ (3 - i))
  if (length(filtered_samples) < 2) {
    return("too few (<2) samples with reliable coverages for the set of contigs")
  }
  filtered_samples_vals <- reshape_filtered(filtered_samples, corrected_covs)

  # Filter contigs by PCA (assume a uniform distribution over the genome and remove duplicates)
  PC1 <- contig_pca(filtered_samples_vals)
  # ksResult <- ks.test(pca$PC1, "punif", min(pca$PC1), max(pca$PC1)) # WHY IS THIS NOT USED???
  range <- select_by_ks_test(sort(PC1$PC1))

  if (range[3] == TRUE) {
    filtered_contigs <- subset(PC1, PC1 >= range[1] & PC1 <= range[2])
  } else {
    return("cannot find a continuous set of contigs in uniform distribution")
  }

  filtered_samples <- subset(filtered_samples, contig %in% rownames(filtered_contigs))
  filtered_corrected <- reshape2::dcast(subset(filtered_samples, select = c("sample", "contig", "correctY")), contig ~ sample, value.var = "correctY")


  # largerClusterContig contains contigs within the range consistent with uniform distribution
  largerClusterContig <- rownames(contigPCAPC1Filtered)

  summeryMeanYSort2 <- subset(summeryMeanYSort, contig %in% largerClusterContig)
  summeryMeanYSortWide <- reshape2::dcast(subset(summeryMeanYSort2, select = c("sample", "contig", "correctY")), contig ~ sample, value.var = "correctY")
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

  browser()
  return(list(samples = SamplesFilteredFinal, correct_ys = summeryMeanYSortFilteredSampleContig, pc1 = contigPCAPC1Filtered, pc1_range = range, samples_y = Samples_filteredY))
}

#' A function for iteration of pipeline until convergence
#' @param Z a matrix of coverages
#' @return a named list
#' \itemize{
#'   \item samples: vector of final filtered samples
#'   \item correct_ys: matrix of sample, contig and corrected coverages
#'   \item pc1: matrix of contig and PC1 values
#'   \item pc1_range: vector of PC1 range
#'   \item samples_y: samples filtered for reliable coverage
#' }
iterate_pipelines <- function(Z) {
  contig <- NULL
  i <- 1

  repeat {
    pipeline <- pipeline(Z, i)
    i <- i + 1
    if (length(pipeline) == 1) {
      return(pipeline)
    }

    ### until convergence
    if ((length(unique(Z$sample)) == length(pipeline$samples)) && (length(unique(Z$contig)) == length(unique(pipeline$pc1$contig)))) {
      return(list(samples = pipeline$samples, correct_ys = pipeline$correct_ys, pc1 = pipeline$pc1, pc1_range = pipeline$pc1_range, samples_y = pipeline$samples_y))
    } else {
      Z <- subset(Z, sample %in% pipeline$samples & contig %in% pipeline$pc1$contig)
    }
  }
}
