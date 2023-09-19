############################################################################
# Copyright (c) 2016-2018 DBEI, UPENN
# All Rights Reserved
# See file LICENSE for details.
############################################################################


#' A convenient function for KS test of uniform distribution
#' @param x a vector without NA
#' @return the p value of KS test
ks <- function(x) {
  ks.result <- ks.test(x, "punif", min(x, na.rm = TRUE), max(x, na.rm = TRUE))
  return(ks.result$p.value)
}

#' A function to remove outlier contigs using KS test
#' @param sortValues a vector of sorted values
#' @return a vector with all values following a uniform distribution
selectAccordingToKSTest <- function(sortValues) {
  len <- length(sortValues)
  if (len < 10) {
    return(c(0, 0, FALSE))
  }

  ksResult <- ks(sortValues)

  if (sortValues[2] - sortValues[1] > sortValues[len] - sortValues[len - 1]) {
    if (ksResult < 0.05) {
      return(selectAccordingToKSTest(sortValues[2:len]))
    } else {
      ks_next <- ks(sortValues[2:len])
      if (ks_next > ksResult + 0.1) {
        return(selectAccordingToKSTest(sortValues[2:len]))
      } else {
        return(c(sortValues[1], sortValues[len], TRUE))
      }
    }
  } else {
    if (ksResult < 0.05) {
      return(selectAccordingToKSTest(sortValues[1:len - 1]))
    } else {
      ks_next <- ks(sortValues[2:len])
      if (ks_next > ksResult + 0.1) {
        return(selectAccordingToKSTest(sortValues[1:len - 1]))
      } else {
        return(c(sortValues[1], sortValues[len], TRUE))
      }
    }
  }
}

#' A convenient function for ordinary linear regression on two vectors
#' @param x first vector
#' @param y second vector
#' @return the coefficient and p value of linear regression
#'
#' @importFrom stats anova
#' @importFrom stats lm
lmColumn <- function(x, y) {
  lmModel <- lm(x ~ y)
  anova_model <- anova(lmModel)
  return(c(lmModel$coefficients[2], anova_model$`Pr(>F)`))
}

#' A helper function for getting the sum and length of a vector
#' @param x the vector
#' @return the vector containing the sum and length of the input
countNA <- function(x) c(sum(x), length(x))

#' A function for sample filtration
#' @param Z a matrix
#' @param cutoffAve threshold of average
#' @param cutoffRatio threshold of ratio
#' @return the coefficient and p value of linear regression
#' Input requirements: 1. have values in more than half of the contigs 2. average log2(cov) > 0 in all these contigs
filterSample <- function(Z, cutoffAve, cutoffRatio) {
  X_summary <- aggregate(correctY ~ (sample), Z, FUN = countNA)

  levelContigX <- length(unique(Z$contig))
  Samples_filtered <- X_summary[X_summary$correctY[, 1] >= cutoffAve * levelContigX & X_summary$correctY[, 2] >= cutoffRatio * levelContigX, ]$sample
  return(Samples_filtered)
}

#' A function for orientation determination
#' @param Z a vector of values
#' @return a subset, where each value has the same majority of orientation
cor_diff <- function(Z) {
  sampleSet1 <- Z[Z$cor > 0, ]$sample
  sampleSet2 <- Z[Z$cor < 0, ]$sample
  if (length(sampleSet1) > length(sampleSet2)) {
    return(sampleSet2)
  } else {
    return(sampleSet1)
  }
}

#' A function for reshape to facilitate PCA, removing all contigs with missing values for designated samples
#' @param Samples_filtered a vector of samples
#' @param Z a matrix of coverage
#' @return a reshaped matrix of coverage
reshapeFiltered <- function(Samples_filtered, Z) {
  X_filter <- subset(Z, sample %in% Samples_filtered, select = c(sample, contig, correctY))
  X_filter_wide <- reshape2::dcast(subset(X_filter, select = c("sample", "contig", "correctY")), contig ~ sample)

  X_filter_wide3 <- apply(subset(X_filter_wide, select = -c(contig)), 1, function(x) length(x[is.na(x)]))
  names(X_filter_wide3) <- X_filter_wide$contig

  X_filter_wide2 <- subset(subset(X_filter_wide, contig %in% names(X_filter_wide3[X_filter_wide3 == 0])), select = -c(contig))

  row.names(X_filter_wide2) <- names(X_filter_wide3[X_filter_wide3 == 0])
  return(X_filter_wide2)
}

#' A function for data frame integration
#' @param x first data frame
#' @param y second data frame
#' @param i 'sample' column
#' @return a data frame with the other column as mean or max of that in the original two
consistTransfer <- function(x, y, i) {
  z <- merge(x, y, by = "sample", all.x = TRUE, all.y = TRUE)
  if (i == 1) {
    estPTRBoth <- apply(subset(z, select = 2:3), 1, function(x) mean(abs(x), na.rm = TRUE))
  } else if (i == 2) {
    estPTRBoth <- apply(subset(z, select = 2:3), 1, function(x) max(abs(x), na.rm = TRUE))
  }
  names(estPTRBoth) <- z$sample
  return(estPTRBoth)
}

#' A function for data frame transfer
#' @param X1 first data frame with six columns
#' @param X2 second data frame with six columns
#' @return a data frame with the same six columns but integrated info
dfTransfer <- function(X1, X2) {
  X12 <- data.frame(
    "sample" = sort(union(X1$sample, X2$sample), method = "shell"),
    "estPTR" = consistTransfer(subset(X1, select = c(sample, estPTR)), subset(X2, select = c(sample, TestPTR2)), 1),
    "coefficient" = consistTransfer(subset(X1, select = c(sample, coefficient)), subset(X2, select = c(sample, coefficient)), 1),
    "pValue" = consistTransfer(subset(X1, select = c(sample, pValue)), subset(X2, select = c(sample, pValue)), 2),
    "cor" = consistTransfer(subset(X1, select = c(sample, cor)), subset(X2, select = c(sample, cor)), 1),
    "correctY" = consistTransfer(subset(X1, select = c(sample, correctY)), subset(X2, select = c(sample, correctY)), 1)
  )
  return(X12)
}

#' A function to test whether the result is reasonable
#' @param a first vector of values
#' @param b second vector of values
#' @return the test result
testReasonable <- function(a, b) {
  c <- c(a, b)
  if (min(c) >= 1) {
    return(min(c) / max(c))
  } else {
    return(0 - length(c[c < 1]))
  }
}

#' A function to return the first dimension of PCA on an input matrix
#' @param X a matrix to undergo PCA
#' @return first dimension of the PCA results
#'
#' @importFrom stats prcomp
contig_pca <- function(X) {
  contigPCA <- prcomp(X) # take first component (PC1)

  return(data.frame("contig" = rownames(contigPCA$x), "PC1" = contigPCA$x[, "PC1"]))
}

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
  range <- selectAccordingToKSTest(sort(pca$PC1))

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
itePipelines <- function(Z) {
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
  return(list(Samples_filtered2, summeryMeanYSortFilteredSampleContig2, contigPCAPC1Filtered2, range2, Samples_filteredY2))
}
