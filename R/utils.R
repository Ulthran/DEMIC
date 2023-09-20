############################################################################
# Copyright (c) 2016-2018 DBEI, UPENN
# All Rights Reserved
# See file LICENSE for details.
############################################################################


#' A convenient function for KS test of uniform distribution
#' @param x a vector without NA
#' @return the p value of KS test
ks <- function(x) {
  ks_result <- ks.test(x, "punif", min(x, na.rm = TRUE), max(x, na.rm = TRUE))

  ks_result$p.value
}

#' A function to remove outlier contigs using KS test
#' @param sort_values a vector of sorted values
#' @return a vector with all values following a uniform distribution
select_by_ks_test <- function(sort_values) {
  len <- length(sort_values)
  if (len < 10) {
    return(c(0, 0, FALSE))
  }

  ks_result <- ks(sort_values)

  if (sort_values[2] - sort_values[1] > sort_values[len] - sort_values[len - 1]) {
    if (ks_result < 0.05) {
      return(select_by_ks_test(sort_values[2:len]))
    } else {
      ks_next <- ks(sort_values[2:len])
      if (ks_next > ks_result + 0.1) {
        return(select_by_ks_test(sort_values[2:len]))
      } else {
        return(c(sort_values[1], sort_values[len], TRUE))
      }
    }
  } else {
    if (ks_result < 0.05) {
      return(select_by_ks_test(sort_values[1:len - 1]))
    } else {
      ks_next <- ks(sort_values[2:len])
      if (ks_next > ks_result + 0.1) {
        return(select_by_ks_test(sort_values[1:len - 1]))
      } else {
        return(c(sort_values[1], sort_values[len], TRUE))
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
