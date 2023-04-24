############################################################################
# Copyright (c) 2016-2018 DBEI, UPENN
# All Rights Reserved
# See file LICENSE for details.
############################################################################

library(logger)
library(stringr)

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
#' @importFrom stats coef
#' @importFrom stats cor
#' @importFrom stats cor.test
#' @importFrom stats ks.test
#' @importFrom stats p.adjust
#' @importFrom logger log_info
pipeline <- function(Y, i) {
  logger::log_info("Starting pipeline run...")

  lmeModel <- lme4::lmer(logCov ~ GC + (1 | sample:contig), data = Y, REML = FALSE)
  summeryMeanY <- aggregate(GC ~ (sample:contig), Y, FUN = "mean")
  summeryMeanY$s_c <- paste(summeryMeanY$sample, summeryMeanY$contig, sep = ":")

  lmeModelCoef <- coef(lmeModel)$`sample:contig`
  lmeModelCoef$s_c <- rownames(lmeModelCoef)

  summeryMeanYSort <- merge(lmeModelCoef, summeryMeanY, by = "s_c")
  summeryMeanYSort$correctY <- summeryMeanYSort$GC.x * mean(summeryMeanYSort$GC.y) + summeryMeanYSort$`(Intercept)` ###

  # remove samples with no coverage for most of contigs
  summeryMeanYSort2 <- summeryMeanYSort

  ### cutoff of filtering samples changes according to parameter i
  ### i=1, cutoffRatio is 0.5; i=2, cutoffRatio is 1 as contig is clean
  i <- 1
  Samples_filteredY <- filterSample(summeryMeanYSort2, 0, 1 / (3 - i))
  if (length(Samples_filteredY) < 2) {
    logger::log_info("Too few (<2) samples with reliable coverages for the set of contigs")
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
    logger::log_info("Cannot find a continuous set of contigs in uniform distribution")
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
#'
#' @importFrom logger log_info
itePipelines <- function(Z) {
  logger::log_info("Starting pipeline iteration...")
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

#' Main function
#'
#' @param input path to the coverage matrix file in csv form (column names: "logCov", "GC", "sample", "contig", "length")
#' @param output path to the output directory (default: "getwd()/output/")
#' @param max_candidate_iter max allowed iterations for estimation of PTR (default: 10)
#' @returns nothing, but outputs file with results
#'
#' @importFrom utils read.csv
#' @importFrom utils write.table
#' @importFrom stats aggregate
#' @importFrom stats prcomp
#' @importFrom logger log_info
#' @importFrom logger log_appender
#' @importFrom logger appender_file
#'
# @examples
# estGrowthRate("tests/testthat/data/all_final_contigs.cov3", "tests/testthat/data/output", 10)
# estGrowthRate("tests/testthat/data/all_final_contigs.cov3")
#'
#' @export
estGrowthRate <- function(input, output, max_candidate_iter) {
  stopifnot(file.exists(input))
  if (missing(output)) {
    logger::log_info(stringr::str_glue("Setting output path to {file.path(getwd(), output)}"))
    output <- file.path(getwd(), output)
  }

  if (file.exists(file.path(output, "log"))) {
    file.remove(file.path(output, "log"))
  }
  logger::log_appender(logger::appender_tee(file.path(output, "log")))

  if (!dir.exists(output)) {
    logger::log_info("Creating output dir")
    dir.create(output)
  }
  if (missing(max_candidate_iter)) {
    logger::log_info("Setting max_candidate_iter to default (10)")
    max_candidate_iter <- 10
  }

  logger::log_info("Starting DEMIC...")

  # Load matrix of .cov3 and rename the heads
  X <- read.csv(input, header = FALSE, stringsAsFactors = TRUE)

  colnames(X) <- c("logCov", "GC", "sample", "contig", "length")

  tag_permu <- 0

  Y <- X

  save.image(file.path(output, "savepoint_demic.RData"))
  logger::log_info(stringr::str_glue("Saved image to {file.path(output, 'savepoint_demic.RData')}"))

  # Attempt the default iteration for contigs
  if (length(levels(X$contig)) >= 20 & length(levels(X$sample)) >= 3) {
    logger::log_info("Attempting default iteration for contigs")
    cor_cutoff <- 0.98
    max_cor <- 0
    for (s2 in 1:3) {
      if (s2 == 2) {
        if (max_cor < 0.9) {
          break
        } else if (max_cor < 0.95) {
          cor_cutoff <- 0.95
        }
      }

      nrm <- floor(length(levels(Y$contig)) / 5)
      set.seed(s2)
      designateR2 <- sample.int(10000, size = length(levels(Y$contig)), replace = FALSE)
      ContigDesignateR2 <- data.frame("Contig" = levels(Y$contig), "number" = designateR2)
      ContigDesignateRSort2 <- ContigDesignateR2[order(ContigDesignateR2[, 2]), ]
      nacontig_id <- NULL
      for (x in 1:4) {
        for (y in (x + 1):5) {
          if (x %in% nacontig_id) {
            break
          }
          if (y %in% nacontig_id) {
            next
          }
          ContigDesignateRemove1 <- ContigDesignateRSort2[(nrm * (x - 1) + 1):(nrm * x), 1]
          ContigDesignateRemove2 <- ContigDesignateRSort2[(nrm * (y - 1) + 1):(nrm * y), 1]

          pipelineX1 <- itePipelines(X[!X$contig %in% ContigDesignateRemove1, ])
          if (length(pipelineX1) == 1) {
            break
          }
          pipelineX2 <- itePipelines(X[!X$contig %in% ContigDesignateRemove2, ])
          if (length(pipelineX2) == 1) {
            nacontig_id <- c(nacontig_id, y)
            next
          }
          Samples_filteredXrm1 <- pipelineX1[[1]]
          summeryMeanYSortFilteredSampleContigXrm1 <- pipelineX1[[2]]
          contigPCAPC1FilteredXrm1 <- pipelineX1[[3]]
          rangeXrm1 <- pipelineX1[[4]]

          Samples_filteredXrm2 <- pipelineX2[[1]]
          summeryMeanYSortFilteredSampleContigXrm2 <- pipelineX2[[2]]
          contigPCAPC1FilteredXrm2 <- pipelineX2[[3]]
          rangeXrm2 <- pipelineX2[[4]]
          if (length(contigPCAPC1FilteredXrm1$contig) - length(intersect(contigPCAPC1FilteredXrm1$contig, contigPCAPC1FilteredXrm2$contig)) < 3 | length(contigPCAPC1FilteredXrm2$contig) - length(intersect(contigPCAPC1FilteredXrm1$contig, contigPCAPC1FilteredXrm2$contig)) < 3) {
            next
          }

          SampleCorrectYWithPC1 <- merge(reshape2::dcast(subset(summeryMeanYSortFilteredSampleContigXrm1, select = c("sample", "contig", "correctY")), contig ~ sample), contigPCAPC1FilteredXrm1)

          lmModelCo <- apply(subset(SampleCorrectYWithPC1, select = -c(contig, PC1)), 2, lmColumn, y = SampleCorrectYWithPC1$PC1)
          cor_model <- apply(subset(SampleCorrectYWithPC1, select = -c(contig, PC1)), 2, function(x) cor.test(SampleCorrectYWithPC1$PC1, x)$estimate)

          estPTRs <- data.frame("estPTR" = 2^abs(lmModelCo[1, ] * (rangeXrm1[1] - rangeXrm1[2])), "coefficient" = lmModelCo[1, ], "pValue" = lmModelCo[2, ], "cor" = cor_model)
          estPTRs$sample <- rownames(estPTRs)
          estPTRsEach1 <- merge(estPTRs, aggregate(correctY ~ sample, summeryMeanYSortFilteredSampleContigXrm1, FUN = "median"), by = "sample")

          SampleCorrectYWithPC1 <- merge(reshape2::dcast(subset(summeryMeanYSortFilteredSampleContigXrm2, select = c("sample", "contig", "correctY")), contig ~ sample), contigPCAPC1FilteredXrm2)

          lmModelCo <- apply(subset(SampleCorrectYWithPC1, select = -c(contig, PC1)), 2, lmColumn, y = SampleCorrectYWithPC1$PC1)
          cor_model <- apply(subset(SampleCorrectYWithPC1, select = -c(contig, PC1)), 2, function(x) cor.test(SampleCorrectYWithPC1$PC1, x)$estimate)

          estPTRs <- data.frame("estPTR" = 2^abs(lmModelCo[1, ] * (rangeXrm2[1] - rangeXrm2[2])), "coefficient" = lmModelCo[1, ], "pValue" = lmModelCo[2, ], "cor" = cor_model)
          estPTRs$sample <- rownames(estPTRs)
          estPTRsEach2 <- merge(estPTRs, aggregate(correctY ~ sample, summeryMeanYSortFilteredSampleContigXrm2, FUN = "median"), by = "sample")

          minor_sample1 <- cor_diff(estPTRsEach1)
          minor_sample2 <- cor_diff(estPTRsEach2)
          if ((length(minor_sample1) > 0 & length(minor_sample2) > 0) | (max(estPTRsEach1$estPTR) < 1.8 & max(estPTRsEach2$estPTR) < 1.8) | (max(estPTRsEach1$estPTR) / min(estPTRsEach1$estPTR) > 5 & max(estPTRsEach2$estPTR) / min(estPTRsEach2$estPTR) > 5)) {
            next
          }

          estPTRsEach12 <- merge(estPTRsEach1, estPTRsEach2, by = "sample")

          if (nrow(estPTRsEach12) > 0.9 * nrow(estPTRsEach1) & nrow(estPTRsEach12) > 0.9 * nrow(estPTRsEach2)) {
            cor_current <- cor(estPTRsEach12$estPTR.x, estPTRsEach12$estPTR.y)
            if (cor_current > max_cor) {
              max_cor <- cor_current
            }

            if (cor(estPTRsEach12$estPTR.x, estPTRsEach12$estPTR.y) > cor_cutoff) {
              tag_permu <- 1
              estPTRsEach12$estPTR <- apply(subset(estPTRsEach12, select = c("estPTR.x", "estPTR.y")), 1, mean)
              estPTRsEach12$coefficient <- apply(subset(estPTRsEach12, select = c("coefficient.x", "coefficient.y")), 1, function(x) mean(abs(x)))
              estPTRsEach12$pValue <- apply(subset(estPTRsEach12, select = c("pValue.x", "pValue.y")), 1, max)
              estPTRsEach12$cor <- apply(subset(estPTRsEach12, select = c("cor.x", "cor.y")), 1, function(x) mean(abs(x)))
              estPTRsEach12$correctY <- apply(subset(estPTRsEach12, select = c("correctY.x", "correctY.y")), 1, mean)

              estPTRs2 <- subset(estPTRsEach12, select = c("sample", "estPTR", "coefficient", "pValue", "cor", "correctY"))
              break
            }
          }
        }
        if (tag_permu == 1) {
          break
        }
      }
      if (tag_permu == 1) {
        break
      }
    }
  }

  # Attempt alternative iteration for samples
  if (tag_permu == 0) {
    logger::log_info("Attempting alternative iteration for samples")
    pipelineY <- itePipelines(Y)
    if (length(pipelineY) == 1) {
      stop("pipelineY failed")
    }
    Samples_filtered1 <- pipelineY[[1]]
    summeryMeanYSortFilteredSampleContig1 <- pipelineY[[2]]
    contigPCAPC1Filtered1 <- pipelineY[[3]]
    range1 <- pipelineY[[4]]
    Samples_filteredY1 <- pipelineY[[5]]

    SampleCorrectYWithPC1 <- merge(reshape2::dcast(subset(summeryMeanYSortFilteredSampleContig1, select = c("sample", "contig", "correctY")), contig ~ sample), contigPCAPC1Filtered1)

    lmModelCo <- apply(subset(SampleCorrectYWithPC1, select = -c(contig, PC1)), 2, lmColumn, y = SampleCorrectYWithPC1$PC1)
    cor_model <- apply(subset(SampleCorrectYWithPC1, select = -c(contig, PC1)), 2, function(x) cor.test(SampleCorrectYWithPC1$PC1, x)$estimate)

    estPTRs <- data.frame("estPTR" = 2^abs(lmModelCo[1, ] * (range1[1] - range1[2])), "coefficient" = lmModelCo[1, ], "pValue" = lmModelCo[2, ], "cor" = cor_model)
    estPTRs$sample <- rownames(estPTRs)
    estPTRs3 <- merge(estPTRs, aggregate(correctY ~ sample, summeryMeanYSortFilteredSampleContig1, FUN = "median"), by = "sample")

    minor_sample3 <- cor_diff(estPTRs3)

    if (length(minor_sample3) == 0 & max(estPTRs3$estPTR) >= 1.8 & max(estPTRs3$estPTR) / min(estPTRs3$estPTR) <= 5) {
      estPTRs2 <- estPTRs3
    } else if ((length(minor_sample3) > 0 | max(estPTRs3$estPTR) < 1.8 | max(estPTRs3$estPTR) / min(estPTRs3$estPTR) > 5) & length(Samples_filteredY1) >= 6) {
      estPTRfinal <- NULL
      for (s in 1:max_candidate_iter) {
        set.seed(s)
        designateR <- sample.int(10000, size = length(Samples_filteredY1), replace = FALSE)
        SampleDesignateR <- data.frame("Sample" = Samples_filteredY1, "number" = designateR)
        SampleDesignateRSort <- SampleDesignateR[order(SampleDesignateR[, 2]), ]
        selectSamples <- NULL
        pipelineX <- NULL
        estPTRsEach <- NULL
        Samples_filteredX <- NULL
        for (q in 0:2) {
          selectSamples[[q + 1]] <- SampleDesignateRSort[(1:length(Samples_filteredY1)) %% 3 == q, ]$Sample

          pipelineX <- itePipelines(X[!X$sample %in% selectSamples[[q + 1]], ])
          if (length(pipelineX) == 1) {
            break
          }
          Samples_filteredX[[q + 1]] <- pipelineX[[1]]
          summeryMeanYSortFilteredSampleContigX <- pipelineX[[2]]
          contigPCAPC1FilteredX <- pipelineX[[3]]
          rangeX <- pipelineX[[4]]

          SampleCorrectYWithPC1 <- merge(reshape2::dcast(subset(summeryMeanYSortFilteredSampleContigX, select = c("sample", "contig", "correctY")), contig ~ sample), contigPCAPC1FilteredX)

          lmModelCo <- apply(subset(SampleCorrectYWithPC1, select = -c(contig, PC1)), 2, lmColumn, y = SampleCorrectYWithPC1$PC1)
          cor_model <- apply(subset(SampleCorrectYWithPC1, select = -c(contig, PC1)), 2, function(x) cor.test(SampleCorrectYWithPC1$PC1, x)$estimate)

          estPTRs <- data.frame("estPTR" = 2^abs(lmModelCo[1, ] * (rangeX[1] - rangeX[2])), "coefficient" = lmModelCo[1, ], "pValue" = lmModelCo[2, ], "cor" = cor_model)
          estPTRs$sample <- rownames(estPTRs)
          estPTRsEach[[q + 1]] <- merge(estPTRs, aggregate(correctY ~ sample, summeryMeanYSortFilteredSampleContig1, FUN = "median"), by = "sample")
        }
        if (length(estPTRsEach) < 3) {
          next
        }
        qmax <- NULL
        rmax <- NULL
        estPTRsIntBest <- NULL
        cormax <- 0
        for (q in 0:1) {
          for (r in (q + 1):2) {
            estPTRsq <- estPTRsEach[[q + 1]]
            estPTRsr <- estPTRsEach[[r + 1]]
            intSamples <- intersect(setdiff(Samples_filteredY1, selectSamples[[q + 1]]), setdiff(Samples_filteredY1, selectSamples[[r + 1]]))
            estPTRsInt <- merge(estPTRsq[estPTRsq$sample %in% intSamples, c("sample", "estPTR")], estPTRsr[estPTRsr$sample %in% intSamples, c("sample", "estPTR")], by = "sample")
            minor_sample_q <- cor_diff(estPTRsq)
            minor_sample_r <- cor_diff(estPTRsr)

            corqr <- cor(estPTRsInt[, 3], estPTRsInt[, 2])
            if (corqr > cormax & length(minor_sample_q) == 0 & length(minor_sample_r) == 0) {
              cormax <- corqr
              qmax <- estPTRsq
              rmax <- estPTRsr
              estPTRsIntBest <- estPTRsInt
            }
          }
        }
        if (cormax > 0.98) {
          rownames(estPTRsIntBest) <- estPTRsIntBest$sample
          estPTRsInt <- subset(estPTRsIntBest, select = -c(sample))

          estPTRsIntPCA <- prcomp(estPTRsInt)

          qmax$TestPTR <- (qmax$estPTR - mean(estPTRsInt$estPTR.x)) / estPTRsIntPCA$rotation[1, 1]
          rmax$TestPTR <- (rmax$estPTR - mean(estPTRsInt$estPTR.y)) / estPTRsIntPCA$rotation[2, 1]
          rmax$TestPTR2 <- rmax$TestPTR * estPTRsIntPCA$rotation[1, 1] + mean(estPTRsInt$estPTR.x)
          qmax$TestPTR2 <- qmax$TestPTR * estPTRsIntPCA$rotation[2, 1] + mean(estPTRsInt$estPTR.y)

          if (testReasonable(qmax$TestPTR2, rmax$estPTR) > testReasonable(rmax$TestPTR2, qmax$estPTR) & testReasonable(qmax$TestPTR2, rmax$estPTR) > 0.2) {
            estPTRfinal <- dfTransfer(rmax, qmax)
            break
          } else if (testReasonable(qmax$TestPTR2, rmax$estPTR) < testReasonable(rmax$TestPTR2, qmax$estPTR) & testReasonable(rmax$TestPTR2, qmax$estPTR) > 0.2) {
            estPTRfinal <- dfTransfer(qmax, rmax)
            break
          } else {
            next
          }
        } else {
          next
        }
      }
      if (nrow(estPTRfinal) < length(Samples_filteredY1)) {
        stop("fail to calculate consistent PTRs or combine PTRs from two subsets of samples!")
      } else {
        estPTRs2 <- estPTRfinal
      }
    } else {
      stop("fail to calculate consistent PTRs!")
    }
  }

  # Output to .eptr
  logger::log_info(stringr::str_glue("Writing output to {file.path(output, 'out.eptr')}"))
  final_output <- file.path(output, "out.eptr")
  write.table(estPTRs2, file = final_output, sep = "\t", quote = FALSE)

  # Edit by scottdaniel25@gmail.com
  save.image(file.path(output, "finished_demic.RData"))
}
