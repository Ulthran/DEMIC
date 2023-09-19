source("R/utils.R")

#' Attempt the default iteration for contigs
#' Requires at least 20 contigs
#'
#' @param X cov3 matrix as a dataframe
#' @return estPTRs matrix on success, NULL otherwise
#'
#' @importFrom stats aggregate prcomp var
contigsPipeline <- function(X) {
  if (length(levels(X$contig)) < 20 | length(levels(X$sample)) < 3) {
    return(NULL)
  }
  cor_cutoff <- 0.98
  max_cor <- 0
  for (s2 in 1:3) {
    if (s2 == 2) {
      if (max_cor < 0.9) {
        return(NULL)
      } else if (max_cor < 0.95) {
        cor_cutoff <- 0.95
      }
    }

    nrm <- floor(length(levels(X$contig)) / 5)
    set.seed(s2)
    designateR2 <- sample.int(10000, size = length(levels(X$contig)), replace = FALSE)
    ContigDesignateR2 <- data.frame("Contig" = levels(X$contig), "number" = designateR2)
    ContigDesignateRSort2 <- ContigDesignateR2[order(ContigDesignateR2[, 2]), ]
    nacontig_id <- NULL
    for (x in 1:4) {
      for (y in (x + 1):5) {
        if (x %in% nacontig_id) {
          return(NULL)
        }
        if (y %in% nacontig_id) {
          next
        }
        ContigDesignateRemove1 <- ContigDesignateRSort2[(nrm * (x - 1) + 1):(nrm * x), 1]
        ContigDesignateRemove2 <- ContigDesignateRSort2[(nrm * (y - 1) + 1):(nrm * y), 1]

        pipelineX1 <- itePipelines(X[!X$contig %in% ContigDesignateRemove1, ])
        if (length(pipelineX1) == 1) {
          return(NULL)
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

        if (nrow(estPTRsEach12) > 0.9 * max(c(nrow(estPTRsEach1), nrow(estPTRsEach2)))) {
          if (var(estPTRsEach12$estPTR.x) == 0 || var(estPTRsEach12$estPTR.y) == 0) {
            estPTRsEach12$estPTR.x <- jitter(estPTRsEach12$estPTR.x)
            estPTRsEach12$estPTR.y <- jitter(estPTRsEach12$estPTR.y)
          }
          cor_current <- cor(estPTRsEach12$estPTR.x, estPTRsEach12$estPTR.y)
          if (is.na(cor_current)) {
            return(NULL)
          }
          if (cor_current > max_cor) {
            max_cor <- cor_current
          }

          if (cor(estPTRsEach12$estPTR.x, estPTRsEach12$estPTR.y) > cor_cutoff) {
            dput(estPTRsEach12)
            estPTRsEach12$estPTR <- apply(subset(estPTRsEach12, select = c("estPTR.x", "estPTR.y")), 1, mean)
            estPTRsEach12$coefficient <- apply(subset(estPTRsEach12, select = c("coefficient.x", "coefficient.y")), 1, function(x) mean(abs(x)))
            estPTRsEach12$pValue <- apply(subset(estPTRsEach12, select = c("pValue.x", "pValue.y")), 1, max)
            estPTRsEach12$cor <- apply(subset(estPTRsEach12, select = c("cor.x", "cor.y")), 1, function(x) mean(abs(x)))
            estPTRsEach12$correctY <- apply(subset(estPTRsEach12, select = c("correctY.x", "correctY.y")), 1, mean)

            estPTRs2 <- subset(estPTRsEach12, select = c("sample", "estPTR", "coefficient", "pValue", "cor", "correctY"))
            return(estPTRs2)
          }
        }
      }
    }
  }
}

#' Attempt alternative iteration for samples
#' Requires at least 3 samples
#'
#' @param X cov3 matrix as a dataframe
#' @param max_candidate_iter max number of tries for samples pipeline iteration
#' @return estPTRs matrix on success, NULL otherwise
#'
#' @importFrom stats prcomp aggregate
samplesPipeline <- function(X, max_candidate_iter) {
  pipelineY <- itePipelines(X)
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

#' Main function
#'
#' @param X dataframe with coverage matrix (column names: "logCov", "GC", "sample", "contig", "length")
#' @param max_candidate_iter max allowed iterations for estimation of PTR (default: 10)
#' @returns dataframe with the estimated PTRs (column names: "estPTR", "coefficient", "pValue", "cor", "correctY")
#'
#' @export
estPTR <- function(X, max_candidate_iter = 10) {
  estPTRs <- suppressWarnings({
    contigsPipeline(X)
  })
  if (is.null(estPTRs)) {
    Y <- X
    estPTRs <- samplesPipeline(Y)
  }

  estPTRs
}
