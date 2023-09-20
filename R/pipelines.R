#' Attempt the default iteration for contigs
#' Requires at least 20 contigs
#'
#' @param X cov3 matrix as a dataframe
#' @return est_PTRs matrix on success, NULL otherwise
#'
#' @importFrom stats aggregate prcomp var
contigs_pipeline <- function(X) {
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
    designate_R <- sample.int(10000, size = length(levels(X$contig)), replace = FALSE)
    contig_designate_R <- data.frame("Contig" = levels(X$contig), "number" = designate_R)
    contig_designate_R_sort <- contig_designate_R[order(contig_designate_R[, 2]), ]
    nacontig_id <- NULL
    for (x in 1:4) {
      for (y in (x + 1):5) {
        if (x %in% nacontig_id) {
          return(NULL)
        }
        if (y %in% nacontig_id) {
          next
        }
        contig_designate_remove_x <- contig_designate_R_sort[(nrm * (x - 1) + 1):(nrm * x), 1]
        contig_designate_remove_y <- contig_designate_R_sort[(nrm * (y - 1) + 1):(nrm * y), 1]

        pipeline_x <- iterate_pipelines(X[!X$contig %in% contig_designate_remove_x, ])
        if (length(pipeline_x) == 1) {
          return(NULL)
        }
        pipeline_y <- iterate_pipelines(X[!X$contig %in% contig_designate_remove_y, ])
        if (length(pipeline_y) == 1) {
          nacontig_id <- c(nacontig_id, y)
          next
        }
        samples_filtered_x <- pipeline_x[[1]]
        sum_mean_sort_sample_contig_x <- pipeline_x[[2]]
        contig_PC1_filtered_x <- pipeline_x[[3]]
        range_x <- pipeline_x[[4]]

        samples_filtered_y <- pipeline_y[[1]]
        sum_mean_sort_sample_contig_y <- pipeline_y[[2]]
        contig_PC1_filtered_y <- pipeline_y[[3]]
        range_y <- pipeline_y[[4]]
        if (length(contig_PC1_filtered_x$contig) - length(intersect(contig_PC1_filtered_x$contig, contig_PC1_filtered_y$contig)) < 3 | length(contig_PC1_filtered_y$contig) - length(intersect(contig_PC1_filtered_x$contig, contig_PC1_filtered_y$contig)) < 3) {
          next
        }

        sample_correct_y_PC1 <- merge(reshape2::dcast(subset(sum_mean_sort_sample_contig_x, select = c("sample", "contig", "correctY")), contig ~ sample), contig_PC1_filtered_x)

        lm_model_co <- apply(subset(sample_correct_y_PC1, select = -c(contig, PC1)), 2, lm_column, y = sample_correct_y_PC1$PC1)
        cor_model <- apply(subset(sample_correct_y_PC1, select = -c(contig, PC1)), 2, function(x) cor.test(sample_correct_y_PC1$PC1, x)$estimate)

        est_PTRs <- data.frame("est_PTR" = 2^abs(lm_model_co[1, ] * (range_x[1] - range_x[2])), "coefficient" = lm_model_co[1, ], "pValue" = lm_model_co[2, ], "cor" = cor_model)
        est_PTRs$sample <- rownames(est_PTRs)
        est_PTRsEach1 <- merge(est_PTRs, aggregate(correctY ~ sample, sum_mean_sort_sample_contig_x, FUN = "median"), by = "sample")

        sample_correct_y_PC1 <- merge(reshape2::dcast(subset(sum_mean_sort_sample_contig_y, select = c("sample", "contig", "correctY")), contig ~ sample), contig_PC1_filtered_y)

        lm_model_co <- apply(subset(sample_correct_y_PC1, select = -c(contig, PC1)), 2, lm_column, y = sample_correct_y_PC1$PC1)
        cor_model <- apply(subset(sample_correct_y_PC1, select = -c(contig, PC1)), 2, function(x) cor.test(sample_correct_y_PC1$PC1, x)$estimate)

        est_PTRs <- data.frame("est_PTR" = 2^abs(lm_model_co[1, ] * (range_y[1] - range_y[2])), "coefficient" = lm_model_co[1, ], "pValue" = lm_model_co[2, ], "cor" = cor_model)
        est_PTRs$sample <- rownames(est_PTRs)
        est_PTRsEach2 <- merge(est_PTRs, aggregate(correctY ~ sample, sum_mean_sort_sample_contig_y, FUN = "median"), by = "sample")

        minor_sample1 <- cor_diff(est_PTRsEach1)
        minor_sample2 <- cor_diff(est_PTRsEach2)
        if ((length(minor_sample1) > 0 & length(minor_sample2) > 0) | (max(est_PTRsEach1$est_PTR) < 1.8 & max(est_PTRsEach2$est_PTR) < 1.8) | (max(est_PTRsEach1$est_PTR) / min(est_PTRsEach1$est_PTR) > 5 & max(est_PTRsEach2$est_PTR) / min(est_PTRsEach2$est_PTR) > 5)) {
          next
        }

        est_PTRsEach12 <- merge(est_PTRsEach1, est_PTRsEach2, by = "sample")

        if (nrow(est_PTRsEach12) > 0.9 * max(c(nrow(est_PTRsEach1), nrow(est_PTRsEach2)))) {
          if (var(est_PTRsEach12$est_PTR.x) == 0 || var(est_PTRsEach12$est_PTR.y) == 0) {
            est_PTRsEach12$est_PTR.x <- jitter(est_PTRsEach12$est_PTR.x)
            est_PTRsEach12$est_PTR.y <- jitter(est_PTRsEach12$est_PTR.y)
          }
          cor_current <- cor(est_PTRsEach12$est_PTR.x, est_PTRsEach12$est_PTR.y)
          if (is.na(cor_current)) {
            return(NULL)
          }
          if (cor_current > max_cor) {
            max_cor <- cor_current
          }

          if (cor(est_PTRsEach12$est_PTR.x, est_PTRsEach12$est_PTR.y) > cor_cutoff) {
            dput(est_PTRsEach12)
            est_PTRsEach12$est_PTR <- apply(subset(est_PTRsEach12, select = c("est_PTR.x", "est_PTR.y")), 1, mean)
            est_PTRsEach12$coefficient <- apply(subset(est_PTRsEach12, select = c("coefficient.x", "coefficient.y")), 1, function(x) mean(abs(x)))
            est_PTRsEach12$pValue <- apply(subset(est_PTRsEach12, select = c("pValue.x", "pValue.y")), 1, max)
            est_PTRsEach12$cor <- apply(subset(est_PTRsEach12, select = c("cor.x", "cor.y")), 1, function(x) mean(abs(x)))
            est_PTRsEach12$correctY <- apply(subset(est_PTRsEach12, select = c("correctY.x", "correctY.y")), 1, mean)

            est_PTRs2 <- subset(est_PTRsEach12, select = c("sample", "est_PTR", "coefficient", "pValue", "cor", "correctY"))
            return(est_PTRs2)
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
#' @return est_PTRs matrix on success, NULL otherwise
#'
#' @importFrom stats prcomp aggregate
samples_pipeline <- function(X, max_candidate_iter) {
  pipelineY <- iterate_pipelines(X)
  if (length(pipelineY) == 1) {
    stop("pipelineY failed")
  }
  Samples_filtered1 <- pipelineY[[1]]
  summeryMeanYSortFilteredSampleContig1 <- pipelineY[[2]]
  contigPCAPC1Filtered1 <- pipelineY[[3]]
  range1 <- pipelineY[[4]]
  Samples_filteredY1 <- pipelineY[[5]]

  sample_correct_y_PC1 <- merge(reshape2::dcast(subset(summeryMeanYSortFilteredSampleContig1, select = c("sample", "contig", "correctY")), contig ~ sample), contigPCAPC1Filtered1)

  lm_model_co <- apply(subset(sample_correct_y_PC1, select = -c(contig, PC1)), 2, lm_column, y = sample_correct_y_PC1$PC1)
  cor_model <- apply(subset(sample_correct_y_PC1, select = -c(contig, PC1)), 2, function(x) cor.test(sample_correct_y_PC1$PC1, x)$estimate)

  est_PTRs <- data.frame("est_PTR" = 2^abs(lm_model_co[1, ] * (range1[1] - range1[2])), "coefficient" = lm_model_co[1, ], "pValue" = lm_model_co[2, ], "cor" = cor_model)
  est_PTRs$sample <- rownames(est_PTRs)
  est_PTRs3 <- merge(est_PTRs, aggregate(correctY ~ sample, summeryMeanYSortFilteredSampleContig1, FUN = "median"), by = "sample")

  minor_sample3 <- cor_diff(est_PTRs3)

  if (length(minor_sample3) == 0 & max(est_PTRs3$est_PTR) >= 1.8 & max(est_PTRs3$est_PTR) / min(est_PTRs3$est_PTR) <= 5) {
    est_PTRs2 <- est_PTRs3
  } else if ((length(minor_sample3) > 0 | max(est_PTRs3$est_PTR) < 1.8 | max(est_PTRs3$est_PTR) / min(est_PTRs3$est_PTR) > 5) & length(Samples_filteredY1) >= 6) {
    est_PTRfinal <- NULL
    for (s in 1:max_candidate_iter) {
      set.seed(s)
      designateR <- sample.int(10000, size = length(Samples_filteredY1), replace = FALSE)
      SampleDesignateR <- data.frame("Sample" = Samples_filteredY1, "number" = designateR)
      SampleDesignateRSort <- SampleDesignateR[order(SampleDesignateR[, 2]), ]
      selectSamples <- NULL
      pipelineX <- NULL
      est_PTRsEach <- NULL
      Samples_filteredX <- NULL
      for (q in 0:2) {
        selectSamples[[q + 1]] <- SampleDesignateRSort[(1:length(Samples_filteredY1)) %% 3 == q, ]$Sample

        pipelineX <- iterate_pipelines(X[!X$sample %in% selectSamples[[q + 1]], ])
        if (length(pipelineX) == 1) {
          break
        }
        Samples_filteredX[[q + 1]] <- pipelineX[[1]]
        summeryMeanYSortFilteredSampleContigX <- pipelineX[[2]]
        contigPCAPC1FilteredX <- pipelineX[[3]]
        rangeX <- pipelineX[[4]]

        sample_correct_y_PC1 <- merge(reshape2::dcast(subset(summeryMeanYSortFilteredSampleContigX, select = c("sample", "contig", "correctY")), contig ~ sample), contigPCAPC1FilteredX)

        lm_model_co <- apply(subset(sample_correct_y_PC1, select = -c(contig, PC1)), 2, lm_column, y = sample_correct_y_PC1$PC1)
        cor_model <- apply(subset(sample_correct_y_PC1, select = -c(contig, PC1)), 2, function(x) cor.test(sample_correct_y_PC1$PC1, x)$estimate)

        est_PTRs <- data.frame("est_PTR" = 2^abs(lm_model_co[1, ] * (rangeX[1] - rangeX[2])), "coefficient" = lm_model_co[1, ], "pValue" = lm_model_co[2, ], "cor" = cor_model)
        est_PTRs$sample <- rownames(est_PTRs)
        est_PTRsEach[[q + 1]] <- merge(est_PTRs, aggregate(correctY ~ sample, summeryMeanYSortFilteredSampleContig1, FUN = "median"), by = "sample")
      }
      if (length(est_PTRsEach) < 3) {
        next
      }
      qmax <- NULL
      rmax <- NULL
      est_PTRsIntBest <- NULL
      cormax <- 0
      for (q in 0:1) {
        for (r in (q + 1):2) {
          est_PTRsq <- est_PTRsEach[[q + 1]]
          est_PTRsr <- est_PTRsEach[[r + 1]]
          intSamples <- intersect(setdiff(Samples_filteredY1, selectSamples[[q + 1]]), setdiff(Samples_filteredY1, selectSamples[[r + 1]]))
          est_PTRsInt <- merge(est_PTRsq[est_PTRsq$sample %in% intSamples, c("sample", "est_PTR")], est_PTRsr[est_PTRsr$sample %in% intSamples, c("sample", "est_PTR")], by = "sample")
          minor_sample_q <- cor_diff(est_PTRsq)
          minor_sample_r <- cor_diff(est_PTRsr)

          corqr <- cor(est_PTRsInt[, 3], est_PTRsInt[, 2])
          if (corqr > cormax & length(minor_sample_q) == 0 & length(minor_sample_r) == 0) {
            cormax <- corqr
            qmax <- est_PTRsq
            rmax <- est_PTRsr
            est_PTRsIntBest <- est_PTRsInt
          }
        }
      }
      if (cormax > 0.98) {
        rownames(est_PTRsIntBest) <- est_PTRsIntBest$sample
        est_PTRsInt <- subset(est_PTRsIntBest, select = -c(sample))

        est_PTRsIntPCA <- prcomp(est_PTRsInt)

        qmax$Test_PTR <- (qmax$est_PTR - mean(est_PTRsInt$est_PTR.x)) / est_PTRsIntPCA$rotation[1, 1]
        rmax$Test_PTR <- (rmax$est_PTR - mean(est_PTRsInt$est_PTR.y)) / est_PTRsIntPCA$rotation[2, 1]
        rmax$Test_PTR2 <- rmax$Test_PTR * est_PTRsIntPCA$rotation[1, 1] + mean(est_PTRsInt$est_PTR.x)
        qmax$Test_PTR2 <- qmax$Test_PTR * est_PTRsIntPCA$rotation[2, 1] + mean(est_PTRsInt$est_PTR.y)

        if (test_reasonable(qmax$Test_PTR2, rmax$est_PTR) > test_reasonable(rmax$Test_PTR2, qmax$est_PTR) & test_reasonable(qmax$Test_PTR2, rmax$est_PTR) > 0.2) {
          est_PTRfinal <- df_transfer(rmax, qmax)
          break
        } else if (test_reasonable(qmax$Test_PTR2, rmax$est_PTR) < test_reasonable(rmax$Test_PTR2, qmax$est_PTR) & test_reasonable(rmax$Test_PTR2, qmax$est_PTR) > 0.2) {
          est_PTRfinal <- df_transfer(qmax, rmax)
          break
        } else {
          next
        }
      } else {
        next
      }
    }
    if (nrow(est_PTRfinal) < length(Samples_filteredY1)) {
      stop("fail to calculate consistent PTRs or combine PTRs from two subsets of samples!")
    } else {
      est_PTRs2 <- est_PTRfinal
    }
  } else {
    stop("fail to calculate consistent PTRs!")
  }
}
