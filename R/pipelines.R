#' Attempt the default iteration for contigs
#' Requires at least 20 contigs
#'
#' @param X cov3 dataframe
#' @return est_ptrs dataframe on success, a message otherwise
#'
#' @importFrom stats aggregate prcomp var
contigs_pipeline <- function(X) {
  cor_cutoff <- 0.98
  max_cor <- 0
  nrm <- floor(length(levels(X$contig)) / 5)

  for (s2 in 1:demic_env$MAX_ITER) {
    contigs <- rand_ordered_contigs(X)
    na_contig_ids <- NULL

    for (x in 1:4) {
      if (x %in% na_contig_ids) {
        return(paste("Found invalid contig (", x, ")"))
      }

      for (y in (x + 1):5) {
        if (y %in% na_contig_ids) {
          next
        }

        contig_subset_x <- contigs[(nrm * (x - 1) + 1):(nrm * x), 1]
        contig_subset_y <- contigs[(nrm * (y - 1) + 1):(nrm * y), 1]

        est_ptrs <- compare_x_y(X, contig_subset_x, contig_subset_y, na_contig_ids, cor_cutoff, max_cor)

        if (length(est_ptrs) == 0) {
          na_contig_ids <- c(na_contig_ids, y)
          next
        } else {
          return(est_ptrs)
        }
      }
    }

    if (max_cor < 0.9) {
      return(paste("Correlation too low (", max_cor, ") (min 0.90)"))
    } else if (max_cor < 0.95) {
      cor_cutoff <- 0.95
    }
  }

  return("Contigs pipeline completed without success")
}

#' Attempt alternative iteration for samples
#' Requires at least 3 samples
#'
#' @param X cov3 dataframe
#' @param max_candidate_iter max number of tries for samples pipeline iteration
#' @param contigs_pipeline_msg message from contigs_pipeline failure
#' @return est_ptrs dataframe
#'
#' @importFrom stats prcomp aggregate
samples_pipeline <- function(X, max_candidate_iter = 10, contigs_pipeline_msg = "") {
  contig <- PC1 <- NULL

  pipelineY <- iterate_pipelines(X)
  if (length(pipelineY) == 1) {
    stop(paste("contigs_pipeline: ", contigs_pipeline_msg, "\nsamples_pipeline: pipeline failed"))
  }
  Samples_filtered1 <- pipelineY[[1]]
  summeryMeanYSortFilteredSampleContig1 <- pipelineY[[2]]
  contigPCAPC1Filtered1 <- pipelineY[[3]]
  range1 <- pipelineY[[4]]
  Samples_filteredY1 <- pipelineY[[5]]

  sample_correct_y_PC1 <- merge(reshape2::dcast(subset(summeryMeanYSortFilteredSampleContig1, select = c("sample", "contig", "correctY")), contig ~ sample, value.var = "correctY"), contigPCAPC1Filtered1)

  lm_model_co <- apply(subset(sample_correct_y_PC1, select = -c(contig, PC1)), 2, lm_column, y = sample_correct_y_PC1$PC1)
  cor_model <- apply(subset(sample_correct_y_PC1, select = -c(contig, PC1)), 2, function(x) cor.test(sample_correct_y_PC1$PC1, x)$estimate)

  est_ptrs <- data.frame("est_ptr" = 2^abs(lm_model_co[1, ] * (range1[1] - range1[2])), "coefficient" = lm_model_co[1, ], "pValue" = lm_model_co[2, ], "cor" = cor_model)
  est_ptrs$sample <- rownames(est_ptrs)
  est_ptrs3 <- merge(est_ptrs, aggregate(correctY ~ sample, summeryMeanYSortFilteredSampleContig1, FUN = "median"), by = "sample")

  minor_sample3 <- cor_diff(est_ptrs3)

  if (length(minor_sample3) == 0 & max(est_ptrs3$est_ptr) >= 1.8 & max(est_ptrs3$est_ptr) / min(est_ptrs3$est_ptr) <= 5) {
    est_ptrs2 <- est_ptrs3
  } else if ((length(minor_sample3) > 0 | max(est_ptrs3$est_ptr) < 1.8 | max(est_ptrs3$est_ptr) / min(est_ptrs3$est_ptr) > 5) & length(Samples_filteredY1) >= 6) {
    est_ptrfinal <- NULL
    for (s in 1:max_candidate_iter) {
      set.seed(s)
      designateR <- sample.int(10000, size = length(Samples_filteredY1), replace = FALSE)
      SampleDesignateR <- data.frame("Sample" = Samples_filteredY1, "number" = designateR)
      SampleDesignateRSort <- SampleDesignateR[order(SampleDesignateR[, 2]), ]
      selectSamples <- NULL
      pipelineX <- NULL
      est_ptrsEach <- NULL
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

        sample_correct_y_PC1 <- merge(reshape2::dcast(subset(summeryMeanYSortFilteredSampleContigX, select = c("sample", "contig", "correctY")), contig ~ sample, value.var = "correctY"), contigPCAPC1FilteredX)

        lm_model_co <- apply(subset(sample_correct_y_PC1, select = -c(contig, PC1)), 2, lm_column, y = sample_correct_y_PC1$PC1)
        cor_model <- apply(subset(sample_correct_y_PC1, select = -c(contig, PC1)), 2, function(x) cor.test(sample_correct_y_PC1$PC1, x)$estimate)

        est_ptrs <- data.frame("est_ptr" = 2^abs(lm_model_co[1, ] * (rangeX[1] - rangeX[2])), "coefficient" = lm_model_co[1, ], "pValue" = lm_model_co[2, ], "cor" = cor_model)
        est_ptrs$sample <- rownames(est_ptrs)
        est_ptrsEach[[q + 1]] <- merge(est_ptrs, aggregate(correctY ~ sample, summeryMeanYSortFilteredSampleContig1, FUN = "median"), by = "sample")
      }
      if (length(est_ptrsEach) < 3) {
        next
      }
      qmax <- NULL
      rmax <- NULL
      est_ptrsIntBest <- NULL
      cormax <- 0
      for (q in 0:1) {
        for (r in (q + 1):2) {
          est_ptrsq <- est_ptrsEach[[q + 1]]
          est_ptrsr <- est_ptrsEach[[r + 1]]
          intSamples <- intersect(setdiff(Samples_filteredY1, selectSamples[[q + 1]]), setdiff(Samples_filteredY1, selectSamples[[r + 1]]))
          est_ptrsInt <- merge(est_ptrsq[est_ptrsq$sample %in% intSamples, c("sample", "est_ptr")], est_ptrsr[est_ptrsr$sample %in% intSamples, c("sample", "est_ptr")], by = "sample")
          minor_sample_q <- cor_diff(est_ptrsq)
          minor_sample_r <- cor_diff(est_ptrsr)

          corqr <- cor(est_ptrsInt[, 3], est_ptrsInt[, 2])
          if (corqr > cormax & length(minor_sample_q) == 0 & length(minor_sample_r) == 0) {
            cormax <- corqr
            qmax <- est_ptrsq
            rmax <- est_ptrsr
            est_ptrsIntBest <- est_ptrsInt
          }
        }
      }
      if (cormax > 0.98) {
        rownames(est_ptrsIntBest) <- est_ptrsIntBest$sample
        est_ptrsInt <- subset(est_ptrsIntBest, select = -c(sample))

        est_ptrsIntPCA <- prcomp(est_ptrsInt)

        qmax$Test_ptr <- (qmax$est_ptr - mean(est_ptrsInt$est_ptr.x)) / est_ptrsIntPCA$rotation[1, 1]
        rmax$Test_ptr <- (rmax$est_ptr - mean(est_ptrsInt$est_ptr.y)) / est_ptrsIntPCA$rotation[2, 1]
        rmax$Test_ptr2 <- rmax$Test_ptr * est_ptrsIntPCA$rotation[1, 1] + mean(est_ptrsInt$est_ptr.x)
        qmax$Test_ptr2 <- qmax$Test_ptr * est_ptrsIntPCA$rotation[2, 1] + mean(est_ptrsInt$est_ptr.y)

        if (test_reasonable(qmax$Test_ptr2, rmax$est_ptr) > test_reasonable(rmax$Test_ptr2, qmax$est_ptr) & test_reasonable(qmax$Test_ptr2, rmax$est_ptr) > 0.2) {
          est_ptrfinal <- df_transfer(rmax, qmax)
          break
        } else if (test_reasonable(qmax$Test_ptr2, rmax$est_ptr) < test_reasonable(rmax$Test_ptr2, qmax$est_ptr) & test_reasonable(rmax$Test_ptr2, qmax$est_ptr) > 0.2) {
          est_ptrfinal <- df_transfer(qmax, rmax)
          break
        } else {
          next
        }
      } else {
        next
      }
    }
    if (nrow(est_ptrfinal) < length(Samples_filteredY1)) {
      stop("fail to calculate consistent PTRs or combine PTRs from two subsets of samples!")
    } else {
      est_ptrs2 <- est_ptrfinal
    }
  } else {
    stop("fail to calculate consistent PTRs!")
  }
}
