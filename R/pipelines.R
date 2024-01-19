#' Tries to estimate based on the whole input dataset
#'
#' @param X cov3 dataframe
#' @return est_ptrs dataframe on success, null otherwise
#' \itemize{
#'  \item est_ptr: estimated PTR values
#'  \item coefficient: coefficient of linear regression
#'  \item pValue: p-value of linear regression
#'  \item cor: correlation coefficient
#'  \item correctY: corrected coverage
#' }
#'
#' @examples
#' est_ptrs_001 <- est_ptr_from_all(ContigCluster1)
#' est_ptrs_001
#'
#' @export
est_ptr_from_all <- function(X) {
  p <- iterate_pipelines(X)
  est_ptrs <- est_ptrs_subset(p)

  if (length(cor_diff(est_ptrs)) > 0) {
    warning("Not all estimates are correlated in the same direction")
    return(NULL)
  }

  rownames(est_ptrs) <- est_ptrs$sample
  est_ptrs <- est_ptrs[, c("est_ptr", "coefficient", "pValue", "cor", "correctY")]
  est_ptrs
}

#' Tries up to max_attempts times to compare each permutation of removing random
#' subsets of contigs/samples from X, and returns the PTR estimate if a valid
#' one comes back from the comparisons
#'
#' Requires a minimum of 2 * num_subsets contigs/samples
#'
#' @param X cov3 dataframe
#' @param max_attempts max number of attempts to find a valid ptr estimate
#' @param num_subsets number of subsets to split contigs/samples into
#' @return est_ptrs dataframe on success, null otherwise
#' \itemize{
#'   \item est_ptr: estimated PTR values
#'   \item coefficient: coefficient of linear regression
#'   \item pValue: p-value of linear regression
#'   \item cor: correlation coefficient
#'   \item correctY: corrected coverage
#' }
#'
#' @examples
#' est_ptrs_001_on_contigs <- est_ptr_on(ContigCluster1, "contig", num_subsets = 5)
#' est_ptrs_001_on_contigs
#'
#' est_ptrs_001_on_samples <- est_ptr_on(ContigCluster1, "sample")
#' est_ptrs_001_on_samples
#'
#' @export
est_ptr_on <- function(X, subset_on, max_attempts = 10, num_subsets = 3, cor_cutoff = 0.98) {
  if (subset_on != "contig" && subset_on != "sample") {
    stop("subset_on must be either 'contig' or 'sample'")
  }
  if (length(unique(X[[subset_on]])) < 2 * num_subsets) {
    warning("Not enough contigs/samples to compare subsampled datasets")
    return(NULL)
  }

  max_cor <- 0

  for (i in 1:max_attempts) {
    unique_values <- sample(unique(X[[subset_on]]))
    subsets <- split(unique_values, cut(seq_along(unique_values), num_subsets, labels = FALSE))

    subset_pipelines <- lapply(subsets, function(z) iterate_pipelines(X[!X[[subset_on]] %in% z, ]))
    subset_est_ptrs <- lapply(subset_pipelines, function(z) est_ptrs_subset(z))

    for (x in 1:(num_subsets-1)) {
      for (y in (x + 1):num_subsets) {
        if (is.null(subset_est_ptrs[[x]]) || is.null(subset_est_ptrs[[y]])) {
          next
        }

        if (subset_on == "contig") {
          # WHY ARE WE KEEPING TRACK OF MAX_COR???
          # Should use it to compare all comparisons in the end instead of just taking the first result over a threshold
          comparison <- compare_contig_subsets(subset_est_ptrs[[x]], subset_est_ptrs[[y]], subset_pipelines[[x]], subset_pipelines[[y]], cor_cutoff, max_cor)
        } else if (subset_on == "sample") {
          comparison <- compare_sample_subsets(subset_est_ptrs[[x]], subset_est_ptrs[[y]], subset_pipelines[[x]], subset_pipelines[[y]], cor_cutoff, max_cor)
        }

        est_ptrs <- comparison$est_ptr
        max_cor <- comparison$max_cor

        if (!is.null(est_ptrs)) {
          rownames(est_ptrs) <- est_ptrs$sample
          est_ptrs <- est_ptrs[, c("est_ptr", "coefficient", "pValue", "cor", "correctY")]
          return(est_ptrs)
        }
      }
    }
  }

  NULL
}




samples_pipeline <- function(X, max_candidate_iter = 10) {
  contig <- PC1 <- NULL
  est_ptrfinal <- NULL

  p <- iterate_pipelines(X)
  Samples_filtered1 <- p[[1]]
  summeryMeanYSortFilteredSampleContig1 <- p[[2]]
  contigPCAPC1Filtered1 <- p[[3]]
  range1 <- p[[4]]
  Samples_filteredY1 <- p[[5]]

  for (s in 1:max_candidate_iter) {
    designateR <- sample.int(10000, size = length(Samples_filtered1), replace = FALSE)
    SampleDesignateR <- data.frame("Sample" = Samples_filtered1, "number" = designateR)
    SampleDesignateRSort <- SampleDesignateR[order(SampleDesignateR[, 2]), ]
    selectSamples <- NULL
    pipelineX <- NULL
    est_ptrsEach <- NULL
    Samples_filteredX <- NULL
    for (q in 0:2) {
      selectSamples[[q + 1]] <- SampleDesignateRSort[(1:length(Samples_filtered1)) %% 3 == q, ]$Sample

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
      browser()

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
}
