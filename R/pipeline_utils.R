#' Verify input data meets pipeline criteria
#'
#' @param X cov3 matrix as a dataframe
#' @return TRUE if the input data meets the criteria, FALSE otherwise
verify_input <- function(X) {
  if (!length(X)) {
    stop("Input data is empty")
  }
  if (length(unique(X$contig)) < demic_env$MIN_CONTIGS) {
    stop(paste("Not enough contigs", length(unique(X$contig)), "<", demic_env$MIN_CONTIGS))
  }
  if (length(unique(X$sample)) < demic_env$MIN_SAMPLES) {
    stop(paste("Not enough samples", length(unique(X$sample)), "<", demic_env$MIN_SAMPLES))
  }
  if (is.null(levels(X$contig))) {
    stop("levels(X$contig) is empty, did you remember to read in the data with strings as factors? e.g. X <- read.csv('/path/to/file.cov3', stringsAsFactors=TRUE)")
  }

  return(TRUE)
}

#' Combine output from the contigs pipeline and samples pipeline
#'
#' @param contigs output from the contigs pipeline
#' @param samples output from the samples pipeline
#' @return a dataframe with the combined estimated PTRs
combine_ests <- function(contigs, samples) {
  est_ptrs <- contigs
  # est_ptrs <- list(contigs=contigs, samples=samples)

  est_ptrs
}

#' Gives a randomly ordered set of contigs
#'
#' @param X input data frame that includes contigs
#' @return contigs randomly ordered
rand_ordered_contigs <- function(X) {
  rands <- sample.int(10000, size = length(levels(X$contig)), replace = FALSE)
  contigs <- data.frame("Contig" = levels(X$contig), "number" = rands)
  contigs <- contigs[order(contigs[, 2]), ]

  contigs
}

#' Get PTR estimates for output of the core pipeline on a subset of data
#'
#' @param p is the pipeline named list
#' @return a dataframe
#' \itemize{
#'   \item sample: sample
#'   \item est_ptr: PTR estimate
#'   \item coefficient: coefficient of linear regression
#'   \item pValue: p-value of linear regression
#'   \item cor: correlation coefficient
#'   \item correctY: corrected coverage
#' }
est_ptrs_subset <- function(p) {
  PC1 <- contig <- NULL

  sample_correct_y_PC1 <- merge(reshape2::dcast(subset(p$correct_ys, select = c("sample", "contig", "correctY")), contig ~ sample, value.var = "correctY"), p$pc1)

  lm_model_co <- apply(subset(sample_correct_y_PC1, select = -c(contig, PC1)), 2, lm_column, y = sample_correct_y_PC1$PC1)
  cor_model <- apply(subset(sample_correct_y_PC1, select = -c(contig, PC1)), 2, function(x) cor.test(sample_correct_y_PC1$PC1, x)$estimate)

  est_ptrs <- data.frame("est_ptr" = 2^abs(lm_model_co[1, ] * (p$pc1_range[1] - p$pc1_range[2])), "coefficient" = lm_model_co[1, ], "pValue" = lm_model_co[2, ], "cor" = cor_model)
  est_ptrs$sample <- rownames(est_ptrs)

  merge(est_ptrs, aggregate(correctY ~ sample, p$correct_ys, FUN = "median"), by = "sample")
}

#' Run mixed linear model with random effect using lme4
#'
#' @param X input data frame
#' @return a dataframe
#'
#' @importFrom lme4 lmer
#' @importFrom stats coef
lme4_model <- function(X) {
  lmeModel <- lmer(log_cov ~ GC_content + (1 | sample:contig), data = X, REML = FALSE)

  lmeModelCoef <- coef(lmeModel)$`sample:contig`
  lmeModelCoef$s_c <- rownames(lmeModelCoef)

  lmeModelCoef
}

#' Compares contig subset x against contig subset y
#'
#' @param X input data frame
#' @param contig_subset_x the first set of contigs
#' @param contig_subset_y the second set of contigs
#' @param na_contig_ids the set of bad contigs discovered so far
#' @param cor_cutoff the correlation cutoff
#' @param max_cor the max correlation
#' @return a dataframe
#' \itemize{
#'   \item sample: sample
#'   \item est_ptr: PTR estimate
#'   \item coefficient: coefficient of linear regression
#'   \item pValue: p-value of linear regression
#'   \item cor: correlation coefficient
#'   \item correctY: corrected coverage
#' }
compare_x_y <- function(X, contig_subset_x, contig_subset_y, na_contig_ids, cor_cutoff, max_cor) {
  pipeline_x <- iterate_pipelines(X[!X$contig %in% contig_subset_x, ])
  if (length(pipeline_x) == 1) {
    return(paste("Pipeline failed (", pipeline_x, ")"))
  }

  pipeline_y <- iterate_pipelines(X[!X$contig %in% contig_subset_y, ])
  if (length(pipeline_y) == 1) {
    return()
  }

  if (length(pipeline_x$pc1$contig) - length(intersect(pipeline_x$pc1$contig, pipeline_y$pc1$contig)) < 3 || length(pipeline_y$pc1$contig) - length(intersect(pipeline_x$pc1$contig, pipeline_y$pc1$contig)) < 3) {
    # Not enough unique contigs in each set
    return()
  }

  est_ptrs_x <- est_ptrs_subset(pipeline_x)
  est_ptrs_y <- est_ptrs_subset(pipeline_y)

  minor_sample1 <- cor_diff(est_ptrs_x)
  minor_sample2 <- cor_diff(est_ptrs_y)
  if ((length(minor_sample1) > 0 & length(minor_sample2) > 0) | (max(est_ptrs_x$est_ptr) < 1.8 & max(est_ptrs_y$est_ptr) < 1.8) | (max(est_ptrs_x$est_ptr) / min(est_ptrs_x$est_ptr) > 5 & max(est_ptrs_y$est_ptr) / min(est_ptrs_y$est_ptr) > 5)) {
    return()
  }

  est_ptrs_x_y <- merge(est_ptrs_x, est_ptrs_y, by = "sample")

  if (nrow(est_ptrs_x_y) > 0.9 * max(c(nrow(est_ptrs_x), nrow(est_ptrs_y)))) {
    if (var(est_ptrs_x_y$est_ptr.x) == 0 || var(est_ptrs_x_y$est_ptr.y) == 0) {
      est_ptrs_x_y$est_ptr.x <- jitter(est_ptrs_x_y$est_ptr.x)
      est_ptrs_x_y$est_ptr.y <- jitter(est_ptrs_x_y$est_ptr.y)
    }
    cor_current <- cor(est_ptrs_x_y$est_ptr.x, est_ptrs_x_y$est_ptr.y)
    if (is.na(cor_current)) {
      return("NA found in correlation calculation")
    }
    if (cor_current > max_cor) {
      max_cor <- cor_current
    }

    if (cor(est_ptrs_x_y$est_ptr.x, est_ptrs_x_y$est_ptr.y) > cor_cutoff) {
      est_ptrs_x_y$est_ptr <- apply(subset(est_ptrs_x_y, select = c("est_ptr.x", "est_ptr.y")), 1, mean)
      est_ptrs_x_y$coefficient <- apply(subset(est_ptrs_x_y, select = c("coefficient.x", "coefficient.y")), 1, function(x) mean(abs(x)))
      est_ptrs_x_y$pValue <- apply(subset(est_ptrs_x_y, select = c("pValue.x", "pValue.y")), 1, max)
      est_ptrs_x_y$cor <- apply(subset(est_ptrs_x_y, select = c("cor.x", "cor.y")), 1, function(x) mean(abs(x)))
      est_ptrs_x_y$correctY <- apply(subset(est_ptrs_x_y, select = c("correctY.x", "correctY.y")), 1, mean)

      subset(est_ptrs_x_y, select = c("sample", "est_ptr", "coefficient", "pValue", "cor", "correctY"))
    }
  }
}
