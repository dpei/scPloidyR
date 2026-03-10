#' @title HMM Core Algorithms
#' @description Core routines for learning reference parameters and fitting the CNV HMM.
#' @name hmm_core
#' @keywords internal
NULL

#' @title Expand Reference CN to Per-Amplicon Vector
#' @description Internal helper to convert scalar/per-chromosome/per-amplicon reference_cn
#'   to a per-amplicon vector.
#' @param reference_cn Numeric value(s) specifying reference cell copy number. Can be:
#'   - Scalar: uniform CN for all amplicons (e.g., 2 or 3)
#'   - Named vector: per-chromosome CN (e.g., c("chr10" = 3, "chr21" = 3))
#'   - Unnamed vector: per-amplicon CN (length must match n_amplicons)
#' @param n_amplicons Integer, number of amplicons.
#' @param amplicon_chr_labels Character/factor vector of chromosome labels for each amplicon.
#'   Required when reference_cn is a named (per-chromosome) vector.
#' @return Numeric vector of length n_amplicons with copy number for each amplicon.
#' @keywords internal
.expand_reference_cn <- function(reference_cn, n_amplicons, amplicon_chr_labels = NULL) {
  # Case 1: Named vector (per-chromosome)
  if (!is.null(names(reference_cn))) {
    if (is.null(amplicon_chr_labels)) {
      stop("amplicon_chr_labels is required when reference_cn is a named (per-chromosome) vector")
    }
    if (length(amplicon_chr_labels) != n_amplicons) {
      stop(sprintf("amplicon_chr_labels length (%d) must match number of amplicons (%d)",
                   length(amplicon_chr_labels), n_amplicons))
    }

    if (!is.numeric(reference_cn) || any(reference_cn < 1) || any(reference_cn > 10)) {
      stop("reference_cn values must be positive numerics between 1 and 10")
    }

    chr_names <- as.character(amplicon_chr_labels)
    reference_cn_per_amplicon <- numeric(n_amplicons)

    for (i in seq_len(n_amplicons)) {
      chr <- chr_names[i]
      if (chr %in% names(reference_cn)) {
        reference_cn_per_amplicon[i] <- reference_cn[chr]
      } else {
        reference_cn_per_amplicon[i] <- 2
      }
    }

    return(reference_cn_per_amplicon)
  }

  # Case 2: Scalar (unnamed, length 1)
  if (length(reference_cn) == 1) {
    if (!is.numeric(reference_cn) || reference_cn < 1 || reference_cn > 10) {
      stop("reference_cn must be a positive numeric value between 1 and 10")
    }
    return(rep(reference_cn, n_amplicons))
  }

  # Case 3: Per-amplicon vector (unnamed, correct length)
  if (length(reference_cn) == n_amplicons) {
    if (!is.numeric(reference_cn) || any(reference_cn < 1) || any(reference_cn > 10)) {
      stop("reference_cn values must be positive numerics between 1 and 10")
    }
    return(reference_cn)
  }

  stop(sprintf(
    "reference_cn must be scalar, per-amplicon vector (length %d), or named per-chromosome vector",
    n_amplicons
  ))
}

#' @title Learn Reference Parameters from Reference Cells
#' @description Learns baseline parameters from reference cells to improve
#'   copy number estimation for test cells. Extracts stable, population-level
#'   estimates of per-amplicon diploid depth baseline, depth variance, BAF sigma,
#'   and per-amplicon heterozygosity rates.
#' @param reference_depth_matrix Numeric matrix (amplicons x reference_cells) of depth values.
#' @param reference_baf_list List of lists, where reference_baf_list[[cell_idx]][[amplicon_idx]]
#'   contains BAF values for that amplicon in that reference cell.
#' @param reference_cn Numeric value(s) specifying reference cell copy number. Can be:
#'   - Scalar: uniform CN for all amplicons (default: 2 for diploid reference)
#'   - Named vector: per-chromosome CN (e.g., c("chr10" = 3) for trisomy 10)
#'   - Unnamed vector: per-amplicon CN (length must match nrow(reference_depth_matrix))
#' @param amplicon_chr_labels Character/factor vector of chromosome labels for each amplicon.
#'   Required when reference_cn is a named (per-chromosome) vector.
#' @return A list containing:
#'   \item{diploid_depth_per_amplicon}{Numeric vector of per-amplicon diploid depth baseline.}
#'   \item{depth_variance}{Numeric, estimated depth variance on normalized scale.}
#'   \item{sigma_maf_init}{Numeric, initial standard deviation for MAF distribution.}
#'   \item{heterozygosity_rate_by_amplicon}{Numeric vector of per-amplicon heterozygosity rates.}
#'   \item{reference_cn_per_amplicon}{Numeric vector of reference CN used for each amplicon.}
#' @keywords internal
learn_reference_parameters <- function(reference_depth_matrix,
                                       reference_baf_list,
                                       reference_cn = 2,
                                       amplicon_chr_labels = NULL) {
  if (!is.matrix(reference_depth_matrix)) {
    stop("reference_depth_matrix must be a numeric matrix (amplicons x reference_cells)")
  }
  if (!is.list(reference_baf_list)) {
    stop("reference_baf_list must be a list of lists")
  }

  n_amplicons <- nrow(reference_depth_matrix)
  n_ref_cells <- ncol(reference_depth_matrix)

  if (length(reference_baf_list) != n_ref_cells) {
    stop("Length of reference_baf_list must match number of columns in reference_depth_matrix")
  }

  reference_cn_per_amplicon <- .expand_reference_cn(
    reference_cn = reference_cn,
    n_amplicons = n_amplicons,
    amplicon_chr_labels = amplicon_chr_labels
  )

  raw_median_per_amplicon <- apply(reference_depth_matrix, 1, stats::median, na.rm = TRUE)
  diploid_depth_per_amplicon <- raw_median_per_amplicon * (2 / reference_cn_per_amplicon)

  normalized_depths <- reference_depth_matrix / diploid_depth_per_amplicon
  depth_variance <- stats::var(as.vector(normalized_depths), na.rm = TRUE)

  het_counts <- numeric(n_amplicons)
  total_counts <- numeric(n_amplicons)

  for (cell_idx in seq_len(n_ref_cells)) {
    cell_baf <- reference_baf_list[[cell_idx]]
    if (!is.list(cell_baf)) {
      next
    }

    for (amp_idx in seq_len(n_amplicons)) {
      if (amp_idx > length(cell_baf)) {
        next
      }

      maf_values <- cell_baf[[amp_idx]]
      if (is.null(maf_values)) {
        next
      }

      maf_values <- pmin(maf_values, 1 - maf_values)
      maf_values <- maf_values[!is.na(maf_values)]
      if (!length(maf_values)) {
        next
      }

      het_counts[amp_idx] <- het_counts[amp_idx] + sum(maf_values > HETEROZYGOSITY_THRESHOLD)
      total_counts[amp_idx] <- total_counts[amp_idx] + length(maf_values)
    }
  }

  heterozygosity_rate_by_amplicon <- ifelse(total_counts > 0,
                                            het_counts / total_counts,
                                            0)
  heterozygosity_rate_by_amplicon <- pmax(0, pmin(1, heterozygosity_rate_by_amplicon))

  all_mafs <- unlist(reference_baf_list)
  all_mafs <- all_mafs[!is.na(all_mafs)]

  if (length(all_mafs) == 0) {
    sigma_maf_init <- 0.08
    warning("No BAF data in reference cells. Using default sigma_maf = 0.08")
  } else {
    het_mafs <- all_mafs[all_mafs > 0.3 & all_mafs < 0.7]

    if (length(het_mafs) >= 10) {
      sigma_maf_init <- stats::sd(het_mafs - 0.5)
      sigma_maf_init <- min(max(sigma_maf_init, MIN_SIGMA_MAF), MAX_SIGMA_MAF)
    } else {
      sigma_maf_init <- stats::sd(all_mafs)
      sigma_maf_init <- min(max(sigma_maf_init, MIN_SIGMA_MAF), MAX_SIGMA_MAF)
    }
  }

  if (any(is.na(diploid_depth_per_amplicon)) || any(diploid_depth_per_amplicon <= 0)) {
    stop("Invalid diploid_depth_per_amplicon computed. Check reference cell data.")
  }
  if (is.na(depth_variance) || depth_variance <= 0) {
    warning("Invalid depth_variance computed. Using default value.")
    depth_variance <- stats::var(as.vector(reference_depth_matrix), na.rm = TRUE)
  }

  return(list(
    diploid_depth_per_amplicon = diploid_depth_per_amplicon,
    depth_variance = depth_variance,
    sigma_maf_init = sigma_maf_init,
    heterozygosity_rate_by_amplicon = heterozygosity_rate_by_amplicon,
    reference_cn_per_amplicon = reference_cn_per_amplicon
  ))
}

#' @title Fit Hidden Markov Model for Copy Number Detection
#' @description Fits a Hidden Markov Model to detect copy number states from depth and BAF
#'   data using EM algorithm for parameter estimation and Viterbi decoding for state prediction.
#' @param depth Numeric vector of depth observations.
#' @param baf Numeric vector or list of B-allele frequencies (0-1).
#' @param reference_params List from learn_reference_parameters() containing:
#'   diploid_depth_per_amplicon, depth_variance, sigma_maf_init,
#'   heterozygosity_rate_by_amplicon.
#' @param states Integer vector of copy number states to model. Defaults to 0:5.
#' @param max_iter Integer, maximum number of EM iterations.
#' @param tol Numeric, convergence tolerance for EM algorithm.
#' @param eps_init Numeric, initial transition probability.
#' @return A list containing:
#'   \item{path}{Integer vector of predicted copy number states.}
#'   \item{mu_depth}{Numeric vector of estimated depth means for each state.}
#'   \item{var_depth}{Numeric, estimated depth variance.}
#'   \item{sigma_maf}{Numeric, estimated standard deviation of MAF distribution.}
#'   \item{epsilon}{Numeric, estimated transition probability.}
#'   \item{logLik}{Numeric, final log-likelihood.}
#'   \item{iter}{Integer, number of EM iterations performed.}
#' @keywords internal
fit_hmm <- function(depth,
                    baf,
                    reference_params,
                    states   = 0:5,
                    max_iter = HMM_DEFAULTS$MAX_ITER,
                    tol      = HMM_DEFAULTS$TOLERANCE,
                    eps_init = HMM_DEFAULTS$EPS_INIT) {

  if (missing(reference_params) || !is.list(reference_params)) {
    stop("reference_params is required and must be a list from learn_reference_parameters()")
  }
  validate_required_fields(
    reference_params,
    c("diploid_depth_per_amplicon", "depth_variance", "sigma_maf_init",
      "heterozygosity_rate_by_amplicon"),
    "reference_params"
  )

  y_raw <- depth
  baf_list <- if (is.list(baf)) baf else as.list(baf)
  n  <- length(y_raw)
  K  <- length(states)
  maf_list <- lapply(baf_list, function(v) {
    if (is.null(v)) {
      return(numeric(0))
    }
    v <- as.numeric(v)
    v <- v[!is.na(v)]
    if (!length(v)) {
      return(numeric(0))
    }
    v <- pmin(v, 1 - v)
    v <- v[!is.na(v)]
    v
  })

  het_rate_by_amplicon <- reference_params$heterozygosity_rate_by_amplicon
  if (length(het_rate_by_amplicon) != n) {
    stop(sprintf(
      "heterozygosity_rate_by_amplicon must have length %d (got %d)",
      n, length(het_rate_by_amplicon)
    ))
  }
  het_rate_by_amplicon <- pmax(0, pmin(1, het_rate_by_amplicon))

  mean_het_rate <- mean(het_rate_by_amplicon, na.rm = TRUE)
  low_het_strength <- 0
  if (is.finite(mean_het_rate) && is.finite(LOW_HET_RATE_THRESHOLD) && LOW_HET_RATE_THRESHOLD > 0) {
    low_het_strength <- (LOW_HET_RATE_THRESHOLD - mean_het_rate) / LOW_HET_RATE_THRESHOLD
    low_het_strength <- max(0, min(1, low_het_strength))
  }

  mu_prior_weight_base <- HMM_DEFAULTS$MU_PRIOR_WEIGHT_LOW_HET
  if (is.null(mu_prior_weight_base) || !is.finite(mu_prior_weight_base) || mu_prior_weight_base <= 0) {
    mu_prior_weight_base <- 10
  }

  diploid_depth_per_amplicon <- reference_params$diploid_depth_per_amplicon
  diploid_depth_per_amplicon <- pmax(diploid_depth_per_amplicon, 1e-3)
  y <- y_raw / diploid_depth_per_amplicon

  expected_maf_matrix <- matrix(NA_real_, nrow = K, ncol = n)
  informative_states_idx <- which(states >= 1)

  genotype_info <- vector("list", K)
  genotype_priors_by_amplicon <- vector("list", K)
  for (k in seq_len(K)) {
    s <- states[k]
    if (s >= 1) {
      genotype_names <- names(get_genotype_prior(s, 0.5))
      expected_mafs <- sapply(genotype_names, get_expected_maf)
      expected_mafs_safe <- expected_mafs

      genotype_info[[k]] <- list(
        names = genotype_names,
        expected_mafs = expected_mafs_safe
      )

      genotype_priors_by_amplicon[[k]] <- vector("list", n)
      for (i in seq_len(n)) {
        priors <- get_genotype_prior(s, het_rate_by_amplicon[i])
        priors <- priors[genotype_names]
        priors[is.na(priors)] <- 0
        if (sum(priors) > 0) {
          priors <- priors / sum(priors)
        }
        genotype_priors_by_amplicon[[k]][[i]] <- priors
        expected_maf_matrix[k, i] <- sum(priors * expected_mafs_safe)
      }
    }
  }

  mu_prior <- states / 2
  mu_prior[states == 0] <- 1e-6
  mu <- mu_prior

  var_y     <- reference_params$depth_variance
  sigma_maf <- reference_params$sigma_maf_init
  eps       <- eps_init

  initial_probs <- rep(1 / K, K)
  log_initial_probs <- log(initial_probs)

  logA <- matrix(log(eps / (K - 1)), K, K)
  diag(logA) <- log(1 - eps)

  log_em <- matrix(0, K, n)

  for (it in seq_len(max_iter)) {

    ## E-step
    for (k in seq_len(K)) {
      s <- states[k]

      for (i in seq_len(n)) {
        log_depth <- log_norm(y[i], mu[k], var_y)

        if (s >= 1 && length(maf_list[[i]]) > 0 && !is.null(genotype_info[[k]])) {
          info <- genotype_info[[k]]
          priors <- genotype_priors_by_amplicon[[k]][[i]]
          n_genotypes <- length(info$names)

          if (is.null(priors) || length(priors) != n_genotypes) {
            priors <- rep(0, n_genotypes)
          }

          log_baf_per_genotype <- rep(-Inf, n_genotypes)

          for (g_idx in seq_len(n_genotypes)) {
            prior <- priors[g_idx]
            if (is.na(prior) || prior <= 0) {
              next
            }

            expected_maf_safe <- info$expected_mafs[g_idx]
            log_baf_likelihood <- sum(log_truncnorm(maf_list[[i]], expected_maf_safe, sigma_maf))
            log_baf_per_genotype[g_idx] <- log(prior) + log_baf_likelihood
          }

          log_baf <- logsumexp(log_baf_per_genotype)
          log_em[k, i] <- log_depth + log_baf
        } else {
          log_em[k, i] <- log_depth
        }
      }
    }

    log_em[is.na(log_em)] <- -1e6

    # Forward algorithm
    logalpha <- matrix(NA_real_, K, n)
    logalpha[, 1] <- log_initial_probs + log_em[, 1]
    for (t in 2:n) {
      for (k in 1:K) {
        logalpha[k, t] <- log_em[k, t] + logsumexp(logalpha[, t - 1] + logA[, k])
      }
    }

    LL <- logsumexp(logalpha[, n])

    # Backward algorithm
    logbeta <- matrix(0, K, n)
    for (t in (n - 1):1) {
      for (k in 1:K) {
        logbeta[k, t] <- logsumexp(logA[k, ] + log_em[, t + 1] + logbeta[, t + 1])
      }
    }

    loggamma <- logalpha + logbeta
    loggamma <- sweep(loggamma, 2, apply(loggamma, 2, logsumexp), FUN = "-")
    gamma <- exp(loggamma)

    n_transitions <- n - 1
    if (n_transitions > 0) {
      logxi <- array(NA_real_, dim = c(K, K, n_transitions))
      for (t in 1:n_transitions) {
        term1 <- matrix(logalpha[, t], K, K, byrow = FALSE)
        term2 <- logA
        term3 <- matrix(log_em[, t + 1] + logbeta[, t + 1], K, K, byrow = TRUE)
        logxi[, , t] <- term1 + term2 + term3
        logxi[, , t] <- logxi[, , t] - logsumexp(logxi[, , t])
      }
      xi <- exp(logxi)
    } else {
      xi <- array(0, dim = c(K, K, 0))
    }

    ## M-step
    mu_num <- rowSums(gamma * y[col(gamma)], na.rm = TRUE)
    mu_den <- rowSums(gamma, na.rm = TRUE)
    mu_new <- mu_num / mu_den

    if (low_het_strength > 0) {
      mu_prior_weight <- mu_prior_weight_base * low_het_strength
      mu_new <- (mu_num + mu_prior_weight * mu_prior) / (mu_den + mu_prior_weight)
    }

    var_new <- sum(gamma * (y[col(gamma)] - mu_new[row(gamma)])^2, na.rm = TRUE) / sum(!is.na(y))

    if (length(informative_states_idx)) {
      variance_num <- 0
      variance_denom <- 0

      for (ii in seq_len(n)) {
        Li <- length(maf_list[[ii]])
        if (Li == 0) next

        for (k_idx in seq_along(informative_states_idx)) {
          k <- informative_states_idx[k_idx]
          expected_maf <- expected_maf_matrix[k, ii]
          if (is.na(expected_maf)) {
            next
          }

          sq_resid <- sum((maf_list[[ii]] - expected_maf)^2)
          variance_num <- variance_num + gamma[k, ii] * sq_resid
          variance_denom <- variance_denom + gamma[k, ii] * Li
        }
      }

      sigma_maf_new <- sqrt(variance_num / variance_denom)
      sigma_maf_new <- max(sigma_maf_new, 0.01, na.rm = TRUE)
      if (is.nan(sigma_maf_new)) sigma_maf_new <- sigma_maf
    }

    if (n_transitions > 0) {
      eps_new <- sum(apply(xi, 3, function(slice) sum(slice) - sum(diag(slice)))) / n_transitions
    } else {
      eps_new <- eps
    }
    eps_new <- min(max(eps_new, 1e-6), 0.5)

    if (it > 1 && max(abs(mu_new - mu), abs(var_new - var_y),
            abs(sigma_maf_new - sigma_maf), abs(eps_new - eps), na.rm = TRUE) < tol) {
      break
    }

    mu        <- mu_new
    var_y     <- var_new
    sigma_maf <- sigma_maf_new
    eps       <- eps_new

    logA <- matrix(log(eps / (K - 1)), K, K)
    diag(logA) <- log(1 - eps)
  }

  # Viterbi decoding
  delta <- matrix(NA_real_, K, n)
  psi   <- matrix(NA_integer_, K, n)
  delta[, 1] <- log_initial_probs + log_em[, 1]
  for (t in 2:n) {
    for (k in 1:K) {
      tmp <- delta[, t - 1] + logA[, k]
      psi[k, t] <- which.max(tmp)
      delta[k, t] <- log_em[k, t] + max(tmp)
    }
  }
  path <- integer(n)
  path[n] <- which.max(delta[, n])
  for (t in (n - 1):1) {
    path[t] <- psi[path[t + 1], t + 1]
  }

  list(
    path        = states[path],
    mu_depth    = mu,
    var_depth   = var_y,
    sigma_maf   = sigma_maf,
    epsilon     = eps,
    logLik      = LL,
    iter        = it
  )
}
