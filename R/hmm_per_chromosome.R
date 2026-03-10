#' @title Per-Chromosome HMM Implementation
#' @description Functions for fitting HMM independently to each chromosome and
#'   parallel multi-cell processing.
#' @name hmm_per_chromosome
NULL

#' @title Fit HMM Independently Per Chromosome
#' @description Fits separate HMM chains to each chromosome, eliminating cross-chromosome
#'   transitions. Each chromosome gets its own parameter estimates (fully independent).
#' @param depth Numeric vector of depth observations (length = n_amplicons).
#' @param baf List of B-allele frequency vectors (length = n_amplicons).
#' @param amplicon_chr_labels Integer or factor vector indicating chromosome for each amplicon.
#' @param reference_params List from learn_reference_parameters().
#' @param states Integer vector of copy number states to model. Defaults to 0:5.
#' @param max_iter Integer, maximum number of EM iterations.
#' @param tol Numeric, convergence tolerance.
#' @param eps_init Numeric, initial transition probability.
#' @return A list containing path, per-chromosome parameters, and convergence info.
#' @keywords internal
fit_hmm_per_chromosome <- function(depth,
                                   baf,
                                   amplicon_chr_labels,
                                   reference_params,
                                   states   = 0:5,
                                   max_iter = HMM_DEFAULTS$MAX_ITER,
                                   tol      = HMM_DEFAULTS$TOLERANCE,
                                   eps_init = HMM_DEFAULTS$EPS_INIT) {

  n_amplicons <- length(depth)
  if (length(baf) != n_amplicons) {
    stop("Length of baf must match length of depth")
  }
  if (length(amplicon_chr_labels) != n_amplicons) {
    stop("Length of amplicon_chr_labels must match length of depth")
  }

  unique_chrs <- unique(amplicon_chr_labels)
  n_chrs <- length(unique_chrs)

  path <- rep(NA_integer_, n_amplicons)
  mu_depth_by_chr <- vector("list", n_chrs)
  var_depth_by_chr <- numeric(n_chrs)
  sigma_maf_by_chr <- numeric(n_chrs)
  epsilon_by_chr <- numeric(n_chrs)
  logLik_by_chr <- numeric(n_chrs)
  iter_by_chr <- integer(n_chrs)
  converged_by_chr <- logical(n_chrs)

  names(mu_depth_by_chr) <- as.character(unique_chrs)
  names(var_depth_by_chr) <- as.character(unique_chrs)
  names(sigma_maf_by_chr) <- as.character(unique_chrs)
  names(epsilon_by_chr) <- as.character(unique_chrs)
  names(logLik_by_chr) <- as.character(unique_chrs)
  names(iter_by_chr) <- as.character(unique_chrs)
  names(converged_by_chr) <- as.character(unique_chrs)

  for (chr_idx in seq_along(unique_chrs)) {
    chr <- unique_chrs[chr_idx]
    amplicon_indices <- which(amplicon_chr_labels == chr)

    depth_chr <- depth[amplicon_indices]
    baf_chr <- baf[amplicon_indices]

    reference_params_chr <- list(
      diploid_depth_per_amplicon = reference_params$diploid_depth_per_amplicon[amplicon_indices],
      depth_variance = reference_params$depth_variance,
      sigma_maf_init = reference_params$sigma_maf_init,
      heterozygosity_rate_by_amplicon = reference_params$heterozygosity_rate_by_amplicon[amplicon_indices]
    )

    tryCatch({
      hmm_result <- fit_hmm(
        depth = depth_chr,
        baf = baf_chr,
        reference_params = reference_params_chr,
        states = states,
        max_iter = max_iter,
        tol = tol,
        eps_init = eps_init
      )

      path[amplicon_indices] <- hmm_result$path
      mu_depth_by_chr[[chr_idx]] <- hmm_result$mu_depth
      var_depth_by_chr[chr_idx] <- hmm_result$var_depth
      sigma_maf_by_chr[chr_idx] <- hmm_result$sigma_maf
      epsilon_by_chr[chr_idx] <- hmm_result$epsilon
      logLik_by_chr[chr_idx] <- hmm_result$logLik
      iter_by_chr[chr_idx] <- hmm_result$iter
      converged_by_chr[chr_idx] <- TRUE

    }, error = function(e) {
      path[amplicon_indices] <<- NA_integer_
      mu_depth_by_chr[[chr_idx]] <<- rep(NA_real_, length(states))
      var_depth_by_chr[chr_idx] <<- NA_real_
      sigma_maf_by_chr[chr_idx] <<- NA_real_
      epsilon_by_chr[chr_idx] <<- NA_real_
      logLik_by_chr[chr_idx] <<- NA_real_
      iter_by_chr[chr_idx] <<- NA_integer_
      converged_by_chr[chr_idx] <<- FALSE

      warning(sprintf("HMM failed for chromosome %s: %s", chr, e$message))
    })
  }

  list(
    path = path,
    mu_depth_by_chr = mu_depth_by_chr,
    var_depth_by_chr = var_depth_by_chr,
    sigma_maf_by_chr = sigma_maf_by_chr,
    epsilon_by_chr = epsilon_by_chr,
    logLik_by_chr = logLik_by_chr,
    iter_by_chr = iter_by_chr,
    converged_by_chr = converged_by_chr
  )
}


#' @title Call Copy Number via Per-Chromosome HMM
#' @description Main entry point for parallel multi-cell per-chromosome HMM copy number
#'   calling. Learns reference parameters from diploid cells, then fits independent
#'   HMM chains per chromosome for each cell in parallel.
#'
#' @param depth_matrix Numeric matrix of depth observations (amplicons x cells).
#' @param baf_data Either a matrix of B-allele frequencies (amplicons x cells), or
#'   a list of lists where baf_data[[cell_idx]] is the BAF list for that cell.
#' @param amplicon_chr_labels Integer or factor vector indicating chromosome for each amplicon.
#' @param reference_cell_indices Integer vector of indices for reference cells.
#' @param reference_cn Numeric value(s) specifying reference cell copy number. Can be:
#'   - Scalar: uniform CN for all amplicons (default: 2 for diploid reference)
#'   - Named vector: per-chromosome CN (e.g., c("chr10" = 3))
#'   - Unnamed vector: per-amplicon CN (length must match nrow(depth_matrix))
#' @param states Integer vector of copy number states to model. Defaults to 0:5.
#' @param max_iter Integer, maximum number of EM iterations.
#' @param tol Numeric, convergence tolerance.
#' @param eps_init Numeric, initial transition probability.
#' @param ncores Integer, number of CPU cores to use. If NULL, uses detectCores() - 1.
#' @param progress Logical, whether to print progress messages. Defaults to TRUE.
#' @param verbose Logical, whether to print detailed per-cell information. Defaults to FALSE.
#'
#' @return A list containing:
#'   \item{cn_calls_amplicon}{Matrix of amplicon-level copy number calls (amplicons x cells).}
#'   \item{mu_depth_by_chr}{List of lists: per-cell, per-chromosome mu_depth vectors.}
#'   \item{var_depth_by_chr}{Matrix (chromosomes x cells) of var_depth values.}
#'   \item{sigma_maf_by_chr}{Matrix (chromosomes x cells) of sigma_maf values.}
#'   \item{epsilon_by_chr}{Matrix (chromosomes x cells) of epsilon values.}
#'   \item{logLik_by_chr}{Matrix (chromosomes x cells) of log-likelihoods.}
#'   \item{iter_by_chr}{Matrix (chromosomes x cells) of iteration counts.}
#'   \item{converged_by_chr}{Matrix (chromosomes x cells) of convergence flags.}
#'   \item{convergence_summary}{List with convergence statistics.}
#'   \item{timing}{List with total_time, mean_time_per_cell, and ncores_used.}
#'
#' @examples
#' \dontrun{
#' # Create test data
#' depth_matrix <- matrix(rnorm(60, 100, 10), nrow = 6, ncol = 10)
#' baf_list <- lapply(1:10, function(i) lapply(1:6, function(j) rnorm(2, 0.5, 0.08)))
#' chr_labels <- rep(1:2, each = 3)
#'
#' result <- call_copy_number(
#'   depth_matrix = depth_matrix,
#'   baf_data = baf_list,
#'   amplicon_chr_labels = chr_labels,
#'   reference_cell_indices = 1:5,
#'   states = 1:5,
#'   ncores = 1
#' )
#' }
#' @export
call_copy_number <- function(depth_matrix,
                            baf_data,
                            amplicon_chr_labels,
                            reference_cell_indices,
                            reference_cn = 2,
                            states   = 0:5,
                            max_iter = HMM_DEFAULTS$MAX_ITER,
                            tol      = HMM_DEFAULTS$TOLERANCE,
                            eps_init = HMM_DEFAULTS$EPS_INIT,
                            ncores   = NULL,
                            progress = TRUE,
                            verbose  = FALSE) {

  if (!is.matrix(depth_matrix)) {
    stop("depth_matrix must be a numeric matrix (amplicons x cells)")
  }
  if (missing(reference_cell_indices) || length(reference_cell_indices) == 0) {
    stop("reference_cell_indices is required and must be a non-empty integer vector")
  }
  if (missing(amplicon_chr_labels)) {
    stop("amplicon_chr_labels is required")
  }

  n_amplicons <- nrow(depth_matrix)
  n_cells <- ncol(depth_matrix)

  if (length(amplicon_chr_labels) != n_amplicons) {
    stop("Length of amplicon_chr_labels must match number of rows in depth_matrix")
  }

  if (progress) {
    cat(sprintf("Starting parallel per-chromosome HMM analysis for %d cells across %d amplicons\n",
                n_cells, n_amplicons))
  }

  if (is.null(ncores)) {
    ncores <- parallel::detectCores() - 1
    ncores <- max(1, ncores)
  }
  ncores <- min(ncores, n_cells)

  if (progress) {
    cat(sprintf("Using %d CPU cores for parallel processing\n", ncores))
  }

  # Convert BAF data to list format if needed
  if (is.matrix(baf_data)) {
    if (progress) cat("Converting BAF matrix to list format...\n")
    baf_list <- lapply(1:n_cells, function(i) {
      as.list(baf_data[, i])
    })
  } else if (is.list(baf_data)) {
    if (length(baf_data) != n_cells) {
      stop("Length of baf_data list must match number of columns in depth_matrix")
    }
    baf_list <- baf_data
  } else {
    stop("baf_data must be either a matrix or a list of lists")
  }

  if (progress) {
    cat(sprintf("Learning reference parameters from %d reference cells...\n",
                length(reference_cell_indices)))
  }

  reference_depth_matrix <- depth_matrix[, reference_cell_indices, drop = FALSE]
  reference_baf_list <- baf_list[reference_cell_indices]

  reference_params <- learn_reference_parameters(
    reference_depth_matrix = reference_depth_matrix,
    reference_baf_list = reference_baf_list,
    reference_cn = reference_cn,
    amplicon_chr_labels = amplicon_chr_labels
  )

  if (progress && verbose) {
    cat("  Reference parameters learned:\n")
    cat(sprintf("    Median diploid depth: %.2f\n",
                stats::median(reference_params$diploid_depth_per_amplicon)))
    cat(sprintf("    Depth variance: %.2f\n", reference_params$depth_variance))
    cat(sprintf("    MAF sigma: %.3f\n", reference_params$sigma_maf_init))
  }

  unique_chrs <- unique(amplicon_chr_labels)
  n_chrs <- length(unique_chrs)

  start_time <- Sys.time()

  # Define worker function
  process_cell <- function(cell_idx) {
    tryCatch({
      depth_vec <- depth_matrix[, cell_idx]
      baf_vec <- baf_list[[cell_idx]]

      result <- fit_hmm_per_chromosome(
        depth = depth_vec,
        baf = baf_vec,
        amplicon_chr_labels = amplicon_chr_labels,
        reference_params = reference_params,
        states = states,
        max_iter = max_iter,
        tol = tol,
        eps_init = eps_init
      )

      return(list(
        success = TRUE,
        cell_idx = cell_idx,
        result = result
      ))
    }, error = function(e) {
      return(list(
        success = FALSE,
        cell_idx = cell_idx,
        error = as.character(e)
      ))
    })
  }

  # Platform-specific parallel processing
  os_type <- .Platform$OS.type

  if (os_type == "unix") {
    if (progress) cat("Using mclapply (fork-based parallelism)\n")

    if (progress && !verbose) {
      results <- vector("list", n_cells)
      chunk_size <- ceiling(n_cells / 10)

      for (chunk_start in seq(1, n_cells, by = chunk_size)) {
        chunk_end <- min(chunk_start + chunk_size - 1, n_cells)
        chunk_indices <- chunk_start:chunk_end

        chunk_results <- parallel::mclapply(
          chunk_indices,
          process_cell,
          mc.cores = ncores
        )

        results[chunk_indices] <- chunk_results

        pct_complete <- round(100 * chunk_end / n_cells)
        cat(sprintf("  Progress: %d%% (%d / %d cells)\n",
                    pct_complete, chunk_end, n_cells))
      }
    } else {
      results <- parallel::mclapply(
        1:n_cells,
        process_cell,
        mc.cores = ncores
      )
    }

  } else {
    if (progress) cat("Using parLapply (PSOCK cluster for Windows)\n")

    cl <- parallel::makeCluster(ncores)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    parallel::clusterExport(cl, c("fit_hmm_per_chromosome", "fit_hmm",
                                   "log_norm", "log_truncnorm",
                                   "logsumexp", "get_genotype_prior",
                                   "get_expected_maf", "learn_reference_parameters"),
                            envir = environment())

    if (progress && !verbose) {
      results <- vector("list", n_cells)
      chunk_size <- ceiling(n_cells / 10)

      for (chunk_start in seq(1, n_cells, by = chunk_size)) {
        chunk_end <- min(chunk_start + chunk_size - 1, n_cells)
        chunk_indices <- chunk_start:chunk_end

        chunk_results <- parallel::parLapply(
          cl,
          chunk_indices,
          process_cell
        )

        results[chunk_indices] <- chunk_results

        pct_complete <- round(100 * chunk_end / n_cells)
        cat(sprintf("  Progress: %d%% (%d / %d cells)\n",
                    pct_complete, chunk_end, n_cells))
      }
    } else {
      results <- parallel::parLapply(cl, 1:n_cells, process_cell)
    }
  }

  end_time <- Sys.time()
  total_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

  errors <- sapply(results, function(x) !x$success)
  if (any(errors)) {
    error_indices <- which(errors)
    warning(sprintf("%d cells failed processing. Cell indices: %s",
                    length(error_indices),
                    paste(utils::head(error_indices, 10), collapse = ", ")))
    if (verbose) {
      for (idx in error_indices) {
        cat(sprintf("  Cell %d error: %s\n", idx, results[[idx]]$error))
      }
    }
  }

  # Extract results
  cn_calls_amplicon <- matrix(NA_integer_, nrow = n_amplicons, ncol = n_cells)
  mu_depth_by_chr <- vector("list", n_cells)
  var_depth_by_chr <- matrix(NA_real_, nrow = n_chrs, ncol = n_cells)
  sigma_maf_by_chr <- matrix(NA_real_, nrow = n_chrs, ncol = n_cells)
  epsilon_by_chr <- matrix(NA_real_, nrow = n_chrs, ncol = n_cells)
  logLik_by_chr <- matrix(NA_real_, nrow = n_chrs, ncol = n_cells)
  iter_by_chr <- matrix(NA_integer_, nrow = n_chrs, ncol = n_cells)
  converged_by_chr <- matrix(NA, nrow = n_chrs, ncol = n_cells)

  if (!is.null(colnames(depth_matrix))) {
    colnames(cn_calls_amplicon) <- colnames(depth_matrix)
    names(mu_depth_by_chr) <- colnames(depth_matrix)
    colnames(var_depth_by_chr) <- colnames(depth_matrix)
    colnames(sigma_maf_by_chr) <- colnames(depth_matrix)
    colnames(epsilon_by_chr) <- colnames(depth_matrix)
    colnames(logLik_by_chr) <- colnames(depth_matrix)
    colnames(iter_by_chr) <- colnames(depth_matrix)
    colnames(converged_by_chr) <- colnames(depth_matrix)
  }

  if (!is.null(rownames(depth_matrix))) {
    rownames(cn_calls_amplicon) <- rownames(depth_matrix)
  }

  rownames(var_depth_by_chr) <- as.character(unique_chrs)
  rownames(sigma_maf_by_chr) <- as.character(unique_chrs)
  rownames(epsilon_by_chr) <- as.character(unique_chrs)
  rownames(logLik_by_chr) <- as.character(unique_chrs)
  rownames(iter_by_chr) <- as.character(unique_chrs)
  rownames(converged_by_chr) <- as.character(unique_chrs)

  for (i in 1:n_cells) {
    if (results[[i]]$success) {
      res <- results[[i]]$result
      cn_calls_amplicon[, i] <- res$path
      mu_depth_by_chr[[i]] <- res$mu_depth_by_chr
      var_depth_by_chr[, i] <- res$var_depth_by_chr
      sigma_maf_by_chr[, i] <- res$sigma_maf_by_chr
      epsilon_by_chr[, i] <- res$epsilon_by_chr
      logLik_by_chr[, i] <- res$logLik_by_chr
      iter_by_chr[, i] <- res$iter_by_chr
      converged_by_chr[, i] <- res$converged_by_chr
    }
  }

  convergence_rate_per_cell <- colMeans(converged_by_chr, na.rm = TRUE)
  convergence_rate_per_cell[is.nan(convergence_rate_per_cell)] <- NA_real_
  mean_convergence_rate <- mean(convergence_rate_per_cell, na.rm = TRUE)
  if (is.nan(mean_convergence_rate)) {
    mean_convergence_rate <- NA_real_
  }

  if (progress) {
    cat(sprintf("\nParallel per-chromosome HMM analysis complete!\n"))
    cat(sprintf("  Total time: %.2f seconds\n", total_time))
    cat(sprintf("  Mean time per cell: %.3f seconds\n", total_time / n_cells))
    cat(sprintf("  Speedup estimate: %.1fx (vs sequential)\n", ncores * 0.85))
    cat(sprintf("  Successful: %d / %d cells\n", sum(!errors), n_cells))
    if (any(errors)) {
      cat(sprintf("  Failed: %d cells\n", sum(errors)))
    }
    cat(sprintf("  Mean convergence rate: %.1f%% (%d/%d chromosomes x %d cells)\n",
                100 * mean_convergence_rate,
                round(mean_convergence_rate * n_chrs), n_chrs, n_cells))
  }

  return(list(
    cn_calls_amplicon = cn_calls_amplicon,
    mu_depth_by_chr = mu_depth_by_chr,
    var_depth_by_chr = var_depth_by_chr,
    sigma_maf_by_chr = sigma_maf_by_chr,
    epsilon_by_chr = epsilon_by_chr,
    logLik_by_chr = logLik_by_chr,
    iter_by_chr = iter_by_chr,
    converged_by_chr = converged_by_chr,
    convergence_summary = list(
      convergence_rate_per_cell = convergence_rate_per_cell,
      mean_convergence_rate = mean_convergence_rate,
      total_chromosomes = n_chrs
    ),
    timing = list(
      total_time = total_time,
      mean_time_per_cell = total_time / n_cells,
      ncores_used = ncores
    )
  ))
}

