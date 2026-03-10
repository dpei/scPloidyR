test_that("fit_hmm_per_chromosome works for single cell", {
  n_chrs <- 3
  n_amplicons_per_chr <- 4
  n_amplicons <- n_chrs * n_amplicons_per_chr
  chr_labels <- rep(1:n_chrs, each = n_amplicons_per_chr)

  set.seed(123)
  depth <- rnorm(n_amplicons, mean = 100, sd = 10)
  baf <- lapply(1:n_amplicons, function(i) rbeta(3, 5, 5))

  reference_params <- list(
    diploid_depth_per_amplicon = rep(100, n_amplicons),
    depth_variance = 0.1,
    sigma_maf_init = 0.08,
    heterozygosity_rate_by_amplicon = rep(1.0, n_amplicons)
  )

  result <- fit_hmm_per_chromosome(
    depth = depth,
    baf = baf,
    amplicon_chr_labels = chr_labels,
    reference_params = reference_params,
    states = 1:5,
    max_iter = 20
  )

  expect_equal(length(result$path), n_amplicons)
  expect_equal(length(result$mu_depth_by_chr), n_chrs)
  expect_true(all(result$converged_by_chr))
  expect_true(all(result$path %in% 1:5))
  expect_true(mean(result$path == 2) > 0.5)
})

test_that("per-chromosome mode respects chromosome boundaries", {
  n_chrs <- 3
  n_amplicons_per_chr <- 6
  n_amplicons <- n_chrs * n_amplicons_per_chr
  chr_labels <- rep(1:n_chrs, each = n_amplicons_per_chr)

  set.seed(456)
  depth <- c(
    rnorm(n_amplicons_per_chr, mean = 40, sd = 3),
    rnorm(n_amplicons_per_chr, mean = 100, sd = 3),
    rnorm(n_amplicons_per_chr, mean = 200, sd = 5)
  )

  baf <- lapply(1:n_amplicons, function(i) numeric(0))

  reference_params <- list(
    diploid_depth_per_amplicon = rep(100, n_amplicons),
    depth_variance = 0.05,
    sigma_maf_init = 0.08,
    heterozygosity_rate_by_amplicon = rep(1.0, n_amplicons)
  )

  result <- fit_hmm_per_chromosome(
    depth = depth,
    baf = baf,
    amplicon_chr_labels = chr_labels,
    reference_params = reference_params,
    states = 1:5,
    max_iter = 50
  )

  expect_true(all(result$converged_by_chr))
  expect_true(all(result$path %in% 1:5))
})

test_that("call_copy_number works for multiple cells", {
  n_chrs <- 2
  n_amplicons_per_chr <- 3
  n_amplicons <- n_chrs * n_amplicons_per_chr
  n_cells <- 4
  chr_labels <- rep(1:n_chrs, each = n_amplicons_per_chr)

  set.seed(789)
  depth_matrix <- matrix(
    rnorm(n_amplicons * n_cells, mean = 100, sd = 10),
    nrow = n_amplicons,
    ncol = n_cells
  )
  colnames(depth_matrix) <- paste0("cell_", 1:n_cells)
  rownames(depth_matrix) <- paste0("amplicon_", 1:n_amplicons)

  baf_list <- lapply(1:n_cells, function(cell) {
    lapply(1:n_amplicons, function(amp) numeric(0))
  })

  result <- call_copy_number(
    depth_matrix = depth_matrix,
    baf_data = baf_list,
    amplicon_chr_labels = chr_labels,
    reference_cell_indices = c(1, 2),
    states = 1:5,
    max_iter = 20,
    ncores = 1,
    progress = FALSE
  )

  expect_equal(dim(result$cn_calls_amplicon), c(n_amplicons, n_cells))
  expect_equal(length(result$mu_depth_by_chr), n_cells)
  expect_equal(dim(result$var_depth_by_chr), c(n_chrs, n_cells))

  convergence_rate <- mean(result$converged_by_chr, na.rm = TRUE)
  expect_true(convergence_rate > 0.5)
})

test_that("convergence summary is robust to cell failures", {
  old_fit <- fit_hmm_per_chromosome
  on.exit({ fit_hmm_per_chromosome <<- old_fit }, add = TRUE)

  fit_hmm_per_chromosome <<- function(...) {
    args <- list(...)
    depth <- args$depth
    if (!is.null(depth) && mean(depth, na.rm = TRUE) > 150) {
      stop("forced per-cell failure")
    }
    old_fit(...)
  }

  n_chrs <- 2
  n_amplicons_per_chr <- 3
  n_amplicons <- n_chrs * n_amplicons_per_chr
  n_cells <- 3
  chr_labels <- rep(1:n_chrs, each = n_amplicons_per_chr)

  set.seed(42)
  depth_matrix <- matrix(rnorm(n_amplicons * n_cells, mean = 100, sd = 5),
                         nrow = n_amplicons, ncol = n_cells)
  depth_matrix[, 3] <- depth_matrix[, 3] + 200

  baf_list <- lapply(1:n_cells, function(cell) {
    lapply(1:n_amplicons, function(amp) numeric(0))
  })

  result <- suppressWarnings(call_copy_number(
    depth_matrix = depth_matrix,
    baf_data = baf_list,
    amplicon_chr_labels = chr_labels,
    reference_cell_indices = c(1, 2),
    states = 1:5,
    max_iter = 5,
    ncores = 1,
    progress = FALSE
  ))

  expect_false(is.nan(result$convergence_summary$mean_convergence_rate))
})

test_that("convergence failures produce valid output", {
  n_chrs <- 2
  n_amplicons_per_chr <- 2
  n_amplicons <- n_chrs * n_amplicons_per_chr
  chr_labels <- rep(1:n_chrs, each = n_amplicons_per_chr)

  depth <- rep(100, n_amplicons)
  baf <- lapply(1:n_amplicons, function(i) numeric(0))

  reference_params <- list(
    diploid_depth_per_amplicon = rep(100, n_amplicons),
    depth_variance = 0.1,
    sigma_maf_init = 0.08,
    heterozygosity_rate_by_amplicon = rep(1.0, n_amplicons)
  )

  result <- suppressWarnings(
    fit_hmm_per_chromosome(
      depth = depth,
      baf = baf,
      amplicon_chr_labels = chr_labels,
      reference_params = reference_params,
      states = 1:5,
      max_iter = 5
    )
  )

  expect_equal(length(result$path), n_amplicons)
  expect_equal(length(result$converged_by_chr), n_chrs)
})

test_that("reference parameters correctly subset per chromosome", {
  n_chrs <- 2
  n_amplicons_per_chr <- 3
  n_amplicons <- n_chrs * n_amplicons_per_chr
  chr_labels <- rep(1:n_chrs, each = n_amplicons_per_chr)

  set.seed(202)
  depth <- rnorm(n_amplicons, mean = 100, sd = 10)
  baf <- lapply(1:n_amplicons, function(i) rbeta(2, 5, 5))

  diploid_depths <- seq(50, 150, length.out = n_amplicons)

  reference_params <- list(
    diploid_depth_per_amplicon = diploid_depths,
    depth_variance = 0.1,
    sigma_maf_init = 0.08,
    heterozygosity_rate_by_amplicon = rep(1.0, n_amplicons)
  )

  result <- fit_hmm_per_chromosome(
    depth = depth,
    baf = baf,
    amplicon_chr_labels = chr_labels,
    reference_params = reference_params,
    states = 1:5,
    max_iter = 20
  )

  expect_true(all(result$converged_by_chr))
  expect_equal(length(result$path), n_amplicons)
})

test_that("output is compatible with aggregate_copy_number", {
  n_chrs <- 2
  n_amplicons_per_chr <- 4
  n_amplicons <- n_chrs * n_amplicons_per_chr
  chr_labels <- rep(1:n_chrs, each = n_amplicons_per_chr)

  set.seed(101)
  depth <- rnorm(n_amplicons, mean = 100, sd = 10)
  baf <- lapply(1:n_amplicons, function(i) rbeta(2, 5, 5))

  reference_params <- list(
    diploid_depth_per_amplicon = rep(100, n_amplicons),
    depth_variance = 0.1,
    sigma_maf_init = 0.08,
    heterozygosity_rate_by_amplicon = rep(1.0, n_amplicons)
  )

  result <- fit_hmm_per_chromosome(
    depth = depth,
    baf = baf,
    amplicon_chr_labels = chr_labels,
    reference_params = reference_params,
    states = 1:5,
    max_iter = 20
  )

  mock_metadata <- data.frame(
    chr = factor(chr_labels),
    amplicon.id = paste0("amp_", 1:n_amplicons)
  )

  path_matrix <- matrix(result$path, ncol = 1)
  colnames(path_matrix) <- "cell_1"
  rownames(path_matrix) <- mock_metadata$amplicon.id

  chr_calls <- aggregate_copy_number(path_matrix, mock_metadata)

  expect_equal(nrow(chr_calls), n_chrs)
  expect_equal(ncol(chr_calls), 1)
})
