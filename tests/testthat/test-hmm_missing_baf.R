test_that("fit_hmm treats NA-only BAF vectors as missing (depth-only emissions)", {
  set.seed(123)

  n_amplicons <- 60
  depth_sd <- 0.07
  depth <- rnorm(n_amplicons, mean = 1.1, sd = depth_sd)

  baf_empty <- replicate(n_amplicons, numeric(0), simplify = FALSE)
  baf_na_only <- replicate(n_amplicons, rep(NA_real_, 3), simplify = FALSE)

  reference_params <- list(
    diploid_depth_per_amplicon = rep(1.0, n_amplicons),
    depth_variance = depth_sd^2,
    sigma_maf_init = 0.08,
    heterozygosity_rate_by_amplicon = rep(0, n_amplicons)
  )

  fit_empty <- fit_hmm(
    depth = depth,
    baf = baf_empty,
    reference_params = reference_params,
    states = 1:5,
    max_iter = 50,
    tol = 1e-6,
    eps_init = 1e-3
  )
  fit_na_only <- fit_hmm(
    depth = depth,
    baf = baf_na_only,
    reference_params = reference_params,
    states = 1:5,
    max_iter = 50,
    tol = 1e-6,
    eps_init = 1e-3
  )

  expect_equal(fit_na_only$path, fit_empty$path)
  expect_gt(mean(fit_na_only$path == 2), 0.95)
})

test_that("fit_hmm_parallel ignores NA BAF matrix entries (no penalty vs depth-only)", {
  set.seed(202)

  n_amplicons <- 50
  n_cells <- 3
  reference_cell_indices <- 1L

  amplicon_baseline <- runif(n_amplicons, min = 80, max = 120)
  depth_matrix <- vapply(1:n_cells, function(i) {
    amplicon_baseline * rnorm(n_amplicons, mean = 1.0, sd = 0.02)
  }, numeric(n_amplicons))
  depth_matrix <- as.matrix(depth_matrix)
  colnames(depth_matrix) <- paste0("cell_", 1:n_cells)
  rownames(depth_matrix) <- paste0("amp_", 1:n_amplicons)

  baf_matrix_na <- matrix(NA_real_, nrow = n_amplicons, ncol = n_cells)
  baf_list_empty <- lapply(1:n_cells, function(i) {
    replicate(n_amplicons, numeric(0), simplify = FALSE)
  })

  amplicon_chr_labels <- rep(1L, n_amplicons)

  res_na <- suppressWarnings(call_copy_number(
    depth_matrix = depth_matrix,
    baf_data = baf_matrix_na,
    amplicon_chr_labels = amplicon_chr_labels,
    reference_cell_indices = reference_cell_indices,
    states = 1:5,
    max_iter = 50,
    tol = 1e-6,
    eps_init = 1e-3,
    ncores = 1,
    progress = FALSE
  ))
  res_empty <- suppressWarnings(call_copy_number(
    depth_matrix = depth_matrix,
    baf_data = baf_list_empty,
    amplicon_chr_labels = amplicon_chr_labels,
    reference_cell_indices = reference_cell_indices,
    states = 1:5,
    max_iter = 50,
    tol = 1e-6,
    eps_init = 1e-3,
    ncores = 1,
    progress = FALSE
  ))

  expect_true(is.matrix(res_na$cn_calls_amplicon))
  expect_equal(dim(res_na$cn_calls_amplicon), c(n_amplicons, n_cells))
  expect_equal(res_na$cn_calls_amplicon, res_empty$cn_calls_amplicon)

  prop_cn2_by_cell <- colMeans(res_na$cn_calls_amplicon == 2)
  expect_true(all(prop_cn2_by_cell > 0.9))
})

test_that("Mixed BAF/NA-only BAF uses depth where BAF is missing", {
  set.seed(404)

  n_amplicons <- 60
  depth_sd <- 0.07
  depth <- rnorm(n_amplicons, mean = 1.1, sd = depth_sd)

  baf_mixed <- vector("list", n_amplicons)
  baf_mixed_zeroed <- vector("list", n_amplicons)
  for (i in 1:n_amplicons) {
    if (i <= n_amplicons / 2) {
      baf_mixed[[i]] <- rnorm(3, mean = 0.5, sd = 0.03)
      baf_mixed_zeroed[[i]] <- baf_mixed[[i]]
    } else {
      baf_mixed[[i]] <- rep(NA_real_, 3)
      baf_mixed_zeroed[[i]] <- numeric(0)
    }
  }

  reference_params <- list(
    diploid_depth_per_amplicon = rep(1.0, n_amplicons),
    depth_variance = depth_sd^2,
    sigma_maf_init = 0.08,
    heterozygosity_rate_by_amplicon = rep(0, n_amplicons)
  )

  fit_mixed <- fit_hmm(
    depth = depth,
    baf = baf_mixed,
    reference_params = reference_params,
    states = 1:5,
    max_iter = 50,
    tol = 1e-6,
    eps_init = 1e-3
  )
  fit_mixed_zeroed <- fit_hmm(
    depth = depth,
    baf = baf_mixed_zeroed,
    reference_params = reference_params,
    states = 1:5,
    max_iter = 50,
    tol = 1e-6,
    eps_init = 1e-3
  )

  expect_equal(fit_mixed$path, fit_mixed_zeroed$path)
  expect_gt(mean(fit_mixed$path == 2), 0.95)
})
