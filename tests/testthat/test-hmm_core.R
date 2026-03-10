test_that("learn_reference_parameters returns correct structure", {
  set.seed(123)
  n_ref_cells <- 10
  n_amplicons <- 20

  amplicon_baselines <- rnorm(n_amplicons, mean = 100, sd = 20)
  ref_depth <- matrix(NA, nrow = n_amplicons, ncol = n_ref_cells)
  for (i in 1:n_ref_cells) {
    ref_depth[, i] <- amplicon_baselines * rnorm(n_amplicons, mean = 1, sd = 0.2)
  }

  ref_baf_list <- lapply(1:n_ref_cells, function(cell) {
    lapply(1:n_amplicons, function(amp) rnorm(3, mean = 0.5, sd = 0.08))
  })

  params <- learn_reference_parameters(ref_depth, ref_baf_list)

  expect_true(is.list(params))
  expect_true("diploid_depth_per_amplicon" %in% names(params))
  expect_true("depth_variance" %in% names(params))
  expect_true("sigma_maf_init" %in% names(params))
  expect_true("heterozygosity_rate_by_amplicon" %in% names(params))
  expect_equal(length(params$diploid_depth_per_amplicon), n_amplicons)
  expect_true(all(params$diploid_depth_per_amplicon > 0))
  expect_true(params$depth_variance > 0)
  expect_true(params$sigma_maf_init > 0.01 & params$sigma_maf_init < 0.2)
  expect_equal(length(params$heterozygosity_rate_by_amplicon), n_amplicons)
  expect_true(all(params$heterozygosity_rate_by_amplicon >= 0))
  expect_true(all(params$heterozygosity_rate_by_amplicon <= 1))
})

test_that("reference depth baseline is accurate", {
  set.seed(456)
  n_ref_cells <- 50
  n_amplicons <- 10

  true_amplicon_depths <- seq(80, 120, length.out = n_amplicons)
  ref_depth <- matrix(NA, nrow = n_amplicons, ncol = n_ref_cells)
  for (i in 1:n_ref_cells) {
    ref_depth[, i] <- true_amplicon_depths * rnorm(n_amplicons, mean = 1, sd = 0.1)
  }

  ref_baf_list <- lapply(1:n_ref_cells, function(cell) {
    lapply(1:n_amplicons, function(amp) rnorm(2, 0.5, 0.05))
  })

  params <- learn_reference_parameters(ref_depth, ref_baf_list)

  relative_error <- abs(params$diploid_depth_per_amplicon - true_amplicon_depths) / true_amplicon_depths
  expect_true(all(relative_error < 0.15))
})

test_that("per-amplicon heterozygosity estimation is accurate", {
  set.seed(789)
  ref_depth <- matrix(100, nrow = 10, ncol = 5)

  ref_baf_high_het <- lapply(1:5, function(cell) {
    lapply(1:10, function(amp) rnorm(5, 0.5, 0.05))
  })
  params_high_het <- learn_reference_parameters(ref_depth, ref_baf_high_het)
  expect_true(all(params_high_het$heterozygosity_rate_by_amplicon > 0.8))

  ref_baf_low_het <- lapply(1:5, function(cell) {
    lapply(1:10, function(amp) rnorm(5, 0, 0.02))
  })
  params_low_het <- learn_reference_parameters(ref_depth, ref_baf_low_het)
  expect_true(all(params_low_het$heterozygosity_rate_by_amplicon < 0.2))
})

test_that("fit_hmm works with reference parameters", {
  set.seed(101112)

  n_amplicons <- 22
  ref_params <- list(
    diploid_depth_per_amplicon = rep(100, n_amplicons),
    depth_variance = 400,
    sigma_maf_init = 0.08,
    heterozygosity_rate_by_amplicon = rep(0.7, n_amplicons)
  )

  depth_vec <- rnorm(n_amplicons, mean = 100, sd = 20)
  baf_list <- lapply(1:n_amplicons, function(i) rnorm(2, 0.5, 0.08))

  result <- fit_hmm(
    depth = depth_vec,
    baf = baf_list,
    reference_params = ref_params,
    states = 1:5
  )

  expect_true(is.list(result))
  expect_true("path" %in% names(result))
  expect_true("mu_depth" %in% names(result))
  expect_equal(length(result$path), n_amplicons)
  expect_true(sum(result$path == 2) > n_amplicons * 0.7)
})

test_that("fit_hmm requires reference_params", {
  expect_error(
    fit_hmm(depth = rep(100, 10), baf = lapply(1:10, function(i) c(0.5))),
    "reference_params"
  )
})
