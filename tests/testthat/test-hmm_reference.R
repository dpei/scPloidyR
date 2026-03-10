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
  n_ref_cells <- 50
  n_amplicons <- 10

  set.seed(456)
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

test_that("call_copy_number works with reference cells", {
  set.seed(131415)

  n_ref <- 20
  n_test <- 10
  n_total <- n_ref + n_test
  n_amplicons <- 10

  depth_matrix <- matrix(NA, nrow = n_amplicons, ncol = n_total)
  depth_matrix[, 1:n_ref] <- matrix(rnorm(n_amplicons * n_ref, 100, 20), nrow = n_amplicons)
  depth_matrix[, (n_ref + 1):n_total] <- matrix(rnorm(n_amplicons * n_test, 150, 25), nrow = n_amplicons)

  baf_data <- lapply(1:n_total, function(cell) {
    lapply(1:n_amplicons, function(amp) {
      if (cell <= n_ref) {
        rnorm(2, 0.5, 0.08)
      } else {
        rnorm(2, 0.33, 0.08)
      }
    })
  })

  amplicon_chr_labels <- rep(1L, n_amplicons)

  result <- call_copy_number(
    depth_matrix = depth_matrix,
    baf_data = baf_data,
    amplicon_chr_labels = amplicon_chr_labels,
    reference_cell_indices = 1:n_ref,
    states = 1:5,
    progress = FALSE
  )

  expect_true(is.matrix(result$cn_calls_amplicon))
  expect_equal(nrow(result$cn_calls_amplicon), n_amplicons)
  expect_equal(ncol(result$cn_calls_amplicon), n_total)

  ref_calls <- result$cn_calls_amplicon[, 1:n_ref]
  expect_true(mean(ref_calls == 2) > 0.6)

  test_calls <- result$cn_calls_amplicon[, (n_ref + 1):n_total]
  expect_true(mean(test_calls) > mean(ref_calls))
})

test_that("functions handle errors appropriately", {
  expect_error(
    fit_hmm(depth = rep(100, 10), baf = lapply(1:10, function(i) c(0.5))),
    "reference_params"
  )

  expect_error(
    call_copy_number(
      depth_matrix = matrix(100, 10, 5),
      baf_data = lapply(1:5, function(i) lapply(1:10, function(j) c(0.5))),
      amplicon_chr_labels = rep(1L, 10)
    ),
    "reference_cell_indices"
  )

  expect_error(
    call_copy_number(
      depth_matrix = matrix(100, 10, 5),
      baf_data = lapply(1:5, function(i) lapply(1:10, function(j) c(0.5))),
      amplicon_chr_labels = rep(1L, 10),
      reference_cell_indices = integer(0)
    ),
    "reference_cell_indices"
  )

  expect_error(
    call_copy_number(
      depth_matrix = matrix(100, 10, 5),
      baf_data = lapply(1:5, function(i) lapply(1:10, function(j) c(0.5))),
      reference_cell_indices = 1:3
    ),
    "amplicon_chr_labels"
  )
})

test_that("reference-aware HMM runs successfully on aneuploid data", {
  set.seed(161718)

  n_ref <- 30
  n_test <- 10
  n_total <- n_ref + n_test
  n_amplicons <- 50

  depth_matrix <- matrix(NA, nrow = n_amplicons, ncol = n_total)
  depth_matrix[, 1:n_ref] <- matrix(rnorm(n_amplicons * n_ref, 100, 10), nrow = n_amplicons)
  depth_matrix[, (n_ref + 1):n_total] <- matrix(rnorm(n_amplicons * n_test, 200, 20), nrow = n_amplicons)

  baf_data <- lapply(1:n_total, function(cell) {
    lapply(1:n_amplicons, function(amp) {
      if (cell <= n_ref) {
        rnorm(5, 0.5, 0.06)
      } else {
        rnorm(5, 0.5, 0.06)
      }
    })
  })

  amplicon_chr_labels <- rep(1L, n_amplicons)

  result <- call_copy_number(
    depth_matrix = depth_matrix,
    baf_data = baf_data,
    amplicon_chr_labels = amplicon_chr_labels,
    reference_cell_indices = 1:n_ref,
    states = 1:5,
    progress = FALSE
  )

  expect_true(is.matrix(result$cn_calls_amplicon))
  expect_equal(dim(result$cn_calls_amplicon), c(n_amplicons, n_total))

  calculate_mode <- function(x) {
    ux <- unique(x[!is.na(x)])
    if (length(ux) == 0) return(NA)
    ux[which.max(tabulate(match(x, ux)))]
  }

  cn_by_cell <- apply(result$cn_calls_amplicon, 2, calculate_mode)
  ref_accuracy <- mean(cn_by_cell[1:n_ref] == 2)

  expect_true(ref_accuracy > 0.8,
              info = "Reference cells should be accurately called as CN=2")
})

test_that("parallel vs sequential produce identical results", {
  set.seed(202501)

  n_ref <- 15
  n_test <- 10
  n_total <- n_ref + n_test
  n_amplicons <- 12

  depth_matrix <- matrix(NA, nrow = n_amplicons, ncol = n_total)
  for (i in 1:n_ref) {
    depth_matrix[, i] <- rnorm(n_amplicons, mean = 100, sd = 15)
  }
  for (i in (n_ref + 1):n_total) {
    depth_matrix[, i] <- rnorm(n_amplicons, mean = 150, sd = 20)
  }

  baf_data <- lapply(1:n_total, function(cell) {
    lapply(1:n_amplicons, function(amp) {
      if (cell <= n_ref) {
        rnorm(3, mean = 0.5, sd = 0.08)
      } else {
        rnorm(3, mean = 0.33, sd = 0.08)
      }
    })
  })

  amplicon_chr_labels <- rep(1L, n_amplicons)

  result_sequential <- call_copy_number(
    depth_matrix = depth_matrix,
    baf_data = baf_data,
    amplicon_chr_labels = amplicon_chr_labels,
    reference_cell_indices = 1:n_ref,
    states = 1:5,
    ncores = 1,
    progress = FALSE,
    verbose = FALSE
  )

  result_parallel <- call_copy_number(
    depth_matrix = depth_matrix,
    baf_data = baf_data,
    amplicon_chr_labels = amplicon_chr_labels,
    reference_cell_indices = 1:n_ref,
    states = 1:5,
    ncores = 2,
    progress = FALSE,
    verbose = FALSE
  )

  expect_equal(dim(result_sequential$cn_calls_amplicon),
               dim(result_parallel$cn_calls_amplicon))

  cn_calls_match <- all(result_sequential$cn_calls_amplicon == result_parallel$cn_calls_amplicon,
                        na.rm = TRUE)
  expect_true(cn_calls_match,
              info = "Sequential and parallel execution should produce identical CN calls")

  expect_equal(result_sequential$var_depth_by_chr, result_parallel$var_depth_by_chr)
  expect_equal(result_sequential$sigma_maf_by_chr, result_parallel$sigma_maf_by_chr)
  expect_equal(result_sequential$epsilon_by_chr, result_parallel$epsilon_by_chr)
})

test_that("multiple parallel runs with same inputs are deterministic", {
  set.seed(202502)

  n_ref <- 10
  n_test <- 5
  n_total <- n_ref + n_test
  n_amplicons <- 8

  depth_matrix <- matrix(rnorm(n_amplicons * n_total, mean = 100, sd = 20),
                         nrow = n_amplicons, ncol = n_total)

  baf_data <- lapply(1:n_total, function(cell) {
    lapply(1:n_amplicons, function(amp) {
      rnorm(2, mean = 0.5, sd = 0.08)
    })
  })

  amplicon_chr_labels <- rep(1L, n_amplicons)

  result_run1 <- call_copy_number(
    depth_matrix = depth_matrix,
    baf_data = baf_data,
    amplicon_chr_labels = amplicon_chr_labels,
    reference_cell_indices = 1:n_ref,
    states = 1:5,
    ncores = 2,
    progress = FALSE
  )

  result_run2 <- call_copy_number(
    depth_matrix = depth_matrix,
    baf_data = baf_data,
    amplicon_chr_labels = amplicon_chr_labels,
    reference_cell_indices = 1:n_ref,
    states = 1:5,
    ncores = 2,
    progress = FALSE
  )

  expect_equal(result_run1$cn_calls_amplicon, result_run2$cn_calls_amplicon,
               info = "Multiple runs with same inputs should produce identical results")
})
