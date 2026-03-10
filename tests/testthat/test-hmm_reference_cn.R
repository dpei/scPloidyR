test_that("default reference_cn=2 maintains backward compatibility", {
  set.seed(101)
  n_ref_cells <- 20
  n_amplicons <- 15

  ref_depth <- matrix(rnorm(n_amplicons * n_ref_cells, mean = 100, sd = 15),
                      nrow = n_amplicons, ncol = n_ref_cells)
  ref_baf_list <- lapply(1:n_ref_cells, function(cell) {
    lapply(1:n_amplicons, function(amp) rnorm(2, 0.5, 0.08))
  })

  params <- learn_reference_parameters(ref_depth, ref_baf_list)

  raw_medians <- apply(ref_depth, 1, median, na.rm = TRUE)
  expect_equal(params$diploid_depth_per_amplicon, raw_medians, tolerance = 1e-10)
  expect_true("reference_cn_per_amplicon" %in% names(params))
  expect_equal(params$reference_cn_per_amplicon, rep(2, n_amplicons))
})

test_that("scalar reference_cn correctly scales diploid baseline", {
  set.seed(202)
  n_ref_cells <- 20
  n_amplicons <- 10

  ref_depth <- matrix(rnorm(n_amplicons * n_ref_cells, mean = 150, sd = 20),
                      nrow = n_amplicons, ncol = n_ref_cells)
  ref_baf_list <- lapply(1:n_ref_cells, function(cell) {
    lapply(1:n_amplicons, function(amp) rnorm(2, 0.33, 0.08))
  })

  params <- learn_reference_parameters(ref_depth, ref_baf_list, reference_cn = 3)

  raw_medians <- apply(ref_depth, 1, median, na.rm = TRUE)
  expected_diploid <- raw_medians * (2 / 3)
  expect_equal(params$diploid_depth_per_amplicon, expected_diploid, tolerance = 1e-10)
  expect_equal(params$reference_cn_per_amplicon, rep(3, n_amplicons))
})

test_that("per-chromosome reference_cn expands correctly", {
  set.seed(303)
  n_ref_cells <- 15
  n_amplicons_per_chr <- 5
  n_chr <- 3
  n_amplicons <- n_amplicons_per_chr * n_chr
  amplicon_chr_labels <- rep(paste0("chr", 1:n_chr), each = n_amplicons_per_chr)

  ref_depth <- matrix(NA, nrow = n_amplicons, ncol = n_ref_cells)
  for (i in 1:n_ref_cells) {
    ref_depth[amplicon_chr_labels == "chr1", i] <- rnorm(n_amplicons_per_chr, 100, 10)
    ref_depth[amplicon_chr_labels == "chr2", i] <- rnorm(n_amplicons_per_chr, 150, 15)
    ref_depth[amplicon_chr_labels == "chr3", i] <- rnorm(n_amplicons_per_chr, 50, 8)
  }

  ref_baf_list <- lapply(1:n_ref_cells, function(cell) {
    lapply(1:n_amplicons, function(amp) rnorm(2, 0.5, 0.08))
  })

  reference_cn <- c("chr1" = 2, "chr2" = 3, "chr3" = 1)

  params <- learn_reference_parameters(
    ref_depth, ref_baf_list,
    reference_cn = reference_cn,
    amplicon_chr_labels = amplicon_chr_labels
  )

  expect_equal(params$reference_cn_per_amplicon[amplicon_chr_labels == "chr1"],
               rep(2, n_amplicons_per_chr))
  expect_equal(params$reference_cn_per_amplicon[amplicon_chr_labels == "chr2"],
               rep(3, n_amplicons_per_chr))
  expect_equal(params$reference_cn_per_amplicon[amplicon_chr_labels == "chr3"],
               rep(1, n_amplicons_per_chr))

  expect_true(all(params$diploid_depth_per_amplicon > 70 &
                  params$diploid_depth_per_amplicon < 130))
})

test_that("per-amplicon reference_cn vector works correctly", {
  set.seed(404)
  n_ref_cells <- 10
  n_amplicons <- 8

  reference_cn_per_amp <- c(2, 2, 3, 3, 1, 1, 4, 4)

  base_diploid <- 100
  ref_depth <- matrix(NA, nrow = n_amplicons, ncol = n_ref_cells)
  for (i in 1:n_ref_cells) {
    ref_depth[, i] <- rnorm(n_amplicons,
                            mean = base_diploid * reference_cn_per_amp / 2,
                            sd = 10)
  }

  ref_baf_list <- lapply(1:n_ref_cells, function(cell) {
    lapply(1:n_amplicons, function(amp) rnorm(2, 0.5, 0.08))
  })

  params <- learn_reference_parameters(ref_depth, ref_baf_list,
                                       reference_cn = reference_cn_per_amp)

  expect_equal(params$reference_cn_per_amplicon, reference_cn_per_amp)
  expect_true(all(abs(params$diploid_depth_per_amplicon - 100) < 30))
})

test_that("invalid reference_cn values produce appropriate errors", {
  set.seed(505)
  n_ref_cells <- 5
  n_amplicons <- 6
  ref_depth <- matrix(100, nrow = n_amplicons, ncol = n_ref_cells)
  ref_baf_list <- lapply(1:n_ref_cells, function(cell) {
    lapply(1:n_amplicons, function(amp) c(0.5))
  })

  expect_error(
    learn_reference_parameters(ref_depth, ref_baf_list, reference_cn = 0),
    "reference_cn"
  )
  expect_error(
    learn_reference_parameters(ref_depth, ref_baf_list, reference_cn = -1),
    "reference_cn"
  )
  expect_error(
    learn_reference_parameters(ref_depth, ref_baf_list, reference_cn = c(2, 2, 3)),
    "reference_cn"
  )
  expect_error(
    learn_reference_parameters(ref_depth, ref_baf_list, reference_cn = c("chr1" = 2)),
    "amplicon_chr_labels"
  )
})

test_that("unspecified chromosomes default to CN=2", {
  set.seed(606)
  n_ref_cells <- 10
  n_amplicons <- 12
  amplicon_chr_labels <- rep(c("chr1", "chr2", "chr3"), each = 4)
  ref_depth <- matrix(rnorm(n_amplicons * n_ref_cells, 100, 15),
                      nrow = n_amplicons, ncol = n_ref_cells)
  ref_baf_list <- lapply(1:n_ref_cells, function(cell) {
    lapply(1:n_amplicons, function(amp) c(0.5))
  })

  reference_cn <- c("chr2" = 3)

  params <- learn_reference_parameters(
    ref_depth, ref_baf_list,
    reference_cn = reference_cn,
    amplicon_chr_labels = amplicon_chr_labels
  )

  expect_equal(params$reference_cn_per_amplicon[amplicon_chr_labels == "chr1"], rep(2, 4))
  expect_equal(params$reference_cn_per_amplicon[amplicon_chr_labels == "chr2"], rep(3, 4))
  expect_equal(params$reference_cn_per_amplicon[amplicon_chr_labels == "chr3"], rep(2, 4))
})

test_that("reference_cn integrates with call_copy_number", {
  set.seed(707)
  n_ref <- 15
  n_test <- 10
  n_total <- n_ref + n_test
  n_amplicons_per_chr <- 4
  n_chr <- 3
  n_amplicons <- n_amplicons_per_chr * n_chr
  amplicon_chr_labels <- rep(paste0("chr", 1:n_chr), each = n_amplicons_per_chr)

  depth_matrix <- matrix(NA, nrow = n_amplicons, ncol = n_total)
  for (i in 1:n_ref) {
    depth_matrix[amplicon_chr_labels == "chr1", i] <- rnorm(n_amplicons_per_chr, 100, 10)
    depth_matrix[amplicon_chr_labels == "chr2", i] <- rnorm(n_amplicons_per_chr, 150, 15)
    depth_matrix[amplicon_chr_labels == "chr3", i] <- rnorm(n_amplicons_per_chr, 100, 10)
  }
  for (i in (n_ref + 1):n_total) {
    depth_matrix[, i] <- rnorm(n_amplicons, 100, 12)
  }

  baf_list <- lapply(1:n_total, function(cell) {
    lapply(1:n_amplicons, function(amp) rnorm(2, 0.5, 0.08))
  })

  reference_cn <- c("chr2" = 3)

  result <- call_copy_number(
    depth_matrix = depth_matrix,
    baf_data = baf_list,
    amplicon_chr_labels = amplicon_chr_labels,
    reference_cell_indices = 1:n_ref,
    reference_cn = reference_cn,
    states = 1:5,
    ncores = 1,
    progress = FALSE
  )

  test_calls <- result$cn_calls_amplicon[, (n_ref + 1):n_total]
  cn2_rate <- mean(test_calls == 2, na.rm = TRUE)

  expect_true(cn2_rate > 0.7,
              info = sprintf("Test cells CN=2 rate: %.2f (expected > 0.7)", cn2_rate))
})
