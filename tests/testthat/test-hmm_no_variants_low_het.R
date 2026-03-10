test_that("learn_reference_parameters returns het_rate=0 when no variants observed", {
  set.seed(20251216)

  n_amplicons <- 40
  n_ref_cells <- 5

  depth_baseline <- runif(n_amplicons, min = 80, max = 120)
  reference_depth_matrix <- vapply(seq_len(n_ref_cells), function(i) {
    depth_baseline * rnorm(n_amplicons, mean = 1.0, sd = 0.03)
  }, numeric(n_amplicons))
  reference_depth_matrix <- as.matrix(reference_depth_matrix)

  reference_baf_list <- lapply(seq_len(n_ref_cells), function(i) {
    replicate(n_amplicons, numeric(0), simplify = FALSE)
  })

  ref <- suppressWarnings(learn_reference_parameters(
    reference_depth_matrix = reference_depth_matrix,
    reference_baf_list = reference_baf_list
  ))

  expect_equal(length(ref$heterozygosity_rate_by_amplicon), n_amplicons)
  expect_true(all(ref$heterozygosity_rate_by_amplicon == 0))
  expect_lte(mean(ref$heterozygosity_rate_by_amplicon), LOW_HET_RATE_THRESHOLD)
})
