test_that("Depth-only HMM calls CN=2 when normalized counts ~1.1", {
  set.seed(123)

  n_amplicons <- 60
  depth_sd <- 0.07
  depth <- rnorm(n_amplicons, mean = 1.1, sd = depth_sd)

  baf <- replicate(n_amplicons, numeric(0), simplify = FALSE)

  reference_params <- list(
    diploid_depth_per_amplicon = rep(1.0, n_amplicons),
    depth_variance = depth_sd^2,
    sigma_maf_init = 0.08,
    heterozygosity_rate_by_amplicon = rep(0, n_amplicons)
  )

  fit <- fit_hmm(
    depth = depth,
    baf = baf,
    reference_params = reference_params,
    states = 1:5,
    max_iter = 50,
    tol = 1e-6,
    eps_init = 1e-3
  )

  expect_true(is.list(fit))
  expect_equal(length(fit$path), n_amplicons)
  prop_cn2 <- mean(fit$path == 2)
  expect_gt(prop_cn2, 0.95)

  tab_states <- sort(table(fit$path), decreasing = TRUE)
  expect_equal(as.integer(names(tab_states)[1]), 2)
})

test_that("Low-het depth-only HMM does not collapse into CN=1", {
  set.seed(999)

  n_amplicons <- 88
  depth <- rep(100, n_amplicons)
  baf <- replicate(n_amplicons, numeric(0), simplify = FALSE)

  reference_params <- list(
    diploid_depth_per_amplicon = rep(100, n_amplicons),
    depth_variance = 0.01,
    sigma_maf_init = 0.08,
    heterozygosity_rate_by_amplicon = rep(0, n_amplicons)
  )

  fit <- fit_hmm(
    depth = depth,
    baf = baf,
    reference_params = reference_params,
    states = 1:5,
    max_iter = 50,
    tol = 1e-6,
    eps_init = 1e-3
  )

  prop_cn2 <- mean(fit$path == 2)
  prop_cn1 <- mean(fit$path == 1)
  expect_gt(prop_cn2, 0.8)
  expect_lt(prop_cn1, 0.2)
})
