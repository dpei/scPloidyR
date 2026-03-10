test_that("Low-het depth-only HMM can enter/exit a short CN=1 block", {
  set.seed(20251216)

  n_amplicons <- 88
  baseline <- 100
  cn2_depth <- baseline
  cn1_depth <- baseline / 2
  noise_sd <- 2

  depth <- rep(cn2_depth, n_amplicons) + rnorm(n_amplicons, mean = 0, sd = noise_sd)
  deletion_idx <- 41:44
  depth[deletion_idx] <- cn1_depth + rnorm(length(deletion_idx), mean = 0, sd = noise_sd)

  baf <- replicate(n_amplicons, numeric(0), simplify = FALSE)

  reference_params <- list(
    diploid_depth_per_amplicon = rep(baseline, n_amplicons),
    depth_variance = (noise_sd / baseline)^2,
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

  prop_cn1_in_block <- mean(fit$path[deletion_idx] == 1)
  prop_cn2_outside <- mean(fit$path[-deletion_idx] == 2)
  prop_cn1_outside <- mean(fit$path[-deletion_idx] == 1)

  expect_gt(prop_cn1_in_block, 0.75)
  expect_gt(prop_cn2_outside, 0.9)
  expect_lt(prop_cn1_outside, 0.1)
})
