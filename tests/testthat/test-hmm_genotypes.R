test_that("low heterozygosity favors balanced homozygotes", {
  priors <- get_genotype_prior(2, het_rate = 0)

  expect_equal(unname(priors["AA"]), unname(priors["BB"]))
  expect_equal(unname(priors["AA"]), 0.5, tolerance = 1e-8)
  expect_equal(unname(priors["AB"]), 0, tolerance = 1e-8)
  expect_equal(sum(priors[c("AA", "BB")]), 1, tolerance = 1e-8)
})

test_that("full heterozygosity recovers hetero templates", {
  priors <- get_genotype_prior(4, het_rate = 1)

  expect_equal(unname(priors["AAAB"]), 1 / 3, tolerance = 1e-8)
  expect_equal(unname(priors["AABB"]), 1 / 3, tolerance = 1e-8)
  expect_equal(unname(priors["ABBB"]), 1 / 3, tolerance = 1e-8)
  expect_equal(unname(priors["AAAA"]), 0)
  expect_equal(unname(priors["BBBB"]), 0)
})

test_that("heterozygous mass increases with het_rate", {
  low <- get_genotype_prior(3, het_rate = 0)
  mid <- get_genotype_prior(3, het_rate = 0.5)
  high <- get_genotype_prior(3, het_rate = 1)

  expect_equal(unname(low["AAB"]), 0)
  expect_lt(unname(low["AAB"]), unname(mid["AAB"]))
  expect_lt(unname(mid["AAB"]), unname(high["AAB"]))
  expect_equal(unname(high["AAA"]), 0)
  expect_equal(unname(high["BBB"]), 0)
})

test_that("copy number one stays unbiased", {
  priors <- get_genotype_prior(1, het_rate = 0.3)

  expect_equal(unname(priors["A"]), 0.5)
  expect_equal(unname(priors["B"]), 0.5)
})

test_that("get_expected_maf computes correctly", {
  expect_equal(get_expected_maf("AB"), 0.5)
  expect_equal(get_expected_maf("AA"), 0)
  expect_equal(get_expected_maf("AAB"), 1 / 3)
  expect_equal(get_expected_maf("AABB"), 0.5)
  expect_equal(get_expected_maf("A"), 0)
})
