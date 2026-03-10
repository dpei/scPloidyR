#' @title HMM Genotype Utilities
#' @description Functions for genotype priors, MAF expectations, and heterozygosity estimation.
#' @name genotypes
#' @keywords internal
NULL

# Hard-coded homozygous and heterozygous genotype templates per CN state.
HOMOZYGOUS_PRIORS <- list(
  `2` = c(AA = 0.5, BB = 0.5),
  `3` = c(AAA = 0.5, BBB = 0.5),
  `4` = c(AAAA = 0.5, BBBB = 0.5),
  `5` = c(AAAAA = 0.5, BBBBB = 0.5)
)

normalize_prior <- function(x) {
  x / sum(x)
}

HETEROZYGOUS_PRIORS <- list(
  `2` = c(AB = 1.0),
  `3` = normalize_prior(c(AAB = 0.40, ABB = 0.40)),
  `4` = normalize_prior(c(AAAB = 0.30, AABB = 0.30, ABBB = 0.30)),
  `5` = normalize_prior(c(AAAAB = 0.23, AAABB = 0.27,
                          AABBB = 0.27, ABBBB = 0.23))
)

#' @title Get Genotype Priors for Copy Number State
#' @description Returns biologically-informed genotype priors for a given CN state,
#'   adapted to observed heterozygosity rate in the data via interpolation
#'   between low-het (uniform) and high-het (historical) templates.
#' @param cn_state Integer, the copy number state (1-5).
#' @param het_rate Numeric, estimated heterozygosity rate from data (0 to 1).
#' @return Named numeric vector of genotype probabilities.
#' @keywords internal
get_genotype_prior <- function(cn_state, het_rate = 0.7) {
  if (cn_state == 1) {
    return(c(A = 0.5, B = 0.5))
  }

  cn_key <- as.character(cn_state)
  homo <- HOMOZYGOUS_PRIORS[[cn_key]]
  heter <- HETEROZYGOUS_PRIORS[[cn_key]]

  if (is.null(homo) || is.null(heter)) {
    return(c(uniform = 1.0))
  }

  het_rate <- max(0, min(1, het_rate))
  het_weight <- het_rate
  homo_weight <- 1 - het_weight

  weights <- c(homo_weight * homo, het_weight * heter)
  total <- sum(weights)

  if (total <= 0) {
    return(c(uniform = 1.0))
  }

  weights / total
}

#' @title Calculate Expected MAF for Genotype
#' @description Calculates the expected minor allele frequency for a given genotype.
#' @param genotype Character, genotype string (e.g., "AA", "AB", "AABB").
#' @return Numeric, expected MAF (between 0 and 0.5).
#' @keywords internal
get_expected_maf <- function(genotype) {
  alleles <- strsplit(genotype, "")[[1]]
  n_A <- sum(alleles == "A")
  n_B <- sum(alleles == "B")
  total <- n_A + n_B

  if (total == 0) {
    return(0)
  }

  min(n_A, n_B) / total
}

#' @title Estimate Heterozygosity Rate from BAF Data
#' @description Estimates the heterozygosity rate by counting variants with MAF > threshold.
#' @param maf_list List of numeric vectors, where each element contains MAF values for one amplicon.
#' @param threshold Numeric, MAF threshold to classify variants as heterozygous.
#' @return Numeric, estimated heterozygosity rate (0 to 1).
#' @keywords internal
estimate_heterozygosity <- function(maf_list, threshold = HETEROZYGOSITY_THRESHOLD) {
  all_mafs <- unlist(maf_list)
  all_mafs <- all_mafs[!is.na(all_mafs)]

  if (length(all_mafs) == 0) {
    return(0)
  }

  het_count <- sum(all_mafs > threshold)
  total_count <- length(all_mafs)

  het_count / total_count
}
