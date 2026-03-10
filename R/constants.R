#' @title HMM Default Parameters
#' @description Default parameters for HMM algorithm tuning.
#' @keywords internal
HMM_DEFAULTS <- list(
  MAX_ITER = 50,
  TOLERANCE = 1e-5,
  EPS_INIT = 1e-3,
  MU_PRIOR_WEIGHT_LOW_HET = 10
)

#' @title Default HMM Core Count
#' @description Default number of CPU cores for HMM per-cell parallelism.
#' @keywords internal
DEFAULT_HMM_CORES <- 12

#' @title Heterozygosity Threshold
#' @description MAF threshold for het/hom classification.
#' @keywords internal
HETEROZYGOSITY_THRESHOLD <- 0.3

#' @title Minimum MAF Sigma
#' @description Lower bound for MAF standard deviation.
#' @keywords internal
MIN_SIGMA_MAF <- 0.05

#' @title Maximum MAF Sigma
#' @description Upper bound for MAF standard deviation.
#' @keywords internal
MAX_SIGMA_MAF <- 0.15

#' @title Low Heterozygosity Rate Threshold
#' @description Threshold for considering a sample "low heterozygosity" when deciding
#'   whether to apply depth-anchored HMM overrides.
#' @keywords internal
LOW_HET_RATE_THRESHOLD <- 0.25
