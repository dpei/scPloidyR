#' @title Log-Density of Normal Distribution
#' @description Calculates the log-density of a normal distribution.
#' @param x Numeric vector of observations.
#' @param mu Numeric, mean parameter.
#' @param var Numeric, variance parameter.
#' @return Numeric vector of log-density values.
#' @keywords internal
log_norm <- function(x, mu, var) stats::dnorm(x, mu, sqrt(var), log = TRUE)

#' @title Log-Density of Truncated Normal Distribution
#' @description Calculates the log-density of a truncated normal distribution on [a, b].
#' @param x Numeric vector of observations.
#' @param mu Numeric, mean parameter.
#' @param sigma Numeric, standard deviation parameter.
#' @param a Numeric, lower bound. Defaults to 0.
#' @param b Numeric, upper bound. Defaults to 0.5.
#' @return Numeric vector of log-density values.
#' @keywords internal
log_truncnorm <- function(x, mu, sigma, a = 0, b = 0.5) {
  log_density <- stats::dnorm(x, mu, sigma, log = TRUE)
  log_normalizer <- log(stats::pnorm(b, mu, sigma) - stats::pnorm(a, mu, sigma))
  log_density - log_normalizer
}

#' @title Log-Sum-Exp Function
#' @description Computes the log-sum-exp in a numerically stable way.
#' @param v Numeric vector of log values.
#' @return Numeric, the log of the sum of exponentials.
#' @keywords internal
logsumexp <- function(v) {
  m <- max(v)
  m + log(sum(exp(v - m)))
}

