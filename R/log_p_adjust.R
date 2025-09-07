
#  Helper Functions -#

#' Methods supported by logp.adjust
#'
#' Character vector of method names accepted by \code{logp.adjust()}.
#' These are derived from p.adjust.methods and should match unless a new method has been added.
#' @seealso [p.adjust.methods()]
#' @format A character vector.
#' @examples
#' logp.adjust.methods
#' @export
logp.adjust.methods <- c("bonferroni","holm","hommel","hochberg","BH","BY","fdr","none")

#- Methods -#

.log_bonferroni <- function(p, n, base){
	pmin(0, log(n,base) + p)
}

.log_holm <- function(p, n, lp, base){
	i <- seq_len(lp)
	o <- order(p)
	ro <- order(o)
	pmin( 0, cummax( log((n + 1L - i), base) + p[o]))[ro]
}

.log_hommel <- function(p, n, lp, base){
  # pad to n with p=1 (log p = 0)
  if (n > lp) p <- c(p, rep.int(0, n - lp))

  # work in sorted order (ascending log p)
  o  <- order(p)
  ro <- order(o)
  ps <- p[o]

  i  <- seq_len(n)
  pa <- rep.int(min(log(n, base) + ps - log(i, base)), n)
  lq <- rep.int(-Inf, n)

  if (n >= 2L) {
    for (j in seq.int(n - 1L, 2L, by = -1L)) {      # explicit countdown
      ij <- seq_len(n - j + 1L)
      i2 <- seq.int(n - j + 2L, n)
      q1 <- min(log(j, base) + ps[i2] - log(seq.int(2L, j), base))
      lq[ij] <- pmin(log(j, base) + ps[ij], q1)
      lq[i2] <- lq[n - j + 1L]
      pa <- pmax(pa, lq)
    }
  }

  res_sorted <- pmax(pa, ps)
  res <- res_sorted[ro]
  if (n > lp) res <- res[seq_len(lp)]
  res
}

.log_hochberg <- function(p, n, lp, base){
	i <- lp:1L
	o <- order(p, decreasing = TRUE)
	ro <- order(o)
	pmin(0, cummin(log((n + 1L - i), base) + p[o]))[ro]
}

.log_BH <- function(p, n, lp, base){
	i <- lp:1L
	o <- order(p, decreasing = TRUE)
	ro <- order(o)
	pmin(0, cummin(log(n,base) - log(i,base) + p[o]))[ro]
}

.log_BY <- function(p, n, lp, base){
	i <- lp:1L
	o <- order(p, decreasing = TRUE)
	ro <- order(o)
	q <- sum(1/(1L:n))
	pmin(0, cummin(log(q,base) + log(n,base) - log(i,base) + p[o]))[ro]
}


#- Dispatcher -#


#' Adjust log P-values for Multiple Comparisons
#'
#' @param p vector of p-values in log space
#' @param method correction method to use
#' @param n number of comparisons, must be greater than or equal to length(p).
#' @param base the log base to perform operations in. Should match the base of the p-values. defaults to natural log.
#'
#' @importFrom stats setNames
#' @return vector of adjusted log p-values
#' @export
#'
logp.adjust <-
	function (p, method = logp.adjust.methods, n = length(p), base = exp(1)){
		# The original code and functions were adapted from the p.adjust function
		method <- match.arg(method)
		if (method == "fdr")
			method <- "BH"
		nm <- names(p)
		p <- as.numeric(p)
		p0 <- setNames(p, nm)
		if (all(nna <- !is.na(p)))
			nna <- TRUE
		else p <- p[nna]
		lp <- length(p)
		stopifnot(n >= lp)
		if (n <= 1)
			return(p0)
		if (n == 2 && method == "hommel")
			method <- "hochberg"

		# The original function defines all functions here
		# These have been defined as their own internal helper functions (above)
		p0[nna] <- switch(
			method,
			bonferroni = .log_bonferroni(p, n, base),
			holm = .log_holm(p, n, lp, base),
			hommel = .log_hommel(p, n, lp, base),
			hochberg = .log_hochberg(p, n, lp, base),
			BH = .log_BH(p, n, lp, base),
			BY = .log_BY(p, n, lp, base),
			none = p
		)
		p0
	}

















