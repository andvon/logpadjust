
#  Helper Functions -#

#' log function factory helper
#'
#' @description
#' Helper function meant to be used internally.
#' Returns a function that represents \code{log(x, base = base)}
#'
#' @param base numeric value to set the base of the log function to. default: exp(1)
#'
#' @return function
#' @export
#'
#' @examples
#' custom_log10 <- .log_factory(10)
#' custom_log10(10)
#'
.log_factory <- function(base = exp(1)){
	if(!is.numeric(base)) stop("The `base` supplied to .log_factory must be numeric.")
	func <- function(x){log(x, base)}
	return(func)
}

#' Methods supported by logp.adjust
#'
#' Character vector of method names accepted by \code{logp.adjust()}.
#' These are derived from p.adjust.methods and should match unless a new method has been added.
#' @seealso [p.adjust.methods]
#' @format A character vector.
#' @examples
#' logp.adjust.methods
#' @export
logp.adjust.methods <- c("bonferroni","holm","hommel","hochberg","BH","BY","fdr","none")

#- Methods -#

#' Multiple Test Correction (Log-Space)
#'
#' These are internal functions that implement multiple test correction methods in **log-space**.
#' They accept log-transformed p-values and return adjusted p-values in
#' log-space using different correction methods.
#'
#' @details
#' The following methods are provided:
#' \itemize{
#'   \item \code{.log_holm} — Holm (1979)
#'   \item \code{.log_hochberg} — Hochberg (1988)
#'   \item \code{.log_hommel} — Hommel (1988)
#'   \item \code{.log_BH} — Benjamini & Hochberg (1995)
#'   \item \code{.log_bonferroni} — Bonferroni correction
#'   \item \code{.log_BY} — Benjamini & Yekutieli (2001)
#' }
#'
#' @param p Numeric vector of p-values in **log-space**.
#' @param n Integer. Total number of tests used for correction.
#' @param base Integer. Logarithmic base to use for transformations.
#' @param lp Integer. Effective length of the p-value vector (Not used in Bonferroni).
#'
#' @return A numeric vector of **adjusted p-values in log-space**.
#'
#' @seealso [stats::p.adjust] for the linear-scale equivalents.
#'
#' @examples
#' pvals <- log(c(0.01, 0.02, 0.05))
#'
#' # Bonferroni
#' .log_bonferroni(pvals, n = length(pvals), base = exp(1))
#'
#' # Other methods (all use the same implementation)
#' .log_BH(pvals, n = length(pvals), lp = length(pvals), base = exp(1))
#'
#' @name log_p_adjust
NULL

#' @rdname log_p_adjust
#' @export
.log_bonferroni <- function(p, n, base){
	lfunc <- .log_factory(base)
	pmin(0, lfunc(n) + p)
}

#' @rdname log_p_adjust
#' @export
.log_holm <- function(p, n, lp, base){
	lfunc <- .log_factory(base)
	i <- seq_len(lp)
	o <- order(p)
	ro <- order(o)
	pmin( 0, cummax( lfunc(n + 1L - i) + p[o]))[ro]
}

#' @rdname log_p_adjust
#' @export
.log_hommel <- function(p, n, lp, base){
	lfunc <- .log_factory(base)
	if (n > lp) p <- c(p, rep.int(0, n - lp))
	i <- seq_len(n)
	o <- order(p)
	p <- p[o]
	ro <- order(o)
	q <- pa <- rep.int(min(lfunc(n) + p - lfunc(i)), n)
	for (j in (n - 1L):2L) {
		ij <- seq_len(n - j + 1L)
		i2 <- (n - j + 2L):n
		q1 <- min(lfunc(j) + p[i2] - lfunc(2L:j))
		q[ij] <- pmin(lfunc(j) + p[ij], q1)
		q[i2] <- q[n - j + 1L]
		pa <- pmax(pa, q)
	}
	pmax(pa, p)[if (lp < n) ro[1L:lp] else ro]
}

#' @rdname log_p_adjust
#' @export
.log_hochberg <- function(p, n, lp, base){
	lfunc <- .log_factory(base)
	i <- lp:1L
	o <- order(p, decreasing = TRUE)
	ro <- order(o)
	pmin(0, cummin(lfunc(n + 1L - i) + p[o]))[ro]
}

#' @rdname log_p_adjust
#' @export
.log_BH <- function(p, n, lp, base){
	lfunc <- .log_factory(base)
	i <- lp:1L
	o <- order(p, decreasing = TRUE)
	ro <- order(o)
	pmin(0, cummin(lfunc(n) - lfunc(i) + p[o]))[ro]
}

#' @rdname log_p_adjust
#' @export
.log_BY <- function(p, n, lp, base){
	lfunc <- .log_factory(base)
	i <- lp:1L
	o <- order(p, decreasing = TRUE)
	ro <- order(o)
	q <- sum(1/(1L:n))
	pmin(0, cummin(lfunc(q) + lfunc(n) - lfunc(i) + p[o]))[ro]
}


#- Dispatcher -#


#' Adjust log P-values for Multiple Comparisons
#'
#' @description
#'
#' \code{logp.adjust()} modifies the code from \code{\link[stats]{p.adjust}} to
#' operate on log-space p-values. Each method should should generate identical results
#' to the original function when converted back to linear.
#'
#' @details
#' This was mainly written for cases when one is working with very small p-values (< 1e-300).
#' In these cases, the value is rounded to 0. When calculating p-values using a log distribution,
#' it is possible to get values that can be lower than this cap. This function should allow users to apply
#' multiple test correction methods to these values without losing potentially valuable information.
#'
#' The original code from \code{\link[stats]{p.adjust}} was used as a starting point and thus borrows heavily from the theory
#' and efforts put into that codebase. Please keep this in mind and cite the original authors!
#' The primary difference is that methods were extracted into their own helper functions for easier interpretation.
#'
#' For the methods themselves, please see the documentation for the helper functions.
#'
#' @references
#' The references here are taken from \code{\link[stats]{p.adjust}}.
#'
#' Benjamini, Y., and Hochberg, Y. (1995).Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society Series B, 57, 289–300.
#' \doi{10.1111/j.2517-6161.1995.tb02031.x}
#'
#' Benjamini, Y., and Yekutieli, D. (2001). The control of the false discovery rate in multiple testing under dependency. Annals of Statistics, 29, 1165–1188.
#' \doi{10.1214/aos/1013699998}
#'
#' Holm, S. (1979). A simple sequentially rejective multiple test procedure. Scandinavian Journal of Statistics, 6, 65–70.
#' \url{https://www.jstor.org/stable/4615733}
#'
#' Hommel, G. (1988). A stagewise rejective multiple test procedure based on a modified Bonferroni test. Biometrika, 75, 383–386.
#' \doi{10.2307/2336190}
#'
#' Hochberg, Y. (1988). A sharper Bonferroni procedure for multiple tests of significance. Biometrika, 75, 800–803.
#' \doi{10.2307/2336325}
#'
#'
#' R Core Team (2023). _R: A Language and Environment for Statistical Computing_. R Foundation for Statistical Computing, Vienna, Austria.
#' \url{https://www.R-project.org/}
#'
#' @param p vector of p-values in log space
#' @param method correction method to use. See: \code{\link{logp.adjust.methods}}
#' @param n number of comparisons, must be greater than or equal to length(p).
#' @param base the log base to perform operations in. Should match the base of the p-values. defaults to natural log.
#' @seealso \code{\link[stats]{p.adjust}}, \code{\link{logp.adjust.methods}}
#' @importFrom stats setNames
#' @return vector of adjusted log p-values
#' @export
#'
#' @examples
#' # Operate on log p-values --------------
#'
#' p <- c(0.001, 0.1, NA)
#' logp.adjust(p = log(p), method = "bonferroni")
#' # Equivalent to
#' log(p.adjust(p, method = "bonferroni"))
#'
#' # Use a different base -----------------
#' logp.adjust(p = log(p, base = 2), method = "bonferroni", base = 2)
#'
#' # Convert back to linear p-values
#' exp( logp.adjust(p = log(p), method = "bonferroni") )
#'
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

















