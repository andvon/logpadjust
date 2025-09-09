
# logpadjust

<!-- badges: start -->
[![R-CMD-check](https://github.com/andvon/logpadjust/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/andvon/logpadjust/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

This package is currently meant to provide access to the function `logp.adjust`,
which is a modified version of `stats::p.adjust` that accepts log-space p-values and returns adjusted p-values in log-space.

While `p.adjust` could be used in the majority of cases by converting values to linear space, certain methods may generate p-values that exceed the minimum value (~1e-300).
This package was written to provide a function that would enable multiple test correction methods on these values without losing information.

## Installation

You can install the development version of logpadjust from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("andvon/logpadjust")
```

## Example

The `logp.adjust` function should be able to be called and used in the same manner as `stats::p.adjust`.

``` r
library(logpadjust)
log_p_values <- log(c(0.005,0.1))
logp.adjust(p = log_p_values,method = "bonferroni")
```

The methods (which should match `stats::p.adjust.methods`) can be viewed with `logp.adjust.methods`.

If p-values were calculated in a space other than the natural log, specifying `base` will update the calculations.

```
log10_p_values <- log10(c(0.005,0.1))
logp.adjust(p = log10_p_values,method = "bonferroni", base = 10)
```

For additional help, please see the package's documentation.




