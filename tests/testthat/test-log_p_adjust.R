
testthat::test_that(
	"linear outputs functionally equivalent to p.adjust",
	{
		set.seed(12345)
		vals <- runif(100)
		for( m in logp.adjust.methods){
			testthat::expect_equal(
				exp(logp.adjust(log(vals),method = m)),
				p.adjust(vals,method = m),
				tolerance = 1e-12
			)
		}
	}
)

testthat::test_that(
	"NA values are handled properly",
	{
		vals <- c(0.005,0.1,NA,0.01)
		for( m in logp.adjust.methods){
			testthat::expect_equal(
				exp(logp.adjust(log(vals),method = m)),
				p.adjust(vals,method = m),
				tolerance = 1e-12
			)
		}
	}
)

testthat::test_that(
	"Custom `n` values are handled properly",
	{
		vals <- c(0.005,0.1,NA,0.01)
		n <- 10
		for( m in logp.adjust.methods){
			testthat::expect_equal(
				exp(logp.adjust(log(vals),method = m,n = n)),
				p.adjust(vals,method = m,n = n),
				tolerance = 1e-12
			)
		}
	}
)

testthat::test_that(
	"Different log spaces handled appropriately",
	{
		vals <- c(0.005,0.1,NA,0.01)
		bases <- 2:10
		for( m in logp.adjust.methods){
			for( base in bases ){
				testthat::expect_equal(
					base^(logp.adjust(log(vals,base = base),method = m,base = base)),
					p.adjust(vals,method = m),
					tolerance = 1e-12
				)
			}
		}
	}
)





