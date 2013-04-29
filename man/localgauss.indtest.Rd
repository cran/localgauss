\name{localgauss.indtest}
\alias{localgauss.indtest}
\title{Pointwise Independence test based on local Gaussian correlation}
\description{
  Routine for testing for local independence based on local Gaussian parameters. It accepts an S3 object produced by \code{localgauss()}, and perfoms a bootstrap-based test with null-hypothesis being that \code{x} and \code{y} are indpendent.  
}
\usage{localgauss.indtest(locobj,R=10,alpha=0.10,seed=1)}
\arguments{
  \item{locobj}{\code{localgauss}-object}
  \item{R}{Number of bootstrap replica}
  \item{alpha}{significance level (note: two sided test)}
  \item{seed}{Random seed in used for bootstrap}
 }
\details{The test is based on producing a null-distribution of local Gaussian correlations were the original data are resampled from their empirical marginal distributions. The bootstrap-based null-distribution is produced for each point specified in \code{xy.mat} in \code{locobj}. An estimated local correlation for the original data significantly larger than the null-distribution is indicated with +1 (returned in the vector \code{test.results}). An estimated local correlation for the original data insignifcant with respect to the null-distribution is indicated with 0. An estimated local correlation for the original datasignificantly smaller than the null-distribution is indicated with -1.
}
\value{S3 object of type \code{localgauss.indtest} containing the fields:
  \item{localgauss}{simply returns \code{locobj}.}
  \item{upper}{Vector containing the 1-alpha/2 quantiles of the null-distributions.}
  \item{lower}{Vector containing the alpha/2 quantiles of the null-distributions.}
  \item{test.results}{Vector containing the test results.}
}

 \references{See Tjoestheim, D. and Hufthammer K. O., 
Local Gaussian correlation: A new measure of dependence, Journal of Econometrics,
172(1),pages 33-48,2013, for a detailed description of local Gaussian correlation
 and Berentsen, G.D. and Tjoestheim D., Recognizing and visualizing departures from independence in bivariate data using local Gaussian correlation,
http://people.uib.no/gbe062/local-gaussian-correlation/ for a description of the local independence test. 
}



 \seealso{
       \code{\link{plot.localgauss.indtest}}.
     }

\examples{
    x=rnorm(n=100)
    y=x^2 + rnorm(n=100)
    lgobj = localgauss(x,y,gsize=8)
    lgind = localgauss.indtest(lgobj)
    plot(lgind)
}

\keyword{localgauss}

