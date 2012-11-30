\name{localgauss}
\alias{localgauss}
\title{local Gaussian parameters}
\description{
  Routine for estimating local Gaussian parameters based on a sample
  from the bivariate distribution under consideration. The routine can either
  estimate local parameters on a grid covering the data controlled by the \code{gsize} and \code{hthresh}
  parameters. Otherwise, local Gaussian parameters can be estimated at coordinates
  specified by the user in \code{xy.mat}. 
}
\usage{localgauss(x,y,b1=1,b2=1,gsize=15,hthresh=0.001,xy.mat=NULL)}
\arguments{
  \item{x,y}{The two data vectors}
  \item{b1,b2}{The bandwidth in the \code{x}-direction and \code{y}-direction, respectively}
  \item{gsize}{The gridsize (only used if \code{xy.mat} is not specified). }
  \item{hthresh}{Gridpoints where a non-parametric density estimate is lower than hthresh are omitted (only used if \code{xy.mat} is not specified). }
  \item{xy.mat}{A M times 2 matrix of points where the local parameters are to be estimated.}

  }
\details{The objective function is maximized using a modified Newton
  method. The user should check whether the field eflag in the
  returned object is zero for all estimates. If not,
  the optimizer has not converged and the estimates should not be trusted.
  For more details, see [Reference to article].
}
\value{
  \item{par.est}{M times 5 matrix of parameter estimates, with columns
    mu1,mu2,sigma1,sigma2,rho.}
  \item{eflag}{M-vector of exitflags from the optimizer. Estimations
    with exit flags other than 0 should not be trusted.}
  \item{hessian}{The negative Hessian of the objective function.}
}

 \references{See Tjoestheim, D. and Hufthammer K. O., 
Local Gaussian correlation: A new measure of dependence, Journal of Econometrics,
172(1),pages 33-48,2013, for a detailed description of local Gaussian correlation.}



 \seealso{
       \code{\link{plot.localgauss}}.
     }

\examples{
    x=rnorm(n=1000)
    y=x^2 + rnorm(n=1000)
    lgobj = localgauss(x,y)
    plot(lgobj)
}

\keyword{localgauss}

