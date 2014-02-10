\name{plot.localgauss}
     \alias{plot.localgauss}
     
     \title{Local Gaussian correlation plot}
     \description{
     Plots estimates of local Gaussian correlation.  
     }
     \usage{\method{plot}{localgauss}(x,...,plot.text=TRUE,plot.points=FALSE,tsize=3,
  lowcol="cyan",highcol="magenta",point.col="black",
  point.size=NULL,xlab="",ylab="",divergent.col.grad=T)
     }
     \arguments{
       \item{x}{S3 object of class "localgauss" produced by the \code{localgauss}-function}
       \item{...}{Not used.}
       \item{plot.text}{If TRUE, the numerical values of the estimated local correlation are added to each tile. }
       \item{plot.points}{If TRUE, the original observations are overlain.}
       \item{tsize}{The font size used if plot.text is TRUE}
       \item{lowcol}{The color used to indicate negative correlation of -1}
       \item{highcol}{The color used to indicate positive correlation of 1}
       \item{point.col}{The colour used for observations points if \code{plot.points} is TRUE.}
       \item{point.size}{The size of observations points if \code{plot.points} 
       is TRUE.}
       \item{xlab,ylab}{The label of x-axis and y-axis, respectively.}
       \item{divergent.col.grad}{If TRUE, a divergent color gradient between lowcol and highcol with 0 as midpoint is used. 
       If FALSE a ordinary color gradient between lowcol and highcol is used.}
     }
     
    \references{Geir Drage Berentsen, Tore Selland Kleppe, Dag Tjostheim,
        Introducing localgauss, an R Package for Estimating and Visualizing Local Gaussian Correlation, Journal of Statistical Software, 56(12), 1-18, 2014, (http://www.jstatsoft.org/v56/i12/). See also Tjoestheim, D. and Hufthammer K. O., 
Local Gaussian correlation: A new measure of dependence, Journal of Econometrics,
172(1),pages 33-48,2013, for a detailed description of local Gaussian correlation.}


     \seealso{
       \code{\link{localgauss}}.
     }
     \examples{
    x=rnorm(n=1000)
    y=x^2 + rnorm(n=1000)
    lgobj = localgauss(x,y)
    plot(lgobj)
    }
    
     
