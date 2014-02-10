\name{plot.localgauss.indtest}
     \alias{plot.localgauss.indtest}
     
     \title{Local Gaussian correlation based independence test plot}
     \description{
     Plots results of local Gaussian correlation based independence test.  
     }
     \usage{\method{plot}{localgauss.indtest}(x,...,plot.points=F,poscol="magenta",
     negcol="cyan",point.col="black",point.size=NULL,xlab="",ylab="")
     }
     \arguments{
       \item{x}{S3 object of class "localgauss.indtest" produced by the \code{localgauss.indtest}-function}
       \item{...}{Not used.}
       \item{plot.points}{If TRUE, the original observations are overlain.}
       \item{poscol}{Color indicating +1 test result}
       \item{negcol}{Color indicating -1 test result}
       \item{point.col}{The colour used for observations points if \code{plot.points} is TRUE.}
       \item{point.size}{The size of observations points if \code{plot.points} 
       is TRUE.}
       \item{xlab,ylab}{The label of x-axis and y-axis, respectively.}
     }
     
    \references{Geir Drage Berentsen, Tore Selland Kleppe, Dag Tjostheim,
        Introducing localgauss, an R Package for Estimating and Visualizing Local Gaussian Correlation, Journal of Statistical Software, 56(12), 1-18, 2014, (http://www.jstatsoft.org/v56/i12/). See also Tjoestheim, D. and Hufthammer K. O., 
Local Gaussian correlation: A new measure of dependence, Journal of Econometrics,
172(1),pages 33-48,2013, for a detailed description of local Gaussian correlation and  Berentsen, G.D. and Tjoestheim D., Recognizing and visualizing departures from independence in bivariate data using local Gaussian correlation,
http://people.uib.no/gbe062/local-gaussian-correlation/ for a description of the local independence test. }


     \seealso{
       \code{\link{localgauss}}.
     }
     \examples{
     x=rnorm(n=100)
     y=x^2 + rnorm(n=100)
     lgobj = localgauss(x,y,gsize=8)
     lgind = localgauss.indtest(lgobj)
     plot(lgind)
    }
    
     
