#######################dependences###########################

# requrired by the plot-function
library(ggplot2)
library(MASS)

######################localgauss###################################
   
  	# x,y:               Datavectors
    # b1,b2;             Bandwidth in x and y direction    	
  	# gsize              Gridsize (used only if xy.mat is NULL)
    # hthresh            Leave out points where the kernel density estimate is lower than this.
    #                    (used only if xy.mat is NULL)
    # xy.mat             Points where the local Gaussian correlation is computed.
    #                    If xy.mat is not NULL, this option overrides the gsize and hthresh options




localgauss <- function(x,y,b1=1,b2=1,gsize=15,hthresh=0.001,xy.mat=NULL){

# work out if default or custom points (i.e. xy.mat is given) are to be used in
# esitmation
if(is.null(xy.mat)){
  #Create xy.mat by selecting points with density>hthresh	
	xy.est=kde2d(x, y, n=gsize, h=c(bandwidth.nrd(x), bandwidth.nrd(y)), lims=c(min(x),max(x),min(y),max(y)))  #density estimation
	xy.est=con2tr(xy.est)                                                                                      #convert list to dataframe
	xy.est=subset(xy.est,xy.est[,3]>hthresh)                                                                   #remove points with to low density
	xy.mat.internal=as.matrix(xy.est[,1:2])     
} else {
  # xy.mat is provided
  xy.mat.internal=xy.mat
}
# handle data in format suitable for the fortran function
n = length(x)
if(is.matrix(xy.mat.internal)){
  nrep = dim(xy.mat.internal)[1];
  x0 = xy.mat.internal[,1];
  y0 = xy.mat.internal[,2];
}
else{
  nrep = 1;
  x0 = as.vector(xy.mat.internal[1]);
  y0 = as.vector(xy.mat.internal[2]);
}
# do not alter function arguments
hx = abs(b1);
hy = abs(b2);
# the fortran function takes separate bandwidths for each row in
# xy.mat
if(length(hx)==1){
  hxx = rep(hx,nrep);
}
else{
  hxx = as.vector(hx);
}

if(length(hy)==1){
  hyy = rep(hy,nrep);
}
else{
  hyy = as.vector(hy);
}
# the fortran function takes an additional argument - rhof -
# that is not used by this package. Essentially, the fortran function
# has the ability to estimate means and std.devs. while keeping correlation
# fixed at some value rhof. Passing rhof that is not a valid correlation causes
# rho to be estimated. However we hide this functionality as it has not
# been tested thoroughly and have not proven useful in practice, and  
# thus a vector of -2.0s is passed to cause rho always to be estimated.

rhofix = rep(-2.0,nrep);

# allocate space for output
res.in = matrix(0.0,nrow=nrep,ncol=7)
hess.in = array(0.0,dim=c(nrep,5,5))
# actual call to the fortran function
f.out <- .Fortran("localgauss",
                  as.integer(n),
                  as.double(as.vector(x)),
                  as.double(as.vector(y)),
                  as.integer(nrep),
                  as.double(as.vector(x0)),
                  as.double(as.vector(y0)),
                  as.double(as.vector(hxx)),
                  as.double(as.vector(hyy)),
		  as.double(as.vector(rhofix)),
                  as.double(res.in),
		  as.double(hess.in),PACKAGE="localgauss")


# format output for easy return
res.out = as.matrix(matrix((f.out[10])[[1]],nrow=nrep,ncol=7));
hess.out = array(f.out[11][[1]],dim=c(nrep,5,5));
par.est= matrix(res.out[,1:5],nrow=nrep,ncol=5);
colnames(par.est) = c("mu_1","mu_2","sig_1","sig_2","rho");
eflag = res.out[,7]
ret = list(x=x,y=y,xy.mat=xy.mat.internal,b1=hx,b2=hy,par.est=par.est,eflag=eflag,hessian=hess.out)
attr(ret,"class") <- "localgauss"
return(ret);
}

######################plot.localgauss###################################

# x           :"localgauss" object from function localgauss
# plot.text   : Should the numerical local correlation value be printed on each tile?
# plot.points : Should the original observation points be overlain?
# tsize       : Font size used if plot.text is true.
# lowcol      : Colour used to indicate negative correlation of −1.
# highcol     : Colour used to indicate positive correlation of +1.
# point.col   : Color used for observations points (if plot.points is true).
# point.size  : Size of observations points (if plot.points is true).
# xlab        : label of x-axis.
# ylab        : label of y-axis.

plot.localgauss<-function(x,..., plot.text=TRUE,plot.points=FALSE,
tsize=3,lowcol="cyan",highcol="magenta",point.col="black",
point.size=NULL,xlab="",ylab=""){
  # to conform with generic/method consistency
  lgobj = x
	# Dummys to trick R CMD check
	x.obs=y.obs=NULL
	
	#Objects to plot.
	dat.obs = data.frame(x.obs=lgobj$x,y.obs=lgobj$y)
	rho = lgobj$par.est[,5]
	x = lgobj$xy.mat[,1]
	y = lgobj$xy.mat[,2]
	
	#Format rho with leading sign.
	rho.label = formatC(rho,digits=2, format="f",flag="+")                      
	#Dataframe for plotting rho.
	d = data.frame(x=x,y=y,rho=rho,rho.label=rho.label,stringsAsFactors=FALSE)  
	
	#Make a ggplot using the geom "tile".
	g = ggplot() + layer(data=d, mapping=aes(x=x, y=y, fill=rho, label=rho.label),geom="tile") # Plot estimated rho.
	# Add a new colour gradient to graph with possible different colours.
	gr = scale_fill_gradient2(midpoint=0,low=lowcol, high=highcol, space="Lab", limits=c(-1,1),breaks=seq(-1, 1, by=0.2))
	g = g + gr 
	
	#Should we overlay the original observations?
	if(plot.points) 
		if(is.null(point.size))
			g = g + layer(data=dat.obs, mapping=aes(x=x.obs,y=y.obs), geom="point", colour=point.col)
		else
			g = g + layer(data=dat.obs, mapping=aes(x=x.obs,y=y.obs), geom="point", colour=point.col,size=point.size)
	
	#Add numerical values of the local correlation to the graph if wanted
	if(plot.text)
		g = g + layer(data=d, mapping=aes(x=x, y=y, label=rho.label),geom="text",size=tsize)
	
	#Add labels
	g = g + scale_x_continuous(name=xlab)
	g = g + scale_y_continuous(name=ylab)
		
	# Return plot
	
	return(g)
	
}


