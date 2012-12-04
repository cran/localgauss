! Module of global variables for easy passing of
! data to the objective function
module globals
implicit none
integer :: gn
real(8),dimension(:),allocatable :: gx
real(8),dimension(:),allocatable :: gy
real(8),dimension(:),allocatable :: gwts
real(8) :: gx0
real(8) :: gy0
real(8) :: ghx
real(8) :: ghy
logical :: cv
real(8) :: fixrho
end module globals

! module implementing the modified newton solver
module newton_solver
implicit none
contains
! This subroutine modifies (if neccesary) the hessian matrix
! and then invert it
subroutine mod_and_inv(n,mat,eflag)
implicit none
integer,intent(in) :: n
real(8),dimension(n,n),intent(inout) :: mat
integer,intent(out) :: eflag
real(8),dimension(n,n) :: matcp,tmp
real(8),dimension(n) :: ww
real(8),dimension(32*n) :: work
integer :: ii

matcp = mat

! eigenvalue-decomposition from LAPACK
call dsyev('V','L',n,matcp,n,ww,work,32*n,eflag)


! product with diagonal matrix
! notice absolute value to produce positive
! definite matrix
if(abs(eflag)==0) then
   do ii=1,n
      tmp(:,ii) = matcp(:,ii)/max(1.0e-12_8,abs(ww(ii)));
   enddo
   ! complete product
   mat = matmul(tmp,transpose(matcp));
endif
end subroutine mod_and_inv


! This is the actual modified newton solver
subroutine n_solver(lfun,dlfun,ddlfun,n,x0,xout,fout,gout,hout,eflag)
use globals
implicit none
! interfaces for function, gradient and hessian
interface 
subroutine lfun(p,pp,ll,eflag)
implicit none
integer,intent(in) :: p
real(8),dimension(p),intent(in) :: pp
real(8),intent(out) :: ll
integer,intent(out) :: eflag
end subroutine lfun
end interface
interface
subroutine dlfun(p,pp,dll,eflag)
implicit none
integer,intent(in) :: p
real(8),dimension(p),intent(in) :: pp
real(8),dimension(p),intent(out) :: dll
integer,intent(out) :: eflag
end subroutine dlfun
end interface
interface
subroutine ddlfun(p,pp,ll,dll,ddll,eflag)
implicit none
integer,intent(in) :: p
real(8),dimension(p),intent(in) :: pp
real(8),intent(out) :: ll
real(8),dimension(p),intent(out) :: dll
real(8),dimension(p,p),intent(out) :: ddll
integer,intent(out) :: eflag
end subroutine ddlfun
end interface
! convergence criteria
real(8),parameter :: TOLX = 1.0e-10_8;
real(8),parameter :: TOLGRAD = 1.0e-8;
real(8),parameter :: c1 = 1.0e-8_8;
integer,parameter :: MAXITER = 200
logical,parameter :: debug = .false.
integer,intent(in) :: n ! length of free variable
real(8),dimension(n),intent(in) :: x0 ! inital condition
real(8),dimension(n),intent(out) :: xout ! final estimate
real(8),intent(out) :: fout ! objective at estimate
real(8),dimension(n),intent(out) :: gout ! gradient at estimate
real(8),dimension(n,n),intent(out) :: hout ! hessian at estimate
integer,intent(out) :: eflag ! exit-flag, 0=convergence
integer :: fef,iter,lsiter,linf
logical :: lsconv
real(8) :: ff,fft,alpha,dirder,alphanew
real(8),dimension(n) :: xc,grad,dir,xct
real(8),dimension(n,1) :: dirm
real(8),dimension(n,n) :: hess,hesscp

xc = x0;
eflag = 16;
do iter=1,MAXITER ! newton iteration master loop
!  evaluate objective function with gradient and hessian
   call ddlfun(n,xc,ff,grad,hess,fef)
   hesscp = hess;

! check for convergence
   if(maxval(abs(grad))<1.0e-12_8) then

      eflag = 0;
      exit;
   endif
   
   dirm(:,1) = -grad;
   ! cholesky 
   call dpotrf('L',n,hess,n,linf)

   ! if not hess>0 do hessian modification
   if(abs(linf)>0) then
      hess = hesscp;
      call mod_and_inv(n,hess,linf);
      if(abs(linf)>0) then
         eflag = 7;
         exit;
      endif
      dir = -matmul(hess,grad);
   else
   ! else, backsubsitution based on cholesky factorization
      call dpotrs('L',n,1,hess,n,dirm,n,linf)
      dir = dirm(:,1);
   endif


  ! check the search direction
   dirder = sum(dir*grad);
   if(dirder>0.0_8) then
      eflag = 4;
      exit;
   endif
   ! Start line search
   alpha = 1.0_8;
   lsconv = .false.
   do lsiter = 1,20
      ! try this point
      xct = xc + alpha*dir;
      call lfun(n,xct,fft,fef)
      ! check first Wolfe condition, if OK exit
      if(fft<ff + c1*alpha*dirder + 100.0_8*epsilon(ff)) then
         xc = xct;
         lsconv = .true.
         exit;
      endif
      ! otherwise 
      ! do quadratic interpolation to find new alpha
      alphanew = 0.5_8*dirder*alpha**2/(ff + dirder*alpha - fft);
      alpha = min(0.9_8*alpha,max(0.1_8*alpha,alphanew))
   enddo
   ! end line search
   ! check if linesearch succeded, if not-exit.
   if(.not. lsconv) then
      eflag = 5;
      exit;
   endif
enddo ! end newton iteration loop

! copy local variables to output
xout = xc;
fout = ff;
gout = grad;
hout = hesscp;

end subroutine n_solver

end module newton_solver

! these are wrapper-functions used to pass to the optimizer
subroutine lfun(p,pp,ll,eflag)
implicit none
integer,intent(in) :: p
real(8),dimension(p),intent(in) :: pp
real(8),intent(out) :: ll
integer,intent(out) :: eflag
call lwrapper(p,pp,ll,eflag)
end subroutine lfun

subroutine dlfun(p,pp,dll,eflag)
implicit none
integer,intent(in) :: p
real(8),dimension(p),intent(in) :: pp
real(8),dimension(p),intent(out) :: dll
integer,intent(out) :: eflag
call dlwrapper(p,pp,dll,eflag)
end subroutine dlfun

subroutine ddlfun(p,pp,ll,dll,ddll,eflag)
implicit none
integer,intent(in) :: p
real(8),dimension(p),intent(in) :: pp
real(8),intent(out) :: ll
real(8),dimension(p),intent(out) :: dll
real(8),dimension(p,p),intent(out) :: ddll
integer,intent(out) :: eflag
call ddlwrapper(p,pp,ll,dll,ddll,eflag)
end subroutine ddlfun


! these functions wrap the original objective function and
! output of the automatic differentiation tool.
! Notice that the global variables are passed to the objective
! function here.
subroutine lwrapper(p,pp,ll,eflag)
use globals
implicit none
integer,intent(in) :: p
real(8),dimension(p),intent(in) :: pp
real(8),intent(out) :: ll
integer,intent(out) :: eflag
call lgobfun(gn,gx,gy,gwts,gx0,gy0,pp,ghx,ghy,ll,cv,fixrho)
ll = - ll;
eflag = 0;
end subroutine lwrapper

subroutine dlwrapper(p,pp,dll,eflag)
use globals
implicit none
integer,intent(in) :: p
real(8),dimension(p),intent(in) :: pp
real(8),dimension(p),intent(out) :: dll
integer,intent(out) :: eflag
real(8),dimension(5,5) :: ppd
real(8) :: ll
integer :: i
ppd = 0.0_8;
do i=1,5
   ppd(i,i) = 1.0_8;
enddo
! this is the gradient function
call lgobfun_dv(gn,gx,gy,gwts,gx0,gy0,pp,ppd,ghx,ghy,ll,dll,cv,fixrho,5)
dll = -dll;
eflag = 0;
end subroutine dlwrapper

subroutine ddlwrapper(p,pp,ll,dll,ddll,eflag)
use globals
implicit none
integer,intent(in) :: p
real(8),dimension(p),intent(in) :: pp
real(8),intent(out) :: ll
real(8),dimension(p),intent(out) :: dll
real(8),dimension(p,p),intent(out) :: ddll
integer,intent(out) :: eflag
real(8),dimension(5,5) :: ppd
real(8),dimension(5,5) :: ppd0
integer :: i

ppd = 0.0_8;
ppd0 = 0.0_8;
do i=1,5
   ppd(i,i) = 1.0_8;
   ppd0(i,i) = 1.0_8;
enddo
! this is the Hessian function
call lgobfun_dv_dv(gn,gx,gy,gwts,gx0,gy0,pp,ppd0,ppd,&
     ghx,ghy,ll,dll,ddll,cv,fixrho,5,5)

ll = -ll;
dll = -dll;
ddll = -ddll;
eflag = 0;
end subroutine ddlwrapper

!-----------------------------------------------------------------------
! subroutine localgauss
!------------------------------------------------------------------------
! This is the routine called by R. Essentially it loops over all the
! optimization problem required
!------------------------------------------------------------------------
subroutine localgauss(n,x,y,nrep,x0,y0,hx,hy,rhof,results,hess)
use globals
!use mle_opt
use newton_solver
implicit none
real(8),parameter :: twopi = 6.283185307179586e+00_8
! same interface used to pass to the optimizer
interface 
subroutine lfun(p,pp,ll,eflag)
implicit none
integer,intent(in) :: p
real(8),dimension(p),intent(in) :: pp
real(8),intent(out) :: ll
integer,intent(out) :: eflag
end subroutine lfun
end interface
interface
subroutine dlfun(p,pp,dll,eflag)
implicit none
integer,intent(in) :: p
real(8),dimension(p),intent(in) :: pp
real(8),dimension(p),intent(out) :: dll
integer,intent(out) :: eflag
end subroutine dlfun
end interface
interface
subroutine ddlfun(p,pp,ll,dll,ddll,eflag)
implicit none
integer,intent(in) :: p
real(8),dimension(p),intent(in) :: pp
real(8),intent(out) :: ll
real(8),dimension(p),intent(out) :: dll
real(8),dimension(p,p),intent(out) :: ddll
integer,intent(out) :: eflag
end subroutine ddlfun
end interface
! end interfaces
integer,intent(in) :: n ! number of data-points
real(8),dimension(n),intent(in) :: x ! x-data
real(8),dimension(n),intent(in) :: y ! y-data
integer,intent(in) :: nrep ! number of estimations
real(8),dimension(nrep),intent(in) :: x0 ! x-coordinate at estimation point
real(8),dimension(nrep),intent(in) :: y0 ! y-coordinate at estimation point
real(8),dimension(nrep),intent(in) :: hx ! bandwidth in x-dir
real(8),dimension(nrep),intent(in) :: hy ! bandwidth in y-dir
real(8),dimension(nrep),intent(in) :: rhof ! vector indicating if rho should be estimated
real(8),dimension(nrep,7),intent(out) :: results ! matrix with results
real(8),dimension(nrep,5,5),intent(out) :: hess ! tensor of hessian matrices
real(8),dimension(5) :: pp,ppout,gradout
real(8),dimension(5,5) :: hessout
real(8) :: meanx,meany,stdx,stdy
real(8) :: fout
integer :: eflag
integer :: rep

! allocate space for the global variables
allocate(gx(n))
allocate(gy(n))
allocate(gwts(n))
! copy data to the global variables
gn = n
gx = x
gy = y
pp = 0.0_8;

! calculations for initial conditions
meanx = sum(x)/(1.0_8*n);
meany = sum(y)/(1.0_8*n);
stdx = sqrt(sum((x-meanx)**2)/(1.0_8*(n-1)));
stdy = sqrt(sum((y-meany)**2)/(1.0_8*(n-1)));

! this is the initial vector passed to the optimizer
! Note, optimization is performed in log-standard-deviations and
! a logit-transformation for rho
pp(1) = meanx;
pp(2) = meany;
pp(3) = log(2.0_8*stdx);
pp(4) = log(2.0_8*stdy);
pp(5) = 0.0_8;

! loop over estimation replica
do rep = 1,nrep
  ! copy data to the global variables
   cv = .true.
   fixrho = rhof(rep);
   gx0 = x0(rep)
   gy0 = y0(rep)
   ghx = hx(rep)
   ghy = hy(rep)
   ! compute the weights
   gwts = exp(-0.5_8*(((x-gx0)/ghx)**2 + ((y-gy0)/ghy)**2))/(twopi*ghx*ghy)

   ! check that at least some weights are non-zero
   if(maxval(gwts)<1.0e-13_8) then
      results(rep,1:6) = -999.0_8;
      results(rep,7) = 22;
   else
      ! call to the optimizer
      call n_solver(lfun,dlfun,ddlfun,5,pp,ppout,fout,gradout,hessout,&
           eflag)
      ! copy output to the results-matrix
      ! means
      results(rep,1:2) = ppout(1:2);
      ! standard-deviations
      results(rep,3:4) = exp(ppout(3:4));
      ! correlations (if fixed, the fixed rho is passed)
      if(abs(fixrho)>=1.0_8) then
         results(rep,5) = -1.0_8 + 2.0_8*exp(ppout(5))/(1.0_8+exp(ppout(5)));
      else
         results(rep,5) = fixrho
      endif
      ! value of objective function
      results(rep,6) = fout;
      ! exif flag from optimizer
      results(rep,7) = eflag
      cv = .false.
      ! compute hessian in original parameters
      call ddlfun(5,results(rep,1:5),fout,gradout,hessout,eflag)
      ! store hessian in hessian-tensor
      hess(rep,:,:) = -hessout;
   endif
enddo
! deallocate dynamic memory for global vectors
deallocate(gwts)
deallocate(gy)
deallocate(gx)
end subroutine localgauss

