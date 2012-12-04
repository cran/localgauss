!----------------------------------------------------------------------
! This file defines the objective function for local gaussian correlation
!----------------------------------------------------------------------
! Throughout, we use the convetion for parameters of a bivariate gaussian
! pars(1) = mu_x
! pars(2) = mu_y
! pars(3) = sigma_x
! pars(4) = sigma_y
! pars(5) = rho
!-------------------------------------------------------------------------

! log-density of a bivariate Gaussian, vectorized in the data
! with parameters fixed
subroutine loggausspdf(n,x,y,pars,res)
implicit none
real(8),parameter :: twopi = 6.283185307179586e+00_8
integer,intent(in) :: n
real(8),dimension(n),intent(in) :: x
real(8),dimension(n),intent(in) :: y
real(8),dimension(5),intent(in) :: pars
real(8),dimension(n),intent(out) :: res
real(8),dimension(n) :: cen1,cen2
real(8) :: t1,f1,f2,f12
t1 = -0.5_8/(1.0_8-pars(5)**2);
f1 = t1/(pars(3)**2);
f2 = t1/(pars(4)**2);
f12 = -2.0*pars(5)*t1/(pars(3)*pars(4));
cen1 = x-pars(1);
cen2 = y-pars(2);
res = -log(twopi*pars(3)*pars(4)*sqrt(1.0_8-pars(5)**2)) &
     + f1*cen1**2 + f2*cen2**2 + f12*cen1*cen2;
end subroutine


! This is the objective function
subroutine lgobfun(n,x,y,wts,x0,y0,pp,hx,hy,ll,cv,fixrho)
implicit none
integer,intent(in) :: n ! number of data points
real(8),dimension(n),intent(in) :: x ! x-data points
real(8),dimension(n),intent(in) :: y ! y-data points
real(8),dimension(n),intent(in) :: wts ! kernel weights of each data point
real(8),intent(in) :: x0 ! x-coordinate where local parameters are computed
real(8),intent(in) :: y0 ! y-coordinate where local parameters are computed
real(8),dimension(5),intent(in) :: pp ! transformed parameters used by optimizer
real(8),intent(in) :: hx ! x-dir bandwidth
real(8),intent(in) :: hy ! y-dir bandwidth
logical,intent(in) :: cv ! should the input parameters be transformed?
real(8),intent(in) :: fixrho ! fixed rho?
real(8),intent(out) :: ll ! objective function out
real(8),dimension(n) :: lgauss
real(8),dimension(5) :: pars2
real(8),dimension(1) :: xtmp,ytmp,restmp
real(8),dimension(5) :: pars
ll = 0.0_8;
! if cv=TRUE; Transform parameters
if(cv) then
   pars(1:2) = pp(1:2);
   pars(3:4) = exp(pp(3:4));
   if(abs(fixrho)<1.0_8) then
      pars(5) = fixrho;
      ! add this to make optimizer work
      ! if rho is fixed
      ll = -0.5_8*pp(5)**2
   else
      pars(5) = -1.0_8 + 2.0_8*exp(pp(5))/(1.0_8+exp(pp(5)));
   endif
! otherwise, pass parameters as they are
else
   pars = pp;
   if(abs(fixrho)<1.0_8) then
      pars(5) = fixrho;
   endif
endif
! compute logdensites for the weighted log-likelihood
call loggausspdf(n,x,y,pars,lgauss)
! weighted log-likelihood
ll = ll + sum(wts*lgauss)/(1.0_8*n);

! work out parameters for the penalty term
pars2(1:2) = pars(1:2);
pars2(3) = sqrt(pars(3)**2 + hx**2);
pars2(4) = sqrt(pars(4)**2 + hy**2);
pars2(5) = pars(5)*pars(3)*pars(4)/(pars2(3)*pars2(4))
xtmp(1) = x0;
ytmp(1) = y0;
call loggausspdf(1,xtmp,ytmp,pars2,restmp);
! add penalty term to objective
ll = ll - exp(restmp(1));
end subroutine lgobfun


