subroutine fevolve(SPR,n,x,y,yprime)
  use ModelParams
  use LambdaGeneral
  use uetc
  implicit none
  type(StringParams) SPR
  integer n
  real(DP) :: x,y(n),yprime(n)
  real(DP) :: a,a2,tau
  real(DP) :: grhob_t,grhor_t,grhoc_t,grhog_t,grhov_t,grho
  real(DP) :: ct,fk
  real(DP) :: adotoa,xl,xv,xv2,xlprime,xvprime
 
  tau = x
  a=y(1)     
  a2=a*a
  grhob_t=grhob/a
  grhoc_t=grhoc/a
  grhor_t=grhornomass/a2
  grhog_t=grhog/a2
  if (w_lam==-1._dl) then
     grhov_t=grhov*a2
  else
     grhov_t=grhov*a**(-1-3*w_lam)
  end if
  grho=grhob_t+grhoc_t+grhor_t+grhog_t+grhov_t
  if (CP%flat) then
      adotoa=sqrt(grho/3.0)
   else
      stop 'Non-flat models not supported'
   end if
   yprime(1)=adotoa*a
   xl=y(2)
   xv=y(3)
   xv2=xv*xv
   ct=(SPR%cr+a*SPR%g*SPR%cm)/(1.0d0+a*SPR%g)
   fk=(SPR%fkr+a*SPR%g*SPR%fkm)/(1.0d0+a*SPR%g)
   xlprime=adotoa*xl*xv2+0.5d0*ct*xv
   xvprime=(1.0d0-xv2)*(fk/xl-2.0d0*adotoa*xv)  
   yprime(2)=xlprime
   yprime(3)=xvprime
end subroutine fevolve
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
function adtauda(a)
  use uetc
  implicit none
  real(DP) adtauda,dtauda,a
  external dtauda
  adtauda = dtauda(a)*a
end function adtauda
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

