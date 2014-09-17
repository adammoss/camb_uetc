subroutine fevolve(SPR,n,x,y,yprime)
  use ModelParams
  use LambdaGeneral
  use uetc
  implicit none
  Type(StringParams) SPR
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
   ct=(SPRa%cr+a*SPRa%g*SPRa%cm)/(1.0d0+a*SPRa%g)
   fk=(SPRa%fkr+a*SPRa%g*SPRa%fkm)/(1.0d0+a*SPRa%g)
   xlprime=adotoa*xl*xv2+0.5d0*ct*xv
   xvprime=(1.0d0-xv2)*(fk/xl-2.0d0*adotoa*xv)  
   yprime(2)=xlprime
   yprime(3)=xvprime
end subroutine fevolve
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
subroutine fevolve4(SPR,n,x,y,yprime)
  use ModelParams
  use LambdaGeneral
  use uetc
  implicit none
  type(StringParams) SPR
  integer n
  real(DP), parameter :: pi1 = 3.1415926535897932384626433832795d0
  real(DP) :: x,y(n),yprime(n)
  real(DP) :: tau,a,a2,grho,adotoa
  real(DP) :: xl,xv,xv2,fk,xlprime,xvprime

      tau=x
      a=y(1)     
      a2=a*a

      grho=grhob/a+grhoc/a+grhog/a2+grhornomass/a2+grhov*a2

      adotoa=dsqrt(grho/3.0d0)
      yprime(1)=adotoa*a
      xl=y(2)
      xv=y(3)
      xv2=xv*xv


      fk=(2.0d0*dsqrt(2.0d0)/pi1)*&
         (1.0d0-8.0d0*xv**6)/(1.0d0+8.0d0*xv**6)

      xlprime=adotoa*xl*xv2+0.5d0*SPRa%cr*xv
      xvprime=(1.0d0-xv2)*(fk/xl-2.0d0*adotoa*xv)

      yprime(2)=xlprime
      yprime(3)=xvprime

end subroutine fevolve4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function adtauda(a)
  use uetc
  implicit none
  real(DP) adtauda,dtauda,a
  external dtauda
  adtauda = dtauda(a)*a
end function adtauda
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

