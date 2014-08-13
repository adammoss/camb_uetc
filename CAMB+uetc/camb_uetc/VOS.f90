  subroutine VOS(nint,ktau_array,xi0,v0,alpha0)
    use ModelParams
    implicit none
    Type(StringParams) :: SPR
    integer, intent(in) :: nint
    real(DP), intent(in) :: ktau_array(nint),xi0,v0,alpha0
    real(DP) :: t0,a0,tau,tauend,tol1,tol2
    real(DP) :: rombint, dtau
    integer, parameter :: l = 3
    real(DP) :: yev(l),yevpr(l),c(24),w(l,9),adot(nint)
    integer :: i,n,m
    real(DP), allocatable :: v(:),xi(:),alpha(:),tphys(:),a(:),t(:)
    real(DP), allocatable :: vpr(:),xipr(:),alphapr(:),tphyspr(:),apr(:),tpr(:)
    external fevolve, rombint, dtauda, adtauda

    n = nint

    allocate(v(n),xi(n),alpha(n),tphys(n),a(n),t(n))
    allocate(vpr(n),xipr(n),alphapr(n),tphyspr(n),apr(n),tpr(n) )

    !Initial time from which the VOS is run
    t0 = 0.1d0*ktau_array(1)
    !Initial scale factor
    a0 = adotrad*t0
    !Evolution parameter for the scale factor
    yev(1) = a0
    !Evolution parameter for the correlation length
    yev(2) = xi0
    !Evolution parameter for the velocity
    yev(3) = v0
    tau = t0
    m = 1
    tol1 = 1.0d-6
    tol2 = 1.0d-4
    !Evolve the sting parameters and their derivatives
    do i = 1,n
    tauend = ktau_array(i)
    call dverk(SPR,l,fevolve,tau,yev,tauend,tol1,m,c,l,w)                                                                            
    call fevolve(SPR,l,tauend,yev,yevpr)
    !Derivative of the scale factor and string parameters
    adot(i) =  yevpr(1)
    a(i) =  yev(1)
    xi(i) = yev(2)
    alpha(i) = 1.0d0+(alpha0-1.0d0)/(adot(i)*tau/yev(1))
    v(i) = yev(3)
    !t as a function of tau
    tphys(i) = rombint(adtauda,0.0d0,yev(1),tol2)
    end do
    !Spline the one-scale functions of time and physical time
    call spline(tphys,t,n,d0lo,d0hi,tpr)
    call spline(t,tphys,n,d0lo,d0hi,tphyspr)
    call spline(t,xi,n,d0lo,d0hi,xipr)
    call spline(t,alpha,n,d0lo,d0hi,alphapr)
    call spline(t,v,n,d0lo,d0hi,vpr)
    call spline(t,a,n,d0lo,d0hi,apr)
  end subroutine

