! STRING MODULE FOR USE IN CAMB/COSMOMC - AM (Nottingham) 
! Uses DVERK integrator
! May 2012
! Computes:
! get_Vk computes V*k, where V is the vector metric potential as in eqn (88)
! of Hu and White, PRD56 (1997) 596
    
module string
  use ModelParams
  use Random
  use am_routines
  implicit none
  integer, parameter :: DP = kind(1.0D0)
  logical :: do_strings = .false.
  logical :: do_string_source = .true.
  real(DP) :: tau_init,taustart_string
  real(DP) :: alfrad,vdev,dksi,xlf,gmu,tkmax      
  integer :: irandomv,ns1,iseed0,iseed1,iseed2
  ! Evolve 
  integer, parameter :: nstep00 = 200
  real(DP) :: evv(nstep00),evvpr(nstep00)
  real(DP) :: cole(nstep00),colepr(nstep00)
  real(DP) :: evalf(nstep00),evalfpr(nstep00)
  real(DP) :: cnorm(nstep00),cnormpr(nstep00)
  real(DP) :: tphys(nstep00),tphyspr(nstep00)
  real(DP) :: tim(nstep00),timpr(nstep00)
  real(DP) :: xaa(nstep00),xaapr(nstep00)
  ! Mixup
  integer, parameter :: ns0 = 801
  real(DP) :: xp3(ns0,nstep00),xp3pr(ns0,nstep00)
  real(DP) :: xd3(ns0,nstep00),xd3pr(ns0,nstep00)
  real(DP) :: ps(ns0,nstep00),pspr(ns0,nstep00)
  real(DP) :: pv(ns0,nstep00),pvpr(ns0,nstep00)
  real(DP) :: pt(ns0,nstep00),ptpr(ns0,nstep00)   
  real(DP) :: x0k(ns0),tf(ns0),xlden(ns0)
  ! Spline
  integer, parameter :: nk_max = 500
  integer, parameter :: nt_fine=10000
  integer :: nt_f(nk_max)
  real(DP) :: ts_f(nk_max,nt_fine)
  real(DP) :: em00(nk_max,nt_fine),em00pr(nk_max,nt_fine)
  real(DP) :: emS(nk_max,nt_fine),emSpr(nk_max,nt_fine)
  real(DP) :: emD(nk_max,nt_fine),emDpr(nk_max,nt_fine)
  real(DP) :: emP(nk_max,nt_fine),emPpr(nk_max,nt_fine)          
  real(DP) :: emV(nk_max,nt_fine),emVpr(nk_max,nt_fine)
  real(DP) :: emT(nk_max,nt_fine),emTpr(nk_max,nt_fine)
  real(DP) :: as(nk_max,nt_fine),aspr(nk_max,nt_fine)
  real(DP) :: Vk(nk_max,nt_fine),Vkpr(nk_max,nt_fine)
  real(DP), parameter :: d0lo=1.0d40
  real(DP), parameter :: d0hi=1.0d40 
  logical :: debug = .false.
  logical, parameter :: random_seed = .false.
  logical, parameter :: random_nr = .true.
  logical :: evolve 
  character(len=100) :: str_output
  logical, parameter :: single_string = .false.
  real(DP) :: cr,cm,g,fkr,fkm
  
  type StrEvolutionVars
     integer wk_ix
     real(DP) wk,emtD,emtP
  end type StrEvolutionVars
 
contains
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  subroutine init_string(lmaxout)
    implicit none
    integer, intent(in) :: lmaxout
    ! Initialize string model parameters	
    logical :: do_string_init = .true. 
    ! the string tension (the most important parameter)
    gmu=1.1d-6   
    ! Whether to evolve string parameters with time
    evolve = .true.
    ! Wiggliness in the radiation era                  
    alfrad=1.9d0
    ! Do segments randomly change velocity around each Hubble time? yes=1, no=0
    irandomv=0
    !  cr,cm,kr,km must be set.
    !cr=0.23d0
    !cm=0.18d0
    !g=300.0d0
    !fkr=0.17d0
    !fkm=0.49d0
    cr=0.23d0
    cm=0.18d0
    g=300.0d0
    fkr=0.17d0
    fkm=0.49d0
    ! vdev is the initial rms string velocity (over all scales)
    !vdev=0.65d0
    vdev = sqrt(fkr/(cr+fkr))
    ! tkmax is the cutoff above which \Theta_D and \Theta_P are set to zero
    tkmax=dble(lmaxout)/2.0d0
    ! dksi is the initial correlation length/initial conformal time
    !dksi=0.13d0      
    dksi = fkr/(2.0*vdev)
    ! ns1 is the number of consolidated string segments
    ns1=500
    if(single_string) ns1 = 1
    ! xlf controls the rate at which strings decay
    !xlf=0.5d0
    xlf = 0.99
    !  tau_init is the earliest conformal time at which one scale model is run   
    tau_init=2.d-2
    ! taustart_string
    taustart_string = 2.2d0
    write(*,*) '************* String model parameters *************'
    !  initialize the random # generator    
    if (random_seed) then
       call InitRandom
    else
       iseed0=-20
       iseed1=25
       iseed2=10
       if (random_nr) then
          write(*,*) 'NR seed:',iseed0
       else
          call InitRandom(iseed1,iseed1)
       end if
    end if
    write(*,*) 'Single string:',single_string
    write(*,*) 'One scale evolution:',evolve
    write(*,*) 'Initial RMS velocity:',vdev
    write(*,*) 'Initial correlation length/tau:',dksi
    write(*,*) 'Initial wiggliness:',alfrad
    write(*,*) 'T_k max:',tkmax
    write(*,*) 'Number of consolidated segments:',ns1
    write(*,*) 'Rate of string decay:',xlf
    write(*,*) '***************************************************'
    do_string_init = .false.
  end subroutine init_string
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  subroutine strings_k(nexp,lmaxout)
    use am_routines
    implicit none
    integer, intent(in) :: nexp,lmaxout
    integer iexp,ik,nk
    integer, parameter :: nk_max_arr = 100
    real(DP) :: tau_test,k_arr(nk_max_arr)
    real(DP) :: emt00,emtS,emtV,emtT,emtD,emtP,emt00dot,emtSdot,V
    real(DP) :: emt00_mean(nk_max_arr),emt00_var(nk_max_arr)
    real(DP) :: emtS_mean(nk_max_arr),emtS_var(nk_max_arr)
    real(DP) :: emtV_mean(nk_max_arr),emtV_var(nk_max_arr)
    real(DP) :: emtT_mean(nk_max_arr),emtT_var(nk_max_arr)
    real(DP) :: emtD_mean(nk_max_arr),emtD_var(nk_max_arr)
    real(DP) :: emtP_mean(nk_max_arr),emtP_var(nk_max_arr)
    real(DP) :: emt00dot_mean(nk_max_arr),emt00dot_var(nk_max_arr)
    real(DP) :: emtSdot_mean(nk_max_arr),emtSdot_var(nk_max_arr)
    real(DP) :: Vk_mean(nk_max_arr),Vk_var(nk_max_arr)
    character(len=100) :: filename
    call init_string(lmaxout)
    call prepare
    !tau_test = 290.0
    tau_test = 100.0
    nk = 25
    write(*,*) 'Doing string k with tau:',tau_test
    filename = trim(str_output)//'emt00.dat'
    open(unit=60,file=filename,status='unknown')
    filename = trim(str_output)//'emtS.dat'
    open(unit=61,file=filename,status='unknown')
    filename = trim(str_output)//'emtV.dat'
    open(unit=62,file=filename,status='unknown')
    filename = trim(str_output)//'emtT.dat'
    open(unit=63,file=filename,status='unknown')
    filename = trim(str_output)//'emtD.dat'
    open(unit=64,file=filename,status='unknown')
    filename = trim(str_output)//'emtP.dat'
    open(unit=65,file=filename,status='unknown')
    filename = trim(str_output)//'emt00dot.dat'
    open(unit=66,file=filename,status='unknown')
    filename = trim(str_output)//'emtSdot.dat'
    open(unit=67,file=filename,status='unknown')  
    filename = trim(str_output)//'Vk.dat'
    open(unit=68,file=filename,status='unknown')
    emt00_mean(:) = 0.0
    emt00_var(:) = 0.0
    emtS_mean(:) = 0.0
    emtS_var(:) = 0.0
    emtV_mean(:) = 0.0
    emtV_var(:) = 0.0
    emtT_mean(:) = 0.0
    emtT_var(:) = 0.0
    emtD_mean(:) = 0.0
    emtD_var(:) = 0.0
    emtP_mean(:) = 0.0
    emtP_var(:) = 0.0
    emt00dot_mean(:) = 0.0
    emt00dot_var(:) = 0.0
    emtSdot_mean(:) = 0.0
    emtSdot_var(:) = 0.0
    Vk_mean(:) = 0.0
    Vk_var(:) = 0.0
    !k_arr(1) = 6.852d-6
    write(*,*) ' K modes:'
    do ik=1,nk 
       k_arr(ik) = 10.0**(-3.0+3.0*(ik-1)/(nk-1))
       write(*,*) ik,k_arr(ik)
    end do
    do iexp=1,nexp
       write(*,*) 'N exp:',iexp
       call mixup
       call strings_spline_arr(nk,k_arr,0.218d0,16504.8d0)
       do ik=1,nk
          emt00 = get_em00(ik,tau_test)
          emtS = get_emS(ik,tau_test)
          emtV = get_emV(ik,tau_test)
          emtT = get_emT(ik,tau_test)
          emtD = get_emD(ik,tau_test)
          emtP = get_emP(ik,tau_test)    
          emt00dot = get_em00_dot(ik,tau_test)
          emtSdot = get_emS_dot(ik,tau_test)
          V = get_Vk(ik,tau_test)
          emt00_mean(ik) = emt00_mean(ik) + emt00
          emt00_var(ik) = emt00_var(ik) + emt00**2.0
          emtS_mean(ik) = emtS_mean(ik) + emtS
          emtS_var(ik) = emtS_var(ik) + emtS**2.0
          emtV_mean(ik) = emtV_mean(ik) + emtV
          emtV_var(ik) = emtV_var(ik) + emtV**2.0
          emtT_mean(ik) = emtT_mean(ik) + emtT
          emtT_var(ik) = emtT_var(ik) + emtT**2.0
          emtD_mean(ik) = emtD_mean(ik) + emtD
          emtD_var(ik) = emtD_var(ik) + emtD**2.0
          emtP_mean(ik) = emtP_mean(ik) + emtP
          emtP_var(ik) = emtP_var(ik) + emtP**2.0
          emt00dot_mean(ik) = emt00dot_mean(ik) + emt00dot
          emt00dot_var(ik) = emt00dot_var(ik) + emt00dot**2.0
          emtSdot_mean(ik) = emtSdot_mean(ik) + emtSdot
          emtSdot_var(ik) = emtSdot_var(ik) + emtSdot**2.0
          Vk_mean(ik) = Vk_mean(ik) + V
          Vk_var(ik) = Vk_var(ik) + V**2.0
       end do
    end do
    do ik=1,nk 
       emt00_mean(ik) =  emt00_mean(ik)/nexp
       emt00_var(ik) =  emt00_var(ik)/nexp
       emtS_mean(ik) =  emtS_mean(ik)/nexp
       emtS_var(ik) =  emtS_var(ik)/nexp
       emtV_mean(ik) =  emtV_mean(ik)/nexp
       emtV_var(ik) =  emtV_var(ik)/nexp
       emtT_mean(ik) =  emtT_mean(ik)/nexp
       emtT_var(ik) =  emtT_var(ik)/nexp
       emtD_mean(ik) =  emtD_mean(ik)/nexp
       emtD_var(ik) =  emtD_var(ik)/nexp
       emtP_mean(ik) =  emtP_mean(ik)/nexp
       emtP_var(ik) =  emtP_var(ik)/nexp 
       emt00dot_mean(ik) =  emt00dot_mean(ik)/nexp
       emt00dot_var(ik) =  emt00dot_var(ik)/nexp
       emtSdot_mean(ik) =  emtSdot_mean(ik)/nexp
       emtSdot_var(ik) =  emtSdot_var(ik)/nexp
       Vk_mean(ik) =  Vk_mean(ik)/nexp
       Vk_var(ik) =  Vk_var(ik)/nexp
       write(60,'(10E15.5)') k_arr(ik),emt00_mean(ik),emt00_var(ik)
       write(61,'(10E15.5)') k_arr(ik),emtS_mean(ik),emtS_var(ik)
       write(62,'(10E15.5)') k_arr(ik),emtV_mean(ik),emtV_var(ik)
       write(63,'(10E15.5)') k_arr(ik),emtT_mean(ik),emtT_var(ik)
       write(64,'(10E15.5)') k_arr(ik),emtD_mean(ik),emtD_var(ik)
       write(65,'(10E15.5)') k_arr(ik),emtP_mean(ik),emtP_var(ik)
       write(66,'(10E15.5)') k_arr(ik),emt00dot_mean(ik),emt00dot_var(ik)
       write(67,'(10E15.5)') k_arr(ik),emtSdot_mean(ik),emtSdot_var(ik)
       write(68,'(10E15.5)') k_arr(ik),Vk_mean(ik),Vk_var(ik)
    end do
    close(60)
    close(61)
    close(62)
    close(63)
    close(64)
    close(65)
    close(66)
    close(67)
    close(68)
  end subroutine strings_k
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  subroutine strings_uetc(nexp,lmaxout)
    use am_routines
    implicit none
    integer, intent(in) :: nexp,lmaxout
    integer iexp,ik,it,jt,nk,nt
    integer, parameter :: nt_max_arr = 400
    integer, parameter :: nk_max_arr = 1
    integer, parameter :: nexp_max_arr = 20000
    real(DP) :: tau_arr(nt_max_arr),k_arr(nk_max_arr)
    real(DP) :: emt00(nexp_max_arr,nk_max_arr,nt_max_arr),&
         emt00_uetc(nk_max_arr,nt_max_arr,nt_max_arr)
    real(DP) :: emtS(nexp_max_arr,nk_max_arr,nt_max_arr),&
         emtS_uetc(nk_max_arr,nt_max_arr,nt_max_arr)
    real(DP) :: emtV(nexp_max_arr,nk_max_arr,nt_max_arr),&
         emtV_uetc(nk_max_arr,nt_max_arr,nt_max_arr)
    real(DP) :: emtT(nexp_max_arr,nk_max_arr,nt_max_arr),&
         emtT_uetc(nk_max_arr,nt_max_arr,nt_max_arr)
    real(DP) :: emtD(nexp_max_arr,nk_max_arr,nt_max_arr),&
         emtD_uetc(nk_max_arr,nt_max_arr,nt_max_arr)
    real(DP) :: emtP(nexp_max_arr,nk_max_arr,nt_max_arr),&
         emtP_uetc(nk_max_arr,nt_max_arr,nt_max_arr)
    real(DP) :: emt00dot(nexp_max_arr,nk_max_arr,nt_max_arr),&
         emt00dot_uetc(nk_max_arr,nt_max_arr,nt_max_arr)
    real(DP) :: emtSdot(nexp_max_arr,nk_max_arr,nt_max_arr),&
         emtSdot_uetc(nk_max_arr,nt_max_arr,nt_max_arr)
    real(DP) emt00S_uetc(nk_max_arr,nt_max_arr,nt_max_arr)
    character(len=100) :: filename,fmt 
    call init_string(lmaxout)
    call prepare
    write(*,*) 'Doing UETC for k:'
    nk = 1
    nt = 400
    filename = trim(str_output)//'k_uetc.dat'
    open(unit=60,file=filename,status='unknown')
    do ik=1,nk 
       k_arr(1) = 0.1
       !k_arr(ik) = 10.0**(-5.0+5.0*(ik-1)/(nk-1))
       write(*,*) ik,k_arr(ik)
       write(60,*) ik,k_arr(ik)
    end do
    close(60)
    filename = trim(str_output)//'tau_uetc.dat'
    open(unit=60,file=filename,status='unknown')
    do it=1,nt 
       tau_arr(it) = 10.0**(-1.0+5.0*(it-1)/(nt-1))
       write(60,*) it,tau_arr(it)
    end do
    close(60)
    emt00(:,:,:) = 0.0
    emtS(:,:,:) = 0.0
    emtV(:,:,:) = 0.0
    emtT(:,:,:) = 0.0
    emtD(:,:,:) = 0.0
    emtP(:,:,:) = 0.0
    emt00dot(:,:,:) = 0.0
    emtSdot(:,:,:) = 0.0
    do iexp=1,nexp
       write(*,*) 'N exp:',iexp
       call mixup
       call strings_spline_arr(nk,k_arr,0.218d0,16504.8d0)
       do ik=1,nk
          do it=1,nt
             emt00(iexp,ik,it) = get_em00(ik,tau_arr(it))
             emtS(iexp,ik,it) = get_emS(ik,tau_arr(it))
             emtV(iexp,ik,it) = get_emV(ik,tau_arr(it))
             emtT(iexp,ik,it) = get_emT(ik,tau_arr(it))
             emtD(iexp,ik,it) = get_emD(ik,tau_arr(it))
             emtP(iexp,ik,it) = get_emP(ik,tau_arr(it)) 
             emt00dot(iexp,ik,it)  = get_em00_dot(ik,tau_arr(it))
             emtSdot(iexp,ik,it)  = get_emS_dot(ik,tau_arr(it))
             !write(90,*) tau_arr(it), emtV(iexp,ik,it)
          end do
       end do
       !stop
    end do 
    emt00_uetc(:,:,:) = 0.0d0
    emtS_uetc(:,:,:) = 0.0d0
    emtV_uetc(:,:,:) = 0.0d0
    emtT_uetc(:,:,:) = 0.0d0
    emtD_uetc(:,:,:) = 0.0d0
    emtP_uetc(:,:,:) = 0.0d0
    emt00dot_uetc(:,:,:) = 0.0d0
    emtSdot_uetc(:,:,:) = 0.0d0
    emt00S_uetc(:,:,:) = 0.0d0
    do iexp=1,nexp 
       do ik=1,nk
          do it=1,nt
             do jt=1,nt
                emt00_uetc(ik,it,jt) = emt00_uetc(ik,it,jt) + &
                     emt00(iexp,ik,it)*emt00(iexp,ik,jt)
                emtS_uetc(ik,it,jt) = emtS_uetc(ik,it,jt) + &
                     emtS(iexp,ik,it)*emtS(iexp,ik,jt)
                emtV_uetc(ik,it,jt) = emtV_uetc(ik,it,jt) + &
                     emtV(iexp,ik,it)*emtV(iexp,ik,jt)
                emtT_uetc(ik,it,jt) = emtT_uetc(ik,it,jt) + &
                     emtT(iexp,ik,it)*emtT(iexp,ik,jt)
                emtD_uetc(ik,it,jt) = emtD_uetc(ik,it,jt) + &
                     emtD(iexp,ik,it)*emtD(iexp,ik,jt)
                emtP_uetc(ik,it,jt) = emtP_uetc(ik,it,jt) + &
                     emtP(iexp,ik,it)*emtP(iexp,ik,jt)
                emt00dot_uetc(ik,it,jt) = emt00dot_uetc(ik,it,jt) + &
                     emt00dot(iexp,ik,it)*emt00dot(iexp,ik,jt)
                emtSdot_uetc(ik,it,jt) = emtSdot_uetc(ik,it,jt) + &
                     emtSdot(iexp,ik,it)*emtSdot(iexp,ik,jt)
                emt00S_uetc(ik,it,jt) = emt00S_uetc(ik,it,jt) + &
                     emt00(iexp,ik,it)*emtS(iexp,ik,jt)
             end do
          end do
       end do
    end do
    emt00_uetc(:,:,:) = emt00_uetc(:,:,:)/nexp
    emtS_uetc(:,:,:) = emtS_uetc(:,:,:)/nexp
    emtV_uetc(:,:,:) = emtV_uetc(:,:,:)/nexp
    emtT_uetc(:,:,:) = emtT_uetc(:,:,:)/nexp
    emtD_uetc(:,:,:) = emtD_uetc(:,:,:)/nexp
    emtP_uetc(:,:,:) = emtP_uetc(:,:,:)/nexp
    emt00dot_uetc(:,:,:) = emt00dot_uetc(:,:,:)/nexp
    emtSdot_uetc(:,:,:) = emtSdot_uetc(:,:,:)/nexp
    emt00S_uetc(:,:,:) = emt00S_uetc(:,:,:)/nexp
    fmt = trim(numcat('(',nt))//'E15.5)' 
    do ik=1,nk
       filename = trim(str_output)//'emt00_uetc_'//trim(IntToStr(ik))//'.dat'
       open(unit=60,file=filename,status='unknown')
       filename = trim(str_output)//'emtS_uetc_'//trim(IntToStr(ik))//'.dat'
       open(unit=61,file=filename,status='unknown')
       filename = trim(str_output)//'emtV_uetc_'//trim(IntToStr(ik))//'.dat'
       open(unit=62,file=filename,status='unknown')
       filename = trim(str_output)//'emtT_uetc_'//trim(IntToStr(ik))//'.dat'
       open(unit=63,file=filename,status='unknown')
       filename = trim(str_output)//'emtD_uetc_'//trim(IntToStr(ik))//'.dat'
       open(unit=64,file=filename,status='unknown')
       filename = trim(str_output)//'emtP_uetc_'//trim(IntToStr(ik))//'.dat'
       open(unit=65,file=filename,status='unknown')
       filename = trim(str_output)//'emt00dot_uetc_'//trim(IntToStr(ik))//'.dat'
       open(unit=66,file=filename,status='unknown')
       filename = trim(str_output)//'emtSdot_uetc_'//trim(IntToStr(ik))//'.dat'
       open(unit=67,file=filename,status='unknown')
       filename = trim(str_output)//'emt00S_uetc_'//trim(IntToStr(ik))//'.dat'
       open(unit=68,file=filename,status='unknown')
       do it=1,nt
          write(60,fmt) emt00_uetc(ik,it,1:nt) 
          write(61,fmt) emtS_uetc(ik,it,1:nt) 
          write(62,fmt) emtV_uetc(ik,it,1:nt) 
          write(63,fmt) emtT_uetc(ik,it,1:nt) 
          write(64,fmt) emtD_uetc(ik,it,1:nt) 
          write(65,fmt) emtP_uetc(ik,it,1:nt) 
          write(66,fmt) emt00dot_uetc(ik,it,1:nt) 
          write(67,fmt) emtSdot_uetc(ik,it,1:nt) 
          write(68,fmt) emt00S_uetc(ik,it,1:nt) 
       end do
       close(60)
       close(61)
       close(62)
       close(63)
       close(64)
       close(65)
       close(66)
       close(67)
       close(68)
    end do
  end subroutine strings_uetc
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  subroutine prepare
    ! One scale model ran between tau_init and tau corresponding to a=10
    implicit none
    type(StrEvolutionVars) SEV
    real(DP) :: tau,taui,tauf,tauend,tol1,tol2,dlnt,tinit,ainit,xl0
    real(DP) :: adot,xl,vel,gam,con,alf
    integer :: it,nvar,ind,j
    integer, parameter :: nw = 3
    real(DP) :: yev(nw),yevpr(nw),c(24),w(nw,9)
    real(DP) dtauda, adtauda, rombint
    external dtauda, adtauda,rombint,fevolve
    !  taui is the earliest conformal time at which one scale model is run   
    taui=tau_init
    !  calculate conformal time at 10 times scale factor today
    tol1=1.d-4
    tauf=rombint(dtauda,0.0d0,10.0d0,tol1)
    !  Set the time grid for the slow functions of time
    !  define tim(it) logarithmically spaced
    dlnt=log(tauf/taui)/dble(nstep00-1)
    do it=1,nstep00
       tim(it)=taui*exp(dble(it-1)*dlnt)
    end do
    tinit=0.1d0*tim(1)      
    ainit=adotrad*tinit
    !  initial correlation length of the network
    xl0=dksi*tinit
    nvar = 3
    yev(1)=adotrad*tinit
    yev(2)=xl0
    yev(3)=vdev
    tau=tinit
    ind=1
    tol1=1.0d-6
    tol2=1.0d-4
    !  evolve and save the paramters of the one-scale model
    !  and their derivatives.
    do  j=1,nstep00
       tauend=tim(j)
       call dverk(SEV,nvar,fevolve,tau,yev,tauend,tol1,ind,c,nw,w)
       call fevolve(SEV,nvar,tauend,yev,yevpr)
       adot=yevpr(1)
       xl=yev(2)
       vel=yev(3)
       ! this is used to store t as a function of tau
       tphys(j)=rombint(adtauda,0.d0,yev(1),tol2)
       gam=adot*tau/yev(1)
       con=alfrad-1.0d0
       alf=1.0d0+con/gam
       !  the scale factor
       xaa(j)=yev(1)     
       if (evolve) then
          !  wiggliness
          evalf(j)=alf    
          !  correlation lenght
          cole(j)=xl    
          !  rms velocity divided by wiggliness
          evv(j)=vel/alf      
       else
          evalf(j) = alfrad
          cole(j) = dksi*tauend
          evv(j) = vdev
       end if
       ! Testing
       write(10,'(10E15.5)') tauend,yev(1),alf,vel,cole(j)/tauend
       ! End testing
    end do
    stop
    ! Spline the one-scale functions of time, also the physical time.      
    call spline(tphys,tim,nstep00,d0lo,d0hi,timpr)
    call spline(tim,tphys,nstep00,d0lo,d0hi,tphyspr)
    call spline(tim,cole,nstep00,d0lo,d0hi,colepr)
    call spline(tim,evalf,nstep00,d0lo,d0hi,evalfpr)      
    call spline(tim,evv,nstep00,d0lo,d0hi,evvpr)     
    call spline(tim,xaa,nstep00,d0lo,d0hi,xaapr)
  end subroutine prepare
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  subroutine mixup
    use am_routines
    implicit none
    real(DP) :: tau,taui,tau2,tauf,dlnt,xl,sum,bnum,off,offdot
    real(DP) :: acteta(ns0),aphi(ns0),apsi(ns0),taum(ns0)
    integer :: m,it,is,nrand,jl,ju,jm,intt
    real(DP) :: u1,u2,y1,y2,cteta,phi,psi,steta,xd1,xd2,xp1,xp2
    real(DP) :: vm,vm2,ovm2,alf
    real(DP) :: rand
    taui = tim(1)
    tf(1) = taui
    if (random_nr) then
       tau2=taui*(1.0d0+ran1_f77(iseed0))
       tauf=tim(nstep00)*(1.0d0+ran1_f77(iseed0))
    else
       tau2=taui*(1.0d0+ranmar())
       tauf=tim(nstep00)*(1.0d0+ranmar())
    end if
    xl = splint_nr(tim,cole,colepr,taui)
    xlden(1)=xl
    dlnt=log(tauf/tau2)/dble(ns1-1)
    ! Testing
    !write(*,*) taui,tau2,tauf,dlnt
    ! End testing
    do m=2,ns1+1
       if (ns1.eq.1) then
          tau=tauf
       else
          tau=tau2*exp(dble(m-2)*dlnt)
       endif
       tf(m)=tau
       xl = splint_nr(tim,cole,colepr,tau)
       xlden(m)=xl
       ! Testing
       !write(*,*) m,tf(m)
       ! End testing
    end do
    do it=1,nstep00
       sum=0.0d0
       tau=tim(it)
       do m=2,ns1
          bnum=1.0d0/xlden(m-1)**3-1.0d0/xlden(m)**3
          call toff(tau,tf(m),off,offdot)
          !*************** Testing *************
          !off=1.0
          !write(*,*) m,bnum
          !************* End testing ***********
          sum=sum+bnum*off        
       end do     
       sum=sum+1.0d0/xlden(ns1+1)**3
       cnorm(it)=1.0d0/(sum*cole(it)**3)
       ! Testing
       !write(*,*) tau,sum,cnorm(it)
       ! End testing
    end do
    !stop
    call spline(tim,cnorm,nstep00,d0lo,d0hi,cnormpr)
    ! Initial time is same for all segments
    taui=tf(1)
    ! final times are also same
    tauf=tf(ns1+1)    
    !  Loop over string segments 
    do m=2,ns1+1
       !  each segment gets a random intitial phase
       if (random_nr) then
          x0k(m)=2.0d0*pi*ran1_f77(iseed0)
       else
          x0k(m)=2.0d0*pi*ranmar()
       end if
       if(irandomv.eq.1) then
          stop 'Not implemented yet'
       else
          nrand=2
          taum(1)=taui
       end if
       taum(nrand)=tauf
       do is=1,nrand-1
          !  Generate directions at random, except for the final time.
          if (random_nr) then
             acteta(is)=2.0d0*ran1_f77(iseed0)-1.0d0
             aphi(is)=2.0d0*pi*ran1_f77(iseed0)
             apsi(is)=2.0d0*pi*ran1_f77(iseed0)
          else
             acteta(is)=2.0d0*ranmar()-1.0d0
             aphi(is)=2.0d0*pi*ranmar()
             apsi(is)=2.0d0*pi*ranmar()
          end if
       enddo
       if(irandomv.eq.1) then
          stop 'Not implemented yet'
       else
          acteta(nrand)=acteta(nrand-1)
          aphi(nrand)=aphi(nrand-1)
          apsi(nrand)=apsi(nrand-1)       
       endif
       if(debug) then
          write(*,*) '*** Random parameters ***'
          write(*,*) 'Tau_f:',taum(1:nrand)
          write(*,*) 'Theta:',acos(acteta(1:nrand))
          write(*,*) 'Phi:',aphi(1:nrand)
          write(*,*) 'Psi:',apsi(1:nrand)
          write(*,*) 'Phase:', x0k(m)
          write(*,*) '************************'
       end if
       do it=1,nstep00
          tau=tim(it)
          jl=0
          ju=nrand+1
          do while (ju-jl.gt.1)
             jm=(ju+jl)/2
             if((taum(nrand).gt.taum(1)).eqv.(tau.ge.taum(jm)))then
                jl=jm
             else
                ju=jm
             endif
          end do
          intt=jl
          u1=(tau-taum(intt))/(taum(intt+1)-taum(intt))
          u2=(taum(intt+1)-tau)/(taum(intt+1)-taum(intt))
          y1=acteta(intt)
          y2=acteta(intt+1)
          cteta=u2*y1+u1*y2
          y1=aphi(intt)
          y2=aphi(intt+1)
          phi=u2*y1+u1*y2	  
          y1=apsi(intt)
          y2=apsi(intt+1)
          psi=u2*y1+u1*y2
          steta=sqrt(1.0d0-cteta**2)
          xp1=steta*sin(phi)
          xp2=-steta*cos(phi)
          xp3(m,it)=cteta
          xd1=cos(psi)*cos(phi)-sin(psi)*sin(phi)*cteta
          xd2=cos(psi)*sin(phi)+sin(psi)*cos(phi)*cteta
          xd3(m,it)=steta*sin(psi)
          !  now the prefactors of S, V and T components
          vm=evv(it)
          vm2=evv(it)**2
          ovm2=1.0d0-vm2
          alf=evalf(it)
          ps(m,it)=0.5d0*(vm2*(3.0d0*xd3(m,it)*xd3(m,it)-1.0d0)-&
               ovm2*(3.0d0*xp3(m,it)*xp3(m,it)-1.0d0)/alf**2)     
          pv(m,it)=vm2*xd1*xd3(m,it)-ovm2*xp1*xp3(m,it)/alf**2
          pt(m,it)=vm2*xd1*xd2-ovm2*xp1*xp2/alf**2  
       end do !it
       call spline(tim,xp3(m,:),nstep00,d0lo,d0hi,xp3pr(m,:))
       call spline(tim,xd3(m,:),nstep00,d0lo,d0hi,xd3pr(m,:))
       call spline(tim,ps(m,:),nstep00,d0lo,d0hi,pspr(m,:))
       call spline(tim,pv(m,:),nstep00,d0lo,d0hi,pvpr(m,:))
       call spline(tim,pt(m,:),nstep00,d0lo,d0hi,ptpr(m,:))
       ! Testing
       !write(*,*) m,splint_nr(tim,xp3(m,:),xp3pr(m,:),2.3d-2)
       !write(*,*) m,splint_nr(tim,xd3(m,:),xd3pr(m,:),2.3d-2)
       ! End testing
    end do !m
  end subroutine mixup
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  ! this subroutine turns the string segments off
  subroutine toff(tau,tauf,off,offdot)
    implicit none
    real(DP), intent(in) :: tau,tauf
    real(DP), intent(out) :: off,offdot
    real(DP) :: var
    var=2.0d0*log(xlf*tauf/tau)/log(xlf)-1.0d0
    if(tau.lt.xlf*tauf) then
       off=1.0d0
       offdot=0.0d0
       return
    endif
    if(tau.ge.tauf) then
       off=0.0d0
       offdot=0.0d0
       return
    endif
    off=0.5d0+0.25d0*(var**3-3.0d0*var)
    offdot=-1.5d0*(var**2-1.0d0)/(tau*log(xlf)) 
  end subroutine toff
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  subroutine strings_spline(SEV,taui,tauf)
    implicit none
    type(StrEvolutionVars) SEV
    real(DP), intent(in) :: taui,tauf
    real(DP) :: dtau,dt,tau,T00,TS,TV,TT,taustart,tauend
    integer :: nvar,ind,i,ik,temp
    integer, parameter :: nw = 1
    real(DP) :: y(nw),yprime(nw),c(24),w(nw,9),tol1
    real(DP) rombint_obj2,sum
    external fevolve,rombint_obj2,fa,dVdtau
    ! Testing
    !write(*,*) taui,tau_init
    ! End testing
    ik = SEV%wk_ix
    ts_f(ik,1)=0.1d0*taui
    if(ts_f(ik,1).lt.tau_init) stop 'One scale model needs an earlier start time'	
    dtau=0.0d0
    temp=1
    do while (ts_f(ik,temp).lt.tauf)
       if(irandomv.eq.1) then
          dtau=0.005d0*ts_f(ik,temp)          
          dt=min(dtau,0.25d0/SEV%wk)
       else
          dtau=0.1d0*ts_f(ik,temp)          
          dt=min(dtau,2.0d0/SEV%wk)      
       endif
       ts_f(ik,temp+1)=ts_f(ik,temp)+dt
       temp=temp+1
    end do
    if (temp.gt.nt_fine) stop 'Need to increase nt_fine'
    nt_f(ik) = temp
    !Testing
    !write(*,*) 'Nt_f:',nt_f(ik),SEV%wk
    ! End testing 
    do i=1,temp
       tau=ts_f(ik,i)
       call strings(SEV%wk,tau,T00,TS,TV,TT) 
       em00(ik,i)=T00
       emS(ik,i)=TS
       emV(ik,i)=TV
       emT(ik,i)=TT
       ! Testing
       !write(*,*) i,T00,TS,TV,TT
       ! End testing
    end do
    ! spline Tmunu
    call spline(ts_f(ik,1:temp),em00(ik,1:temp),temp,d0lo,d0hi,em00pr(ik,1:temp))
    call spline(ts_f(ik,1:temp),emS(ik,1:temp),temp,d0lo,d0hi,emSpr(ik,1:temp))
    call spline(ts_f(ik,1:temp),emV(ik,1:temp),temp,d0lo,d0hi,emVpr(ik,1:temp))
    call spline(ts_f(ik,1:temp),emT(ik,1:temp),temp,d0lo,d0hi,emTpr(ik,1:temp))
    ! find the other 2 scalar mode components from conservation     
    call emtconserve(SEV)
    call spline(ts_f(ik,1:temp),emD(ik,1:temp),nt_f(ik),d0lo,d0hi,emDpr(ik,1:temp)) 
    call spline(ts_f(ik,1:temp),emP(ik,1:temp),nt_f(ik),d0lo,d0hi,emPpr(ik,1:temp)) 
    ind=1
    nvar=1
    tol1=1.0d-8
    taustart=ts_f(ik,1)
    tau=taustart
    as(ik,1)=tau*adotrad
    y(1)=as(ik,1)                       
    do i=2,temp            
       tauend=ts_f(ik,i)
       call dverk(SEV,nvar,fa,tau,y,tauend,tol1,ind,c,nw,w)
       as(ik,i) = y(1)
    end do
    call spline(ts_f(ik,1:temp),as(ik,1:temp),temp,d0lo,d0hi,aspr(ik,1:temp))   
    Vk(ik,1) = 0.0d0
    sum = 0.0
    do i=2,temp     
       taustart=ts_f(ik,i-1)
       tauend=ts_f(ik,i)
       sum = sum + rombint_obj2(ik,dVdtau,taustart,tauend,tol1)
       Vk(ik,i) = -(8.0*pi/get_a(ik,tauend)**2)*sum
     end do
     call spline(ts_f(ik,1:temp),Vk(ik,1:temp),temp,d0lo,d0hi,Vkpr(ik,1:temp))
  end subroutine strings_spline
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  subroutine strings_spline_arr(nk,k_arr,taui,tauf)
    implicit none
    type(StrEvolutionVars) SEV
    integer, intent(in) :: nk
    real(DP), intent(in) :: k_arr(nk),taui,tauf
    integer :: ik
    if (nk.gt.nk_max) stop 'Need to increase nk_max'
    !$OMP PARAllEl DO DEFAUlT(SHARED),SCHEDUlE(DYNAMIC) &
    !$OMP & PRIVATE(ik,SEV)
    do ik=1,nk
       SEV%wk_ix = ik
       SEV%wk = k_arr(ik)
       call strings_spline(SEV,taui,tauf)
       ! Testing
       !write(*,*) ik,SEV%wk_ix,SEV%wk,get_em00(ik,290.0d0),get_emD(ik,290.0d0)
       ! End testing
    end do
    !$OMP END PARAllEl DO
  end subroutine strings_spline_arr
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  subroutine strings(wk,tau,T00,TS,TV,TT)
    implicit none
    real(DP), intent(in) :: wk,tau
    real(DP), intent(out) :: T00,TS,TV,TT
    real(DP) :: vm,xl,alf,snorm,osv
    real(DP) :: xpr3,xdt3,xps,xpv,xpt
    real(DP) :: off,offdot,fnorm,ampl
    real(DP) :: var,phase,bet,sinV,cosV,sinP,cosP,trig,edenre
    real(DP) :: x(nstep00),y(nstep00),ypr(nstep00)
    integer :: m
    !  Components of string stress-energy and the derivatives will be calculated 
    !  for the wave vector k=wk and conformal time tau 
    !  Other two scalar components are evaluated from the conservation
    !  equations in the subroutine fderivs.
    T00=0.0d0
    TS=0.0d0
    TT=0.0d0
    TV=0.0d0  
    vm = splint_nr(tim,evv,evvpr,tau)
    xl = splint_nr(tim,cole,colepr,tau)
    alf = splint_nr(tim,evalf,evalfpr,tau)
    snorm = splint_nr(tim,cnorm,cnormpr,tau)
    osv = 1.0d0/sqrt(1.0d0-vm*vm)
    do m=2,ns1+1
       x(:) = tim(:)
       y(:) = xp3(m,:)
       ypr(:) = xp3pr(m,:)
       xpr3 = splint_nr(x,y,ypr,tau)
       y(:) = xd3(m,:)
       ypr(:) = xd3pr(m,:)
       xdt3 = splint_nr(x,y,ypr,tau)
       xps = splint_nr(tim,ps(m,:),pspr(m,:),tau)
       xpv = splint_nr(tim,pv(m,:),pvpr(m,:),tau)
       xpt = splint_nr(tim,pt(m,:),ptpr(m,:),tau)
       call toff(tau,tf(m),off,offdot)
       fnorm=sqrt(snorm/xlden(m-1)**3-&
            snorm/xlden(m)**3)
       if(m.eq.ns1+1) then     
          fnorm=sqrt(snorm/xlden(m)**3)
       endif
       !*************** Testing *************
       if(single_string) then
          off = 1.0 
          fnorm = 1.0
       end if
       !write(*,*) m
       !write(*,'(10E15.5)') tau,fnorm
       !************* End testing ***********
       ampl=fnorm*off
       if(ampl.le.1.0d-13) continue
       var=0.5d0*wk*xpr3*xl	  
       phase=x0k(m)+wk*xdt3*vm*tau	 
       bet=0.5d0*wk*xpr3
       sinV=sin(var)
       cosV=cos(var)
       sinP=sin(phase)
       cosP=cos(phase)
       trig=alf*osv*sinV*cosP/bet 
       ! The factor of sqrt(2) in the equations below is compensating
       ! for the fact that we evaluate only the real part of the Fourier
       ! transform of the stress-energy tensor. 
       edenre=sqrt(2.0d0)*(ampl*trig)
       T00=T00+edenre	  
       TS=TS+xps*edenre
       TV=TV+xpv*edenre	  	  
       TT=TT+xpt*edenre
       ! Testing
       !write(*,*) m,edenre,xps,xpv,xpt
       ! End testing
    end do
  end subroutine strings
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  subroutine emtconserve(SEV)
    implicit none
    type(StrEvolutionVars) SEV
    real(DP) :: wk,tau,taustart,tauend,a,a2,tol1
    integer :: nvar,ind,j,ik
    integer, parameter :: nw = 2
    real(DP) :: y(nw),yprime(nw),c(24),w(nw,9)
    external fstrings
    ik=SEV%wk_ix
    wk=SEV%wk
    taustart=ts_f(ik,1)
    tau=taustart
    a=tau*adotrad
    a2=a*a
    ! Initial values of emtP and emt D set to 0   
    emP(ik,1)=0.0d0
    emD(ik,1)=0.0d0
    y(2)=0.0d0
    y(1)=a                
    ind=1
    nvar=2	
    if (wk.lt.0.2) then 
       tol1=1.0d-8
    else
       tol1=1.0d-6
    endif        
!  Begin timestep loop.
    do j=2,nt_f(ik)              
       tauend=ts_f(ik,j)
       if(wk*tau.le.tkmax) then
          call dverk(SEV,nvar,fstrings,tau,y,tauend,tol1,ind,c,nw,w)
          call fstrings(SEV,nvar,tau,y,yprime)
          emP(ik,j)=SEV%emtP
          emD(ik,j)=SEV%emtD
       else
          emP(ik,j)=0.0d0
          emD(ik,j)=0.0d0	       
       endif
    end do
  end subroutine emtconserve
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  function get_em00(ik,tau) 
    implicit none
    integer :: ik
    real(DP) :: get_em00,tau
    get_em00 = splint_nr(ts_f(ik,1:nt_f(ik)),em00(ik,1:nt_f(ik)),em00pr(ik,1:nt_f(ik)),tau)
  end function get_em00
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  function get_em00_dot(ik,tau) 
    implicit none
    integer :: ik
    real(DP) :: get_em00_dot,tau
    get_em00_dot = splint_deriv_nr(ts_f(ik,1:nt_f(ik)),em00(ik,1:nt_f(ik)),em00pr(ik,1:nt_f(ik)),tau)
  end function get_em00_dot
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  function get_emS(ik,tau) 
    implicit none
    integer :: ik
    real(DP) :: get_emS,tau
    get_emS = splint_nr(ts_f(ik,1:nt_f(ik)),emS(ik,1:nt_f(ik)),emSpr(ik,1:nt_f(ik)),tau)
  end function get_emS
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  function get_emS_dot(ik,tau) 
    implicit none
    integer :: ik
    real(DP) :: get_emS_dot,tau
    get_emS_dot = splint_deriv_nr(ts_f(ik,1:nt_f(ik)),emS(ik,1:nt_f(ik)),emSpr(ik,1:nt_f(ik)),tau)
  end function get_emS_dot
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  function get_emV(ik,tau) 
    implicit none
    integer :: ik
    real(DP) :: get_emV,tau
    get_emV = splint_nr(ts_f(ik,1:nt_f(ik)),emV(ik,1:nt_f(ik)),emVpr(ik,1:nt_f(ik)),tau)
  end function get_emV
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  function get_emT(ik,tau) 
    implicit none
    integer :: ik
    real(DP) :: get_emT,tau
    get_emT = splint_nr(ts_f(ik,1:nt_f(ik)),emT(ik,1:nt_f(ik)),emTpr(ik,1:nt_f(ik)),tau)
  end function get_emT
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  function get_emD(ik,tau) 
    implicit none
    integer :: ik
    real(DP) :: get_emD,tau
    get_emD = splint_nr(ts_f(ik,1:nt_f(ik)),emD(ik,1:nt_f(ik)),emDpr(ik,1:nt_f(ik)),tau)
  end function get_emD
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  function get_emP(ik,tau) 
    implicit none
    integer :: ik
    real(DP) :: get_emP,tau
    get_emP = splint_nr(ts_f(ik,1:nt_f(ik)),emP(ik,1:nt_f(ik)),emPpr(ik,1:nt_f(ik)),tau)
  end function get_emP
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  function get_a(ik,tau) 
    implicit none
    integer :: ik
    real(DP) :: get_a,tau
    get_a = splint_nr(ts_f(ik,1:nt_f(ik)),as(ik,1:nt_f(ik)),aspr(ik,1:nt_f(ik)),tau)
  end function get_a
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  function get_Vk(ik,tau) 
    implicit none
    integer :: ik
    real(DP) :: get_Vk,tau
    get_Vk = splint_nr(ts_f(ik,1:nt_f(ik)),Vk(ik,1:nt_f(ik)),Vkpr(ik,1:nt_f(ik)),tau)
  end function get_Vk
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
end module string

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function adtauda(a)
  use string
  implicit none
  real(DP) adtauda,dtauda,a
  external dtauda
  adtauda = dtauda(a)*a
end function adtauda
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
subroutine fstrings(SEV,n,x,y,yprime)
  use ModelParams
  use LambdaGeneral
  use string
  implicit none
  type(StrEvolutionVars) SEV
  integer n,ik
  real(DP) :: x,y(n),yprime(n)
  real(DP) :: a,a2,ak2,tau,adotoa
  real(DP) :: grhob_t,grhor_t,grhoc_t,grhog_t,grhov_t,grho
  real(DP) :: emt00,emtS,emt00dot,emtSdot,emtDdot,emtD,emtP
  tau = x
  a=y(1)
  a2=a*a
  ak2=SEV%wk*SEV%wk
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
   ik = SEV%wk_ix
   emt00 = get_em00(ik,tau)   
   emtS = get_emS(ik,tau)  
   emt00dot = get_em00_dot(ik,tau) 
   emtSdot = get_emS(ik,tau) 
   emtD=y(2)
   yprime(1)=adotoa*a
   emtP=(emtD-emt00dot)/adotoa-emt00 
!  String conservation equation
   emtDdot=-2.0d0*adotoa*emtD-(ak2/3.0d0)*(emtP+2.0d0*emtS)
   yprime(2)=emtDdot   
   SEV%emtD = emtD
   SEV%emtP = emtP
end subroutine fstrings
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
subroutine fevolve(SEV,n,x,y,yprime)
  use ModelParams
  use LambdaGeneral
  use string
  implicit none
  type(StrEvolutionVars) SEV
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
   ct=(cr+a*g*cm)/(1.0d0+a*g)
   fk=(fkr+a*g*fkm)/(1.0d0+a*g)
   xlprime=adotoa*xl*xv2+0.5d0*ct*xv
   xvprime=(1.0d0-xv2)*(fk/xl-2.0d0*adotoa*xv)  
   yprime(2)=xlprime
   yprime(3)=xvprime
end subroutine fevolve
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
subroutine fa(SEV,n,x,y,yprime)
  use ModelParams
  use LambdaGeneral
  use string
  implicit none
  type(StrEvolutionVars) SEV
  integer n,temp,ik
  real(DP) :: x,y(n),yprime(n)
  real(DP) :: a,a2,tau,adotoa
  real(DP) :: grhob_t,grhor_t,grhoc_t,grhog_t,grhov_t,grho
  logical :: nomatter =.false.
  tau = x
  a=y(1)
  a2=a*a
  if (nomatter) then
     grhov = grhov+grhoc+grhob
     grhoc = 0.0
     grhob = 0.0
  end if
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
end subroutine fa
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function dVdtau(ik,tau)
  use string
  implicit none
  integer ik
  real(DP) dVdtau,tau
  dVdtau = get_emV(ik,tau)*get_a(ik,tau)**2
end function dVdtau

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function rombint_obj2(obj,f,a,b,tol, maxit)
        use Precision
!  Rombint returns the integral from a to b of using Romberg integration.
!  The method converges provided that f(x) is continuous in (a,b).
!  f must be real(dl) and must be declared external in the calling
!  routine.  tol indicates the desired relative accuracy in the integral.
!
        implicit none
        integer, intent(in), optional :: maxit
        integer :: MAXITER=20
        integer, parameter :: MAXJ=5
        dimension g(MAXJ+1)
        integer obj !dummy
        real(dl) f
        external f
        real(dl) :: rombint_obj2
        real(dl), intent(in) :: a,b,tol
        integer :: nint, i, k, jmax, j
        real(dl) :: h, gmax, error, g, g0, g1, fourj
!
        if (present(maxit)) then
            MaxIter = maxit
        end if
        h=0.5d0*(b-a)
        gmax=h*(f(obj,a)+f(obj,b))
        g(1)=gmax
        nint=1
        error=1.0d20
        i=0
10        i=i+1
          if (i.gt.MAXITER.or.(i.gt.5.and.abs(error).lt.tol)) &
            go to 40
!  Calculate next trapezoidal rule approximation to integral.
          g0=0._dl
            do 20 k=1,nint
            g0=g0+f(obj,a+(k+k-1)*h)
20        continue
          g0=0.5d0*g(1)+h*g0
          h=0.5d0*h
          nint=nint+nint
          jmax=min(i,MAXJ)
          fourj=1._dl
            do 30 j=1,jmax
!  Use Richardson extrapolation.
            fourj=4._dl*fourj
            g1=g0+(g0-g(j))/(fourj-1._dl)
            g(j)=g0
            g0=g1
30        continue
          if (abs(g0).gt.tol) then
            error=1._dl-gmax/g0
          else
            error=gmax
          end if
          gmax=g0
          g(jmax+1)=g0
        go to 10
40      rombint_obj2=g0
        if (i.gt.MAXITER.and.abs(error).gt.tol)  then
          write(*,*) 'Warning: Rombint failed to converge; '
          write (*,*)'integral, error, tol:', rombint_obj2,error, tol
        end if
        
        end function rombint_obj2
