! Code to compute analytic form of USM string UETC
! AM, Nottingham, July 2012
! Uses NAG libraries to compute Sine integral and spherical Bessel function 
! Uses LAPACK routine for eigenvector decomposition
! Latest version Oct 30 2012

module uetc 
  USE nag_library, ONLY : nag_wp,s13adf,s17def,s14aaf
  use am_routines
  implicit none
  integer, parameter :: DP = kind(1.0D0)
  integer, parameter :: QP = kind(1.0Q0)
  real(DP), parameter :: pi = 3.1415926535897932384626433832795d0
  ! Whether to use gridded I1 and I4
  logical :: use_grid = .true.
  real(DP), parameter :: taustart_string = 1.0d0
  logical :: do_strings = .false.
  logical :: do_string_source = .false.
  integer :: scaling_option = 2
  integer, parameter :: min_terms = 10
  real, parameter :: scale_terms = 3.0
  real(DP), parameter :: xmin = 0.05d0
  real(DP), parameter :: xmax = 20.0d0
  real(DP), parameter :: etcmin = 0.001d0
  real(DP), parameter :: xapr = 1.0d0
  real(DP) :: weighting = 0.0d0
  real(DP) :: ktau_min,ktau_max,dktau,dx_grid
  integer :: nktau
  real(DP), allocatable :: ktau_array(:),log_ktau_array(:)    
  real(DP), allocatable :: ss00_evec(:,:,:),ss00_eval(:,:),ss00_evec_pr(:,:,:)
  real(DP), allocatable :: ss_evec(:,:,:),ss_eval(:,:),ss_evec_pr(:,:,:)
  real(DP), allocatable :: vv_evec(:,:,:),vv_eval(:,:),vv_evec_pr(:,:,:)
  real(DP), allocatable :: tt_evec(:,:,:),tt_eval(:,:),tt_evec_pr(:,:,:) 
  Type StringParams
     real(DP) :: xi,alpha,mu,v,L
  end type StringParams 
  real(DP), parameter :: d0lo=1.0d40
  real(DP), parameter :: d0hi=1.0d40 
  character(len=100) :: str_output
  character(len=100) :: uetc_file='test_uetc'
  integer :: nint
  !------- Grid parameters -------
  logical :: do_grid = .false.
  integer :: nx_grid = 800
  real(DP) :: x_min_grid = 1.0d-5
  real(DP) :: x_max_grid = 1.0d4
  !-------------------------------
  logical :: do_init_grid = .true.
  real(DP), allocatable :: x_array_grid(:),log_x_array_grid(:)  
  real, allocatable :: I1_grid(:,:,:),I4_grid(:,:,:)
  character(len=100) :: grid_root='/home/adammoss/work/code/string/grid/'
  Type(StringParams) :: SPR_rad,SPR_mat
  logical :: output_uetc = .false.
  integer :: nmodes
  logical :: do_diagonalize = .false.
  logical :: no_evolve
  real(DP) :: tau_s
  integer :: uetc_feedback = 1
  logical :: do_num = .true.
contains
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  subroutine init_string(SPR)
    implicit none
    Type(StringParams) :: SPR
    ! Initialize string model parameters	
    ! the string tension (the most important parameter)
    !SPR%mu=1.1d-6 
    SPR%mu = 1.0d0
    ! Wiggliness                   
    SPR%alpha=1.9d0
    ! v is the initial rms string velocity (over all scales)
    SPR%v=0.65d0
    ! dksi is the initial correlation length/initial conformal time
    SPR%xi=0.13d0    
    ! xlf 
    SPR%L = 0.95d0
  end subroutine init_string
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine write_string_params(SPR)
    Type(StringParams), intent(in) :: SPR
    write(*,*) '************* String model parameters *************'
    write(*,*) 'String tension:',SPR%mu
    write(*,*) 'Initial RMS velocity:',SPR%v
    write(*,*) 'Initial correlation length/tau:',SPR%xi
    write(*,*) 'Initial wiggliness:',SPR%alpha
    write(*,*) 'Rate of string decay:',SPR%L
    write(*,*) '***************************************************'
  end subroutine write_string_params
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine deallocate_string()
    implicit none
    deallocate(ktau_array,log_ktau_array)
    deallocate(ss00_evec,ss00_eval,ss00_evec_pr)
    deallocate(ss_evec,ss_eval,ss_evec_pr)
    deallocate(vv_evec,vv_eval,vv_evec_pr)
    deallocate(tt_evec,tt_eval,tt_evec_pr)
  end subroutine deallocate_string
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine get_SPR(tau,SPR)
    implicit none
    real(DP), intent(in) :: tau
    Type(StringParams) :: SPR
    real(DP) :: g,ft
    g = 0.002d0
    ft = (1.0d0-tanh(g*(tau-tau_s)))/(1.0d0+tanh(g*tau_s))
    SPR%mu = SPR_rad%mu
    SPR%xi = SPR_rad%xi*ft + SPR_mat%xi*(1.0d0-ft)
    SPR%alpha = SPR_rad%alpha
    SPR%v = SPR_rad%v
    SPR%L = SPR_rad%L
  end subroutine get_SPR
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine diagonalize_uetc(nk,k_arr)
    ! Eigenvectors of UETC times weighting factor as function of ktau. This is achieved by setting k=1/Mpc in the UETC calculations
    use AMLutils
    implicit none
    integer, intent(in) :: nk
    real(DP), intent(in) :: k_arr(1:nk)
    real(DP) :: tau1,tau2,uetc_val
    real(DP), allocatable :: ss00array(:,:),ssarray(:,:),sscrossarray(:,:)
    real(DP), allocatable :: vvarray(:,:),ttarray(:,:)
    real(DP), allocatable :: scalar_evec_temp(:,:),scalar_eval_temp(:)
    real(DP), allocatable :: ss00array_recon(:,:),ssarray_recon(:,:)
    real(DP), allocatable :: vv_evec_temp(:,:),vv_eval_temp(:),vvarray_recon(:,:) 
    real(DP), allocatable :: tt_evec_temp(:,:),tt_eval_temp(:),ttarray_recon(:,:)
    real(DP), allocatable :: ss00_evec_temp(:,:),ss00_eval_temp(:)
    real(DP), allocatable :: ss_evec_temp(:,:),ss_eval_temp(:)
    real(DP), allocatable :: wtau(:,:)
    integer i,j,k,ii,status,n
    real(DP) :: evec,devec,sevec(2),dsevec(2),uetc_arr(5)
    real(DP) :: tau1test,tau2test,x1test,x2test,rhotest,I1_val,I4_val
    Type(StringParams) :: SPR

    if(no_evolve) then
       nint = 1
    else
       nint = nk
    end if

    if (uetc_feedback.gt.0) then
       write(*,*) 'Computing UETC'
       write(*,*) 'N inter: ',nint
       write(*,*) 'Output UETC:',output_uetc
       write(*,*) 'k tau weighting:',weighting
       write(*,*) 'UETC n grid:',nktau
       write(*,*) 'ktau_min:',ktau_min 
       write(*,*) 'ktau_max:',ktau_max 
    end if
    n = nktau

    allocate(ktau_array(n),log_ktau_array(n))

    allocate(ss00_evec(nint,nmodes,n),ss00_eval(nint,nmodes),ss00_evec_pr(nint,nmodes,n))
    allocate(ss_evec(nint,nmodes,n),ss_eval(nint,nmodes),ss_evec_pr(nint,nmodes,n))
    allocate(vv_evec(nint,nmodes,n),vv_eval(nint,nmodes),vv_evec_pr(nint,nmodes,n))
    allocate(tt_evec(nint,nmodes,n),tt_eval(nint,nmodes),tt_evec_pr(nint,nmodes,n))

    allocate(ss00array(n,n),ssarray(n,n),sscrossarray(n,n))
    allocate(vvarray(n,n),ttarray(n,n))
    allocate(scalar_evec_temp(2*n,2*n),scalar_eval_temp(2*n))
    allocate(ss00array_recon(n,n),ssarray_recon(n,n))
    allocate(vv_evec_temp(n,n),vv_eval_temp(n))
    allocate(vvarray_recon(n,n),ttarray_recon(n,n))
    allocate(tt_evec_temp(n,n),tt_eval_temp(n))
    allocate(ss00_evec_temp(n,n),ss00_eval_temp(n))
    allocate(ss_evec_temp(n,n),ss_eval_temp(n))
    allocate(wtau(n,n))

    dktau = log(ktau_max/ktau_min)/dble(n-1)
    do i=1,n	    
       ktau_array(i)=ktau_min*exp(dktau*(i-1))
       log_ktau_array(i)=log(ktau_array(i))
       !write(*,*) i,ktau_array(i)
    end do

    if (uetc_feedback.gt.0) then
       call write_string_params(SPR_rad)
       call write_string_params(SPR_mat)
    end if

    if (do_init_grid .and. use_grid) then
       if(do_grid) then
          call init_grid()
       else
          call read_grid()
       end if
       ! Testing 
       !x1test = 100.0d0
       !x2test = 20.0d0
       !rhotest = 10.0d0
       !call get_grid(x1test,x2test,rhotest,I1_val,I4_val,status)
       !write(*,*) I1_int_num(x1test,x2test,rhotest)
       !write(*,*) I1_val
       !write(*,*) I4_int_num(x1test,x2test,rhotest)
       !write(*,*) I4_val
       ! End testing  
    end if

    if(.false.) then
       ! Testing
       tau1test = 50.0d0
       tau2test = 40.0d0
       call get_uetc(tau1test,tau2test,1.0d0,uetc_arr)
       write(*,*) uetc_arr(:)
       stop
    end if

    do ii=1,nint 
       if (uetc_feedback.gt.0) write(*,*) 'K mode:',ii
       do i=1,n
          !$OMP PARAllEl DO DEFAUlT(SHARED),SCHEDUlE(DYNAMIC) &
          !$OMP & PRIVATE(j,tau1,tau2,uetc_arr)
          do j=1,n 
             if (j.gt.i) cycle
             if(no_evolve) then
                tau1 = ktau_array(i)
                tau2 = ktau_array(j)
                call get_uetc(tau1,tau2,1.0d0,uetc_arr)
             else
                tau1 = ktau_array(i)/k_arr(ii)
                tau2 = ktau_array(j)/k_arr(ii)
                call get_uetc(tau1,tau2,k_arr(ii),uetc_arr)
             end if
             wtau(i,j) = (ktau_array(i)*ktau_array(j))**weighting
             ss00array(i,j) = sqrt(tau1*tau2)*uetc_arr(1)
             ssarray(i,j) = sqrt(tau1*tau2)*uetc_arr(2)
             vvarray(i,j) = sqrt(tau1*tau2)*uetc_arr(3)		  
             ttarray(i,j) = sqrt(tau1*tau2)*uetc_arr(4)
             sscrossarray(i,j) = sqrt(tau1*tau2)*uetc_arr(5)
          end do
          !$OMP END PARAllEl DO
       end do
       do i=1,n
          do j=1,n 
             if (j.gt.i) ss00array(i,j)=ss00array(j,i)
             if (j.gt.i) ssarray(i,j)=ssarray(j,i)
             if (j.gt.i) vvarray(i,j)=vvarray(j,i)
             if (j.gt.i) ttarray(i,j)=ttarray(j,i)
             if (j.gt.i) sscrossarray(i,j)=sscrossarray(j,i)
             if (j.gt.i) wtau(i,j)=wtau(j,i)
          end do
       end do
       if(output_uetc .and. ii.eq.1) then 
          open(unit=50,file =trim(uetc_file)//'_ktau.dat')
          do i=1,n	    
             write(50,'(1E15.5)') ktau_array(i)
          end do
          close(50)
          open(unit=50,file = trim(uetc_file)//'_ss00.dat')
          open(unit=51,file = trim(uetc_file)//'_ss.dat')	
          open(unit=52,file = trim(uetc_file)//'_vv.dat')
          open(unit=53,file = trim(uetc_file)//'_tt.dat')
          open(unit=54,file = trim(uetc_file)//'_sscross.dat')
          do i=1,n
             write(50,'(500E15.5)') ss00array(i,:)
             write(51,'(5000E15.5)') ssarray(i,:)
             write(52,'(5000E15.5)') vvarray(i,:)
             write(53,'(5000E15.5)') ttarray(i,:)
             write(54,'(5000E15.5)') sscrossarray(i,:)
          end do
          close(50)
          close(51)
          close(52)
          close(53)
          close(54)
       end if
       do i=1,n
          do j=1,n
             vv_evec_temp(i,j) = wtau(i,j)*vvarray(i,j)
             tt_evec_temp(i,j) = wtau(i,j)*ttarray(i,j)
             scalar_evec_temp(i,j)=wtau(i,j)*ss00array(i,j)
             scalar_evec_temp(i+n,j+n)=wtau(i,j)*ssarray(i,j) 
             scalar_evec_temp(i,j+n)=wtau(i,j)*sscrossarray(i,j)
             scalar_evec_temp(j+n,i)=wtau(i,j)*scalar_evec_temp(i,j+n)         
             ss00_evec_temp(i,j) = wtau(i,j)*ss00array(i,j)
             ss_evec_temp(i,j) = wtau(i,j)*ssarray(i,j)
             !!write(*,*) ss00_evec_temp(i,j),ss_evec_temp(i,j),vv_evec_temp(i,j),tt_evec_temp(i,j),wtau(i,j)
          end do
       end do
       !write(*,*) 'Diagonalizing UETC'
       call Matrix_Diagonalize(scalar_evec_temp,scalar_eval_temp,2*n)   
       call Matrix_Diagonalize(vv_evec_temp,vv_eval_temp,n)
       call Matrix_Diagonalize(tt_evec_temp,tt_eval_temp,n)
       do i=1,nmodes
          ss00_eval(ii,i) = scalar_eval_temp(2*n-i+1)
          ss00_evec(ii,i,:) = scalar_evec_temp(1:n,2*n-i+1)
          call spline_nr(log(ktau_array),ss00_evec(ii,i,:),d0lo,d0hi,ss00_evec_pr(ii,i,:))
          ss_eval(ii,i) = ss00_eval(ii,i)
          ss_evec(ii,i,:) = scalar_evec_temp(n+1:2*n,2*n-i+1)
          call spline_nr(log(ktau_array),ss_evec(ii,i,:),d0lo,d0hi,ss_evec_pr(ii,i,:))
          vv_eval(ii,i) = vv_eval_temp(n-i+1)
          vv_evec(ii,i,:) = vv_evec_temp(:,n-i+1)
          call spline_nr(log(ktau_array),vv_evec(ii,i,:),d0lo,d0hi,vv_evec_pr(ii,i,:))
          tt_eval(ii,i) = tt_eval_temp(n-i+1)
          tt_evec(ii,i,:) = tt_evec_temp(:,n-i+1)
          call spline_nr(log(ktau_array),tt_evec(ii,i,:),d0lo,d0hi,tt_evec_pr(ii,i,:))
       end do

       if(output_uetc.and. ii.eq.1) then
          ! Test reconstructing UETC from a number of eigenmodes
          ss00array_recon(:,:) = 0.0d0
          ssarray_recon(:,:) = 0.0d0
          vvarray_recon(:,:) = 0.0d0
          ttarray_recon(:,:) = 0.0d0
          do k=1,nmodes 
             do i=1,n
                do j=1,n
                   ss00array_recon(i,j) = ss00array_recon(i,j) + ss00_eval(ii,k)*ss00_evec(ii,k,i)*ss00_evec(ii,k,j)/wtau(i,j)
                   ssarray_recon(i,j) = ssarray_recon(i,j) + ss_eval(ii,k)*ss_evec(ii,k,i)*ss_evec(ii,k,j)/wtau(i,j)
                   vvarray_recon(i,j) = vvarray_recon(i,j) + vv_eval(ii,k)*vv_evec(ii,k,i)*vv_evec(ii,k,j)/wtau(i,j)
                   ttarray_recon(i,j) = ttarray_recon(i,j) + tt_eval(ii,k)*tt_evec(ii,k,i)*tt_evec(ii,k,j)/wtau(i,j)
                end do
             end do
          end do
          open(unit=50,file = trim(uetc_file)//'_ss00_recon.dat')
          open(unit=51,file = trim(uetc_file)//'_ss_recon.dat')
          open(unit=52,file = trim(uetc_file)//'_vv_recon.dat')
          open(unit=53,file = trim(uetc_file)//'_tt_recon.dat')
          do i=1,n
             write(50,'(5000E15.5)') ss00array_recon(i,:)
             write(51,'(5000E15.5)') ssarray_recon(i,:)
             write(52,'(5000E15.5)') vvarray_recon(i,:)
             write(53,'(5000E15.5)') ttarray_recon(i,:)
          end do
          close(50)
          close(51)
          close(52)
          close(53)
          open(unit=50,file = trim(uetc_file)//'_'//trim(IntToStr(ii))//'_evec.dat')
          do i=1,n
             write(50,'(5000E15.5)') vv_evec(ii,i,:)
          end do
          close(50)
       end if
    end do

    deallocate(ss00array,ssarray,sscrossarray)
    deallocate(vvarray,ttarray)
    deallocate(scalar_evec_temp,scalar_eval_temp)
    deallocate(ss00array_recon,ssarray_recon)
    deallocate(vv_evec_temp,vv_eval_temp)
    deallocate(vvarray_recon,ttarray_recon)
    deallocate(tt_evec_temp,tt_eval_temp)
    deallocate(ss00_evec_temp,ss00_eval_temp)
    deallocate(ss_evec_temp,ss_eval_temp)
    deallocate(wtau)

    do_diagonalize = .false.
    if (uetc_feedback.gt.0) write(*,*) 'Done computing UETC'

  end subroutine diagonalize_uetc
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  subroutine init_grid()
    ! Initialize I1 and I4 grids 
    implicit none
    real(DP) :: x1,x2,rho
    integer :: i,j,k
    integer status,unit,blocksize,bitpix,naxis,naxes(3)
    integer group,fpixel,nelements
    character (len=100) filename
    logical simple,extend

    write(*,*) 'Initializing I1 and I4 grids'
    allocate(x_array_grid(nx_grid),log_x_array_grid(nx_grid))
    allocate(I1_grid(nx_grid,nx_grid,nx_grid),I4_grid(nx_grid,nx_grid,nx_grid))
    dx_grid = log(x_max_grid/x_min_grid)/dble(nx_grid-1)
    do i=1,nx_grid	    
       x_array_grid(i)=x_min_grid*exp(dx_grid*(i-1))
       log_x_array_grid(i)=log(x_array_grid(i))
       !write(*,*) i,x_array_grid(i)
    end do

    do i=1,nx_grid
       write(*,*) i
       do j=1,nx_grid
          !$OMP PARAllEl DO DEFAUlT(SHARED),SCHEDUlE(DYNAMIC) &
          !$OMP & PRIVATE(k,x1,x2,rho)
          do k=1,j
             rho = x_array_grid(i)
             x1 = x_array_grid(j)
             x2 = x_array_grid(k)
             I1_grid(j,k,i) = I1_int_num(x1,x2,rho)
             if (k.ne.j) I1_grid(k,j,i) = I1_grid(j,k,i)
             I4_grid(j,k,i) = I4_int_num(x1,x2,rho)
             if (k.ne.j) I4_grid(k,j,i) = I4_grid(j,k,i)
          end do
          !$OMP END PARAllEl DO
       end do
    end do

    status = 0 
    filename = trim(grid_root)//'I1_grid.fits'
    call delete_file(filename,status)
    call ftgiou(unit,status)
    blocksize=1
    call ftinit(unit,filename,blocksize,status)
    simple=.true.
    bitpix=-32
    extend=.true.
    naxis=3
    naxes(1:3)=nx_grid
    call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
    group=1
    fpixel=1
    nelements=naxes(1)**3
    call ftppre(unit,group,fpixel,nelements,I1_grid,status)
    call ftpkyd(unit,'x_min_grid',x_min_grid,5,'Minimum grid size in x=k \xi \tau',status)
    call ftpkyd(unit,'x_max_grid',x_max_grid,5,'Maximum grid size in x=k \xi \tau',status)
    call ftclos(unit,status)
    call ftfiou(unit,status)

    status = 0 
    filename = trim(grid_root)//'I4_grid.fits'
    call delete_file(filename,status)
    call ftgiou(unit,status)
    blocksize=1
    call ftinit(unit,filename,blocksize,status)
    simple=.true.
    bitpix=-32
    extend=.true.
    naxis=3
    naxes(1:3)=nx_grid
    call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
    group=1
    fpixel=1
    nelements=naxes(1)**3
    call ftppre(unit,group,fpixel,nelements,I4_grid,status)
    call ftpkyd(unit,'x_min_grid',x_min_grid,5,'Minimum grid size in x=k \xi \tau',status)
    call ftpkyd(unit,'x_max_grid',x_max_grid,5,'Maximum grid size in x=k \xi \tau',status)
    call ftclos(unit,status)
    call ftfiou(unit,status)

    do_init_grid =.false.
    write(*,*) 'Done'

  end subroutine init_grid
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  subroutine read_grid()
    ! Read I1 and I4 grids
    implicit none
    integer status,unit,readwrite,blocksize,naxes(3),nfound
    integer group,fpixel,nelements,i
    logical anyf
    character(len=100) filename
    character(len=80) comment
    real(DP) keyval

    write(*,*) 'Reading I1 and I4 grids'
    filename=trim(grid_root)//'I1_grid.fits'
    status=0
    call ftgiou(unit,status)
    readwrite=0
    call ftopen(unit,filename,readwrite,blocksize,status)
    call ftgknj(unit,'NAXIS',1,3,naxes,nfound,status)
    call ftgkyd(unit,'x_min_grid',x_min_grid,comment,status)
    call ftgkyd(unit,'x_max_grid',x_max_grid,comment,status)
    nx_grid = naxes(1)
    allocate(x_array_grid(nx_grid),log_x_array_grid(nx_grid))
    allocate(I1_grid(nx_grid,nx_grid,nx_grid),I4_grid(nx_grid,nx_grid,nx_grid))
    dx_grid = log(x_max_grid/x_min_grid)/dble(nx_grid-1)
    do i=1,nx_grid
       x_array_grid(i)=x_min_grid*exp(dx_grid*(i-1))
       log_x_array_grid(i)=log(x_array_grid(i))
       !write(*,*) i,x_array_grid(i),log_x_array_grid(i)
    end do
    nelements=naxes(1)**3
    group=1
    fpixel=1
    call ftgpve(unit,group,fpixel,nelements,0,I1_grid,anyf,status)
    call ftclos(unit,status)
    call ftfiou(unit,status)

    filename=trim(grid_root)//'I4_grid.fits'
    status=0
    call ftgiou(unit,status)
    readwrite=0
    call ftopen(unit,filename,readwrite,blocksize,status)
    group=1
    fpixel=1
    call ftgpve(unit,group,fpixel,nelements,0,I4_grid,anyf,status)
    call ftclos(unit,status)
    call ftfiou(unit,status)

    do_init_grid =.false.
    write(*,*) 'Done'

  end subroutine read_grid
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  subroutine get_grid(x1,x2,rho,I1_val,I4_val,status)
    ! Use trilinear interpolation to obtain I1 and I4 
    implicit none
    real(DP), intent(in) :: x1,x2,rho
    real(DP), intent(out) :: I1_val,I4_val
    integer, intent(out) :: status
    integer :: nx1,nx2,nrho
    integer :: nrho_1,nrho_2,nx1_1,nx1_2,nx2_1,nx2_2
    real(DP) :: log_x1,log_x2,log_rho   
    real(DP) :: c00,c10,c01,c11,c0,c1

    status = 0
    I1_val = 0.0d0
    I4_val = 0.0d0

    log_x1 = log(x1)
    log_x2 = log(x2)
    log_rho = log(rho)

    nrho_1 = floor((log_rho-log_x_array_grid(1))/dx_grid+1)
    if (nrho_1.eq.0) then 
       !write(*,*) 'Warning, nrho_1=0 for rho:',rho
       status = 1
       return
    end if
    nrho_2 = nrho_1 + 1
    nx1_1 = floor((log_x1-log_x_array_grid(1))/dx_grid+1)
    if (nx1_1.eq.0) then 
       !write(*,*) 'Warning, nx1_1=0 for x1:',x1
       status = 1
       return
    end if
    nx1_2 = nx1_1 + 1
    nx2_1 = floor((log_x2-log_x_array_grid(1))/dx_grid+1)
    if (nx2_1.eq.0) then 
       !write(*,*) 'Warning, nx2_1=0 for x2:',x2
       status = 1
       return
    end if
    nx2_2 = nx2_1 + 1

    c00 =  (I1_grid(nx1_1,nx2_1,nrho_1)*(log_x_array_grid(nrho_2)-log_rho) + &
            I1_grid(nx1_1,nx2_1,nrho_2)*(log_rho-log_x_array_grid(nrho_1)))/dx_grid
    c10 =  (I1_grid(nx1_2,nx2_1,nrho_1)*(log_x_array_grid(nrho_2)-log_rho) + &
            I1_grid(nx1_2,nx2_1,nrho_2)*(log_rho-log_x_array_grid(nrho_1)))/dx_grid
    c01 =  (I1_grid(nx1_1,nx2_2,nrho_1)*(log_x_array_grid(nrho_2)-log_rho) + &
            I1_grid(nx1_1,nx2_2,nrho_2)*(log_rho-log_x_array_grid(nrho_1)))/dx_grid
    c11 =  (I1_grid(nx1_2,nx2_2,nrho_1)*(log_x_array_grid(nrho_2)-log_rho) + &
            I1_grid(nx1_2,nx2_2,nrho_2)*(log_rho-log_x_array_grid(nrho_1)))/dx_grid
    c0 = (c00*(log_x_array_grid(nx1_2)-log_x1) + &
          c10*(log_x1-log_x_array_grid(nx1_1)))/dx_grid
    c1 = (c01*(log_x_array_grid(nx1_2)-log_x1) + &
          c11*(log_x1-log_x_array_grid(nx1_1)))/dx_grid
    I1_val = (c0*(log_x_array_grid(nx2_2)-log_x2) + &
              c1*(log_x2-log_x_array_grid(nx2_1)))/dx_grid

    c00 =  (I4_grid(nx1_1,nx2_1,nrho_1)*(log_x_array_grid(nrho_2)-log_rho) + &
            I4_grid(nx1_1,nx2_1,nrho_2)*(log_rho-log_x_array_grid(nrho_1)))/dx_grid
    c10 =  (I4_grid(nx1_2,nx2_1,nrho_1)*(log_x_array_grid(nrho_2)-log_rho) + &
            I4_grid(nx1_2,nx2_1,nrho_2)*(log_rho-log_x_array_grid(nrho_1)))/dx_grid
    c01 =  (I4_grid(nx1_1,nx2_2,nrho_1)*(log_x_array_grid(nrho_2)-log_rho) + &
            I4_grid(nx1_1,nx2_2,nrho_2)*(log_rho-log_x_array_grid(nrho_1)))/dx_grid
    c11 =  (I4_grid(nx1_2,nx2_2,nrho_1)*(log_x_array_grid(nrho_2)-log_rho) + &
            I4_grid(nx1_2,nx2_2,nrho_2)*(log_rho-log_x_array_grid(nrho_1)))/dx_grid
    c0 = (c00*(log_x_array_grid(nx1_2)-log_x1) + &
          c10*(log_x1-log_x_array_grid(nx1_1)))/dx_grid
    c1 = (c01*(log_x_array_grid(nx1_2)-log_x1) + &
          c11*(log_x1-log_x_array_grid(nx1_1)))/dx_grid
    I4_val = (c0*(log_x_array_grid(nx2_2)-log_x2) + &
              c1*(log_x2-log_x_array_grid(nx2_1)))/dx_grid

  end subroutine get_grid
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  subroutine get_scalar_evec(ii,i,ktau,evec,devec,warning) 
    implicit none
    integer, intent(in) :: ii,i
    real(DP), intent(in) :: ktau
    real(DP), intent(out) :: evec(2),devec(2)
    logical, intent(in), optional :: warning
    real(DP) :: log_ktau,a,b,dely
    integer :: index,index_p1
    if((ktau.lt.ktau_min).or.(ktau.gt.ktau_max)) then
       if (present(warning)) write(*,*) 'Warning, get_ss_evec out of range, setting to zero with ktau=',ktau
       evec(2) = 0.0d0
       devec(2) = 0.0d0
       return
    end if
    log_ktau = log(ktau)    
    index=int(log(ktau/ktau_min)/dktau+1.0d0)
    index_p1=index+1
    a=(log_ktau_array(index_p1)-log_ktau)/dktau
    b=(log_ktau-log_ktau_array(index))/dktau
    evec(1)=a*ss00_evec(ii,i,index)+b*ss00_evec(ii,i,index_p1)+((a**3-a)*ss00_evec_pr(ii,i,index)+&
         (b**3-b)*ss00_evec_pr(ii,i,index_p1))*(dktau**2)/6.0d0
    evec(2)=a*ss_evec(ii,i,index)+b*ss_evec(ii,i,index_p1)+((a**3-a)*ss_evec_pr(ii,i,index)+&
         (b**3-b)*ss_evec_pr(ii,i,index_p1))*(dktau**2)/6.0d0
    dely=ss00_evec(ii,i,index_p1)-ss00_evec(ii,i,index)
    devec(1)=dely/dktau-(3.d0*a**2-1.0d0)*dktau*ss00_evec_pr(ii,i,index)/6.0d0 + &
         (3.d0*b**2-1.0d0)*dktau*ss00_evec_pr(ii,i,index_p1)/6.0d0 
    dely=ss_evec(ii,i,index_p1)-ss_evec(ii,i,index)
    devec(2)=dely/dktau-(3.d0*a**2-1.0d0)*dktau*ss_evec_pr(ii,i,index)/6.0d0 + &
         (3.d0*b**2-1.0d0)*dktau*ss_evec_pr(ii,i,index_p1)/6.0d0 
  end subroutine get_scalar_evec
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  subroutine get_vv_evec(ii,i,ktau,evec,warning) 
    implicit none
    integer, intent(in) :: ii,i
    real(DP), intent(in) :: ktau
    real(DP), intent(out) :: evec
    logical, intent(in), optional :: warning
    real(DP) :: log_ktau,a,b
    integer :: index,index_p1
    if((ktau.lt.ktau_min).or.(ktau.gt.ktau_max)) then
       if (present(warning)) write(*,*) 'Warning, get_vv_evec out of range, setting to zero with ktau=',ktau
       evec = 0.0d0 
       return
    end if
    log_ktau = log(ktau)
    index=int(log(ktau/ktau_min)/dktau+1.0d0)
    index_p1=index+1
    a=(log_ktau_array(index_p1)-log_ktau)/dktau
    b=(log_ktau-log_ktau_array(index))/dktau
    evec=a*vv_evec(ii,i,index)+b*vv_evec(ii,i,index_p1)+((a**3-a)*vv_evec_pr(ii,i,index)+&
         (b**3-b)*vv_evec_pr(ii,i,index_p1))*(dktau**2)/6.0d0
  end subroutine get_vv_evec
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  subroutine get_tt_evec(ii,i,ktau,evec,warning) 
    implicit none
    integer, intent(in) :: ii,i
    real(DP), intent(in) :: ktau
    real(DP), intent(out) :: evec
    logical, intent(in), optional :: warning
    real(DP) :: log_ktau,a,b
    integer :: index,index_p1
    if((ktau.lt.ktau_min).or.(ktau.gt.ktau_max)) then
       if (present(warning)) write(*,*) 'Warning, get_tt_evec out of range, setting to zero with ktau=',ktau
       evec = 0.0d0 
       return
    end if
    log_ktau = log(ktau)
    index=int(log(ktau/ktau_min)/dktau+1.0d0)
    index_p1=index+1
    a=(log_ktau_array(index_p1)-log_ktau)/dktau
    b=(log_ktau-log_ktau_array(index))/dktau
    evec=a*tt_evec(ii,i,index)+b*tt_evec(ii,i,index_p1)+((a**3-a)*tt_evec_pr(ii,i,index)+&
         (b**3-b)*tt_evec_pr(ii,i,index_p1))*(dktau**2)/6.0d0 
  end subroutine get_tt_evec
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine get_uetc(tau1,tau2,k,uetc_val)
    implicit none
    real(DP), intent(in) :: tau1,tau2,k
    real(DP), intent(out) :: uetc_val(5)
    real(DP) :: x,x1,x2,xp,xm,rho,sum
    real(DP) :: I1,I2,I3,I4,I5,I6
    real(DP) :: a1,a2,a3,a4,a5,a6,C
    integer i,nterms,status
    real(DP) :: xi1,xi2,alpha,mu,v,L
    Type(StringParams) :: SPR1,SPR2

    call get_SPR(tau1,SPR1)
    call get_SPR(tau1,SPR2)

    xi1 = SPR1%xi
    xi2 = SPR2%xi
    alpha = SPR1%alpha
    mu = SPR1%mu
    v=SPR1%v
    L=SPR1%L

    x1 = k*tau1*xi1
    x2 = k*tau2*xi2
    xp=(x1+x2)/2.0d0
    xm=(x1-x2)/2.0d0
    rho = k*v*abs(tau1-tau2)

    if (uetc_feedback.gt.1) then 
       write(*,*) 'x1:',x1
       write(*,*) 'x2:',x2
       write(*,*) 'rho:',rho
    end if

    uetc_val(:) = 0.0d0

    if ((x1.le.xmin) .and. (x2.le.xmin)) then
       ! Use form for small x1 and x2
       if (uetc_feedback.gt.1) write(*,*) 'Using small x'
       uetc_val(1) = mu**2.0d0*x1*x2*alpha**2/(k**2.0d0*(1.0d0-v**2.0d0))
       uetc_val(2) = mu**2.0d0*x1*x2*(1.0d0+v**2*(alpha**2-2.0d0)+v**4*(1.0d0-alpha**2+alpha**4))/&
            (5.0d0*k**2.0d0*(1.0d0-v**2.0d0)*alpha**2)
       uetc_val(3) = mu**2*x1*x2*(1.0d0+v**2*(alpha**2-2.0d0)+v**4*(1.0d0-alpha**2+alpha**4))/&
            (15.0d0*k**2*(1.0d0-v**2)*alpha**2)
       uetc_val(4) = mu**2*x1*x2*(1.0d0+v**2*(alpha**2-2.0d0)+v**4*(1.0d0-alpha**2+alpha**4))/&
            (15.0d0*k**2*(1.0d0-v**2)*alpha**2)
       uetc_val(5) = mu**2.0d0*x1*x2*((x1**2+x2**2)*(2.0d0+v**2*(alpha**2-2.0d0))-12.0d0*rho**2*(1.0d0+&
            v**2*(2.0d0*alpha**2-1.0d0)))/(360*k**2.0d0*(1.0d0-v**2.0d0))
       uetc_val(1:5) = scaling_factor(tau1,tau2,xi1,xi2,L)*uetc_val(1:5)
       return
    end if

    if (abs(x1-x2).le.etcmin) then
       ! Use ETC expression
       if (uetc_feedback.gt.1) write(*,*) 'Using ETC'
       x=(x1+x2)/2.0d0
       uetc_val(1) = 2.0d0*alpha**2*mu**2*(-1.0d0+cos(x)+x*sine_integral(x))/(k**2.0d0*(1.0d0-v**2.0d0))
       uetc_val(2) = (mu**2*((8*(-18 + x**2) + 8*(-2 + alpha**2)*v**2*(-18 + x**2) + &
               v**4*(8*(-18 + x**2) - 8*alpha**2*(-18 + x**2) + alpha**4*(-54 + 11*x**2)))*cos(x) + &
               (-32*(1 + (-2 + alpha**2)*v**2 + (1 - alpha**2 + alpha**4)*v**4)*x**3 + &
               3*(-8*(-6 + x**2) - 8*(-2 + alpha**2)*v**2*(-6 + x**2) + &
               v**4*(-8*(-6 + x**2) + 8*alpha**2*(-6 + x**2) + alpha**4*(18 + x**2)))*sin(x))/x + &
               (8 + 8*(-2 + alpha**2)*v**2 + (8 - 8*alpha**2 + 11*alpha**4)*v**4)*x**3*sine_integral(x)))/ &
               (16.*alpha**2*k**2*(1 - v**2)*x**2)
       uetc_val(3) = (mu**2*((2.0d0*(8.0d0+8.0d0*(-2.0d0+alpha**2)*v**2 + &
            (8.0d0-8.0d0*alpha**2+3*alpha**4)*v**4)*(x**3 + 3.0d0*x*cos(x) - 3.0d0*sin(x)))/(3.0d0*x**3) &
            + alpha**4*v**4*(-2.0d0+cos(x)+sin(x)/(x)+x*sine_integral(x))))/(8.0d0*alpha**2*k**2*(1-v**2))
       uetc_val(4) = (mu**2*(3*(8 + 8*(-2 + alpha**2)*v**2 + (8 - 8*alpha**2 + 3*alpha**4)*v**4)*(-2 + x**2)*cos(x) + &
            (64*(-1 + v**2)*(1 + (-1 + alpha**2)*v**2)*x**3 - &
            3*(-8*(2 + x**2) - 8*(-2 + alpha**2)*v**2*(2 + x**2) + &
            v**4*(-8*(2 + x**2) + 8*alpha**2*(2 + x**2) + alpha**4*(-6 + 5*x**2)))*sin(x))/x + &
            3*(8 + 8*(-2 + alpha**2)*v**2 + (8 - 8*alpha**2 + 3*alpha**4)*v**4)*x**3*sine_integral(x)))/ &
            (96.*alpha**2*k**2*(1 - v**2)*x**2) 
       uetc_val(5) = (mu**2*(2 + (-2 + alpha**2)*v**2)*(-4 + cos(x) + (3*sin(x))/x + x*sine_integral(x)))/(2.*k**2*(1 - v**2))
       uetc_val(1:5) = scaling_factor(tau1,tau2,xi1,xi2,L)*uetc_val(1:5)
       return
    end if

    if (abs(x1-x2).ge.xapr) then
       !Use Ed's approximation
       if (uetc_feedback.gt.1) write(*,*) 'Using Eds Approximation'
       xp=(x1+x2)/2.0d0
       xm=(x1-x2)/2.0d0
       I1 = I1_int_a(min(x1,x2),rho)
       I4 = I4_int_a(min(x1,x2),rho)

    else

    nterms = max(min_terms,int(scale_terms*xp))
    if(uetc_feedback.gt.1) write(*,*) 'N terms:',nterms

    ! Combined integral
    if(use_grid) then
       call get_grid(x1,x2,rho,I1,I4,status)
       if (status.eq.1) then
          I1 = I1_int_num(x1,x2,rho)
          I4 = I4_int_num(x1,x2,rho)
       end if
    else
       
       if (do_num.or.(x1.gt.xmax).or.(x2.gt.xmax)) then 
          if (uetc_feedback.gt.1) write(*,*) 'Exceeds xmax, integrating I1 and I4'
          I1 = I1_int_num(x1,x2,rho)
          I4 = I4_int_num(x1,x2,rho)
       else
          I1 = I1_int(xm,rho,nterms)-I1_int(xp,rho,nterms)
          I4 = I4_int(xm,rho,nterms)-I4_int(xp,rho,nterms)
       end if
    end if
    end if
    I2 = I2_int(xm,rho)-I2_int(xp,rho)
    I3 = I3_int(xm,rho)-I3_int(xp,rho)
    I5 = I5_int(xm,rho)-I5_int(xp,rho)
    I6 = I6_int(xm,rho)-I6_int(xp,rho)

    if (uetc_feedback.gt.1) then
       write(*,*) 'Components:'
       write(*,'(2E15.5)') I1
       write(*,'(2E15.5)') I2
       write(*,'(2E15.5)') I3
       write(*,'(2E15.5)') I4
       write(*,'(2E15.5)') I5
       write(*,'(2E15.5)') I6
    endif

    a1 = 2.0d0*alpha**2
    sum = a1*I1
    uetc_val(1) = sum*mu**2*scaling_factor(tau1,tau2,xi1,xi2,L)/(k**2*(1.0d0-v**2))
    
    a1 = (-27*alpha**4*v**4 + rho**2*(1 + (-1 + 2*alpha**2)*v**2)**2)/(2.*alpha**2*rho**2)
    a2 = (3*(9*alpha**4*v**4 + rho**2*(1 - 2*v**2 - (-1 + alpha**4)*v**4)))/(2.*alpha**2*rho**2)
    a3 = (-9*(1 + (-1 + alpha**2)*v**2)**2)/(2.*alpha**2)
    a4 = 3*v**2*(-1 + (1 + alpha**2*(-2 + 9/rho**2))*v**2)
    a5 = 3*v**2 + (-3 + alpha**2*(6 - 27/rho**2))*v**4
    a6 = 9*(v**2 + (-1 + alpha**2)*v**4)
    sum = a1*I1 + a2*I2 + a3*I3 + a4*I4 + a5*I5 + a6*I6
    uetc_val(2) = sum*mu**2*scaling_factor(tau1,tau2,xi1,xi2,L)/(k**2*(1.0d0-v**2))

    C = 1.0d0+v**2*(alpha**2-1.0d0)
    a1 = 3.0d0*v**4*alpha**2/rho**2
    a2 = -3.0d0*v**4*alpha**2/rho**2
    a3 = C**2/alpha**2
    a4 = -v**4*alpha**2*(6.0d0/rho**2-1.0d0)
    a5 = v**4*alpha**2*(6.0d0/rho**2-1.0d0)
    a6 = -2.0d0*v**2*C
    sum = a1*I1 + a2*I2 + a3*I3 + a4*I4 + a5*I5 + a6*I6
    uetc_val(3) = sum*mu**2*scaling_factor(tau1,tau2,xi1,xi2,L)/(k**2*(1.0d0-v**2))

    a1 = (-3*alpha**4*v**4 + rho**2*(-1 + v**2)**2)/(4.*alpha**2*rho**2)
    a2 = (3*alpha**4*v**4 + rho**2*(-1 + 2*v**2 + (-1 + alpha**4)*v**4))/(4.*alpha**2*rho**2)
    a3 = -(1 + (-1 + alpha**2)*v**2)**2/(4.*alpha**2)
    a4 = (v**2*(1 + (-1 + (3*alpha**2)/rho**2)*v**2))/2.
    a5 = (-3*alpha**2*v**4 + rho**2*v**2*(-1 + v**2))/(2.*rho**2)
    a6 = (v**2 + (-1 + alpha**2)*v**4)/2.
    sum = a1*I1 + a2*I2 + a3*I3 + a4*I4 + a5*I5 + a6*I6
    uetc_val(4) = sum*mu**2*scaling_factor(tau1,tau2,xi1,xi2,L)/(k**2*(1.0d0-v**2))

    a1 = 1.0d0 + (-1.0d0 + 2*alpha**2)*v**2
    a2 = -3.0d0 - 3.0d0*(-1.0d0 + alpha**2)*v**2
    a3 = 0.0d0
    a4 = -3*alpha**2*v**2
    a5 = 3*alpha**2*v**2
    a6 = 0.0d0
    sum = a1*I1 + a2*I2 + a3*I3 + a4*I4 + a5*I5 + a6*I6
    uetc_val(5) = sum*mu**2*scaling_factor(tau1,tau2,xi1,xi2,L)/(k**2*(1.0d0-v**2))

  end subroutine get_uetc 
 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function I1_int(x,rho,n)
    implicit none
    real(DP) I1_int
    real(DP) x,rho,term
    integer i,n
    I1_int = 0.0d0
    do i=1,n
       term = 1.0d0/factorial(i)*(rho/(2.0d0*i-1.0d0))*(-x**2/(2.0d0*rho))**i*spher_bessel(i-1,rho)
       I1_int = I1_int + term
    end do
  end function I1_int
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 function I1_int_a(x,rho)
    implicit none
    real(DP) I1_int_a
    real(DP) x,rho,J0_rho
    J0_rho = BESSEL_JN(0,rho)
    I1_int_a = (pi*x/2.0d0)*J0_rho
 end function I1_int_a
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 function I2_int(x,rho)
    implicit none
    real(DP) I2_int
    real(DP) x,rho,px,rpx,srpx
    px = rho**2+x**2
    rpx = sqrt(px)
    srpx = sin(rpx)
    I2_int = srpx/rpx
  end function I2_int
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 function I3_int(x,rho)
    implicit none
    real(DP) I3_int
    real(DP) x,rho,px,rpx,srpx,crpx
    px = rho**2+x**2
    rpx = sqrt(px)
    srpx = sin(rpx)
    crpx = cos(rpx)
    I3_int = crpx/px*(1.0d0-3.0d0*x**2/px)+srpx/rpx*(1.0d0-(1.0d0+x**2)/px+3.0d0*x**2/px**2)
  end function I3_int
 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function I4_int(x,rho,n)
    implicit none
    real(DP) I4_int
    real(DP) x,rho,term
    integer i,n
    I4_int = cos(x)/rho**2
    do i=1,n 
       term = - 1.0d0/factorial(i)*(1.0d0/(2.0d0*i-1.0d0))*(-x**2/(2.0d0*rho))**i*spher_bessel(i-2,rho)
       I4_int = I4_int + term
    end do
  end function I4_int
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 function I4_int_a(x,rho)
   implicit none
   real(DP) I4_int_a
   real(DP) x,rho,J1_rho
   J1_rho = BESSEL_JN(1,rho)
   I4_int_a = (pi*x*J1_rho)/(2.0d0*rho)
 end function I4_int_a
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 function I5_int(x,rho)
    implicit none
    real(DP) I5_int
    real(DP) x,rho,px,rpx,srpx,crpx
    px = rho**2+x**2
    rpx = sqrt(px)
    crpx = cos(rpx)
    I5_int = (cos(x)-crpx)/rho**2
end function I5_int
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 function I6_int(x,rho)
    implicit none
    real(DP) I6_int
    real(DP) x,rho,px,rpx,srpx,crpx
    px = rho**2+x**2
    rpx = sqrt(px)
    srpx = sin(rpx)
    crpx = cos(rpx)
    I6_int = (srpx/rpx-crpx)/px
end function I6_int
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function scaling_factor(tau1,tau2,xi1,xi2,L)
    implicit none
    real(DP) scaling_factor
    real(DP), intent(in) :: tau1,tau2,xi1,xi2,L
    select case (scaling_option)
    case(1)
       scaling_factor = 1.0d0
    case(2)
       scaling_factor = 1.0d0/(max(xi1*tau1,xi2*tau2))**3
    case(3)
       stop 'Not implemented yet'
    case default
       write(*,*) 'Setting scaling factor to 1'
       scaling_factor = 1.0d0
    end select
  end function scaling_factor
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  function spher_bessel(n,x)
    implicit none
    real(DP) spher_bessel
    real(DP) x
    integer n
    real(DP), parameter :: pi = 3.1415926535897932384626433832795d0
    COMPLEX (KIND=nag_wp)           :: z,cy(1)                      
    REAL (KIND=nag_wp)              :: fnu                     
    INTEGER                         :: nz,ifail      
    CHARACTER (1)                   :: scal     
    if (n.eq.-1) then
        spher_bessel=cos(x)/x
       return
    end if
    fnu=0.5d0+n
    z = dcmplx(x,0.0d0)
    scal='u'
    ifail=0
    CALL s17def(fnu,z,1,scal,cy,nz,ifail)
    if (ifail.ne.0) stop
    spher_bessel=sqrt(pi/(2.0d0*x))*dreal(cy(1))
  end function spher_bessel
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function sine_integral(x)
    implicit none
    real(DP) sine_integral
    real(DP) x
    integer ifail
    ifail=0
    sine_integral = s13adf(x,ifail)
    !call cisi(x,ci,si)
    !sine_integral = si
    if (ifail.ne.0) stop
  end function sine_integral
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function factorial(n)
    implicit none
    real(DP) factorial,x
    integer n,ifail
    x=1.0d0+n
    factorial = s14aaf(x,ifail)
    if (ifail.ne.0) stop
  end function factorial
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine Matrix_Diagonalize(M,diag,n)
  ! Compute eigenvectors and eigenvalues
  ! Eigenvectors are outputted as (M,:,1), M(:,2) ...
  ! corresponding to eigenvalues diag(1), diag(2), ...
  implicit none
  integer, intent(in) :: n
  real(DP), intent(inout):: m(n,n)
  real(DP), intent(out) :: diag(n)
  integer ierr, tmpsize
  real(DP), allocatable, dimension(:) :: tmp
  tmpsize =  max( (ILAENV_wrap(1,'DSYTRD','U',n,n,n,n)+2)*N,max(1,3*n-1))  !3*n**2
  allocate(tmp(tmpsize));
  call DSYEV('V','U',n,m,n,diag,tmp,tmpsize,ierr) !evalues and vectors of symmetric matrix
  if (ierr /= 0) stop 'Error in Matrix_Diagonalize'
  deallocate(tmp)
end subroutine Matrix_Diagonalize
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function ILAENV_wrap(i,S1,S2,a,b,c,d)
  implicit none
  integer ILAENV_wrap
  integer, intent(in) :: i,a,b,c,d
  character(LEN=*), intent(in) :: S1, S2
  integer, external :: ILAENV
  ILAENV_wrap =  ILAENV(i,S1,S2,a,b,c,d)
end  function ILAENV_wrap
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function I1_int_num(x1,x2,rho)
  implicit none
  real(DP) I1_int_num,x1,x2,rho
  real(DP) I1_integrand,rombint_mod
  external I1_integrand,rombint_mod
  real(DP) :: pio2 = 1.5707963267948966192d0
  real(DP) :: eps = 1.0d-8
  I1_int_num = 2.0d0*rombint_mod(x1,x2,rho,I1_integrand,0.0d0+eps,pio2-eps,1.0d-8)
end function I1_int_num
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function I4_int_num(x1,x2,rho)
  implicit none
  real(DP) I4_int_num,x1,x2,rho
  real(DP) I4_integrand,rombint_mod
  external I4_integrand,rombint_mod
  real(DP) :: pio2 = 1.5707963267948966192d0
  real(DP) :: eps = 1.0d-8 
  I4_int_num = 2.0d0*rombint_mod(x1,x2,rho,I4_integrand,0.0d0+eps,pio2-eps,1.0d-8)
end function I4_int_num
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine delete_file(filename,status)
!
!*******************************************************************************
!
!! DELETE_FILE deletes a FITS file.
!
  integer status,unit,blocksize
  character (len =*) filename
!
!  Simply return if status is greater than zero.
!
  if (status > 0) then
    return
  end if
!
!  Get an unused Logical Unit Number to use to open the FITS file
!
  call ftgiou ( unit, status )
!
!  Try to open the file, to see if it exists
!
  call ftopen ( unit, filename, 1, blocksize, status )
  if ( status == 0 ) then
!
!  File was opened;  so now delete it 
!
    call ftdelt(unit,status)
  else if (status == 103)then
!
!  File doesn't exist, so just reset status to zero and clear errors
!
    status=0
    call ftcmsg
  else
!
!  There was some other error opening the file; delete the file anyway
!
      status=0
      call ftcmsg
      call ftdelt(unit,status)
  end if
!
!  Free the unit number for later reuse.
!
  call ftfiou(unit, status)
  return
end subroutine delete_file
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end module uetc

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function I1_integrand(x1,x2,rho,theta)
  use Precision
  implicit none
  real(dl) I1_integrand,theta,x1,x2,rho
  I1_integrand = sin(theta)/cos(theta)**2*sin(x1/2.0d0*cos(theta))*sin(x2/2.0d0*cos(theta))* &
       BESSEL_J0(rho*sin(theta))
end function I1_integrand
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function I4_integrand(x1,x2,rho,theta)
  use Precision
  implicit none
  real(dl) I4_integrand,theta,x1,x2,rho
  I4_integrand = sin(theta)/cos(theta)**2*sin(x1/2.0d0*cos(theta))*sin(x2/2.0d0*cos(theta))* &
       BESSEL_J1(rho*sin(theta))/(rho*sin(theta))
end function I4_integrand
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 function rombint_mod(x1,x2,x3,f,a,b,tol, maxit)
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
        real(dl) x1,x2,x3 !dummy
        real(dl) f
        external f
        real(dl) :: rombint_mod
        real(dl), intent(in) :: a,b,tol
        integer :: nint, i, k, jmax, j
        real(dl) :: h, gmax, error, g, g0, g1, fourj
!

        if (present(maxit)) then
            MaxIter = maxit
        end if
        h=0.5d0*(b-a)
        gmax=h*(f(x1,x2,x3,a)+f(x1,x2,x3,b))
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
            g0=g0+f(x1,x2,x3,a+(k+k-1)*h)
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
40      rombint_mod=g0
        if (i.gt.MAXITER.and.abs(error).gt.tol)  then
          write(*,*) 'Warning: Rombint failed to converge; '
          write (*,*)'integral, error, tol:', rombint_mod,error, tol
        end if
        
        end function rombint_mod
