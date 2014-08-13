 !Interface module for CAMB. Call CAMB_GetResults to do the work.

     module CAMB
         use Precision
         use ModelParams
         use ModelData
         use Transfer
         use GaugeInterface
         use InitialPower
         use Reionization
         use Recombination
         use lensing
         use uetc
         implicit none
         
         Type CAMBdata
            Type (ClTransferData) :: ClTransScal,ClTransTens,ClTransVec
            Type (MatterTransferData) :: MTrans
            Type (CAMBparams) :: Params
         end Type CAMBdata

!         public CAMB_GetTransfers, CAMB_GetResults, CAMB_GetCls, CAMB_SetDefParams, & 
!                CAMB_ValidateParams, CAMB_GetAge,CAMB_InitCAMBdata, 
     contains

       subroutine CAMB_GetTransfers(Params, OutData, error)
        use CAMBmain
        use lensing
        type(CAMBparams) :: Params
        type (CAMBdata)  :: OutData
        integer :: error !Zero if OK
  
      !Set internal types from OutData so it always 'owns' the memory, prevent leaks
      
        MT =  OutData%MTrans

        CTransScal = OutData%ClTransScal  
        CTransVec  = OutData%ClTransVec  
        CTransTens = OutData%ClTransTens
  

        call CAMB_GetResults(Params, error)
        OutData%Params = Params
        OutData%MTrans = MT
        OutData%ClTransScal = CTransScal
        OutData%ClTransVec  = CTransVec
        OutData%ClTransTens = CTransTens
  
       end subroutine CAMB_GetTransfers

       
       ! AJM
       subroutine CAMB_GetTransfers_string(Params,OutData,plot_uetc,error)
         use CAMBmain
         use lensing
         implicit none
         type (CAMBparams) :: Params
         type (CAMBdata) :: OutData
         integer :: error 
         logical plot_uetc
         MT = OutData%MTrans
         CTransScal = OutData%ClTransScal
         CTransVec = OutData%ClTransVec
         CTransTens = OutData%ClTransTens
         call CAMB_GetResults_string(Params,plot_uetc,error)
         OutData%Params = Params
         OutData%MTrans = MT
         OutData%ClTransScal = CTransScal
         OutData%ClTransVec = CTransVec
         OutData%ClTransTens = CTransTens
       end subroutine CAMB_GetTransfers_string

       subroutine CAMB_InitCAMBdata(Dat)
        type (CAMBdata) :: Dat
 
!Comment these out to try to avoid intel bugs with status deallocating uninitialized pointers   
        call Ranges_Nullify(Dat%ClTransScal%q)
        call Ranges_Nullify(Dat%ClTransVec%q)
        call Ranges_Nullify(Dat%ClTransTens%q)

        nullify(Dat%ClTransScal%Delta_p_l_k)
        nullify(Dat%ClTransVec%Delta_p_l_k)
        nullify(Dat%ClTransTens%Delta_p_l_k)
        nullify(Dat%MTrans%sigma_8,Dat%MTrans%TransferData,Dat%MTrans%q_trans)
             
       end subroutine CAMB_InitCAMBdata


       subroutine CAMB_FreeCAMBdata(Dat)
            type (CAMBdata) :: Dat

            call Free_ClTransfer(Dat%ClTransScal)
            call Free_ClTransfer(Dat%ClTransVec)
            call Free_ClTransfer(Dat%ClTransTens)
            call Transfer_Free(Dat%MTrans)

       end subroutine CAMB_FreeCAMBdata


       subroutine CAMB_TransfersToPowers(CData)
        use CAMBmain
        use lensing
        type (CAMBdata) :: CData

        CP = CData%Params
        call InitializePowers(CP%InitPower,CP%curv)
        if (global_error_flag/=0) return
        if (CData%Params%WantCls) then
          call ClTransferToCl(CData%ClTransScal,CData%ClTransTens, CData%ClTransvec)  
          if (CP%DoLensing .and. global_error_flag==0) call lens_Cls
          if (global_error_flag/=0) return
          if (CP%OutputNormalization == outCOBE) call COBEnormalize
        end if
        if (CData%Params%WantTransfer) call Transfer_Get_sigma8(Cdata%MTrans,8._dl)
     
       end subroutine CAMB_TransfersToPowers
  


       !Call this routine with a set of parameters to generate the results you want.
       subroutine CAMB_GetResults(Params, error)
        use CAMBmain
        use lensing
        use Bispectrum
        use Errors
        type(CAMBparams) :: Params
        integer, optional :: error !Zero if OK
        type(CAMBparams) P
        logical :: separate = .true. !whether to do P_k in separate call or not
        logical :: InReionization
        
        if (Params%DoLensing .and. Params%NonLinear==NonLinear_Lens) separate = .false.
        InReionization = Params%Reion%Reionization
        global_error_flag = 0
        
         if (Params%WantCls .and. Params%WantScalars) then
          P = Params
          if (separate) then
          P%WantTransfer = .false.
          P%Transfer%high_precision = .false.
          end if
          P%WantTensors = .false.
          P%WantVectors = .false.
          call CAMBParams_Set(P) 
          if (global_error_flag==0) call cmbmain
          if (global_error_flag/=0) then
            if (present(error)) error =global_error_flag
            return 
          end if
          call_again = .true.
          !Need to store CP%flat etc, but must keep original P_k settings
          CP%Transfer%high_precision = Params%Transfer%high_precision
          CP%WantTransfer = Params%WantTransfer
          CP%WantTensors = Params%WantTensors
          CP%WantVectors = Params%WantVectors
          CP%Transfer%num_redshifts = Params%Transfer%num_redshifts
          Params = CP            
         end if
 
         if (Params%WantCls .and. Params%WantTensors) then
          P=Params
          P%WantTransfer = .false.
          P%Transfer%high_precision = .false.
          P%WantScalars = .false.
          P%WantVectors = .false.
          call CAMBParams_Set(P)  
          if (global_error_flag==0) call cmbmain
          if (global_error_flag/=0) then
            if (present(error)) error =global_error_flag
            return 
          end if
           call_again = .true.
           CP%Transfer%high_precision = Params%Transfer%high_precision
           CP%WantTransfer = Params%WantTransfer
           CP%WantScalars = Params%WantScalars
           CP%WantVectors = Params%WantVectors
           CP%Transfer%num_redshifts = Params%Transfer%num_redshifts
           Params = CP            
   
         end if
 
         if (Params%WantCls .and. Params%WantVectors) then
          P=Params
          P%WantTransfer = .false.
          P%Transfer%high_precision = .false.
          P%WantScalars = .false.
          P%WantTensors = .false.
          call CAMBParams_Set(P)  
          if (global_error_flag==0) call cmbmain
          if (global_error_flag/=0) then
            if (present(error)) error =global_error_flag
            return 
          end if
          call_again = .true.
          CP%Transfer%high_precision = Params%Transfer%high_precision
          CP%WantTransfer = Params%WantTransfer
          CP%WantTensors = Params%WantTensors
          CP%WantScalars = Params%WantScalars
          CP%Transfer%num_redshifts = Params%Transfer%num_redshifts
          Params = CP            
   
         end if

         if (Params%WantTransfer .and. &
          .not. (Params%WantCls .and. Params%WantScalars .and. .not. separate)) then
          P=Params
          P%WantCls = .false.
          P%WantScalars = .false.
          P%WantTensors = .false.
          P%WantVectors = .false.
          call CAMBParams_Set(P)  
          if (global_error_flag==0) call cmbmain
          if (global_error_flag/=0) then
            if (present(error)) error =global_error_flag
            return 
          end if
          !Need to store num redshifts etc
          CP%WantScalars = Params%WantScalars
          CP%WantCls =  Params%WantCls
          CP%WantTensors = Params%WantTensors
          CP%WantVectors = Params%WantVectors
          CP%Reion%Reionization = InReionization
          Params = CP            

         end if

        call_again = .false.
   

        if (.not. CP%OnlyTransfers) then

         if (CP%WantCls .and. CP%OutputNormalization == outCOBE) call COBEnormalize

         if (CP%DoLensing .and. global_error_flag==0) then
           call lens_Cls 
         end if
         
         if (do_bispectrum .and. global_error_flag==0) call GetBispectrum(CTransScal) 

        end if

        end subroutine CAMB_GetResults


        !Return real (NOT double precision) arrays of the computed CMB  Cls
        !Output is l(l+1)C_l/2pi
        !If GC_Conventions = .false. use E-B conventions (as the rest of CAMB does)
        subroutine CAMB_GetCls(Cls, lmax, in, GC_conventions)
          integer, intent(IN) :: lmax, in
          logical, intent(IN) :: GC_conventions
          real, intent(OUT) :: Cls(2:lmax,1:4)
          integer l
      
          Cls = 0
          do l=2, lmax
             if (CP%WantScalars .and. l<= CP%Max_l) then
             Cls(l,1:2) = Cl_scalar(l, in,  C_Temp:C_E)
             Cls(l,4) = Cl_scalar(l, in,  C_Cross)
             end if
             if (CP%WantTensors .and. l <= CP%Max_l_tensor) then
                Cls(l,1:4) = Cls(l,1:4) + Cl_tensor(l, in,  CT_Temp:CT_Cross)
             end if
          end do
          if (GC_conventions) then
             Cls(:,2:3) = Cls(:,2:3)/2
             Cls(:,4)   = Cls(:,4)/sqrt(2.0)             
          end if
 
        end subroutine CAMB_GetCls

        function CAMB_GetAge(P)
           !Return age in gigayears, returns -1 on error
           use constants
           type(CAMBparams), intent(in) :: P
           real(dl) CAMB_GetAge
           real(dl) atol,a1,a2, dtda, rombint
!           real(dl), parameter :: Mpc = 3.085678e22_dl, &
!                 c = 2.99792458e8_dl, Gyr=3.1556926e16
           integer error
           external dtda,rombint


           call  CAMBParams_Set(P, error, .false.)

           if (error/=0) then
            CAMB_GetAge = -1
           else

           atol = 1d-4
           a1=0
           a2=1
           CAMB_GetAge = rombint(dtda,a1,a2,atol)*Mpc/c/Gyr
           end if
    
         end function CAMB_GetAge


        function CAMB_GetZreFromTau(P, tau)
           type(CAMBparams) :: P
           real(dl) tau
           real(dl) CAMB_GetZreFromTau
           integer error

            P%Reion%use_optical_depth = .true.
            P%Reion%optical_depth = tau
            call CAMBParams_Set(P,error)
            
            CAMB_GetZreFromTau = CP%Reion%redshift

        end function CAMB_GetZreFromTau

      
        subroutine CAMB_SetDefParams(P)
            use Bispectrum
            type(CAMBparams), intent(out) :: P

            P%WantTransfer= .false.
            P%WantCls = .true.

            P%omegab  = .045
            P%omegac  = 0.255
            P%omegav  = 0.7
            P%omegan  = 0
            P%H0      = 65

            P%TCMB    = 2.726
            P%YHe     = 0.24
            P%Num_Nu_massless =3.04
            P%Num_Nu_massive  =0
            P%Nu_mass_splittings = .false.
            P%Nu_mass_eigenstates = 0
           
            P%Scalar_initial_condition =initial_adiabatic
            P%NonLinear = NonLinear_none
            
            call SetDefPowerParams(P%InitPower)

            call Recombination_SetDefParams(P%Recomb)
          
            call Reionization_SetDefParams(P%Reion)
            
            P%Transfer%high_precision=.false.
    
            P%OutputNormalization = outNone

            P%WantScalars = .true.
            P%WantVectors = .false.
            P%WantTensors = .false.
            P%want_zstar = .true.  !!JH 
            P%want_zdrag = .true.  !!JH  

            P%Max_l=1500
            P%Max_eta_k=3000
            P%Max_l_tensor=400
            P%Max_eta_k_tensor=800
            !Set up transfer just enough to get sigma_8 OK
            P%Transfer%kmax=0.9  
            P%Transfer%k_per_logint=0
            P%Transfer%num_redshifts=1
            P%Transfer%redshifts=0

            P%AccuratePolarization = .true.
            P%AccurateReionization = .false.
            P%AccurateBB = .false.

            P%DoLensing = .false.

            P%MassiveNuMethod = Nu_best
            P%OnlyTransfers = .false.

         end subroutine CAMB_SetDefParams


         !Stop with error is not good
         function CAMB_ValidateParams(P) result(OK)
            type(CAMBparams), intent(in) :: P
            logical OK

             OK = .true.
             if (.not. P%WantTransfer .and. .not. P%WantCls) then
                OK = .false.
                write(*,*) 'There is nothing to do! Do transfer functions or Cls.'
             end if

             if (P%h0 < 20._dl.or.P%h0 > 100._dl) then
               OK = .false.
               write(*,*) '  Warning: H0 has units of km/s/Mpc. You have:', P%h0
            end if
             if (P%tcmb < 2.7d0.or.P%tcmb > 2.8d0) then
                write(*,*) '  Warning: Tcmb has units of K.  Your have:', P%tcmb
             end if

             if (P%yhe < 0.2d0.or.P%yhe > 0.8d0) then
                OK = .false.
                write(*,*) &
                     '  Warning: YHe is the Helium fraction of baryons.', &
                     '  Your have:', P%yhe
             end if
             if (P%Num_Nu_massive < 0.or.P%Num_Nu_massive > 3.1) then
                OK = .false.
                write(*,*) &
                     'Warning: Num_Nu_massive is strange:',P%Num_Nu_massive 
              end if
             if (P%Num_Nu_massless < 0.or.P%Num_Nu_massless > 3.1) then
                OK = .false.
                write(*,*) &
                     'Warning: Num_nu_massless is strange:', P%Num_Nu_massless
              end if
             if (P%Num_Nu_massive < 1 .and. P%omegan > 0.0) then
                OK = .false.
                write(*,*) &
                     'Warning: You have omega_neutrino > 0, but no massive species'
              end if


             if (P%omegab<0.001 .or. P%omegac<0 .or. P%omegab>1 .or. P%omegac>3) then
                OK = .false.
                write(*,*) 'Your matter densities are strange'
             end if

             if (P%WantScalars .and. P%Max_eta_k < P%Max_l .or.  &
                  P%WantTensors .and. P%Max_eta_k_tensor < P%Max_l_tensor) then
                OK = .false.
                write(*,*) 'You need Max_eta_k larger than Max_l to get good results'
             end if
             
             call Reionization_Validate(P%Reion, OK)
             call Recombination_Validate(P%Recomb, OK)

             if (P%WantTransfer) then
              if (P%transfer%num_redshifts > max_transfer_redshifts .or. P%transfer%num_redshifts<1) then
                OK = .false.
                write(*,*) 'Maximum ',  max_transfer_redshifts, &
                     'redshifts. You have: ', P%transfer%num_redshifts 
              end if
              if (P%transfer%kmax < 0.01 .or. P%transfer%kmax > 50000 .or. &
                     P%transfer%k_per_logint>0 .and.  P%transfer%k_per_logint <1) then
                 OK = .false.
                 write(*,*) 'Strange transfer function settings.'
              end if
              if (P%transfer%num_redshifts > max_transfer_redshifts .or. P%transfer%num_redshifts<1) then
                OK = .false.
                write(*,*) 'Maximum ',  max_transfer_redshifts, &
                     'redshifts. You have: ', P%transfer%num_redshifts 
              end if


             end if

         end function CAMB_ValidateParams

         subroutine CAMB_cleanup
          use ThermoData
          use SpherBessels
          use ModelData
          use Transfer

            !Free memory
           call ThermoData_Free
           call Bessels_Free
           call ModelData_Free  
           call Transfer_Free(MT)

         end subroutine CAMB_cleanup

!AJM ------------------------------------------------------------------------------

  subroutine CAMB_GetResults_string(Params, plot_uetc,error)
        use CAMBmain
        use lensing
        use Bispectrum
        use Errors
        use timer
        type(CAMBparams) :: Params
        logical, intent(in) :: plot_uetc
        integer, optional :: error !Zero if OK
        type(CAMBparams) P
        logical :: separate = .true. !whether to do P_k in separate call or not
        logical :: InReionization
        integer :: last_c,im,il
        real(DP) :: scale_factor,vec_factor,tens_factor
        CHARACTER(LEN=8) :: ELAPSED_TIME
        real(DP) evec
        real(DP) :: scalar_power(5000),vector_power(5000),tensor_power(5000)
        logical :: output_power = .false.

        CALL TIMER_START

        if (Params%DoLensing .and. Params%NonLinear==NonLinear_Lens) separate = .false.
        InReionization = Params%Reion%Reionization
        global_error_flag = 0

        Params%Max_l_tensor = Params%Max_l
        Params%Max_eta_k_tensor = Params%Max_eta_k 

        ! Set to white noise and unit amplitude
        Params%InitPower%an(1) = 4.0
        Params%InitPower%ant(1) = 3.0
        Params%InitPower%rat(1) = 1.0
        Params%InitPower%ScalarPowerAmp(1) = 1.0
        Params%InitPower%k_0_scalar = 1.0
        Params%InitPower%k_0_tensor = 1.0
        Params%Scalar_initial_condition = 6
        nmodes = Params%nmodes

        !call CAMBParams_Set(P) 
        !call InitVars
        !call SetkValuesForSources 
        do_string_source = .true.

!        call init_string(SPR_rad)
!        call init_string(SPR_mat)
!        ! Multiply final spectrum by mu^2
!        SPR_rad%mu = 1.0! Params%mu
!        SPR_rad%xi = Params%xi_rad
!        SPR_rad%alpha = Params%alpha_rad
!        SPR_rad%v = Params%v_rad 
!        SPR_rad%L = Params%L
!        SPR_mat%mu = 1.0! Params%mu
!        SPR_mat%xi = Params%xi_mat
!        SPR_mat%alpha = Params%alpha_mat
!        SPR_mat%v = Params%v_mat 
!        SPR_mat%L = Params%L
        ktau_min = Params%uetc_kt_min
        ktau_max = Params%uetc_kt_max
        nktau = Params%uetc_n
        output_uetc = plot_uetc
!        if ((SPR_rad%mu.eq.SPR_mat%mu).and.(SPR_rad%xi.eq.SPR_mat%xi).and.(SPR_rad%v.eq.SPR_mat%v)&
!             .and.(SPR_rad%L.eq.SPR_mat%L).and.(SPR_rad%alpha.eq.SPR_mat%alpha)) then
!           no_evolve=.true.
!        else
!           no_evolve=.false.
!        end if

        ! Note on normalization 
        ! C_l defined as C_l=4pi \int dk/k \Delta_zeta^2 Delta_l(k)^2 where
        ! Delta_zeta^2=k^3P_zeta(k)/(2pi^2) is power per log interval of curvature 
        ! perturbation
        scale_factor = 1.0/(2.0*pi**2)
        ! AJM changed vector factor
        vec_factor = 8.0
        tens_factor = 16.0

        scalar_power(:) = 0.0d0
        vector_power(:) = 0.0d0
        tensor_power(:) = 0.0d0

        if(output_power) open(unit=80,file=trim(uetc_file)//'_power.dat')

        do_diagonalize = .true.

        do im=1,nmodes
           
           write(*,*) 'Mode:',im   

           P = Params
           P%WantTransfer = .false.
           P%Transfer%high_precision = .false.
           P%WantScalars = .true.
           P%WantTensors = .false.
           P%WantVectors = .false.
           !write(*,*) 'Doing scalars'
           call CAMBParams_Set(P) 
           if (global_error_flag==0) call cmbmain(im)
           call_again = .true.
           CP%Transfer%high_precision = Params%Transfer%high_precision
           CP%WantTransfer = Params%WantTransfer
           CP%WantScalars = Params%WantScalars
           CP%WantVectors = Params%WantVectors
           CP%Transfer%num_redshifts = Params%Transfer%num_redshifts
           Params = CP      

           last_C=min(C_PhiTemp,C_last)
        
           P = Params
           P%WantTransfer = .false.
           P%Transfer%high_precision = .false.
           P%WantScalars = .false.
           P%WantTensors = .true.
           P%WantVectors = .false.
           !write(*,*) 'Doing tensors'
           call CAMBParams_Set(P)
           if (global_error_flag==0) call cmbmain(im)
           call_again = .true.
           CP%Transfer%high_precision = Params%Transfer%high_precision
           CP%WantTransfer = Params%WantTransfer
           CP%WantScalars = Params%WantScalars
           CP%WantVectors = Params%WantVectors
           CP%Transfer%num_redshifts = Params%Transfer%num_redshifts
           Params = CP        
           
           P = Params
           P%WantTransfer = .false.
           P%Transfer%high_precision = .false.
           P%WantScalars = .false.
           P%WantTensors = .false.
           P%WantVectors = .true.
           !write(*,*) 'Doing vectors'
           call CAMBParams_Set(P)
           if (global_error_flag==0) call cmbmain(im)
           call_again = .true.
           CP%Transfer%high_precision = Params%Transfer%high_precision
           CP%WantTransfer = Params%WantTransfer
           CP%WantScalars = Params%WantScalars
           CP%WantVectors = Params%WantVectors
           CP%Transfer%num_redshifts = Params%Transfer%num_redshifts
           Params = CP         
  
           Cl_scalar_nmodes(:, :, :,im) =  Cl_scalar(:, :, :)*scale_factor
           Cl_vector_nmodes(:, :, :,im) =  Cl_vector(:, :, :)*scale_factor*vec_factor
           Cl_tensor_nmodes(:, :, :,im) =  Cl_tensor(:, :, :)*scale_factor*tens_factor
           Cl_scalar_mean(:, :, :) = Cl_scalar_mean(:, :, :) + Cl_scalar_nmodes(:, :, :,im)
           Cl_vector_mean(:, :, :) = Cl_vector_mean(:, :, :) + Cl_vector_nmodes(:, :, :,im)
           Cl_tensor_mean(:, :, :) = Cl_tensor_mean(:, :, :) + Cl_tensor_nmodes(:, :, :,im)

           do il=2,P%Max_l
              scalar_power(im) = scalar_power(im) + Cl_scalar_nmodes(il,1,1,im)/il*2.0*pi
              vector_power(im) = vector_power(im) + Cl_vector_nmodes(il,1,1,im)/il*2.0*pi
              tensor_power(im) = tensor_power(im) + Cl_tensor_nmodes(il,1,1,im)/il*2.0*pi
           end do
           if (im.gt.0) then
              scalar_power(im) = scalar_power(im) + scalar_power(im-1)
              vector_power(im) = vector_power(im) + vector_power(im-1)
              tensor_power(im) = tensor_power(im) + tensor_power(im-1)
           end if

           if (output_power) write(80,'(1I5,3E15.5)') im,scalar_power(im),vector_power(im),tensor_power(im)

        end do

        if(output_power) close(80)

        call_again = .false.
        Cl_scalar_nmodes(:,:,:,:) = 0.0d0
        Cl_vector_nmodes(:,:,:,:) = 0.0d0
        Cl_tensor_nmodes(:,:,:,:) = 0.0d0
        call deallocate_string()

        CP%WantScalars = .true.
        CP%WantTensors = .true.
        CP%WantVectors = .true.

        CALL TIMER_STOP
        CALL TIMER_ELAPSED(ELAPSED_TIME)
        PRINT *, ' TIME ELAPSED = '//ELAPSED_TIME

        end subroutine CAMB_GetResults_string


!---------------------------------------------------------------------------------


  end module CAMB


  function dtda(a)
          use Precision
          implicit none
          real(dl) dtda,dtauda,a
          external dtauda
          
          dtda= dtauda(a)*a
  end function

        

