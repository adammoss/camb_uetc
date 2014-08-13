       module LambdaGeneral
         use precision
         implicit none
          
         real(dl)  :: w_lam = -1 !p/rho for the dark energy (assumed constant) 
         real(dl) :: cs2_lam = 1_dl 
          !comoving sound speed. Always exactly 1 for quintessence 
          !(otherwise assumed constant, though this is almost certainly unrealistic)

         logical :: w_perturb = .true.

       end module LambdaGeneral


    
!Return OmegaK - modify this if you add extra fluid components
        function GetOmegak()
        use precision
        use ModelParams
        real(dl)  GetOmegak
         GetOmegak = 1 - (CP%omegab+CP%omegac+CP%omegav+CP%omegan) 
          
        end function GetOmegak
  
  
       subroutine init_background
         !This is only called once per model, and is a good point to do any extra initialization.
         !It is called before first call to dtauda, but after
         !massive neutrinos are initialized and after GetOmegak
       end  subroutine init_background


!Background evolution
        function dtauda(a)
         !get d tau / d a
        use precision
        use ModelParams
        use MassiveNu
        use LambdaGeneral
        implicit none
        real(dl) dtauda
        real(dl), intent(IN) :: a
        real(dl) rhonu,grhoa2, a2
        integer nu_i

        a2=a**2

!  8*pi*G*rho*a**4.
        grhoa2=grhok*a2+(grhoc+grhob)*a+grhog+grhornomass
         if (w_lam == -1._dl) then
           grhoa2=grhoa2+grhov*a2**2
         else
           grhoa2=grhoa2+grhov*a**(1-3*w_lam)
         end if
        if (CP%Num_Nu_massive /= 0) then
!Get massive neutrino density relative to massless
           do nu_i = 1, CP%nu_mass_eigenstates
            call Nu_rho(a*nu_masses(nu_i),rhonu)
            grhoa2=grhoa2+rhonu*grhormass(nu_i)
           end do
        end if

        dtauda=sqrt(3/grhoa2)
     
        end function dtauda
