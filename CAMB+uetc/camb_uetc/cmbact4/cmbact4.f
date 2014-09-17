c CMBACT 4.0: evaluates CMB and matter spectra sourced by active 
c sources such as cosmic strings.
c--------------------------------------------------------------------
c Levon Pogosian, July 2013
c E-mail: levon@sfu.ca
c Latest updates at http://www.sfu.ca/~levon/cmbact.html
c--------------------------------------------------------------------
c Based on CMBFAST written by Uros Seljak and Mattias Zaldarriaga
c http://www.cmbfast.org
c -------------------------------------------------------------------
c Requires cmbfast.inc, jlgen.dat
c To create jlgen.dat compile and run jlgen.f
c when prompted, enter lmax=3000, akmax=6000 and name the output file jlgen.dat
c If other jlgen parameters are desired you will need to edit 
c subroutine initlval in this code.

      implicit double precision(a-h,o-z)
      parameter (lmax0=3300,nnmax=1,nk0=200)
      real*8 akout(nk0),trsum(nk0)
      common/trfun/nkount
      real*8 xclts(lmax0,nnmax),xcltt(lmax0,nnmax),xcltv(lmax0,nnmax)
      real*8 xclcs(lmax0,nnmax),xclct(lmax0,nnmax),xclcv(lmax0,nnmax)
      real*8 xcles(lmax0,nnmax),xclet(lmax0,nnmax),xclev(lmax0,nnmax)
      real*8 xclbv(lmax0,nnmax),xclbt(lmax0,nnmax)
      common/clout/xclts,xcltt,xcltv,
     &             xclcs,xclct,xclcv,
     &             xcles,xclet,xclev,
     &             xclbv,xclbt,
     &             akout,trsum
      common/strparam/tau_init,alf,vdev,dksi,ctilde,xlf,ns1
      common/seed/iseed
      common/volcom/tkmax     

      open(unit=65,file='../output/data/alpha1/cmbact/cl_tt_no_100.d',
     &	            status='unknown',form='formatted')
     
      open(unit=66,file='../output/data/alpha1/cmbact/cl_te_no_100.d',
     &	            status='unknown',form='formatted')
     
      open(unit=67,file='../output/data/alpha1/cmbact/cl_ee_no_100.d',
     &	            status='unknown',form='formatted')

      open(unit=70,file='../output/data/alpha1/cmbact/cl_bb_no_100.d',
     &	            status='unknown',form='formatted')

      open(unit=80,file='../output/data/alpha1/cmbact/pk_lin_no_100.d',
     &	            status='unknown',form='formatted')

c      write(*,*)'hello'

c set parameters:
c the string tension (the most important parameter)
       gmu=1.1d-6

c       gmu=0.0d0

c the number of string network realizations (>100 is good)
      nexp=100
      
c if string tension is zero, evaluate inflationary spectrum      
      if (gmu.eq.0.0d0) nexp=1

c maximum l to output (<= number used in jlgen.dat, currently 3000)
      lmaxout=3000

c Omega_B h^2
      ob=0.02249d0
      
c Omega_M h^2      
      om=0.112d0
      
c Hubble const      
      h=0.704d0
      
c Constant effective wiggliness parameter                  
      alf=1.0d0

c vdev is the initial rms string velocity (over all scales)
      vdev=0.65d0
      
c tkmax is the cutoff above which \Theta_D and \Theta_P are set to zero
      tkmax=dble(lmaxout)/2.0d0

c dksi is the initial correlation length/initial conformal time
      dksi=0.13d0

c ctilde is the loop chopping efficiency (the only parameter in VOS)
      ctilde=0.23d0
      
c ns1 is the number of consolidated string segments
      ns1=200

c xlf controls the rate at which strings decay
      xlf=0.95d0
      
c  initialize the random # generator      
      iseed=-19

c  tau_init is the earliest conformal time at which one scale model is run   
      tau_init=2.d-2

      do ik=1,nk0
      akout(ik)=0.0d0
      enddo

      call stringcmb(nexp,ob,om,h,gmu)

      do l=2,lmaxout      
      write(65,'(1I5,3E15.5)')
     & l,xclts(l,nnmax),xcltv(l,nnmax),xcltt(l,nnmax)
      write(66,'(1I5,3E15.5)')
     & l,xclcs(l,nnmax),xclcv(l,nnmax),xclct(l,nnmax)
      write(67,'(1I5,3E15.5)')
     & l,xcles(l,nnmax),xclev(l,nnmax),xclet(l,nnmax)
      write(70,'(1I5,2E15.5)')
     & l,xclbv(l,nnmax),xclbt(l,nnmax)
      enddo

      do ik=1,nkount      
      write(80,*)akout(ik),trsum(ik)      
      enddo    
    
      stop
      end

	subroutine stringcmb(nexp,ob,om,hub,gmu)
	implicit double precision(a-h,o-z)
      parameter (lmax0=3300,nnmax=1,lmax=20+lmax0/50)
	real*8 ztf(nnmax)

      common/strparam/tau_init,alf,vdev,dksi,ctilde,xlf,ns1
c
cLP Output arrays are scalar, vector and tensor temperature anisotropy:
c   clts, cltv, cltt and the power spectrum: trsum

        parameter (nk0=200)
        real*8 akout(nk0),trsum(nk0),trA(nk0)
        real*8 xclts(lmax0,nnmax),xcltt(lmax0,nnmax)
        real*8 xclcs(lmax0,nnmax),xclct(lmax0,nnmax)
        real*8 xclet(lmax0,nnmax),xclbt(lmax0,nnmax)
        real*8 xclev(lmax0,nnmax),xclbv(lmax0,nnmax)
        real*8 xcles(lmax0,nnmax),xcltv(lmax0,nnmax)

        real*8 xclcv(lmax0,nnmax)
        real*8 xclbtv(lmax0,nnmax),xclbtt(lmax0,nnmax)
        real*8 xclebv(lmax0,nnmax),xclebt(lmax0,nnmax)

      common/clout/xclts,xcltt,xcltv,
     &             xclcs,xclct,xclcv,
     &             xcles,xclet,xclev,
     &             xclbv,xclbt,
     &             akout,trsum

        real*8 cl2av(lmax0,nnmax),clerr(lmax0,nnmax)
        real*8 cltsA(lmax0,nnmax),clttA(lmax0,nnmax)
        real*8 cltvA(lmax0,nnmax),clesA(lmax0,nnmax)
        real*8 clcsA(lmax0,nnmax),clctA(lmax0,nnmax)
        real*8 cletA(lmax0,nnmax),clbtA(lmax0,nnmax)
        real*8 clevA(lmax0,nnmax),clbvA(lmax0,nnmax)

        real*8 clcvA(lmax0,nnmax)
        real*8 clbtvA(lmax0,nnmax),clbttA(lmax0,nnmax)
        real*8 clebvA(lmax0,nnmax),clebtA(lmax0,nnmax)



	integer l(lmax),l0
	common /lvalues1/ l,l0,lmo
	save /lvalues1/
	common /lingerinc/ omegab,omegac,omegav,omegan,omegak,h0
     &                     ,tcmb,yhe,nnur,nnunr
	common /genparm/ grhom,grhog,grhor,adotrad,taurst,dtaurec
	common /initcase/ initfl
      common /reionization/zri,taurist,zristp,tauristp,rif

      common /initialps/ an(nnmax),nn
      common /tensor/ ant(nnmax),rat(nnmax),itflag,lmaxt
      common /transfer/akmaxt,ztf,nlnkt,ict,ntf
	common /transfer2/ictpres
c	
	common/trfun/nkount
	
   
c	write(*,*)'CMB (0), transfer functions (1) or both (2):'
c	read(*,*)ict
	ict=2

	if (ict.ne.0) then
c	   write(*,*)'If you want a high precision (1%)in' 
c           write(*,*)'the transfer function at small scales'
c	   write(*,*)'enter 1. (The precision in the Cls is' 
c           write(*,*)'not affected by this choice but the code'
c	   write(*,*)'will become slower)'
c	   read(*,*)ictpres
       ictpres=0
	end if
	if (ict.ne.1) then
	   call initlval
	endif

       if (ict.ne.0) then
c        write(*,*)'Enter number and redshifts of the output: (1,10)'
c        write(*,*)'If more than one tf is requested the redshifts'
c        write(*,*)'have to be specified in decreasing order.'
c        read(*,*)ntf,(ztf(i),i=1,ntf)
         ntf=1
         ztf(1)=0.1d0
       endif
1      continue
c  Read initial parameters.
c     write(*,*)
c     2  'Enter Omega_b, Omega_c, Omega_v, Omega_nu (e.g. .05 .95 0 0)'
c	read(*,*) omegab,omegac,omegav,omegan

      omegab=ob/hub**2
      omegac=om/hub**2-omegab
      omegav=1.0d0-omegab-omegac
      omegan=0.0d0

c      write(*,*)omegab,omegac,omegav,hub
c
	omega=omegab+omegac+omegav+omegan
	omegak=1.0d0-omega
	if (abs(omegak).gt.1.d10) then
	  write(*,*) 'Currently works only for a flat universes'
	  write(*,*) '  You request Omega_curvature=',omegak
	  write(*,*) '  Illegal choice.  Try again:'
	  goto 1
	end if
        if (abs(omegak).gt.1.d10.and.omegan.gt.0) then
         write(*,*) 'Currently Omega_nu>0 works only in flat universe'
         write(*,*) ' You request Omega_nu,Omega_curv=',omegan,omegak
         write(*,*) '  Illegal choice.  Try again:'
         goto 1
        end if

2      continue
c       write(*,*)
c     2    'Enter  H0, Tcmb, Y_He,N_nu(massive) (e.g. 50 2.726 0.24 0)'
c	read(*,*) h0,tcmb,yhe,nnunr
c
       h0=hub*100.0d0
       tcmb=2.725d0
       yhe=0.24d0
       nnunr=0
c
	if (h0.lt.25.d0.or.h0.gt.100.d0) then
	  write(*,*)
     2      '  Warning: H0 has units of km/s/Mpc.  Your value is weird.'
	end if
	if (tcmb.lt.2.7d0.or.tcmb.gt.2.8d0) then
	  write(*,*)
     2      '  Warning: Tcmb has units of K.  Your value is weird.'
	end if

	if (yhe.lt.0.2d0.or.yhe.gt.0.3d0) then
	  write(*,*)
     2      '  Warning: Y_He is the Helium fraction of baryons.',
     3      '  Your value is weird.'
	end if
        if (nnunr.lt.0.or.nnunr.gt.3) then
          write(*,*)
     2      'Warning: N_nu(massive) should be 0, 1, 2, or 3'
          write(*,*) '  Illegal choice.  Try again:'
          go to 2
        end if

c	write (*,*) 'Enter 0 for no reionization'
c	write (*,*) 'Enter 1 for specified optical depth to lss (xe=1)'
c	write (*,*) 'Enter 2 for specified redshift and xe'
c	read (*,*) riflag
      riflag=1
c
c  Evaluate general cosmology constants.
        nnur=3-nnunr
c  grho gives the contribution to the expansion rate from: (g) photons,
c  (r) one flavor of relativistic neutrino (2 degrees of freedom),
c  (m) nonrelativistic matter (for Omega=1).  grho is actually
c  8*pi*G*rho/c^2 at a=1, with units of Mpc**(-2).
c  a=tau(Mpc)*adotrad, with a=1 today, assuming 3 neutrinos.
c  (Used only to set the initial conformal time.)
	grhom=3.3379d-11*h0*h0
	grhog=1.4952d-13*tcmb**4
	grhor=3.3957d-14*tcmb**4

c  adotrad gives the relation a(tau) in the radiation era:
	adotrad=2.8948d-7*tcmb*tcmb

	zri=0.0d0
	rif=0.0d0

	if (riflag.eq.1) then
c	   write (*,*) 'Enter optical depth to lss'
c	   read (*,*) optdlss
         optdlss=0.09d0
	   rif=1.0d0
c Calculating when reionization starts.
           call reiopar(optdlss,zri,zristp,rif)
	end if
c         write(*,*)optdlss,zri,zristp,rif

	if (riflag.eq.2) then
c	   write (*,*)'Enter redshift, ionization fraction(e.g. 50 0.2)'
c	   read (*,*) zri,rif

         zri=20.0d0
         rif=0.2d0

	   zristp=0.07d0*zri-1.0d0
	   if (zristp.lt.0.0) zristp=0.0d0
	end if

	if (ict.ne.1) then
c	      if (abs(omegak).gt.1.0d-3.and.itflag.ne.1) then
c		 write(*,*) 'tensor modes only allowed in flat models'
c	       itflag=
c	    else
c	      write (*,*) 'Enter 0 for scalar alone,'
c	      write(*,*) '1 for tensor+scalar+vector, 2 for tensors alone'
c             wriet(*,*)' or 3 for vectors alone'
c	      read(*,*)itflag
            if(gmu.eq.0.0d0) then
	     itflag=0
	    else
	     itflag=1
	    endif	
	    
c	   endif
	   if (itflag.ne.2.and.itflag.ne.3) then	 
c	      write(*,*)'number and values of scal. spec. ind. n(1,1)'
c	      read(*,*)nn,(an(i),i=1,nn)
cLP  Spectral index of 4 will give white noise intitial power spectrum
            nn=1
	    if(gmu.eq.0.0d0) then
             an(1)=0.955d0
	    else
	     an(1)=4.0d0
	    endif 
	   end if
	   if (itflag.eq.1) then 
c	      write(*,*)
c     2 'Tensor spectral index given by nt=ns-1 (0) or different (1)'
c	      read(*,*)itn
           itn=0
cLP Tensor and Vector spectral index are the same and set to give white noise.
	      if (itn.eq.0) then
		 do 123 in=1,nn
		    ant(in)=an(in)-1
123		 continue
	      else
		 write(*,*)'values of tensor spectral indexes:'
		 read(*,*)(ant(in),in=1,nn)
	      endif
	   endif
	   if (itflag.eq.2.or.itflag.eq.3) then
c	      write(*,*)'number and values of tens. spec. ind. (1,0)'
c	      read(*,*)nn,(ant(in),in=1,nn)
            nn=1
            ant(1)=3.0d0
	   end if
c
	else
	   itflag=0
	   nn=0
	end if

c	if (itflag.ne.2.and.itflag.ne.3) then
c	   if (abs(omegak).gt.1.0d-3) then
c	      write (*,*) 'CURRENTLY ONLY ISENTROPIC INITIAL CONDITIONS'
c	      write(*,*) 'ARE  AVAILABLE FOR OPEN MODELS'
c	      initfl=1
c	      write (*,*) 'Enter initial conditions'
c	      write (*,*) '1= Isentropic (adiabatic)'
c	      write (*,*) '2= Isocurvature CDM'
c	      write (*,*) '3= Isocurvature baryon'
c	      write (*,*) '4= Isocurvature seed conditions'
c	      read (*,*) initfl
c	   else
c	      write (*,*) 'Enter initial conditions'
c	      write (*,*) '1= Isentropic (adiabatic)'
c	      write (*,*) '2= Isocurvature CDM'
c	      write (*,*) '3= Isocurvature baryon'
c	      write (*,*) '4= Isocurvature seed conditions'
c	      read (*,*) initfl
cLP  The intitial conditions are the isocurvature seed conditions with
c    all intitial perturbations set to zero. Initial conditions are not
c    so important for this particular model if you start early enough.
            if(gmu.eq.0.0d0) then
	     initfl=1
	    else
             initfl=4
	    endif
c	   end if
c	else
c	   initfl=1
c	end if
c
	  if (ict.ne.1) then
cLP jlgen.dat is read in initjl	  
	     call initjl
	  end if
     
cLP If vector part is requested precalculate certain numerical factors
      if(itflag.ne.0) then
      call savetime
      endif
c
cLP Initialize all arrays used in the string model and evaluate the
c   parameters of the one-scale model:

      if (gmu.ne.0.0d0) call prepare

	do ik=1,nk0
        trsum(ik)=0.0d0
        enddo

	do in=1,nnmax
	   do il=1,lmax0
	      clerr(il,in)=0.0d0
	      cl2av(il,in)=0.0d0
         xclts(il,in)=0.0d0
         xcles(il,in)=0.0d0
	 xcltt(il,in)=0.0d0
         xclet(il,in)=0.0d0
	 xclbt(il,in)=0.0d0
         xclev(il,in)=0.0d0
	 xclbv(il,in)=0.0d0
         xclcs(il,in)=0.0d0
	 xclct(il,in)=0.0d0
	 xcltv(il,in)=0.0d0
         xclcv(il,in)=0.0d0
	 xclbtv(il,in)=0.0d0
	 xclbtt(il,in)=0.0d0
	 xclebv(il,in)=0.0d0
	 xclebt(il,in)=0.0d0
	   end do
	end do


      do 100 iexp=1,nexp
      
c      write(*,*)iexp
c
cLP  First generate all random parameters of this particular realization
c    of the string network.
      if (gmu.ne.0.0d0) call mixup
c
	do in=1,nnmax
	   do il=1,lmax0
	    cltsA(il,in)=0.0d0
            clesA(il,in)=0.0d0
	    clttA(il,in)=0.0d0
	    cltvA(il,in)=0.0d0
            clcsA(il,in)=0.0d0
            clctA(il,in)=0.0d0
	    cletA(il,in)=0.0d0
	    clbtA(il,in)=0.0d0
            clevA(il,in)=0.0d0
	    clbvA(il,in)=0.0d0
            clcvA(il,in)=0.0d0
	    clbtvA(il,in)=0.0d0
	    clbttA(il,in)=0.0d0
            clebvA(il,in)=0.0d0
	    clebtA(il,in)=0.0d0
	   end do
	end do
c
cLP cmbflat calculates CMB anisotropy for this particular
c   realization of the string network. The string stress-energy
c   tensor is calculated by subroutine strings_splined called from cmbflat

      call cmbflat(cltsA,clttA,cltvA,clesA,clevA,clbvA,
     &             cletA,clbtA,clcsA,clctA,clcvA,clbtvA,
     &             clbttA,clebvA,clebtA,trA,akout)


c
cLP calculate the average power and the average square of power
c   to find the error bar.
      do ik=1,nkount
      trsum(ik)=trsum(ik)+trA(ik)/dble(nexp)
      enddo


     	do in=1,nnmax
	  do il=2,lmax0
	  
c	write(70+iexp,*)il,cltsA(il,in)
	  
cLP  Calculate the average Cl's and
c    the average of cl^2 to find the error bar
         cl2av(il,in)=cl2av(il,in)+
     & (cltsA(il,in)+cltvA(il,in)+clttA(il,in))**2/dble(nexp)

       xclts(il,in)=xclts(il,in)+cltsA(il,in)/dble(nexp)
       xcles(il,in)=xcles(il,in)+clesA(il,in)/dble(nexp)
       xcltt(il,in)=xcltt(il,in)+clttA(il,in)/dble(nexp)
       xcltv(il,in)=xcltv(il,in)+cltvA(il,in)/dble(nexp)
       xclet(il,in)=xclet(il,in)+cletA(il,in)/dble(nexp)
       xclbt(il,in)=xclbt(il,in)+clbtA(il,in)/dble(nexp)
       xclev(il,in)=xclev(il,in)+clevA(il,in)/dble(nexp)
       xclbv(il,in)=xclbv(il,in)+clbvA(il,in)/dble(nexp)
       xclcs(il,in)=xclcs(il,in)+clcsA(il,in)/dble(nexp)
       xclct(il,in)=xclct(il,in)+clctA(il,in)/dble(nexp)
       xclcv(il,in)=xclcv(il,in)+clcvA(il,in)/dble(nexp)
       xclbtv(il,in)=xclbtv(il,in)+clbtvA(il,in)/dble(nexp)
       xclbtt(il,in)=xclbtt(il,in)+clbttA(il,in)/dble(nexp)
       xclebv(il,in)=xclebv(il,in)+clebvA(il,in)/dble(nexp)
       xclebt(il,in)=xclebt(il,in)+clebtA(il,in)/dble(nexp)
       
       
c	write(74+iexp,*)il,xclts(il,in)       

        end do
	end do

100   continue



cLP  if igmu=0, then normalize to COBE, if igmu=1 need a value of G\mu
      if(gmu.ne.0.0d0) igmu=1

      call COBEnormalize(igmu,gmu,xclts,xcltt,xcltv,cl2av,
     & clerr,xcles,
     & xclev,xclbv,xclet,xclbt,xclcs,xclct,xclcv,xclbtv,xclbtt,
     & xclebv,xclebt,trsum,akout)


777   continue

      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine prepare
      implicit double precision (a-h,o-z)
      parameter(nstep00=200)

      common/picom/pi,twopi3

      common/strparam/tau_init,alf,vdev,dksi,ctilde,xlf,ns1

      common/evolved/evv(nstep00),evvpr(nstep00)
     &              ,cole(nstep00),colepr(nstep00)
c     &              ,evalf(nstep00),evalfpr(nstep00)
     &              ,cnorm(nstep00),cnormpr(nstep00)
     &              ,tphys(nstep00),tphyspr(nstep00)     
     &              ,tim(nstep00),timpr(nstep00)
     &              ,xaa(nstep00),xaapr(nstep00),ntime

      common /lingerinc/ omegab,omegac,omegav,omegan,omegak,h0
     &                     ,tcmb,yhe,nnur,nnunr
      common /genparm/ grhom,grhog,grhor,adotrad,taurst,dtaurec

      parameter (nw=3)
      parameter (d0hi=1.0d40,d0lo=1.0d40)

      logical do_evolve
      real*8 yev(nw),yevpr(nw)
      dimension c(24),w(nw,9)
              double precision dtauda,adtauda
      external fevolve
      external dtauda
      external adtauda

      twopi3=(2.0d0*pi)**3
      pi=4.0d0*datan(1.0d0)

      do_evolve = .false.
       
c  taui is the earliest conformal time at which one scale model is run   
      taui=tau_init

c  calculate conformal time at 10 times scale factor today
      tol1=1.d-4
      tauf=rombint(dtauda,0.0d0,10.0d0,tol1)

c  Set the time grid for the slow functions of time
c  define tim(it) logarithmically spaced
      ntime=200

      dlnt=log(tauf/taui)/dble(ntime-1)
      do 10 it=1,ntime
      tim(it)=taui*exp(dble(it-1)*dlnt)
10    continue

      tinit=0.1d0*tim(1)      
      ainit=adotrad*tinit

c  initial correlation length of the network
      xl0=dksi*tinit

      nvar=3
      yev(1)=adotrad*tinit
      yev(2)=xl0
      yev(3)=vdev
      tau=tinit
      ind=1
      tol1=1.0d-6
      tol2=1.0d-4
c  evolve and save the paramters of the one-scale model
c  and their derivatives.

      do 20 j=1,ntime
      tauend=tim(j)

        if (do_evolve) then
          call dverk(nvar,fevolve,tau,yev,tauend,tol1,ind,c,nw,w)
          call fevolve(nvar,tauend,yev,yevpr)
       
          adot=yevpr(1)
          xl=yev(2)
          vel=yev(3)
          else
          xl = xl0
          vel = vdev
        endif

c this is used to store t as a function of tau

      tphys(j)=rombint(adtauda,0.d0,yev(1),tol2)

c  the scale factor
      xaa(j)=yev(1)     

c  wiggliness
c      evalf(j)=alf    
c
c  correlation lenght
      cole(j)=xl
      
c  rms velocity
      evv(j)=vel        

c      write(58,'(4E15.5)')xaa(j),cole(j)/tim(j),
c     & cole(j)*xaa(j)/tphys(j),evv(j)

20    continue

c Spline the one-scale functions of time, also the physical time.      

      call spline(tphys,tim,ntime,d0lo,d0hi,timpr)
      call spline(tim,tphys,ntime,d0lo,d0hi,tphyspr)
      call spline(tim,cole,ntime,d0lo,d0hi,colepr)
c      call spline(tim,evalf,ntime,d0lo,d0hi,evalfpr)      
      call spline(tim,evv,ntime,d0lo,d0hi,evvpr)     
      call spline(tim,xaa,ntime,d0lo,d0hi,xaapr)
      
      return
      end

      subroutine fevolve(n,x,y,yprime)
c this subroutine evaluates the one-scale model equations
      implicit double precision (a-h,o-z)
      common/picom/pi,twopi3
      dimension y(n),yprime(n)

      common /lingerinc/ omegab,omegac,omegav,omegan,omegak,h0
     &                     ,tcmb,yhe,nnur,nnunr
      common /genparm/ grhom,grhog,grhor,adotrad,taurst,dtaurec

      common/strparam/tau_init,alf,vdev,dksi,ctilde,xlf,ns1

      tau=x
      a=y(1)     
      a2=a*a

      rhonu=0.0d0
      grho=grhom*(omegac+omegab)/a+(grhog+grhor*(nnur+nnunr*rhonu))/a2
     2      +grhom*omegav*a2

      rho=grho/a2
      addot=grhom*(0.5d0*(omegac+omegab)+2.0d0*omegav*a**3)/3.0d0

      adotoa=sqrt(grho/3.0d0)
      yprime(1)=adotoa*a
      xl=y(2)
      xv=y(3)
      xv2=xv*xv

c  cr is the only free parameter

      fk=(2.0d0*sqrt(2.0d0)/pi)*
     & (1.0d0-8.0d0*xv**6)/(1.0d0+8.0d0*xv**6)

      xlprime=adotoa*xl*xv2+0.5d0*ctilde*xv
      xvprime=(1.0d0-xv2)*(fk/xl-2.0d0*adotoa*xv)

      yprime(2)=xlprime
      yprime(3)=xvprime

      return
      end

      subroutine mixup
      implicit double precision (a-h,o-z)
      parameter(nstep00=200,ns0=801)

      common/seed/iseed
      common/picom/pi,twopi3
      parameter (d0hi=1.0d40,d0lo=1.0d40)
      common/strparam/tau_init,alf,vdev,dksi,ctilde,xlf,ns1

      common/evolved/evv(nstep00),evvpr(nstep00)
     &              ,cole(nstep00),colepr(nstep00)
c     &              ,evalf(nstep00),evalfpr(nstep00)
     &              ,cnorm(nstep00),cnormpr(nstep00)
     &              ,tphys(nstep00),tphyspr(nstep00)     
     &              ,tim(nstep00),timpr(nstep00)
     &              ,xaa(nstep00),xaapr(nstep00),ntime

      common/mixcom/xp3(ns0,nstep00),xp3pr(ns0,nstep00)
     &             ,xd3(ns0,nstep00),xd3pr(ns0,nstep00)
     &             ,ps(ns0,nstep00),pspr(ns0,nstep00)
     &             ,pv(ns0,nstep00),pvpr(ns0,nstep00)
     &             ,pt(ns0,nstep00),ptpr(ns0,nstep00)   
     &             ,x0k(ns0),tf(ns0),xlden(ns0)

c      write(*,*)iseed	

c   Set the times at which the string segments decay: tf()

      taui=tim(1)
      tf(1)=taui

	call splint1(tim,cole,colepr,ntime,taui,xl)
      xlden(1)=xl

      tau2=tf(1)*(1.0d0+ran1(iseed))
      tauf=tim(ntime)*(1.0d0+ran1(iseed))      

      dlnt=log(tauf/tau2)/dble(ns1-1)

      do 15 m=2,ns1+1

      tau=tau2*exp(dble(m-2)*dlnt)
      tf(m)=tau
	call splint1(tim,cole,colepr,ntime,tau,xl)
      xlden(m)=xl

15    continue

c  normalize the string density to 1/l**3
c  and save the normalization and it's derivative.
      do 30 it=1,ntime
      sum=0.0d0
      tau=tim(it)

      do 40 m=2,ns1
      bnum=1.0d0/xlden(m-1)**3-
     &                  1.0d0/xlden(m)**3
      call toff(tau,tf(m),xlf,off,offdot)
      sum=sum+bnum*off
40    continue

      sum=sum+1.0d0/xlden(ns1+1)**3

      cnorm(it)=1.0d0/(sum*cole(it)**3)

30    continue

      call spline(tim,cnorm,ntime,d0lo,d0hi,cnormpr)

c------------------------------------

c Initial time is same for all segments
      taui=tf(1)

c final times are also same
      tauf=tf(ns1+1)
     
c  Loop over string segments 
      do 300 m=2,ns1+1

c  each segment gets a random intitial phase
      x0k(m)=2.0d0*pi*ran1(iseed)

c  Generate directions at random
      cteta=2.0d0*ran1(iseed)-1.0d0
      phi=2.0d0*pi*ran1(iseed)
      psi=2.0d0*pi*ran1(iseed)   

c For each segment, pre-calculate all slowly varying functions of time

      do it=1,ntime
      tau=tim(it)

      steta=sqrt(1.0d0-cteta**2)

      xp1=steta*sin(phi)
      xp2=-steta*cos(phi)
      xp3(m,it)=cteta
      xd1=cos(psi)*cos(phi)-sin(psi)*sin(phi)*cteta
      xd2=cos(psi)*sin(phi)+sin(psi)*cos(phi)*cteta
      xd3(m,it)=steta*sin(psi)    
      
c  now the prefactors of S, V and T components
      vm=evv(it)
      vm2=evv(it)**2
      ovm2=1.0d0-vm2
c      alf=evalf(it)

      ps(m,it)=0.5d0*(vm2*(3.0d0*xd3(m,it)*xd3(m,it)-1.0d0)-
     #       ovm2*(3.0d0*xp3(m,it)*xp3(m,it)-1.0d0)/alf**2)     
   
      pv(m,it)=vm2*xd1*xd3(m,it)-ovm2*xp1*xp3(m,it)/alf**2
      pt(m,it)=vm2*xd1*xd2-ovm2*xp1*xp2/alf**2           
	  
	  enddo

	  call spline(tim,xp3(m,1),ntime,d0lo,d0hi,xp3pr(m,1))
	  call spline(tim,xd3(m,1),ntime,d0lo,d0hi,xd3pr(m,1))
	  call spline(tim,ps(m,1),ntime,d0lo,d0hi,pspr(m,1))
	  call spline(tim,pv(m,1),ntime,d0lo,d0hi,pvpr(m,1))
	  call spline(tim,pt(m,1),ntime,d0lo,d0hi,ptpr(m,1))	  
	 	 
300   continue  

      return
      end

      subroutine strings(wk,tau,T00,TS,TV,TT)
c  Components of string stress-energy and the derivatives will be calculated 
c  for the wave vector k=wk and conformal time tau
      implicit double precision (a-h,o-z)
      parameter(nstep00=200,ns0=801)

      common/picom/pi,twopi3

      common/strparam/tau_init,alf,vdev,dksi,ctilde,xlf,ns1

      common/evolved/evv(nstep00),evvpr(nstep00)
     &              ,cole(nstep00),colepr(nstep00)
c     &              ,evalf(nstep00),evalfpr(nstep00)
     &              ,cnorm(nstep00),cnormpr(nstep00)
     &              ,tphys(nstep00),tphyspr(nstep00)     
     &              ,tim(nstep00),timpr(nstep00)
     &              ,xaa(nstep00),xaapr(nstep00),ntime

      common/mixcom/xp3(ns0,nstep00),xp3pr(ns0,nstep00)
     &             ,xd3(ns0,nstep00),xd3pr(ns0,nstep00)
     &             ,ps(ns0,nstep00),pspr(ns0,nstep00)
     &             ,pv(ns0,nstep00),pvpr(ns0,nstep00)
     &             ,pt(ns0,nstep00),ptpr(ns0,nstep00)   
     &             ,x0k(ns0),tf(ns0),xlden(ns0)
	 

c  Components of string stress-energy and the derivatives will be calculated 
c  for the wave vector k=wk and conformal time tau 
c
c  Other two scalar components are evaluated from the conservation
c  equations in the subroutine fderivs.


      T00=0.0d0
c      T00dot=0.0d0
      TS=0.0d0
c      TSdot=0.0d0
      TT=0.0d0
      TV=0.0d0  

c parameters required by the cubic spline
      klo=1
      khi=ntime
111   if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(tim(k).gt.tau)then
          khi=k
        else
          klo=k
        endif
      goto 111
      endif
      h=tim(khi)-tim(klo)
      if (h.eq.0.d0) then
      write(*,*)'bad xa input in splint'
      stop
      endif
      a=(tim(khi)-tau)/h
      b=(tau-tim(klo))/h
c end of cubic spline prep

c find v,l,alpha and their derivatives

      vm=a*evv(klo)+b*evv(khi)
     &    +((a**3-a)*evvpr(klo)+(b**3-b)*evvpr(khi))*(h**2)/6.d0

c      vmdot=(evv(khi)-evv(klo))/h
c     &       -(3.d0*a**2-1.0d0)*h*evvpr(klo)/6.0d0
c     &       +(3.d0*b**2-1.0d0)*h*evvpr(khi)/6.0d0

      xl=a*cole(klo)+b*cole(khi)
     &    +((a**3-a)*colepr(klo)+(b**3-b)*colepr(khi))*(h**2)/6.d0

c      xldot=(cole(khi)-cole(klo))/h
c     &       -(3.d0*a**2-1.0d0)*h*colepr(klo)/6.0d0
c     &       +(3.d0*b**2-1.0d0)*h*colepr(khi)/6.0d0

c      alf=a*evalf(klo)+b*evalf(khi)
c     &    +((a**3-a)*evalfpr(klo)+(b**3-b)*evalfpr(khi))*(h**2)/6.d0

c      alfdot=(evalf(khi)-evalf(klo))/h
c     &       -(3.d0*a**2-1.0d0)*h*evalfpr(klo)/6.0d0
c     &       +(3.d0*b**2-1.0d0)*h*evalfpr(khi)/6.0d0

c       snorm=1.0d0
c       snormdot=0.0d0

      snorm=a*cnorm(klo)+b*cnorm(khi)
     &    +((a**3-a)*cnormpr(klo)+(b**3-b)*cnormpr(khi))*(h**2)/6.d0

c      snormdot=(cnorm(khi)-cnorm(klo))/h
c     &       -(3.d0*a**2-1.0d0)*h*cnormpr(klo)/6.0d0
c     &       +(3.d0*b**2-1.0d0)*h*cnormpr(khi)/6.0d0

      
	  osv=1.0d0/sqrt(1.0d0-vm*vm)
c	  osvdot=osv**3*vm*vmdot
      
      do 100 m=2,ns1+1

c evaluate the splined functions and their derivatives
      xpr3=a*xp3(m,klo)+b*xp3(m,khi)
     &    +((a**3-a)*xp3pr(m,klo)+(b**3-b)*xp3pr(m,khi))*(h**2)/6.d0

c      xpr3dot=(xp3(m,khi)-xp3(m,klo))/h
c     &       -(3.d0*a**2-1.0d0)*h*xp3pr(m,klo)/6.0d0
c     &       +(3.d0*b**2-1.0d0)*h*xp3pr(m,khi)/6.0d0 

      xdt3=a*xd3(m,klo)+b*xd3(m,khi)
     &    +((a**3-a)*xd3pr(m,klo)+(b**3-b)*xd3pr(m,khi))*(h**2)/6.d0

c      xdt3dot=(xd3(m,khi)-xd3(m,klo))/h
c     &       -(3.d0*a**2-1.0d0)*h*xd3pr(m,klo)/6.0d0
c     &       +(3.d0*b**2-1.0d0)*h*xd3pr(m,khi)/6.0d0 

      xps=a*ps(m,klo)+b*ps(m,khi)
     &    +((a**3-a)*pspr(m,klo)+(b**3-b)*pspr(m,khi))*(h**2)/6.d0

c      xpsdot=(ps(m,khi)-ps(m,klo))/h
c     &       -(3.d0*a**2-1.0d0)*h*pspr(m,klo)/6.0d0
c     &       +(3.d0*b**2-1.0d0)*h*pspr(m,khi)/6.0d0

      xpv=a*pv(m,klo)+b*pv(m,khi)
     &    +((a**3-a)*pvpr(m,klo)+(b**3-b)*pvpr(m,khi))*(h**2)/6.d0

      xpt=a*pt(m,klo)+b*pt(m,khi)
     &    +((a**3-a)*ptpr(m,klo)+(b**3-b)*ptpr(m,khi))*(h**2)/6.d0
c end spline evaluations

c      if(m.eq.15) then
c      write(*,*)tau,xpr3,xpr3dot,xdt3,xdt3dot
c      endif

      call toff(tau,tf(m),xlf,off,offdot)

      fnorm=sqrt(snorm/xlden(m-1)**3-
     &                         snorm/xlden(m)**3)

c     fnormdot=snormdot*0.5d0*fnorm/snorm


      if(m.eq.ns1+1) then     
      fnorm=sqrt(snorm/xlden(m)**3)
      endif

      ampl=fnorm*off

      if(ampl.le.1.0d-13) goto 100
c      ampldot=fnormdot*off+fnorm*offdot

      var=0.5d0*wk*xpr3*xl	  
c      vardot=0.5d0*wk*(xpr3*xldot+xpr3dot*xl)
	  
      phase=x0k(m)+wk*xdt3*vm*tau	  
c      phasedot=wk*(xdt3dot*vm*tau+xdt3*vmdot*tau+xdt3*vm)
	  
      bet=0.5d0*wk*xpr3
c      betdot=0.5d0*wk*xpr3dot

      sinV=sin(var)
      cosV=cos(var)
      sinP=sin(phase)
      cosP=cos(phase)

      trig=alf*osv*sinV*cosP/bet

c      trigdot=alfdot*osv*sinV*cosP/bet
c     &	 +alf*osvdot*sinV*cosP/bet
c     &       +alf*osv*(cosV*vardot)*cosP/bet
c     &       -alf*osv*sinV*(sinP*phasedot)/bet
c     &       -alf*osv*sinV*cosP*(betdot/bet**2)
	  
	  							
c The factor of sqrt(2) in the equations below is compensating
c for the fact that we evaluate only the real part of the Fourier
c transform of the stress-energy tensor. 

      edenre=sqrt(2.0d0)*(ampl*trig)
c      edot=sqrt(2.0d0)*(ampl*trigdot+ampldot*trig)

      T00=T00+edenre
c      T00dot=T00dot+edot
	  
      TS=TS+xps*edenre
c      TSdot=TSdot+xps*edot+xpsdot*edenre

      TV=TV+xpv*edenre
	  	  
      TT=TT+xpt*edenre

100   continue


      return
      end
      
      subroutine strings_spline(wk,taui,tauf)
      implicit double precision (a-h,o-z)  
      common/strparam/tau_init,alf,vdev,dksi,ctilde,xlf,ns1    
      parameter (d0hi=1.0d40,d0lo=1.0d40)
      parameter (nt_fine=5000)
      common/splined_strings/ts_f(nt_fine)
     &                      ,em00(nt_fine),em00pr(nt_fine)
     &                      ,emS(nt_fine),emSpr(nt_fine)
     &                      ,emD(nt_fine),emDpr(nt_fine)
     &                      ,emP(nt_fine),emPpr(nt_fine)          
     &                      ,emV(nt_fine),emVpr(nt_fine)
     &                      ,emT(nt_fine),emTpr(nt_fine)	    
     &                      ,nt_f
c      parameter(nstep00=200,ns0=801)
c      common/evolved/evv(nstep00),evvpr(nstep00)
c     &              ,cole(nstep00),colepr(nstep00)
c     &              ,evalf(nstep00),evalfpr(nstep00)
c     &              ,cnorm(nstep00),cnormpr(nstep00)
c     &              ,tphys(nstep00),tphyspr(nstep00)     
c     &              ,tim(nstep00),timpr(nstep00)
c     &              ,xaa(nstep00),xaapr(nstep00),ntime     

    
      ts_f(1)=0.1d0*taui

      if(ts_f(1).lt.tau_init) then
      write(*,*)'one scale model needs an earlier start time'
      stop
      endif

      dtau=0.0d0
      nt_f=1	

1     continue    

    
      dtau=0.1d0*ts_f(nt_f)          
      dt=min(dtau,2.0d0/wk)      
    
      ts_f(nt_f+1)=ts_f(nt_f)+dt

      nt_f=nt_f+1
      
      if(ts_f(nt_f).lt.tauf) goto 1
       

c	do i=1,nt_f
c	write(*,*)i,ts_f(i)
c	enddo

                
c      nt_f=int((tauf-taui)/dt)+1
      
c      write(*,*)wk,nt_f


      do i=1,nt_f

c      ts_f(i)=taui+dble(i-1)*dt
c      ts_f(i)=taui*exp(dble(i-1)*dlnt)	     
      
      tau=ts_f(i)
      
      call strings(wk,tau,T00,TS,TV,TT)
      
      em00(i)=T00
      emS(i)=TS
      emV(i)=TV
      emT(i)=TT
           
      enddo      

c spline Tmunu

      call spline(ts_f,em00,nt_f,d0lo,d0hi,em00pr)
      call spline(ts_f,emS,nt_f,d0lo,d0hi,emSpr)
      call spline(ts_f,emV,nt_f,d0lo,d0hi,emVpr)
      call spline(ts_f,emT,nt_f,d0lo,d0hi,emTpr)            

c find the other 2 scalar mode components from conservation           
      call emtconserve(wk)
      
      call spline(ts_f,emD,nt_f,d0lo,d0hi,emDpr) 
      call spline(ts_f,emP,nt_f,d0lo,d0hi,emPpr)           	    
      
      return
      end
      
      
      subroutine scalar_source(tau,t00,t00dot,ts,tsdot,td,tp)
      implicit double precision (a-h,o-z)  
      
      parameter (d0hi=1.0d40,d0lo=1.0d40)
      parameter (nt_fine=5000)
      common/splined_strings/ts_f(nt_fine)
     &                      ,em00(nt_fine),em00pr(nt_fine)
     &                      ,emS(nt_fine),emSpr(nt_fine)
     &                      ,emD(nt_fine),emDpr(nt_fine)
     &                      ,emP(nt_fine),emPpr(nt_fine)          
     &                      ,emV(nt_fine),emVpr(nt_fine)
     &                      ,emT(nt_fine),emTpr(nt_fine)	    
     &                      ,nt_f
     
     
c parameters required by the cubic spline
      klo=1
      khi=nt_f
111   if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(ts_f(k).gt.tau)then
          khi=k
        else
          klo=k
        endif
      goto 111
      endif
      h=ts_f(khi)-ts_f(klo)
      if (h.eq.0.d0) then
      write(*,*)'bad xa input in splint'
      stop
      endif
      a=(ts_f(khi)-tau)/h
      b=(tau-ts_f(klo))/h
c end of cubic spline prep


      t00=a*em00(klo)+b*em00(khi)
     &    +((a**3-a)*em00pr(klo)+(b**3-b)*em00pr(khi))*(h**2)/6.d0

      t00dot=(em00(khi)-em00(klo))/h
     &       -(3.d0*a**2-1.0d0)*h*em00pr(klo)/6.0d0
     &       +(3.d0*b**2-1.0d0)*h*em00pr(khi)/6.0d0     
     
      ts=a*emS(klo)+b*emS(khi)
     &    +((a**3-a)*emSpr(klo)+(b**3-b)*emSpr(khi))*(h**2)/6.d0    
     
      tsdot=(emS(khi)-emS(klo))/h
     &       -(3.d0*a**2-1.0d0)*h*emSpr(klo)/6.0d0
     &       +(3.d0*b**2-1.0d0)*h*emSpr(khi)/6.0d0
     
      td=a*emD(klo)+b*emD(khi)
     &    +((a**3-a)*emDpr(klo)+(b**3-b)*emDpr(khi))*(h**2)/6.d0  
     
      tp=a*emP(klo)+b*emP(khi)
     &    +((a**3-a)*emPpr(klo)+(b**3-b)*emPpr(khi))*(h**2)/6.d0                   
     
     
      return
      end        

c emt conservation
        subroutine emtconserve(qk)
c  Find the rest of the scalar components from energy conservation.        

	implicit double precision(a-h,o-z)
        parameter (nnmax=1)
     
	common /lingerinc/ omegab,omegac,omegav,omegan,omegak,h0
     &                     ,tcmb,yhe,nnur,nnunr     
     
	common/fstcom/ak,emtD,emtP
        common/volcom/tkmax
	
	common /genparm/ grhom,grhog,grhor,adotrad,taurst,dtaurec

c
c  The metric perturbation variables are computed algebraically.
c>>>>>>>>>>>>>>>>>>>>
c added 1 to the number of scalar variables        
        parameter (nvar0=2)

c<<<<<<<<<<<<<<<<<<<<      

c Scalar perturbations

	real*8 y(nvar0),yprime(nvar0)

      parameter (nt_fine=5000)
      common/splined_strings/ts_f(nt_fine)
     &                      ,em00(nt_fine),em00pr(nt_fine)
     &                      ,emS(nt_fine),emSpr(nt_fine)
     &                      ,emD(nt_fine),emDpr(nt_fine)
     &                      ,emP(nt_fine),emPpr(nt_fine)          
     &                      ,emV(nt_fine),emVpr(nt_fine)
     &                      ,emT(nt_fine),emTpr(nt_fine)	    
     &                      ,nt_f

c
c  dverk  parameters.
	parameter (nw=nvar0)
        dimension c(24),w(nw,9)
c
	external fstrings
c>>>>>>>>>>>>>>>>>
        ak=qk
        taustart=ts_f(1)
c
        tau=taustart
	a=tau*adotrad
	a2=a*a
c Initial values of emtP and emt D set to 0
        
	emP(1)=0.0d0
	emD(1)=0.0d0

c
        y(2)=0.0d0
	y(1)=a                
        ind=1
        nvar=2
	
                     if (ak.lt.0.2) then 
                        tol1=1.0d-8
                     else
                        tol1=1.0d-6
                     endif	
                
c  Begin timestep loop.

                do 10 j=2,nt_f
                
                   tauend=ts_f(j)

             if(ak*tau.le.tkmax) then

                      call dverk(nvar,fstrings,tau,y,tauend,tol1
     &                     ,ind,c,nw,w)

                      call fstrings(nvar,tau,y,yprime)

                      emP(j)=emtP
                      emD(j)=emtD 
             else
                      emP(j)=0.0d0
                      emD(j)=0.0d0	     
	     
	     endif		      
		      		      		      		          
10             end do

c
       return
       end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine fstrings(n,x,y,yprime)
c  Evaluate the time derivatives of the perturbations.
c
	implicit double precision (a-h,o-z)

	dimension y(n),yprime(n)

      parameter (nt_fine=5000)
      common/splined_strings/ts_f(nt_fine)
     &                      ,em00(nt_fine),em00pr(nt_fine)
     &                      ,emS(nt_fine),emSpr(nt_fine)
     &                      ,emD(nt_fine),emDpr(nt_fine)
     &                      ,emP(nt_fine),emPpr(nt_fine)          
     &                      ,emV(nt_fine),emVpr(nt_fine)
     &                      ,emT(nt_fine),emTpr(nt_fine)	    
     &                      ,nt_f

      
        common/picom/pi,twopi3

	common /lingerinc/ omegab,omegac,omegav,omegan,omegak,h0
     &                     ,tcmb,yhe,nnur,nnunr

	common/fstcom/ak,emtD,emtP
	
	common /genparm/ grhom,grhog,grhor,adotrad,taurst,dtaurec

	common /out1/ adotoa,hdot,dgshear,rhonu,shearnu
c
c     inverse curvature radius squared
c      hc=2.998d5/h0
c      curv=-omegak/hc/hc
c
	tau=x
	a=y(1)
	a2=a*a
	ak2=ak*ak
c
c  STRING ENERGY MOMENTUM TENSOR emt_\mu_\nu

       call splint1(ts_f,em00,em00pr,nt_f,tau,emt00)  
       call splint1(ts_f,emS,emSpr,nt_f,tau,emtS)  
       call splint2(ts_f,em00,em00pr,nt_f,tau,emt00dot)
       call splint2(ts_f,emS,emSpr,nt_f,tau,emtSdot)       

       emtD=y(2)

c
c  8*pi*G*rho*a**2 and 8*pi*G*P*a**2.
c	grho=grhom*(omegac+omegab)/a
c     1      +(grhog+grhor*(annur+annunr*rhonu))/a2
c     2      +grhom*omegav*a2+grhom*omegak
          
      grho=grhom*(omegac+omegab)/a+(grhog+grhor*(nnur+nnunr*rhonu))/a2
     2      +grhom*omegav*a2
     
	adotoa=sqrt(grho/3.0d0)
	yprime(1)=adotoa*a
 
      emtP=(emtD-emt00dot)/adotoa-emt00
      
      
c  String conservation equation

c      emtDdot=-2.0d0*adotoa*emtD-(ak**2/3.0d0)*(emtP
c     1            +2.0d0*(1.0d0-3.0d0*curv/ak2)*emtS)
     
      emtDdot=-2.0d0*adotoa*emtD-(ak2/3.0d0)*(emtP+2.0d0*emtS)
     
       yprime(2)=emtDdot       


	return
	end
      

c this subroutine turns the string segments off
      subroutine toff(tau,tauf,xlf,off,offdot)
      implicit double precision (a-h,o-z)
      var=2.0d0*log(xlf*tauf/tau)/log(xlf)-1.0d0
      if(tau.lt.xlf*tauf) then
      off=1.0d0
      offdot=0.0d0
      goto 1
      endif
      if(tau.ge.tauf) then
      off=0.0d0
      offdot=0.0d0
      goto 1
      endif
      off=0.5d0+0.25d0*(var**3-3.0d0*var)
      offdot=-1.5d0*(var**2-1.0d0)/(tau*log(xlf))
1     continue
      return
      end

      FUNCTION gasdev(idum)
      INTEGER idum
      real*8 gasdev
CU    USES ran1
      INTEGER iset
      real*8 fac,gset,rsq,v1,v2,ran1
      SAVE iset,gset
      DATA iset/0/
      if (iset.eq.0) then
1       v1=2.d0*ran1(idum)-1.
        v2=2.d0*ran1(idum)-1.
        rsq=v1**2+v2**2
        if(rsq.ge.1.d0.or.rsq.eq.0.d0)goto 1
        fac=sqrt(-2.d0*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 0@.1Y..

      function ran1(idum)
c RANDOM NUMBER GENERATOR FROM NUMERICAL RECIPES, 2ND EDITION 
      implicit double precision (a-h,o-z)     
      integer idum,ia,im,iq,ir,ntab,ndiv

      parameter (ia=16807,im=2147483647,am=1.d0/im,iq=127773,ir=2836,
     * ntab=32,ndiv=1+(im-1)/ntab,eps=1.2d-7,rnmx=1.d0-eps)
      integer j,k,iv(ntab),iy
      save iv,iy
      data iv /ntab*0/, iy /0/

      if (idum.le.0.or.iy.eq.0) then
         idum=max(-idum,1)
         do j=ntab+8,1,-1
            k=idum/iq
            idum=ia*(idum-k*iq)-ir*k
            if (idum.lt.0) idum=idum+im
            if (j.le.ntab) iv(j)=idum
         enddo
         iy=iv(1)
      endif
      k=idum/iq
      idum=ia*(idum-k*iq)-ir*k
      if (idum.lt.0) idum=idum+im
      j=1+iy/ndiv
      iy=iv(j)
      iv(j)=idum
      ran1=min(am*iy,rnmx)
c      print *, ran1,idum*am
      return
      end

      subroutine cmbflat(clts,cltt,cltv,cles,clev,clbv,
     &           clet,clbt,clcs,clct,clcv,clbtv,clbtt, 
     &           clebv,clebt,trA,ak0)

c  This program integrates the linearized equations of general relativity,
c  the Boltzmann equations and the fluid equations for scalar perturbations
c  of a spatially flat Friedmann-Robertson-Walker universe in synchronous
c  gauge.  The time variable is conformal time dtau=dt/a(t) and the spatial
c  dependence is Fourier transformed.  ak is the comoving wavevector;
c  comoving distances are x=r/a(t), with a(t)=1 today.  The units of both
c  length and time are Mpc.
c
c  Presently restricted to a flat (K=0) universe.
c
c  All parameters are passed in common statements.
	implicit double precision(a-h,o-z)
        parameter (nnmax=1)
	integer  initfl
	common /initcase/ initfl
	common /lingerinc/ omegab,omegac,omegav,omegan,omegak,h0
     &                     ,tcmb,yhe,nnur,nnunr
	common /cosmoparm/ ak,ak2,amnu,lmaxg,lmaxnr,lmaxnu
     &                     ,nqmax,iq0,iq1,iq2
	common /genparm/ grhom,grhog,grhor,adotrad,taurst,dtaurec
        common /timesteps/ atau0,dtau1,dtau2
        common /reionization/zri,taurist,zristp,tauristp,rif
        common /reionization2/ j2ri1,nri,nri0

        common /initialps/ an(nnmax),nn
        common /tensor/ ant(nnmax),rat(nnmax),itflag,lmaxt
        common /vector/ lmaxv
        common /transfer/akmaxt,ztf,nlnkt,ict,ntf
	common /transfer2/ictpres

c
c  Dependent integration variables: a,phi,delta_c,theta_c,delta_b,theta_b,
c  2*(lmaxg+1) photon moments (total intensity plus polarization),
c  (lmaxnr+1) massless neutrino moments, (lmaxnu+1) massive neutrino
c  moments, nqmax momenta for each l.
c
c  For the tensor modes me will calculate a, ht, htprime, and
c  2*(lmaxt+1) photon moments (total intensity plus polarization).
c
c  For vector modes the integration variables are a, hv, vvb,
c  3*(lmaxv+1) photon moments (total intensity plus 2 polarization
c  states, (lmaxv+1) massless neutrino moments
c
c  The metric perturbation variables are computed algebraically.
        parameter (lmax0=8,lmaxnr0=25,lmaxnu0=25,nqmax0=15,lmaxt0=10)
        parameter (lmaxv0=8)
        parameter (lm1=8,lm2=7,lm3=4,nq1=15)
c lmx0 MUST be grater than all of the previous lmax0.
	parameter (lmx0=30)
        parameter (nk0=200)
        parameter (nstep0=6000)
c>>>>>>>>>>>>>>>>>>>>
c in old version used to add 1 to the number of scalar variables
c      parameter (nvar0=1+7+2*(lmax0+1)+(lmaxnr0+1)+nqmax0*(lmaxnu0+1))
c      common/maxvcom/maxvar
        parameter (nvar0=7+2*(lmax0+1)+(lmaxnr0+1)+nqmax0*(lmaxnu0+1))
c<<<<<<<<<<<<<<<<<<<<
    	parameter (nvar0t=2*(lmaxt0+1)+2+1)
    	parameter (nvar0v=2*lmaxv0+2*(lmaxv0-1)+3)
        parameter (nkmax=4000)
        parameter (l0max=3300,lmax=20+l0max/50,ketamax0=2*l0max+126)

c Transfer functions

        real*8 tautf(nnmax),ztf(nnmax)
c
c Scalar perturbations
      dimension trA(nk0)
	real*8 y(nvar0),yprime(nvar0)
        real*8 d(nk0,nstep0),dpr(nk0,nstep0)
        real*8 dp(nk0,nstep0),dppr(nk0,nstep0)
        real*8 dl(lmax),dl2(lmax),dl3(lmax)
        real*8 dpl2(lmax)
        real*8 s2(nstep0),sp2(nstep0)
        real*8 cl(lmax,nnmax),clpr(lmax,nnmax)
        real*8 cpl(lmax,nnmax),cplpr(lmax,nnmax)
        real*8 ccl(lmax,nnmax),cclpr(lmax,nnmax)
c
c Tensor Perturbations
	real*8 yt(nvar0t),ytprime(nvar0t)
        real*8 dt(nk0,nstep0),dtpr(nk0,nstep0)
        real*8 dte(nk0,nstep0),dtepr(nk0,nstep0)
        real*8 dtb(nk0,nstep0),dtbpr(nk0,nstep0)
        real*8 dtl(lmax),dtl2(lmax),dtl3(lmax)
        real*8 dtel2(lmax),dtbl2(lmax)
        real*8 st2(nstep0),ste2(nstep0)
        real*8 stb2(nstep0)
        real*8 ctl(lmax,nnmax),ctlpr(lmax,nnmax)
        real*8 ctel(lmax,nnmax),ctelpr(lmax,nnmax)
        real*8 ctbl(lmax,nnmax),ctblpr(lmax,nnmax)
        real*8 ctcl(lmax,nnmax),ctclpr(lmax,nnmax)
        real*8 ctbtl(lmax,nnmax),ctbtlpr(lmax,nnmax)
        real*8 ctebl(lmax,nnmax),cteblpr(lmax,nnmax)

c
c Vector perturbations
	  real*8 yv(nvar0v),yvprime(nvar0v)

        real*8 dv(nk0,nstep0),dvpr(nk0,nstep0)
        real*8 dve(nk0,nstep0),dvepr(nk0,nstep0)
        real*8 dvb(nk0,nstep0),dvbpr(nk0,nstep0)

        real*8 dvl(lmax),dvl2(lmax),dvl3(lmax)
        real*8 dvel(lmax),dvel2(lmax),dvel3(lmax)
        real*8 dvbl(lmax),dvbl2(lmax),dvbl3(lmax)

        real*8 sv2(nstep0),sve2(nstep0),svb2(nstep0)

        real*8 cvl(lmax,nnmax),cvlpr(lmax,nnmax)
        real*8 cvel(lmax,nnmax),cvelpr(lmax,nnmax)
        real*8 cvbl(lmax,nnmax),cvblpr(lmax,nnmax)

        real*8 cvcl(lmax,nnmax),cvclpr(lmax,nnmax)
        real*8 cvebl(lmax,nnmax),cveblpr(lmax,nnmax)
        real*8 cvbtl(lmax,nnmax),cvbtlpr(lmax,nnmax)

c
c Output arrays temperature: clts, cltt, cltv ; e spectra:
c cles, clet ; b perturbation(only tensor contrubution):
c clbt ; cross correlation: clcs, clct.

        real*8 clts(l0max,nnmax),cltt(l0max,nnmax)
        real*8 cltv(l0max,nnmax)
        real*8 cles(l0max,nnmax),clet(l0max,nnmax)
        real*8 clev(l0max,nnmax),clbv(l0max,nnmax)
        real*8 clbt(l0max,nnmax)
        real*8 clcs(l0max,nnmax),clct(l0max,nnmax)
        real*8 clebv(l0max,nnmax),clebt(l0max,nnmax)
        real*8 clcv(l0max,nnmax)
        real*8 clbtv(l0max,nnmax),clbtt(l0max,nnmax)

c
        real*8 atau0(nstep0),dtau1(nstep0),dtau2(nstep0)
        real*8 ak0(nk0)
        real*8 ajl(ketamax0,lmax),ajlpr(ketamax0,lmax)
        real*8 xx(ketamax0), dxx(ketamax0)
	integer mxx(nstep0)
	real*8 ak1(nkmax),dak1(nkmax)
        parameter (d0hi=1.0d40,d0lo=1.0d40)
        integer l(lmax)
	real*8 xl(lmax)

        common /lvalues1/ l,l0,lmo
        common /lvalues2/ akmax0
        save /lvalues1/
        save /lvalues2/

        common /jlgen/ ajl,ajlpr,xx,dxx,kmax
        save /jlgen/

c
	common /epsilon/ epsw
c
c
c  const=7*pi**4/120.
	parameter (const=5.68219698d0,zeta3=1.20205690d0)
c
        parameter (fourpi=4.0d0*3.14159265d0,xlimmin=10.0d0)
c  dverk  parameters.
	parameter (nw=nvar0,nwt=nvar0t,nwv=nvar0v)
        dimension c(24),w(nw,9),wt(nwt,9),wv(nwv,9)
c
c
	dimension denl(lmx0),dlfdlq(nqmax0)
	common /store/ denl,dlfdlq
      common/picom/pi,twopi3
      common/trfun/nkount
c
        double precision dtauda
	external fderivs,fderivst,fderivsv,dtauda

c zst is the value of z where the program stops the calculation.
        if ((abs(omegab+omegac-1.0d0).lt.1.0d-3).and.(zri.eq.0.0).
     & and.(itflag.eq.0).and.(h0.gt.40.0)) then
c>>>>>>>>>>>>>>>>>
cLP previously zst=10.0d0     
         zst=0.0d0
        else
         zst=0.0d0
        endif
        if (ict.ne.0) zst=min(zst,ztf(ntf))
c
	if (itflag.ne.2.and.itflag.ne.3) then
           lmaxg=lmax0
           if (ictpres.eq.1) then
              lmaxnr=lmaxnr0
           else
              lmaxnr=lm2
           end if
	end if
        nstep=nstep0
	if (itflag.ne.0.and.itflag.ne.3) then
	   lmaxt=lmaxt0
        else
           lmaxt=0
        end if
        if (itflag.ne.0.and.itflag.ne.2) then
	   lmaxv=lmaxv0
	   
        else
           lmaxv=0
        end if
c 

	do 133 j=1,lmx0
	  denl(j)=1.0d0/dble((2*j+1))
133	end do
c
c
c  Initialize neutrino mass and storage.
	if (nnunr.eq.0.or.omegan.eq.0.0) then
	  amnu=0.0d0
	  nqmax=0
          lmaxnu=0
	else
c  amnu=m_nu*c**2/(k_B*T_nu0).
c  The 1.0e-18 is to prevent certain compilers (e.g. Sun's) from
c  thinking that a divide by zero will occur.
	  amnu=const/(1.5d0*zeta3)*grhom/grhor*omegan
     &                            /dble((nnunr+1.0e-18))

          if (ictpres.eq.1) then
             nqmax=nqmax0
             lmaxnu=lmaxnu0
          else
             nqmax=nq1
             lmaxnu=lm3
          end if

	  dq=1.0d0	
	  do 134 i=1,nqmax
	    q=i*dq-0.5d0
	    expq=exp(-q)
	    dlfdlq(i)=-q/(1.0d0+expq)
134	  end do
	end if
c Calculate number of equations
	if (itflag.ne.2.and.itflag.ne.3) then
	   iq0=11+2*lmaxg+lmaxnr
	   iq1=iq0+nqmax
	   iq2=iq1+nqmax
c	   nvar=1+7+2*(lmaxg+1)+(lmaxnr+1)+nqmax*(lmaxnu+1)
	   nvar=7+2*(lmaxg+1)+(lmaxnr+1)+nqmax*(lmaxnu+1)
	else
	   nvar=0
	end if
	if (itflag.ne.0.and.itflag.ne.3) then
	   nvart=(2*lmaxt+5)
	else
	   nvart=0
	end if
	if (itflag.ne.0.and.itflag.ne.2) then
	   nvarv=(2*lmaxv+2*(lmaxv-1)+3)
	else
	   nvarv=0
	end if	
c  Initialize massive neutrinos.
	if ((nnunr.ne.0).or.(omegan.ne.0.0)) then
	   arel=1.d-3/amnu
	   call initnu1(amnu)
	end if

c  Time today
        tol1=1.d-6
        tau0=rombint(dtauda,0.0d0,1.0d0,tol1)
	epsw=100.0d0/tau0
c
c
c Maximum and minimum k-values. dlnk is the logarithmic spacing
c used for low k.
        akmax=akmax0/tau0
        akmin=0.15/tau0
	
        if (initfl.eq.4) then
c>>>>>>>>>>>>>>>>
cLP For strings I changed dlnk from 0.1 to 0.04
           dlnk=0.04d0
c<<<<<<<<<<<<<<<<
        else
           dlnk=0.1d0
        endif
c
c
c  Timesteps during recombination (tentative, the actual
c  timestep is the minimum between this value and taurst/40,
c  where taurst is the time when recombination starts.
c>>>>>>>>>>>>>>>>>>>>>>> canged dtaurec and dlntau0
	if (akmax.ne.0) then
	  if (initfl.eq.4) then
	    dtaurec=1.5d0/akmax
	   else
	    dtaurec=2.0d0/akmax
	   endif
	else
	   dtaurec=0.0d0
	end if
c  Stoping time
        taumax=rombint(dtauda,0.0d0,1.0d0/(zst+1.0d0),tol1)

c Timesteps after recombination, exponentially separated.
         if (initfl.eq.4) then
	   dlntau0=0.07d-1
	  else
	   dlntau0=0.1d-1
	  endif

	if (itflag.ne.0) dlntau0=0.5d0*dlntau0
	
c Time at which reionization takes place
	if (zri.ne.0.0) then
	   taurist=rombint(dtauda,0.0d0,1.0d0/(1+zri),tol1)
	   tauristp=rombint(dtauda,0.0d0,1.0d0/(1+zristp),tol1)
	else
	   taurist=tau0
	   tauristp=tau0
	end if
c
        taumin=.001d0/akmax
c
	taumin=min(taumin,0.1d0)

c

	if (amnu.ne.0.0d0) then
	   taumin=min(taumin,arel/adotrad)
	end if
c  Initialize baryon temperature and ionization fractions vs. time.
c  This subroutine also fixes the timesteps where the sources are
c  saved in order to do the integration.
        call finithermo(taumin,taumax,tau0,taurend,dlntau0,n1,nstep)

c
c  Calculating the times for the outputs of the transfer functions.
c
        if (ict.ne.0) then
           do 150 itf=1,ntf
              tautf(itf)=rombint(dtauda,0.0d0
     &                        ,1.0d0/(ztf(itf)+1.0d0),tol1)
              tautf(itf)=min(tautf(itf),atau0(nstep))
150       enddo
        endif

c
         nriend=0
c Integration will only be carried out after z=10 for low
c k, k < k10. If there is reionization the boundary will not be
c z=10 but tauristp
c	if (zst.ne.10.0d0) then
c	   t10=rombint(dtauda,0.0d0,1.0d0/11.0d0,tol1)
c	   ak10=250.0d0/tau0
c	   n10=n1+int(log(t10/taurend)/dlntau0)
	   if (zri.ne.0.0) then
	      t10=tau0
	      nriend=j2ri1+nri+nri0
	      t10=max(t10,atau0(nriend))
	      n10=nriend+int(log(t10/atau0(nriend))/dlntau0)
	   end if
c	else
	   t10=tau0
	   ak10=akmax
	   n10=nstep+1
c	end if

c
c  Calculation of the CMB sources.
c

       nk=0
       if (ict.ne.1) then
c     set k values for which the sources for the anisotropy and
c     polarization will be calculated. For low values of k we
c     use a logarithmic spacing.
c
	  if (zri.ne.0.0) then
	     dlnk0=2.0d0*dlnk
	  else
	     dlnk0=5.0d0*dlnk
	  end if

          dkn1=0.8d0/taurst
          dkn2=1.5d0/taurst
          nk1=int(log(dkn1/akmin/dlnk0)/dlnk0)+1
          if (akmax.gt.(1500.0d0/tau0)) then
             nk2=int((1500.0d0/tau0-akmin*exp((nk1-1)*dlnk0))
     &                /dkn1)+nk1+1
             nk=int((akmax-1500.0d0/tau0)/dkn2)+nk2+1
          else
             nk=int((akmax-akmin*exp((nk1-1)*dlnk0))/dkn1)+nk1+1
             nk2=nk+1
          end if
          if (nk.gt.nk0) then
             write(*,*)
     2            'Sorry, the arrays were dimensioned for a maximum of'
             write(*,*) nk0, 'k modes.'
             write(*,*)'The model you requested needs',nk
             write(*,*)'Please make the arrays bigger by making '
             write(*,*)'nk0 bigger where it appears'
             stop
          end if

          do 25 ik=1,nk
             if (ik.le.nk1) then
                ak0(ik)=akmin*exp((ik-1)*dlnk0)
             else
                if (ik.gt.nk2) then
                   ak0(ik)=ak0(nk2)+dble(ik-nk2)*dkn2
                else
                   ak0(ik)=ak0(nk1)+dble(ik-nk1)*dkn1
                end if
             end if
25       continue

c  Loop over wavenumbers.
          nkount=nk

          do 30 ik=1,nk	  

             ak=ak0(ik)
	     ak2=ak*ak
	     
c	   write(*,*)nstep,ik,ak	     
	     
c LP strings can have later start time, initfl=4 is for strings
          if(initfl.ne.4) then
	  
c  Begin when wave is far outside horizon.
c  Conformal time (in Mpc) in the radiation era, for photons plus 3 species
c  of relativistic neutrinos.
             taustart=0.001d0/ak
c  Make sure to start early in the radiation era. We are not including
c neutrinos as sources of the tensor modes.
             taustart=min(taustart,0.1d0)
             taustart0=taustart
c  Start when massive neutrinos are strongly relativistic.
             if (amnu.ne.0.0d0) then
                arel=1.d-3/amnu
                taustart=min(taustart0,arel/adotrad)
             end if

           else 

c            maxvar=nvar 
c LP in defect models one doesn't need to set conditions on superhorizon scales	    
	    taustart=0.01d0*atau0(2)
	    taustart0=taustart

c LP precalculate strings on an appropriate time grid	    
	    akx=ak
	    tau_i=0.1d0*taustart
	    tau_f=1.1d0*tau0
	    call strings_spline(akx,tau_i,tau_f)
	    
	   endif 
	    
             if (itflag.ne.2.and.itflag.ne.3) then
                call finitial(y,taustart)
                tau=taustart
                ind=1
c  Begin timestep loop.
c
c d contains the sources for the anisotropy and dp 
c for the polarization. t means tensor.

                     if (ak.lt.0.2d0) then 
                        tol1=1.d-5
                     else
                        tol1=1.d-4 
                        if(initfl.eq.4) tol1=1.d-2
                     endif

                if (ict.ne.0) itf=1

c                write(*,*)'s',tol1
                do 10 j=2,nstep								  

                   tauend=atau0(j)
c                    
                   if (ak.gt.ak10.and.tauend.gt.t10.and.
     2(ict.eq.0.or.tau.gt.tautf(ntf))) then
c     
                      d(ik,j)=0.0d0
                      dp(ik,j)=0.0d0
                   else

                      call dverk(nvar,fderivs,tau,y,tauend,tol1
     &                     ,ind,c,nw,w)
                      call fderivs(nvar,tau,y,yprime)

                      call foutput(nvar,y,yprime,j,tau0,tau
     &                                       ,d(ik,j),dp(ik,j))

c     Calculation of transfer functions.
101                  continue
                      if (ict.ne.0.and.itf.le.ntf) then 
                         if (j.lt.nstep.and.tauend.lt.tautf(itf).
     2                        and.atau0(j+1).gt.tautf(itf)) then
                            call dverk(nvar,fderivs,tau,y,tautf(itf)
     2                           ,tol1,ind,c,nw,w)
                         endif
c output transfer functions for this k-values.                     
                         if (abs(tau-tautf(itf)).lt.1.0d-5) then
cLP   delta_CDM**2/k**4 is saved.
c     it will get multiplied by k**4 and the proper normalization
c     in subroutine COBEnormalize
                            trA(ik)=(y(4)/ak2)**2
                            itf=itf+1
                            if (j.lt.nstep.and.itf.le.ntf.and.
     2                           atau0(j+1).gt.tautf(itf)) goto 101
                         endif
                      endif
                   end if
10             end do
                d(ik,nstep)=0.0d0
                dp(ik,nstep)=0.0d0
             end if
          indic=1
	  

c time loop for tensors if requested.
c
c Tensors will only be calculated if k*tau< stsp.
c>>>>>>>>>>>>>>>>>>>>>>>
cLP I changed stpt from 50 to 1.e8
c
             stpt=1.0d8
c
c
	if (itflag.ne.0.and.itflag.ne.3) then
	   call finitialt(yt,taustart0)
	   tau=taustart0
	   ind=1

c  Begin timestep loop.
c dt contains the sources for the anisotropy and dtp
c for the polarization. t means tensor.
               if (ak.lt.0.2) then 
                tol1=1.0d-5
               else
                tol1=1.0d-4
                if(initfl.eq.4) tol1=1.d-3
               endif
c                write(*,*)'t',tol1
	   do 11 j=2,nstep-1	   
	    
		    tauend=atau0(j)
	      if ((ak*tauend).gt.stpt) then
		 dt(ik,j)=0.0d0
		 dte(ik,j)=0.0d0
		 dtb(ik,j)=0.0d0
	      else
	       call dverk(nvart,fderivst,tau,yt,tauend,tol1,
     &                                          ind,c,nwt,wt)
               call fderivst(nvart,tau,yt,ytprime)
	       call foutputt(nvart,yt,ytprime,j,tau0,tau,
     &                           dt(ik,j),dte(ik,j),dtb(ik,j))
	      end if
11	   end do
           dt(ik,nstep)=0.0d0
           dte(ik,nstep)=0.0d0
           dtb(ik,nstep)=0.0d0
	end if
c
c>>>>>>>>>>>>>>>>>>>
cLP  start the time loop for vectors
	if (itflag.ne.0.and.itflag.ne.2) then
	   call finitialv(yv,taustart0)
	   tau=taustart0
	   ind=1
                     if (ak.lt.0.2d0) then 
                        tol1=1.d-4
                     else
                        tol1=1.d-4
                        if(initfl.eq.4) tol1=1.d-2
                     endif
c  Begin timestep loop.
c dv contains the sources for the anisotropy 
c v means vector.
c                write(*,*)'v',tol1
	   do 911 j=2,nstep-1	   
	   
	      tauend=atau0(j)
	      if ((ak*tauend).gt.1.d8) then
		 dv(ik,j)=0.0d0
                 dve(ik,j)=0.0d0
                 dvb(ik,j)=0.0d0
	      else

	       call dverk(nvarv,fderivsv,tau,yv,tauend,tol1,
     &                                          ind,c,nwv,wv)
          call fderivsv(nvarv,tau,yv,yvprime)
	  call foutputv(nvarv,yv,yvprime,j,tau0,tau,
     &                        dv(ik,j),dve(ik,j),dvb(ik,j))
c
	      end if
911	   end do
           dv(ik,nstep)=0.0d0
           dve(ik,nstep)=0.0d0
           dvb(ik,nstep)=0.0d0

	end if	
c
c
c
30     end do
      endif
 
c
c if CMB calculations are requested, calculate the Cl by 
c integrating the sources over time and over k.
c
c
	if (ict.ne.1) then
c get the interpolation matrix for the sources to interpolate them
c for other k-values, scalar case.
        if (itflag.ne.2.and.itflag.ne.3) then	    
	   do 100 i=1,nstep
	      call spline(ak0,d(1,i),nk,d0lo,d0hi,dpr(1,i))
	      call spline(ak0,dp(1,i),nk,d0lo,d0hi,dppr(1,i))
100	   continue
	   do 105 in=1,nn
	      do 106 j=1,l0
		 cl(j,in)=0.0d0
		 cpl(j,in)=0.0d0
		 ccl(j,in)=0.0d0
106	      continue
105	   continue
	end if
c
c get the interpolation matrix for the tensor sources.
	if (itflag.ne.0.and.itflag.ne.3) then
	   do 110 i=1,nstep
	      call spline(ak0,dt(1,i),nk,d0lo,d0hi,dtpr(1,i))
	      call spline(ak0,dte(1,i),nk,d0lo,d0hi,dtepr(1,i))
	      call spline(ak0,dtb(1,i),nk,d0lo,d0hi,dtbpr(1,i))
110	   continue
	   do 115 in=1,nn
	      do 116 j=1,l0
		 ctl(j,in)=0.0d0
		 ctel(j,in)=0.0d0
		 ctbl(j,in)=0.0d0
		 ctcl(j,in)=0.0d0
                 ctbtl(j,in)=0.0d0
		 ctebl(j,in)=0.0d0
116	      continue
115	   continue
	end if
c
c get the interpolation matrix for the vector sources.
	if (itflag.ne.0.and.itflag.ne.2) then
	   do 120 i=1,nstep
	      call spline(ak0,dv(1,i),nk,d0lo,d0hi,dvpr(1,i))
              call spline(ak0,dve(1,i),nk,d0lo,d0hi,dvepr(1,i))
              call spline(ak0,dvb(1,i),nk,d0lo,d0hi,dvbpr(1,i))
120	   continue
	   do 125 in=1,nn
	      do 126 j=1,l0
		 cvl(j,in)=0.0d0
                 cvel(j,in)=0.0d0
                 cvbl(j,in)=0.0d0
                 cvcl(j,in)=0.0d0
                 cvbtl(j,in)=0.0d0
                 cvebl(j,in)=0.0d0
126	      continue
125	   continue
	end if
c 
c Fixing the # of k for the integration.
	dk=2.5d0/tau0
	dk0=1.5d0/tau0

	no=500
        dlnk1=0.1d0
	no1=int(log(10.0d0*dk0/akmin)/dlnk1)+1
	if (akmax.gt.(no*dk0)) then 
	   nko=int((akmax-no*dk0)/dk)+no
	else
	   no=int((akmax-10.0d0*dk0)/dk0)+no1
	   nko=no
	end if

        if (nko.gt.nkmax) then
           write(*,*)
     2          'Sorry, the arrays were dimensioned for a maximum of'
           write(*,*) nkmax, 'k modes.'
           write(*,*)'The model you requested needs',nk 
           write(*,*)'Please make the arrays bigger by making '
           write(*,*)'nkmax bigger where it appears'
           stop
        end if
	do 198 k=1,nko
	   if (k.le.no) then
	      if (k.le.no1) then
		 ak1(k)=10.0d0*dk0*exp(-(no1-k)*dlnk1)
                 dak1(k)=ak1(k)*dlnk1
	      else
		 ak1(k)=ak1(no1)+(k-no1)*dk0
                 dak1(k)=dk0
              end if
           else
              ak1(k)=ak1(no)+(k-no)*dk
	      dak1(k)=dk
           end if
198	continue
	do 199 k=2,(nko-1)
199	continue
	dak1(1)=0.5d0*dak1(1)
	dak1(no1)=0.5d0*(dak1(no1)+dk0)
	dak1(nko)=0.5d0*dak1(nko)
c
	klo=1
	khi=2
c Begin k-loop
       do 200 k=1,nko

          akt=ak1(k)
c finding position of k in table ak0 to do the interpolation.
211	  continue
	  if ((akt.gt.ak0(klo+1)).and.(klo.lt.(nk-1))) then
	     klo=klo+1
	     khi=klo+1
	     goto 211
	  end if      
	  ho=ak0(khi)-ak0(klo)
	  a0=(ak0(khi)-akt)/ho
	  b0=(akt-ak0(klo))/ho

	  nstps=0
	  nstpt=0
	  nstpv=0

	  if (itflag.ne.2.and.itflag.ne.3) then
	     do 250 j=1,l0
		dl(j)=0.0d0
		dl2(j)=0.0d0
		dl3(j)=0.0d0
		dpl2(j)=0.0d0
250	     continue
c
c  Interpolating the source as a function of time for the present
c  wavelength.
	     if (akt.lt.ak10) then
		nstps=nstep-1
	     else
		nstps=n10
	     end if
		   
	     s2(1)=0.0d0
	     sp2(1)=0.0d0
	     do 304 i=2,nstps
		s2(i)=a0*d(klo,i)+b0*d(khi,i)+((a0**3-a0)*dpr(klo,i)
     &                   +(b0**3-b0)*dpr(khi,i))*ho*ho/6.d0 
		sp2(i)=a0*dp(klo,i)+b0*dp(khi,i)+((a0**3-a0)*dppr(klo,i)
     &                   +(b0**3-b0)*dppr(khi,i))*ho*ho/6.d0
304	     continue
	     s2(nstps+1)=0.0d0
	     sp2(nstps+1)=0.0d0
	  end if

c If tensors wanted
	 if (itflag.ne.0.and.itflag.ne.3) then
	    do 260 j=1,l0
	       dtl(j)=0.0d0
	       dtl2(j)=0.0d0
	       dtl3(j)=0.0d0
	       dtel2(j)=0.0d0
	       dtbl2(j)=0.0d0
260	    continue
c
c  Interpolating the tensor source as a function of time for the present
c  wavelength.
	    st2(1)=0.0d0
	    ste2(1)=0.0d0
	    stb2(1)=0.0d0
	    nstpt=2
	    do 306 i=2,(nstep-2)
	       xf=akt*(tau0-atau0(i))
	       if (((akt*atau0(i)).lt.stpt).and.(xf.gt.1.0d-8)) then
		  nstpt=i
		  st2(i)=a0*dt(klo,i)+b0*dt(khi,i)+((a0**3-a0)
     &                *dtpr(klo,i)+(b0**3-b0)*dtpr(khi,i))*ho*ho/6.d0 
		  ste2(i)=a0*dte(klo,i)+b0*dte(khi,i)+((a0**3-a0)
     &               *dtepr(klo,i)+(b0**3-b0)*dtepr(khi,i))*ho*ho/6.d0
		  stb2(i)=a0*dtb(klo,i)+b0*dtb(khi,i)+((a0**3-a0)
     &               *dtbpr(klo,i)+(b0**3-b0)*dtbpr(khi,i))*ho*ho/6.d0
	       else
		  st2(i)=0.0d0
		  ste2(i)=0.0d0
		  stb2(i)=0.0d0
	       end if
306	    end do
	    nstpt=max(nstpt,n1)
	    st2(nstpt+1)=0.0d0
	    ste2(nstpt+1)=0.0d0
	    stb2(nstpt+1)=0.0d0
	 end if

c>>>>>>>>>>>>>>>>>>>>>
c If vectors wanted
	 if (itflag.ne.0.and.itflag.ne.2) then
	    do 360 j=1,l0
	       dvl(j)=0.0d0
	       dvl2(j)=0.0d0
	       dvl3(j)=0.0d0
c
               dvel(j)=0.0d0
	       dvel2(j)=0.0d0
	       dvel3(j)=0.0d0
c
               dvbl(j)=0.0d0
	       dvbl2(j)=0.0d0
	       dvbl3(j)=0.0d0

360	    continue
c
c  Interpolating the tensor source as a function of time for the present
c  wavelength.
	    sv2(1)=0.0d0
	    nstpv=2
	    do 406 i=2,(nstep-2)
	       xf=akt*(tau0-atau0(i))
	       if (((akt*atau0(i)).lt.stpt).and.(xf.gt.1.0d-8)) then
		  nstpv=i
		  sv2(i)=a0*dv(klo,i)+b0*dv(khi,i)+((a0**3-a0)
     &              *dvpr(klo,i)+(b0**3-b0)*dvpr(khi,i))*ho*ho/6.d0
                  sve2(i)=a0*dve(klo,i)+b0*dve(khi,i)+((a0**3-a0)
     &              *dvepr(klo,i)+(b0**3-b0)*dvepr(khi,i))*ho*ho/6.d0
                  svb2(i)=a0*dvb(klo,i)+b0*dvb(khi,i)+((a0**3-a0)
     &              *dvbpr(klo,i)+(b0**3-b0)*dvbpr(khi,i))*ho*ho/6.d0
	       else
		  sv2(i)=0.0d0
                  sve2(i)=0.0d0
                  svb2(i)=0.0d0
	       end if
406	    end do
	    nstpv=max(nstpv,n1)
	    sv2(nstpv+1)=0.0d0
            sve2(nstpv+1)=0.0d0
            svb2(nstpv+1)=0.0d0
	 end if
c
c
         nstpt=max(nstpv,nstpt)

c Findind the position in the xx table for the x correponding to each
c timestep
	  do 405 i=1,max(nstps,nstpt)+1
	     xf=abs(akt*(tau0-atau0(i)))
	     if (xf.le.25.0) then
		if (xf.le.5.0) then
		   de2=10.0d0*xf+1.0d0
		else
		   de2=(xf-5.0d0)*5.0d0+51.0d0
		end if
	     else
		de2=(xf-25.0d0)+151.0d0
	     end if
	     mxx(i)=int(de2)
405	  end do

c
c Begin l and  time-loop to integrate scalar perturbations.
c
c Determining ranges of integration

	 if (itflag.ne.2.and.itflag.ne.3) then
	    do 310 j=1,l0
	       xlim=0.05d0*l(j)
	       xlim=max(xlim,xlimmin)
	       xlim=l(j)-xlim
	       xlmax1=80.0d0*l(j)
	       xlmax2=dble(min(2*l(j),kmax))
	       tmin=tau0-xlmax1/akt
	       tmax=tau0-xlim/akt
               tmax=min(tau0,tmax)
               tmin=max(atau0(2),tmin)
c
	       if (tmax.lt.atau0(2)) goto 400

	       if (zri.eq.0.0) then

		  if (tmin.lt.taurend) then
		     nstart1=2
		  else
		     nstart1=n1+int(log(tmin/taurend)/dlntau0)
		  end if
		  if (tmax.lt.taurend) then
		     nstop1=n1
		     nstop1a=n1
		  else
		     nstop1=n1+int(log(tmax/taurend)/dlntau0)
		     nstop1=min(nstop1,nstps)
		     if ((akt*dtau1(nstop1)).gt.1) then
			nstop1a=max(n1,n1-int(log(akt*dlntau0*taurend)
     &                          /dlntau0))
		     else
			nstop1a=nstop1
		     end if
		  end if

	       else

		  if (tmin.lt.taurend) then
		     nstart1=2
		     nstart2=j2ri1
		  else
		     if (tmin.lt.atau0(nriend)) then
			nstart1=min(n1+int(log(tmin/taurend)/dlntau0),
     &                    j2ri1)
			nstart2=j2ri1
		     else
			nstart1=nstep+1
			nstart2=nriend+int(log(tmin/atau0(nriend)
     &                      /dlntau0))
		     end if
		  end if

		  if (tmax.lt.taurend) then
		     nstop1=n1
		     nstop2=0
		  else
		     if (tmax.lt.atau0(j2ri1)) then
			nstop1=n1+int(log(tmax/taurend)/dlntau0)
			nstop2=0
		     else
			nstop1=j2ri1-1
			nstop2=max(nriend,nriend+int(log(tmax
     &                                  /atau0(nriend))/dlntau0))

			nstop2=min(nstop2,nstps)
			if ((akt*dtau1(nstop2)).gt.1) then
			   nstop2a=max(nriend,nriend-int(log(akt
     &                        *dlntau0*atau0(nriend))/dlntau0))
			else
			   nstop2a=nstop2
			end if

		     end if

		     nstop1=min(nstop1,j2ri1-1)
		     if ((akt*dtau1(nstop1)).gt.1) then
			nstop1a=max(n1,n1-int(log(akt
     &                             *dlntau0*taurend)/dlntau0))
		     else
			nstop1a=nstop1
		     end if

		  end if
	       end if

c Integration before reionization.
c
c Interpolating jls at points where the
c sources are recorded.
	       do i=nstart1,nstop1a
		  xf=akt*(tau0-atau0(i))
		  m2=mxx(i)
		  h2=xx(m2+1)-xx(m2)
		  a2=(xx(m2+1)-xf)/h2
		  b2=(xf-xx(m2))/h2
		  ajl0=a2*ajl(m2,j)+b2*ajl(m2+1,j)+((a2**3-a2)
     &     *ajlpr(m2,j)+(b2**3-b2)*ajlpr(m2+1,j))*(h2*h2)/6.d0
		  dl2(j)=dl2(j)+s2(i)*ajl0*dtau2(i)
		  dpl2(j)=dpl2(j)+sp2(i)*ajl0*dtau2(i)
	       end do

c Interpolating sources at points where the
c jls are recorded.

	       do i=nstop1a+1,nstop1
		  xf=akt*(tau0-atau0(i))
		  m2=mxx(i)
		  dtau3=dtau1(i)*abs(s2(i)/(s2(i+1)-s2(i)+1.0d-10))
		  if ((xf.lt.xlmax2).or.((akt*dtau3).lt.1.0)) then
		     xi=xf-akt*dtau1(i)
		     m1=mxx(i+1)
		     ddt1=s2(i)
		     ddt2=s2(i+1)
		     do lx=m1+1,m2
			x=xx(lx)
			ddt=(ddt1-ddt2)*(x-xi)/(xf-xi)+ddt2
			dl3(j)=dl3(j)+ajl(lx,j)*ddt*dxx(lx)
		     end do
		  end if
	       end do

c Integration after reionization
c Interpolating jls at points where the
c sources are recorded.

	       if (zri.ne.0.0) then
		  do i=nstart2,nstop2a
		     xf=akt*(tau0-atau0(i))
		     m2=mxx(i)
		     h2=xx(m2+1)-xx(m2)
		     a2=(xx(m2+1)-xf)/h2
		     b2=(xf-xx(m2))/h2
		     ajl0=a2*ajl(m2,j)+b2*ajl(m2+1,j)+((a2**3-a2)*
     &  	         ajlpr(m2,j)+(b2**3-b2)*ajlpr(m2+1,j))
     &                                      *(h2*h2)/6.d0
		     dl2(j)=dl2(j)+s2(i)*ajl0*dtau2(i)
		     dpl2(j)=dpl2(j)+sp2(i)*ajl0*dtau2(i)
		  end do

c Interpolating sources at points where the
c jls are recorded.

		  do i=nstop2a+1,nstop2
		     xi=xf-akt*dtau1(i)
		     xf=akt*(tau0-atau0(i))
		     m2=mxx(i)
		     dtau3=dtau1(i)*abs(s2(i)/(s2(i+1)-s2(i)+1.0d-10))
		     if ((xf.lt.xlmax2).or.((akt*dtau3).lt.1)) then
			m1=mxx(i+1)
			ddt1=s2(i)
			ddt2=s2(i+1)
			do lx=m1+1,m2
			   x=xx(lx)
			   ddt=(ddt1-ddt2)*(x-xi)/(xf-xi)+ddt2
			   dl3(j)=dl3(j)+ajl(lx,j)*ddt*dxx(lx)
			end do
		     end if
		  end do
	       end if
	    
	       dl(j)=dl2(j)+dl3(j)/akt
310	    continue

400	    continue
c

c Adding to calculate the integral over k.
c Scalar case
c
	    do in=1,nn
	       do j=1,l0
                  call powersflat(akt,in,apowers)
		  ckj=apowers*dl(j)*dl(j)*dak1(k)
		  cl(j,in)=cl(j,in)+ckj
		  cpkj=apowers*dpl2(j)*dpl2(j)*dak1(k)
		  cpl(j,in)=cpl(j,in)+cpkj
		  cckj=apowers*dl(j)*dpl2(j)*dak1(k)
		  ccl(j,in)=ccl(j,in)+cckj
	       end do
	    end do
c
	 end if
c
c Begin l and  time-loop to integrate tensor perturbations.
c
c Finding the ranges of integration.

	 if (itflag.ne.0.and.itflag.ne.3) then
	    do 530 j=1,l0
	       xlim=0.05d0*l(j)
	       xlim=max(xlim,xlimmin)
	       xlim=l(j)-xlim
	       xlmax1=80.0d0*l(j)
c	       xlmax2=dble(min(2*l(j),kmax))
	       tmin=tau0-xlmax1/akt
	       tmax=tau0-xlim/akt
               tmax=min(tau0,tmax)
               tmin=max(atau0(2),tmin)

c
	       if (tmax.lt.atau0(2)) goto 590

	       if (zri.eq.0.0) then

		  if (tmin.lt.taurend) then
		     nstart1=2
		  else
		     nstart1=n1+int(log(tmin/taurend)/dlntau0)
		  end if
		  if (tmax.lt.taurend) then
		     nstop1=n1
		     nstop1a=n1
		  else
		     nstop1=n1+int(log(tmax/taurend)/dlntau0)
		     nstop1=min(nstop1,nstpt)
		     if ((akt*dtau1(nstop1)).gt.1) then
			nstop1a=max(n1,n1-int(log(akt*dlntau0*taurend)
     &                          /dlntau0))
		     else
			nstop1a=nstop1
		     end if
		  end if

	       else

		  if (tmin.lt.taurend) then
		     nstart1=2
		     nstart2=j2ri1
		  else
		     if (tmin.lt.atau0(nriend)) then
			nstart1=min(n1+int(log(tmin/taurend)/dlntau0),
     &                    j2ri1)
			nstart2=j2ri1
		     else
			nstart1=nstep+1
			nstart2=nriend+int(log(tmin/atau0(nriend)
     &                      /dlntau0))
		     end if
		  end if

		  if (tmax.lt.taurend) then
		     nstop1=n1
		     nstop2=0
		  else
		     if (tmax.lt.atau0(j2ri1)) then
			nstop1=n1+int(log(tmax/taurend)/dlntau0)
			nstop2=0
		     else
			nstop1=j2ri1-1
			nstop2=max(nriend,nriend+int(log(tmax
     &                                  /atau0(nriend))/dlntau0))

			nstop2=min(nstop2,nstpt)    
			if ((akt*dtau1(nstop2)).gt.1) then
			   nstop2a=max(nriend,nriend-int(log(akt
     &                        *dlntau0*atau0(nriend))/dlntau0))
			else
			   nstop2a=nstop2
			end if

		     end if

		     nstop1=min(nstop1,j2ri1-1)
		     if ((akt*dtau1(nstop1)).gt.1) then
			nstop1a=max(n1,n1-int(log(akt
     &                             *dlntau0*taurend)/dlntau0))
		     else
			nstop1a=nstop1
		     end if

		  end if
	       end if

c Integration before reionization.
c
c Interpolating jls at points where the
c sources are recorded.

	       do i=nstart1,nstop1a
		  xf=akt*(tau0-atau0(i))
		  m2=mxx(i)
		  h2=xx(m2+1)-xx(m2)
		  a2=(xx(m2+1)-xf)/h2
		  b2=(xf-xx(m2))/h2
		  ajl0=a2*ajl(m2,j)+b2*ajl(m2+1,j)+((a2**3-a2)
     &                *ajlpr(m2,j)+(b2**3-b2)*ajlpr(m2+1,j))
     &                                        *(h2*h2)/6.d0


		  dtl2(j)=dtl2(j)+st2(i)*ajl0*dtau2(i)
		  dtel2(j)=dtel2(j)+ste2(i)*ajl0*dtau2(i)
		  dtbl2(j)=dtbl2(j)+stb2(i)*ajl0*dtau2(i)
	       end do

c Interpolating sources at points where the
c jls are recorded.

	       do i=nstop1a+1,nstop1
		  xf=akt*(tau0-atau0(i))
		  m2=mxx(i)
		  xi=xf-akt*dtau1(i)
		  m1=mxx(i+1)
		  dtdt1=st2(i)
		  dtdt2=st2(i+1)
		  do lx=m1+2,m2+1
		     x=xx(lx)
		     dtdt=(dtdt1-dtdt2)*(x-xi)/(xf-xi)+dtdt2
		     dtl3(j)=dtl3(j)+ajl(lx,j)*dtdt*dxx(lx)
 		  end do
	       end do

c Integration after reionization
	       if (zri.ne.0.0) then
c Interpolating jls at points where the
c sources are recorded.
c
		  do i=nstart2,nstop2a
		     xf=akt*(tau0-atau0(i))
		     m2=mxx(i)
		     h2=xx(m2+1)-xx(m2)
		     a2=(xx(m2+1)-xf)/h2
		     b2=(xf-xx(m2))/h2
		     ajl0=a2*ajl(m2,j)+b2*ajl(m2+1,j)+((a2**3-a2)*
     &  	          ajlpr(m2,j)+(b2**3-b2)
     &                    *ajlpr(m2+1,j))*(h2*h2)/6.d0
		     dtl2(j)=dtl2(j)+st2(i)*ajl0*dtau2(i)
		     dtel2(j)=dtel2(j)+ste2(i)*ajl0*dtau2(i)
		     dtbl2(j)=dtbl2(j)+stb2(i)*ajl0*dtau2(i)
		  end do

c Interpolating sources at points where the
c jls are recorded.

		  do i=nstop2a+1,nstop2
		     xi=xf-akt*dtau1(i)
		     xf=akt*(tau0-atau0(i))
		     m2=mxx(i)
		     m1=mxx(i+1)
		     dtdt1=st2(i)
		     dtdt2=st2(i+1)
		     do lx=m1+2,m2+1
			x=xx(lx)
			dtdt=(dtdt1-dtdt2)*(x-xi)/(xf-xi)+dtdt2
			dtl3(j)=dtl3(j)+ajl(lx,j)*dtdt*dxx(lx)
		     end do
		  end do
	       end if

	       dtl(j)=dtl2(j)+dtl3(j)/akt
530	    continue

590	    continue


c Tensor case.
c
c Adding to calculate the integral over k.
c Tensor case
	    do in=1,nn
	       do j=1,l0
                  call powertflat(akt,in,apowert)

		  ctkj=apowert*dtl(j)*dtl(j)*dak1(k)
		  ctl(j,in)=ctl(j,in)+ctkj

		  ctekj=apowert*dtel2(j)*dtel2(j)*dak1(k)
		  ctel(j,in)=ctel(j,in)+ctekj

		  ctbkj=apowert*dtbl2(j)*dtbl2(j)*dak1(k)
		  ctbl(j,in)=ctbl(j,in)+ctbkj

		  ctckj=apowert*dtl(j)*dtel2(j)*dak1(k)
		  ctcl(j,in)=ctcl(j,in)+ctckj

                  ctbtkj=apowert*dtbl2(j)*dtl(j)*dak1(k)
		  ctbtl(j,in)=ctbtl(j,in)+ctbtkj

                  ctebkj=apowert*dtel2(j)*dtbl2(j)*dak1(k)
		  ctebl(j,in)=ctebl(j,in)+ctebkj

	       end do
	    end do
c
	 end if
c
c Begin l and  time-loop to integrate vector perturbations.
c
c Finding the ranges of integration.

	 if (itflag.ne.0.and.itflag.ne.2) then
	    do 730 j=1,l0
	       xlim=0.05d0*l(j)
	       xlim=max(xlim,xlimmin)
	       xlim=l(j)-xlim
	       xlmax1=80.0d0*l(j)
c	       xlmax2=dble(min(2*l(j),kmax))
	       tmin=tau0-xlmax1/akt
	       tmax=tau0-xlim/akt
               tmax=min(tau0,tmax)
               tmin=max(atau0(2),tmin)

c
	       if (tmax.lt.atau0(2)) goto 790

	       if (zri.eq.0.0) then

		  if (tmin.lt.taurend) then
		     nstart1=2
		  else
		     nstart1=n1+int(log(tmin/taurend)/dlntau0)
		  end if
		  if (tmax.lt.taurend) then
		     nstop1=n1
		     nstop1a=n1
		  else
		     nstop1=n1+int(log(tmax/taurend)/dlntau0)
		     nstop1=min(nstop1,nstpv)
		     if ((akt*dtau1(nstop1)).gt.1) then
			nstop1a=max(n1,n1-int(log(akt*dlntau0*taurend)
     &                          /dlntau0))
		     else
			nstop1a=nstop1
		     end if
		  end if

	       else

		  if (tmin.lt.taurend) then
		     nstart1=2
		     nstart2=j2ri1
		  else
		     if (tmin.lt.atau0(nriend)) then
			nstart1=min(n1+int(log(tmin/taurend)/dlntau0),
     &                    j2ri1)
			nstart2=j2ri1
		     else
			nstart1=nstep+1
			nstart2=nriend+int(log(tmin/atau0(nriend)
     &                      /dlntau0))
		     end if
		  end if

		  if (tmax.lt.taurend) then
		     nstop1=n1
		     nstop2=0
		  else
		     if (tmax.lt.atau0(j2ri1)) then
			nstop1=n1+int(log(tmax/taurend)/dlntau0)
			nstop2=0
		     else
			nstop1=j2ri1-1
			nstop2=max(nriend,nriend+int(log(tmax
     &                                  /atau0(nriend))/dlntau0))

			nstop2=min(nstop2,nstpv)
			if ((akt*dtau1(nstop2)).gt.1) then
			   nstop2a=max(nriend,nriend-int(log(akt
     &                        *dlntau0*atau0(nriend))/dlntau0))
			else
			   nstop2a=nstop2
			end if

		     end if

		     nstop1=min(nstop1,j2ri1-1)
		     if ((akt*dtau1(nstop1)).gt.1) then
			nstop1a=max(n1,n1-int(log(akt
     &                             *dlntau0*taurend)/dlntau0))
		     else
			nstop1a=nstop1
		     end if

		  end if
	       end if

c Integration before reionization.
c
c Interpolating jls at points where the
c sources are recorded.

	       do i=nstart1,nstop1a
		  xf=akt*(tau0-atau0(i))
		  m2=mxx(i)
		  h2=xx(m2+1)-xx(m2)
		  a2=(xx(m2+1)-xf)/h2
		  b2=(xf-xx(m2))/h2
		  ajl0=a2*ajl(m2,j)+b2*ajl(m2+1,j)+((a2**3-a2)
     &                *ajlpr(m2,j)+(b2**3-b2)*ajlpr(m2+1,j))
     &                                        *(h2*h2)/6.d0


		  dvl2(j)=dvl2(j)+sv2(i)*ajl0*dtau2(i)
                  dvel2(j)=dvel2(j)+sve2(i)*ajl0*dtau2(i)
                  dvbl2(j)=dvbl2(j)+svb2(i)*ajl0*dtau2(i)
	       end do

c Interpolating sources at points where the
c jls are recorded.

	       do i=nstop1a+1,nstop1
		  xf=akt*(tau0-atau0(i))
		  m2=mxx(i)
		  xi=xf-akt*dtau1(i)
		  m1=mxx(i+1)
		  dtdt1=sv2(i)
		  dtdt2=sv2(i+1)

                  dtdt1e=sve2(i)
		  dtdt2e=sve2(i+1)

                  dtdt1b=svb2(i)
		  dtdt2b=svb2(i+1)

		  do lx=m1+2,m2+1
		     x=xx(lx)

		     dtdt=(dtdt1-dtdt2)*(x-xi)/(xf-xi)+dtdt2
		     dvl3(j)=dvl3(j)+ajl(lx,j)*dtdt*dxx(lx)

                     dtdte=(dtdt1e-dtdt2e)*(x-xi)/(xf-xi)+dtdt2e
		     dvel3(j)=dvel3(j)+ajl(lx,j)*dtdte*dxx(lx)

                     dtdtb=(dtdt1b-dtdt2b)*(x-xi)/(xf-xi)+dtdt2b
		     dvbl3(j)=dvbl3(j)+ajl(lx,j)*dtdtb*dxx(lx)
 		  end do
	       end do

c Integration after reionization
	       if (zri.ne.0.0) then
c Interpolating jls at points where the
c sources are recorded.
c
		  do i=nstart2,nstop2a
		     xf=akt*(tau0-atau0(i))
		     m2=mxx(i)
		     h2=xx(m2+1)-xx(m2)
		     a2=(xx(m2+1)-xf)/h2
		     b2=(xf-xx(m2))/h2
		     ajl0=a2*ajl(m2,j)+b2*ajl(m2+1,j)+((a2**3-a2)*
     &  	          ajlpr(m2,j)+(b2**3-b2)
     &                    *ajlpr(m2+1,j))*(h2*h2)/6.d0
		     dvl2(j)=dvl2(j)+sv2(i)*ajl0*dtau2(i)
                     dvel2(j)=dvel2(j)+sve2(i)*ajl0*dtau2(i)
                     dvbl2(j)=dvbl2(j)+svb2(i)*ajl0*dtau2(i)
		  end do

c Interpolating sources at points where the
c jls are recorded.

		  do i=nstop2a+1,nstop2
		     xi=xf-akt*dtau1(i)
		     xf=akt*(tau0-atau0(i))
		     m2=mxx(i)
		     m1=mxx(i+1)

                     dtdt1=sv2(i)
		     dtdt2=sv2(i+1)

                     dtdt1e=sve2(i)
		     dtdt2e=sve2(i+1)

                     dtdt1b=svb2(i)
		     dtdt2b=svb2(i+1)

		     do lx=m1+2,m2+1
			x=xx(lx)
			dtdt=(dtdt1-dtdt2)*(x-xi)/(xf-xi)+dtdt2
			dvl3(j)=dvl3(j)+ajl(lx,j)*dtdt*dxx(lx)

                     dtdte=(dtdt1e-dtdt2e)*(x-xi)/(xf-xi)+dtdt2e
		     dvel3(j)=dvel3(j)+ajl(lx,j)*dtdte*dxx(lx)

                     dtdtb=(dtdt1b-dtdt2b)*(x-xi)/(xf-xi)+dtdt2b
		     dvbl3(j)=dvbl3(j)+ajl(lx,j)*dtdtb*dxx(lx)
		     end do
		  end do
	       end if
c
	       dvl(j)=dvl2(j)+dvl3(j)/akt
               dvel(j)=dvel2(j)+dvel3(j)/akt
               dvbl(j)=dvbl2(j)+dvbl3(j)/akt
730	    continue

790	    continue
c
c
c
c Adding to calculate the integral over k.
c Vector case
	    do in=1,nn
	       do j=1,l0
           call powertflat(akt,in,apowert)
           cvkj=apowert*dvl(j)*dvl(j)*dak1(k)
           cvl(j,in)=cvl(j,in)+cvkj

           cvekj=apowert*dvel(j)*dvel(j)*dak1(k)
           cvel(j,in)=cvel(j,in)+cvekj

           cvckj=apowert*dvl(j)*dvel(j)*dak1(k)
           cvcl(j,in)=cvcl(j,in)+cvckj

           cvbkj=apowert*dvbl(j)*dvbl(j)*dak1(k)
           cvbl(j,in)=cvbl(j,in)+cvbkj

           cvbtkj=apowert*dvbl(j)*dvl(j)*dak1(k)
           cvbtl(j,in)=cvbtl(j,in)+cvbtkj

           cvebkj=apowert*dvel(j)*dvbl(j)*dak1(k)
           cvebl(j,in)=cvebl(j,in)+cvebkj

	       end do
	    end do
c
	 end if
c
c
200	continue
c
c
c
c Final calculations for CMB output.
c
	do j=1,l0
	   xl(j)=dble(l(j))
	end do

c Scalar case
	if (itflag.ne.2.and.itflag.ne.3) then
	   do 600 in=1,nn
	      do 620 j=1,l0
                 ctnorm=dble(l(j)*(l(j)-1.0d0)*(l(j)+1)*(l(j)+2))

ccccc    cl(j,in)=8.0d0*cl(j,in)*dble(l(j)*(l(j)+1))/fourpi

                 cl(j,in)=fourpi*fourpi*cl(j,in)*dble(l(j)*(l(j)+1))

		 cpl(j,in)=fourpi*fourpi*ctnorm*cpl(j,in)
     &                                         *dble(l(j)*(l(j)+1))
		 ccl(j,in)=fourpi*fourpi*sqrt(ctnorm)
     &                               *ccl(j,in)*dble(l(j)*(l(j)+1))
 620	      continue
 600	   continue
c Making the interpolation tables to get other l-values.
           llo=1
           cllo=1.0d30
           clhi=1.0d30
           do 630 in=1,nn
              call spline(xl,cl(1,in),l0,cllo,clhi,clpr(1,in))
              call spline(xl,cpl(1,in),l0,cllo,clhi,cplpr(1,in))
              call spline(xl,ccl(1,in),l0,cllo,clhi,cclpr(1,in))
 630       continue
        end if
c
c Tensor Case
	if (itflag.ne.0.and.itflag.ne.3) then
c Normalization
	   do 610 j=1,l0
	      ctnorm=dble(l(j)*(l(j)-1.0d0)*(l(j)+1)*(l(j)+2))
            do 640 in=1,nn
	      ctl(j,in)=fourpi*fourpi*ctnorm*ctl(j,in)
     &                                  *dble(l(j)*(l(j)+1))

	      ctel(j,in)=fourpi*fourpi*ctel(j,in)*dble(l(j)*(l(j)+1))

	      ctbl(j,in)=fourpi*fourpi*ctbl(j,in)*dble(l(j)*(l(j)+1))

	      ctcl(j,in)=fourpi*fourpi*sqrt(ctnorm)*ctcl(j,in)
     &                                           *dble(l(j)*(l(j)+1))

              ctebl(j,in)=fourpi*fourpi*ctebl(j,in)*dble(l(j)*(l(j)+1))

	      ctbtl(j,in)=fourpi*fourpi*sqrt(ctnorm)*ctbtl(j,in)
     &                                           *dble(l(j)*(l(j)+1))

 640	   continue
 610	   continue
c Making the interpolation tables to get other l-values.
	   cllo=1.0d30
	   clhi=1.0d30
           do 650 in=1,nn
	    call spline(xl,ctl(1,in),l0,cllo,clhi,ctlpr(1,in))
	    call spline(xl,ctel(1,in),l0,cllo,clhi,ctelpr(1,in))
	    call spline(xl,ctbl(1,in),l0,cllo,clhi,ctblpr(1,in))
	    call spline(xl,ctcl(1,in),l0,cllo,clhi,ctclpr(1,in))
            call spline(xl,ctbtl(1,in),l0,cllo,clhi,ctbtlpr(1,in))
            call spline(xl,ctebl(1,in),l0,cllo,clhi,cteblpr(1,in))
650        continue
	end if
c
c Vector Case
	if (itflag.ne.0.and.itflag.ne.2) then
c Normalization
	   do 810 j=1,l0
            do 840 in=1,nn
            cvnorm=dble((l(j)-1)*(l(j)+2))
       cvl(j,in)=fourpi*fourpi*cvl(j,in)*dble(l(j)*(l(j)+1))**2

       cvel(j,in)=
     &  fourpi*fourpi*cvnorm*cvel(j,in)*dble(l(j)*(l(j)+1))

       cvbl(j,in)=
     &  fourpi*fourpi*cvnorm*cvbl(j,in)*dble(l(j)*(l(j)+1))

       cvcl(j,in)=
     &  fourpi*fourpi*sqrt(cvnorm)*cvcl(j,in)
     &  *dble(l(j)*(l(j)+1))**(3.d0/2.d0)

       cvbtl(j,in)=
     &  fourpi*fourpi*cvnorm*cvbtl(j,in)*dble(l(j)*(l(j)+1))

       cvebl(j,in)=
     &  fourpi*fourpi*cvnorm*cvebl(j,in)*dble(l(j)*(l(j)+1))

c If computing cross-correlators, remember to include a factor of sqrt(2)!!

840	   continue
810	   continue
c Making the interpolation tables to get other l-values.
	   cllo=1.0d30
	   clhi=1.0d30
           do 850 in=1,nn
	    call spline(xl,cvl(1,in),l0,cllo,clhi,cvlpr(1,in))
            call spline(xl,cvel(1,in),l0,cllo,clhi,cvelpr(1,in))
            call spline(xl,cvbl(1,in),l0,cllo,clhi,cvblpr(1,in))
            call spline(xl,cvcl(1,in),l0,cllo,clhi,cvclpr(1,in))
            call spline(xl,cvbtl(1,in),l0,cllo,clhi,cvbtlpr(1,in))
            call spline(xl,cvebl(1,in),l0,cllo,clhi,cveblpr(1,in))
850        continue
	end if
c
c Calculating Cls for every l.

        do 700 in=1,nn
	   llo=1
	   do 710 il=2,l(l0)
	      xi=il
	      if ((xi.gt.xl(llo+1)).and.(llo.lt.l0)) then
		 llo=llo+1
	      end if
	      lhi=llo+1
	      ho=xl(lhi)-xl(llo)
	      a0=(xl(lhi)-xi)/ho
	      b0=(xi-xl(llo))/ho
	      if (itflag.ne.2.and.itflag.ne.3) then
		 clint=a0*cl(llo,in)+b0*cl(lhi,in)+((a0**3-a0)*
     &              clpr(llo,in)+(b0**3-b0)*clpr(lhi,in))*ho*ho /6.d0
		 cplint=a0*cpl(llo,in)+b0*cpl(lhi,in)+((a0**3-a0)
     &            *cplpr(llo,in)+(b0**3-b0)*cplpr(lhi,in))*ho*ho /6.d0
		 cclint=a0*ccl(llo,in)+b0*ccl(lhi,in)+((a0**3-a0)
     &            *cclpr(llo,in)+(b0**3-b0)*cclpr(lhi,in))*ho*ho /6.d0
	      end if
              if (itflag.ne.0.and.itflag.ne.3) then
		ctlint=a0*ctl(llo,in)+b0*ctl(lhi,in)+((a0**3-a0)
     &            *ctlpr(llo,in)+(b0**3-b0)*ctlpr(lhi,in))*ho*ho/6.d0
		ctelint=a0*ctel(llo,in)+b0*ctel(lhi,in)+((a0**3-a0)
     &            *ctelpr(llo,in)+(b0**3-b0)*ctelpr(lhi,in))*ho*ho/6.d0
		ctblint=a0*ctbl(llo,in)+b0*ctbl(lhi,in)+((a0**3-a0)
     &            *ctblpr(llo,in)+(b0**3-b0)*ctblpr(lhi,in))*ho*ho/6.d0
		ctclint=a0*ctcl(llo,in)+b0*ctcl(lhi,in)+((a0**3-a0)
     &            *ctclpr(llo,in)+(b0**3-b0)*ctclpr(lhi,in))*ho*ho/6.d0
               ctbtlint=a0*ctbtl(llo,in)+b0*ctbtl(lhi,in)+((a0**3-a0)
     &          *ctbtlpr(llo,in)+(b0**3-b0)*ctbtlpr(lhi,in))*ho*ho/6.d0
               cteblint=a0*ctebl(llo,in)+b0*ctebl(lhi,in)+((a0**3-a0)
     &          *cteblpr(llo,in)+(b0**3-b0)*cteblpr(lhi,in))*ho*ho/6.d0
	     end if
              if (itflag.ne.0.and.itflag.ne.2) then
		cvlint=a0*cvl(llo,in)+b0*cvl(lhi,in)+((a0**3-a0)
     &            *cvlpr(llo,in)+(b0**3-b0)*cvlpr(lhi,in))*ho*ho/6.d0
       cvelint=a0*cvel(llo,in)+b0*cvel(lhi,in)+((a0**3-a0)
     &          *cvelpr(llo,in)+(b0**3-b0)*cvelpr(lhi,in))*ho*ho/6.d0
       cvblint=a0*cvbl(llo,in)+b0*cvbl(lhi,in)+((a0**3-a0)
     &          *cvblpr(llo,in)+(b0**3-b0)*cvblpr(lhi,in))*ho*ho/6.d0
       cvbtlint=a0*cvbtl(llo,in)+b0*cvbtl(lhi,in)+((a0**3-a0)
     &         *cvbtlpr(llo,in)+(b0**3-b0)*cvbtlpr(lhi,in))*ho*ho/6.d0
       cveblint=a0*cvebl(llo,in)+b0*cvebl(lhi,in)+((a0**3-a0)
     &         *cveblpr(llo,in)+(b0**3-b0)*cveblpr(lhi,in))*ho*ho/6.d0
       cvclint=a0*cvcl(llo,in)+b0*cvcl(lhi,in)+((a0**3-a0)
     &          *cvclpr(llo,in)+(b0**3-b0)*cvclpr(lhi,in))*ho*ho/6.d0
	     end if
	     if (itflag.ne.2.and.itflag.ne.3) then
                clts(il,in)=clint
                cles(il,in)=cplint
                clcs(il,in)=cclint
	     end if
	     if (itflag.ne.0.and.itflag.ne.3) then
                cltt(il,in)=ctlint
                clet(il,in)=ctelint
                clbt(il,in)=ctblint
                clct(il,in)=ctclint
                clebt(il,in)=cteblint
                clbtt(il,in)=ctbtlint
	     end if
	     if (itflag.ne.0.and.itflag.ne.2) then
                cltv(il,in)=cvlint
                clev(il,in)=cvelint
                clbv(il,in)=cvblint
                clcv(il,in)=cvclint
                clebv(il,in)=cveblint
                clbtv(il,in)=cvbtlint
	     end if
 710      continue
          if (itflag.ne.2.and.itflag.ne.3) then
             clts(1,in)=cl(1,in)
          end if
          if (itflag.ne.0.and.itflag.ne.3) then
             cltt(1,in)=ctl(1,in)
          end if
          if (itflag.ne.0.and.itflag.ne.2) then
             cltv(1,in)=cvl(1,in)
          end if
 700      continue
c
       end if

       return
       end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine fderivs(n,x,y,yprime)
c  Evaluate the time derivatives of the perturbations.
c
	implicit double precision (a-h,o-z)
        integer initfl
	dimension y(n),yprime(n)

c      common/volcom/tkmax
      common/picom/pi,twopi3
c      common/maxvcom/maxvar   
  
	common /lingerinc/ omegab,omegac,omegav,omegan,omegak,h0
     &                     ,tcmb,yhe,nnur,nnunr
	common /cosmoparm/ ak,ak2,amnu,lmax,lmaxnr,lmaxnu,
     &                     nqmax,iq0,iq1,iq2
	common /genparm/ grhom,grhog,grhor,adotrad,taurst,dtaurec
	common /initcase/ initfl
	common /out1/ adotoa,hdot,dgshear,rhonu,shearnu
c
	parameter (lmaxnu0=25,nqmax0=15)
	parameter (lmx0=30)
	parameter (ep0=1.0d-2)
c  Internal variables.
	dimension denl(lmx0),dlfdlq(nqmax0),akv(nqmax0)
	common /store/ denl,dlfdlq
c
	common /epsilon/ epsw
c
c ep is used to stop the tight coupling approximation.
	if (ak.gt.epsw) then
	   ep=ep0
	else
	   ep=0.5d0*ep0
	end if
	tau=x
	a=y(1)
c	ahdot=y(2)
	eta=y(3)

        if (initfl.eq.4) then
	 call scalar_source(tau,emt00,emt00dot,emtS,emtSdot,emtD,emtP)
	else
	 emt00=0.0d0
	 emt00dot=0.0d0
	 emtS=0.0d0
	 emtSdot=0.0d0
	 emtD=0.0d0
	 emtP=0.0d0
        endif
c  CDM.
	deltac=y(4)
	thetac=y(5)
c  Baryons.
	deltab=y(6)
	thetab=y(7)
c  Photons.
	deltag=y(8)
	thetag=y(9)
	shearg=y(10)/2.0d0
c  Polarization term.
	polter=y(10)+y(9+lmax)+y(11+lmax)
c  Massless neutrinos.
	deltar=y(10+2*lmax)
	thetar=y(11+2*lmax)
	shearr=y(12+2*lmax)/2.0d0

	a2=a*a
	call thermo(tau,xcs2,xopac)

c Tight Coupling parameters
	tcp=0.0d0
	tcp1=ak/xopac
	tcp2=1.0d0/(xopac*tau)
	if ((tcp1.gt.ep).or.(tcp2.gt.ep)) then
	   tcp=1.0d0
	end if   
c  Photon mass density over baryon mass density.
	photbar=grhog/(grhom*omegab*a)
	pb43=4.0d0/3.0d0*photbar
c  Compute expansion rate.
	if (amnu.eq.0.0d0) then
	   rhonu=1.0d0
	   pnu=1.0d0/3.0d0
	   drhonu=0.0d0
	   fnu=0.0d0
	   dpnu=0.0d0
	   shearnu=0.0d0
	else
	   call nu1(a,rhonu,pnu)
	   call nu2(a,drhonu,fnu,dpnu,shearnu,y(iq0),y(iq1),y(iq2))
	end if
c  8*pi*G*rho*a**2 and 8*pi*G*P*a**2.
	grho=grhom*(omegac+omegab)/a+(grhog+grhor*(nnur+nnunr*rhonu))/a2
     2      +grhom*omegav*a2
	adotoa=sqrt(grho/3.0d0)
	yprime(1)=adotoa*a
	gpres=((grhog+grhor*nnur)/3.0d0+grhor*nnunr*pnu)/a2
     2      -grhom*omegav*a2
c  Evaluate metric and massive neutrino perturbations.
c	deltan=drhonu/rhonu
c	thetan=ak*fnu/(rhonu+pnu)
c
cLP  used to CALCULATE THE REST OF STRING COMPONENTS    
c      emtP=(emtD-emt00dot)/adotoa-emt00
c  String conservation equation
c      if(ak*tau.gt.tkmax) then
c        yprime(maxvar)=0.0d0
c       else
c      emtDdot=-2.0d0*adotoa*emtD-(ak**2/3.0d0)*(emtP+2.0d0*emtS)
c       yprime(maxvar)=emtDdot
c      endif
c
c
c  8*pi*G*delta_rho*a**2 and 8*pi*G*delta_P*a**2.
	dgrho=grhom*(omegac*deltac+omegab*deltab)/a
     2    +(grhog*deltag+grhor*(nnur*deltar+nnunr*drhonu))/a2
     &    +8.0d0*pi*emt00
	dgpres=(grhog*deltag+grhor*nnur*deltar)/a2/3.0d0
     2    +grhor*nnunr*dpnu/a2
     &    +8.0d0*pi*emtP/3.0d0

	dahdotdtau=-(dgrho+3.0d0*dgpres)*a
	yprime(2)=dahdotdtau
c  Force energy conservation.
	hdot=(2.0d0*ak2*eta+dgrho)/adotoa
c  8*pi*G*(rho+P)*theta*a**2.
	dgtheta=grhom*(omegac*thetac+omegab*thetab)/a
     2    +4.0d0/3.0d0*(grhog*thetag+nnur*grhor*thetar)/a2
     3    +nnunr*grhor*ak*fnu/a2
     &    -8.0d0*pi*emtD
	etadot=0.5d0*dgtheta/ak2
	yprime(3)=etadot
c  8*pi*G*(rho+P)*sigma*a**2.
	dgshear=4.0d0/3.0d0*(grhog*shearg+nnur*grhor*shearr)/a2
     2    +nnunr*grhor*shearnu/a2
     &    -16.0d0*pi*emtS/3.0d0
c  CDM equations of motion.
	deltacdot=-thetac-0.5d0*hdot
	yprime(4)=deltacdot
	thetacdot=-adotoa*thetac
	yprime(5)=thetacdot
c  Baryon equations of motion.
	deltabdot=-thetab-0.5d0*hdot
	yprime(6)=deltabdot
c  Need photon perturbation for first-order correction to tightly-coupled
c  baryon-photon approximation.
	deltagdot=4.0d0/3.0d0*(-thetag-0.5d0*hdot)
	drag=xopac*(thetag-thetab)
	if (tcp.eq.1) then
c  Treat baryons and photons as uncoupled.
	  thetabdot=-adotoa*thetab+ak2*xcs2*deltab+pb43*drag
	else
c  Treat baryons and photons as tightly coupled.
c  Zeroth-order approximation to baryon velocity.
	  thetabdot=(-adotoa*thetab+ak2*xcs2*deltab
     2  +ak2*pb43*(0.25d0*deltag-shearg))/(1.0d0+pb43)
c  (\ddot a)/a.
	  adotdota=0.5d0*(adotoa*adotoa-gpres)
c  First-order approximation to baryon-photon slip, thetabdot-thetagdot.
	  slip=2.0d0*pb43/(1.0d0+pb43)*adotoa*(thetab-thetag)
     2     +1.0d0/xopac*(-adotdota*thetab-adotoa*ak2*0.5d0*deltag
     3     +ak2*(xcs2*deltabdot-0.25d0*deltagdot))/(1.0d0+pb43)
c  First-order approximation to baryon velocity.
	  thetabdot=thetabdot+pb43/(1.0d0+pb43)*slip
	end if
	yprime(7)=thetabdot
c  Photon total intensity and polarization equations of motion.
	yprime(8)=deltagdot
	thetagdot=(-thetabdot-adotoa*thetab+ak2*xcs2*deltab)/pb43
     2   +ak2*(0.25d0*deltag-shearg)
	yprime(9)=thetagdot
	if (tcp.eq.1) then
c  Treat baryons and photons as uncoupled.
	  yprime(10)=8.0d0/15.0d0*thetag-0.6d0*ak*y(11)-xopac*y(10)
     2        +4.0d0/15.0d0*hdot+8.0d0/5.0d0*etadot+0.1d0*xopac*polter
c  Polarization equations for l = 0, 1, 2.
	  yprime(9+lmax)=-ak*y(10+lmax)-xopac*y(9+lmax)+0.5d0*xopac*polter
	  yprime(10+lmax)=ak/3.0d0*(y(9+lmax)-2.0d0*y(11+lmax))
     2             -xopac*y(10+lmax)
	  yprime(11+lmax)=ak*(0.4d0*y(10+lmax)-0.6d0*y(12+lmax))
     2             -xopac*y(11+lmax)+0.1d0*xopac*polter
	    do 10 l=3,lmax-1
	    yprime(8+l)=ak*denl(l)*(l*y(7+l)-(l+1)*y(9+l))-xopac*y(8+l)
	    yprime(9+lmax+l)=ak*denl(l)*(l*y(8+lmax+l)-(l+1)*
     2                 y(10+lmax+l))-xopac*y(9+lmax+l)
10	  continue
	else
c  Treat baryons and photons as tightly coupled (with no polarization).
	  yprime(10)=0.0d0
	  yprime(9+lmax)=0.0d0
	  yprime(10+lmax)=0.0d0
	  yprime(11+lmax)=0.0d0
	    do 15 l=3,lmax-1
	    yprime(8+l)=0.0d0
	    yprime(9+lmax+l)=0.0d0
15	  continue
	end if
c  Truncate moment expansion.
c	yprime(8+lmax)=ak*lmax*y(7+lmax)/(2*lmax+1)-xopac*y(8+lmax)
c	yprime(9+2*lmax)=ak*lmax*y(8+2*lmax)/(2*lmax+1)-xopac*y(9+2*lmax)
	yprime(8+lmax)=ak*y(7+lmax)-(lmax+1)/tau*y(8+lmax)
     2                -xopac*y(8+lmax)
	yprime(9+2*lmax)=ak*y(8+2*lmax)-(lmax+1)/tau*y(9+2*lmax)
     2                -xopac*y(9+2*lmax)
c  Massless neutrino equations of motion.
	deltardot=4.0d0/3.0d0*(-thetar-0.5d0*hdot)
	yprime(10+2*lmax)=deltardot
	thetardot=ak2*(0.25d0*deltar-shearr)
	yprime(11+2*lmax)=thetardot
	yprime(12+2*lmax)=8.0d0/15.0d0*thetar-0.6d0*ak*y(13+2*lmax)
     2             +4.0d0/15.0d0*hdot+8.0d0/5.0d0*etadot
	  do 20 l=3,lmaxnr-1
	  yprime(10+2*lmax+l)=ak*denl(l)*(l*y(9+2*lmax+l)
     2                    -(l+1)*y(11+2*lmax+l))
20	continue
c  Truncate moment expansion.
c	yprime(10+2*lmax+lmaxnr)=ak*lmax*y(9+2*lmax+lmaxnr)/(2*lmaxnr+1)
	yprime(10+2*lmax+lmaxnr)=ak*y(9+2*lmax+lmaxnr)-(lmaxnr+1)/tau
     7  *y(10+2*lmax+lmaxnr)
c  Massive neutrino equations of motion.
	if (nqmax.eq.0) return
c	dq=qmax/nqmax
	dq=1.0d0
	  do i=1,nqmax
	  q=i*dq-0.5d0
	  aq=a*amnu/q
	  v=1.0d0/sqrt(1.0d0+aq*aq)
	  akv(i)=ak*v
	end do
c  l = 0, 1, 2,lmaxnu.
	  do 30 i=1,nqmax
	  ind=iq0+i-1
	  yprime(ind)=-akv(i)*y(ind+nqmax)+hdot*dlfdlq(i)/6.0d0
	  ind=iq1+i-1
	  yprime(ind)=akv(i)*(y(ind-nqmax)-2*y(ind+nqmax))/3
	  ind=iq2+i-1
	  yprime(ind)=akv(i)*(2*y(ind-nqmax)-3*y(ind+nqmax))/5
     2           -(hdot/15.0d0+2.0d0/5.0d0*etadot)*dlfdlq(i)
	  ind=10+2*lmax+lmaxnr+i+lmaxnu*nqmax
c  Truncate moment expansion.
c	  yprime(ind)=akv*lmaxnu*y(ind-nqmax)/(2*lmaxnu+1)
	  yprime(ind)=akv(i)*y(ind-nqmax)-(lmaxnu+1)/tau*y(ind)
30	continue
	  do 50 l=3,lmaxnu-1
	    do 40 i=1,nqmax
	    ind=10+2*lmax+lmaxnr+i+l*nqmax
          yprime(ind)=akv(i)*denl(l)*(l*y(ind-nqmax)-(l+1)*y(ind+nqmax))
40	  continue
50	continue
c
	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine foutput(n,y,yprime,j,tau0,tau,d,dp)
        implicit double precision (a-h,o-z)
        integer n,j
c COMMON picom,lingerinc,genparm
	parameter (nstep0=6000)
        common/picom/pi,twopi3
	
	integer initfl
	common /initcase/ initfl

	common /lingerinc/ omegab,omegac,omegav,omegan,omegak,h0
     &                     ,tcmb,yhe,nnur,nnunr
	common /cosmoparm/ ak,ak2,amnu,lmaxg,lmaxnr,lmaxnu,
     &                     nqmax,iq0,iq1,iq2
	common /genparm/ grhom,grhog,grhor,adotrad,taurst,dtaurec
	common /out1/ adotoa,hdot,dgshear,rhonu,shearnu
c
        dimension y(n),yprime(n)
        dimension vis(nstep0),dvis(nstep0),ddvis(nstep0)
        dimension opac(nstep0),dopac(nstep0),expmmu(nstep0)
        common /visib/ vis,dvis,ddvis,expmmu,opac,dopac        
        save /visib/


        if (initfl.eq.4) then
	 call scalar_source(tau,emt00,emt00dot,emtS,emtSdot,emtD,emtP)
	else
	 emt00=0.0d0
	 emt00dot=0.0d0
	 emtS=0.0d0
	 emtSdot=0.0d0
	 emtD=0.0d0
	 emtP=0.0d0
        endif 

c       call splint1(ts_f,emS,emSpr,nt_f,tau,emtS)
c       call splint2(ts_f,emS,emSpr,nt_f,tau,emtSdot)
c
        x=ak*(tau0-tau)
	a=y(1)
	a2=a*a
c	ahdot=y(2)
	eta=y(3)
	etadot=yprime(3)
        alpha=(hdot+6*etadot)/(2.0d0*ak2)
        alphadot=-3*dgshear/(2.0d0*ak2)+eta-2.0d0*adotoa*alpha
c  Baryons.
c	deltab=y(6)
	thetab=y(7)
        thetabdot=yprime(7)
c  Photons.
	deltag=y(8)
c	thetag=y(9)
c	shearg=y(10)/2.0d0
	thetagdot=yprime(9)
        sheargdot=yprime(10)/2.0d0
c  Polarization term.
	polter=y(10)+y(9+lmaxg)+y(11+lmaxg)
c        coupl=8.0d0*(thetag+ak2*alpha)/15.0d0-ak*0.6d0
c     2      *(y(11)+y(10+lmaxg)+y(12+lmaxg))
        coupldot=8.0d0*(thetagdot+ak2*alphadot)/15.0d0
     2      -ak*0.6d0*(yprime(11)+yprime(10+lmaxg)+yprime(12+lmaxg))
        polterdot=yprime(10)+yprime(9+lmaxg)+yprime(11+lmaxg)
        polterddot=coupldot-0.3d0*(dopac(j)*polter+opac(j)*polterdot)
c  Massless neutrinos.
c	deltar=y(10+2*lmaxg)
c	thetar=y(11+2*lmaxg)
c	shearr=y(12+2*lmaxg)/2.0d0
	shearrdot=yprime(12+2*lmaxg)/2.0d0
c  Second derivative of expansion rate
	if (amnu.eq.0.0) then
	   rhonudot=0.0d0
	   shearnudot=0.0d0
	else
	   call nuder(a,adotoa,rhonu,rhonudot,shearnudot,
     2     		y(iq2),yprime(iq2))
	end if  
	grhodot=(-grhom*(omegac+omegab)/a-2.0d0
     2      *(grhog+grhor*(nnur+nnunr*rhonu))/a2
     3      +2*grhom*omegav*a2)*adotoa
     4      +grhor*nnunr*rhonudot/a2
        adotoadot=grhodot/(6*adotoa)
c  Derivative of the shear
	dgsheardot=4.0d0/3.0d0*(grhog*sheargdot+nnur*grhor
     2  *shearrdot)/a2-2.0d0*adotoa*(dgshear+16.0d0*pi*emtS/3.0d0)
     3  +nnunr*grhor*shearnudot/a2
     &  -16.0d0*pi*emtSdot/3.0d0        
c  Calculation of the sources
        alphaddot=-3*dgsheardot/(2.0d0*ak2)+etadot
     2      -2.0d0*adotoadot*alpha-2.0d0*adotoa*alphadot
c
         s1=etadot+alphaddot
         s2=2*alphadot
c
        d=expmmu(j)*s1+vis(j)*(0.25d0*deltag
     2     +s2+polter/16+thetabdot/ak2+3.0d0/16.0d0
     3     *polterddot/ak2)+dvis(j)*(alpha+thetab/ak2+3.0d0/8.0d0
     4     *polterdot/ak2)+ddvis(j)*3.0d0/16.0d0*polter/ak2

        if (x.gt.0.0d0) then
           dp=vis(j)*3.0d0/16.0d0*polter/x**2
        else
           dp=0.0d0
        end if
c       
	return
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine fderivst(n,x,y,yprime)
c  Evaluate the time derivatives of the perturbations.
c COMMON volcom,picom,genparm,lingerinc,tensor,store
	implicit double precision (a-h,o-z)
	dimension y(n),yprime(n)
c
	integer initfl
	common /initcase/ initfl
	
        parameter (nnmax=1)
	common /lingerinc/ omegab,omegac,omegav,omegan,omegak,h0
     &                     ,tcmb,yhe,nnur,nnunr
	common /cosmoparm/ ak,ak2,amnu,lmax,lmaxnr,lmaxnu,
     &                     nqmax,iq0,iq1,iq2
c>>>>>>>>>>>>>>>>
      common/volcom/tkmax
      common/picom/pi,twopi3 
      
      parameter (nt_fine=5000)
      common/splined_strings/ts_f(nt_fine)
     &                      ,em00(nt_fine),em00pr(nt_fine)
     &                      ,emS(nt_fine),emSpr(nt_fine)
     &                      ,emD(nt_fine),emDpr(nt_fine)
     &                      ,emP(nt_fine),emPpr(nt_fine)          
     &                      ,emV(nt_fine),emVpr(nt_fine)
     &                      ,emT(nt_fine),emTpr(nt_fine)	    
     &                      ,nt_f    
      
c<<<<<<<<<<<<<<<<                 
	common /genparm/ grhom,grhog,grhor,adotrad,taurst,dtaurec
        common /tensor/ ant(nnmax),rat(nnmax),itflag,lmaxt
c        common /reionization/zri,taurist,zristp,tauristp,rif
c        common /reionization2/ j2ri1,nri,nri0
c
	parameter (nqmax0=15)
	parameter(lmx0=30)
	parameter (ep0=1.0d-2)
c
	common /epsilon/ epsw
c
c  Internal variables.
	dimension denl(lmx0),dlfdlq(nqmax0)
	common /store/ denl,dlfdlq
c

c ep is used to stop the tight coupling approximation.
	if (ak.gt.0.06d0*epsw) then
	   ep=ep0
	else
           ep=1.17d0*ep0
	end if

	tau=x
	a=y(1)
c
	a2=a*a
	call thermo(tau,xcs2,xopac)
c
c Tight Coupling parameters
	tcp=0.0d0
	tcp1=ak/xopac
	tcp2=1.0d0/(xopac*tau)
	if ((tcp1.gt.ep).or.(tcp2.gt.ep)) then
	   tcp=1.0d0
	end if   
        
c  Compute expansion rate.
	if (amnu.eq.0) then
	   rhonu=1.0d0
	   pnu=1.0d0/3.0d0
	else
	   call nu1(a,rhonu,pnu)
	end if
c  8*pi*G*rho*a**2 and 8*pi*G*P*a**2.
	grho=grhom*(omegac+omegab)/a+(grhog+grhor*(nnur+nnunr*rhonu))/a2
     2      +grhom*omegav*a2
	adotoa=sqrt(grho/3.0d0)
	yprime(1)=adotoa*a

c Tensors
	ht=y(2)
	htpr=y(3)
	yprime(2)=htpr
c>>>>>>>>>>>>>>>>>>      
c        if(ak*tau.gt.tkmax) then
c        emtT=0.0d0
c       else

        if (initfl.eq.4) then
       call splint1(ts_f,emT,emTpr,nt_f,tau,emtT)
	else
	 emtT=0.0d0
        endif 

       call splint1(ts_f,emT,emTpr,nt_f,tau,emtT)

c      endif
c      	
	htdpr=-2*adotoa*htpr-ak2*ht+8.0d0*pi*emtT
c<<<<<<<<<<<<<<<<<<	
	yprime(3)=htdpr
c
c Photon perturbations
	ind1=4
	ind2=ind1+lmaxt+1
	psie=y(ind1)/10.0d0+y(ind1+2)/7.0d0+3d0*y(ind1+4)/70.0d0
     &       -3.0d0*y(ind2)/5.0d0+
     &      6.0d0*y(ind2+2)/7.0d0-3.0d0*y(ind2+4)/70.0d0
	if (tcp.eq.1) then
c no tight coupling approx
	   yprime(ind1)=-ak*y(ind1+1)-xopac*y(ind1)+xopac*psie-htpr
	   yprime(ind2)=-ak*y(ind2+1)-xopac*y(ind2)-xopac*psie
c l=1...lmaxt
	   do 70 l=1,lmaxt-1
	      yprime(ind1+l)=ak*denl(l)*(l*y(ind1-1+l)-(l+1)
     &   *y(ind1+1+l))-xopac*y(ind1+l)
	      yprime(ind2+l)=ak*denl(l)*(l*y(ind2-1+l)-(l+1)
     &   *y(ind2+1+l))-xopac*y(ind2+l)
 70	   continue
c
c Truncate moment expansion
	   yprime(ind1+lmaxt)=ak*y(ind1-1+lmaxt)-(lmaxt+1)/tau
     2	   *y(ind1+lmaxt)-xopac*y(ind1+lmaxt)
	   yprime(ind2+lmaxt)=ak*y(ind2-1+lmaxt)-(lmaxt+1)/tau
     2	   *y(ind2+lmaxt)-xopac*y(ind2+lmaxt)
	else
           deltat0=-4.0d0*htpr/xopac/3.0d0
           deltap0=-deltat0/4.0d0
           y(ind1)=deltat0
           y(ind2)=deltap0
	   do 80 l=0,lmaxt
	      yprime(ind1+l)=0.0d0
	      yprime(ind2+l)=0.0d0
 80	   continue
	end if
	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine foutputt(n,y,ypr,j,tau0,tau,dt,dte,dtb)

        implicit double precision (a-h,o-z)
        integer j,n
        parameter (nnmax=1)
c COMMON tensor
	parameter (nstep0=6000)
	common /cosmoparm/ ak,ak2,amnu,lmaxg,lmaxnr,lmaxnu,
     &                     nqmax,iq0,iq1,iq2
        common /tensor/ ant(nnmax),rat(nnmax),itflag,lmaxt
c
        dimension y(n),ypr(n)
        dimension vis(nstep0),dvis(nstep0),ddvis(nstep0)
        dimension opac(nstep0),dopac(nstep0),expmmu(nstep0)
        common /visib/ vis,dvis,ddvis,expmmu,opac,dopac 
c                       
        save /visib/
c     
c Tensors
c
        x=ak*(tau0-tau)
        x2=x*x
        ind1=4
        ind2=ind1+lmaxt+1
        htpr=y(3)
        htdpr=ypr(3)
        psie=y(ind1)/10.0d0+y(ind1+2)/7.0d0+3.0d0
     & *y(ind1+4)/70.0d0-3.d0*y(ind2)/5.0d0+6.0d0*y(ind2+2)
     & /7.0d0-3.0d0*y(ind2+4)/70.0d0
        psiedot=ypr(ind1)/10.0d0+ypr(ind1+2)/7.0d0+3.0d0
     & *ypr(ind1+4)/70.0d0-3.d0*ypr(ind2)/5.0d0+6.0d0*ypr(ind2+2)
     & /7.0d0-3.0d0*ypr(ind2+4)/70.0d0
        psieddot=-0.3d0*(opac(j)*psiedot+dopac(j)*psie)
     & -0.1d0*htdpr-ak*(3.0d0*ypr(ind1+1)/70.0d0+ypr(ind1+3)/15.0d0
     & +ypr(ind1+5)/42.0d0-33.0d0*ypr(ind2+1)/35.0d0
     & +8.0d0*ypr(ind2+3)/15.0d0-ypr(ind2+5)/42.0d0)
c

        if (x.gt.0.0d0) then
           dt=(-expmmu(j)*htpr+vis(j)*psie)/x2

           dte=vis(j)*(psie-psieddot/ak2-6.0d0*psie/x2
     &     -4.0d0*psiedot/ak/x)-dvis(j)*(4.0d0*psie/x/ak
     &     +2.0d0*psiedot/ak2)-ddvis(j)*psie/ak2

           dtb=2.0d0*(vis(j)*(2.0d0*psie/x+psiedot/ak)
     &     +dvis(j)*psie/ak)
        else
           dt=0.0d0
           dte=0.0d0
           dtb=0.0d0
        end if     
c        
        dte=-dte
        dtb=-dtb
	return
        end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine fderivsv(n,x,y,yprime)
c  Evaluate the time derivatives of the vector perturbations.
c
	implicit double precision (a-h,o-z)
	dimension y(n),yprime(n)
c COMMON volcom,picom,genparm,lingerinc,tensor,store,vector
c COMMON vectime
        parameter (nnmax=1)
	common /lingerinc/ omegab,omegac,omegav,omegan,omegak,h0
     &                     ,tcmb,yhe,nnur,nnunr
	common /cosmoparm/ ak,ak2,amnu,lmax,lmaxnr,lmaxnu,
     &                     nqmax,iq0,iq1,iq2
c>>>>>>>>>>>>>>>>
      common/volcom/tkmax
      common/picom/pi,twopi3

	integer initfl
	common /initcase/ initfl
      
      parameter (nt_fine=5000)
      common/splined_strings/ts_f(nt_fine)
     &                      ,em00(nt_fine),em00pr(nt_fine)
     &                      ,emS(nt_fine),emSpr(nt_fine)
     &                      ,emD(nt_fine),emDpr(nt_fine)
     &                      ,emP(nt_fine),emPpr(nt_fine)          
     &                      ,emV(nt_fine),emVpr(nt_fine)
     &                      ,emT(nt_fine),emTpr(nt_fine)	    
     &                      ,nt_f    
       
c<<<<<<<<<<<<<<<<                 
	common /genparm/ grhom,grhog,grhor,adotrad,taurst,dtaurec
        common /tensor/ ant(nnmax),rat(nnmax),itflag,lmaxt
        common /vector/ lmaxv
c        common /reionization/zri,taurist,zristp,tauristp,rif
c        common /reionization2/ j2ri1,nri,nri0
      parameter (lmaxv0=8)
      common/vectime/dl2m1(lmaxv0),dl2p2l(lmaxv0),dkap(lmaxv0),
     &               dkap1(lmaxv0),dlp1(lmaxv0) 
c
	parameter (nqmax0=15)
	parameter(lmx0=30)
	parameter (ep0=1.0d-2)
c
	common /epsilon/ epsw
c
c  Internal variables.
	dimension denl(lmx0),dlfdlq(nqmax0)
	common /store/ denl,dlfdlq
c

c ep is used to stop the tight coupling approximation.
	if (ak.gt.0.01d0*epsw) then
	   ep=ep0
	else
c           ep=1.17d0*ep0
        ep=0.5d0*ep0
	end if

	tau=x
	a=y(1)
c
	a2=a*a
	call thermo(tau,xcs2,xopac)
c
c Tight Coupling parameters
	tcp=0.0d0
	tcp1=ak/xopac
	tcp2=1.0d0/(xopac*tau)
	if ((tcp1.gt.ep).or.(tcp2.gt.ep)) then
	   tcp=1.0d0
	end if
        
c  Compute expansion rate.
	if (amnu.eq.0) then
	   rhonu=1.0d0
	   pnu=1.0d0/3.0d0
	else
	   call nu1(a,rhonu,pnu)
	end if
c  8*pi*G*rho*a**2 and 8*pi*G*P*a**2.
	grho=grhom*(omegac+omegab)/a+(grhog+grhor*(nnur+nnunr*rhonu))/a2
     2      +grhom*omegav*a2
	adotoa=sqrt(grho/3.0d0)
	yprime(1)=adotoa*a
c
        gpresg=(grhog/3.0d0)/a2
        gpresn=(grhor*nnur/3.0d0)/a2
c
	ind1=3
	ind2=ind1+lmaxv
	ind3=ind2+lmaxv-1
	ind4=ind3+lmaxv-1
c Vectors
	hv=y(2)
	vb=y(3)
c>>>>>>>>>>>>>>>>>>
cLP vector component of the stress energy
c      if(ak*tau.gt.tkmax) then
c        emtV=0.0d0
c       else

       if (initfl.eq.4) then
       call splint1(ts_f,emV,emVpr,nt_f,tau,emtV)
	else
	 emtV=0.0d0
       endif 
       
c      endif
c<<<<<<<<<<<<<<<<<<
        shg=8.0d0*dsqrt(3.0d0)*y(ind1+2)
        shn=8.0d0*dsqrt(3.0d0)*y(ind2+2)

cLP The sqrt(2) multiplying emtV accounts for left and right handed modes. 
c This means there is no factor of 2 in the final summation of vector C_l 

	hvpr=-2.0d0*adotoa*hv+8.0d0*pi*(2.0d0*dsqrt(2.0d0)*emtV)/ak
     &       -(gpresg*shg+gpresn*shn)/ak
c
c  Photon mass density over baryon mass density.
	photbar=grhog/(grhom*omegab*a)
	pb43=4.0d0/3.0d0*photbar

	vbpr=hvpr-adotoa*(vb-hv)+pb43*xopac*(3.0d0*y(ind1+1)-vb)
c
        yprime(2)=hvpr
	yprime(3)=vbpr
c
c  Polarization term, polv
        polv=0.5d0*(y(ind1+2)-dsqrt(6.0d0)*y(ind3+2))
c
	if (tcp.eq.1) then
c no tight coupling approx
c ----------------Photon perturbations--------------------------
c 1st and 2nd moments
      yprime(ind1+1)=-ak*y(ind1+2)/dsqrt(3.0d0)
     &               -xopac*y(ind1+1)+(xopac*vb+hvpr)/3.0d0
      yprime(ind1+2)=0.2d0*ak*(y(ind1+1)*dsqrt(3.0d0)
     &               -dsqrt(8.0d0)*y(ind1+3))
     &               -xopac*y(ind1+2)+xopac*polv/5.0d0
c l=3...lmaxv-1
	   do 71 l=3,lmaxv-1
c
      yprime(ind1+l)=
     & ak*denl(l)*(dl2m1(l)*y(ind1-1+l)
     & -dl2p2l(l)*y(ind1+1+l))
     & -xopac*y(ind1+l)
71	   continue
c
c Truncate photon moment expansion
      yprime(ind1+lmaxv)=ak*y(ind1-1+lmaxv)
     &     -(lmaxv+1)*y(ind1+lmaxv)/tau
     &	   -xopac*y(ind1+lmaxv)
c ---------------------------------------------------------------
c ----------------Massless Neutrinos-----------------------------
c 1st moment
      yprime(ind2+1)=-ak*y(ind2+2)/dsqrt(3.0d0)+hvpr/3.0d0
c l=2...lmaxv-1
	   do 72 l=2,lmaxv-1
c
      yprime(ind2+l)=
     & ak*denl(l)*(dl2m1(l)*y(ind2-1+l)
     & -dl2p2l(l)*y(ind2+1+l))
72	   continue
c
c Truncate neutrino moment expansion
      yprime(ind2+lmaxv)=ak*y(ind2-1+lmaxv)
     &     -(lmaxv+1)*y(ind2+lmaxv)/tau
c ---------------------------------------------------------------
c --------------- E polarization --------------------------------
c 2nd moment
      yprime(ind3+2)=-ak*(y(ind4+2)/3.0d0+
     &                    dsqrt(40.0d0/9.0d0)*y(ind3+3)*denl(2))
     &               -xopac*y(ind3+2)-xopac*dsqrt(6.0d0)*polv/5.0d0
c l=3...lmaxv-1
	   do 73 l=3,lmaxv-1
c
      yprime(ind3+l)=ak*(dkap(l)*y(ind3+l-1)*denl(l)-
     &                   2.0d0*y(ind4+l)*dlp1(l)-
     &                  dkap1(l)*y(ind3+l+1)*denl(l))
     &               -xopac*y(ind3+l)
73	   continue
c
c Truncate E moment expansion
      yprime(ind3+lmaxv)=ak*y(ind3-1+lmaxv)
     &     -(lmaxv+1)*y(ind3+lmaxv)/tau
     &	   -xopac*y(ind3+lmaxv)
c ---------------------------------------------------------------
c --------------- B polarization --------------------------------
c l=2...lmaxv-1
	   do 74 l=2,lmaxv-1
c
      yprime(ind4+l)=ak*(dkap(l)*y(ind4+l-1)*denl(l)+
     &                   2.0d0*y(ind3+l)*dlp1(l)-
     &                  dkap1(l)*y(ind4+l+1)*denl(l))
     &               -xopac*y(ind4+l)
74	   continue
c
c Truncate B moment expansion
      yprime(ind4+lmaxv)=ak*y(ind4-1+lmaxv)
     &     -(lmaxv+1)*y(ind4+lmaxv)/tau
     &	   -xopac*y(ind4+lmaxv)
c ---------------------------------------------------------------
	else
	   nv=2*lmaxv+2*(lmaxv-1)+3
	   do 80 l=3,nv
	      yprime(l)=0.0d0
80	   continue
           yprime(ind1+1)=hvpr/3.0d0
c
	end if
	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine foutputv(n,y,ypr,j,tau0,tau,dv,dve,dvb)
c
        implicit double precision (a-h,o-z)
        integer j,n
        parameter (nnmax=1)
c COMMON tensor,vector
	parameter (nstep0=6000)
	common /cosmoparm/ ak,ak2,amnu,lmaxg,lmaxnr,lmaxnu,
     &                     nqmax,iq0,iq1,iq2
        common /tensor/ ant(nnmax),rat(nnmax),itflag,lmaxt
        common /vector/ lmaxv
c
        dimension y(n),ypr(n)
        dimension vis(nstep0),dvis(nstep0),ddvis(nstep0)
        dimension opac(nstep0),dopac(nstep0),expmmu(nstep0)
        common /visib/ vis,dvis,ddvis,expmmu,opac,dopac
        save /visib/
c
c Vectors
c
        x=ak*(tau0-tau)
        x2=x*x
	ind1=3
	ind2=ind1+lmaxv
	ind3=ind2+lmaxv-1
	ind4=ind3+lmaxv-1
        hv=y(2)
        hvpr=ypr(2)
        vb=y(3)
        polv=0.5d0*(y(ind1+2)-dsqrt(6.0d0)*y(ind3+2))
        polvpr=0.5d0*(ypr(ind1+2)-dsqrt(6.0d0)*ypr(ind3+2))
c
c

        if (x.gt.0.0d0) then

         dv=expmmu(j)*hvpr/(dsqrt(2.0d0)*x)
     &       +vis(j)*(vb/(dsqrt(2.0d0)*x)+
     &                polvpr*dsqrt(1.5d0)/(ak*x))
     &       +dvis(j)*polv*dsqrt(1.5d0)/(ak*x)

c LP Eliminated minus sign in front of dve and dvb, that
c was in Whu and White and in previous cmbact, to
c conform to conventinos used for scalar modes. This
c affects the sign of vector TE correlation, which should
c be correct now

           dve=(sqrt(6.0d0)/2.0d0)*(2.0d0*vis(j)*polv/x**2
     &         + (vis(j)*polvpr + dvis(j)*polv)/(ak*x)  )


           dvb=(sqrt(6.0d0)/2.0d0)*vis(j)*polv/x
c
        else
           dv=0.0d0
           dve=0.0d0
           dvb=0.0d0
        end if

	return
        end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine finithermo(taumin,taumax,tau0,taurend,
     2      dlntau0,n1,nstep)
c  Compute and save unperturbed baryon temperature and ionization fraction
c  as a function of time.  With nthermo=10000, xe(tau) has a relative 
c accuine foutpuracy (numerical integration precision) better than 1.e-5.

	implicit double precision (a-h,o-z)
c COMMON lingerinc,genparm,reionization
	common /lingerinc/ omegab,omegac,omegav,omegan,omegak,h0
     &                     ,tcmb,yhe,nnur,nnunr
	common /genparm/ grhom,grhog,grhor,adotrad,taurst,dtaurec
c
	parameter (barssc0=9.1820d-14)
	parameter (nthermo=10000)
	parameter (nstep0=6000)

	dimension tb(nthermo),cs2(nthermo),xe(nthermo)
	dimension dcs2(nthermo)
        dimension dotmu(nthermo),sdotmu(nthermo),emmu(nthermo)
        dimension demmu(nthermo),ddotmu(nthermo)
        dimension dddotmu(nthermo),ddddotmu(nthermo)
        dimension vis(nstep0),dvis(nstep0),ddvis(nstep0)
        dimension opac(nstep0),dopac(nstep0),expmmu(nstep0)
        real*8 atau0(nstep0),dtau1(nstep0),dtau2(nstep0)
        common/thermod/dotmu,ddotmu,cs2,dcs2,tauminn,dlntau
        common /visib/ vis,dvis,ddvis,expmmu,opac,dopac
        common /timesteps/ atau0,dtau1,dtau2
        common /reionization/zri,taurist,zristp,tauristp,rif
        common /reionization2/ j2ri1,nri,nri0
c	  save/thermod/
        save /visib/

        ncount=0
	thomc0=5.0577d-8*tcmb**4
        akthom=(2.3048d-9)*(1-yhe)*omegab*h0*h0
	tauminn=0.05d0*taumin
	dlntau=log(tau0/tauminn)/(nthermo-1)
c
c  Initial conditions: assume radiation-dominated universe.
	tau01=tauminn
	adot0=adotrad
	a0=adotrad*tauminn
        a02=a0*a0
c  Assume that any entropy generation occurs before tauminn.
c  This gives wrong temperature before pair annihilation, but
c  the error is harmless.
	tb(1)=tcmb/a0
	xe0=1.0d0
	x1=0.0d0
	x2=1.0d0
	xe(1)=xe0+0.25d0*yhe/(1.0d0-yhe)*(x1+2*x2)
	barssc=barssc0*(1.d0-0.75d0*yhe+(1.d0-yhe)*xe(1))
	cs2(1)=4.0d0/3.0d0*barssc*tb(1)
        dotmu(1)=xe(1)*akthom/a02
        sdotmu(1)=0
c
	  do 10 i=2,nthermo
	  tau=tauminn*exp((i-1)*dlntau)
	  dtau=tau-tau01
c  Integrate Friedmann equation using inverse trapezoidal rule.
	  a=a0+adot0*dtau
	  a2=a*a
	  call nu1(a,rhonu,pnu)
	grho=grhom*(omegac+omegab)/a+(grhog+grhor*(nnur+nnunr*rhonu))/a2
     2      +grhom*omegav*a2
	  adot=sqrt(grho/3.0d0)*a
	  a=a0+2.0d0*dtau/(1.0d0/adot0+1.0d0/adot)
c  Baryon temperature evolution: adiabatic except for Thomson cooling.
c  Use  quadrature solution.
	  tg0=tcmb/a0
	  ahalf=0.5d0*(a0+a)
	  adothalf=0.5d0*(adot0+adot)
c  fe=number of free electrons divided by total number of free baryon
c  particles (e+p+H+He).  Evaluate at timstep i-1 for convenience; if
c  more accuracy is required (unlikely) then this can be iterated with
c  the solution of the ionization equation.
	  fe=(1.d0-yhe)*xe(i-1)/(1.d0-0.75d0*yhe+(1.d0-yhe)*xe(i-1))
	  thomc=thomc0*fe/adothalf/ahalf**3
	  etc=exp(-thomc*(a-a0))
	  a2t=a0*a0*(tb(i-1)-tg0)*etc-tcmb/thomc*(1.d0-etc)
	  tb(i)=tcmb/a+a2t/(a*a)
c  Integrate ionization equation.
	  tbhalf=0.5d0*(tb(i-1)+tb(i))
	  call ionize(tbhalf,ahalf,adothalf,dtau,xe0)
	  call ionhe(tb(i),a,xe0,x1,x2)
c If there is re-ionization, smoothly increase xe to the 
c requested value.
          if ((zri.ne.0.0).and.(tau.gt.(9.0d0*taurist/10.0d0))) then
             if(ncount.eq.0) then
                ncount=i-1
             end if   
             xod=150.0d0*(tau-taurist)/taurist
             if (xod.gt.100) then
                tgh=1
             else
                tgh=(exp(xod)-exp(-xod))/(exp(xod)+exp(-xod))
             end if
             xe(i)=(rif-xe(ncount))*(tgh+1.0d0)/2.0d0+xe(ncount)
          else
             xe(i)=xe0+0.25d0*yhe/(1.0d0-yhe)*(x1+2*x2)
          end if
c  Baryon sound speed squared (over c**2).
	  dtbdla=-2.0d0*tb(i)-thomc*adothalf/adot*(a*tb(i)-tcmb)
	  barssc=barssc0*(1.d0-0.75d0*yhe+(1.d0-yhe)*xe(i))
	  cs2(i)=barssc*tb(i)*(1-dtbdla/tb(i)/3.0d0)
c Calculation of the visibility function
          dotmu(i)=xe(i)*akthom/a2
          if (tau.lt.0.001) then
             sdotmu(i)=0
             go to 15
          end if   
         sdotmu(i)=sdotmu(i-1)+2.0d0*dtau/(1.0d0/dotmu(i)
     2        +1.0d0/dotmu(i-1))
c
 15       a0=a
	  tau01=tau
	  adot0=adot
 10	continue
        do 20 j1=1,nthermo
           emmu(j1)=exp(sdotmu(j1)-sdotmu(nthermo))+1.0d-30
 20     continue   
c
        iv=0
        vfi=0.0d0
c Getting the starting and finishing times for decoupling.
	if (ncount.eq.0) then
	   cf1=1.0d0
	   ns=nthermo
	   else
	      cf1=exp(sdotmu(nthermo)-sdotmu(ncount))
	      ns=ncount
	   end if
        do 30 j1=1,ns
           vfi=vfi+emmu(j1)*dotmu(j1)*cf1
           if ((iv.eq.0).and.(vfi.gt.1.0d-5)) then
              taurst=9.0d0/10.0d0*tauminn*exp(dble(j1-1)*dlntau)
              iv=1
           end if
           if ((iv.eq.1).and.(vfi.gt.0.99)) then
              taurend1=1.5d0*(tauminn*exp(dble(j1-1)*dlntau))
	      taurend1=max(taurend1,taurend1*sqrt(2500.0d0/
     &  	   (omegac+omegab)/h0**2))  
              iv=2
           end if
 30        continue
	   if (iv.ne.2) then
	      taurend1=1.5d0*(tauminn*exp(dble(ncount-1)*dlntau))
	   end if
c Calculating the timesteps during recombination.

	   if (dtaurec.ne.0.0) then
	      dtaurec=min(dtaurec,taurst/40.0d0)
	   else
	      dtaurec=taurst/40.0d0
	   end if

           taurend2=dtaurec/dlntau0
           taurend=max(taurend1,taurend2)
           taurend=min(taurend,9.0d0*taurist/10.0d0)
	   n1=int((taurend-taurst)/dtaurec+1)
           dtaurec=(taurend-taurst)/(n1-1)
           n1=n1+1
c Calculating the timesteps after recombination (logarithmic
c outside re-ionization scattering surface).
           nstep=int(log(taumax/taurend)/dlntau0)
           dlntau0=log(taumax/taurend)/nstep
           nstep=nstep+n1


c Adjusting if there is reionization.
c There will be nri0 points to sample the quick rise in 
c the free electron density. After that, timesteps of length
c dtauri until tauristp.
c
           nri0=50
           if (zri.ne.0) then
              j2ri1=int(log(9.0d0*taurist/taurend/10.0d0)
     &                             /dlntau0+n1)
              j2ri2=int(log(21.0d0*taurist/taurend/20.0d0)
     &                             /dlntau0+n1)
              j2ri3=int(log(tauristp/taurend)/dlntau0+n1)
              dtri0=taurend*(exp(dlntau0*(j2ri2-n1))-
     2            exp(dlntau0*(j2ri1-n1)))
              dtri=taurend*(exp(dlntau0*(j2ri3-n1))-
     2            exp(dlntau0*(j2ri2-n1)))
	      dtauri0=dtri0/nri0
	      dtauri=1.5d0*dtaurec
	      dtauri=min(dtauri,dtri/10.0d0)
	      nri=int(dtri/dtauri)
	      dtauri=dtri/nri
              nstep=nstep+nri0+nri+j2ri1-j2ri3
              else
                 j2ri1=0
                 j2ri2=0
	   end if
c
	   if (nstep.gt.nstep0) then
	      write(*,*) 
     2   	'Sorry, the arrays were dimensioned for a maximum of'
              write(*,*)nstep0, 'timesteps.'
              write(*,*)'The model you requested needs'
              write(*,*)nstep,'Please make the arrays bigger by making'
              write(*,*)'nstep0 bigger where it appears'
	      stop
	   end if

	call splini
	call splder(cs2,dcs2,nthermo)
        call splder(dotmu,ddotmu,nthermo)
        call splder(ddotmu,dddotmu,nthermo)
        call splder(dddotmu,ddddotmu,nthermo)
        call splder(emmu,demmu,nthermo)
c
c Saving tau values and the quantities needed to calculate
c the derivatives of the visibility function appearing in the sources.
        do 40 j2=2,nstep
        if (j2.le.n1) then
           tau=taurst+dble(j2-2)*dtaurec
           else
              if ((zri.eq.0).or.(j2.le.j2ri1)) then
                 tau=taurend*exp(dlntau0*dble(j2-n1))
                 else
                    if (j2.lt.(j2ri1+nri+nri0)) then
		       if (j2.le.(j2ri1+nri0)) then
                       tau=atau0(j2ri1)+dtauri0*dble(j2-j2ri1)
		       else
			  tau=atau0(j2ri1+nri0)+dtauri
     &                          *dble(j2-j2ri1-nri0)
		       end if
                       else
                          tau=taurend*exp(dlntau0*dble(j2-j2ri1-nri0
     2                        -nri+j2ri3-n1))
                    end if
              end if      
        end if 
        atau0(j2)=tau
c Cubic-spline interpolation.
           d=log(tau/tauminn)/dlntau+1.0d0
           i=int(d)
           d=d-i
           if (i.lt.nthermo) then
	  opac(j2)=dotmu(i)+d*(ddotmu(i)+d*(3.0d0*(dotmu(i+1)-dotmu(i))
     2         -2.0d0*ddotmu(i)-ddotmu(i+1)+d*(ddotmu(i)+ddotmu(i+1)
     3         +2.0d0*(dotmu(i)-dotmu(i+1)))))
	  dopac(j2)=(ddotmu(i)+d*(dddotmu(i)+d*(3.0d0*(ddotmu(i+1)
     2         -ddotmu(i))-2.0d0*dddotmu(i)-dddotmu(i+1)+d*(dddotmu(i)
     3         +dddotmu(i+1)+2.0d0*(ddotmu(i)-ddotmu(i+1))))))/(tau
     4         *dlntau)
	  ddopac=(dddotmu(i)+d*(ddddotmu(i)+d*(3.0d0*(dddotmu(i+1)
     2         -dddotmu(i))-2.0d0*ddddotmu(i)-ddddotmu(i+1)
     3         +d*(ddddotmu(i)+ddddotmu(i+1)+2.0d0*(dddotmu(i)
     4         -dddotmu(i+1)))))-(dlntau**2)*tau*dopac(j2))
     5         /(tau*dlntau)**2
	  expmmu(j2)=emmu(i)+d*(demmu(i)+d*(3.0d0*(emmu(i+1)-emmu(i))
     2         -2.0d0*demmu(i)-demmu(i+1)+d*(demmu(i)+demmu(i+1)
     3         +2.0d0*(emmu(i)-emmu(i+1)))))

          vis(j2)=opac(j2)*expmmu(j2)
          dvis(j2)=expmmu(j2)*(opac(j2)**2+dopac(j2))
          ddvis(j2)=expmmu(j2)*(opac(j2)**3+3*opac(j2)*dopac(j2)+ddopac)
          else
          opac(j2)=dotmu(nthermo)
          dopac(j2)=ddotmu(nthermo)
          ddopac=dddotmu(nthermo)
          expmmu(j2)=emmu(nthermo)
          vis(j2)=opac(j2)*expmmu(j2)
          dvis(j2)=expmmu(j2)*(opac(j2)**2+dopac(j2))
          ddvis(j2)=expmmu(j2)*(opac(j2)**3+3.0d0*opac(j2)
     2        *dopac(j2)+ddopac)
c
          end if
  40       continue
          atau0(1)=0.0d0
          atau0(nstep)=min(tau0,atau0(nstep))
c saving the length of the timesteps for the time integration.
          do 50 j2=2,nstep-1
             dtau1(j2)=atau0(j2+1)-atau0(j2)
             dtau2(j2)=abs(atau0(j2+1)-atau0(j2-1))/2.0d0
 50       continue

          dtau1(nstep)=dtau1(nstep-1)
          dtau2(2)=dtau1(2)/2.0d0
          dtau2(nstep)=(atau0(nstep)-atau0(nstep-1))/2.0d0
	return
	end	

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine finitial(y,tau)
c  Initial conditions.
	implicit double precision (a-h,o-z)
        integer initfl

      common/volcom/tkmax
      common/picom/pi,twopi3
c      common/maxvcom/maxvar
c
	common /lingerinc/ omegab,omegac,omegav,omegan,omegak,h0
     &                     ,tcmb,yhe,nnur,nnunr
	common /cosmoparm/ ak,ak2,amnu,lmax,lmaxnr,lmaxnu,
     &                     nqmax,iq0,iq1,iq2
	common /genparm/ grhom,grhog,grhor,adotrad,taurst,dtaurec
	common /initcase/ initfl
c
	parameter (lmax0=8,lmaxnr0=25,lmaxnu0=25,nqmax0=15)
	
c    	parameter (nvar0=1+7+2*(lmax0+1)+(lmaxnr0+1)+nqmax0*(lmaxnu0+1))
        parameter (nvar0=7+2*(lmax0+1)+(lmaxnr0+1)+nqmax0*(lmaxnu0+1))
	dimension y(nvar0)
c
	a=tau*adotrad
	a2=a*a
	call nu1(a,rhonu,pnu)
c  8*pi*G*rho*a**2 and 8*pi*G*P*a**2.
	grho=grhom*(omegac+omegab)/a+(grhog+grhor*(nnur+nnunr*rhonu))/a2
     2      +grhom*omegav*a2
c	adotoa=sqrt(grho/3.0d0)
	gpres=((grhog+grhor*nnur)/3.0d0+grhor*nnunr*pnu)/a2
     2       -grhom*omegav*a2
	s=grho+gpres
	fracnu=grhor*4.0d0/3.0d0*(nnur+nnunr)/a2/s
c  Use yrad=rho_matter/rho_rad to correct initial conditions for
c  matter+radiation.
	yrad=grhom*(omegac+omegab)*a/(grhog+grhor*(nnur+nnunr*rhonu))
c
c  Choose one of the following four cases for initial conditions, or
c  add your own.  Comment out the other cases.
c
	if (initfl.eq.1) then
c-------------------------------------------------------------------------------
c  First case.
c  Isentropic ("adiabatic") initial conditions.
	psi=-1.0d0
	C=(15.0d0+4.0d0*fracnu)/20.0d0*psi
	akt2=(ak*tau)**2
	h=C*akt2*(1.0d0-0.2d0*yrad)
	eta=2.0d0*C-(5.0d0+4.0d0*fracnu)/6.0d0/(15.0d0+4.0d0*fracnu)*
     2              C*akt2*(1.0d0-yrad/3.0d0)
	f1=(23.0d0+4.0d0*fracnu)/(15.0d0+4.0d0*fracnu)
	deltac=-0.5d0*h
	deltag=-2.0d0/3.0d0*h*(1.0d0-akt2/36.0d0)
	deltab=0.75d0*deltag
	deltar=-2.0d0/3.0d0*h*(1.0d0-akt2/36.0d0*f1)
	thetac=0.0d0
	thetag=-C/18.0d0*akt2*akt2/tau
	thetab=thetag
	thetar=f1*thetag
	shearr=4.0d0/15.0d0*ak2/s*psi*(1.0d0+7.0d0/36.0d0*yrad)
	ahdot=2.0d0*C*ak2*tau*a*(1.0d0-0.3d0*yrad)
c-------------------------------------------------------------------------------
	else if (initfl.eq.2) then
c  Second case.
c  Isocurvature CDM initial conditions: perturb only CDM as a --> 0.
	delta0=1.0d0
	h=delta0*yrad*(1.0d0/(1.0d0+omegab/omegac)-0.5d0*yrad)
	deltac=delta0-0.5d0*h
c  Compensate perturbation with everything else.
	deltag=-2.0d0/3.0d0*h
	deltab=0.75d0*deltag
	deltar=deltag
	thetac=0.0d0
	thetag=-h/12.0d0*ak2*tau
	thetab=thetag
	thetar=thetag
	shearr=0.0d0
	ahdot=adotrad*h*(1.0d0-0.5d0*yrad)
	eta=-h/6.0d0
c-------------------------------------------------------------------------------
	else if (initfl.eq.3) then
c  Third case.
c  Isocurvature baryon initial conditions: perturb only baryons as a->0.
	delta0=1.0d0
	h=delta0*yrad*(1.0d0/(1.0d0+omegac/omegab)-0.5d0*yrad)
	deltab=delta0-0.5d0*h
c  Compensate perturbation with everything else.
	deltac=-0.5d0*h
	deltag=-2.0d0/3.0d0*h
	deltar=deltag
	thetac=0.0d0
	thetag=-h/12.0d0*ak2*tau
	thetab=thetag
	thetar=thetag
	shearr=0.0d0
	ahdot=adotrad*h*(1.0d0-0.5d0*yrad)
	eta=-h/6.0d0
c-------------------------------------------------------------------------------
	else if (initfl.eq.4) then
c  Fourth case.
c  Isocurvature seed initial conditions:everything is unperturned as a->0
	delta0=1.0d0
	h=0.0d0
c  Compensate perturbation with everything else.
	deltab=-0.5d0*h
	deltac=-0.5d0*h
	deltag=-2.0d0/3.0d0*h
	deltar=deltag
	thetac=0.0d0
	thetag=-h/12.0d0*ak2*tau
	thetab=thetag
	thetar=thetag
	shearr=0.0d0
	ahdot=adotrad*h*(1.0d0-0.5d0*yrad)
	eta=-h/6.0d0
c-------------------------------------------------------------------------------
	else
	  write(*,*) 'initfl must equal 1-4! initfl=',initfl
	  stop
	end if
c
	deltan=deltar
	thetan=thetar
c
c      y(maxvar)=0.0d0
c
	y(1)=a
	y(2)=ahdot
	y(3)=eta
c  CDM.
	y(4)=deltac
	y(5)=thetac
c  Baryons.
	y(6)=deltab
	y(7)=thetab
c  Photons (total intensity and polarization).
	y(8)=deltag
	y(9)=thetag
c	shearg=0.0d0
	y(9+lmax)=0.0d0
	y(10+lmax)=0.0d0
	  do 10 l=2,lmax
	  y(8+l)=0.0d0
	  y(9+lmax+l)=0.0d0
10	continue
c  Massless neutrinos.
	y(10+2*lmax)=deltar
	y(11+2*lmax)=thetar
	y(12+2*lmax)=shearr*2.0d0
	  do 20 l=3,lmaxnr
	  y(10+2*lmax+l)=0.0d0
20	continue
c  Massive neutrinos.
	if (nqmax.eq.0) go to 50
c	dq=qmax/nqmax
	dq=1.0d0
	  do 40 i=1,nqmax
	  q=i*dq-0.5d0
	  aq=a*amnu/q
	  v=1.0d0/sqrt(1.0d0+aq*aq)
	  akv=ak*v
	  expq=exp(-q)
	  dlfdlq=-q/(1.0d0+expq)
	  y(iq0+i-1)=-0.25d0*dlfdlq*deltan
c  Divide by v to get first-order correction for neutrino mass.
	  y(iq1+i-1)=-dlfdlq*thetan/akv/3.0d0
	  y(iq2+i-1)=-0.5d0*dlfdlq*shearr
	    do 30 l=3,lmaxnu
	    ind=10+2*lmax+lmaxnr+i+l*nqmax
	    y(ind)=0.0d0
30	  continue
40	continue
c  Check energy constraint equation.
50	call nu2(a,drhonu,fnu,dpnu,shearnu,y(iq0),y(iq1),y(iq2))
	deltan=drhonu/rhonu
	thetan=ak*fnu/(rhonu+pnu)
c	shearn=shearnu/(rhonu+pnu)
	dgrho=grhom*(omegac*deltac+omegab*deltab)/a
     2    +(grhog*deltag+grhor*(nnur*deltar+nnunr*drhonu))/a2
c  Add a seed if desired.
c	if (initfl.eq.4) dgrho=dgrho+grhom/a
c	econ=(adotoa*ahdot/a-2.0d0*ak2*eta-dgrho)/grho
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine finitialt(y,tau)
c  Initial conditions.
	implicit double precision (a-h,o-z)
c COMMON genparm,tensor
        parameter (nnmax=1)
	common /genparm/ grhom,grhog,grhor,adotrad,taurst,dtaurec
        common /tensor/ ant(nnmax),rat(nnmax),itflag,lmaxt
c
	parameter (lmaxt0=10)
    	parameter (nvar0t=2*(lmaxt0+1)+2+1)
c
	dimension y(nvar0t)
c
	a=tau*adotrad
c
	y(1)=a
c Tensor modes
	y(2)=0.0d0
	y(3)=0.0d0
	ind1=4
	ind2=ind1+lmaxt
	do l=0,lmaxt
	   y(ind1+l)=0.0d0
	   y(ind2+l)=0.0d0
	end do
	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine finitialv(y,tau)
c  Initial conditions.
	implicit double precision (a-h,o-z)
c COMMON genparm,vector,tensor
        parameter (nnmax=1)
	common /genparm/ grhom,grhog,grhor,adotrad,taurst,dtaurec
        common /tensor/ ant(nnmax),rat(nnmax),itflag,lmaxt
        common /vector/ lmaxv
c
	parameter (lmaxv0=8)
    	parameter (nvar0v=2*lmaxv0+2*(lmaxv0-1)+2+1)
c
	dimension y(nvar0v)
c
	a=tau*adotrad
c
	y(1)=a
c Vector modes
	y(2)=0.0d0
	y(3)=0.0d0
c	
        nv=2*lmaxv+2*(lmaxv-1)+3
c        
	do l=4,nv
	   y(l)=0.0d0
	end do
c
	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine initjl

c This subroutine reads the jl files from disk and 
c initializes other variables needed in CMBFAST.
  
      implicit double precision(a-h,o-z)
      parameter (l0max=3300,lmax=20+l0max/50,ketamax=2*l0max+126)
      parameter (d0hi=1.0d40,d0lo=1.0d40,xlimmin=10.0d0)
c COMMON lvalues1,lvalues2,jlgen
      real*8 x
      real*4 aj
      real*8 ajl(ketamax,lmax),ajlpr(ketamax,lmax)
      real*8 xx(ketamax), dxx(ketamax)
      integer l(lmax),l0
c      character*80 filename

      common /lvalues1/ l,l0,lmo
      common /lvalues2/ akmax0
      save /lvalues1/
      save /lvalues2/

      common /jlgen/ ajl,ajlpr,xx,dxx,kmax
      save /jlgen/

c      write(*,*)'Enter input filename for jl'
c      read(*,'(a80)')filename
       open(unit=10,file='jlgen.dat'
     &	   ,status='unknown',form='unformatted')
      rewind 10
      read(10)lmofile	
      read(10)kmaxfile

      if ((lmo.gt.lmofile).or.(akmax0.gt.kmaxfile)) then
         write(*,*)'You have entered a lmax and/or kmax'
         write(*,*)'inconsistent with those in the file'
         write(*,*)lmofile,kmaxfile
         write(*,*)'You will have to start again'
         stop
      end if

      kmaxfile=kmaxfile-25+151
      kmax=int(akmax0-25+151)


c
c reading  j_l 
c remember to create jl.dat with jlgen.f first using 
c correct lmax and ketamax

      do 40 i=1,kmaxfile
         if (i.le.151) then
            if (i.le.51) then
               xx(i)=dble((i-1))/10.0d0
            else
               xx(i)=dble((i-51))/5.0d0+5.0d0 
            end if
         else
            xx(i)=dble((i-151))+25.0d0
         end if
 40   continue

      do 50 j=1,l0
         do 60 i=1,kmaxfile
            x=xx(i)
            xlim=0.05d0*dble(l(j))
            xlim=max(xlim,xlimmin)
            xlim=l(j)-xlim
            if (x.gt.xlim) then
               read(10)aj
               ajl(i,j)=dble(aj)
            else
               ajl(i,j)=0.0d0
            end if
 60      continue
 50   continue
      close(10)

      do 70 i=2,(kmaxfile-1)
         dxx(i)=(xx(i+1)-xx(i-1))/2.0d0
 70   continue
      dxx(1)=xx(2)/2.0d0
      dxx(kmaxfile)=(xx(kmaxfile)-xx(kmaxfile-1))/2.0d0
c
c get the interpolation matrix for bessel functions

      do 80 j=1,l0
         call spline(xx,ajl(1,j),kmaxfile,d0lo,d0hi,ajlpr(1,j))
 80   continue

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine powersflat(ak,in,apower)

c This subroutine computes the power spectra
c for mode ak of the scalar perturbations. 
c Now it is set to a power law.

      implicit double precision(a-h,o-z)
      parameter (nnmax=1)
c COMMON initialps
      common /initialps/ an(nnmax),nn

      apower=exp((an(in)-2.0d0)*log(ak))

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine powertflat(ak,in,apower)

c This subroutine computes the power spectra
c for mode ak of the scalar perturbations. 
c Now it is set to a power law.
  
      implicit double precision(a-h,o-z)
      parameter (nnmax=1)
c COMMON tensor
      common /tensor/ ant(nnmax),rat(nnmax),itflag,lmaxt

      apower=exp((ant(in)-1.0d0)*log(ak))

      return
      end


	subroutine COBEnormalize(igmu,gmu,clts,cltt,cltv,cl2av,
     & clerr,cles,clev,clbv,clet,clbt,clcs,clct,clcv,clbtv,
     & clbtt,clebv,clebt,trsum,akout)
c
	implicit double precision(a-h,o-z)
c COMMON lvalues1,tensor,vector,transfer,initialps,trfun
c COMMON lingerinc
	parameter (lmax0=3300,lmax=20+lmax0/50)
	parameter (sigr0=8.0d0,nnmax=1)
      parameter (fourpi=4.0d0*3.14159265d0)
c
	real*8 ztf(nnmax)
c
        real*8 clts(lmax0,nnmax),cltt(lmax0,nnmax)
        real*8 cltv(lmax0,nnmax),clerr(lmax0,nnmax)
        real*8 cl2av(lmax0,nnmax)
        real*8 cles(lmax0,nnmax),clet(lmax0,nnmax)
        real*8 clbt(lmax0,nnmax)
        real*8 clcs(lmax0,nnmax),clct(lmax0,nnmax)
        real*8 clev(lmax0,nnmax),clbv(lmax0,nnmax)
        real*8 clcv(lmax0,nnmax)
        real*8 clbtv(lmax0,nnmax),clbtt(lmax0,nnmax)
        real*8 clebv(lmax0,nnmax),clebt(lmax0,nnmax)
c
	integer l(lmax),l0
	real*8 xl(lmax)
c
	common /lvalues1/ l,l0,lmo
	save /lvalues1/
c
	common/trfun/nkount
      parameter (nk0=200)
      dimension akout(nk0),trsum(nk0)
c
        common /initialps/ an(nnmax),nn
        common /tensor/ ant(nnmax),rat(nnmax),itflag,lmaxt
        common /vector/ lmaxv
        common /transfer/akmaxt,ztf,nlnkt,ict,ntf
	common /lingerinc/ omegab,omegac,omegav,omegan,omegak,h0
     &                     ,tcmb,yhe,nnur,nnunr
c
	xlog10=log(10.0d0)
	xlnh=log(h0/100.0d0)
	h=h0/100.0d0

	do j=1,l0
	   xl(j)=dble(l(j))
	end do
c Curvature radius
        if (omegak.ne.0.0d0) then
         hc=2.998d5/h0
         curv=-omegak/hc/hc
         r=1.0d0/sqrt(abs(curv))
        endif
c
c COBE normalization
c fit the spectrum to a quadratic around C_10 with equal weights in logl

	do 700 in=1,nn

       if(igmu.eq.1) then

         c10=gmu**2/((0.5d0*fourpi)**3*1.1d-9)

       else


	   c10=clts(10,in)+cltt(10,in)+cltv(10,in)

	   d1=(clts(int(xl(2)),in)+cltt(int(xl(2)),in)
     &	   +cltv(int(xl(2)),in))/c10-1.0d0
	   d2=(clts(int(xl(3)),in)+cltt(int(xl(3)),in)
     &	   +cltv(int(xl(3)),in))/c10-1.0d0
	   d3=(clts(int(xl(5)),in)+cltt(int(xl(5)),in)
     &	   +cltv(int(xl(5)),in))/c10-1.0d0
	   d4=(clts(int(xl(7)),in)+cltt(int(xl(7)),in)
     &	   +cltv(int(xl(7)),in))/c10-1.0d0
	   d5=(clts(int(xl(10)),in)+cltt(int(xl(10)),in)
     &	   +cltv(int(xl(10)),in))/c10-1.0d0
	   d6=(clts(int(xl(11)),in)+cltt(int(xl(11)),in)
     &	   +cltv(int(xl(11)),in))/c10-1.0d0
	   d7=(clts(int(xl(12)),in)+cltt(int(xl(12)),in)
     &	   +cltv(int(xl(12)),in))/c10-1.0d0

	   x1=log(xl(2))/xlog10-1.0d0
           x2=log(xL(3))/xlog10-1.0d0
           x3=log(xl(5))/xlog10-1.0d0
           x4=log(xl(7))/xlog10-1.0d0
           x5=log(xl(10))/xlog10-1.0d0
           x6=log(xl(11))/xlog10-1.0d0
           x7=log(xl(12))/xlog10-1.0d0
           sy=x1*d1+x2*d2+x3*d3+x4*d4+x5*d5+x6*d6+x7*d7
           s=x1*x1+x2*x2+x3*x3+x4*x4+x5*x5+x6*x6+x7*x7
           sx=x1**3+x2**3+x3**3+x4**3+x5**3+x6**3+x7**3
           sxy=x1**2*d1+x2**2*d2+x3**2*d3+x4**2*d4+
     2        x5**2*d5+x6**2*d6+x7**2*d7
           sxx=x1**4+x2**4+x3**4+x4**4+x5**4+x6**4+x7**4
           delt=s*sxx-sx*sx
           d1pr=(sxx*sy-sx*sxy)/delt
           d1ppr=2.0d0*(s*sxy-sx*sy)/delt

c Bunn and White fitting formula
           c10=(0.64575d0+0.02282d0*d1pr+0.01391d0*d1pr*d1pr
     2     -0.01819d0*d1ppr-0.00646d0*d1pr*d1ppr
     2     +0.00103d0*d1ppr*d1ppr)/c10
c logl
           xlogl=-0.01669d0+1.19895d0*d1pr-0.83527d0*d1pr*d1pr
     2           -0.43541d0*d1ppr-0.03421d0*d1pr*d1ppr
     2           +0.01049d0*d1ppr*d1ppr
c	   write(19,*)'COBE Likelihood relative to flat=',exp(xlogl)

       endif

c density power spectrum normalization;
c Delta2=4pik^3P(k)=d2norm*k^{n+3}*tf^2
c correct for h*k instead of k

           if (ict.eq.2) then

ccc             d2norm=c10*1.1d-9*fourpi*h**3/(0.5d0*fourpi)**3

        d2norm=c10*1.1d-9*fourpi*h**3

c
           do ik=1,nkount
             ak=akout(ik)
             trsum(ik)=trsum(ik)*d2norm*ak**4
             akout(ik)=ak/h
           enddo
           endif
c C_l normalization; output l(l+1)C_l/twopi
ccccccc write(19,*)'norm',c10,sqrt(c10*1.1d-9)

c          write(19,*)'norm',c10,sqrt((0.5d0*fourpi)**3*c10*1.1d-9)

          c10=c10*2.2d-9/fourpi
	   do il=2,l(l0)
	   clsum2=(clts(il,in)+cltv(il,in)+cltt(il,in))**2
           if((cl2av(il,in)-clsum2).le.0.0d0) goto 699
           clerr(il,in)=sqrt(cl2av(il,in)-clsum2)*c10
c
699      continue
            clts(il,in)=clts(il,in)*c10*2.725d6**2
            cles(il,in)=cles(il,in)*c10*2.725d6**2

            cltv(il,in)=cltv(il,in)*c10*2.725d6**2
            clev(il,in)=clev(il,in)*c10*2.725d6**2
            clbv(il,in)=clbv(il,in)*c10*2.725d6**2
            clcv(il,in)=clcv(il,in)*c10*2.725d6**2
            clebv(il,in)=clebv(il,in)*c10*2.725d6**2
            clbtv(il,in)=clbtv(il,in)*c10*2.725d6**2
c
            cltt(il,in)=cltt(il,in)*c10*2.725d6**2
            clet(il,in)=clet(il,in)*c10*2.725d6**2
            clbt(il,in)=clbt(il,in)*c10*2.725d6**2
            clebt(il,in)=clebt(il,in)*c10*2.725d6**2
            clbtt(il,in)=clbtt(il,in)*c10*2.725d6**2
c
            clcs(il,in)=clcs(il,in)*c10*2.725d6**2
            clct(il,in)=clct(il,in)*c10*2.725d6**2

           end do
700       continue
c
	   return
	   end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine savetime
      implicit double precision (a-h,o-z)
      parameter (lmaxv0=8)
c COMMON vectime
      common/vectime/dl2m1(lmaxv0),dl2p2l(lmaxv0),dkap(lmaxv0),
     &               dkap1(lmaxv0),dlp1(lmaxv0)
c
      lend=lmaxv0
      do 10 l=1,lend
c
      rl=dble(l)
      rl1=dble(l+1)
c
      dl2m1(l)=dsqrt(rl**2-1.0d0)
      dl2p2l(l)=dsqrt(rl**2+2.0d0*rl)
c
      dkap(l)=dsqrt((rl**2-1.0d0)*(rl**2-4.0d0))/rl
      dkap1(l)=dsqrt((rl1**2-1.0d0)*(rl1**2-4.0d0))/rl1
c
      dlp1(l)=1.0d0/(rl*rl1)
c
10    continue
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine initlval

c This subroutines initializes l arrays.

	implicit double precision(a-h,o-z)
	parameter (l0max=3300,lmax=20+l0max/50)

	integer l(lmax),l0
c COMMON lvalues1,lvalues2
	common /lvalues1/ l,l0,lmo
	common /lvalues2/ akmax0
	save /lvalues1/
	save /lvalues2/


c	write(*,*)'Value of lmax, ketamax (e.g. 1500 3000)'
c	write(*,*)'Remember to be consistent with the file'
c	write(*,*)'in the flat case.'
c	read(*,*)lmo,akmax0
      lmo=3000
      akmax0=6000.0d0

c have l's in l(j); from 2-10,12,15,20,30,40,60,80,100, every 50 beyond
c This are l's which will be calculated, the rest will be interpolated
	lind=1
	do 22 lvar=2,10
	   l(lind)=lvar
	   lind=lind+1
 22	continue
	l(lind)=12
	lind=lind+1
	l(lind)=15
	lind=lind+1
	l(lind)=20
	lind=lind+1
	l(lind)=30
	lind=lind+1
	l(lind)=40

	lind=lind+1
	l(lind)=50
	lind=lind+1
	l(lind)=60
	lind=lind+1
	l(lind)=70
	lind=lind+1
	l(lind)=80
	lind=lind+1
	l(lind)=90
	lind=lind+1
	l(lind)=110
	lind=lind+1
	l(lind)=130
	do 24 lvar=150,lmo,50
	   lind=lind+1
	   l(lind)=lvar
 24	continue
	l0=lind


	return 
	end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine thermo(tau,xcs2,xopac)
c  Compute unperturbed baryon temperature, sound speed squared,
c  and ionization fraction by interpolating pre-computed tables.

	implicit double precision (a-h,o-z)

	parameter (nthermo=10000)
	dimension cs2(nthermo),dcs2(nthermo)
      dimension dotmu(nthermo),ddotmu(nthermo)
	common/thermod/dotmu,ddotmu,cs2,dcs2,tauminn,dlntau
c	save/thermod/

	d=log(tau/tauminn)/dlntau+1.0d0

	i=int(d)
	d=d-i
	if (i.lt.1) then
c  Linear interpolation if out of bounds (should not occur).
	  xcs2=cs2(1)+(d-1)*dcs2(1)
	  xopac=dotmu(1)+(d-1)*ddotmu(1)
	else if (i.ge.nthermo) then
	  xcs2=cs2(nthermo)+(d-nthermo)*dcs2(nthermo)
	  xopac=dotmu(nthermo)+(d-nthermo)*ddotmu(nthermo)
	else
c  Cubic spline interpolation.
	  xcs2=cs2(i)+d*(dcs2(i)+d*(3.0d0*(cs2(i+1)-cs2(i))
     2         -2.0d0*dcs2(i)-dcs2(i+1)+d*(dcs2(i)+dcs2(i+1)
     3         +2.0d0*(cs2(i)-cs2(i+1)))))
	  xopac=dotmu(i)+d*(ddotmu(i)+d*(3.0d0*(dotmu(i+1)-dotmu(i))
     2         -2.0d0*ddotmu(i)-ddotmu(i+1)+d*(ddotmu(i)+ddotmu(i+1)
     3         +2.0d0*(dotmu(i)-dotmu(i+1)))))

       end if
c
	return
	end	
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine ionize(tempb,a,adot,dtau,xe)
c  Integrate the ionization fraction xe for hydrogen semi-implicitly
c  from tau to tau+dtau, treating tempb, a, and adot as constants
c  during this interval.
c
	implicit double precision (a-h,o-z)
c COMMON lingerinc
	common /lingerinc/ omegab,omegac,omegav,omegan,omegak,h0,
     &                     tcmb,yhe,nnur,nnunr
c  Ionization temperature and coefficient.
	parameter (tion=1.5789d5,beta0=43.082d0)
c  Two-photon decay rate (in 1/Mpc).
	parameter (dec2g=8.468d14)
c  Switch for implicit (switch=1.0) or semi-implicit (switch=0.5) scheme.
	parameter (switch=0.5d0)
c
c  Recombination coefficient (in sqrt(K)/Mpc).
	alpha0=2.3866d-6*(1-yhe)*omegab*h0*h0
c  Coefficient for correction of radiative decay (dimensionless).
	crec=8.0138d-26*(1-yhe)*omegab*h0*h0
c  Recombination and ionization rates.
	phi2=0.448d0*log(tion/tempb)
	phi2=max(phi2,0.0d0)
	alpha=alpha0/sqrt(tempb)*phi2/(a*a*a)
	beta=tempb*phi2*exp(beta0-tion/tempb)
c  Peebles' correction factor.
	if (tempb.lt.200.0d0) then
	  cpeebles=1.0d0
	else
	  cp1=crec*dec2g*(1.0d0-xe)/(a*adot)
	  cp2=crec*tempb*phi2*exp(beta0-0.25d0*tion/tempb)*
     2        (1.0d0-xe)/(a*adot)
	  cpeebles=(1.0d0+cp1)/(1.0d0+cp1+cp2)
	end if
c  Integrate dxe=bb*(1-xe)-aa*xe*xe by averaging rhs at current tau
c  (fraction 1-switch) and future tau (fraction switch).
	aa=a*dtau*alpha*cpeebles
	bb=a*dtau*beta*cpeebles
	b1=1.0d0+switch*bb
	bbxe=bb+xe-(1.0d0-switch)*(bb*xe+aa*xe*xe)
	rat=switch*aa*bbxe/(b1*b1)
c  Prevent roundoff error.
	if (rat.lt.1.0d-6) then
	  xe=bbxe/b1*(1.0d0-rat)
	else
	  xe=b1/(2.0d0*switch*aa)*(sqrt(4.0d0*rat+1.0d0)-1.0d0)
	end if
c
	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine ionhe(tempb,a,x0,x1,x2)
c  Compute the helium ionization fractions using the Saha equation.
c  x0 is the hydrogen ionization fraction n(H+)/n(H) (input),
c  x1 is the helium first ionization fraction n(He+)/n(He)
c    (input and output), and
c  x2 is the helium second ionization fraction n(He++)/n(He)
c    (input and output).
c
	implicit double precision (a-h,o-z)
c COMMON lingerinc
	common /lingerinc/ omegab,omegac,omegav,omegan,omegak,h0,
     &                     tcmb,yhe,nnur,nnunr
c  Ionization temperatures.
	parameter (tion1=2.855d5,tion2=6.313d5)
c
c  Constant for electron partition function per baryon.
	b0=2.150d24/((1.0d0-yhe)*omegab*h0*h0)
c  Electron partition function per baryon.
	b=b0*a*a*a*tempb*sqrt(tempb)
c  Dimensionless right-hand sides in Saha equations.
	if (abs(tion1/tempb).lt.100) then
	   r1=4.0d0*b*exp(-tion1/tempb)
	else
	   r1=0.0d0
	end if
	if (abs(tion2/tempb).lt.150) then
	   r2=b*exp(-tion2/tempb)
	else
	   r2=0.0d0
	end if
c
c  Solve coupled equations iteratively.
	c=0.25d0*yhe/(1.0d0-yhe)
	err=1.0d0
	niter=0
10	niter=niter+1
	if (err.lt.1.0d-12) return
	  xe=x0+c*(x1+2*x2)
	  x2new=r1*r2/(r1*r2+xe*r1+xe*xe)
	  x1=xe*r1/(r1*r2+xe*r1+xe*xe)
	  err=abs(x2new-x2)
	  x2=x2new
	  go to 10
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine nu1(a,rhonu,pnu)
c  Compute massive neutrino density and pressure in units of the mean
c  density of one flavor of massless neutrinos.  Use cubic splines to
c  interpolate from a table.
c
	implicit double precision (a-h,o-z)
c
	common /cosmoparm/ ak,ak2,amnu,lmax,lmaxnr,lmaxnu,
     &                            nqmax,iq0,iq1,iq2
	parameter (nrhopn=10000)
	dimension r1(nrhopn),p1(nrhopn)
	dimension dr1(nrhopn),dp1(nrhopn),ddr1(nrhopn)
	common /nu1d/ amin,dlna,r1,p1,dr1,dp1,ddr1
	save /nu1d/

	if (amnu.eq.0.0d0) then
	  rhonu=1.0d0
	  pnu=1.0d0/3.0d0
	  return
	end if

	d=log(a/amin)/dlna+1.0d0
	i=int(d)
	d=d-i
	if (i.lt.1) then
c  Use linear interpolation, bounded by results for massless neutrinos.
	  rhonu=r1(1)+(d-1)*dr1(1)
	  pnu=p1(1)+(d-1)*dp1(1)
	  rhonu=min(exp(rhonu),1.0d0)
	  pnu=min(exp(pnu),0.3333333333d0)
	else if (i.ge.nrhopn) then
c  This should not happen, unless the user evolves to z<0!
	  rhonu=r1(nrhopn)+(d+i-nrhopn)*dr1(nrhopn)
	  pnu=p1(nrhopn)+(d+i-nrhopn)*dp1(nrhopn)
	  rhonu=exp(rhonu)
	  pnu=exp(pnu)
	else
c  Cubic spline interpolation.
	  rhonu=r1(i)+d*(dr1(i)+d*(3.0d0*(r1(i+1)-r1(i))-2.0d0*dr1(i)
     2         -dr1(i+1)+d*(dr1(i)+dr1(i+1)+2.0d0*(r1(i)-r1(i+1)))))
	  pnu=p1(i)+d*(dp1(i)+d*(3.0d0*(p1(i+1)-p1(i))-2.0d0*dp1(i)
     2         -dp1(i+1)+d*(dp1(i)+dp1(i+1)+2.0d0*(p1(i)-p1(i+1)))))
	  rhonu=exp(rhonu)
	  pnu=exp(pnu)
	end if

c
	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine initnu1(amnu)
c  Initialize interpolation tables for massive neutrinos.
c  Use cubic splines interpolation of log rhonu and pnu vs. log a.
	implicit double precision (a-h,o-z)
c
	parameter (nrhopn=10000,nqmax0=15)
	dimension r1(nrhopn),p1(nrhopn)
	dimension dr1(nrhopn),dp1(nrhopn)
	dimension ddr1(nrhopn)
	dimension qdn(nqmax0)
	common /nu1d/ amin,dlna,r1,p1,dr1,dp1,ddr1
	save /nu1d/
	common /fermi/ qdn
	save /fermi/
c
	amin=1.0d-9
	dlna=-log(amin)/(nrhopn-1)
c
c  Check whether correct interpolation table already exists on disk.
	open(12,err=10,file='nu1.dat',status='old',form='unformatted')
	read(12,end=10) amnu1
	err=abs(amnu1/amnu-1.0d0)
	if (err.gt.1.0e-8) then
	  go to 10
	else
	  read(12,end=10) (r1(j),p1(j),j=1,nrhopn)
	  close(12)
c  Yes, interpolation table was okay.
	  go to 20
	end if

c  No, correct interpolation tables aren't on disk: compute and save.
10	close(12)
	do i=1,nrhopn
	  a=amin*exp((i-1)*dlna)
	  call ninu1(a,rhonu,pnu)
	  r1(i)=log(rhonu)
	  p1(i)=log(pnu)
	end do

20	call splini
	call splder(r1,dr1,nrhopn)
	call splder(p1,dp1,nrhopn)
	call splder(dr1,ddr1,nrhopn)

	open(12,file='nu1.dat',status='unknown',form='unformatted')
	rewind 12
	write(12) amnu
	write(12) (r1(j),p1(j),j=1,nrhopn)
	close(12)
c Calculate qdn
	dq=1.0d0
	  do 30 i=1,nqmax0
	  q=i-0.5d0
c	  aq=a*amnu/q
c	  v=1.0d0/sqrt(1.0d0+aq*aq)
	  qdn(i)=dq*q*q*q/(exp(q)+1.0d0)
30	continue

	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine ninu1(a,rhonu,pnu)
c  Compute the density and pressure of one flavor of massive neutrinos,
c  in units of the mean density of one flavor of massless neutrinos.
c  Relative accuracy is better than 1.e-10 for all a.
c
	implicit double precision (a-h,o-z)
c
	common /cosmoparm/ ak,ak2,amnu,lmax,lmaxnr,lmaxnu,
     &                              nqmax,iq0,iq1,iq2
c
	parameter (qmax=30.0d0,nq=1000,nq1=nq+1)
	dimension dum1(nq1),dum2(nq1)
c  const=7*pi**4/120.
	parameter (const=5.68219698d0)
c
	if (amnu.eq.0.0d0) then
	  rhonu=1.0d0
	  pnu=1.0d0/3.0d0
	  return
	end if
c
c  q is the comoving momentum in units of k_B*T_nu0/c.
c  Integrate up to qmax and then use asymptotic expansion for remainder.
	dq=qmax/nq
	dum1(1)=0.0d0
	dum2(1)=0.0d0
	  do 10 i=1,nq
	  q=i*dq
	  aq=a*amnu/q
	  v=1.0d0/sqrt(1.0d0+aq*aq)
	  qdn=dq*q*q*q/(exp(q)+1.0d0)
	  dum1(i+1)=qdn/v
	  dum2(i+1)=qdn*v
10	continue
	call splint(dum1,rhonu,nq1)
	call splint(dum2,pnu,nq1)
c  Apply asymptotic corrrection for q>qmax and normalize by relativistic
c  energy density.
	rhonu=(rhonu+dum1(nq1)/dq)/const
	pnu=(pnu+dum2(nq1)/dq)/const/3.0d0
c
	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine nu2(a,drhonu,fnu,dpnu,shearnu,psi0,psi1,psi2)
c  Compute the perturbations of density, energy flux, pressure, and
c  shear stress of one flavor of massive neutrinos, in units of the mean
c  density of one flavor of massless neutrinos, by integrating over 
c  momentum.
c
	implicit double precision (a-h,o-z)
c
	common /cosmoparm/ ak,ak2,amnu,lmax,lmaxnr,lmaxnu,
     &                        nqmax,iq0,iq1,iq2
	common /fermi/ qdn
	save /fermi/
c
	parameter (nqmax0=15,qmax=nqmax0-0.5d0)
	dimension psi0(nqmax0),psi1(nqmax0),psi2(nqmax0)
c  const=7*pi**4/120.
	parameter (const=5.68219698d0)
c
	dimension g0(4),g1(nqmax0+1),g2(nqmax0+1)
	dimension g3(nqmax0+1),g4(nqmax0+1)
	dimension qdn(nqmax0)
c
	if (nqmax.eq.0) then
	  drhonu=0.0d0
	  fnu=0.0d0
	  dpnu=0.0d0
	  shearnu=0.0d0
	  return
	end if
c
c  q is the comoving momentum in units of k_B*T_nu0/c.
	g1(1)=0.0d0
	g2(1)=0.0d0
	g3(1)=0.0d0
	g4(1)=0.0d0
	do 40 iq=2,(nqmax0+1)
	    q=iq-1.5d0
	    aq=a*amnu/q
	    v=1.0d0/sqrt(1.0d0+aq*aq)
	    g1(iq)=qdn(iq-1)*psi0(iq-1)/v
	    g2(iq)=qdn(iq-1)*psi0(iq-1)*v
	    g3(iq)=qdn(iq-1)*psi1(iq-1)
	    g4(iq)=qdn(iq-1)*psi2(iq-1)*v
40	continue
	call splint(g1,g0(1),nqmax0+1)
	call splint(g2,g0(2),nqmax0+1)
	call splint(g3,g0(3),nqmax0+1)
	call splint(g4,g0(4),nqmax0+1)
 	gf1=g1(nqmax0+1)
	gf2=g2(nqmax0+1)
	gf3=g3(nqmax0+1)
	gf4=g4(nqmax0+1)
	drhonu=(g0(1)+gf1*2.0d0/qmax)/const
	fnu=(g0(3)+gf3*2.0d0/qmax)/const
	dpnu=(g0(2)+gf2*2.0d0/qmax)/const/3.0d0
	shearnu=(g0(4)+gf4*2.0d0/qmax)/const*2.0d0/3.0d0
c
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine nuder(a,adotoa,rhonu,rhonudot,shearnudot,
     2                  psi2,psi2dot)
c  Compute the time derivative of the mean density in massive neutrinos 
c  and the shear perturbation.
c
	implicit double precision (a-h,o-z)
c
	common /cosmoparm/ ak,ak2,amnu,lmax,lmaxnr,
     &                    lmaxnu,nqmax,iq0,iq1,iq2
	common /fermi/ qdn
c
	parameter (nrhopn=10000)
	dimension r1(nrhopn),p1(nrhopn)
	dimension dr1(nrhopn),dp1(nrhopn),ddr1(nrhopn)
	common /nu1d/ amin,dlna,r1,p1,dr1,dp1,ddr1
c
	parameter (nqmax0=15,qmax=nqmax0-0.5d0)
	dimension psi2(nqmax0),psi2dot(nqmax0)
c
	parameter (const=5.68219698d0)
c
	dimension g1(nqmax0+1)
	dimension qdn(nqmax0)
c
	save /nu1d/
	save /fermi/
c
	if (nqmax.eq.0) then
	  rhonudot=0.0d0
	  shearnudot=0.0d0
	  return
	end if
c
c  q is the comoving momentum in units of k_B*T_nu0/c.
	g1(1)=0.0d0
	do 40 iq=2,(nqmax0+1)
	    q=iq-1.5d0
	    aq=a*amnu/q
	    aqdot=aq*adotoa
	    v=1.0d0/sqrt(1.0d0+aq*aq)
	    vdot=-aq*aqdot/(1.0d0+aq*aq)**1.5d0
	    g1(iq)=qdn(iq-1)*(psi2dot(iq-1)*v+psi2(iq-1)*vdot)
40	continue
	call splint(g1,g0,nqmax0+1)
	gf1=g1(nqmax0+1)
	shearnudot=(g0+gf1*2.0d0/qmax)/const*2.0d0/3.0d0
c
c
	d=log(a/amin)/dlna+1.0d0
	i=int(d)
	d=d-i
	if (i.lt.1) then
c  Use linear interpolation
	  rhonudot=dr1(1)+(d-1)*ddr1(1)
	else if (i.gt.nrhopn) then
c  This should not happen, unless the user evolves to z<0!
	  rhonudot=dr1(nrhopn)+(d+i-nrhopn)*ddr1(nrhopn)
	else

c  Cubic spline interpolation for rhonudot.
	  rhonudot=dr1(i)+d*(ddr1(i)+d*(3.0d0*(dr1(i+1)-dr1(i))
     2         -2.0d0*ddr1(i)-ddr1(i+1)+d*(ddr1(i)+ddr1(i+1)
     3         +2.0d0*(dr1(i)-dr1(i+1)))))
	end if
	rhonudot=rhonu*adotoa*rhonudot/dlna
	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	function dtauda(a)
	implicit double precision (a-h,o-z)
c COMMON lingerinc,genparm
	common /lingerinc/ omegab,omegac,omegav,omegan,omegak,h0,
     &                     tcmb,yhe,nnur,nnunr
	common /genparm/ grhom,grhog,grhor,adotrad,taurst,dtaurec
        double precision dtauda
c  Compute expansion rate.
	if (a.lt.1.0d-9) then
	   rhonu=1.0d0
	   pnu=1.0d0/3.0d0
	else
	   call nu1(a,rhonu,pnu)
	end if
c  8*pi*G*rho*a**4.
	grho2=grhom*(omegac+omegab)*a+(grhog+grhor*(nnur+nnunr*rhonu))
     2       +grhom*omegav*a**4+grhom*omegak*a**2
	dtauda=sqrt(3.0d0/grho2)
c
	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function adtauda(a)
        implicit double precision (a-h,o-z)
c COMMON lingerinc,genparm        
c
      common /lingerinc/ omegab,omegac,omegav,omegan,omegak,h0,
     &                     tcmb,yhe,nnur,nnunr
      common /genparm/ grhom,grhog,grhor,adotrad,taurst,dtaurec
              double precision adtauda
c
c  Compute expansion rate.
c  8*pi*G*rho*a**4.
        rhonu=0.0d0
        grho2=grhom*(omegac+omegab)*a+(grhog+grhor*(nnur+nnunr*rhonu))
     2       +grhom*omegav*a**4
        adtauda=a*sqrt(3.0d0/grho2)
c
        return
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine reiopar(optdlss,zri,zristp,rif)

c This subroutine calculates the redshift of reionization
c from an optdlss and sets the reionization fraction rif=1.
      implicit double precision(a-h,o-z)
c COMMON lingerinc
      common /lingerinc/ omegab,omegac,omegav,omegan,omegak,h0
     &                     ,tcmb,yhe,nnur,nnunr

       zri=0.0d0
      rif=0.0d0
      if (optdlss.ne.0) then
         akthom=(2.3048d-9)*(1-yhe)*omegab*h0*h0
         rif=1.0d0

c Calculating the redshift of reionation the parameter
c needed by CMBFAST.
         na=1
         da=0.00001d0
         optd=0.0d0
 5       a=1-na*da
         optd=optd+da*rif*akthom*dtauda(a)/a**2
         if (optd.lt.optdlss) then
            na=na+1
            goto 5
         end if
         zri=1.0d0/a-1.0d0
         zristp=0.07d0*zri-1.0d0
         if (zristp.lt.0.0) zristp=0.0d0
c       TESTING
c        write(*,*)'zristp=',zristp
c       Problems in choosing zristp may lead to inaccuracies,
c       so it is better to always choose zristp=0, this leads
c       to a small increase in computation time. (Bug reported
c       by Andrew Liddle) July 2000.
         zristp=0.0

      end if
      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine splder(y,dy,n)
c  Splder fits a cubic spline to y and returns the first derivatives at
c  the grid points in dy.  Dy is equivalent to a 4th-order Pade
c  difference formula for dy/di.
c
	implicit double precision (a-h,o-z)
	dimension f(20001),g(20001),y(n),dy(n)
	common /spl/ g
	save /spl/ 
c
	n1=n-1
	if (n1.gt.20000) write(*,*) 'Spline array overflow!!! n1=',
     2      n1,'>20000'
c  Quartic fit to dy/di at boundaries, assuming d3y/di3=0.
	f(1)=(-10.0d0*y(1)+15.0d0*y(2)-6.0d0*y(3)+y(4))/6.0d0
	f(n)=(10.0d0*y(n)-15.0d0*y(n1)+6.0d0*y(n-2)-y(n-3))/6.0d0
c  Solve the tridiagonal system
c  dy(i-1)+4*dy(i)+dy(i+1)=3*(y(i+1)-y(i-1)), i=2,3,...,n1,
c  with dy(1)=f(1), dy(n)=f(n).
	  do 10 i=2,n1
	  f(i)=g(i)*(3.0d0*(y(i+1)-y(i-1))-f(i-1))
10	continue
	dy(n)=f(n)
	  do 20 i=n1,1,-1
	  dy(i)=f(i)-g(i)*dy(i+1)
20	continue
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine splini
c  Splini must be called before splder to initialize array g in common.
c
	implicit double precision (a-h,o-z)
	dimension g(20001)
	common /spl/ g
	save /spl/
c
	g(1)=0.0d0
	  do 10 i=2,20001
	  g(i)=1.0d0/(4.0d0-g(i-1))
10	continue
	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	function rombint(f,a,b,tol)
c  Rombint returns the integral from a to b of using Romberg integration.
c  The method converges provided that f(x) is continuous in (a,b).  
c  f must be double precision and must be declared external in the calling
c  routine.  tol indicates the desired relative accuracy in the integral.

	parameter (MAXITER=18,MAXJ=5)
	implicit double precision (a-h,o-z)
	dimension g(MAXJ+1)
	double precision f
	external f
c
	h=0.5d0*(b-a)
	gmax=h*(f(a)+f(b))
	g(1)=gmax
	nint=1
	error=1.0d20
	i=0
10	  i=i+1
	  if (i.gt.MAXITER.or.(i.gt.5.and.abs(error).lt.tol))
     2      go to 40
c  Calculate next trapezoidal rule approximation to integral.
	  g0=0.0d0
	    do 20 k=1,nint
	    g0=g0+f(a+(k+k-1)*h)
20	  continue
	  g0=0.5d0*g(1)+h*g0
	  h=0.5d0*h
	  nint=nint+nint
	  jmax=min(i,MAXJ)
	  fourj=1.0d0
	    do 30 j=1,jmax
c  Use Richardson extrapolation.
	    fourj=4.0d0*fourj
	    g1=g0+(g0-g(j))/(fourj-1.0d0)
	    g(j)=g0
	    g0=g1
30	  continue
	  if (abs(g0).gt.tol) then
	    error=1.0d0-gmax/g0
	  else
	    error=gmax
	  end if
	  gmax=g0
	  g(jmax+1)=g0
	go to 10
40	rombint=g0
	if (i.gt.MAXITER.and.abs(error).gt.tol)
     2    write(*,*) 'Rombint failed to converge; integral, error=',
     3    rombint,error
	return
	end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE splint1(xa,ya,y2a,n,x,y)
      implicit double precision(a-h,o-z)
      INTEGER n
      double precision x,y,xa(n),y2a(n),ya(n)
      INTEGER k,khi,klo
      double precision a,b,h
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.d0) then
      write(*,*)'bad xa input in splint'
      stop
      endif
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**
     *2)/6.d0
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 0@.1Y..

      SUBROUTINE splint2(xa,ya,y2a,n,x,ypr)
      implicit double precision(a-h,o-z)
      INTEGER n
      double precision xa(n),y2a(n),ya(n),x,dely
      INTEGER k,khi,klo
      double precision a,b,h
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) then
      write(*,*)'bad xa input in splint'
      stop
      endif
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      dely=ya(khi)-ya(klo)
      ypr=dely/h-(3.d0*a**2-1.0d0)*h*y2a(klo)/6.0d0 + 
     & (3.d0*b**2-1.0d0)*h*y2a(khi)/6.0d0 

      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 0@.1Y..


      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      INTEGER n,NMAX
      real*8 yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=50010)
      INTEGER i,k
      real*8 p,qn,sig,un,u(NMAX)
      if (yp1.gt..99d30) then
        y2(1)=0.d0
        u(1)=0.d0
      else
        y2(1)=-0.5d0
        u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.d0
        y2(i)=(sig-1.d0)/p
        u(i)=(6.d0*((y(i+1)-y(i))/(x(i+
     *1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*
     *u(i-1))/p
11    continue
      if (ypn.gt..99d30) then
        qn=0.d0
        un=0.d0
      else
        qn=0.5d0
        un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
C  (C) Copr. 1986-92 Numerical Recipes Software =$j*m,).
      END

c

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
	subroutine splint(y,z,n)
c  Splint integrates a cubic spline, providing the ouput value
c  z = integral from 1 to n of s(i)di, where s(i) is the spline fit
c  to y(i).

	implicit double precision (a-h,o-z)
	dimension y(n)

	n1=n-1
c  Cubic fit to dy/di at boundaries.
	dy1=(-11.0d0*y(1)+18.0d0*y(2)-9.0d0*y(3)+2.0d0*y(4))/6.0d0
	dy1=0.0d0
	dyn=(11.0d0*y(n)-18.0d0*y(n1)+9.0d0*y(n-2)-2.0d0*y(n-3))/6.0d0

	z=0.5d0*(y(1)+y(n))+(dy1-dyn)/12.0d0
	  do 10 i=2,n1
	  z=z+y(i)
10	continue
	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine outtransf(n,y,curv,ik,itf)

	implicit double precision(a-h,o-z)

        dimension y(n)
	parameter (nnmax=1)
c COMMON lingerinc
	common /cosmoparm/ beta,beta2,amnu,lmaxg,lmaxnr,
     &                              lmaxnu,nqmax,iq0,iq1,iq2
        common /lingerinc/ omegab,omegac,omegav,omegan,omegak,h0,
     &                     tcmb,yhe,nnur,nnunr


	tfc=y(4)/beta2
	tfb=y(6)/beta2
	tfg=y(8)/beta2
	tfn=y(10+2*lmaxg)/beta2

	if (amnu.ne.0.0d0) then
           a=y(1)
           call nu1(a,rhonu,pnu)
           call nu2(a,drhonu,fnu,dpnu,shearnu,y(iq0),y(iq1),y(iq2))
           deltan=drhonu/rhonu
           tfnm=deltan/beta2
	endif
c	if (amnu.eq.0.0d0) then 
c output k in units of Mpc/h
c           write(12+itf,'(5E14.6)')beta/(h0/100.d0),tfc,tfb
c	else
c	   write(12+itf,'(6E14.6)')beta/(h0/100.d0),tfc,tfb
c	endif
	
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine out(a,b,c,n,file)

c This subroutine just writes vectors a and b 
c of lenght n to file      

        real*8 a(n),b(n),c(n)
        integer n,file,i

        do i=2,n
           write(file,*) i,a(i),b(i),c(i)
        end do

        return
        end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine dverk (n, fcn, x, y, xend, tol, ind, c, nw, w)
      integer n, ind, nw, k
      double precision x, y(n), xend, tol, c(*), w(nw,9), temp
c
c***********************************************************************
c                                                                      *
c note added 11/14/85.                                                 *
c                                                                      *
c if you discover any errors in this subroutine, please contact        *
c                                                                      *
c        kenneth r. jackson                                            *
c        department of computer science                                *
c        university of toronto                                         *
c        toronto, ontario,                                             *
c        canada   m5s 1a4                                              *
c                                                                      *
c        phone: 416-978-7075                                           *
c                                                                      *
c        electronic mail:                                              *
c        uucp:   {cornell,decvax,ihnp4,linus,uw-beaver}!utcsri!krj     *
c        csnet:  krj@toronto                                           *
c        arpa:   krj.toronto@csnet-relay                               *
c        bitnet: krj%toronto@csnet-relay.arpa                          *
c                                                                      *
c dverk is written in fortran 66.                                      *
c                                                                      *
c the constants dwarf and rreb -- c(10) and c(11), respectively -- are *
c set for a  vax  in  double  precision.  they  should  be  reset,  as *
c described below, if this program is run on another machine.          *
c                                                                      *
c the c array is declared in this subroutine to have one element only, *
c although  more  elements  are  referenced  in this subroutine.  this *
c causes some compilers to issue warning messages.  there is,  though, *
c no  error  provided  c is declared sufficiently large in the calling *
c program, as described below.                                         *
c                                                                      *
c the following external statement  for  fcn  was  added  to  avoid  a *
c warning  message  from  the  unix  f77 compiler.  the original dverk *
c comments and code follow it.                                         *
c                                                                      *
c***********************************************************************
c
      external fcn
c
c***********************************************************************
c                                                                      *
c     purpose - this is a runge-kutta  subroutine  based  on  verner's *
c fifth and sixth order pair of formulas for finding approximations to *
c the solution of  a  system  of  first  order  ordinary  differential *
c equations  with  initial  conditions. it attempts to keep the global *
c error proportional to  a  tolerance  specified  by  the  user.  (the *
c proportionality  depends  on the kind of error control that is used, *
c as well as the differential equation and the range of integration.)  *
c                                                                      *
c     various options are available to the user,  including  different *
c kinds  of  error control, restrictions on step sizes, and interrupts *
c which permit the user to examine the state of the  calculation  (and *
c perhaps make modifications) during intermediate stages.              *
c                                                                      *
c     the program is efficient for non-stiff systems.  however, a good *
c variable-order-adams  method  will probably be more efficient if the *
c function evaluations are very costly.  such a method would  also  be *
c more suitable if one wanted to obtain a large number of intermediate *
c solution values by interpolation, as might be the case  for  example *
c with graphical output.                                               *
c                                                                      *
c                                    hull-enright-jackson   1/10/76    *
c                                                                      *
c***********************************************************************
c                                                                      *
c     use - the user must specify each of the following                *
c                                                                      *
c     n  number of equations                                           *
c                                                                      *
c   fcn  name of subroutine for evaluating functions - the  subroutine *
c           itself must also be provided by the user - it should be of *
c           the following form                                         *
c              subroutine fcn(n, x, y, yprime)                         *
c              integer n                                               *
c              double precision x, y(n), yprime(n)                     *
c                      *** etc ***                                     *
c           and it should evaluate yprime, given n, x and y            *
c                                                                      *
c     x  independent variable - initial value supplied by user         *
c                                                                      *
c     y  dependent variable - initial values of components y(1), y(2), *
c           ..., y(n) supplied by user                                 *
c                                                                      *
c  xend  value of x to which integration is to be carried out - it may *
c           be less than the initial value of x                        *
c                                                                      *
c   tol  tolerance - the subroutine attempts to control a norm of  the *
c           local  error  in  such  a  way  that  the  global error is *
c           proportional to tol. in some problems there will be enough *
c           damping  of  errors, as well as some cancellation, so that *
c           the global error will be less than tol. alternatively, the *
c           control   can   be  viewed  as  attempting  to  provide  a *
c           calculated value of y at xend which is the exact  solution *
c           to  the  problem y' = f(x,y) + e(x) where the norm of e(x) *
c           is proportional to tol.  (the norm  is  a  max  norm  with *
c           weights  that  depend on the error control strategy chosen *
c           by the user.  the default weight for the k-th component is *
c           1/max(1,abs(y(k))),  which therefore provides a mixture of *
c           absolute and relative error control.)                      *
c                                                                      *
c   ind  indicator - on initial entry ind must be set equal to  either *
c           1  or  2. if the user does not wish to use any options, he *
c           should set ind to 1 - all that remains for the user to  do *
c           then  is  to  declare c and w, and to specify nw. the user *
c           may also  select  various  options  on  initial  entry  by *
c           setting ind = 2 and initializing the first 9 components of *
c           c as described in the next section.  he may also  re-enter *
c           the  subroutine  with ind = 3 as mentioned again below. in *
c           any event, the subroutine returns with ind equal to        *
c              3 after a normal return                                 *
c              4, 5, or 6 after an interrupt (see options c(8), c(9))  *
c              -1, -2, or -3 after an error condition (see below)      *
c                                                                      *
c     c  communications vector - the dimension must be greater than or *
c           equal to 24, unless option c(1) = 4 or 5 is used, in which *
c           case the dimension must be greater than or equal to n+30   *
c                                                                      *
c    nw  first dimension of workspace w -  must  be  greater  than  or *
c           equal to n                                                 *
c                                                                      *
c     w  workspace matrix - first dimension must be nw and second must *
c           be greater than or equal to 9                              *
c                                                                      *
c     the subroutine  will  normally  return  with  ind  =  3,  having *
c replaced the initial values of x and y with, respectively, the value *
c of xend and an approximation to y at xend.  the  subroutine  can  be *
c called  repeatedly  with new values of xend without having to change *
c any other argument.  however, changes in tol, or any of the  options *
c described below, may also be made on such a re-entry if desired.     *
c                                                                      *
c     three error returns are also possible, in which  case  x  and  y *
c will be the most recently accepted values -                          *
c     with ind = -3 the subroutine was unable  to  satisfy  the  error *
c        requirement  with a particular step-size that is less than or *
c        equal to hmin, which may mean that tol is too small           *
c     with ind = -2 the value of hmin  is  greater  than  hmax,  which *
c        probably  means  that the requested tol (which is used in the *
c        calculation of hmin) is too small                             *
c     with ind = -1 the allowed maximum number of fcn evaluations  has *
c        been  exceeded,  but  this  can only occur if option c(7), as *
c        described in the next section, has been used                  *
c                                                                      *
c     there are several circumstances that will cause the calculations *
c to  be  terminated,  along with output of information that will help *
c the user determine the cause of  the  trouble.  these  circumstances *
c involve  entry with illegal or inconsistent values of the arguments, *
c such as attempting a normal  re-entry  without  first  changing  the *
c value of xend, or attempting to re-enter with ind less than zero.    *
c                                                                      *
c***********************************************************************
c                                                                      *
c     options - if the subroutine is entered with ind = 1, the first 9 *
c components of the communications vector are initialized to zero, and *
c the subroutine uses only default values  for  each  option.  if  the *
c subroutine  is  entered  with ind = 2, the user must specify each of *
c these 9 components - normally he would first set them all  to  zero, *
c and  then  make  non-zero  those  that  correspond to the particular *
c options he wishes to select. in any event, options may be changed on *
c re-entry  to  the  subroutine  -  but if the user changes any of the *
c options, or tol, in the course of a calculation he should be careful *
c about  how  such changes affect the subroutine - it may be better to *
c restart with ind = 1 or 2. (components 10 to 24 of c are used by the *
c program  -  the information is available to the user, but should not *
c normally be changed by him.)                                         *
c                                                                      *
c  c(1)  error control indicator - the norm of the local error is  the *
c           max  norm  of  the  weighted  error  estimate  vector, the *
c           weights being determined according to the value of c(1) -  *
c              if c(1)=1 the weights are 1 (absolute error control)    *
c              if c(1)=2 the weights are 1/abs(y(k))  (relative  error *
c                 control)                                             *
c              if c(1)=3 the  weights  are  1/max(abs(c(2)),abs(y(k))) *
c                 (relative  error  control,  unless abs(y(k)) is less *
c                 than the floor value, abs(c(2)) )                    *
c              if c(1)=4 the weights are 1/max(abs(c(k+30)),abs(y(k))) *
c                 (here individual floor values are used)              *
c              if c(1)=5 the weights are 1/abs(c(k+30))                *
c              for all other values of c(1), including  c(1) = 0,  the *
c                 default  values  of  the  weights  are  taken  to be *
c                 1/max(1,abs(y(k))), as mentioned earlier             *
c           (in the two cases c(1) = 4 or 5 the user must declare  the *
c           dimension of c to be at least n+30 and must initialize the *
c           components c(31), c(32), ..., c(n+30).)                    *
c                                                                      *
c  c(2)  floor value - used when the indicator c(1) has the value 3    *
c                                                                      *
c  c(3)  hmin specification - if not zero, the subroutine chooses hmin *
c           to be abs(c(3)) - otherwise it uses the default value      *
c              10*max(dwarf,rreb*max(weighted norm y/tol,abs(x))),     *
c           where dwarf is a very small positive  machine  number  and *
c           rreb is the relative roundoff error bound                  *
c                                                                      *
c  c(4)  hstart specification - if not zero, the subroutine  will  use *
c           an  initial  hmag equal to abs(c(4)), except of course for *
c           the restrictions imposed by hmin and hmax  -  otherwise it *
c           uses the default value of hmax*(tol)**(1/6)                *
c                                                                      *
c  c(5)  scale specification - this is intended to be a measure of the *
c           scale of the problem - larger values of scale tend to make *
c           the method more reliable, first  by  possibly  restricting *
c           hmax  (as  described  below) and second, by tightening the *
c           acceptance requirement - if c(5) is zero, a default  value *
c           of  1  is  used.  for  linear  homogeneous  problems  with *
c           constant coefficients, an appropriate value for scale is a *
c           norm  of  the  associated  matrix.  for other problems, an *
c           approximation to  an  average  value  of  a  norm  of  the *
c           jacobian along the trajectory may be appropriate           *
c                                                                      *
c  c(6)  hmax specification - four cases are possible                  *
c           if c(6).ne.0 and c(5).ne.0, hmax is taken to be            *
c              min(abs(c(6)),2/abs(c(5)))                              *
c           if c(6).ne.0 and c(5).eq.0, hmax is taken to be  abs(c(6)) *
c           if c(6).eq.0 and c(5).ne.0, hmax is taken to be            *
c              2/abs(c(5))                                             *
c           if c(6).eq.0 and c(5).eq.0, hmax is given a default  value *
c              of 2                                                    *
c                                                                      *
c  c(7)  maximum number of function evaluations  -  if  not  zero,  an *
c           error  return with ind = -1 will be caused when the number *
c           of function evaluations exceeds abs(c(7))                  *
c                                                                      *
c  c(8)  interrupt number  1  -  if  not  zero,  the  subroutine  will *
c           interrupt   the  calculations  after  it  has  chosen  its *
c           preliminary value of hmag, and just before choosing htrial *
c           and  xtrial  in  preparation for taking a step (htrial may *
c           differ from hmag in sign, and may  require  adjustment  if *
c           xend  is  near) - the subroutine returns with ind = 4, and *
c           will resume calculation at the point  of  interruption  if *
c           re-entered with ind = 4                                    *
c                                                                      *
c  c(9)  interrupt number  2  -  if  not  zero,  the  subroutine  will *
c           interrupt   the  calculations  immediately  after  it  has *
c           decided whether or not to accept the result  of  the  most *
c           recent  trial step, with ind = 5 if it plans to accept, or *
c           ind = 6 if it plans to reject -  y(*)  is  the  previously *
c           accepted  result, while w(*,9) is the newly computed trial *
c           value, and w(*,2) is the unweighted error estimate vector. *
c           the  subroutine  will  resume calculations at the point of *
c           interruption on re-entry with ind = 5 or 6. (the user  may *
c           change ind in this case if he wishes, for example to force *
c           acceptance of a step that would otherwise be rejected,  or *
c           vice versa. he can also restart with ind = 1 or 2.)        *
c                                                                      *
c***********************************************************************
c                                                                      *
c  summary of the components of the communications vector              *
c                                                                      *
c     prescribed at the option       determined by the program         *
c           of the user                                                *
c                                                                      *
c                                    c(10) rreb(rel roundoff err bnd)  *
c     c(1) error control indicator   c(11) dwarf (very small mach no)  *
c     c(2) floor value               c(12) weighted norm y             *
c     c(3) hmin specification        c(13) hmin                        *
c     c(4) hstart specification      c(14) hmag                        *
c     c(5) scale specification       c(15) scale                       *
c     c(6) hmax specification        c(16) hmax                        *
c     c(7) max no of fcn evals       c(17) xtrial                      *
c     c(8) interrupt no 1            c(18) htrial                      *
c     c(9) interrupt no 2            c(19) est                         *
c                                    c(20) previous xend               *
c                                    c(21) flag for xend               *
c                                    c(22) no of successful steps      *
c                                    c(23) no of successive failures   *
c                                    c(24) no of fcn evals             *
c                                                                      *
c  if c(1) = 4 or 5, c(31), c(32), ... c(n+30) are floor values        *
c                                                                      *
c***********************************************************************
c                                                                      *
c  an overview of the program                                          *
c                                                                      *
c     begin initialization, parameter checking, interrupt re-entries   *
c  ......abort if ind out of range 1 to 6                              *
c  .     cases - initial entry, normal re-entry, interrupt re-entries  *
c  .     case 1 - initial entry (ind .eq. 1 or 2)                      *
c  v........abort if n.gt.nw or tol.le.0                               *
c  .        if initial entry without options (ind .eq. 1)              *
c  .           set c(1) to c(9) equal to zero                          *
c  .        else initial entry with options (ind .eq. 2)               *
c  .           make c(1) to c(9) non-negative                          *
c  .           make floor values non-negative if they are to be used   *
c  .        end if                                                     *
c  .        initialize rreb, dwarf, prev xend, flag, counts            *
c  .     case 2 - normal re-entry (ind .eq. 3)                         *
c  .........abort if xend reached, and either x changed or xend not    *
c  .        re-initialize flag                                         *
c  .     case 3 - re-entry following an interrupt (ind .eq. 4 to 6)    *
c  v        transfer control to the appropriate re-entry point.......  *
c  .     end cases                                                  .  *
c  .  end initialization, etc.                                      .  *
c  .                                                                v  *
c  .  loop through the following 4 stages, once for each trial step .  *
c  .     stage 1 - prepare                                          .  *
c***********error return (with ind=-1) if no of fcn evals too great .  *
c  .        calc slope (adding 1 to no of fcn evals) if ind .ne. 6  .  *
c  .        calc hmin, scale, hmax                                  .  *
c***********error return (with ind=-2) if hmin .gt. hmax            .  *
c  .        calc preliminary hmag                                   .  *
c***********interrupt no 1 (with ind=4) if requested.......re-entry.v  *
c  .        calc hmag, xtrial and htrial                            .  *
c  .     end stage 1                                                .  *
c  v     stage 2 - calc ytrial (adding 7 to no of fcn evals)        .  *
c  .     stage 3 - calc the error estimate                          .  *
c  .     stage 4 - make decisions                                   .  *
c  .        set ind=5 if step acceptable, else set ind=6            .  *
c***********interrupt no 2 if requested....................re-entry.v  *
c  .        if step accepted (ind .eq. 5)                              *
c  .           update x, y from xtrial, ytrial                         *
c  .           add 1 to no of successful steps                         *
c  .           set no of successive failures to zero                   *
c**************return(with ind=3, xend saved, flag set) if x .eq. xend *
c  .        else step not accepted (ind .eq. 6)                        *
c  .           add 1 to no of successive failures                      *
c**************error return (with ind=-3) if hmag .le. hmin            *
c  .        end if                                                     *
c  .     end stage 4                                                   *
c  .  end loop                                                         *
c  .                                                                   *
c  begin abort action                                                  *
c     output appropriate  message  about  stopping  the  calculations, *
c        along with values of ind, n, nw, tol, hmin,  hmax,  x,  xend, *
c        previous xend,  no of  successful  steps,  no  of  successive *
c        failures, no of fcn evals, and the components of y            *
c     stop                                                             *
c  end abort action                                                    *
c                                                                      *
c***********************************************************************
c
c     ******************************************************************
c     * begin initialization, parameter checking, interrupt re-entries *
c     ******************************************************************
c
c  ......abort if ind out of range 1 to 6
         if (ind.lt.1 .or. ind.gt.6) go to 500
c
c        cases - initial entry, normal re-entry, interrupt re-entries
         go to (5, 5, 45, 1111, 2222, 2222), ind
c        case 1 - initial entry (ind .eq. 1 or 2)
c  .........abort if n.gt.nw or tol.le.0
    5       if (n.gt.nw .or. tol.le.0.d0) go to 500
            if (ind.eq. 2) go to 15
c              initial entry without options (ind .eq. 1)
c              set c(1) to c(9) equal to 0
               do 10 k = 1, 9
                  c(k) = 0.d0
   10          continue
               go to 35
   15       continue
c              initial entry with options (ind .eq. 2)
c              make c(1) to c(9) non-negative
               do 20 k = 1, 9
                  c(k) = dabs(c(k))
   20          continue
c              make floor values non-negative if they are to be used
               if (c(1).ne.4.d0 .and. c(1).ne.5.d0) go to 30
                  do 25 k = 1, n
                     c(k+30) = dabs(c(k+30))
   25             continue
   30          continue
   35       continue
c           initialize rreb, dwarf, prev xend, flag, counts
            c(10) = 2.d0**(-56)
            c(11) = 1.d-35
c           set previous xend initially to initial value of x
            c(20) = x
            do 40 k = 21, 24
               c(k) = 0.d0
   40       continue
            go to 50
c        case 2 - normal re-entry (ind .eq. 3)
c  .........abort if xend reached, and either x changed or xend not
   45       if (c(21).ne.0.d0 .and.
     +                        (x.ne.c(20) .or. xend.eq.c(20))) go to 500
c           re-initialize flag
            c(21) = 0.d0
            go to 50
c        case 3 - re-entry following an interrupt (ind .eq. 4 to 6)
c           transfer control to the appropriate re-entry point..........
c           this has already been handled by the computed go to        .
c        end cases                                                     v
   50    continue
c
c     end initialization, etc.
c
c     ******************************************************************
c     * loop through the following 4 stages, once for each trial  step *
c     * until the occurrence of one of the following                   *
c     *    (a) the normal return (with ind .eq. 3) on reaching xend in *
c     *        stage 4                                                 *
c     *    (b) an error return (with ind .lt. 0) in stage 1 or stage 4 *
c     *    (c) an interrupt return (with ind  .eq.  4,  5  or  6),  if *
c     *        requested, in stage 1 or stage 4                        *
c     ******************************************************************
c
99999 continue
c
c        ***************************************************************
c        * stage 1 - prepare - do calculations of  hmin,  hmax,  etc., *
c        * and some parameter  checking,  and  end  up  with  suitable *
c        * values of hmag, xtrial and htrial in preparation for taking *
c        * an integration step.                                        *
c        ***************************************************************
c
c***********error return (with ind=-1) if no of fcn evals too great
            if (c(7).eq.0.d0 .or. c(24).lt.c(7)) go to 100
               ind = -1
               return
  100       continue
c
c           calculate slope (adding 1 to no of fcn evals) if ind .ne. 6
            if (ind .eq. 6) go to 105
               call fcn(n, x, y, w(1,1))
               c(24) = c(24) + 1.d0
  105       continue
c
c           calculate hmin - use default unless value prescribed
            c(13) = c(3)
            if (c(3) .ne. 0.d0) go to 165
c              calculate default value of hmin
c              first calculate weighted norm y - c(12) - as specified
c              by the error control indicator c(1)
               temp = 0.d0
               if (c(1) .ne. 1.d0) go to 115
c                 absolute error control - weights are 1
                  do 110 k = 1, n
                     temp = dmax1(temp, dabs(y(k)))
  110             continue
                  c(12) = temp
                  go to 160
  115          if (c(1) .ne. 2.d0) go to 120
c                 relative error control - weights are 1/dabs(y(k)) so
c                 weighted norm y is 1
                  c(12) = 1.d0
                  go to 160
  120          if (c(1) .ne. 3.d0) go to 130
c                 weights are 1/max(c(2),abs(y(k)))
                  do 125 k = 1, n
                     temp = dmax1(temp, dabs(y(k))/c(2))
  125             continue
                  c(12) = dmin1(temp, 1.d0)
                  go to 160
  130          if (c(1) .ne. 4.d0) go to 140
c                 weights are 1/max(c(k+30),abs(y(k)))
                  do 135 k = 1, n
                     temp = dmax1(temp, dabs(y(k))/c(k+30))
  135             continue
                  c(12) = dmin1(temp, 1.d0)
                  go to 160
  140          if (c(1) .ne. 5.d0) go to 150
c                 weights are 1/c(k+30)
                  do 145 k = 1, n
                     temp = dmax1(temp, dabs(y(k))/c(k+30))
  145             continue
                  c(12) = temp
                  go to 160
  150          continue
c                 default case - weights are 1/max(1,abs(y(k)))
                  do 155 k = 1, n
                     temp = dmax1(temp, dabs(y(k)))
  155             continue
                  c(12) = dmin1(temp, 1.d0)
  160          continue
               c(13) = 10.d0*dmax1(c(11),c(10)*dmax1(c(12)/tol,dabs(x)))
  165       continue
c
c           calculate scale - use default unless value prescribed
            c(15) = c(5)
            if (c(5) .eq. 0.d0) c(15) = 1.d0
c
c           calculate hmax - consider 4 cases
c           case 1 both hmax and scale prescribed
               if (c(6).ne.0.d0 .and. c(5).ne.0.d0)
     +                                    c(16) = dmin1(c(6), 2.d0/c(5))
c           case 2 - hmax prescribed, but scale not
               if (c(6).ne.0.d0 .and. c(5).eq.0.d0) c(16) = c(6)
c           case 3 - hmax not prescribed, but scale is
               if (c(6).eq.0.d0 .and. c(5).ne.0.d0) c(16) = 2.d0/c(5)
c           case 4 - neither hmax nor scale is provided
               if (c(6).eq.0.d0 .and. c(5).eq.0.d0) c(16) = 2.d0
c
c***********error return (with ind=-2) if hmin .gt. hmax
            if (c(13) .le. c(16)) go to 170
               ind = -2
               return
  170       continue
c
c           calculate preliminary hmag - consider 3 cases
            if (ind .gt. 2) go to 175
c           case 1 - initial entry - use prescribed value of hstart, if
c              any, else default
               c(14) = c(4)
               if (c(4) .eq. 0.d0) c(14) = c(16)*tol**(1.0d0/6.0d0)
               go to 185
  175       if (c(23) .gt. 1.d0) go to 180
c           case 2 - after a successful step, or at most  one  failure,
c              use min(2, .9*(tol/est)**(1/6))*hmag, but avoid possible
c              overflow. then avoid reduction by more than half.
               temp = 2.d0*c(14)
               if (tol .lt. (2.d0/.9d0)**6*c(19))
     +                 temp = .9d0*(tol/c(19))**(1.0d0/6.0d0)*c(14)
               c(14) = dmax1(temp, .5d0*c(14))
               go to 185
  180       continue
c           case 3 - after two or more successive failures
               c(14) = .5d0*c(14)
  185       continue
c
c           check against hmax
            c(14) = dmin1(c(14), c(16))
c
c           check against hmin
            c(14) = dmax1(c(14), c(13))
c
c***********interrupt no 1 (with ind=4) if requested
            if (c(8) .eq. 0.d0) go to 1111
               ind = 4
               return
c           resume here on re-entry with ind .eq. 4   ........re-entry..
 1111       continue
c
c           calculate hmag, xtrial - depending on preliminary hmag, xend
            if (c(14) .ge. dabs(xend - x)) go to 190
c              do not step more than half way to xend
               c(14) = dmin1(c(14), .5d0*dabs(xend - x))
               c(17) = x + dsign(c(14), xend - x)
               go to 195
  190       continue
c              hit xend exactly
               c(14) = dabs(xend - x)
               c(17) = xend
  195       continue
c
c           calculate htrial
            c(18) = c(17) - x
c
c        end stage 1
c
c        ***************************************************************
c        * stage 2 - calculate ytrial (adding 7 to no of  fcn  evals). *
c        * w(*,2), ... w(*,8)  hold  intermediate  results  needed  in *
c        * stage 3. w(*,9) is temporary storage until finally it holds *
c        * ytrial.                                                     *
c        ***************************************************************
c
            temp = c(18)/1398169080000.d0
c
            do 200 k = 1, n
               w(k,9) = y(k) + temp*w(k,1)*233028180000.d0
  200       continue
            call fcn(n, x + c(18)/6.d0, w(1,9), w(1,2))
c
            do 205 k = 1, n
               w(k,9) = y(k) + temp*(   w(k,1)*74569017600.d0
     +                                + w(k,2)*298276070400.d0  )
  205       continue
            call fcn(n, x + c(18)*(4.d0/15.d0), w(1,9), w(1,3))
c
            do 210 k = 1, n
               w(k,9) = y(k) + temp*(   w(k,1)*1165140900000.d0
     +                                - w(k,2)*3728450880000.d0
     +                                + w(k,3)*3495422700000.d0 )
  210       continue
            call fcn(n, x + c(18)*(2.d0/3.d0), w(1,9), w(1,4))
c
            do 215 k = 1, n
               w(k,9) = y(k) + temp*( - w(k,1)*3604654659375.d0
     +                                + w(k,2)*12816549900000.d0
     +                                - w(k,3)*9284716546875.d0
     +                                + w(k,4)*1237962206250.d0 )
  215       continue
            call fcn(n, x + c(18)*(5.d0/6.d0), w(1,9), w(1,5))
c
            do 220 k = 1, n
               w(k,9) = y(k) + temp*(   w(k,1)*3355605792000.d0
     +                                - w(k,2)*11185352640000.d0
     +                                + w(k,3)*9172628850000.d0
     +                                - w(k,4)*427218330000.d0
     +                                + w(k,5)*482505408000.d0  )
  220       continue
            call fcn(n, x + c(18), w(1,9), w(1,6))
c
            do 225 k = 1, n
               w(k,9) = y(k) + temp*( - w(k,1)*770204740536.d0
     +                                + w(k,2)*2311639545600.d0
     +                                - w(k,3)*1322092233000.d0
     +                                - w(k,4)*453006781920.d0
     +                                + w(k,5)*326875481856.d0  )
  225       continue
            call fcn(n, x + c(18)/15.d0, w(1,9), w(1,7))
c
            do 230 k = 1, n
               w(k,9) = y(k) + temp*(   w(k,1)*2845924389000.d0
     +                                - w(k,2)*9754668000000.d0
     +                                + w(k,3)*7897110375000.d0
     +                                - w(k,4)*192082660000.d0
     +                                + w(k,5)*400298976000.d0
     +                                + w(k,7)*201586000000.d0  )
  230       continue
            call fcn(n, x + c(18), w(1,9), w(1,8))
c
c           calculate ytrial, the extrapolated approximation and store
c              in w(*,9)
            do 235 k = 1, n
               w(k,9) = y(k) + temp*(   w(k,1)*104862681000.d0
     +                                + w(k,3)*545186250000.d0
     +                                + w(k,4)*446637345000.d0
     +                                + w(k,5)*188806464000.d0
     +                                + w(k,7)*15076875000.d0
     +                                + w(k,8)*97599465000.d0   )
  235       continue
c
c           add 7 to the no of fcn evals
            c(24) = c(24) + 7.d0
c
c        end stage 2
c
c        ***************************************************************
c        * stage 3 - calculate the error estimate est. first calculate *
c        * the  unweighted  absolute  error  estimate vector (per unit *
c        * step) for the unextrapolated approximation and store it  in *
c        * w(*,2).  then  calculate the weighted max norm of w(*,2) as *
c        * specified by the error  control  indicator  c(1).  finally, *
c        * modify  this result to produce est, the error estimate (per *
c        * unit step) for the extrapolated approximation ytrial.       *
c        ***************************************************************
c
c           calculate the unweighted absolute error estimate vector
            do 300 k = 1, n
               w(k,2) = (   w(k,1)*8738556750.d0
     +                    + w(k,3)*9735468750.d0
     +                    - w(k,4)*9709507500.d0
     +                    + w(k,5)*8582112000.d0
     +                    + w(k,6)*95329710000.d0
     +                    - w(k,7)*15076875000.d0
     +                    - w(k,8)*97599465000.d0)/1398169080000.d0
  300       continue
c
c           calculate the weighted max norm of w(*,2) as specified by
c           the error control indicator c(1)
            temp = 0.d0
            if (c(1) .ne. 1.d0) go to 310
c              absolute error control
               do 305 k = 1, n
                  temp = dmax1(temp,dabs(w(k,2)))
  305          continue
               go to 360
  310       if (c(1) .ne. 2.d0) go to 320
c              relative error control
               do 315 k = 1, n
                  temp = dmax1(temp, dabs(w(k,2)/y(k)))
  315          continue
               go to 360
  320       if (c(1) .ne. 3.d0) go to 330
c              weights are 1/max(c(2),abs(y(k)))
               do 325 k = 1, n
                  temp = dmax1(temp, dabs(w(k,2))
     +                             / dmax1(c(2), dabs(y(k))) )
  325          continue
               go to 360
  330       if (c(1) .ne. 4.d0) go to 340
c              weights are 1/max(c(k+30),abs(y(k)))
               do 335 k = 1, n
                  temp = dmax1(temp, dabs(w(k,2))
     +                             / dmax1(c(k+30), dabs(y(k))) )
  335          continue
               go to 360
  340       if (c(1) .ne. 5.d0) go to 350
c              weights are 1/c(k+30)
               do 345 k = 1, n
                  temp = dmax1(temp, dabs(w(k,2)/c(k+30)))
  345          continue
               go to 360
  350       continue
c              default case - weights are 1/max(1,abs(y(k)))
               do 355 k = 1, n
                  temp = dmax1(temp, dabs(w(k,2))
     +                             / dmax1(1.d0, dabs(y(k))) )
  355          continue
  360       continue
c
c           calculate est - (the weighted max norm of w(*,2))*hmag*scale
c              - est is intended to be a measure of the error  per  unit
c              step in ytrial
            c(19) = temp*c(14)*c(15)
c
c        end stage 3
c
c        ***************************************************************
c        * stage 4 - make decisions.                                   *
c        ***************************************************************
c
c           set ind=5 if step acceptable, else set ind=6
            ind = 5
            if (c(19) .gt. tol) ind = 6
c
c***********interrupt no 2 if requested
            if (c(9) .eq. 0.d0) go to 2222
               return
c           resume here on re-entry with ind .eq. 5 or 6   ...re-entry..
 2222       continue
c
            if (ind .eq. 6) go to 410
c              step accepted (ind .eq. 5), so update x, y from xtrial,
c                 ytrial, add 1 to the no of successful steps, and set
c                 the no of successive failures to zero
               x = c(17)
               do 400 k = 1, n
                  y(k) = w(k,9)
  400          continue
               c(22) = c(22) + 1.d0
               c(23) = 0.d0
c**************return(with ind=3, xend saved, flag set) if x .eq. xend
               if (x .ne. xend) go to 405
                  ind = 3
                  c(20) = xend
                  c(21) = 1.d0
                  return
  405          continue
               go to 420
  410       continue
c              step not accepted (ind .eq. 6), so add 1 to the no of
c                 successive failures
               c(23) = c(23) + 1.d0
c**************error return (with ind=-3) if hmag .le. hmin
               if (c(14) .gt. c(13)) go to 415
                  ind = -3
                  return
  415          continue
  420       continue
c
c        end stage 4
c
      go to 99999
c     end loop
c
c  begin abort action
  500 continue
c
      write(6,505) ind, tol, x, n, c(13), xend, nw, c(16), c(20),
     +      c(22), c(23), c(24), (y(k), k = 1, n)
  505 format( /// 1h0, 58hcomputation stopped in dverk with the followin
     +g values -
     +   / 1h0, 5hind =, i4, 5x, 6htol  =, 1pd13.6, 5x, 11hx         =,
     +          1pd22.15
     +   / 1h , 5hn   =, i4, 5x, 6hhmin =, 1pd13.6, 5x, 11hxend      =,
     +          1pd22.15
     +   / 1h , 5hnw  =, i4, 5x, 6hhmax =, 1pd13.6, 5x, 11hprev xend =,
     +          1pd22.15
     +   / 1h0, 14x, 27hno of successful steps    =, 0pf8.0
     +   / 1h , 14x, 27hno of successive failures =, 0pf8.0
     +   / 1h , 14x, 27hno of function evals      =, 0pf8.0
     +   / 1h0, 23hthe components of y are
     +   // (1h , 1p5d24.15)                                           )
c
      stop
c
c  end abort action
c
      end
