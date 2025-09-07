c  	constants

      	double precision pi, m_n, h_par, k_b, avog
      	double precision he4mol, gam, rgas
      	parameter(	pi		= 3.141592653589,
     +		 	m_n            	= 0.1045407751,	! meV/(Angstrom/psec)^2
     +		 	h_par    	= 0.6582122,	! meV*psec
     +		 	k_b      	= 0.08617,      ! Boltzman const (meV/K)
     +		 	avog		= 0.6022,       ! Avodadro's no. (10^24 atom/mol)
     +		 	he4mol		= 4.0,		! helium4 (gm/mol)
     +		 	gam		= 1.6666666667,	! adiabatic const
     +		 	rgas		= 8.314)	! gas constant (J/mol*K)

 
      double precision gy, ambBx, ambBy, ambBz,icscale,ocscale,piscale
      parameter(gy     =18.324718)   ! 1/(mGaus*sec)

!!!!!   set Ambient Field which is added to all coil fields in low-field region (Bfield.f)
      parameter( ambBx = 0.1d-3, ambBy=0.d0, ambBz=0.1d-3,icscale=0.0d0,ocscale=0.0d0,piscale=0.0d0)
      double precision gradBxx,gradBxy,gradBzx,gradBzy   ! Gauss/cm  B=ambBx + gradBsx*x + gradBxy*y
      parameter(gradBxx=0.0d-3,gradBxy=0.0d-3,gradBzx=0.0d-3,gradBzy=0.0d-3)     
     
      	double precision k_0, bb
      	parameter(	k_0   		=0.75,         	! 1/Ang
c     +		 	bb    		=3.07729087d-15)!scattering length (m)
     +		 	bb    		=0.326)!scattering length (1E-12 cm)

****	Number of Events   *****
      integer*4, parameter :: nevts = dint(1.0d5) 
      
****	Choose if limiting q<0.1 1/Ang and doing weighting or not
       	logical qlimflg
	double precision qlimit
c	parameter(	qlimflg		=.true.)   	! limit q<0.1 (1/Ang) with weighting
	parameter(	qlimflg		=.false.)  	! use full S(q,w) q<4.8 (1/Ang)
c	parameter(	qlimit		=0.2d0)		! limit for q-values
	parameter(	qlimit		=4.8d0)		! limit for q-values

c     ***Choice of using MCSTAS generated PSM output on NGC****
       logical PSMngcflg, befilterflg
c	   parameter (PSMngcflg = .true. )  !  xy,lambda,angles all done from MCSTAS histograms
	   parameter (PSMngcflg = .false.)  !  need to make choices below for e,xy,lamba dists
	   parameter (befilterflg = .true.) !  removes lambda<3.5angstr
c	   parameter (befilterflg = .false.) !  removes lambda<3.5angstr
cccccccc choices below only matter if PSMngcflg = .false.  ccccccc
****	Lambda choice from lansce spectrum	****
        logical lansceflg
c       parameter (	lansceflg	=.true.)      	!yes
        parameter (	lansceflg	=.false.)     	!no

cccccccc choices below only matter if PSMngcflg and lansceflg = .false.  ccccccc
****	Lambda choice either from distribution or Lambda=5Ang	****
	logical ngcflg
	integer ngcflgnum, lamflg
!	parameter(	lamflg		=0)   		!  lambda from distribution
	parameter(	lamflg		=1)  		!  lambda = 6Angst
!	parameter(	lamflg		=2)  		!  lambda flat
!	parameter(	ngcflg		=.true. )  	!  lambda from NGC
!	parameter(	ngcflg		num = 1 )  	!  lambda from NGC NO filter
!	parameter(	ngcflg		num = 2 )  	!  lambda from NGC Be filter
	parameter(	ngcflg		=.false.)  	!  lambda from NG6
c     ***xy Distribution choice****
	   integer xydist
c	   parameter (xydist=0)   !  Flat xy distribution
!	   parameter (xydist=1)   !  Gradient xy distribution
!	   parameter (xydist=2)   !  NSRII PSM Beam Image xy distribution
!	   parameter (xydist=3)   !  INverted parabola xy distribution
	   parameter (xydist=4)   !  point xy distribution
c     ***Angular Distribution choice****
	   integer angdist
c	   parameter (angdist=0)   !  Isotropic Angular distribution
!	   parameter (angdist=1)   !  Flat Angular distribution
!	   parameter (angdist=2)   !  Cosine Angular distribution
!	   parameter (angdist=3)   !  Angular distribution from J. Cook sim of NGC
	   parameter (angdist=4)   !  Test so that neutrons go straight down beamline
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	
****	Scattering Weighting numerator for limited q<0.1 range	***
	double precision wgtA
	logical sqwhbarflg
c  	removed, set in code now
	parameter (	wgtA 	= 0.00436) 		!uncommented
c	parameter (	wgtA 	= 0.436)  		!artificially increased
       parameter (sqwhbarflg=.true.)   ! adjusts Sqw with hbar factors to get good fit to experimental Sommers cross sections vs. E
c      parameter (sqwhbarflg=.false.)   ! NO adjustment to Sqw model, so not good fit to Sommers sig(E)

      logical isoscatflg
      parameter( isoscatflg=.false.)  !  usual case
c      parameter( isoscatflg=.true.)   ! for testing purposes

c     *****PA weighting of Neutrons*****
     	logical PAflg
c       parameter (PAflg=.true.)  ! PA weight included
       	parameter (PAflg=.false.) ! PA weight = 1

c     *****Pi Coil ON/OFF******
	logical picoilflg
	double precision pilen, pibfield
	parameter (picoilflg=.true.)  	! Pi Coil ON
c	parameter (picoilflg=.false.) 	! Pi Coil Off
c        parameter (pilen = 5.0)   	! no filter!  not using length since shows up in nSpinSim and Bfield.f and moves downstream components
!		parameter (pibfield = 2.1) !  no filter
		parameter (pilen = 4.97)   !  with Be
		parameter (pibfield = 2.6) !  with Be  !  2.27G for 6angstr through 4.97cm for 180deg


c     *****Septum for inguide and incoil*****
        logical sptmflg
	parameter (sptmflg = .true.)  ! septum in place
c	parameter (sptmflg = .false.) ! no septum

c     *****Air or Helium Gas in Beam*****
        logical hebeam
	parameter (hebeam = .true.)   !He gas in beam
c	paramter  (hebeam = .false.)  !air in beam

c	  *****Bfield Gradients****
	    integer bzgradnum,strongField
		character bzfile1*20,bzfile2*20,bzfile3*20,bzfile4*20

c  Chose which field setting below
c		parameter (bzgradnum=0)  	! NO Bz gradient run
		parameter (bzgradnum=1) 	!  sBrot FIELD MAP
c		parameter (bzgradnum=2) 	!  Large Field
c		parameter (bzgradnum=3) 	!  Bzgrad small (0.616mG/m)
c		parameter (bzgradnum=4) 	!  for He gastar, Bzgrad chg with tarstate - !set HEGAS!
		parameter (strongField=0)  	! Spin is not rotated outside low field region
c		parameter (strongField=1)  	! Strong magnetic field is rotating thtot
		
c file names determined by choice above
		parameter (bzfile1='bfield.dat') ! Large Field gradient
c		parameter (bzfile2='Bz_LargeFld.dat') ! Large Field
c		parameter (bzfile3='Bz_run1154.dat') ! Bzgrad small (0.616mG/m)
c		parameter (bzfile4='Bz_run1449.dat') ! for He gastar, Bzgrad chg with tarstate


**********Helium Target properties*********
      double precision T 
c      parameter (T  = 2.10)     ! helium temp. (K)
c      parameter (T  = 4.25)     ! helium temp. (K)
       parameter (T  = 4.24)     ! He temp. (K) for Murad new S(q)
c
      double precision tarx1, tarx2, tarx3, tarx4
	  double precision tary1, tary2, tarz
	  double precision ndenhegas
	  logical hegastar, airtar

c  !!two of the following four should be commented out!!
c  !! for liquid He, both should be false below!!	  
      	parameter (hegastar = .false.)  ! false if liquid He target or air
!      	parameter (hegastar = .true.)   ! true if gas He target
      	parameter (airtar   = .false.)  ! false if liquid He target or He gas
!      	parameter (airtar   = .true.)   ! true if air in target

c    4He gas number density at 273K/1atm is 2.688e19/cc  -KG
c    4He gas number density at 57K/1atm is 1.287e20/cc (ideal gas law) -KG
c    4He gas number density at 4.25K/1atm is 1.727e21/cc (ideal gas law) -KG
c	parameter (ndenhegas = 0.0000269 )  	!! Room Temp He numb. density (1e24 atm/cc)
c   	parameter (ndenhegas = 0.0001287 )  	!! He numb. density (1e24/cc) -KG
c	parameter (ndenhegas = 0.001727 )  	!! He numb. density (1e24/cc)-KG
	parameter (ndenhegas = 0.002484 )  	!! 4K He Gas numb. density (1e24 atm/cc) -KG
c	parameter (ndenhegas = 0.0000 )  	!!

      double precision m_he
      parameter (m_he = 0.418) !4He mass -mev/(Angstrom/psec) from atomic mass ratio to Al -KG

      double precision gtindexx, gtindexy 
      parameter(gtindexx = 0.0, gtindexy = 0.0)  ! rad/Ang
      parameter	(tarx1 = -5.0,    ! cm
     +           tarx2 = -0.15,         ! cm
     +           tarx3 = 0.15,          ! cm
     +           tarx4 = 5.0,           ! cm
     +           tary1 = -5.0,          ! cm
     +           tary2 = 5.0,           ! cm
     +           tarz  = 41.10)          ! cm
c     +           tarz  = 10.110)          ! cm   JUST FOR TEST
******************************************
c     *********Aluminum Window properties**********
      double precision ndenAl, m_al,ndenCu,m_cu
      parameter (ndenAl  = 0.06026d0)     ! Al numb. density (1e24 atm/cc)
c      parameter (ndenAl  = 0.00d0)     ! Al numb. density (1e24 atm/cc)
      parameter (m_al = 2.82d0)                ! Al mass - meV/(Angstrom/psec)^2
      parameter (ndenCu = 0.08492d0)	  !Number density of Cu(1e24/cc)
c      parameter (ndenCu = 0.d0)	  !Number density of Cu(1e24/cc)
      parameter (m_cu = 6.6416d0)     ! meV/(Angstrom/psec)^2
cccccc

      	double precision talx1, talx2, talx3, talx4, taly1, taly2
      	double precision talz1, talz2, talz3, talz4, talz5, talz6, talz7
      	double precision Alentlen, pAlentrad, nAlentrad
      	parameter(talx1 = -6.0,		! (cm)
     +       	talx2 = 6.0,          	! (cm)
     +       	talx3 = -6.0,          	! (cm)
     +         	talx4 = 6.0,           	! (cm)
     +         	taly1 = -6.0,          	! (cm)
     +         	taly2 = 6.0,           	! (cm)
     +         	talz1 = 0.03,          	! (cm)
     +         	talz2 = 0.04,          	! (cm)
     +	 	talz3 = 0.05,		! (cm)
     +		talz4 = 0.20,		! (cm)
     +         	talz5 = 0.032,          	! (cm)
     +         	talz6 = 0.38,          	! (cm)
     +		Alentlen = 0.6858,
     +    	pAlentrad = 10.414,
     +	 	nAlentrad = -10.414)       !cm
************************************
c
*******Air Gaps***********
      	double precision m_air, ndenair, ndenairtar
	parameter (ndenair = 0.000054)  	!! air numb. density (1e24 atm/cc)
c	parameter (ndenair = 0.0)  		!! air numb. density (1e24 atm/cc)
ccccccc
c	parameter (ndenairtar = 0.000054)  	!! air numb. density (1e24 atm/cc)
c	parameter (ndenairtar = 0.00019)  	!! 77K N numb. density (1e24 atm/cc)
	parameter (ndenairtar = 0.00)  		!! air numb. density (1e24 atm/cc)
ccccccc
      	parameter (m_air = 1.57)     		!! air mass - meV/(Angstrom/psec)^2

      	double precision tairx1,tairx2,tairy1,tairy2,tairz0,tairz1,tairz2
	double precision tairz3,tairz4,tairz5,tairz6, tairx0,tairy0
	double precision tairz7,tairz8,tairz9,tairz10,tairz11,tairz12
	double precision tairz13,tairz14,tairz15,tairz16,tairz17
      	parameter(tairx0 = 5.0,    	! normally set to size of PSM
     +		tairy0  = 5.0,		! (cm)
     +		tairx1  = -6.0,      	! (cm)
     +		tairx2  = 6.0,         	! (cm)
     +        	tairy1  = -6.0,       	! (cm)
     +         	tairy2  = 6.0,         	! (cm)
!    +         	tairz0  = 150.0,      	! (cm)gap BEFORE PSM     -- with Be Filter
     +         	tairz0  = 0.0,         	! (cm)gap BEFORE PSM       -- NO Be Filter
     +         	tairz1  = 4.10,        	! (cm)gap after PSM
     +	       	tairz2  = 1.15,    	! (cm)gap after guide extensions 
     +         	tairz3  = 1.95,       	! (cm)gap after input coil
     +         	tairz4  = 2.86,       	! (cm)gap after endcap
     +         	tairz5  = 3.67, 	! (cm)gap after vaccum entrance window
     +         	tairz6  = 3.95,        	! (cm)gap after LN shield entrance window
     +         	tairz7  = 1.52,        	! (cm)gap after LHe entrance window
     +         	tairz8  = 1.50,        	! (cm)gap after 1st target
c     +         	tairz8  = 100.50,        	! (cm)gap after 1st target  !!!! TEST TO MAKE SOLID ANGLES VERY DIFF
     +         	tairz9  = 4.97,        	! (cm)gap after pi coil upstream
     +         	tairz10 = 1.37,      	! (cm)gap after pi coil downstream
     +         	tairz11 = 1.40,        	! (cm)gap after 2nd target
     +         	tairz12 = 6.08,        	! (cm)gap after LHe exit window
     +         	tairz13 = 19.95,      	! (cm)gap after LN shield exit window
c     +         	tairz13 = 0.95,      	! (cm)gap after LN shield exit window   !!TEST TO MAKE SOLIDANGLES VERY DIFF
     +         	tairz14 = 0.57,        	! (cm)gap after vacuum exit window
     +         	tairz15 = 3.36,        	! (cm)gap after output coil endcap
     +         	tairz16 = 3.41,        	! (cm)gap after output coil
     +         	tairz17 = 7.00)		! (cm)gap after super mirror analyzer
*******Air Gaps***********
c
c	  **BEAM GEOMETRY*******************************
c
c     ******Initial x, y-range and guide length***
	  double precision gl,xin,yin,ng6indexx,ng6indexy,ng6index
	parameter(	xin = 5.5,           	! starting x range (E-W) after PSM
     +           	yin = 5.5,            	! Height of PSM is 4.5cm  -KG    
     +           	gl = 53.34,           	! Guide length in z-direction (cm)
c     +           	ng6index = 0.002)    	! NG6 index rad/ang
c     +           	ng6index = 0.0032)    	! NG6 index rad/ang
c
     +           	ng6indexx = 0.0035, 	! Greater divergence after PSM -KG 
     +          	ng6indexy = 0.010)     	! Greater divergence after PSM -KG 
c     +          	ng6indexx = 0.0000001, 	! for no divergence -rcm
c     +          	ng6indexy = 0.0000001) 	! for no divergence -rcm 
c

c    *******Input Guide********
      	double precision gx1, gx2, gx3, gx4, gy1, gy2
	logical reflectflg
c	parameter (reflectflg = .true.)   ! true reflectivity = 1.0
      	parameter (reflectflg = .false.)  ! false reflectivity varies
      	parameter(	gx1 = -5.0,		!cm  not -2.20 -KG 
     +           	gx2 = -0.15,          	!cm, 3mm septum -rcm 
     +           	gx3 = 0.15,         	!cm, 3mm septum -rcm
     +           	gx4 = 5.0,          	!cm  
     +           	gy1 = -5.0,         	!cm  not -2.25 -KG 
     +          	gy2 = 5.0)          	!cm not 2.25 -KG   


c    *******Input Coil*********
   	double precision ix1, ix2, ix3, ix4, iy1, iy2, iz0, ginx, giny
c	parameter(	ginx = 0.00,      	!Input critical index (rad/Ang)
c     +           	giny = 0.00,      	!Input critical index (rad/Ang)
	parameter(	ginx = 0.00346,      	!Input critical index (rad/Ang)
     +           	giny = 0.00346,      	!Input critical index (rad/Ang)
     +	         	ix1 = -5.0,         	!not -2.2 -KG  
     +           	ix2 = -0.15,        	!cm, 3mm septum -rcm   
     +           	ix3 = 0.15,         	!cm, 3mm septum -rcm
     +           	ix4 = 5.0,          	!cm 
     +           	iy1 = -5.0,         	!cm  not -2.25 -KG 
     +           	iy2 = 5.0,     		!cm  not 2.25 -KG
     +          	iz0 = 202.64)          	!cm
c
c    ****************************
c
c    *******Output Coil**********
      double precision ox1, ox2, ox3, ox4, oy1, oy2, oz, gotx, goty
c      parameter (gotx = 0.00,	! Output critical index (rad/Ang)
c     +           goty = 0.00,       	! Output critical index (rad/Ang)
      parameter (gotx = 0.00346,	! Output critical index (rad/Ang)
     +           goty = 0.00346,       	! Output critical index (rad/Ang)
     +           ox1 = -5.0,       	! cm  not -3.3255 -KG
     +           ox2 = -0.15,       	! cm
     +           ox3 = 0.15,        	! cm
     +           ox4 = 5.0,       	! cm  not 3.3255 -KG
     +           oy1 = -5.0,      	! cm  not -2.575 -KG
     +           oy2 = 5.0,        	! cm  not 2.575 -KG
     +           oz  = 127.64)         	! cm
c     +           oz  = 20.)         	! cm   !!!!TEST FOR SOLID ANGLE TEST
c
c    ******Analyzing Super Mirror (ASM)****
      double precision sx1, sx2, sy1, sy2, sz, detrad, gsm
      parameter	(sx1 = -5.0,        	! cm 
     +           sx2 = 5.0,         	! cm
     +           sy1 = -5.0,        	! cm
     +           sy2 = 5.0,         	! cm
     +           sz  = 42.00,        	! cm
     +           detrad = 7.1,       	! cm
     +           gsm = 0.003)        	! (Pyrex glass+ Ti + Ga) rad/Ang



c  *****ion detector****
      	double precision winrad, ionz1, ionz2, ionz3, ionz4
      	double precision partpreshe, partpresar

	parameter(winrad = 9.144)		!cm should be 9.144

!   NGC run
	parameter(
     +  ionz1 = 2.54, 		!cm
     +  ionz2 = 8.68, 		!cm
     +  ionz3 = 16.11,		!cm
     +  ionz4 = 24.66  )	!cm   

!  from PRC paper NSRII run with old ion detector
!0, 2.8, 9.2, 16.6, 25.3
!	parameter(
!     +  ionz1 = 2.8, 		!cm
!     +  ionz2 = 9.2, 		!cm
!     +  ionz3 = 16.6,		!cm
!     +  ionz4 = 25.3  )	!cm   

!  ???NSRII run with old ion detector????
!0, 2.385, 6.195, 13.815, 25.245
!	parameter(
!	ionz1 = 2.385, 		!cm
!	ionz2 = 6.195, 		!cm
!	ionz3 = 13.815,		!cm
!	ionz4 = 25.245	)  	!cm   

!  ????diff NSRII run with old ion detector
!0, 2.91, 10.17, 19.32, 25.25
!	parameter(
!	ionz1 = 2.91, 		!cm
!	ionz2 = 10.17, 		!cm
!	ionz3 = 19.32,		!cm
!	ionz4 = 25.25  )	!cm   

!       for equal percent spacing NGC,5Angstroms
!	parameter(
!	ionz1 = 1.2, 		!cm
!	ionz2 = 3.0 ,		!cm
!	ionz3 = 6.0,		!cm
!	ionz4 = 24.66  )	!cm   

	parameter(
c     +	partpreshe 	= .37,
     +  partpreshe 	= .22,   	! NSR II
!	partpresar  	= 0.001,
     +  partpresar  	= 3 )


c 	Axial offset for bfield 
	double precision xBfieldOff,yBfieldOff
	parameter(xBfieldOff = 0)
	parameter(yBfieldOff = 0)


	logical hfflg
	parameter(hfflg=.false.)


