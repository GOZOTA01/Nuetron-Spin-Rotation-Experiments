	subroutine sBrot(E,XI,YI,ZI,THETAII,PHII,XF,YF,ZF,ThetaF,Phif,
     &			ztargstart,myfieldfile,P,R,a,fieldfile2)
c
c	Bec Sep2020   dz, bgrad, etc haven't been used in couple years...removing them	
c	Edited May 2018 by Sam Infanger
c	 April 2016, BEC, modified totdep_phic and bfield to allow choice of reference 
c            axis for theta and phi of inital spin and allow input of B-field from a file.
c	     Defining step by 1/10 of 2pi around in a Larmor precession, so omegaL*dt = PI/5
c            Then the linear step, dz, is given by dz=v*dt=v*PI/5/omegaL=PI*v/(5*gyro*B)
c
c	    same as totdep_phi except now can put another collimator in
c	    middle of beamline.  That is instead of just having starting
c	    and stopping X and Y values, now you can make sure that a 
c	    given trajectory passes through some XM and YM at ZM.
C
C		Written 5/22/91 by W. Scott Wilburn, TUNL
C		Based on the FORTRAN program SPINX by D. Bowman, LANL
C		Modified for small offaxis bfield calculation via subroutine
C		BFIELD.FOR by Bret Crawford
C
C		****************************************************************
C		Must be linked with an external subroutine BFIELD(R,B) with 
C		a 3 component R*8 giving the position at which to calculate the
C		magnetic field (cm) and B a 3 component R*8 giving the 
C		components of the B field returned by the subroutine (G).
C		****************************************************************
C
C		This program calculates the spin precession of neutrons
C		traveling through a gap in two solenoids.  They can travel at
C		a given radial distance a along the z-axis or at frandom\
C		trajectories.  The initial
C		two-component spinor wave function is calculated from the
C		initial theta and phi given.  Each neutron is transported
C		along its chosen trajectory in small enough steps that the
C		spin precession can be treated as infinitesimal.  A rotation
C		matrix is generated for each step and the product of all
C		the rotation matrices is computed along the way.  At the end
C		of the trajectory, a rotation matrix representing the entire
C		precession has been generated.  The step size is determined
C		with a convergence criteria which compares the result of
C		taking a whole step with that of taking two half steps.  If
C		the components of the rotation matrix from the two step sizes
C		differ by more than specified amount, the step size is cut in
C		half and we try again.  At the end of each trajectory, the 
C		final wave function is calculated by operating the initial
C		wave function with the rotation operator.  Spin projection
C		expectation values are calculated along the x,y, and z axes and
C		written to ouput as are the final angles, theta and phi.  
C		The output contains both these final results
C		and the results at each z-value curing the trip.  The program
C		TOTDEPOL_PLOT.FOR uses the output of this program to plot the 
C		amount of depolarization as a function of axial distance.  The
C		input file is BFIELD.IN.
C		
C		Input variables:
C
C			E	Neutron energy (eV) See subroutine Rot
C			PHII	Initial azimuthal spin angle (mrad)
C			THETAI	Initial polar spin angle (mrad)
C			THETAII Initial polar spin angle (mrad)
C			ZF	Final z position (cm)
C			ZI	Initial z position (cm)
C			XI	Initial x position (cm)
C			XF	Final X position (cm)
C			YI	Initial Y position (cm)
C			YF	Final Y position (cm)
C			ZM	Middle Z position, used for the collimator (cm)
C			YM 	Middle Y position, used for the collimator (cm)
C
C
C		Output variables:
C
C			P		Spin projections on cartesian axes
C			PHIF		Azimuthal angle of final spin (mrad)
C			RF		Final position vector (cm)
C			RI		Initial position vector (cm)
C			THETAF		Polar angle of final spin (mrad)
C
C		Internal Variables:
C
C			CVG		Indicates whether convergence has
C						occured in calculating rotation
C						matrix
C			D		Step size in fraction of total distance
C						to be travelled
C			DELTA_T		Difference in component of rotation
C						matrices calculated using step
C						D and D/2
C			DX		Magnitude of step (cm)
C			J		Counting variable
C			K		Counting variable
C			L		Counting variable
C			F		Fraction of total distance travelled
C			PSIF		Final spinor wavefunction
C			PSII		Initial spinor wavefunction
C			R		Position vector (cm)
C			R_B		Position vector at which to calculate
C						B field for full step (cm)
C			RD		Position difference vector, RF-RI (cm)
C			RL_B		Position vector at which to calculate
C						B field for lower half step (cm)
C			RU_B		Position vector at which to calculate
C						B field for upper half step (cm)
C			S		Spin direction vector
C			T		Rotation matrix for a full step
C			TL		Rotation matrix for lower half step
C			TT		Matrix product TU*TL
C			TU		Rotation matrix for upper half step
C			U		Rotation matrix for entire trajectory
C
C		Constants:
C
C			stepFrac	Initial stepsize in fraction of total distance
C					to be moved (1)
C			EPS		Convergence criteria on roation matrices for
C					determining stepsize
C			I		Imaginary unit
C			PI		Pi
C			X		Unit vector in x direction
C			Y		Unit vector in y direction
C			Z		Unit vector in z direction
C
		IMPLICIT NONE
		include 'const.inc'
C
		LOGICAL*1 CVG
		logical readflg,ocflg

c		common/bfld/readflg,db,zgrad2,zgrad3
		common/bfld/readflg,ocflg

		integer*4 j,k,l,n,a
c		REAL*8 myz(1000), mythetax(1000), mythetay(1000), mythetaz(1000)
C
c		INTEGER*4 J,K,L,N,v,g,q,w,zgrad2(100000),zgrad3(100000),   !SMI
c     >                    ce,le						!SMI
                REAL*8 env, rsqrt						!SMI
		REAL*8 D,stepFrac,DELTA_T,DX,E,EPS,F
     >			,R(3),R_B(3),RD(3),RF(3),RI(3),RL_B(3),RU_B(3)
     		REAL*8 YI,YF,XI,XF,X(3),Y(3),Z(3),ZF,ZI,zgrad(100000)	!SMI
     >			,ZM,YM,midY,Brot(3),Pstart
     >                  ,Poavg(3),Poerr(3),Po(3),db(100000)		!SMI
C    		REAL*8 C,YF,XI,XF,X(3),Y(3),Z(3),ZF,ZI	!SMI
c     >			,ZM,YM,midY,Brot(3),Po(3)		!SMI
C
c		REAL*8 P(3),PHII,S(3),PSUM(3),PSQSUM(3),THETAI,PHIF,THETAF,
c     >			depol,THETAFDEV,PHIFDEV,PDEV(3),depoldev,PAVE(3),
c     >			THETAII,PINIT(3),meandev,ztargstart,BGRAD(100000), !SMI
c     >			dz(100000),GRAD(100000) !SMI
		REAL*8 P(3),PHII,S(3),PSUM(3),PSQSUM(3),THETAI,PHIF,THETAF,
     >			THETAII,ztargstart
C
		COMPLEX*16 I
C
		COMPLEX*16 PSIF(2),PSII(2),Ttot(2,2),TL(2,2),TT(2,2)
		COMPLEX*16 TU(2,2),U(2,2)
C
		Character comments*80,filei*20,fileo*20,AXIS*1,myfieldfile*80,fieldfile2*80
		

		SAVE /bfld/
   		open (action='write', file='thetax.txt', unit=8, status='replace')
		open (action='write', file='thetay.txt', unit=1, status='replace')
		open (action='write', file='thetaz.txt', unit=2, status='replace')
		open (action='write', file='quadrants.txt', unit=11, status='replace')
		open (action='write', file='Bfieldsx.txt', unit=18, status='replace')
		open (action='write', file='Bfieldsy.txt', unit=28, status='replace')
		open (action='write', file='Bfieldsz.txt', unit=38, status='replace')

c		Open(unit=99,file=fieldfile2,status='old',action='read')

C
C		Read in Location and Geometry
C
		

		DATA EPS	/1.0D-5/   ! orig 1E-5
		DATA I		/(0.0,1.0)/
		DATA X		/1.0,0.0,0.0/
		DATA Y		/0.0,1.0,0.0/
		DATA Z		/0.0,0.0,1.0/
		DATA N		/1/
		DATA AXIS	/'Z'/
		!june6
		!write(6,*) E,XI,YI,ZI,THETAII,PHII,XF,YF,ZF,ztargstart
c  take an intial step that is stepFrac, which is fraction of complete step thru object
		if ((ZF-ZI).lt.1.0) then  ! if total zx is <1cm, take bigger fractional step
c			!june13
			stepFrac = 1.0d0/3.0d0     ! orig 1/3
		else
			stepFrac = 0.2d0
		end if
c
c     for the output coil field, need to take small steps all the time, so start with small step
c		stepFrac=1.0d-4

100		FORMAT (' Initial z position (cm)            ',F7.2)
110		FORMAT (' Final z position (cm)              ',F7.2)
120		FORMAT (' Initial X position (cm)            ',F7.2)
130		FORMAT (' Final X position (cm)              ',F7.2)
134		FORMAT (' Initial Y position (cm)            ',F7.2)
136		FORMAT (' Final Y position (cm)              ',F7.2)
137		FORMAT (' Middle Z position (cm)             ',F7.2)
138		FORMAT (' Middle Y position (cm)             ',F7.2)
140		FORMAT (' Neutron energy (eV)                 ',1pE9.2)
150		FORMAT (' Initial spin angle theta (degrees)   ',F7.2)
160		FORMAT (' Initial spin angle phi (degrees)     ',F7.2)
164		FORMAT (' Initial spin angle ref. axis	       ',a1)
161		FORMAT (' Earth Magnetic Field (Gauss)	       ',F10.7)
162		FORMAT (' Angle of Earth Field (degrees)       ',F7.2)
170		FORMAT (/,9x'Z',11x,'Px',11x,'Py',11x,'Pz'
     >                   ,11x,'Bx',11x,'By',11x,'Bz')
C
C
cccc  this line just for skipping rotations	
cccc       goto 201

		THETAI=THETAII/1000.0D0
		PHII=PHII/1000.0D0

c	write(6,*) ' theta ',thetaii, ' phi',phii

C		Calculate initial spin wavefunction, PSII, from initial theta
C		and phi -- where theta is from axis of initial polarization x,y, or z
c               and phi is the azimuthal angle, so if in input file theta=90, phi=90 with X 
C               as the polarization axis, then the polarization would be along the z-axis,
c               This is done so one can specify 0,0,Y for instance, so polarized along the Y-direction and
c               then let the code choose random phi directions, that is, random phases for the y-polarization
c
c	If (AXIS.eq.'X') then
c	write(6,*) ' X '
c		PSII(1)=(DCOS(THETAI/2.0D0)*DCOS(PHII/2.0D0)+
c     >                  DSIN(THETAI/2.0D0)*DSIN(PHII/2.0D0)+
c     >		    I*( DSIN(THETAI/2.0D0)*DCOS(PHII/2.0D0)+
c     >              DCOS(THETAI/2.0D0)*DSIN(PHII/2.0D0)))/DSQRT(2.0D0)
c		PSII(2)=(DCOS(THETAI/2.0D0)*DCOS(PHII/2.0D0)-
c     >                  DSIN(THETAI/2.0D0)*DSIN(PHII/2.0D0)+
c     >		    I*( DCOS(THETAI/2.0D0)*DSIN(PHII/2.0D0)-
c     >              DSIN(THETAI/2.0D0)*DCOS(PHII/2.0D0)))/DSQRT(2.0D0)
c	else if (AXIS.eq.'Y') then
c	write(6,*) ' Y '
c	  	PSII(1)=(DCOS(PHII/2.0D0)*
c     >			(DCOS(THETAI/2.0D0)+DSIN(THETAI/2.0D0))-
c     >			I*DSIN(PHII/2.0D0)*
c     >			(DSIN(THETAI/2.0D0)-DCOS(THETAI/2.0D0)))/DSQRT(2.0D0)
c	  	PSII(2)=(DSIN(PHII/2.0D0)*
c     >			(DCOS(THETAI/2.0D0)+DSIN(THETAI/2.0D0))+
c     >			I*DCOS(PHII/2.0D0)*
c     >			(DSIN(THETAI/2.0D0)-DCOS(THETAI/2.0D0)))/DSQRT(2.0D0)
c	else if (AXIS.eq.'Z') then
c	write(6,*) ' Z '
		PSII(1)=DCOS(THETAI/2.0D0)*DCOS(PHII/2.0D0)+
     >			I*DCOS(THETAI/2.0D0)*DSIN(PHII/2.0D0)
		PSII(2)=DSIN(THETAI/2.0D0)*DCOS(PHII/2.0D0)-
     >			I*DSIN(THETAI/2.0D0)*DSIN(PHII/2.0D0)
c	else
c		write(6,*) ' oops, no axis specified in input '
c	end if
c		write(6,*) thetai,phii
c		write(6,'(a10,f7.3,"  + ",f7.3," I")') 'PSII 1',psii(1)
c		write(6,'(a10,f7.3,"  + ",f7.3," I")') 'PSII 2',psii(2)

		CALL PROJ(PSII,X,P(1))
		CALL PROJ(PSII,Y,P(2))
		CALL PROJ(PSII,Z,P(3))
C		write(2,*)
C		write(2,77) ' Initial Px ',P(1)
c		write(2,77) ' Initial Py ',P(2)
c		write(2,77) ' Initial Pz ',P(3)
c		write(2,*)
c		write(6,*)
c		write(6,77) ' Initial Px ',P(1)
c		write(6,77) ' Initial Py ',P(2)
c		write(6,77) ' Initial Pz ',P(3)
c		write(6,*)



77       format(2x,a13,f11.9)
C		WRITE (2,170)
C		Initialize variables.  PSUM are the spin projection
C		expectation values for each axis summed over all trajectories
C		calculated.  PSQSUM are the sums of the squares of these
C		values.  They are used at the end to calculate averages and
C		standard deviations.
C
C		Calculate N random trajectories, propagate a neutron spin
C		wavefunction along each trajectory, and calculate expectation
C		values of the spin projections on the x, y, and z axes.
C

		CALL PROJ(PSII,X,Po(1))
		CALL PROJ(PSII,Y,Po(2))
		CALL PROJ(PSII,Z,Po(3))

C
C		Choose a random trajectory for a neutron starting at
C		z=ZI and ending at z=ZF.  Starting X values are
C		randomly chosen in the range -XI<=X<=XI ditto for Y.
C		Final X values are
C		chosen by picking random trajectories with a cosine
C		probability distribution such that the pass through the
C		region -XF<=X<=XF.  Ditto for Y.  RI and
C		RF are the initial and final vectors so chosen and RD is
C		RF-RI.
C
99		CONTINUE

C
              	RI(1) = XI + xBfieldOff
		RI(2) = YI + yBfieldOff
		RI(3) = ZI - ztargstart
		!june18
		!WRITE(6,*)'cursora ', RI(3),ZI,ZF,ztargstart
C
              	RF(1) = XF + xBfieldOff
               	RF(2) = YF + yBfieldOff
		RF(3) = ZF - ztargstart

		DO K=1,3
			RD(K)=RF(K)-RI(K)   ! RD is the total distance needed to be traveled
		END DO
			



c		ZM = (ZF - ZI)/2 + ZI
c		YM = (YF - YI)/2 + YI
C		Check that the trajectory passes through middle collimators
c		midY = RI(2) + ((RF(2)-RI(2))/(RF(3)-RI(3)))*(ZM-RI(3))
c		If (ABS(midY).LE.ABS(YM)) then
c			GOTO 111
c		ELSE
c			GOTO 99
c		ENDIF
111		CONTINUe
			
C
			
C
C		Transport a neutron spin along the trajectory from RI to
C		RF starting with initial spin wave function PSII.  PSIF
C		is the calculated final spin wave function.
C		Set initial values of position (R), fraction of
C		trajectory travelled (F), and fractional incremental distance(D).
C               RD is the total distance needed to be traveled.
C
		DO K=1,3
			R(K)=RI(K)
		END DO
C
		D=stepFrac
		F=0.0
C
c		write(6,*) r(1),r(2),r(3)
C		Initialize U, a 2x2 matrix corresponding to the
C		cumulative rotation of the neutron spin
C
		U(1,1)=(1.0,0.0)
		U(1,2)=(0.0,0.0)
		U(2,1)=(0.0,0.0)
		U(2,2)=(1.0,0.0)
		DO WHILE (F.LT.1.0)
C
C			Check to see if we are at the last step and it
C			will take us past the end of the trajectory.  If
C			so, adjust the last step to end at the right
C			point.
C
			IF ((F+D).GT.1.0) THEN
				D=1.0D0-F
			END IF
C
C			Calculate the rotation matrix for a whole step
C			(T) and for the interval taken as two half steps
C			(TL and TU).  If T differs from the product of
C			TL and TU by more than EPS, it is not a good
C			approximation to approximate the step with an
C			infintesimal rotation, and we cut the step size
C			in half and try again.
C
			DO K=1,3
				R_B(K)=R(K)+(D/2.0D0)*RD(K)
				RL_B(K)=R(K)+(D/4.0D0)*RD(K)
				RU_B(K)=R(K)+(0.75D0*D)*RD(K)
			END DO

			DX=D*DSQRT(RD(1)**2+RD(2)**2+RD(3)**2)
			!cursor
c			WRITE(6,*) ; write(6,*) 'calling rot full', R_B(3)
			CALL ROT(R_B,DX,E,Ttot,Brot,myfieldfile,fieldfile2)
c			write(6,*) 'calling rot halved'
			CALL ROT(RL_B,DX/2.0,E,TL,Brot,myfieldfile,fieldfile2)
			!write(6,*) 'calling rot halved'
			CALL ROT(RU_B,DX/2.0,E,TU,Brot,myfieldfile,fieldfile2)

C			Calculate TT and the difference between it and
C			T, DELTA_T
C
			CALL MULT(TL,TU,TT)
c			write(6,*)'fail'
			CVG=.TRUE.
C
			DO K=1,2
				DO L=1,2
					DELTA_T=ABS(Ttot(K,L)-TT(K,L))
c					write(6,*) 'stepFrac, d, delta_t ',stepFrac,D,Delta_T
					IF (DELTA_T.GT.EPS) THEN
						CVG=.FALSE.
					END IF
				END DO
			END DO
C
C			If the convergence criteria is not met, cut the
C			step size in half.  If it is met, increment the
C			fraction of the distance travelled, F, by the
C			step taken, D, multiply the cumulative rotation
C			matrix, U, by the rotation performed in this
C			step, T, and set the step size back to its
C			original value.
C
			IF (CVG.EQV..FALSE.) then
				D=D/2.0D0
c					write(6,*) 'it didnt converge'
c					write(6,'(a6,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3)') 'No CVG', 
c     &    				R(3), D, d*(zf-zi), F 
			ELSE 					
				F=F+D
				!june19
				!WRITE(6,*) 'cursorb ',R(3),D,RD(3)
				DO K=1,3
					R(K)=R(K)+D*RD(K)
				END DO
				!june19
				!WRITE(6,*) 'cursorc ',R(3),D,RD(3)
				CALL MULT(U,Ttot,U)		
c				le = floor(R(3))+20			!SMI
c				env = ZF-ZI				!SMI
c				write(6,*) R(3),le,zgrad2(le)
c				write(6,*) zgrad3(zgrad2(le)),db(zgrad3(zgrad2(le)))
				D=stepFrac
c				write(6,'(a6,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3)') 'CVG', 
c     &    			R(3), D, d*(zf-zi), F 

c        			D=RSQRT(dabs((9.0D-4*RSQRT(2*E/m_n))/				!SMI
c     &  			(gy*db(zgrad3(zgrad2(le))))))/env 				!SMI

c				if(D.gt.1.0D0) then						!SMI
c					D=5.0d-1
c				end if
c				write(6,*) d

C       			From here to 250 is to see what's
C       			happening to the neutron's precession during its flight through field.

                              	DO K=1,2
                                  	PSIF(K) = (0.0,0.0)
                                END DO
                                DO K=1,2
                                        DO L=1,2
                                        	PSIF(K)=PSIF(K)+U(K,L)*PSII(L)
                                        END DO
                                END DO
	      		        CALL PROJ(PSIF,X,P(1))
          		      	CALL PROJ(PSIF,Y,P(2))
         		      	CALL PROJ(PSIF,Z,P(3))

				

c				write(6,*) r(1),r(2),r(3)
c            			WRITE(6,250) R(3),P(1),P(2),P(3),
c     >				Brot(1),Brot(2),Brot(3)

cccccccc 			uncomment to write out the rotation angles in each diredcion for each step in z-direction   Ccccc
c            			WRITE(6,250) R(3),asin(2.0*P(1)-1.0),
c     &   			asin(2.0*P(2)-1.0),asin(2.0*P(3)-1),Brot(1),Brot(2),Brot(3)
c	write(6,*) ' new z value ',R(3), asin(2.0*P(1)-1.0)
	write(8,*) R(3), asin(2.0*P(1)-1.0)
	write(1,*) R(3), asin(2.0*P(2)-1.0)
	write(2,*) R(3), asin(2.0*P(3)-1.0)
	write(18,*) R(3), Brot(1)
	write(28,*) R(3), Brot(2)
	write(38,*) R(3), Brot(3)

		write(11,*) XI,YI
		
	
cccccccccccccccc
cccccccc 			uncomment to write out the projection operators each diredcion for each step in z-direction   Ccccc
				!june
cc   to screen				WRITE(6,250) R(1),R(2),R(3),P(1),P(2),P(3),Brot(1),Brot(2),Brot(3)
c				if (ocflg) then
c					if ((R(1)>0).and.(R(2)>0)) WRITE(170,250) R(1),R(2),R(3),P(1),P(2),P(3),Brot(1),Brot(2),Brot(3)
c					if ((R(1)<0).and.(R(2)>0)) WRITE(171,250) R(1),R(2),R(3),P(1),P(2),P(3),Brot(1),Brot(2),Brot(3)
c					if ((R(1)<0).and.(R(2)<0)) WRITE(172,250) R(1),R(2),R(3),P(1),P(2),P(3),Brot(1),Brot(2),Brot(3)
c					if ((R(1)>0).and.(R(2)<0)) WRITE(173,250) R(1),R(2),R(3),P(1),P(2),P(3),Brot(1),Brot(2),Brot(3)
c				else
c					if ((R(1)>0).and.(R(2)>0)) WRITE(174,250) R(1),R(2),R(3),P(1),P(2),P(3),Brot(1),Brot(2),Brot(3)
c					if ((R(1)<0).and.(R(2)>0)) WRITE(175,250) R(1),R(2),R(3),P(1),P(2),P(3),Brot(1),Brot(2),Brot(3)
c					if ((R(1)<0).and.(R(2)<0)) WRITE(176,250) R(1),R(2),R(3),P(1),P(2),P(3),Brot(1),Brot(2),Brot(3)
c					if ((R(1)>0).and.(R(2)<0)) WRITE(177,250) R(1),R(2),R(3),P(1),P(2),P(3),Brot(1),Brot(2),Brot(3)
c				end if
cccccccccccccccc
c				WRITE(6,*) R(3),D,F
c				write(6,'(a6,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3)') 'RotStep',
c     &    			R(3), D, d*(zf-zi), F
				!june18
c				if((R(3).GT.46.8995).and.(47.1496.GT.R(3)))then
c					DO K=1,3
c						!WRITE(6,*) R_B(K)
c					END DO
c					!WRITE(6,*)DX,E,Ttot
c				end if

  250                           FORMAT(2X,9(2X,F12.8))
			END IF

		END DO
c	write(15,*) R(1), R(2), asin(2.0*P(1)-1.0)
c	write(16,*) R(1), R(2), asin(2.0*P(2)-1.0)
c	write(17,*) R(1), R(2), asin(2.0*P(3)-1.0)

	
C		Calculate the average theta and phi for final spin wavefunction,
C		with standard deviations
C
C		Get unit vector for spin direction.
C
		DO K=1,3
			S(K)=2.0D0*P(K)-1.0D0
		END DO
C


		THETAF=ACOS(S(3))
		
C
		IF ((S(1).EQ.0.0).AND.(S(2).EQ.0.0)) THEN
			PHIF=0.0
			
		ELSE
			IF (S(2).GE.0.0) THEN
				if ((S(1)**2+S(2)**2).le.0.0d0) 
     &                   	write(6,*) 'xyproj<0',(S(1)**2+S(2)**2)
				PHIF=DACOS(S(1)/DSQRT(S(1)**2+S(2)**2))
			ELSE
				PHIF=2.0D0*PI-DACOS(S(1)/DSQRT(S(1)**2+S(2)**2))
			END IF
		END IF


		PHII = PHII*1000D0    ! convert angles back to mrad for main code
		PHIF = PHIF*1000D0
		THETAI = THETAI*1000D0
		THETAF = THETAF*1000D0
		!june6
		!write(6,*)	'output ' ThetaF,PhiF
		if (isnan(THETAF).or.isnan(PHIF)) write(6,*) ' isnan sBrot R',
     &          R(1),R(2),R(3)

c		write(6,*) 'end of sBrot fine',RF(3)

C
C		Write final results to output file
C
C		WRite (2,400)
C		WRITE (2,401) Po(1)
C		WRITE (2,402) Po(2)
C		WRITE (2,403) Po(3)
C		WRite (2,404)
C		WRITE (2,405) P(1)
C		WRITE (2,410) P(2)
C		WRITE (2,420) P(3)
C		WRITE (2,430) THETAF
C		WRITE (2,440) PHIF

c		WRite (6,400)
c		WRITE (6,401) Po(1)
c		WRITE (6,402) Po(2)
c		WRITE (6,403) Po(3)
c		WRite (6,404)
c		WRITE (6,405) P(1)
c		WRITE (6,410) P(2)
c		WRITE (6,420) P(3)
c		WRITE (6,430) THETAF
c		WRITE (6,440) PHII





  400		FORMAT (//,' Starting Averages')
  401		Format (/,'    PX = ',F6.4)
  402		FORMAT ('    PY = ',F6.4)
  403		FORMAT ('    PZ = ',F6.4,//)
  404		FORMAT (//,' Final Results')
  405		Format (/,'    PX = ',F6.4)
  410		FORMAT ('    PY = ',F6.4)
  420		FORMAT ('    PZ = ',F6.4,//)
  430		FORMAT (' THETA = ',f9.4)
  440		FORMAT ('   PHI = ',f9.4,/)


C	close(2)
	close(8)
	close(1)
	close(2)
	close(11)
	close(18)
	close(28)
	close(38)

	IF (a == 0) then
	call execute_command_line('gnuplot -p ' // 'Allplot.plt')
c	call execute_command_line('gnuplot -p ' // 'Bfieldsplots.plt')
c	call execute_command_line('python3 ' // 'plotBfields.py')
	end if
c	call execute_command_line('gnuplot -p ' // 'scatterplot.plt')
c	call execute_command_line('gnuplot -p ' // 'plot.plt')
c	call execute_command_line('gnuplot -p ' // 'plot2.plt')
c	call execute_command_line('gnuplot -p ' // 'plot3.plt')
	goto 202 
 201    continue
	THETAF = -1234567.
	PHIF = -1234567.
 202    continue
	END
C
C
C
	SUBROUTINE MULT(U1,U2,U3)
C
C		Multiplies two spin rotation operators to obtain a third,
C		U3=U2*U1
C
C		Input variables:
C
C			U1	First matrix to be multiplied
C			U2	Second matrix to be multiplied
C
C		Output variable:
C
C			U3	Matrix product U2*U1
C
C		Internal variables:
C
C			K	Counting variable
C			L	Counting variable
C			M	Counting variable
C			UT	Matrix product U2*U1
C
		IMPLICIT NONE
C
		INTEGER K,L,M
C
		COMPLEX*16 U1(2,2),U2(2,2),U3(2,2),UT(2,2)
C
C		Initialize UT, which will temporarily be the product.  This
C		allows U3 to be the same variable as U1 or U2, if desired
C
		DO K=1,2
			DO L=1,2
				UT(K,L)=(0.0,0.0)
			END DO
		END DO
C
C		Calculate UT=U2*U1
C
		DO K=1,2
			DO L=1,2
				DO M=1,2
					UT(K,L)=UT(K,L)+U2(K,M)*U1(M,L)
				END DO
			END DO
		END DO
C
C		Set U3=UT
C
		DO K=1,2
			DO L=1,2
				U3(K,L)=UT(K,L)
			END DO
		END DO
C
		RETURN
	END
C
C
C
	SUBROUTINE ROT(R1,DX1,E1,T1,B,myfieldfile,fieldfile2)
C
C		Constructs the approximate rotation operator for an incremental
C		step by assuming it is infinitesmal
C
C		Input variables:
C
C			DX1	Length of step (cm)
C			E1	Energy of neutron (eV)
C			R1	Position at which field is to be calculated (cm)
C					(center of step)
C
C		Output variable:
C
C			T1	Rotation matrix for step
C
C		Internal variables:
C
C			B	B field vector at middle of step (G)
C			BMAG	Magnitude of B field (G)
C			BUN	Unit vector in direction of B field
C			K	Counting variable
C			L	Counting variable
C			M	Counting variable
C			THETA	Precession angle in rotation matrix (rad.)
C
C		Constants:
C
C			E0	Neutron rest energy (eV)
C			I	Imaginary number, i
C			OMEGA	Precession of 1 eV neutron after travelling 1
C					cm in 1 G field
C			SI	2x2 identity matrix
C			SIGMA	3x2x2 matrix, Pauli spin vector
C			SX	\
C			SY	 }Pauli spin matrices
C			SZ	/
C
		IMPLICIT NONE
C
		INTEGER K,L,M,P
C
		REAL*8 DX1,E0,E1,OMEGA,R1(3),PI,E2,R2(3)
C
		REAL*8 B(3),B2(3),BMAG,BUN(3),THETA,BMAGAVG
C
		COMPLEX*16 I,SI(2,2),SIGMA(2,2,3),SX(2,2),SY(2,2),SZ(2,2)
C
		COMPLEX*16 T1(2,2)
		CHARACTER myfieldfile*80,fieldfile2*80
C
		EQUIVALENCE (SIGMA(1,1,1),SX),(SIGMA(1,1,2),SY),(SIGMA(1,1,3),
     >			SZ)	
c		write(6,*)myfieldfile
C
C
		DATA E0		/939.573D6/
		DATA I		/(0.0,1.0)/
		DATA OMEGA	/-1.325D-2/
		DATA SI		/(1.0,0.0),(0.0,0.0),(0.0,0.0),(1.0,0.0)/
		DATA SX		/(0.0,0.0),(1.0,0.0),(1.0,0.0),(0.0,0.0)/
		DATA SY		/(0.0,0.0),(0.0,-1.0),(0.0,1.0),(0.0,0.0)/ ! 11, 21, 12, 22
		DATA SZ		/(1.0,0.0),(0.0,0.0),(0.0,0.0),(-1.0,0.0)/
		DATA PI		/3.14159/
C
C		Calculate magnetic field vector

		CALL BFIELD(R1,B,myfieldfile)
		R2(1)=R1(1)
		R2(2)=R1(2)
		R2(3)=R1(3)-260



		if ((R2(3).gt.-195).and.(R2(3).lt.-4)) then
C			write(6,*) ' calling Bfield2 with file ', fieldfile2, R2
			CALL BFIELD2(R2,B2,fieldfile2)
		else
C			write(6,*) ' skipping ASM field R2(3)=',R2
			B2(1)=0 ; B2(2)=0; B2(3)=0
		end if
		B(1)=B2(1)+B(1)
		B(2)=B2(2)+B(2)
		B(3)=B2(3)+B(3)

		!write(6,*)'brosorB',B(2)
		!write(6,*)
C		Calculate magnitude of B
C
		BMAG=DSQRT(B(1)**2+B(2)**2+B(3)**2)
c		write(6,*)'Yack', B(1),B(2),B(3)
C
C		Calculate unit vector in direction of B in  Cartesian
C		Coordinates.
C
			IF (dabs(BMAG).lt.1.0d-7) THEN
				BUN(1)=0.0;BUN(2)=1.0;BUN(3)=0.0
				BMAG=1.0d-7
			else if (dabs(BMAG).gt.1.0d7) THEN
				write(6,*) '!! Very Larg Bmag !!, z=',R1(3)
				BUN(1)=0.0;BUN(2)=1.0;BUN(3)=0.0
c				BMAG=1.0d-7
			else if (isnan(BMAG)) THEN
			write(6,*) '!! NaN Bmag !!,',R1(1),R1(2),R1(3),
     &				B(1),B(2),B(3)
				BUN(1)=0.0;BUN(2)=1.0;BUN(3)=0.0
c				BMAG=1.0d-7
			ELSE
				BUN(1)=B(1)/BMAG
				BUN(2)=B(2)/BMAG
				BUN(3)=B(3)/BMAG
			END IF

		E2 = E1/1000.0

C		Calculate precession of spin during step
C
		THETA=OMEGA*DSQRT((1.0D0/E0)+(1.0D0/E2))*DX1*BMAG

C
C
C		Initialize T1, rotation matrix
C
		DO K=1,2
			DO L=1,2
				T1(K,L)=(0.0,0.0)
			END DO
		END DO
C
C		Calculate T1=cos(THETA/2)+isin(THETA/2)(SIGMA*BUN)
C
		DO K=1,2
			DO L=1,2
				T1(K,L)=T1(K,L)+DCOS(THETA/2.0D0)*SI(K,L)
				DO M=1,3
				   T1(K,L)=T1(K,L)+I*DSIN(THETA/2.0D0)*BUN(M)
     >						*SIGMA(K,L,M)
				END DO
			END DO
		END DO
C
		RETURN
	END
C
C
C
	SUBROUTINE PROJ (PSI1,N1,P1)
C
C		Calculates the projection of the spin wavefunction PSI1 on the
C		normal vector N1.  P1 is the result
C
C		Input variables:
C
C			N1	Normal vector
C			PSI1	Spinor wavefunction
C
C		Output variable:
C
C			P1	Projection of PSI1 on N1
C
C		Internal variables:
C
C			I	Counting variable
C			J	Counting variable
C			K	Counting variable
C			PHI	PSI1 operated on with  projection operator
C
C		Constants:
C
C			SI	2x2 identity matrix
C			SIGMA	3x2x2 matrix, Pauli spin vector
C			SX	\
C			SY	 }Pauli spin matrices
C			SZ	/
C
		IMPLICIT NONE
C
		INTEGER I,J,K
C
		REAL*8 N1(3)
C
		REAL*8 P1
C
		COMPLEX*16 SI(2,2),SIGMA(2,2,3),SX(2,2),SY(2,2),SZ(2,2)
C
		COMPLEX*16 WPHI(2),PSI1(2)
 
		EQUIVALENCE (SIGMA(1,1,1),SX),(SIGMA(1,1,2),SY),(SIGMA(1,1,3),
     >			SZ)
C
		DATA SI	/(1.0,0.0),(0.0,0.0),(0.0,0.0),(1.0,0.0)/
		DATA SX	/(0.0,0.0),(1.0,0.0),(1.0,0.0),(0.0,0.0)/
		DATA SY	/(0.0,0.0),(0.0,-1.0),(0.0,1.0),(0.0,0.0)/ ! ij=11, 21, 12, 22
		DATA SZ	/(1.0,0.0),(0.0,0.0),(0.0,0.0),(-1.0,0.0)/
C
C		Initialize PHI, result of operating on PSI1 with projection
C		operator
C
		WPHI(1)=(0.0,0.0)
		WPHI(2)=(0.0,0.0)
C

		DO I=1,2
			DO J=1,2
				WPHI(I)=WPHI(I)+0.5D0*SI(I,J)*PSI1(J)  

				DO K=1,3
					WPHI(I)=WPHI(I)+0.5D0*SIGMA(I,J,K)*N1(K)
     >						*PSI1(J)
				END DO
			END DO
		END DO



	
C
C		Initialize P1, the expectation value of the projection of PSI1
C		along N1
C
		P1=0.0D0
C
C		Calculate P1=<PSI1|PHI>
C
		DO I=1,2
			P1=P1+DBLE(CONJG(PSI1(I))*WPHI(I))
		END DO

		RETURN
	END
