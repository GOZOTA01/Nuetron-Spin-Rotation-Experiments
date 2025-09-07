	SUBROUTINE BFIELD2(R2,TOTB1,fieldfile)
	
C               edited May 2018 by Sam Infanger to optimize speed
c		edited June 2017 Peter Yergeau to do tri-linear interpolation of field values
c		edited April 2016 BEC
c			edited to fake NSR input/output coil transitions for calculated field and
c                       to allow B-field to be read from a file instead of calculated
C		Written 5/25/94 by Bret Crawford, TUNL
C
C		Calculates the magnetic field due to gaps in solenoids and
C		and the Earth's field in all three dimensions.   For
C		now the calculatiion is a small off axis approximation
C		obtained by assuming the axial field is independent of radial
C		distance and the radial field is then found from the zero
C		divergence of the magnetic field.  The axis of the solenoid
C		is along the z direction.  The first
C		solenoid ends at z=0 and the second starts at z=d...
C		The input file is bFIELD.IN.
C
C		Input variables:
C
C			R2	Position at which field is calculated (cm)

C
C			B1	Three D magnetic field due to one gap (Gauss)
C
C		Output variable:
C
C			TOTB1	Calculated magnetic field (G)
C
C
C		Constants:
C
C			PI	pi
C			mu	mu zero - free space permeability (10^4 10^2 H/m
C					i.e., Gauss centimeters per amp)
C
		IMPLICIT NONE

		include 'const.inc'

		integer*4 	i,j,ntot(5),fi
		REAL*8		R2(3),NI(40),a(40),B1b(3),R3(3)
		REAL*8		B1(3),mu,TOTB1(3),length(40),wid,zshift

		real*8		xx2(5,500000),yy2(5,500000),zz2(5,500000),
     >			Bx2(5,500000),By2(5,500000),Bz2(5,500000),Binterp(3)
		integer 	nin000,nin100,nin001,nin101,nin010,nin110
		integer 	nin011,nin111
		double precision xmin2(5),xmax2(5),xstep2(5),ymin2(5),ymax2(5),ystep2(5),zmin2(5),zmax2(5),zstep2(5)
		double precision rindx,rindy,rindz,nx2(5),ny2(5),nz2(5),icz0,icz1,pciz0,pciz1,pcoz0,pcoz1,ocz0,ocz1
		double precision mypibfield,mypiscale
C
		Character comments*80,fieldfile*80
		logical fileflg2,found,readflg,readflg2,ocflg
		 
c		common/bfld/readflg,db,zgrad2,zgrad3
		common/bfld/readflg,readflg2,ocflg

		SAVE xx2,yy2,zz2,Bx2,By2,Bz2
		SAVE xmin2,xmax2,ymin2,ymax2,zmin2,zmax2,xstep2,ystep2,zstep2
		SAVE nx2,ny2,nz2,fileflg2
		SAVE /bfld/
		
	        found=.false.
c		write(6,*)'bfield'
C		Calculate pi
C
cbec		PI=DACOS(-1.0D0)
C
C		Calculate mu -- mu zero in Guass centimeters per amp
C
		mu= 4.0D-1*PI
C 
C		Initialize B-field components
C
	Do j=1,3
		TOTB1(j) = 0.0
		Binterp(j) = 0.0
	End do		
C

cccc   start and end z-positions for field maps
	

C	write(6,*) 'at bfield stuff'
c	fieldfile='VertPlateBfield2.txt'

	if (.not.(readflg2)) then
	Open(unit=89,file=fieldfile,status='old',action='read')
	    readflg2 = .true.
c
C	    write(6,'(a40,a30)') ' Reading Bfield2 input file...',fieldfile
cc	    write(6,'(a40,5(1x,a20),a14,3(f5.2))') ' Reading Bfield input file...',fieldfile
cc	    write(3,'(a40,5(1x,a20),a14,3(f5.2))') ' Reading Bfield input file...',fieldfile
C
C		Read in Number of coils and Earth field
C

	    i=1  ! just reading one file
	    read(89,'(a)') comments
	    read(89,'(a)') comments
	    read(89,'(a)') comments
	    read(89,*) xmin2(i),xmax2(i),xstep2(i),ymin2(i),ymax2(i),ystep2(i),zmin2(i),zmax2(i),zstep2(i)
	    read(89,'(a)') comments
	
C		write(6,*) 'I got to the do loop'
c	      if (i.eq.4) write(98,'(6(2x,f11.4))') xmin(i),xmax(i),xstep(i),ymin(i),ymax(i),ystep(i),zmin(i),zmax(i),zstep(i)
c	      if (i.eq.5) write(98,'(6(2x,f11.4))') xmin(i),xmax(i),xstep(i),ymin(i),ymax(i),ystep(i),zmin(i),zmax(i),zstep(i)

ccccccccccc   Read in data arrays, get total number of lines, find number of z's and number in the xy groups
	    do j=1,100000    ! will jump out of loop when at end of file
	      read(89,*,end=100) xx2(i,j),yy2(i,j),zz2(i,j),Bx2(i,j),By2(i,j),Bz2(i,j)
c	      if ((i.eq.4).or.(i.eq.5)) write(98,'(6(2x,f11.4))') xx2(i,j),yy(i,j),zz(i,j),Bx(i,j),By(i,j),Bz(i,j)
	    end do
	    
100	continue
	   ntot(i)=j-1


	   backspace(unit=89)
	   CLOSE(unit=89)


C	   write(6,*) ' Done reading Field Map'
c	   write(6,*) ; write(6,*) ' check on stepping '

	   nx2(i) = int((xmax2(i)-xmin2(i))/xstep2(i))+1
	   ny2(i) = int((ymax2(i)-ymin2(i))/ystep2(i))+1
	   nz2(i) = int((zmax2(i)-zmin2(i))/zstep2(i))+1

c	   write(6,*) 'nx ',ny(i)
c	   write(6,*) 'ny ',nx(i)
c	   write(6,*) 'nz ',nz(i)
	

        end if   !  end of file reads and storing data in arrays for future neutrons

	fi = 1  
      
cccccccccccccccccccccccccccccccc   Find locations of 8 grid points surrounging point of interest  cccccccccccccccccccccc

cc


c 	write(6,'(8(2x,f9.2))') icz0,icz1,pciz0,pciz1,pcoz0,pcoz1,ocz0,ocz1    
C	write(6,*)  ' Bfield.f  rvalue ',r2(1),r2(2),r2(3)


	R3(1)=R2(1);R3(2)=R2(2);R3(3)=R2(3) !SMI    creating copy r-value
c	write(6,*)  r3(1),r3(2),r3(3)

c	write(6,*) ' nx,ny,nz ',nx,ny,nz
c	write(6,*) ' ntot ',ntot,' nx*ny*nz ',nx*ny*nz
	

	fi = 1 
	rindx = int((R3(1)-xmin2(fi))/xstep2(fi))+1  !   int takes floor so this should be indx for <= value
	
	rindy = int((R3(2)-ymin2(fi))/ystep2(fi))+1
c	write(6,*) rindx
	rindz = int((R3(3)-zshift-zmin2(fi))/zstep2(fi))+1
c	write(6,*) rindx

	nin000 = rindx + (rindy-1)*nx2(1) + (rindz-1)*nx2(1)*ny2(1) 	! data array pos for xyz point < field point
	nin100 = rindx+1 + (rindy-1)*nx2(1) + (rindz-1)*nx2(1)*ny2(1) ! data array pos for grid point one > in x
	nin010 = rindx + (rindy)*nx2(1) + (rindz-1)*nx2(1)*ny2(1)  	! data array pos for grid point one > in y
	nin001 = rindx + (rindy-1)*nx2(1) + (rindz)*nx2(1)*ny2(1)  	! data array pos for grid point one > in z
	nin110 = rindx+1 + (rindy)*nx2(1) + (rindz-1)*nx2(1)*ny2(1)   ! data array pos for grid point one > in x and y
	nin101 = rindx+1 + (rindy-1)*nx2(1) + (rindz)*nx2(1)*ny2(1)   ! data array pos for grid point one > in x and z
	nin011 = rindx + (rindy)*nx2(1) + (rindz)*nx2(1)*ny2(1)  	 ! data array pos for grid point one > in y and z
	nin111 = rindx+1 + (rindy)*nx2(1) + (rindz)*nx2(1)*ny2(1) 	 ! data array pos for grid point one > in x,y, and z

c	if(R3(3).gt.70.0) write(6,*) '   R point ',R3(1),R3(2),R3(3),' file# ',fi
c	if(R3(3).gt.70.0) write(6,*) ' zshift',zshift,' xmin ',xmin(fi)
c	if(R3(3).gt.70.0) write(6,*) ' rindex x,y,z ',rindx,rindy,rindz
c	if(R3(3).gt.70.0) write(6,*) ' (rindz-1)*nx*ny ',(rindz-1)*nx*ny 
c	if(R3(3).gt.70.0) write(6,*) ' (rindy-1)*nx  ',(rindy-1)*nx , '   rindx ', rindx
c	if(R3(3).gt.70.0) write(6,*) '   nin000  ',nin000, ' nin100 ',nin100,' nin010 ',nin010, ' nin001 ',nin001
c	if(R3(3).gt.70.0) write(6,*) '   nin110  ',nin110, ' nin101 ',nin101,' nin011 ',nin011, ' nin111 ',nin111
c	if(R3(3).gt.70.0) write(6,*) '  Bz   ', Bz(fi,nin000),Bz(fi,nin100),Bz(fi,nin001),
c     &  Bz(fi,nin101),Bz(fi,nin010),Bz(fi,nin110),Bz(fi,nin011),
c     &	Bz(fi,nin111)
c        if(R3(3).gt.70.0) write(6,*)
	if(nin100.eq.0)then
	nin100=nin000
	end if


	if(nin010.eq.0)then
	nin010=nin100
	end if


	if(nin001.eq.0)then
	nin001=nin000
	end if
	

	if(nin110.eq.0)then
	nin110=nin010
	if(nin110.eq.0.and.nin010.eq.0)then
	nin110=nin100
	end if
	end if
	
	if(nin101.eq.0)then
	nin101=nin001
	if(nin101.eq.0.and.nin001.eq.0) then
	nin101=nin100
	end if
	end if

	if(nin011.eq.0)then
	nin011=nin001
	if(nin011.eq.0.and.nin001.eq.0)then
	nin011=nin010
	end if
	end if

	if(nin111.eq.0) then
	nin111=nin101
	if(nin111.eq.0.and.nin101.eq.0)then
	nin111=nin011
	if(nin111.eq.0.and.nin101.eq.0.and.nin011.eq.0)then
	nin111=nin110
	end if
	end if
	end if
	
cccccc  Find interpolated Field values  ccccccccc

c	write(6,*) R2(1),R2(2),R2(3),bx(nin000),by(nin000),bz(nin000)

	call trilinearinterp(R3(1),R3(2),R3(3),xx2(fi,nin000),
     & 	yy2(fi,nin000),zz2(fi,nin000),xx2(fi,nin111),yy2(fi,nin111),
     &  zz2(fi,nin111),Bx2(fi,nin000),Bx2(fi,nin100),Bx2(fi,nin001),
     &  Bx2(fi,nin101),Bx2(fi,nin010),Bx2(fi,nin110),Bx2(fi,nin011),
     &	Bx2(fi,nin111),Binterp(1))
c	write(6,*) nin000,bx(fi,nin000),nin111,bx(fi,nin111)
		B1(1) = Binterp(1)
c	write(6,*)R2(3),nin000,nin111
	call trilinearinterp(R3(1),R3(2),R3(3),xx2(fi,nin000),
     & 	yy2(fi,nin000),zz2(fi,nin000),xx2(fi,nin111),yy2(fi,nin111),
     &  zz2(fi,nin111),By2(fi,nin000),By2(fi,nin100),By2(fi,nin001),
     &  By2(fi,nin101),By2(fi,nin010),By2(fi,nin110),By2(fi,nin011),
     & 	By2(fi,nin111),Binterp(2))

		B1(2) = Binterp(2)
!		write(6,*)'fieldsor1',By(fi,nin000),By(fi,nin100),By(fi,nin001),
!     &  By(fi,nin101),By(fi,nin010),By(fi,nin110),By(fi,nin011),
!     & 	By(fi,nin111),binterp(2),R3(3)
!		write(6,*)nin000,nin100,nin001,
!     &  nin101,nin010,nin110,nin011,
!     & 	nin111

	call trilinearinterp(R3(1),R3(2),R3(3),xx2(fi,nin000),
     & 	yy2(fi,nin000),zz2(fi,nin000),xx2(fi,nin111),yy2(fi,nin111),
     &  zz2(fi,nin111),Bz2(fi,nin000),Bz2(fi,nin100),Bz2(fi,nin001),
     &  Bz2(fi,nin101),Bz2(fi,nin010),Bz2(fi,nin110),Bz2(fi,nin011),
     &	Bz2(fi,nin111),Binterp(3))

		B1(3) = Binterp(3)	

			found=.true.
c		write(6,*) R2(3),B1(1),B1(2),'end sam'

127	format(8x,a20,1x,i8,1xf8.2,1x,f8.2,1x,f8.2)

111	continue
	if (.not.(found)) write(6,'(a40,1x,f8.2,1x,f8.2,1x,f8.2)') 
     >      ' Could not find B-values at R=',R2(1),R2(2),R2(3)


	  TOTB1(1) = B1(1) 
	  TOTB1(2) = B1(2) 
	  TOTB1(3) = B1(3)
c
c	if ((R3(3).ge.icz0).and.(R3(3).lt.icz1)) then
c	  TOTB1(1) = icscale*B1(1) + ambBx + gradBxx*R2(1)+gradBxy*R2(2)  ! so gradBxx is in mGauss/cm
c	  TOTB1(2) = icscale*B1(2) + ambBy
c	  TOTB1(3) = icscale*B1(3) + ambBz + gradBzx*R2(1)+gradBzy*R2(2)
c	elseif ((R3(3).ge.ocz0).and.(R3(3).lt.ocz1)) then
c	  TOTB1(1) = ocscale*B1(1) + ambBx + gradBxx*R2(1)+gradBxy*R2(2)
c	  TOTB1(2) = ocscale*B1(2) + ambBy
c	  TOTB1(3) = ocscale*B1(3) + ambBz + gradBzx*R2(1)+gradBzy*R2(2)
c	elseif (((R3(3).ge.pcoz0).and.(R3(3).lt.pcoz1)).or.((R3(3).ge.pciz0).and.(R3(3).lt.pciz1))) then
c	  TOTB1(1) = mypiscale*B1(1) + ambBx + gradBxx*R2(1)+gradBxy*R2(2)
c	  TOTB1(2) = mypiscale*B1(2) + ambBy
c	  TOTB1(3) = mypiscale*B1(3) + ambBz + gradBzx*R2(1)+gradBzy*R2(2)
c	else
c	  TOTB1(1) = B1(1) + ambBx + gradBxx*R2(1)+gradBxy*R2(2)
c	  TOTB1(2) = B1(2) + ambBy
c	  TOTB1(3) = B1(3) + ambBz + gradBzx*R2(1)+gradBzy*R2(2)
c	end if	

c	write(6,*) 'end of bfield.f fine',R2(3)
c	write(6,*)'Tot', TotB1

		RETURN
	END
