	PROGRAM Outcoil
		implicit none
		double precision E, XI, YI, ZI, XF, YF, ZF, THETAI, PHII, THETAF, PHIF, ztarg, PI
		double precision xmin, xmax,xstep,ymin,ymax,ystep,zmin,zmax,zstep,P(3),stepsize,t(5)
		INTEGER i,j,numsteps,ix,iy,a,numneutrons,datetime(8)
		character fieldfile*80, comments*50, date*12, time*12, zone*12, fieldfile2*80
		CHARACTER day*2,month*2,hour*2,minute*2,second*2
		CHARACTER filename1*80,filename2*80,filename3*80	
c		REAL*8, DIMENSION(:), ALLOCATABLE::xgrid(:,:), ygrid(:,:), zgrid(:,:)
c		REAL*8, DIMENSION(:), ALLOCATABLE::Ngridx(:,:),Ngridy(:,:),Ngridz(:,:),avgthetax(:,:),avgthetay(:,:)
		REAL*8 xgrid(200,200), ygrid(200,200), zgrid(200,200)
		REAL*8 Ngridx(200,200),Ngridy(200,200),Ngridz(200,200),avgthetax(200,200),avgthetay(200,200)
		REAL*8  R(3),xx,yy,xygrid(2500,2)	
c		REAL*8, DIMENSION(:), ALLOCATABLE:: avgthetaz(:,:)	
		REAL*8 avgthetaz(200,200)	
		character(len=128) :: fmtString
		parameter (PI=3.14159)

		THETAI= 1000.d0*PI/2
		PHII = 0

		CALL Date_and_time(date,time,zone,datetime)

		day= CHAR(INT(datetime(3)/10)+48)//CHAR(mod(datetime(3),10)+48)
    		month= CHAR(int(datetime(2)/10)+48)//CHAR(MOD(datetime(2),10)+48)
      		hour= CHAR(INT(datetime(5)/10)+48)//CHAR(mod(datetime(5),10)+48)
      		minute= CHAR(INT(datetime(6)/10)+48)//CHAR(mod(datetime(6),10)+48)
      		second= CHAR(INT(datetime(7)/10)+48)//CHAR(mod(datetime(7),10)+48)	!SMI

		filename1='thetax'//month//day//'-'//hour//minute//'.txt'

		fieldfile = 'OutCoilFullZASMNeg.dat'
c		fieldfile = ADJUSTL(fieldfile)
c		fieldfile = ADJUSTR(fieldfile)
		
		fieldfile2 = 'MagneticBox.txt'

		Open(unit=9,file=fieldfile,status='old',action='read')

c	   	write(6,'(a40,5(1x,a20),a14,3(f5.2))') ' Reading Bfield input file...',fieldfile
	    	write(3,'(a40,5(1x,a20),a14,3(f5.2))') ' Reading Bfield input file...',fieldfile


		open (action='write', file='finalthetax.txt', unit=15, status='replace')
		open (action='write', file='finalthetay.txt', unit=16, status='replace')
		open (action='write', file='finalthetaz.txt', unit=17, status='replace')
c		open (action='write', file='thetagridx.txt', unit=20, status='replace')
c		open (action='write', file='thetagridy.txt', unit=24, status='replace')
c		open (action='write', file='thetagridz.txt', unit=23, status='replace')
		open (action='write', file='xygrid.txt', unit=21, status='replace')

C
C		Read in Number of coils and Earth field
C

	    	read(9,'(a)') comments
	    	read(9,'(a)') comments
	    	read(9,'(a)') comments
	    	read(9,*) xmin,xmax,xstep,ymin,ymax,ystep,zmin,zmax,zstep
	    	read(9,'(a)') comments

		close(unit=9)
		
		ZI = zmin
		ZF= zmax
		write(6,*) ' Running from Zi= ',ZI,' to Zf= ',ZF
C		Assigning random numbers to the initial x and y values
C		call RANDOM_NUMBER(XI)
C		XI = (XI*10) - 5

C		call RANDOM_NUMBER(YI)
C		YI = (YI*10) - 5
	
		ztarg = 0
		write(6,*) "Input Energy "
		read(5,*) E   
		
		write(6,*) "0=Quadrant Simulation, 1=Numsteps Simulation, 2=Random Simulation"
		read(5,*) a
		
	  	IF (a == 0) then
		XI = 2.5
		XF = 2.5
		YI = 2.5
		YF = 2.5
	do i=1,2
		do j=1,2
		call sBrot(E,XI,YI,ZI,THETAI,PHII,XF,YF,ZF,THETAF,THETAF, ztarg,fieldfile,P,R,a,fieldfile2)
cf		write(6,*) 'sBrot ran'
		YI=(YI)*(-1)
		YF=YI*(-1)
		end do
		XI=XI*(-1)
		XF=XF*(-1)
	end do	

	ELSE IF (a == 1) then
	XI = xmin
	XF = xmin
	YI = ymin
	YF = ymin

	write(6,*) "How many Steps? "
	read(5,*) numsteps 

	filename1='thetagridx'//month//day//'-'//hour//minute//'.txt'
	filename2='thetagridy'//month//day//'-'//hour//minute//'.txt'
	filename3='thetagridz'//month//day//'-'//hour//minute//'.txt'
	open (file=filename1, unit=31, status='new')
	open (file=filename2, unit=32, status='new')
	open (file=filename3, unit=33, status='new')

	write(31,*) 'numsteps =', numsteps
	write(32,*) 'numsteps =', numsteps
	write(33,*) 'numsteps =', numsteps

	write(6,*) filename1
	write(6,*) filename2
	write(6,*) filename3

	stepsize = ((xmax-XI)/(numsteps-1))

c	allocate(xgrid(numsteps,numsteps))
c	allocate(ygrid(numsteps,numsteps))
c	allocate(zgrid(numsteps,numsteps))
	
	do iy=1,numsteps
	do ix=1,numsteps
		call sBrot(E,XI,YI,ZI,THETAI,PHII,XF,YF,ZF,THETAF,THETAF, ztarg,fieldfile,P,R,a,fieldfile2)
		
c		write(24,*) R(1), R(2), asin(2.0*P(1)-1.0)
c		write(20,*) R(1), R(2), asin(2.0*P(2)-1.0)
c		write(23,*) R(1), R(2), asin(2.0*P(3)-1.0)
		
		xgrid(ix,iy) =asin(2.0*P(1)-1.0)
		ygrid(ix,iy) =asin(2.0*P(2)-1.0)
		zgrid(ix,iy) =asin(2.0*P(3)-1.0)
		
c		write(24,*) asin(2.0*P(1)-1.0)
c		write(20,*) asin(2.0*P(2)-1.0)
c		write(23,*) asin(2.0*P(3)-1.0)

		xygrid(ix,1)=XI

		XI=XI+stepsize
		XF=XF+stepsize
	end do
	xygrid(iy,2)=YI
	write(21,*) xygrid(iy, 1),xygrid(iy, 2)

	write(31,*) ygrid(1:numsteps,iy)
	write(32,*) xgrid(1:numsteps,iy)
	write(33,*) zgrid(1:numsteps,iy)

	YI=YI+stepsize
	YF=YF+stepsize
	XI=xmin
	XF=xmin
	end do

	ELSE IF (a == 2) then


	write(6,*) "How many Neutrons? "
	read(5,*) numneutrons

	write(6,*) "How many Bins? "
	read(5,*) numsteps



	filename1='thetagridx'//month//day//'-'//hour//':'//minute//'.txt'
	filename2='thetagridy'//month//day//'-'//hour//':'//minute//'.txt'
	filename3='thetagridz'//month//day//'-'//hour//':'//minute//'.txt'
	
	open (file=filename1, unit=31, status='new')
	open (file=filename2, unit=32, status='new')
	open (file=filename3, unit=33, status='new')

	write(31,*) 'numsteps =', numsteps, ',nuetrons=', numneutrons 
	write(32,*) 'numsteps =', numsteps, ',nuetrons=', numneutrons 
	write(33,*) 'numsteps =', numsteps, ',nuetrons=', numneutrons 


	write(6,*) filename1
	write(6,*) filename2
	write(6,*) filename3
	write(fmtString,*) numsteps
c	write(6,*) numsteps, fmtString
	fmtString = "("//trim(adjustl(fmtString))//"(1pe11.4,1x))"
	
c	allocate(xgrid(numsteps,numsteps))
c	allocate(ygrid(numsteps,numsteps))
c	allocate(zgrid(numsteps,numsteps))
c
c	allocate(Ngridx(numsteps,numsteps))
c	allocate(Ngridy(numsteps,numsteps))
c	allocate(Ngridz(numsteps,numsteps))
c	
c	allocate(avgthetax(numsteps,numsteps))
c	allocate(avgthetay(numsteps,numsteps))
c	allocate(avgthetaz(numsteps,numsteps))

	do iy=1,numsteps
		do ix=1,numsteps
			xgrid(ix,iy) = 0.0d0
			ygrid(ix,iy) = 0.0d0
			zgrid(ix,iy) = 0.0d0
			Ngridx(ix,iy) = 0
			Ngridy(ix,iy) = 0
			Ngridz(ix,iy) = 0
		end do
	end do
	  
	do i=1,numneutrons
	
		call RANDOM_NUMBER(t)
c		write(6,*)t
		XI = t(1)*(xmax-xmin)+xmin
		YI = t(2)*(ymax-ymin)+ymin
		XF = t(3)*(xmax-xmin)+xmin
		YF = t(4)*(ymax-ymin)+ymin
		call sBrot(E,XI,YI,ZI,THETAI,PHII,XF,YF,ZF,THETAF,THETAF, ztarg,fieldfile,P,R,a,fieldfile2)
c		write(6,*) i,XI,XF,YI,YF
		ix = int(((XF-xmin)/(xmax-xmin))*numsteps)  
		iy = int(((YF-ymin)/(ymax-ymin))*numsteps)  
		xgrid(ix,iy) = xgrid(ix,iy)+asin(2.0*P(1)-1.0)
		ygrid(ix,iy) = ygrid(ix,iy)+asin(2.0*P(2)-1.0)
		zgrid(ix,iy) = zgrid(ix,iy)+asin(2.0*P(3)-1.0)

		Ngridx(ix,iy) = Ngridx(ix,iy) + 1 
		Ngridy(ix,iy) = Ngridy(ix,iy) + 1 
		Ngridz(ix,iy) = Ngridz(ix,iy) + 1 

	end do

c	write(6,*) ' Done with neutrons '

	do iy=1,numsteps
		do ix=1,numsteps
			if (dabs(Ngridy(ix,iy)).gt.1.0d-10) then
				avgthetay(ix,iy) = ygrid(ix,iy)/ Ngridy(ix,iy) 
			else
				avgthetay(ix,iy) = 0
			END if

			if (dabs(Ngridx(ix,iy)).gt.1.0d-10) then
				avgthetax(ix,iy) = xgrid(ix,iy)/ Ngridx(ix,iy) 
			else
				avgthetax(ix,iy) = 0
			END if

			if (dabs(Ngridz(ix,iy)).gt.1.0d-10) then
				avgthetaz(ix,iy) = zgrid(ix,iy)/ Ngridz(ix,iy) 
			else
				avgthetaz(ix,iy) = 0
			END if

c		if (dabs(Ngridx(ix,iy)).gt.1.0d-10) avgthetax(ix,iy) = xgrid(ix,iy)/ Ngridx(ix,iy)
c		if (dabs(Ngridz(ix,iy)).gt.1.0d-10) avgthetaz(ix,iy) = zgrid(ix,iy)/ Ngridz(ix,iy)

		end do

		write(31,fmtString) avgthetay(1:numsteps,iy)
		write(32,fmtString) avgthetax(1:numsteps,iy)
		write(33,fmtString) avgthetaz(1:numsteps,iy)
	end do

		do i=1,numsteps
			xygrid(i,1)=(i-1)*(xmax-xmin)/(numsteps-1) + xmin			
			xygrid(i,2)=(i-1)*(ymax-ymin)/(numsteps-1) + ymin
			write(21, *) xygrid(i, 1),xygrid(i, 2)
		end do

	END IF



C	write(15,*) R(1), R(2), asin(2.0*P(1)-1.0)
C	write(16,*) R(1), R(2), asin(2.0*P(2)-1.0)
C	write(17,*) R(1), R(2), asin(2.0*P(3)-1.0)

c	do i=1,100
c		call RANDOM_NUMBER (XI)
c		XI = (XI*10) - 5 
c		XF=XI
c		call RANDOM_NUMBER(YI)
c		YI = (YI*10) - 5
c		YF=YI
c		call sBrot(E,XI,YI,ZI,THETAI,PHII,XF,YF,ZF,THETAF,THETAF, ztarg,fieldfile,fieldfile2)
c	end do
	
	close(15)
	close(16)
	close(17)
	close(9)
c	deallocate(xgrid,ygrid,zgrid,Ngridx,Ngridy,Ngridz,avgthetax,avgthetay,avgthetaz)

C	call execute_command_line('gnuplot -p ' // 'scatterplotx.plt')
c	call execute_command_line('gnuplot -p ' // 'scatterploty.plt')
c	call execute_command_line('gnuplot -p ' // 'scatterplotz.plt')
	IF ((a == 1).or.(a == 2)) then
	call execute_command_line('python3 ' // 'heatmap.py')
	END IF

		end


