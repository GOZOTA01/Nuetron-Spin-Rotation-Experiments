	subroutine trilinearinterp(x,y,z,x0,y0,z0,x1,y1,z1,V000,V100,V001,
     &		V101,V010,V110,V011,V111,tot)
c	Edited by Sam Infanger 6/13/18 to fix things
c	Written by Peter Yergeau, Finished:7/10/17	
c	
c	A method of interpolation in three dimensions to be used in Bfield.f
c	so that we can interpolate points that are not in our field map.
c	
c	Written using the method presented by the wikipedia page: Trilinear Interpolation
c

	implicit none
	
	real*8 x,x0,x1,xd,y,y0,y1,yd,z,z0,z1,zd 
	real*8 c00,c01,c10,c11
	real*8 V000,V100,V001,V101,V010,V110,V011,V111 
	real*8 g0,g1,tot 

c	VARIABLE DEFINITIONS
c	x is the position we are interpolating for, x0 is the nearest neighbor below and x1 is the nearest neighbor above, xd is the fractional difference, same applies to y and z
c	c00 is the x value of the field at y0,z0, c01 is the value of the x field at y0,z1, etc., V is the function of the field at the x,y,z positions for each c
c	g0 is the y value of the field at c00 and c10, g1 is the y value of the field at c01 and c11, tot is the z value of the field at g0 and g1 which is the predicted value at x,y,z
c
	xd = 0.0
	yd = 0.0
	zd = 0.0
c	write(6,*) 'x,x0,x1',x,x0,x1
c	write(6,*) 'y,y0,y1',y,y0,y1
c	write(6,*) 'z,z0,z1',z,z0,z1
c	write(6,*) (z - z0)/(z1 - z0)

	if (dabs((x1-x0)/x0).lt.1.0d-8) then
c	  write(6,*) 'x,x0,x1',x,x0,x1
	  xd=0.0d0
	else
	  xd = (x - x0)/(x1 - x0)
	end if
	if (dabs((y1-y0)/y0).lt.1.0d-8) then
c	  write(6,*) 'y,y0,y1',y,y0,y1
	  yd = 0.0d0	
	else
	  yd = (y - y0)/(y1 - y0)
	end if
	if (dabs((z1-z0)/z0).lt.1.0d-8) then
c	  write(6,*) 'z,z0,z1',z,z0,z1
	  zd=0.0d0
	else
	  zd = (z - z0)/(z1 - z0)
	end if
	if((xd.gt.1).or.isnan(xd))then	!SMI
	xd=1
	end if
	if((yd.gt.1).or.isnan(yd))then	!SMI
	yd=1
	end if
	if((zd.gt.1).or.isnan(zd))then	!SMI
	zd=1
	end if
c	write(6,*) z,y,y0,y1,yd,zd,'start'
c	write(6,*) 'Valete',v000,v010,v110,v100,v001,v011,v101,v111 
	c00 = V000*(1-xd) + V100*xd
	c01 = V001*(1-xd) + V101*xd
	c10 = V010*(1-xd) + V110*xd
	c11 = V011*(1-xd) + V111*xd

	g0 = c00*(1-yd) + c10*yd
	g1 = c01*(1-yd) + c11*yd

	tot = g0*(1-zd) + g1*zd
c	write(6,*) c00,c01,c10,c11,g0,g1,tot,'end'
c	write(6,*) 'x',x,x0,x1,xd
c	write(6,*) 'y',y,y0,y1,yd
c	write(6,*) 'z',z,z0,z1,zd
c	write(6,*) tot
	return
	end
	
	
	
