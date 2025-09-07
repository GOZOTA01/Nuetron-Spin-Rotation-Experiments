FC	      = gfortran
GC            = gfortran 

FFLAGS	      = -m64 -ffixed-line-length-132

LDFLAGS	      = 

MAIN	      = OC

LIBDIR			=/usr/lib/

OBJS	      = Outcoil.o sBrot.o Bfield.o trilinearinterp.o Bfield2.o

none:	        main

all:	        OC

sim_ray:        OC

main: $(OBJS)
	$(FC) -o $(MAIN)  $(OBJS) $(LIBS) $(FFLAGS) $(LDFLAGS)

clean:
	rm -f *.o OC