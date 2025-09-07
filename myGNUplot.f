	! example.f90
	program main
    	implicit none
    	character(len=*), parameter :: OUT_FILE = 'data.txt' ! Output file.
    	character(len=*), parameter :: PLT_FILE = 'plot.plt' ! Gnuplot file.
    	integer,          parameter :: N        = 30         ! Number of values.

    	integer :: i, fu
    	real    :: x(N), y(N)

    	x(1) = 0.0
    	y(1) = 0.0

    	do i = 2, N
        x(i) = 0.1 * i
        y(i) = x(i)**2
    	end do

    	open (action='write', file=OUT_FILE, newunit=fu, status='replace')

    	do i = 1, N
        write (fu, *) x(i), y(i)
    	end do

    	close (fu)
    	call execute_command_line('gnuplot -p ' // PLT_FILE)
	end program main
