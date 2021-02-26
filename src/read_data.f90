! this subroutine reads the needed data 
! grid as two 1D-arrays and the horizontal and vertical velocity-data as 2D-arrays

subroutine read_data2D()
    use parameter_module, only: x_input, y_input, U, V, nx_input, ny_input, filename
    
    
    !open(12,file = 'dimdg')
    open(12,file = 'dim')
	!open(12,file = 'dimdat')
	read(12,*) nx_input, ny_input
    close(12)
    
    
    print*, 'number of x- and y-values = ', nx_input,ny_input

    allocate(U(nx_input,ny_input),V(nx_input,ny_input),x_input(nx_input),y_input(ny_input)) 
    
    !read x-grid
	!open(12,file = 'lon')
    open(12,file = 'xgrid')
	!open(12,file = 'xdat')
    do i = 1,nx_input
        read(12,*) x_input(i)
    end do
    close(12)    

    ! read y-grid
	!open(12,file = 'lat')
    open(12,file = 'ygrid')
	!open(12,file = 'ydat')
    do i = 1,ny_input
        read(12,*) y_input(i)
    end do
    close(12)    

    ! read horizontal velocity
	!open(12,file = 'U')
    open(12,file = 'veloX')
	!open(12,file = 'Udat')
    do i = 1,nx_input
        read(12,*) (U(i,j), j=1,ny_input)
    end do
    close(12)

    ! read vertical velocity
	!open(12,file = 'V')
    open(12,file = 'veloY')
	!open(12,file = 'Vdat')
    do i = 1,nx_input
        read(12,*) (V(i,j), j=1,ny_input)
    end do
    close(12)
    
end subroutine
