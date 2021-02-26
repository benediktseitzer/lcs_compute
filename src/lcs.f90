program FTLE2D
    use parameter_module
    implicit none

    write(*,*) ''
    write(*,*) ''
    write(*,*) '###############################'
    write(*,*) '        2D lcs code            '
    write(*,*) '###############################'
    write(*,*) ''
    write(*,*) ''
    
    call input_std()
    ! ######## 2D - lcs ########
    !read data
    call read_data2D()
	!FTLE and Cauchy-Green
    call gridadvect2D(euler,forward)
	!variational theory
	call variational()
	!write data
    call write_output2D()
    
    deallocate(u,v,x,y,x_input,y_input)
 	deallocate(gridY0,gridY1,gridY2)
	deallocate(gridX0,gridX1,gridX2)
end program
