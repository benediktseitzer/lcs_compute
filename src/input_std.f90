subroutine input_std
	use parameter_module, ONLY: filename,t,dt,euler,forward,x_res,y_res,lf,lmin,stepsize
	implicit none
	character             :: INTEGRATION_DIRECTION *20
	character             :: INTEGRATION_METHOD *20
	
	namelist / input_velocity / filename
	namelist / resolution /  x_res, y_res
	namelist / integration / INTEGRATION_METHOD, INTEGRATION_DIRECTION, t, dt
	namelist / variational / stepsize, lf, lmin

	open(77,file='Par',delim='APOSTROPHE')
	
	read(77,nml=INPUT_VELOCITY)
	read(77,nml=RESOLUTION)
	read(77,nml=INTEGRATION)
	read(77,nml=VARIATIONAL)
	
	close(77)

	if (trim(INTEGRATION_METHOD)=="EULER") then
		euler = .true.
	elseif (trim(INTEGRATION_METHOD)=="RK4") then
		euler = .false.
	else
		write(*,*) '  INTEGRATION METHOD IS UNDEFINED'
		stop
	end if

	if (trim(INTEGRATION_DIRECTION)=="FORWARDS") then
		forward = .true.
	elseif (trim(INTEGRATION_DIRECTION)=="BACKWARDS") then
		forward = .false.
	else
		write(*,*) ' INTEGRATION DIRECTION IS UNDEFINED'
	end if

	write(*,*) ' FILENAME: ',filename
	write(*,*) ''

end subroutine
