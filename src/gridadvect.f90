! this subroutine creates the new grid and advects the tracers.
! after this is done, the CG-tensor can be computed
! with the CG-tensor, the FTLE-Field and also the LCS from their variational theory can be calculated... hopefully!

subroutine gridadvect2D(euler,forward)
    use parameter_module, only: nx_input,ny_input,x_input,y_input,x_res,y_res,gridX0,gridX1,gridX2,gridY0,gridY1,gridY2,DP,U,V,nx,ny,x,y,sigma,t,dt,steps,rank,alpha,lambda1,lambda2,eigv1_x,eigv1_y,eigv2_x,eigv2_y,div_eig
    implicit none
    real(DP)                :: min_x,max_x,min_y,max_y
    integer                 :: i,j,it,tr_out
    real(DP)                :: timefac,t1,t_span,tstep,max_vel,dx,dy,start,ende
    integer,allocatable     :: TracerLeftDomain(:,:),calcFTLE(:,:)
    logical,intent(in)      :: euler, forward
        
    where(isnan(U)) U = 0.0_DP
    where(isnan(V)) V = 0.0_DP

    min_x=minval(x_input); min_y=minval(y_input); max_x=maxval(x_input); max_y=maxval(y_input)
    
    ! ##############################################################    
    ! ######## create additional grid ########
    ! ##############################################################

    if (x_res /= 1.0_DP) then
        nx = ceiling(nx_input/x_res)
        dx = (max_x-min_x)/(nx-1)
        allocate(x(nx))
        do i = 1,nx
            x(i) = min_x + ((i-1)*dx)
		end do
    else
        allocate(x(nx_input))     
        x = x_input
        nx = nx_input
    end if

    if (y_res /= 1.0_DP) then
        ny = ceiling(ny_input/y_res)
        dy = (max_y-min_y)/(ny-1)
        allocate(y(ny))
        do i = 1,ny
            y(i) = min_y + ((i-1)*dy)
        end do
    else
        allocate(y(ny_input))
        y = y_input
        ny = ny_input
    end if
    
    allocate(gridX0(nx,ny),gridX1(nx,ny),gridX2(nx,ny))
    allocate(gridY0(nx,ny),gridY1(nx,ny),gridY2(nx,ny))

    do j = 1,ny
        do i = 1,nx
            gridX0(i,j) = x(i)  
 			gridY0(i,j) = y(j)
        end do
    end do    
    


    ! set initial positions of tracers
    gridX1 = gridX0
    gridY1 = gridY0
    
    print*, ''
    print*, 'resolution = ', nx,ny
    print*, 'Number of tracers =', nx*ny
    print*, ''
    
    max_vel = max(maxval(U),maxval(V))
    print*, ''
    print*, 'maximum velocity =', max_vel
    print*, ' -> dt_min:', min( (( max_x-min_x)/(nx-1))/max_vel , ( (max_y-min_y)/(ny-1))/max_vel )
    print*, ''    
    
    if (rank == 0) then
        print*, ''
        if (euler) then
            print*, 'euler time-stepping'
        else
            print*, 'rk4 time-stepping'
        end if    
         
        if (forward) then
            print*, 'forward advection'
        else
            print*, 'backward advection'
        end if
        
        print*, '    dt =',dt, '    t =', t 
        print*, ''
    end if 
 
   
    
    ! ##############################################################
    ! ######## advection of tracergrid ########
    ! ##############################################################
    allocate(sigma(nx,ny),calcFTLE(nx,ny),TracerLeftDomain(nx,ny))
    allocate(alpha(nx,ny),lambda1(nx,ny),lambda2(nx,ny),eigv1_x(nx,ny),eigv1_y(nx,ny),eigv2_x(nx,ny),eigv2_y(nx,ny),div_eig(nx,ny))
	sigma = 0.0_DP; calcFTLE = 1; TracerLeftDomain = 0; gridX2 = 0.0_DP; gridY2 = 0.0_DP
    alpha = 0.0_DP; lambda1 = 0.0_DP; lambda2 = 0.0_DP; div_eig = 0.0_DP
	
	timeFac=1.0_DP
    steps=ceiling(t/dt)
    
    do it = 1,steps

        print*, ' ## steps:',it,'/',steps,'integration time:', it*dt, 'of', steps*dt
        t1 = it*dt
        tstep = dt*timeFac
        t_span = t1
        ! decide which time integration scheme to use:
        if (euler) then
        	print*, 'nicht implementiert!'
		else
		    call rk42D(tstep)
        end if

        !is tracer still in domain? 
        do j = 1,ny
            do i = 1,nx
                if (TracerLeftDomain(i,j) == 0) then 
                    ! TRACER INSIDE: calculate trajectory. 
                    if (((gridX2(i,j)-min_x)*(gridX2(i,j)-max_x) >= 0.0_DP)  .or. ((gridY2(i,j)-min_y)*(gridY2(i,j)-max_y) >= 0.0_DP)) then
                        ! TRACER ADVECTED OUT OF DOMAIN: calculate FTLE at this point and the adjacent points.
                        TracerLeftDomain(i,j) = 1
                        
                        if (calcFTLE(i,j) == 1 ) then
                            call calculate_FTLE2D(i,j,t_span)
                            calcFTLE(i,j) = 0
                            
                            if (i > 1) then
                                if (calcFTLE(i-1,j) == 1 ) then
                                    call calculate_FTLE2D(i-1,j,t_span)
                                    !calcFTLE(i-1,j) = 0
                                end if
                            end if
                            
                            if (i < nx) then
                                if ( calcFTLE(i+1,j) == 1 ) then
                                    call calculate_FTLE2D(i+1,j,t_span)
                                    !calcFTLE(i+1,j) = 0
                                end if
                            end if

                            if (j > 1 ) then
                                if ( calcFTLE(i,j-1) == 1 ) then 
                                    call calculate_FTLE2D(i,j-1,t_span)
                                    !calcFTLE(i,j-1) = 0
                                end if
                            end if

                            if (j < ny) then
                                if ( calcFTLE(i,j+1) == 1 ) then
                                    call calculate_FTLE2D(i,j+1,t_span)
                                    !calcFTLE(i,j+1) = 0
                                end if
                            end if

                        end if   
                     
                    end if

                else
                    ! TRACER OUTSIDE: do not calculate trajectory. keep tracers position!
                    gridX2(i,j) = gridX1(i,j)
                    gridY2(i,j) = gridY1(i,j)
                end if
                
                
            end do ! i = 1,nx
        end do ! j = 1,ny
        
        ! update all tracerpositions for next integration-step
        do j = 1,ny
            do i = 1,nx
                if (TracerLeftDomain(i,j) == 0) then
                    gridX1(i,j) = gridX2(i,j)
                    gridY1(i,j) = gridY2(i,j)
                end if
            end do
        end do
        
    end do ! it= 1,steps

    ! positions after last time-step
    do j = 1,ny
        do i = 1,nx
            if (calcFTLE(i,j) == 1 ) then
                call calculate_FTLE2D(i,j,t_span)
                calcFTLE(i,j) = 0
            end if
        end do
    end do
	
	print*, ''
	print*, '############ FTLE is ready ############'
	print*, ''
	print*, 'nx =', nx, 'ny =', ny
	print*, ''	

    deallocate(TracerLeftDomain,calcFTLE)

end subroutine
