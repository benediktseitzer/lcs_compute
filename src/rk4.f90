! rk4 for advecting particles in velocity-field
subroutine rk42D(tstep_sec)
    use parameter_module, ONLY: GridX1,GridY1,GridX2,GridY2,DP,pi,forward,x_input,y_input,U,V,nx,ny,nx_input,ny_input
    use interp_module, ONLY: interp_bilinear,polint2d,polint1d,interp_poly,spline,splint,splie2,splin2,interp_spline
    use omp_lib
    implicit none
    real(DP),intent(in) :: tstep_sec
    real(DP)        :: gridX1_new(nx,ny), gridY1_new(nx,ny)
    real(DP)        :: rkU1(nx,ny),rkU2(nx,ny),rkU3(nx,ny),rkU4(nx,ny),start,ende
    real(DP)        :: rkV1(nx,ny),rkV2(nx,ny),rkV3(nx,ny),rkV4(nx,ny)
    integer         :: i,j
    character*6     :: interpolation

    ! type 'spline', 'polyno', 'bicubi' or 'biline' for different interpolations.
    interpolation = 'biline'
	
	
    rku1=0.0_DP; rku2=0.0_DP; rku3=0.0_DP; rku4=0.0_DP
    rkv1=0.0_DP; rkv2=0.0_DP; rkv3=0.0_DP; rkv4=0.0_DP


    if (forward) then
        ! forward integration
        !***********************************
        ! 1. step
        !***********************************
 		
		start = omp_get_wtime()
		   	
		!$omp parallel do
        do j=1,ny
            do i=1,nx
                if (interpolation == 'spline') then
					rkU1(i,j) = interp_spline(x_input,y_input,nx_input,ny_input,gridX1(i,j),gridY1(i,j),U)
					rkV1(i,j) = interp_spline(x_input,y_input,nx_input,ny_input,gridX1(i,j),gridY1(i,j),V)
				elseif(interpolation == 'polyno') then
					rkU1(i,j) = interp_poly(x_input,y_input,nx_input,ny_input,gridX1(i,j),gridY1(i,j),U)
					rkV1(i,j) = interp_poly(x_input,y_input,nx_input,ny_input,gridX1(i,j),gridY1(i,j),V)        
				elseif(interpolation == 'biline') then
					rkU1(i,j) = interp_bilinear(x_input,y_input,nx_input,ny_input,gridX1(i,j),gridY1(i,j),U)
					rkV1(i,j) = interp_bilinear(x_input,y_input,nx_input,ny_input,gridX1(i,j),gridY1(i,j),V)
				end if
            end do
        end do
    	!$omp end parallel do
        where (isnan(rku1)) rku1=0.0_DP
        where (isnan(rkv1)) rkv1=0.0_DP

        !new coordinates
		!$omp parallel do
        do j=1,ny
            do i=1,nx
                gridX1_new(i,j)=gridX1(i,j)+((tstep_sec*rkU1(i,j))/2.0_DP)
                gridY1_new(i,j)=gridY1(i,j)+((tstep_sec*rkV1(i,j))/2.0_DP)
            end do
        end do
    	!$omp end parallel do
        print*, 'first RK4-step computed'

        !***********************************
        ! 2.step
        !***********************************
		!$omp parallel do
        do j=1,ny
            do i=1,nx
                if (interpolation == 'spline') then
                    rkU2(i,j) = interp_spline(x_input,y_input,nx_input,ny_input,gridX1_new(i,j),gridY1_new(i,j),U)
                    rkV2(i,j) = interp_spline(x_input,y_input,nx_input,ny_input,gridX1_new(i,j),gridY1_new(i,j),V)
                elseif(interpolation == 'polyno') then
                    rkU2(i,j) = interp_poly(x_input,y_input,nx_input,ny_input,gridX1_new(i,j),gridY1_new(i,j),U)
                    rkV2(i,j) = interp_poly(x_input,y_input,nx_input,ny_input,gridX1_new(i,j),gridY1_new(i,j),V)
                elseif(interpolation == 'biline') then
                    rkU2(i,j) = interp_bilinear(x_input,y_input,nx_input,ny_input,gridX1_new(i,j),gridY1_new(i,j),U)
                    rkV2(i,j) = interp_bilinear(x_input,y_input,nx_input,ny_input,gridX1_new(i,j),gridY1_new(i,j),V)
                end if
            end do
        end do    
 		!$omp end parallel do
        where (isnan(rku2)) rku2 = 0.0_DP
        where (isnan(rkv2)) rkv2 = 0.0_DP

        !new coordinates
		!$omp parallel do
        do j=1,ny
            do i=1,nx
                gridX1_new(i,j)=gridX1(i,j)+((tstep_sec*rkU2(i,j))/2.0_DP)
                gridY1_new(i,j)=gridY1(i,j)+((tstep_sec*rkV2(i,j))/2.0_DP)
            enddo
        enddo
    	!$omp end parallel do
        print*, 'second RK4-step computed'
    
        !***********************************
        ! 3. step
        !***********************************
		!$omp parallel do
        do j=1,ny
            do i=1,nx
                if (interpolation == 'spline') then
                    rkU3(i,j) = interp_spline(x_input,y_input,nx_input,ny_input,gridX1_new(i,j),gridY1_new(i,j),U)
                    rkV3(i,j) = interp_spline(x_input,y_input,nx_input,ny_input,gridX1_new(i,j),gridY1_new(i,j),V)
                elseif (interpolation == 'polyno') then
                    rkU3(i,j) = interp_poly(x_input,y_input,nx_input,ny_input,gridX1_new(i,j),gridY1_new(i,j),U)
                    rkV3(i,j) = interp_poly(x_input,y_input,nx_input,ny_input,gridX1_new(i,j),gridY1_new(i,j),V)
                elseif (interpolation == 'biline') then
                    rkU3(i,j) = interp_bilinear(x_input,y_input,nx_input,ny_input,gridX1_new(i,j),gridY1_new(i,j),U)
                    rkV3(i,j) = interp_bilinear(x_input,y_input,nx_input,ny_input,gridX1_new(i,j),gridY1_new(i,j),V)
                end if
            enddo
        enddo    
  		!$omp end parallel do
        where (isnan(rku3)) rku3 = 0.0_DP
        where (isnan(rkv3)) rkv3 = 0.0_DP

        !new coordinates
		!$omp parallel do
        do j=1,ny
            do i=1,nx
                gridX1_new(i,j)=gridX1(i,j)+(tstep_sec*rkU3(i,j))
                gridY1_new(i,j)=gridY1(i,j)+(tstep_sec*rkV3(i,j))
            enddo
        enddo
    	!$omp end parallel do
        print*, 'third RK4-step computed'
    
        !***********************************
        ! 4. step
        !***********************************
		!$omp parallel do
        do j=1,ny
            do i=1,nx
                if (interpolation == 'spline') then
                    rkU4(i,j) = interp_spline(x_input,y_input,nx_input,ny_input,gridX1_new(i,j),gridY1_new(i,j),U)
                    rkV4(i,j) = interp_spline(x_input,y_input,nx_input,ny_input,gridX1_new(i,j),gridY1_new(i,j),V)
                elseif (interpolation == 'polyno') then
                    rkU4(i,j) = interp_poly(x_input,y_input,nx_input,ny_input,gridX1_new(i,j),gridY1_new(i,j),U)
                    rkV4(i,j) = interp_poly(x_input,y_input,nx_input,ny_input,gridX1_new(i,j),gridY1_new(i,j),V)
                elseif (interpolation == 'biline') then
                    rkU4(i,j) = interp_bilinear(x_input,y_input,nx_input,ny_input,gridX1_new(i,j),gridY1_new(i,j),U)
                    rkV4(i,j) = interp_bilinear(x_input,y_input,nx_input,ny_input,gridX1_new(i,j),gridY1_new(i,j),V)
                end if
            enddo
        enddo
    	!$omp end parallel do
        where (isnan(rku4)) rku4 = 0.0_DP
        where (isnan(rkv4)) rkv4 = 0.0_DP

        print*, 'fourth RK4-step computed'
    
        ! advect tracer
		!$omp parallel do
        do j = 1,ny
            do i = 1,nx
                gridX2(i,j)= gridX1(i,j) + (tstep_sec*((rkU1(i,j)+2.0_DP*rkU2(i,j)+ 2.0_DP* rkU3(i,j) +rkU4(i,j))/6.0_DP))
                gridY2(i,j)= gridY1(i,j) + (tstep_sec*((rkV1(i,j)+2.0_DP*rkV2(i,j)+ 2.0_DP* rkV3(i,j) +rkV4(i,j))/6.0_DP))
            end do    
        end do
       	!$omp end parallel do
        
		ende = omp_get_wtime()
		print*, '					cpu-time =		', ende-start

		print*, ''
        print*, 'tracer advected'
        print*, ''
        
        
    else 
        
        !backwards integration
        !***********************************
        ! 1. step
        !***********************************
        
		
        start = omp_get_wtime()
 		       
        !$omp parallel do 
        do j=1,ny
            do i=1,nx
                if (interpolation == 'spline') then
                    rkU1(i,j) = interp_spline(x_input,y_input,nx_input,ny_input,gridX1(i,j),gridY1(i,j),U)
                    rkV1(i,j) = interp_spline(x_input,y_input,nx_input,ny_input,gridX1(i,j),gridY1(i,j),V)
                elseif(interpolation == 'polyno') then
                    rkU1(i,j) = interp_poly(x_input,y_input,nx_input,ny_input,gridX1(i,j),gridY1(i,j),U)
                    rkV1(i,j) = interp_poly(x_input,y_input,nx_input,ny_input,gridX1(i,j),gridY1(i,j),V)        
                elseif(interpolation == 'biline') then
                    rkU1(i,j) = interp_bilinear(x_input,y_input,nx_input,ny_input,gridX1(i,j),gridY1(i,j),U)
                    rkV1(i,j) = interp_bilinear(x_input,y_input,nx_input,ny_input,gridX1(i,j),gridY1(i,j),V)
                end if
            enddo
        enddo    
        !$omp end parallel do

        where (isnan(rku1)) rku1=0.0_DP
        where (isnan(rkv1)) rkv1=0.0_DP
    
        !new coordinates

        !$omp parallel do
        do j=1,ny
            do i=1,nx
                gridX1_new(i,j)=gridX1(i,j)-((tstep_sec*rkU1(i,j))/2.0_DP)
                gridY1_new(i,j)=gridY1(i,j)-((tstep_sec*rkV1(i,j))/2.0_DP)
            end do
        end do
        !$omp end parallel do
        
        print*, 'first RK4-step computed'
    
        !***********************************
        ! 2.step
        !***********************************

        !$omp parallel do 
        do j=1,ny
            do i=1,nx
                if (interpolation == 'spline') then
                    rkU2(i,j) = interp_spline(x_input,y_input,nx_input,ny_input,gridX1_new(i,j),gridY1_new(i,j),U)
                    rkV2(i,j) = interp_spline(x_input,y_input,nx_input,ny_input,gridX1_new(i,j),gridY1_new(i,j),V)
                elseif(interpolation == 'polyno') then
                    rkU2(i,j) = interp_poly(x_input,y_input,nx_input,ny_input,gridX1_new(i,j),gridY1_new(i,j),U)
                    rkV2(i,j) = interp_poly(x_input,y_input,nx_input,ny_input,gridX1_new(i,j),gridY1_new(i,j),V)
                elseif(interpolation == 'biline') then
                    rkU2(i,j) = interp_bilinear(x_input,y_input,nx_input,ny_input,gridX1_new(i,j),gridY1_new(i,j),U)
                    rkV2(i,j) = interp_bilinear(x_input,y_input,nx_input,ny_input,gridX1_new(i,j),gridY1_new(i,j),V)
                end if
            enddo
        enddo    
        !$omp end parallel do
        
        where (isnan(rku2)) rku2 = 0.0_DP
        where (isnan(rkv2)) rkv2 = 0.0_DP

        !new coordinates
        !$omp parallel do
        do j=1,ny
            do i=1,nx
                gridX1_new(i,j)=gridX1(i,j)-((tstep_sec*rkU2(i,j))/2.0_DP)
                gridY1_new(i,j)=gridY1(i,j)-((tstep_sec*rkV2(i,j))/2.0_DP)
            enddo
        enddo
        !$omp end parallel do
        
        print*, 'second RK4-step computed'
    
        !***********************************
        ! 3. step
        !***********************************

        !$omp parallel do
        do j=1,ny
            do i=1,nx
                if (interpolation == 'spline') then
                    rkU3(i,j) = interp_spline(x_input,y_input,nx_input,ny_input,gridX1_new(i,j),gridY1_new(i,j),U)
                    rkV3(i,j) = interp_spline(x_input,y_input,nx_input,ny_input,gridX1_new(i,j),gridY1_new(i,j),V)
                elseif (interpolation == 'polyno') then
                    rkU3(i,j) = interp_poly(x_input,y_input,nx_input,ny_input,gridX1_new(i,j),gridY1_new(i,j),U)
                    rkV3(i,j) = interp_poly(x_input,y_input,nx_input,ny_input,gridX1_new(i,j),gridY1_new(i,j),V)
                elseif (interpolation == 'biline') then
                    rkU3(i,j) = interp_bilinear(x_input,y_input,nx_input,ny_input,gridX1_new(i,j),gridY1_new(i,j),U)
                    rkV3(i,j) = interp_bilinear(x_input,y_input,nx_input,ny_input,gridX1_new(i,j),gridY1_new(i,j),V)
                end if
            enddo
        enddo    
        !$omp end parallel do

        where (isnan(rku3)) rku3 = 0.0_DP
        where (isnan(rkv3)) rkv3 = 0.0_DP
        
        !new coordinates
        !$omp parallel do
        do j=1,ny
            do i=1,nx
                gridX1_new(i,j)=gridX1(i,j)-(tstep_sec*rkU3(i,j))
                gridY1_new(i,j)=gridY1(i,j)-(tstep_sec*rkV3(i,j))
            enddo
        enddo
        !$omp end parallel do
    
        print*, 'third RK4-step computed'
    
        !***********************************
        ! 4. step
        !***********************************

        !$omp parallel do 
        do j=1,ny
            do i=1,nx
                if (interpolation == 'spline') then
                    rkU4(i,j) = interp_spline(x_input,y_input,nx_input,ny_input,gridX1_new(i,j),gridY1_new(i,j),U)
                    rkV4(i,j) = interp_spline(x_input,y_input,nx_input,ny_input,gridX1_new(i,j),gridY1_new(i,j),V)
                elseif (interpolation == 'polyno') then
                    rkU4(i,j) = interp_poly(x_input,y_input,nx_input,ny_input,gridX1_new(i,j),gridY1_new(i,j),U)
                    rkV4(i,j) = interp_poly(x_input,y_input,nx_input,ny_input,gridX1_new(i,j),gridY1_new(i,j),V)
                elseif (interpolation == 'biline') then
                    rkU4(i,j) = interp_bilinear(x_input,y_input,nx_input,ny_input,gridX1_new(i,j),gridY1_new(i,j),U)
                    rkV4(i,j) = interp_bilinear(x_input,y_input,nx_input,ny_input,gridX1_new(i,j),gridY1_new(i,j),V)
                end if
            enddo
        enddo    
        !$omp end parallel do

        where (isnan(rku4)) rku4 = 0.0_DP
        where (isnan(rkv4)) rkv4 = 0.0_DP

        print*, 'fourth RK4-step computed'
    
        ! advect tracer
    
        !$omp parallel do
        do j=1,ny
            do i=1,nx
                gridX2(i,j)= gridX1(i,j) - (tstep_sec*((rkU1(i,j)+2.0_DP*rkU2(i,j)+ 2.0_DP* rkU3(i,j) +rkU4(i,j))/6.0_DP))
                gridY2(i,j)= gridY1(i,j) - (tstep_sec*((rkV1(i,j)+2.0_DP*rkV2(i,j)+ 2.0_DP* rkV3(i,j) +rkV4(i,j))/6.0_DP))
            enddo   
        enddo
        !$omp end parallel do

        ende = omp_get_wtime()
        print*, '					cpu-time =		', ende-start

        print*, ''
        print*, 'tracer advected'
        print*, ''
    
    endif ! if (forward) then

end subroutine
