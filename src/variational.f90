! this will solve the algorithm for the variational theory. (Farazmand and Haller [2012])
! used functions and subroutines can be found in variational_module.f90 and interp_module.f90
subroutine variational()
	use parameter_module, only: nx,ny,x,y,colR,Ndim,gridX0,gridY0,DP,sigma,alpha,lambda1,lambda2,eigv1_x,eigv1_y,eigv2_x,eigv2_y,Rx,Ry,lf,lmin,stepsize,div_eig
	use variational_module
	use interp_module
	implicit none
	integer					:: i,j,k,n,ix,iy,ixx,iyy,abontherun,conditionAB(nx,ny),counter,isdiscontinuous,disconti(nx,ny)
	integer,allocatable		:: Nurow(:)
	real(DP)				:: L,conditionB(nx,ny),XMAX,XMIN,YMAX,YMIN,foldx,foldy,fnewx,fnewy,fzwx,fzwy,pi
	real(DP)				:: r0(2),r(2)
	real(DP),allocatable	:: length(:), lambda2_mean(:)

	XMAX = maxval(gridX0); XMIN = minval(gridX0)
	YMAX = maxval(gridY0); YMIN = minval(gridY0)
	pi = 4.0_DP * atan(1.0_DP)
	! first, determine subgrid G_0 where conditions (A) and (B) are still valid	
	conditionAB = 0; k = 0; conditionB = 0

	do iy = 1,ny
		do ix = 1,nx
			
			if (lambda1(ix,iy) /= lambda2(ix,iy) .and. lambda2(ix,iy) > 1 .and. alpha(ix,iy)>=0.5_DP) then
				conditionB(ix,iy) = calc_conditionB(ix,iy)
				if (conditionB(ix,iy) <= 0.0_DP) then
					conditionAB(ix,iy) = 1
					!if (ix==nx/3 .or. ix==nx*2/3 .or. iy==ny/3 .or. iy==ny*2/3) then
					!if (ix==2 .or. ix==nx/3 .or. ix==nx*2/3 .or. ix==nx-2 .or. iy==2 .or. iy==ny/3 .or. iy==ny*2/3 .or. iy==ny-2) then
					if ( ix==nx/5 .or. ix==nx*2/5 .or. ix==nx*3/5 .or. ix==nx*4/5 .or. iy==ny/5 .or. iy==ny*2/5 .or. iy==ny*3/5 .or. iy==ny*4/5 ) then
						conditionAB(ix,iy) = 2
						k = k+1
					end if
				end if
			end if
		end do
	end do
	
	colR = k; Ndim = 3000
	k = 0
	print*, 'colR =', colR
	
	allocate(Rx(Ndim,colR),Ry(Ndim,colR))
	allocate(length(colR),lambda2_mean(colR),Nurow(colR))	
	
	!##################################
	!heart of darkness:
 	do iy = 1,ny
 		do ix = 1,nx
 			
 			if (conditionAB(ix,iy) == 2) then
 				k = k + 1; L = 0.0_DP; n = 2
 						
 				r0(1) = gridX0(ix,iy); r0(2) = gridY0(ix,iy)
 				Rx(1,k) = r0(1); Ry(1,k) = r0(2)
 				foldx = alpha(ix,iy) * eigv1_x(ix,iy)
 				foldy = alpha(ix,iy) * eigv1_y(ix,iy)
 
 				do while(n<=Ndim .and. L < lf .and. XMAX>r0(1) .and. XMIN<r0(1) .and. YMAX>r0(2) .and. YMIN<r0(2) .and. interp_spline(x,y,nx,ny,r0(1),r0(2),alpha)>=0.5)

					!first compute alpha*eigv1_i
 					if (n>2) then
						fnewx = interp_bilinear(x,y,nx,ny,r0(1),r0(2),alpha) * interp_bilinear(x,y,nx,ny,r0(1),r0(2),eigv1_x)
						fnewy = interp_bilinear(x,y,nx,ny,r0(1),r0(2),alpha) * interp_bilinear(x,y,nx,ny,r0(1),r0(2),eigv1_y)
					else
 						fnewx = foldx
 						fnewy = foldy
 					end if
 					
					!reverse if necessary
					call reverse(fnewx,fnewy,foldx,foldy)
 					
					!solver returning r(1) and r(2)
 					call ode_vec(fnewx,fnewy,r,r0,stepsize,disconti)
 					
 					!check if trajectory is still in G_0 and probably add length of new part to L for while-condition
 					abontherun = calc_ontherun(r,conditionB)
 					if (abontherun == 0) then
 						L = L + sqrt((Rx(n-1,k)-r(1))**2.0_DP+(Ry(n-1,k)-r(2))**2.0_DP)
 					else
 						L = 0.0_DP
 					end if
 					
 					Rx(n,k) = r(1); Ry(n,k) = r(2)
 					r0 = r; foldx = fnewx; foldy = fnewy; n = n+1
 				end do !while
 
 				! save the number of rows for each trajectory
 				Nurow(k) = n
 				
 				! compute length of trajectory and mean of lambda2 over trajectory
 				if (n>3) then
 					length(k) = calc_length(Nurow(k),k)
 					lambda2_mean(k)  = calc_mean(Nurow(k),k)/length(k)
 				end if				
 				
				if (n>3) then
 					! print some stuff 
 					print*, ''
 					print*, ''
 					print*, 'length =',length(k)
 					print*, 'lambda2_mean =', lambda2_mean(k)				 
 					print*, 'L = ', L
 					print*, '                       colR    =', colR
 					print*, '                       n           =', n
 					print*, '####----> trajectory', k,'of',colR, ' computed <----####'
 				end if
			end if
 		end do
 	end do				
 

	print*, ''
	print*, '###############################################'
	print*, ''
	print*, '######## all trajectories are computed ########'
	print*, ''
 	print*, '###############################################'
 	print*, ''

	counter = 0
 	! if the trajectory is worth it... save it!
	open(40, file = 'output_lcs')
 	do j = 1,colR
 		if (length(j)>lmin .and. lambda2_mean(j)>=lambda2_mean(j-1) .and. lambda2_mean(j)>=lambda2_mean(j+1)) then
			do i = 1,Nurow(j)-1
 				 write(40,*) Rx(i,j), Ry(i,j)
 			end do
 			counter = counter + 1
 			write(40,*) '' 
 		end if
 	end do
	close(40)
	
	!zum checken wie cond (A) und (B) verteilt sind:
	lambda1 = conditionAB
	
	print*, ''
	print*, '###############################################'
	print*, ''
	print*, '########', counter, 'LCS FOUND!	########'
	print*, ''
	print*, '###############################################'
	print*, ''
end subroutine
