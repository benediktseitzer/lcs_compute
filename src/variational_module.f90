!here are all the functions and subroutines important for computing LCS using the variational theory
module variational_module
	implicit none
	
	public				:: calc_conditionB,calc_ontherun

	contains
	


	!##################################################################################################
	!##################################################################################################
	!this guy calculates the euclidean inner product and the derivatives for condition B
	real(DP) function calc_conditionB(i,j)
		use parameter_module, only: DP,nx,ny,lambda2,eigv2_x,eigv2_y,gridX0,gridY0
		implicit none
		
		integer			:: i,j
		real(DP)		:: derx,dery
		character*2		:: order
		
		! order of derivative for laplace-operator
		order = 'o4'
		
		if (((i-1)*(i-nx)<0) .and. ((j-1)*(j-ny)<0)) then
			! derivatives in second argument of product
			if (order=='o2') then
				! Order 2
				derx = (lambda2(i+1,j)-2.0_DP*lambda2(i,j)+lambda2(i-1,j))/(gridX0(i+1,j)-gridX0(i-1,j))**2.0_DP
				dery = (lambda2(i,j+1)-2.0_DP*lambda2(i,j)+lambda2(i,j-1))/(gridY0(i,j+1)-gridY0(i,j-1))**2.0_DP
			elseif (order=='o4') then
				! Order 4
				derx = (-lambda2(i+2,j)+16.0_DP*lambda2(i+1,j)-30.0_DP*lambda2(i,j)+16.0_DP*lambda2(i-1,j)-lambda2(i-2,j))/(12.0_DP*(gridX0(i+2,j)-gridX0(i-2,j))**2.0_DP)
				dery = (-lambda2(i,j+2)+16.0_DP*lambda2(i,j+1)-30.0_DP*lambda2(i,j)+16.0_DP*lambda2(i,j-1)-lambda2(i,j-2))/(12.0_DP*(gridY0(i,j+2)-gridY0(i,j-2))**2.0_DP)
			end if

			!euclidean inner product
			calc_conditionB = eigv2_x(i,j)*(derx+dery)*eigv2_x(i,j) + eigv2_y(i,j)*(derx+dery)*eigv2_y(i,j)
		end if
	end function
	
	
	
	!##################################################################################################
	!##################################################################################################
	!checks condition (A) and (B) of Theorem 1 on the run to see if integrated particle is still okay!
	integer function calc_ontherun(r,conditionB)
		use parameter_module, only: DP,x,y,nx,ny,lambda2,lambda1,alpha
		use interp_module
		implicit none
		
		integer			:: i,j,k
		real(DP)		:: l1,l2,conditionB(nx,ny),cB,a
		real(DP)		:: r(2)
		
		cB = 0.0_DP; l1 = 0.0_DP; l2 = 0.0_DP
		
		cB = interp_spline(x,y,nx,ny,r(1),r(2),conditionB)
		a = interp_spline(x,y,nx,ny,r(1),r(2),alpha)
		
		if (cB<=0.0_DP .and. a>= 0.75_DP) then
			calc_ontherun = 1
		else
			calc_ontherun = 0
		end if
	end function
	
	
	
	!##################################################################################################
	!##################################################################################################
	!whereami
	subroutine whereami(isdiscontinuous,r0,disconti)
		use parameter_module, only: DP,x_input,y_input,nx,ny,eigv1_x,eigv1_y
		implicit none		

		integer			:: ix1,iy1,ix2,iy2,ix3,iy3,ix4,iy4
		integer			:: isdiscontinuous,disconti(nx,ny)
		real(DP)		:: r0(2)
		real(DP)		:: min_x,max_x,min_y,max_y,dx,dy,pi
		real(DP)		:: phi1,phi2,phi3,phi4

		min_x=minval(x_input); min_y=minval(y_input); max_x=maxval(x_input); max_y=maxval(y_input)
		dx = (max_x-min_x)/(nx-1)
		dy = (max_y-min_y)/(ny-1)
		pi = 4.0_DP*atan(1.0_DP)
		isdiscontinuous = 0		
		
		
		!corner1 - upper right:
		ix1 = ceiling((r0(1)-min_x)/dx)+1
		iy1 = ceiling((r0(2)-min_y))+1
		if (disconti(ix1,iy1) == 1) then
			isdiscontinuous = 1
		end if			
		!corner2 - upper left:
		ix2 = ix1-1
		iy2 = iy1
		if (disconti(ix2,iy2) == 1) then
			isdiscontinuous = 1
		end if
		!corner3 - lower left:
		ix3 = ix1-1
		iy3 = iy1-1
		if (disconti(ix3,iy3) == 1) then
			isdiscontinuous = 1
		end if
		!corner4 - lower right:
		ix4 = ix1
		iy4 = iy1-1
		if (disconti(ix4,iy4) ==1) then
			isdiscontinuous = 1
		end if
	end subroutine
	
	
	
	!##################################################################################################
	!##################################################################################################
	!reverses the sign at orientational discontinuities
	subroutine reverse(fnx,fny,fox,foy)
		use parameter_module, only: DP,div_eig
		implicit none

		integer			:: i,j,nx,ny
		real(DP) 		:: pi,fnx,fny,fox,foy,phin,phio,critang,angle

		pi = 4.0_DP*atan(1.0_DP)
		critang = 90.0_DP		

		!compute angles
		if (fnx>=0.0_DP .and. fny>=0.0_DP) then			! I.
			phin = datan2(fny,fnx)*180.0_DP/pi
			phin = phin
		else if (fnx<=0.0_DP .and. fny>=0.0_DP) then		! II.
			phin = datan2(fny,fnx)*180.0_DP/pi
			phin = phin
		else if (fnx<=0.0_DP .and. fny<=0.0_DP) then		! III.
			phin = datan2(fny,fnx)*180.0_DP/pi
			phin = phin + 360.0_DP
		else if (fnx>=0.0_DP .and. fny<=0.0_DP) then		! IV.
			phin = datan2(fny,fnx)*180.0_DP/pi
			phin = phin + 360.0_DP
		end if

		if (fox>=0.0_DP .and. foy>=0.0_DP) then
			phio = datan2(foy,fox)*180.0_DP/pi
			phio = phio
			!print*, 'I.', phio,phin
		else if (fox<=0.0_DP .and. foy>=0.0_DP) then
			phio = datan2(foy,fox)*180.0_DP/pi
			phio = phio 
			!print*, 'II.', phio, phin
		else if (fox<=0.0_DP .and. foy<=0.0_DP) then
			phio = datan2(foy,fox)*180.0_DP/pi
			phio = phio + 360.0_DP
			!print*, 'III.', phio, phin
		else if (fox>=0.0_DP .and. foy<=0.0_DP) then
			phio = datan2(foy,fox)*180.0_DP/pi
			phio = phio + 360.0_DP
			!print*, 'IV.', phio, phin
		end if
		!print*, ''
		
		!ask and reverse
		if (phio<=critang) then												! I.
			if (phin<=phio+critang .or. phin>=360.0_DP-(critang-phio)) then
			else
				fnx = -1.0_DP*fnx; fny = -1.0_DP*fny
			end if
		else if (phio>critang .and. phio<(360.0_DP-critang)) then 			! II. und III.
			if (phin>phio+critang .or. phin<phio-critang) then
				fnx = -1.0_DP*fnx; fny = -1.0_DP*fny
			end if
		else if (phio>=360.0_DP-critang) then								! IV.
			if (phin>=phio-critang .or. phin<=abs(critang-(360.0_DP-phio))) then
			else
				fnx = -1.0_DP*fnx; fny = -1.0_DP*fny
			end if
		end if		
	end subroutine
	


	!##################################################################################################
	!##################################################################################################
	!compute alpha*xi
	subroutine compute_axi(isdiscontinuous,r0,foldx,foldy,fnewx,fnewy)
		use parameter_module, only: DP,x,y,nx,ny,alpha,eigv1_x,eigv1_y
		use interp_module
		implicit none
		
		integer			:: isdiscontinuous
		real(DP)		:: foldx,foldy,fnewx,fnewy,r0(2)
		
		if (isdiscontinuous == 1) then
			fnewx = foldx
			fnewy = foldy
		else
			fnewx = interp_bilinear(x,y,nx,ny,r0(1),r0(2),alpha) * interp_bilinear(x,y,nx,ny,r0(1),r0(2),eigv1_x)
			fnewy = interp_bilinear(x,y,nx,ny,r0(1),r0(2),alpha) * interp_bilinear(x,y,nx,ny,r0(1),r0(2),eigv1_y)
		end if
	end subroutine



	!##################################################################################################
	!##################################################################################################
	!solve ode: r'=xi_1(r)
	!we use a rk4 to calculate the trajectories
	subroutine ode_vec(fx,fy,r,r0,stepsize,disconti)
		use parameter_module, only: DP,nx,ny,x,y,eigv1_x,eigv1_y,gridX0,gridY0,alpha,div_eig
		use interp_module
		implicit none
		
		integer			:: i,j,k,isdiscontinuous,disconti(nx,ny)
		real(DP)		:: r(2),r0(2),rzw(2),fx,fy
		real(DP) 		:: rkx1,rkx2,rkx3,rkx4,rky1,rky2,rky3,rky4,stepsize
		
		rkx1 = 0.0_DP; rkx2 = 0.0_DP; rkx3 = 0.0_DP; rkx4 = 0.0_DP
		rky1 = 0.0_DP; rky2 = 0.0_DP; rky3 = 0.0_DP; rky4 = 0.0_DP
		rzw = 0.0_DP; isdiscontinuous = 0
	
		!compute rk-coeffs:
		rkx1 = fx
		rky1 = fy
		
		rzw(1) = r0(1) + ((stepsize * rkx1)/2.0_DP)
		rzw(2) = r0(2) + ((stepsize * rky1)/2.0_DP)
		!call whereami(isdiscontinuous,rzw,disconti)
		rkx2 = interp_bilinear(x,y,nx,ny,r0(1),r0(2),alpha) * interp_bilinear(x,y,nx,ny,r0(1),r0(2),eigv1_x)
		rky2 = interp_bilinear(x,y,nx,ny,r0(1),r0(2),alpha) * interp_bilinear(x,y,nx,ny,r0(1),r0(2),eigv1_y)
		!call compute_axi(isdiscontinuous,rzw,rkx1,rky1,rkx2,rky2)
		call reverse(rkx2,rky2,rkx1,rky1)

		rzw(1) = r0(1) + ((stepsize * rkx2)/2.0_DP)
		rzw(2) = r0(2) + ((stepsize * rky2)/2.0_DP)
		
		!call whereami(isdiscontinuous,rzw,disconti)
		rkx3 = interp_bilinear(x,y,nx,ny,r0(1),r0(2),alpha) * interp_bilinear(x,y,nx,ny,r0(1),r0(2),eigv1_x)
		rky3 = interp_bilinear(x,y,nx,ny,r0(1),r0(2),alpha) * interp_bilinear(x,y,nx,ny,r0(1),r0(2),eigv1_y)
		!call compute_axi(isdiscontinuous,rzw,rkx2,rky2,rkx3,rky3)
		call reverse(rkx3,rky3,rkx2,rky2)
			
		rzw(1) = r0(1) + (stepsize * rkx3)
		rzw(2) = r0(2) + (stepsize * rky3)
		
		!call whereami(isdiscontinuous,rzw,disconti)
		rkx4 = interp_bilinear(x,y,nx,ny,r0(1),r0(2),alpha) * interp_bilinear(x,y,nx,ny,r0(1),r0(2),eigv1_x)
		rky4 = interp_bilinear(x,y,nx,ny,r0(1),r0(2),alpha) * interp_bilinear(x,y,nx,ny,r0(1),r0(2),eigv1_y)
		!call compute_axi(isdiscontinuous,rzw,rkx3,rky3,rkx4,rky4)
		call reverse(rkx4,rky4,rkx3,rky3)

		r(1) = r0(1) + (stepsize * (rkx1+2.0_DP*rkx2+2.0_DP*rkx3+rkx4)/6.0_DP)
		r(2) = r0(2) + (stepsize * (rky1+2.0_DP*rky2+2.0_DP*rky3+rky4)/6.0_DP)

	end subroutine
	
	
	
	!##################################################################################################
	!##################################################################################################
	!needed to compute length of LCS
	real(DP) function calc_length(row,column)
		use parameter_module, only: DP,nx,ny,x,y,Ndim,colR,Rx,Ry
		use interp_module
		implicit none
		
		integer			:: i,j,k,column,row
		real(DP)		:: leeength
		
		leeength = 0.0_DP
		
		do i = 1,row-2
			leeength = leeength + sqrt((abs(Rx(i+1,column)-Rx(i,column)))**2.0_DP + (abs(Ry(i+1,column)-Ry(i,column)))**2.0_DP)
		end do
		calc_length = leeength

	end function
	
	
	
	!##################################################################################################
	!##################################################################################################
	!needed to copmute mean of lambda2 on LCS
	real(DP) function calc_mean(row,column)
		use parameter_module, only: DP,nx,ny,x,y,Ndim,colR,lambda2,Rx,Ry
		use interp_module
		implicit none
		
		integer			:: i,j,k,row,column
		real(DP)		:: sum_mean
		
		sum_mean = 0.0_DP
		do  i = 1,row-2
			sum_mean = sum_mean + interp_bilinear(x,y,nx,ny,Rx(i,column),Ry(i,column),lambda2)
		end do
		calc_mean = sum_mean		
	end function
	
end module
