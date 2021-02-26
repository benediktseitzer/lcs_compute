!calculates C-G strain tensor and FTLE
!eig_v2 corresponds to larger eigenvalue lambda2 and eig_v1 to lambda1
subroutine calculate_FTLE2D(ix,iy,t_span)
    use parameter_module, only: nx,ny,gridX1,gridX0,gridY1,gridY0,DP,sigma,alpha,lambda1,lambda2,eigv1_x,eigv1_y,eigv2_x,eigv2_y,div_eig
    use lapack95
    implicit none
    integer, intent(in)     :: ix, iy
    real(DP), intent(in)    :: t_span
    real(DP)                :: A(2,2), CG(2,2)
    real(DP)                :: eig_r(2),eig_i(2),eig_v(2,2),eig_v1(2),eig_v2(2)
  
  
    ! computing Cauchy-Green strain tensor (CG):
    if ( ((ix-1)*(ix-nx)<0) .and. ((iy - 1)*(iy - ny)<0) ) then
        A(1,1)= (gridX1(ix+1,iy) - gridX1(ix-1,iy) ) / (gridX0(ix+1,iy) - gridX0(ix-1,iy) )
        A(1,2)= (gridX1(ix,iy+1) - gridX1(ix,iy-1) ) / (gridY0(ix,iy+1) - gridY0(ix,iy-1) )
        A(2,1)= (gridY1(ix+1,iy) - gridY1(ix-1,iy) ) / (gridX0(ix+1,iy) - gridX0(ix-1,iy) ) 
        A(2,2)= (gridY1(ix,iy+1) - gridY1(ix,iy-1) ) / (gridY0(ix,iy+1) - gridY0(ix,iy-1) )

        ! CG strain tensor:
        CG = matmul(transpose(A),A)
        
		! compute eigenvalues and eigenvectors of CG strain tensor:
        call geev(CG, eig_r, eig_i,eig_v)
		lambda2(ix,iy) = maxval(eig_r)
		lambda1(ix,iy) = minval(eig_r)		
		
		! sort eigenvectors to corresponding eigenvalues:
		if (lambda1(ix,iy) == eig_r(1)) then
			eig_v1 = eig_v(:,1)
			eig_v2 = eig_v(:,2)
		else
			eig_v1 = eig_v(:,2)
			eig_v2 = eig_v(:,1)
		end if
		
		eigv1_x(ix,iy) = eig_v1(1)
		eigv1_y(ix,iy) = eig_v1(2)
		
		eigv2_x(ix,iy) = eig_v2(1)
		eigv2_y(ix,iy) = eig_v2(2)
		
		! calculate divergence of xi_1
		div_eig(ix,iy) = (eigv1_x(ix+1,iy-eigv1_x(ix-1,iy)))/(2.0_DP*(gridX0(ix+1,iy)-gridX0(ix-1,iy))) + (eigv1_y(ix,iy+1)-eigv1_y(ix,iy-1))/(2.0_DP*(gridY0(ix,iy+1)-gridY0(ix,iy-1)))
	
		! calculate norm for strainline integration:
		alpha(ix,iy) = ((lambda2(ix,iy)-lambda1(ix,iy))/(lambda2(ix,iy)+lambda1(ix,iy)))**2.0_DP
		
		! compute FTLE:
		sigma(ix,iy) = log(lambda2(ix,iy))/(2.0_DP*t_span)
    else
        sigma(ix,iy) = 0.0_DP
    end if
    
end subroutine
