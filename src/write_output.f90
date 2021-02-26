subroutine write_output2D()
    use parameter_module, ONLY: DP,U,V,sigma,alpha,lambda1,lambda2,nx,ny,colR,rank,x,y,eigv1_x,eigv1_y,eigv2_x,eigv2_y,Rx,Ry,x,y,gridX1,gridX0,gridY1,gridY0,div_eig
    implicit none
    character*25    :: name
    integer     :: i,j,k
    real(8)     :: pi,AA,BB,CSB,B,ZMAX,ZMIN,CS,XMAX,XMIN



    !##############################################################################
	!##############################################################################
    print*, ''
    print*, 'nlen,nx,ny', nx*ny, nx, ny
    print*, ''
     
    do j = 1,ny
        do i = 1,nx
            if (sigma(i,j) /= sigma(i,j)) then
                sigma(i,j) = 0.0_DP
            end if
        end do
    end do
	!##############################################################################
	!##############################################################################



	!##############################################################################
	!##############################################################################
	!FTLE in columns
    open(40, file = 'output_sigma_plot')
    do j = 1,ny
        do i = 1,nx
            write (40,*) x(i),y(j),sigma(i,j)
        end do
    end do  
    close(40)
	
	!FTLE matrix
    open(40, file= 'output_sigma')
    !do i=1,nx
    do i=1,nx
	    write (40,'(<ny>E15.7)') sigma(i,:)
    enddo
    close(40)

	!Norm nach Tchon et al
	open(40, file = 'output_alpha_plot')
	do j = 1,ny
		do i = 1,nx
			write(40,*) x(i),y(j),alpha(i,j)	
		end do
	end do	
	close(40)

	!alpha matrix
    open(40, file= 'output_alpha')
    !do i=1,nx
    do i=1,nx
	    write (40,'(<ny>E15.7)') alpha(i,:)
    enddo
    close(40)
	
	print*, 'soon...' 
	
	!Xgrid
	open(40, file= 'output_xgrid')
	do j=1,ny
		do i=1,nx
			write (40,*) gridX1(i,j)
		end do
	end do
	close(40)
	
	!Ygrid
	open(40, file= 'output_ygrid')
	do j=1,ny
		do i=1,nx
			write (40,*) gridY1(i,j)
		end do
	end do
	close(40)

	!Eigenvalue (larger)
	open(40, file = 'output_lambda2_plot')
	do j = 1,ny
		do i = 1,nx
			write(40,*) x(i),y(j),lambda2(i,j)
		end do
	end do
	close(40)
	
	!Eigenvalue (larger) matrix
    open(40, file= 'output_lambda2')
    !do i=1,nx
    do i=1,nx
	    write (40,'(<ny>E15.7)') lambda2(i,:)
    enddo
    close(40)

	!Eigenvalue (smaller)
	open(40, file = 'output_initcond_plot')
	do j = 1,ny
		do i = 1,nx
			write(40,*) x(i),y(j),lambda1(i,j)
		end do
	end do
	close(40)	
	
	!initcond matrix
    open(40, file= 'output_initcond')
    !do i=1,nx
    do i=1,nx
	    write (40,'(<ny>E15.7)') lambda1(i,:)
    enddo
    close(40)

	!Eigenvector to lambda2 (x-coordinate)
	open(40,file = 'output_eigv1_x_plot')
	do j = 1,ny
		do i = 1,nx
			write (40,*) x(i),y(j), eigv1_x(i,j)
			!write(40,*) x(i),y(j),U(i,j)
		end do
	end do
	close(40)
	
	print*, 'be patient...'

	!Eigenvector to lambda2 (y-coordinate)
	open(40, file = 'output_eigv1_y_plot')
	do j = 1,ny
		do i = 1,nx
			write(40,*) x(i),y(j), eigv1_y(i,j)
			!write(40,*) x(i),y(j),V(i,j)
		end do
	end do
	close(40)


	print*, 'all the data is written!'
	print*, 'i guess that is it'

	deallocate(Rx,Ry,lambda1,lambda2,eigv1_x,eigv1_y,eigv2_x,eigv2_y,alpha)
end subroutine
