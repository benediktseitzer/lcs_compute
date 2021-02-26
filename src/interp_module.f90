! interpolation routines are taken from numerical recipes in FORTRAN 
module interp_module
    implicit none

    public :: polint1d,polint2d,locate_in_grid,interp_bilinear,spline,splint,splie2,splin2,interp_spline



    contains


    !##########################################################################
 	!##########################     bicubic spline    #########################
	!##########################################################################
	
    subroutine spline(x,y,n,yp1,ypn,y2)
        use parameter_module, ONLY:DP
        implicit none
        integer     :: n
        real(DP)    :: yp1,ypn,x(n),y(n),y2(n)
        integer     :: i,k
        real(DP)    :: p,qn,sig,un,u(n)
    
        if (yp1 > .99d30) then
            y2(1) = 0.0d0
            u(1) = 0.0d0
        else
            y2(1) = -0.5d0
            u(1) = (3.0d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
        end if
    
        do i = 2,n-1
            sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
            p = sig*y2(i-1)+2.0d0
            y2(i) = (sig - 1.0d0)/p
            u(i) = (6.0d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
        end do  !do i =2,n-1
    
        if (ypn > .99d30) then
            qn = 0.0d0
            un = 0.0d0
        else
            qn = 0.5d0
            un = (3.0d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
        end if
    
        y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.0d0)
    
        do k = n-1,1,-1
            y2(k) = y2(k)*y2(k+1)+u(k)
        end do
    
        return
    end subroutine
 
    subroutine splint(xa,ya,y2a,n,x,y)
        use parameter_module, ONLY:DP  
        implicit none
        integer     :: n
        real(DP)    :: x,y,xa(n),y2a(n),ya(n)
        integer     :: k,khi,klo
        real(DP)    :: a,b,h
        klo = 1
        khi = n
    
        do while (khi-klo > 1)
            k = (khi+klo)/2
            if (xa(k) > x) then
                khi = k
            else
                klo = k
            end if
        end do  !do while (khi-klo > 1)
    
        h = xa(khi)-xa(klo)
    
        if (h == 0.d0) pause
    
        a = (xa(khi)-x)/h
        b = (x- xa(klo))/h  
        y = a*ya(klo)+b*ya(khi)+((a**3.0d0-a)*y2a(klo)+(b**3.0d0-b)*y2a(khi))*(h**2.0d0)/6.0d0
    
        return
    end subroutine
 
    subroutine splie2(x1a,x2a,ya,m,n,y2a)
        use parameter_module, ONLY:DP
        use omp_lib
        implicit none
        integer     :: m,n,NN
        real(DP)    :: x1a(m),x2a(n),y2a(m,n),ya(m,n),ende,start
        parameter(NN=100)  
        integer     :: j,k
        real(DP)    :: y2tmp(n),ytmp(n)
        
        do j = 1,m
            do k = 1,n
                ytmp(k) = ya(j,k)
            end do  !k = 1,n
            
            call spline(x2a,ytmp,n,1.d30,1.d30,y2tmp)

            do k = 1,n
                y2a(j,k) = y2tmp(k)
            end do  !k = 1,n
        end do  !j = 1,m

        return
    end subroutine !splie2


    subroutine splin2(x1a,x2a,ya,y2a,m,n,x1,x2,y)
        use parameter_module, ONLY:DP
        implicit none
        integer     :: m,n,NN
        real(DP)    :: x1,x2,y,x1a(m),x2a(n),y2a(m,n),ya(m,n),ende,start
        parameter (NN=100)
        integer     :: j,k
        real(DP)    :: y2tmp(n),ytmp(n),yytmp(n)
    

        do j = 1,m
            do k = 1,n
                ytmp(k) = ya(j,k)
                y2tmp(k) = y2a(j,k)
            end do  !k = 1,n
        
            call splint(x2a,ytmp,y2tmp,n,x2,yytmp(j))
        
        end do  !j = 1,m
    
    
        call spline(x1a,yytmp,m,1.d30,1.d30,y2tmp)
        call splint(x1a,yytmp,y2tmp,m,x1,y)
    

        return
    end subroutine
 	

	!##########################################################################
	!###########################      polinom       ###########################
	!##########################################################################

    SUBROUTINE polint1d(xa,ya,n,x,y,dy)
        use parameter_module, ONLY:DP
        implicit none
        INTEGER :: n
        REAL(DP), intent(IN) :: xa(n),ya(n)
        REAL(DP), intent(IN) :: x
        INTEGER :: i,m,ns
        REAL(DP),intent(OUT) :: dy,y
        integer,PARAMETER :: nmax=10
        REAL(DP):: den,dif,dift,ho,hp,w,C(nmax),d(nmax)
        ns=1
        dif=abs(x-xa(1))

        do i=1,n 
            dift=abs(x-xa(i))
            if (dift.lt.dif) then
                ns=i
                dif=dift
            endif
            c(i)=ya(i)
            d(i)=ya(i)
        enddo
        
        y=ya(ns)
        ns=ns-1
        
        do m=1,n-1
            do i=1,n-m
                ho=xa(i)-x
                hp=xa(i+m)-x
                w=c(i+1)-d(i)
                den=ho-hp
            !       IF(den.EQ.0.) PAUSE
                
                if (den == 0.) then
                    den = 0.001
                end if
                
                den=w/den
                d(i)=hp*den
                c(i)=ho*den
            enddo
            
            if (2.0_DP * ns .lt. n-m ) then
                dy=c(ns+1)
            else
                dy=d(ns)
                ns=ns-1
            endif
            y=y+dy
        enddo
        return
    end subroutine polint1d


    subroutine polint2d(x1a,x2a,ya,m,n,x1,x2,y,dy)
        use parameter_module, ONLY:DP
        implicit none
        integer, intent(in)   :: m,n
        real(DP) ,intent(in)  :: x1a(m),x2a(n),ya(m,n)
        real(DP), intent(in)  :: x1,x2
        real(DP), intent(out) :: y,dy
        integer:: j
        real(DP) :: ymtmp(m),yntmp(n)
        
        do j=1,m
            yntmp=ya(j,:)
            call polint1d(x2a,yntmp,n,x2,ymtmp(j),dy)
        enddo
        
        call polint1d(x1a,ymtmp,m,x1,y,dy)
        
    end subroutine polint2d


	!##########################################################################
	!#########################  locate tracer  ################################
	!##########################################################################

    integer function locate_in_grid(gridvals,n,x)
        use parameter_module, ONLY:DP
        implicit none
        integer,intent(in)  :: n
        real(DP),intent(in) :: x,gridvals(n)
        integer  :: jl, jm ,ju
        logical  :: ascnd

        jl=0; ju=n+1
        ascnd = (gridvals(n) >= gridvals(1))
        do
            if (ju-jl <= 1) exit 
                jm=(ju+jl)/2 
            if (ascnd .eqv. (x >= gridvals(jm))) then
                jl=jm 
            else
                ju=jm 
            end if
        end do
        
        if (x == gridvals(1)) then 
            locate_in_grid=1
        else if (x == gridvals(n)) then
            locate_in_grid=n-1
        else
            locate_in_grid=jl
        end if

    end function
	!##########################################################################
	!######################  calling functions  ###############################
	!##########################################################################
	
	! spline
    real(DP) recursive function interp_spline(coordsx,coordsy,nx,ny,x,y,A)
        use parameter_module, only:DP
        implicit none
        real(DP),intent(in) :: x,y
        integer,intent(in)  :: nx,ny
        real(DP),intent(in) :: coordsx(nx),coordsy(ny)
        real(DP),intent(in) :: A(nx,ny)
        integer     :: ix,iy,ki,kj
        real(DP)        :: y2a(nx,ny),ende,start
    
        !locate in grid
        ix = locate_in_grid(coordsx,nx,x)
        iy = locate_in_grid(coordsy,ny,y)
    
        
        !precompute derivative
        call splie2(coordsx,coordsy,A,nx,ny,y2a)
        
        !start 2d-spline-interpolation
        call splin2(coordsx,coordsy,A,y2a,nx,ny,x,y,interp_spline)

    end function interp_spline

	!polynomial
    real(DP) function interp_poly(coordsx,coordsy,nx,ny,x,y,A)
        use parameter_module, only: DP
        implicit none
        real(DP),intent(in) :: x,y
        integer,intent(in)  :: nx,ny
        real(DP),intent(in) :: coordsx(nx),coordsy(ny)
        real(DP),intent(in) :: A(nx,ny)
        real(DP)        :: dy
        integer     :: ix,iy,ki,kj,m_order

        m_order = 3     ! Ordnung-1 der Interpolation
        
        ix = locate_in_grid(coordsx,nx,x)
        iy = locate_in_grid(coordsy,ny,y)
        
        ki = min(max(ix-(m_order-1)/2,1),nx+1-m_order)
        kj = min(max(iy-(m_order-1)/2,1),ny+1-m_order)
        
        call polint2d(coordsx,coordsy,A(ki,kj),m_order,m_order,x,y,interp_poly,dy)
        
    end function interp_poly

	! bilinear
    real(DP) function interp_bilinear(coordsx ,coordsy ,nx, ny, x, y, A)
        use parameter_module, only: DP
        implicit none
        real(DP),intent(in)   :: x,y
        integer, intent(in)   :: nx, ny
        real(DP),intent(in)   :: coordsx(nx), coordsy(ny)
        real(DP), intent(in)  :: A(nx,ny)
        integer               :: ix,iy
        real(DP)              :: t,u,tt,uu
        ! point on new grid x,y
        ! coordinates of old grid coordsx, coordsy
        ! field: A

        ix = locate_in_grid(coordsx,nx,x)
        iy = locate_in_grid(coordsy,ny,y)

        t = ( x-coordsx(ix) ) / ( coordsx(ix+1)-coordsx(ix) )
        u = ( y-coordsy(iy) ) / ( coordsy(iy+1)-coordsy(iy) )

        tt = 1.0_DP - t
        uu = 1.0_DP - u

        interp_bilinear = tt*uu*A(ix,iy) + t*uu*A(ix+1,iy) + t*u*A(ix+1,iy+1) + tt*u*A(ix,iy+1)
        
    end function

end module
