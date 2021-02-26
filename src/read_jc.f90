program read_jc
	implicit none
	real(8), allocatable	:: vx(:),vy(:),vx_zw(:),vy_zw(:),T(:),p(:),vx_ma(:,:),vy_ma(:,:)
	real(8), allocatable	:: xgit(:),ygit(:),git(:)
	integer		:: IERR,JCOPEN,JCRINFO,JCDREAD,JCCLOSE,ICREAD
	integer		:: i,j,k,iddc,numx,numy,numz,nlen,niter,quality,time_write
	real			:: ra,rac,le,time,dt
	real(8)		:: pi,AA,BB,CSB,B,YMAX,YMIN,CS,XMAX,XMIN
	character*7		:: inputfile
	character*2		:: acc
	integer*4		:: iu,magic
	
	inputfile = 'c__data'	! filename
	acc = 'r'		! mode 'read'
	iu = 10		! input 
	magic = 76		! 
	pi = 4.0d0 * atan(1.0d0)
  
	print*, 'welcher Zeitschritt'
	read(*,*), time_write
	
	!Datein zum Schreiben oeffnen
	open(40, file = 'xgrid')
	open(50, file = 'ygrid')
	
	!open c__data
	IERR=JCOPEN(iu, magic, trim(inputfile), acc)
	print*, 'IERR =', IERR  
	
	!read header info
	IF(IERR.NE.0) PRINT*,'ERROR ICOPEN : IERR=',IERR
	IERR=JCRINFO(iu,ra,rac,iddc,le,numx,numy,numz)
	IF(IERR.NE.0) PRINT*,'ERROR ICRINFO : IERR=',IERR
	print*,'Ra=',ra,'NX=',numx,'NZ=',numy
	print*, 'iddc =', iddc
	
	!allocate vectors
	nlen=numx*numy
	allocate(vx(nlen),vy(nlen),vx_zw(nlen),vy_zw(nlen))
	allocate(T(nlen),p(nlen))
	allocate(git(2*nlen))
	allocate(xgit(nlen),ygit(nlen))
	allocate(vx_ma(numx,numy),vy_ma(numx,numy))
	
	!read grid
	IERR=ICREAD(iu,niter,time,dt,git,2*nlen)
	print*, IERR
	IF(IERR.NE.0) PRINT*,'ERROR ICREAD GITTER : IERR=',IERR
	
	!write grid
	XMAX = 3.0d0
	XMIN = 0.0d0
	YMAX = 1.0d0
	YMIN = 0.0d0
	B = 1.0d0                                
	CSB = COS(PI*B)
	AA = 0.5D0*(YMAX+YMIN)
	BB = 0.5D0*(YMAX-YMIN)/CSB
  
  	!selber stützstellen schreiben...  
	do i = 1,numx
		xgit(i) = (i-1)*(XMAX)/float(numx-1)
		write(40,*) xgit(i)
	end do
  
	do i = 1,numy
		!    ygrid(i) = (YMIN*(numy+1-i)+YMAX*(i-1))/numy
		CS = COS(PI*(B+(i-1)*(1.0D0-2.0D0*B)/(numy-1)))
		ygit(i) = AA-BB*CS
		write(50,*) ygit(i)
	end do 
  
	!Dateien zum schreiben öffnen
	open(20,file = 'veloX')
	open(30,file = 'veloY')
	
	!write for every time-step
	do while ((IERR.ne.(-1)))

		!T
		IERR=JCDREAD(iu,niter,time,dt,p,numx,numy,numz,quality)
		IF(IERR.NE.0) PRINT*,'ERROR JCDREAD T: IERR=',IERR  
		
		!vx
		IERR=JCDREAD(iu,niter,time,dt,vx,numx,numy,numz,quality)
		IF(IERR.NE.0) PRINT*,'ERROR JCDREAD vx: IERR=',IERR

		!vy
		IERR=JCDREAD(iu,niter,time,dt,vy,numx,numy,numz,quality)
		IF (IERR.NE.0) PRINT*,'ERROR JCDREAD vy: IERR=',IERR

		!P
		IERR=JCDREAD(iu,niter,time,dt,T,numx,numy,numz,quality)
		IF(IERR.NE.0) PRINT*,'ERROR JCDREAD p: IERR=',IERR  
		
		print*, ''
		print*, ''
		print*, 'Schritt =', niter,time,dt,nlen, numx, numy
		print*, ''
		print*, ''		

		!Zeitschritt lesen
		if (niter == time_write) then
			print*, 'Einlesen'
			print*, ''
			do i = 1,numx
				do j = 1,numy
					k = (j+(i-1)*numy)
					vx_zw(k) = vx((i-1)*numy+j)
					vy_zw(k) = vy((i-1)*numy+j)
				end do
			end do
			exit
		endif
	enddo

  
	!zum LCS berechnen	
	do j = 1,numy
		do i = 1,numx
			k = (j+(i-1)*numy)
			vx_ma(i,j) = vx_zw(k)
			vy_ma(i,j) = vy_zw(k)
		end do
	end do
  
	do i = 1,numx
		write(20,'(<numy>E15.7)') vx_ma(i,:)
		write(30,'(<numy>E15.7)') vy_ma(i,:)
	end do

	print*, numx,numy, time_write
	open(11,file = 'dim')
	write(11,*) numx,numy
	close(11)
	
	IERR=JCCLOSE(iu)
	close(20)
	close(30)
  end program
  
