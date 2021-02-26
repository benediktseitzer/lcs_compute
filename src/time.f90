program time
	implicit none
	integer		:: i,j,N,Ns
	
	Ns = 35530
	N = 35540
	
	open(11,file = 'time')
	
	do i = Ns,N,10
		
		write(11,*) i
		write(11,*) i
		
	end do
	
	close(11)
end program
