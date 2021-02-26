! this is the module where important parameters are saved
module parameter_module
    implicit none
    integer, parameter      :: SP = kind(1.0)
    integer, parameter      :: DP = selected_real_kind(2*precision(1.0_SP))
    real(DP), parameter     :: PI = 3.141592653589793238462643383279502884_DP

    ! velocity fields
    real(DP),allocatable    :: U(:,:),V(:,:)    
    
    ! coordinates
    real(DP), allocatable   :: x_input(:),y_input(:),x(:),y(:)
    
    ! grids and FTLE-Field
    real(DP), allocatable   :: sigma(:,:), gridX0(:,:),gridX1(:,:),gridX2(:,:),gridY0(:,:),gridY1(:,:),gridY2(:,:)    
    real(DP), allocatable	:: alpha(:,:), lambda1(:,:), lambda2(:,:), eigv1_x(:,:), eigv1_y(:,:), eigv2_x(:,:), eigv2_y(:,:),div_eig(:,:)
	real(DP), allocatable	:: Rx(:,:),Ry(:,:)
	real(DP)                :: dx,dy,dz
    integer                 :: nx_input, ny_input, nx, ny, colR, Ndim
    real(DP)                :: x_res, y_res, lf, lmin, stepsize
    
        
    !timestep and time intervall
    real(DP)                :: dt,t
    integer                 :: steps    

    ! euler or runge-kutta?
    logical                 :: euler

    ! integration direction: backward or forward?
    logical                 :: forward

    ! mpi
    integer                 :: rank,numprocs

    character*20            :: filename 

end module
