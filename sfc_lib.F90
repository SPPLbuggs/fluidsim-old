! *** Dielectric Surface Charge Module ***
    module sfc_lib
    use props
    implicit none
    
    real(8), allocatable :: sig(:,:)
    
    ! variables
    public  :: sig
    !private ::
    
    ! subroutines
    public  :: sfc_step, sfc_init
    !private :: 
    
    contains

! *** Surface Charge Timestepping ***
    subroutine sfc_step(g, jx)
    type(grid), intent(in) :: g
    real(8), intent(in) :: jx(:,:)
    integer :: j
    
    do j = 2, g%by + 1
        if (g%type_x(1,j-1) == -3) &
            sig(1,j) = sig(1,j) - g%dt * jx(2,j)
        if (g%type_x(g%bx,j-1) == 3) &
            sig(2,j) = sig(2,j) - g%dt * jx(g%bx+2,j)
    end do
    
    end subroutine

! *** Surface Charge Initialization ***
    subroutine sfc_init(g)
    type(grid), intent(inout) :: g
    
    allocate(sig(2, g%by+2))
    
    sig = 0
    
    end subroutine
    end module
