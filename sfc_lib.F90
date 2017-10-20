! *** Dielectric Surface Charge Module ***
    module sfc_lib
    use props
    use ptcl_lib
    use ptcl_props
    implicit none
    
    real(8), allocatable :: sig(:)
    
    ! variables
    public  :: sig
    !private ::
    
    ! subroutines
    public  :: sfc_step, sfc_init
    !private :: 
    
    contains

! *** Surface Charge Timestepping ***
    subroutine sfc_step(g, phi)
    type(grid), intent(in) :: g
    real(8), intent(in) :: phi(:,:)
    integer :: i, j
    real(8) :: a, fluxi, fluxe, Ey, Te, mue, ve
    
    if (ry == py-1) then
        j = g%by+1
        do i = 2, g%bx + 1
            Ey = -(phi(i,j+1) - phi(i,j)) / g%dy(j)
            
            if (Ey > 0) then
                a = 1
            else
                a = 0
            end if
            
            Te  = get_Te(nt(i,j,2), ne(i,j,2))
            mue = get_mue(Te)
            ve  = sqrt((16.0 * e * phi0 * Te) / (3.0 * pi * me)) * t0 / x0
            
            fluxi = a * mui * Ey * ni(i,j,2) + 0.25 * vi * ni(i,j,2)
            fluxe = (a - 1) * mue * Ey * ne(i,j,2) + 0.25 * ve * ne(i,j,2)
            sig(i) = sig(i) + g%dt * (fluxi - fluxe)
        end do
    
        call MPI_Send(sig(g%bx+1), 1, etype, east, 9, comm, ierr)
        call MPI_Send(sig(2),      1, etype, west, 9, comm, ierr)
        call MPI_Recv(sig(g%bx+2), 1, etype, east, 9, comm, stat, ierr)
        call MPI_Recv(sig(1),      1, etype, west, 9, comm, stat, ierr)
    end if
    
    call MPI_Barrier(comm, ierr)
    end subroutine

! *** Surface Charge Initialization ***
    subroutine sfc_init(g)
    type(grid), intent(inout) :: g
    
    allocate(sig(g%bx+2))
    
    sig = 0
    
    end subroutine
    end module
