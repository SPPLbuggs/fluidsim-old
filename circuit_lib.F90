! *** External Circuit Module ***
    module circ_lib
    use props
    use ptcl_lib
    use ptcl_props
    implicit none
    
    real(8) :: Vmax, Vd, Id, Vsrc, Cap, Res
    
    ! variables
    public  :: Vd, Id
    private :: Vmax, Vsrc, Cap, Res
    
    ! subroutines
    public  :: circ_step, circ_init
    !private :: 
    
    contains

! *** External Circuit Timestepping ***
    subroutine circ_step(g, phi)
    type(grid), intent(in) :: g
    real(8), intent(in) :: phi(:,:)
    integer :: i, j
    real(8) :: a, fluxi, fluxe, Ex, Te, mue, ve, Rtemp
    
    if (g%t < 1d0 / 16d0) then
        Vsrc = Vmax * sin(8d0 * pi * g%t)**3
    else
        Vsrc = Vmax
    end if
    
    Id = 0
    if (rx == 0) then
        i = 2
        do j = 2, g%by + 1
            Ex = -(phi(i,j) - phi(i-1,j)) / g%dx(i-1)
            
            if (Ex > 0) then
                a = 1
            else
                a = 0
            end if
            
            Te  = get_Te(nt(i,j,2), ne(i,j,2))
            mue = get_mue(Te)
            ve  = sqrt((16.0 * e * phi0 * Te) / (3.0 * pi * me)) * t0 / x0
            
            ! Flux at i - 1/2
            fluxi = (1 - a) * mui * Ex * ni(i,j,2) - 0.25 * vi * ni(i,j,2)
            
            fluxe = - a * mue * Ex * ne(i,j,2) &
                    - 0.25 * ve * ne(i,j,2) &
                    - gam * fluxi
            
            if (g%ny > 1) then
                if (cyl) then
                    Id = Id + (fluxi + fluxe) * g%dy(j-1) * g%r(j-1) * 2 * pi
                else
                    Id = Id + (fluxi + fluxe) * g%dy(j-1)
                end if
            else
                Id = Id + (fluxi + fluxe) * g%w**2 * pi
            end if
        end do
    end if
    
    call MPI_Allreduce(MPI_In_Place, Id, 1, etype, MPI_Sum, comm, ierr)
    
    if (g%t < 7d-2) then
        Rtemp = min(1e4 * e / (phi0 * t0), Res)
    else
        Rtemp = Res
    end if
    
    Vd = Vd + g%dt / Cap * (Id - (Vd - Vsrc)/Rtemp)
    end subroutine

! *** External Circuit Initialization ***
    subroutine circ_init(Vset, R0)
    real(8), intent(in) :: Vset, R0
    
    Vmax = Vset
    Vsrc = 0
    Vd = 0
    Cap = 1e-12 * phi0 / e
    Res = R0 * e / (phi0 * t0)
    
    end subroutine
    end module
