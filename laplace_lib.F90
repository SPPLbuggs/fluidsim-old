    module lapl_lib
    use props
    use petsc_lib
    use ptcl_lib
    implicit none
    
    real(8), allocatable :: phi_pl(:,:)
    real(8) :: ampl
    integer :: rftype
    logical :: assem = .true.
    
    ! variables
    public  :: phi_pl
    private :: ampl, rftype, assem
    
    ! subroutines and functions
    public  :: lapl_solve, lapl_init
    private :: lapl, jacob, upd_ep

    contains
    
! *** Electromagnetic Timestepping Routine ***
    subroutine lapl_solve(g)
    type(grid), intent(in) :: g
    integer :: j, conv_reason
    
    ! Update electrode potential
    call upd_ep(g%t, phiL)
    
    ! Update boundary conditions:
    if (rx == 0) then
        do j = 2, g%by+1
            if (g%type_x(1,j-1) == -2) then
                phi_pl(1,j) = phiL
            else
                phi_pl(1,j) = phi_pl(2,j)
            end if
        end do
    end if
    
    if (rx == px-1) then
        do j = 2, g%by+1
            if (g%type_x(g%bx,j-1) == 2) then
                phi_pl(g%bx+2,j) = phiR
            else
                phi_pl(g%bx+2,j) = phi_pl(g%bx+1,j)
            end if
        end do
    end if
    
    ! Assemble jacobian and RHS
    call assem_Ab(g, phi_pl, phi_eval)
    if (assem) then
        call MatSetOption(A, Mat_New_Nonzero_Locations, PETSc_False, ierr)
        assem = .false.
    end if

    ! Solve system:
    call KSPSetOperators(ksp, A, A, ierr)
    call KSPSolve(ksp, b, x, ierr)
    
    call KSPGetConvergedReason(ksp, conv_reason, ierr)
    if ((my_id == 0) .and. (conv_reason .ne. 2)) write(*,3) conv_reason
    3 format('KSP did not converge; Reason = ',i0)
    
    ! Update variables with solution:
    call upd_soln(g, phi_pl)
    
    if (ry == 0) phi_pl(:,1) = phi_pl(:,2)
    if (ry == py-1) phi_pl(:,g%by+2) = phi_pl(:,g%by+1)
    
    end subroutine

! *** Initialization ***
    subroutine lapl_init(g, vl, intype)
    type(grid), intent(inout) :: g
    real(8), intent(in) :: vl
    integer, intent(in) :: intype
    
    allocate(phi_pl(g%bx+2, g%by+2))
    phi_pl = 0
    
    rftype = intype
    ampl = vl
    
    ! Create PETSc objects
    call petsc_create(g)
    end subroutine

! *** Evaluate phi equation ***
    subroutine phi_eval(g, i, j, n, m, phi, b_temp)
    type(grid), intent(in) :: g
    integer, intent(in) :: i, j, n, m
    real(8), intent(in) :: phi(n,m)
    real(8), intent(out) :: b_temp(1)
    
    call lapl(g, i, j, phi(:,:), b_temp(1))
    
    end subroutine

! *** Laplace's Equation ***
    subroutine lapl(g, i, j, phi, b_temp)
    type(grid), intent(in) :: g
    integer, intent(in) :: i, j
    real(8), intent(in) :: phi(:,:)
    real(8), intent(out):: b_temp
    real(8):: dphi_dx = 0, dphi_dy = 0, term_s = 0, &
              dfluxi_dx = 0, dfluxi_dy = 0, fluxi_x(2), fluxi_y(2), &
              dfluxe_dx = 0, dfluxe_dy = 0, fluxe_x(2), fluxe_y(2), &
              fluxt_x(2), fluxt_y(2)
    
    ! Calculate particle fluxes
    call calc_fluxi(g, i, j, phi, fluxi_x, fluxi_y)
    call calc_fluxe(g, i, j, phi, fluxe_x, fluxe_y, fluxt_x, fluxt_y)
    
    if (g%nx > 1) then
        ! Left boundary
        if (g%type_x(i-1,j-1) == -1) then
            dphi_dx = (phi(i+1,j) - phi(i,j)) / g%dx(i)
        
        ! Righ boundary
        else if (g%type_x(i-1,j-1) == 1) then
            dphi_dx = (phi(i-1,j) - phi(i,j)) / g%dx(i-1)
        
        ! Domain center
        else
            dphi_dx = (phi(i+1,j) - phi(i,j)) / g%dx(i) &
                      - (phi(i,j) - phi(i-1,j)) / g%dx(i-1)
        end if
        
        dfluxi_dx = (fluxi_x(2) - fluxi_x(1)) / g%dlx(i - 1)
        dfluxe_dx = (fluxe_x(2) - fluxe_x(1)) / g%dlx(i - 1)
    end if
    
    if (g%ny > 1) then
        ! Left boundary
        if (g%type_y(i-1,j-1) == -1) then
            if (cyl) then
                dphi_dy = (phi(i,j+1) - phi(i,j)) / g%dy(i) &
                          * (g%r(j+1) + g%r(j)) / 2.0
            else
                dphi_dy = (phi(i,j+1) - phi(i,j)) / g%dy(i)
            end if
        
        ! Right boundary
        else if (g%type_y(i-1,j-1) == 1) then
            if (cyl) then
                dphi_dy = (phi(i,j-1) - phi(i,j)) / g%dy(j-1) &
                          * (g%r(j) + g%r(j-1)) / 2.0
            else        
                dphi_dy = (phi(i,j-1) - phi(i,j)) / g%dy(j-1)
            end if
        
        ! Domain center
        else
            if (cyl) then
                dphi_dy = (phi(i,j+1) - phi(i,j)) / g%dy(j)     &
                          * (g%r(j+1) + g%r(j)) / 2.0           &
                          - (phi(i,j) - phi(i,j-1)) / g%dy(j-1) &
                          * (g%r(j) + g%r(j-1)) / 2.0
            else
                dphi_dy = (phi(i,j+1) - phi(i,j)) / g%dy(j) &
                          - (phi(i,j) - phi(i,j-1)) / g%dy(j-1)
            end if
        end if
        
        dphi_dy = dphi_dy / g%dly(j-1)
        if (cyl) dphi_dy = dphi_dy / g%r(j)
        
        if (cyl) then
            dfluxi_dy = (fluxi_y(2) * (g%r(j+1) + g%r(j)) / 2.0   &
                         - fluxi_y(1) * (g%r(j) + g%r(j-1)) / 2.0) &
                         / g%dly(j-1) / g%r(j)
            dfluxe_dy = (fluxe_y(2) * (g%r(j+1) + g%r(j)) / 2.0   &
                         - fluxe_y(1) * (g%r(j) + g%r(j-1)) / 2.0) &
                         / g%dly(j-1) / g%r(j)
         else
            dfluxi_dy = (fluxi_y(2) - fluxi_y(1)) / g%dly(j - 1)
            dfluxe_dy = (fluxe_y(2) - fluxe_y(1)) / g%dly(j - 1)
        end if
    end if
    
    ! Source term: e/epsilon*(n_i-n_e + dt*(del . flux_e - del . flux_i))
    term_s = -e / (eps0 * phi0 * x0) * (ni(i,j,2) - ne(i,j,2) + &
             g%dt * (dfluxe_dx + dfluxe_dy - dfluxi_dx - dfluxe_dx))
    
    b_temp = dphi_dx + (dphi_dy - term_s) * g%dlx(i-1)
    
    if (isnan(term_s)) then
        write(*,*) 'Error. NaN detected in Laplace source term. Stop.'
        call MPI_Abort(comm, 1, ierr)
    end if

    if (isnan(b_temp)) then
        write(*,*) 'Error. NaN detected in Laplace residual. Stop.'
        call MPI_Abort(comm, 1, ierr)
    end if

    end subroutine

! *** Update Electrode Potential ***
    subroutine upd_ep( t, p )
    real(8), intent(in)  :: t
    real(8), intent(out) :: p
    
    if (rftype == 1) then
            p  = ampl * sin( 2 * pi * t ) ! 1 MHz
            
    else if (rftype == 2) then
        if (t < 2.5d-1) then
            p = ampl * sin(4 * pi * t)**3
        else
            p = 0
        end if
        
    else if (rftype == 3) then
        if (t < 5d-1) then
            p = ampl * sin(4 * pi * t)**3 * (0.75 - t) / 0.626681d0
        else
            p = 0
        end if
        
    else
        if (t < 1d0 / 16d0) then
            p = ampl * sin(8d0 * pi * t)
        else
            p = ampl
        end if
        
    end if
    end subroutine
    end module

