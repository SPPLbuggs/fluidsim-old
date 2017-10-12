    module lapl_lib
    use props
    use petsc_lib
    use ptcl_lib
    implicit none
    
    real(8), allocatable :: phif(:,:,:)
    real(8) :: ampl
    integer :: rftype
    
    ! variables
    public  :: phif
    private :: ampl, rftype
    
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
                phif(1,j,1) = phiL
            else
                phif(1,j,1) = phif(2,j,1)
            end if
        end do
    end if
    
    if (rx == px-1) then
        do j = 2, g%by+1
            if (g%type_x(g%bx,j-1) == 2) then
                phif(g%bx+2,j,1) = phiR
            else
                phif(g%bx+2,j,1) = phif(g%bx+1,j,1)
            end if
        end do
    end if
    
    ! Assemble new residual b:
    call assem_b(g, g%node, 1, phif, phif_eval)
    
    !call view
    
    ! Solve system:
    call KSPSolve(ksp, b, b, ierr)
    
    call KSPGetConvergedReason(ksp, conv_reason, ierr)
    if ((my_id == 0) .and. (conv_reason .ne. 2)) write(*,3) conv_reason
3 format('Converged not based on relative tolerance. Reason = ',i0)
    
    ! Update Phi with new values:
    call upd_soln(g, g%node, 1, phif)
    
    if (ry == 0) phif(:,1,1) = phif(:,2,1)
    if (ry == py-1) phif(:,g%by+2,1) = phif(:,g%by+1,1)
    
    end subroutine

! *** Initialization ***
    subroutine lapl_init(g, vl, intype)
    type(grid), intent(inout) :: g
    real(8), intent(in) :: vl
    integer, intent(in) :: intype
    
    allocate(phif(g%bx+2, g%by+2,1))
    phif = 0
    
    rftype = intype
    ampl = vl
    
    ! Create PETSc objects and assemble
    call petsc_create(g%nloc, g%nglob, 1)
    
    call assem_Ab(g, g%node, 1, phif, phif_eval)
    end subroutine

! *** Evaluate phi equation ***
    subroutine phif_eval(g, i, j, n, m, dof, phi, b_temp)
    type(grid), intent(in) :: g
    integer, intent(in) :: i, j, n, m, dof
    real(8), intent(in) :: phi(n,m,dof)
    real(8), intent(out) :: b_temp(dof)
    
    call lapl(g, i, j, phi(:,:,1), b_temp(1))
    
    end subroutine

! *** Laplace's Equation ***
    subroutine lapl(g, i, j, phi, b_temp)
    type(grid), intent(in) :: g
    integer, intent(in) :: i, j
    real(8), intent(in) :: phi(:,:)
    real(8), intent(out):: b_temp
    real(8):: dphi_dx = 0, dphi_dy = 0, term_s = 0
    
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
        
    end if
    
    ! Source term: e/epsilon*(n_i-n_e + dt*(del . flux_e - del . flux_i))
    term_s = e / (eps0 * phi0 * x0) * rho(i,j)
    
    b_temp = dphi_dx + (dphi_dy + term_s) * g%dlx(i-1)
    
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

