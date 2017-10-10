    module lapl_lib
    use props
    use ptcl_lib
    implicit none
    
    real(8), allocatable :: phi(:,:)
    real(8) :: ampl
    integer :: rftype
    
    ! variables
    public  :: phi
    private :: ampl, rftype
    
    ! subroutines and functions
    public  :: lapl_solve, lapl_init
    private :: lapl, jacob, upd_ep

    contains
    
! *** Electromagnetic Timestepping Routine ***
    subroutine lapl_solve(g, rho)
    type(grid), intent(in) :: g
    real(8), intent(in) :: rho(:,:)
    integer :: i, j, node, conv_reason
    real(8) :: soln, b_temp
    
    ! Update electrode potential
    call upd_ep(g%t, phiL)
    
    ! Update boundary conditions
    if (rx == 0) then
        do j = 2, g%by+1
            if (g%type_x(1,j-1) == -2) then
                phi(1,j) = phiL
            else
                phi(1,j) = phi(2,j)
            end if
        end do
    end if
    
    if (rx == px-1) then
        do j = 2, g%by+1
            if (g%type_x(g%bx,j-1) == 2) then
                phi(g%bx+2,j) = phiR
            else
                phi(g%bx+2,j) = phi(g%bx+1,j)
            end if
        end do
    end if
    
    ! Assemble b for next timestep
    do j = 2, g%by+1
        do i = 2, g%bx+1
            node = g%nnode(i,j)
            b_temp = 0
            
            call lapl(g, i, j, rho(i,j), b_temp)
            call VecSetValues(b, 1, node, -b_temp, Insert_Values, ierr)
        end do
    end do
    
    call VecAssemblyBegin(b, ierr)
    call VecAssemblyEnd(b, ierr)
    
    !call view
    
    call KSPSolve(ksp, b, b, ierr)
    
    call KSPGetConvergedReason(ksp, conv_reason, ierr)
    if ((my_id == 0) .and. (conv_reason .ne. 2)) write(*,3) conv_reason
3 format('Converged not based on relative tolerance. Reason = ',i0)
    
    ! Update Phi with new values
    do j = 2, g%by+1
        do i = 2, g%bx+1
            node = g%nnode(i,j)
            call VecGetValues(b, 1, node, soln, ierr)
            phi(i, j) = phi(i, j) + soln
        end do
    end do
    
    if (ry == 0) phi(:,1) = phi(:,2)
    if (ry == py-1) phi(:,g%by+2) = phi(:,g%by+1)
    
    call comm_real(g%bx, g%by, phi)
    
    end subroutine

! *** Initialization ***
    subroutine lapl_init(g, vl, intype)
    type(grid), intent(inout) :: g
    real(8), intent(in) :: vl
    integer, intent(in) :: intype
    integer :: i, j, node, cols(5)
    real(8) :: A_temp(1,5), b_temp
    
    allocate(phi(g%bx+2, g%by+2))
    phi = 0
    
    rftype = intype
    ampl = vl
    
    ! Assemble A
    do j = 2, g%by+1
        do i = 2, g%bx+1
            b_temp =  0
            A_temp =  0
            cols   = -1
            
            node = g%nnode(i,j)
            call lapl(g, i, j, 0d0, b_temp)
            call jacob(g, i, j, b_temp, cols, A_temp)
            
            call MatSetValues(A, 1, node, 5, cols, A_temp, Insert_Values, ierr)
        end do
    end do

    call MatAssemblyBegin(A, Mat_Final_Assembly, ierr)
    call MatAssemblyEnd(A, Mat_Final_Assembly, ierr)
    !call MatSetOption(A, Mat_Symmetric, PETSC_True, ierr);
    
    call KSPCreate(comm, ksp, ierr)
    call KSPSetOperators(ksp, A, A, ierr)
    call KSPSetFromOptions(ksp, ierr)
    
    end subroutine


! *** Laplace's Equation ***
    subroutine lapl(g, i, j, rho, b_temp)
    type(grid), intent(in) :: g
    integer, intent(in) :: i, j
    real(8), intent(in) :: rho
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
    term_s = e / (eps0 * phi0 * x0) * rho
    
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

! *** Numerical Jacobian ***
    subroutine jacob(g, i_loc, j_loc, b_temp, cols, A_temp)
    type(grid), intent(in) :: g
    integer, intent(in):: i_loc, j_loc
    integer, intent(inout):: cols(5)
    real(8), intent(in):: b_temp
    real(8), intent(inout):: A_temp(1,5)
    logical:: zero_perturb
    real(8):: perturb, b_pert, temp
    integer:: I, J, K, width, k_start, k_stop
    integer, dimension(5,2):: stencil

    ! initialize
    temp = 0
    Perturb = 1e-4
    width = 0
    
    k_start = 1
    k_stop  = 5
    if (g%nx .eq. 1) k_start = 3
    if (g%ny .eq. 1) k_stop  = 3
    
    do K = k_start, k_stop
        if ((K .eq. 1) .and. (g%type_x(i_loc-1,j_loc-1) .ge. 0)) then
            width = width + 1
            stencil(width,1) = -1
            stencil(width,2) =  0
        else if ((K .eq. 2) .and. (g%type_x(i_loc-1,j_loc-1) .le. 0)) then
            width = width + 1
            stencil(width,1) = 1
            stencil(width,2) = 0
        else if (K .eq. 3) then
            width = width + 1
            stencil(width,1) = 0
            stencil(width,2) = 0
        else if ((K .eq. 4) .and. (g%type_y(i_loc-1,j_loc-1) .le. 0)) then
            width = width + 1
            stencil(width,1) = 0
            stencil(width,2) = 1
        else if ((K .eq. 5) .and. (g%type_y(i_loc-1,j_loc-1) .ge. 0)) then
            width = width + 1
            stencil(width,1) =  0
            stencil(width,2) = -1
        end if
    end do
    
    do k = 1, width
        I = i_loc + stencil(k,1)
        J = j_loc + stencil(k,2)
        
        zero_perturb = .false.
        temp = phi(I,J)
        if (Abs(phi(I,J)) > 0) then
            phi(I,J) = phi(I,J) + &
                phi(I,J)*perturb
        else
            zero_perturb = .true.
            phi(I,J) = perturb
        end if
        
        call lapl(g, i_loc, j_loc, 0d0, b_pert)
        
        if (.not. zero_perturb) then
            phi(I,J) = temp
        else
            phi(I,J) = 1
        end if

        cols(k) = g%nnode(i,j)
        A_temp(1,k) = (b_pert - b_temp) / (phi(I,J) * perturb)
        
        if (Zero_Perturb) then
            phi(I,J) = temp
        end if
    end do
    end subroutine

! *** Update Electrode Potential ***
    subroutine upd_ep( t, p )
    real(8), intent(in)  :: t
    real(8), intent(out) :: p
    
    if (rftype == 1) then
            p  = ampl * sin( 2 * pi * t )
            
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

