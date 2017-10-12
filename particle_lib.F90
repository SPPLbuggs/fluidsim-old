    module ptcl_lib
    use props
    use ptcl_props
    implicit none
    
    real(8), allocatable :: ki(:,:,:), ni(:,:,:), fluxi_x(:,:), fluxi_y(:,:), &
                            ke(:,:,:), ne(:,:,:), fluxe_x(:,:), fluxe_y(:,:), &
                            km(:,:,:), nm(:,:,:), fluxm_x(:,:), fluxm_y(:,:), &
                            kt(:,:,:), nt(:,:,:), fluxt_x(:,:), fluxt_y(:,:), &
                            rho(:,:),  jx(:,:), jy(:,:)
    real(8) :: err_prev = 1, tau_m
    
    ! variables
    public  :: ni, ne, nm, nt, rho, jx, jy, tau_m
    private :: ki, ke, km, kt, fluxi_x, fluxi_y, fluxe_x, fluxe_y,&
               fluxm_x, fluxm_y, fluxt_x, fluxt_y, err_prev
    
    ! subroutines
    public  :: p_step, p_init
    private :: continuity, get_flux
    
    contains

! *** Particle Timestepping ***
    subroutine p_step(g, phi)
    type(grid), intent(inout) :: g
    real(8), intent(in) :: phi(:,:)
    integer :: i, j, stage
    real(8) :: err_cur, scfac, err_ni, err_ne, err_nt, err_nm
    
    stage = 0
    tau_m = 1
    do
        stage = stage + 1
        if (stage > 5) exit
        
        ! Update particle fluxes
        call upd_fluxi(g, phi)
        call upd_fluxe(g, phi)
        call upd_fluxm(g)
        
        ! Update particle densities
        do j = 2, g%by+1
            do i = 2, g%bx+1
                call continuity(g, i, j, stage, phi)
            end do
        end do
        
        call rk_step(g, ni, ki, stage, g%dt, err_ni, 4)
        call rk_step(g, ne, ke, stage, g%dt, err_ne, 4)
        call rk_step(g, nt, kt, stage, g%dt, err_nt, 4)
        call rk_step(g, nm, km, stage, g%dt, err_nm, 4)
        
        if (stage == 5) then
            err_cur = (err_ni**2 + err_ne**2 + &
                       err_nm**2 + err_nt**2)**0.5
            scfac = 0.8 * err_cur**(-0.7 / 4.) &
                    * err_prev**( 0.4 / 4.)
            scfac = min(2.5d0, max(3d-1, scfac))
            g%dt = scfac * g%dt
            
            g%dt = min(max(g%dt, 1d-9), tau_m/2.0)
            
            call MPI_Bcast(g%dt, 1, etype, 0, comm, ierr)
           
            if (g%dt <= 1.001d-8) then
                write(*,*) 'minimum timestep reached; finishing simulation'
                write(*,'(es10.2)') err_cur
                stop
            end if
            
            if (err_cur .le. 1.0) then
                err_prev = err_cur
            else
                stage = 0
            end if
        end if
    end do
    
    rho = (ni(:,:,1) - ne(:,:,1)) + rho
    ni(:,:,3) = ni(:,:,2)
    ne(:,:,3) = ne(:,:,2)
    nm(:,:,3) = nm(:,:,2)
    nt(:,:,3) = nt(:,:,2)
    
    ni(:,:,2) = ni(:,:,1)
    ne(:,:,2) = ne(:,:,1)
    nm(:,:,2) = nm(:,:,1)
    nt(:,:,2) = nt(:,:,1)
    
    end subroutine

! *** Particle Initialization ***
    subroutine p_init(g)
    type(grid), intent(inout) :: g
    
    allocate( ki(g%bx+2, g%by+2, 5), ni(g%bx+2, g%by+2, 3), &
              fluxi_x(g%bx+1,g%by), fluxi_y(g%bx,g%by+1), &
              ke(g%bx+2, g%by+2, 5), ne(g%bx+2, g%by+2, 3), &
              fluxe_x(g%bx+1,g%by), fluxe_y(g%bx,g%by+1), &
              kt(g%bx+2, g%by+2, 5), nt(g%bx+2, g%by+2, 3), &
              fluxt_x(g%bx+1,g%by), fluxt_y(g%bx,g%by+1), &
              km(g%bx+2, g%by+2, 5), nm(g%bx+2, g%by+2, 3), &
              fluxm_x(g%bx+1,g%by), fluxm_y(g%bx,g%by+1), &
              rho(g%bx+2, g%by+2),   jx(g%bx+3, g%by+2), jy(g%bx+2, g%by+3) )
    
    ki = 0
    ke = 0
    kt = 0
    km = 0
    ni = n_init
    ne = n_init
    nt = n_init / phi0
    nm = n_init
    rho = 0
    jx = 0
    jy = 0
    end subroutine

! *** Particle Continuity ***
    subroutine continuity(g, i, j, stage, phi)
    type(grid), intent(in) :: g
    integer, intent(in)    :: i, j, stage
    real(8), intent(in)    :: phi(:,:)
    real(8) :: Ex = 0, Ey = 0, Te, &
               dfluxi_dx = 0, dfluxi_dy = 0, fluxi_xm = 0, fluxi_ym = 0, &
               dfluxe_dx = 0, dfluxe_dy = 0, fluxe_xm = 0, fluxe_ym = 0, &
               dfluxt_dx = 0, dfluxt_dy = 0, &
               dfluxm_dx = 0, dfluxm_dy = 0, &
               term_sie, term_st1, term_st2, term_st3, term_sm, &
               k_ir, k_ex, k_sc, k_si, nu
    
    ! X-dir flux gradients
    dfluxi_dx = (fluxi_x(i,j-1) - fluxi_x(i-1,j-1)) / g%dlx(i-1)
    dfluxe_dx = (fluxe_x(i,j-1) - fluxe_x(i-1,j-1)) / g%dlx(i-1)
    dfluxt_dx = (fluxt_x(i,j-1) - fluxt_x(i-1,j-1)) / g%dlx(i-1)
    dfluxm_dx = (fluxm_x(i,j-1) - fluxm_x(i-1,j-1)) / g%dlx(i-1)
    
    ! X-dir midpoint fluxes
    fluxi_xm = 0.5 * (fluxi_x(i,j-1) + fluxi_x(i-1,j-1))
    fluxe_xm = 0.5 * (fluxe_x(i,j-1) + fluxe_x(i-1,j-1))
    Ex = -((phi(i+1,j) - phi(i,j)) / g%dx(i) &
         +(phi(i,j) - phi(i-1,j)) / g%dx(i-1) &
         ) / 2.0
    
    if (g%ny > 1) then
        
        ! Y-dir flux gradients
        if (cyl) then
            dfluxi_dy = (fluxi_y(i-1,j) * (g%r(j+1) + g%r(j)) / 2.0   &
                         - fluxi_y(i-1,j-1) * (g%r(j) + g%r(j-1)) / 2.0) &
                         / g%dly(j-1) / g%r(j)
            dfluxe_dy = (fluxe_y(i-1,j) * (g%r(j+1) + g%r(j)) / 2.0   &
                         - fluxe_y(i-1,j-1) * (g%r(j) + g%r(j-1)) / 2.0) &
                         / g%dly(j-1) / g%r(j)
            dfluxt_dy = (fluxt_y(i-1,j) * (g%r(j+1) + g%r(j)) / 2.0   &
                         - fluxt_y(i-1,j-1) * (g%r(j) + g%r(j-1)) / 2.0) &
                         / g%dly(j-1) / g%r(j)
            dfluxm_dy = (fluxm_y(i-1,j) * (g%r(j+1) + g%r(j)) / 2.0   &
                         - fluxm_y(i-1,j-1) * (g%r(j) + g%r(j-1)) / 2.0) &
                         / g%dly(j-1) / g%r(j)
        else
            dfluxi_dy = (fluxi_y(i-1,j) - fluxi_y(i-1,j-1)) / g%dly(j-1)
            dfluxe_dy = (fluxe_y(i-1,j) - fluxe_y(i-1,j-1)) / g%dly(j-1)
            dfluxt_dy = (fluxt_y(i-1,j) - fluxt_y(i-1,j-1)) / g%dly(j-1)
            dfluxm_dy = (fluxm_y(i-1,j) - fluxm_y(i-1,j-1)) / g%dly(j-1)
        end if
        
        ! Y-dir midpoint fluxes
        fluxi_ym = 0.5 * (fluxi_y(i-1,j) + fluxi_y(i-1,j-1))
        fluxe_ym = 0.5 * (fluxe_y(i-1,j) + fluxe_y(i-1,j-1))
        Ey = -((phi(i,j+1) - phi(i,j)) / g%dy(j) &
             +(phi(i,j) - phi(i,j-1)) / g%dy(j-1) &
             ) / 2.0
    end if
    
    ! rates and coefficients
    Te = get_Te(nt(i,j,2), ne(i,j,2))
    k_ir = get_k_ir(Te)
    k_sc = get_k_sc(Te)
    k_si = get_k_si(Te)
    k_ex = get_k_ex(Te)
    nu   = get_nu(Te)

    ! evaluate source terms
    term_sie =   k_ir * ninf * ne(i,j,2) &
               - beta * ni(i,j,2) * ne(i,j,2) &
               + k_si * nm(i,j,2) * ne(i,j,2) &
               + k_mp * nm(i,j,2)**2
    
    term_sm =   k_ex * ninf * ne(i,j,2) &
              - k_si * nm(i,j,2) * ne(i,j,2) &
              - k_sc * nm(i,j,2) * ne(i,j,2) &
              - k_r  * nm(i,j,2) * ne(i,j,2) &
              - 2d0  * k_mp    * nm(i,j,2)**2 &
              - k_2q * ninf * nm(i,j,2) &
              - k_3q * ninf**2 * nm(i,j,2)
    
    ! -e flux_e . E
    term_st1 = - ( fluxe_xm * Ex + fluxe_ym * Ey )

    ! -me/mg nu_e (Te - Tg)
    term_st2 = nt(i,j,2) * nu * me/mi

    ! reactions
    term_st3 =   h_ir * k_ir * ninf * ne(i,j,2) &
                + h_si * k_si * nm(i,j,2) * ne(i,j,2) &
                + h_ex * k_ex * ninf * ne(i,j,2) &
                + h_sc * k_sc * nm(i,j,2) * ne(i,j,2)

    ! evaluate expression
    ki(i,j,stage) = -dfluxi_dx - dfluxi_dy + term_sie
    ke(i,j,stage) = -dfluxe_dx - dfluxe_dy + term_sie
    kt(i,j,stage) = -dfluxt_dx - dfluxt_dy + term_st1 - term_st2 - term_st3
    km(i,j,stage) = -dfluxm_dx - dfluxm_dy + term_sm
    
    if (stage == 5) then
        tau_m = min(tau_m, eps0 * phi0 * x0 / e  &
                / (get_mue(Te) * ne(i,j,2) + mui * ni(i,j,2)))
        rho(i,j) = g%dt * (dfluxe_dx - dfluxi_dx + dfluxe_dy - dfluxi_dy)
        jx(i,j) = fluxi_xm - fluxe_xm
        jy(i,j) = fluxi_ym - fluxe_ym
    end if
    
    if (isnan(ki(i,j,stage))) then
        write(*,*) 'i,j,stage ', i,j,stage
        write(*,*) 'ki(i,j,stage)',ki(i,j,stage)
        call MPI_Finalize(ierr)
    end if
    if (isnan(ke(i,j,stage))) then
        write(*,*) 'i,j,stage ', i,j,stage
        write(*,*) 'ke(i,j,stage)',ke(i,j,stage)
        call MPI_Finalize(ierr)
    end if
    if (isnan(kt(i,j,stage))) then
        write(*,*) 'i,j,stage ', i,j,stage
        write(*,*) 'kt(i,j,stage)',kt(i,j,stage)
        call MPI_Finalize(ierr)
    end if
    if (isnan(km(i,j,stage))) then
        write(*,*) 'i,j,stage ', i,j,stage
        write(*,*) 'km(i,j,stage)',km(i,j,stage)
        call MPI_Finalize(ierr)
    end if
    
    end subroutine

! *** Update Ion Flux ***
    subroutine upd_fluxi(g, phi)
    type(grid), intent(in) :: g
    real(8), intent(in)    :: phi(:,:)
    integer :: i, j
    real(8) :: a, Ex, Ey
    
    ! X-dir Flux at i + 1/2
    do j = 1, g%by
        do i = 2, g%bx
            Ex = -(phi(i+1,j+1) - phi(i,j+1)) / g%dx(i)
            
            call get_flux(fluxi_x(i,j), Ex, g%dx(i), 1, mui, Di, &
                          ni(i,j+1,2), ni(i+1,j+1,2))
        end do
    end do
    
    ! Y-dir Flux at j + 1/2
    if (g%ny > 1) then
        do j = 2, g%by
            do i = 1, g%bx
                Ey = -(phi(i+1,j+1) - phi(i+1,j)) / g%dy(j)
                call get_flux(fluxi_y(i,j), Ey, g%dy(j), 1, mui, Di, &
                              ni(i+1,j,2), ni(i+1,j+1,2))
            end do
        end do
    end if
    
    ! X-Dir Boundary
    do j = 1, g%by
        ! Left boundary
        i = 1
        Ex = -(phi(i+1,j+1) - phi(i,j+1)) / g%dx(i)
        
        ! Case: electrode
        if (g%type_x(i,j) == -2) then
            if (-Ex > 0) then
                a = 1
            else
                a = 0
            end if
            
            fluxi_x(i,j) = a * mui * Ex * ni(i+1,j+1,2) &
                           - 0.25 * vi * ni(i+1,j+1,2)
                           
        ! Case: vacuum
        else if (g%type_x(i,j) == -1) then
            fluxi_x(i,j) = 0
        
        ! Case: interior domain (parallel only)
        else
            call get_flux(fluxi_x(i,j), Ex, g%dx(i), 1, mui, Di, &
                          ni(i,j+1,2), ni(i+1,j+1,2))
        end if
        
        ! Right boundary
        i = g%bx+1
        Ex = -(phi(i+1,j+1) - phi(i,j+1)) / g%dx(i)
        
        ! Case: electrode
        if (g%type_x(i-1,j) == 2) then
            if (Ex > 0) then
                a = 1
            else
                a = 0
            end if
            
            fluxi_x(i,j) = a * mui * Ex * ni(i,j+1,2) &
                           + 0.25 * vi * ni(i,j+1,2)
                          

        ! Case: vacuum
        else if (g%type_x(i-1,j) == 1) then
            fluxi_x(i,j) = 0
            
        ! Case: interior domain (parallel only)
        else
            call get_flux(fluxi_x(i,j), Ex, g%dx(i), 1, mui, Di, &
                          ni(i,j+1,2), ni(i+1,j+1,2))
        end if
    end do
    
    ! Y-Dir Boundary
    if (g%ny > 1) then
        do i = 1, g%bx
            ! Left boundary
            j = 1
            
            ! Case: vacuum
            if (g%type_y(i,j) == -1) then
                fluxi_y(i,j) = 0
            
            ! Case: interior domain (parallel only)
            else
                Ey = -(phi(i+1,j+1) - phi(i+1,j)) / g%dy(j)
                call get_flux(fluxi_y(i,j), Ey, g%dy(j), 1, mui, Di, &
                              ni(i+1,j,2), ni(i+1,j+1,2))
            end if
            
            ! Right boundary
            j = g%by+1
            
            ! Case: vacuum
            if (g%type_y(i,j-1) == 1) then
                fluxi_y(i,j) = 0

            ! Case: interior domain (parallel only)
            else
                Ey = -(phi(i+1,j+1) - phi(i+1,j)) / g%dy(j)
                call get_flux(fluxi_y(i,j), Ey, g%dy(j), 1, mui, Di, &
                              ni(i+1,j,2), ni(i+1,j+1,2))
            end if
        end do
    end if
    end subroutine

! *** Update Electron Flux ***
    subroutine upd_fluxe(g, phi)
    type(grid), intent(in) :: g
    real(8), intent(in)    :: phi(:,:)
    integer :: i, j
    real(8) :: a, Te(2), mue, mut, ve, Ex, Ey
    
    ! X-dir Flux at i + 1/2
    do j = 1, g%by
        do i = 2, g%bx
            ! rates and coefficients
            Te(1) = get_Te(nt(i,j+1,2), ne(i,j+1,2))
            Te(2) = get_Te(nt(i+1,j+1,2), ne(i+1,j+1,2))
            mue = 5d-1 * (get_mue(Te(1)) + get_mue(Te(2)))
            mut = 5d-1 * (get_mut(Te(1)) + get_mut(Te(2)))
            
            Ex = -(phi(i+1,j+1) - phi(i,j+1)) / g%dx(i)
            
            call get_fluxe(fluxe_x(i,j), Ex, mue, g%dx(i), &
                           ne(i:i+1,j+1,2), Te)
            call get_fluxt(fluxt_x(i,j), fluxe_x(i,j), mut, g%dx(i), &
                           nt(i:i+1,j+1,2), Te)
        end do
    end do
    
    ! Y-dir Flux at j + 1/2
    if (g%ny > 1) then
        do j = 2, g%by
            do i = 1, g%bx
                ! rates and coefficients
                Te(1) = get_Te(nt(i+1,j,2), ne(i+1,j,2))
                Te(2) = get_Te(nt(i+1,j+1,2), ne(i+1,j+1,2))
                mue = 5d-1 * (get_mue(Te(1)) + get_mue(Te(2)))
                mut = 5d-1 * (get_mut(Te(1)) + get_mut(Te(2)))
                
                Ey = -(phi(i+1,j+1) - phi(i+1,j)) / g%dy(j)
                
                call get_fluxe(fluxe_y(i,j), Ey, mue, g%dy(j), &
                               ne(i+1,j:j+1,2), Te)
                call get_fluxt(fluxt_y(i,j), fluxe_y(i,j), mut, g%dy(j), &
                               nt(i+1,j:j+1,2), Te)
            end do
        end do
    end if
    
    ! X-Dir Boundary
    do j = 1, g%by
        ! Left boundary
        i = 1
        Ex = -(phi(i+1,j+1) - phi(i,j+1)) / g%dx(i)
        
        ! Case: electrode
        if (g%type_x(i,j) == -2) then
            if (Ex > 0) then
                a = 1
            else
                a = 0
            end if
            
            Te(1) = get_Te(nt(i+1,j+1,2), ne(i+1,j+1,2))
            mue = get_mue(Te(1))
            mut = get_mut(Te(1))
            ve = sqrt( (16.0 * e * phi0 * Te(1)) / (3.0 * pi * me)) * t0 / x0
            
            fluxe_x(i,j) = - a * mue * Ex * ne(i+1,j+1,2) &
                           - 0.25 * ve * ne(i+1,j+1,2) &
                           - gam * fluxi_x(i,j)
            
            fluxt_x(i,j) = - a * mut * Ex * nt(i+1,j+1,2) &
                           - 1.0/3.0 * ve * nt(i+1,j+1,2) &
                           - gam * Te(1) * fluxi_x(i,j)

        ! Case: vacuum
        else if (g%type_x(i,j) == -1) then
            fluxe_x(i,j) = 0
            fluxt_x(i,j) = 0
        
        ! Case: interior domain (parallel only)
        else
            ! rates and coefficients
            Te(1) = get_Te(nt(i,j+1,2), ne(i,j+1,2))
            Te(2) = get_Te(nt(i+1,j+1,2), ne(i+1,j+1,2))
            mue = 5d-1 * (get_mue(Te(1)) + get_mue(Te(2)))
            mut = 5d-1 * (get_mut(Te(1)) + get_mut(Te(2)))
            
            Ex = -(phi(i+1,j+1) - phi(i,j+1)) / g%dx(i)
            
            call get_fluxe(fluxe_x(i,j), Ex, mue, g%dx(i), &
                           ne(i:i+1,j+1,2), Te)
            call get_fluxt(fluxt_x(i,j), fluxe_x(i,j), mut, g%dx(i), &
                           nt(i:i+1,j+1,2), Te)
        end if
        
        ! Right boundary
        i = g%bx+1
        Ex = -(phi(i+1,j+1) - phi(i,j+1)) / g%dx(i)
        
        ! Case: electrode
        if (g%type_x(i-1,j) == 2) then
            if (-Ex > 0) then
                a = 1
            else
                a = 0
            end if
            
            Te(1) = get_Te(nt(i,j+1,2), ne(i,j+1,2))
            mue = get_mue(Te(1))
            mut = get_mut(Te(1))
            ve = sqrt( (16.0 * e * phi0 * Te(1)) / (3.0 * pi * me)) * t0 / x0
            
            fluxe_x(i,j) = - a * mue * Ex * ne(i,j+1,2) &
                           + 0.25 * ve * ne(i,j+1,2) &
                           - gam * fluxi_x(i,j)

            fluxt_x(i,j) = - a * mut * Ex * nt(i,j+1,2) &
                           + 1.0/3.0 * ve * nt(i,j+1,2) &
                           - gam * Te(1) * fluxi_x(i,j)

        ! Case: vacuum
        else if (g%type_x(i-1,j) == 1) then
            fluxe_x(i,j) = 0
            fluxt_x(i,j) = 0
            
        ! Case: interior domain (parallel only)
        else
            ! rates and coefficients
            Te(1) = get_Te(nt(i,j+1,2), ne(i,j+1,2))
            Te(2) = get_Te(nt(i+1,j+1,2), ne(i+1,j+1,2))
            mue = 5d-1 * (get_mue(Te(1)) + get_mue(Te(2)))
            mut = 5d-1 * (get_mut(Te(1)) + get_mut(Te(2)))
            
            Ex = -(phi(i+1,j+1) - phi(i,j+1)) / g%dx(i)
            
            call get_fluxe(fluxe_x(i,j), Ex, mue, g%dx(i), &
                           ne(i:i+1,j+1,2), Te)
            call get_fluxt(fluxt_x(i,j), fluxe_x(i,j), mut, g%dx(i), &
                           nt(i:i+1,j+1,2), Te)
        end if
    end do
    
    ! Y-Dir Boundary
    if (g%ny > 1) then
        do i = 1, g%bx
            ! Left boundary
            j = 1
            
            ! Case: vacuum
            if (g%type_y(i,j) == -1) then
                fluxe_y(i,j) = 0
                fluxt_y(i,j) = 0
            
            ! Case: interior domain (parallel only)
            else
                ! rates and coefficients
                Te(1) = get_Te(nt(i+1,j,2), ne(i+1,j,2))
                Te(2) = get_Te(nt(i+1,j+1,2), ne(i+1,j+1,2))
                mue = 5d-1 * (get_mue(Te(1)) + get_mue(Te(2)))
                mut = 5d-1 * (get_mut(Te(1)) + get_mut(Te(2)))
                
                Ey = -(phi(i+1,j+1) - phi(i+1,j)) / g%dy(j)
                
                call get_fluxe(fluxe_y(i,j), Ey, mue, g%dy(j), &
                               ne(i+1,j:j+1,2), Te)
                call get_fluxt(fluxt_y(i,j), fluxe_y(i,j), mut, g%dy(j), &
                               nt(i+1,j:j+1,2), Te)
            end if
            
            ! Right boundary
            j = g%by+1
            
            ! Case: vacuum
            if (g%type_y(i,j-1) == 1) then
                fluxe_y(i,j) = 0
                fluxt_y(i,j) = 0

            ! Case: interior domain (parallel only)
            else
                ! rates and coefficients
                Te(1) = get_Te(nt(i+1,j,2), ne(i+1,j,2))
                Te(2) = get_Te(nt(i+1,j+1,2), ne(i+1,j+1,2))
                mue = 5d-1 * (get_mue(Te(1)) + get_mue(Te(2)))
                mut = 5d-1 * (get_mut(Te(1)) + get_mut(Te(2)))
                
                Ey = -(phi(i+1,j+1) - phi(i+1,j)) / g%dy(j)
                
                call get_fluxe(fluxe_y(i,j), Ey, mue, g%dy(j), &
                               ne(i+1,j:j+1,2), Te)
                call get_fluxt(fluxt_y(i,j), fluxe_y(i,j), mut, g%dy(j), &
                               nt(i+1,j:j+1,2), Te)
            end if
        end do
    end if
    end subroutine

! *** Update Metastable Flux ***
    subroutine upd_fluxm(g)
    type(grid), intent(in) :: g
    integer :: i, j
    
    ! X-dir Flux at i + 1/2
    do j = 1, g%by
        fluxm_x(:,j) = -Dm * (nm(2:,j+1,2) - nm(:g%bx+1,j+1,2)) / g%dx
    end do
    
    ! Y-dir Flux at j + 1/2
    if (g%ny > 1) then
        do i = 1, g%bx
            fluxm_y(i,:) = -Dm * (nm(i+1,2:,2) - nm(i+1,:g%by+1,2)) / g%dy
        end do
    end if
    
    ! X-Dir Boundary
    do j = 1, g%by
        ! Left boundary
        i = 1
        
        ! Case: electrode
        if (g%type_x(i,j) == -2) then
            fluxm_x(i,j) = -0.25 * vi * nm(i+1,j+1,2)

        ! Case: vacuum
        else if (g%type_x(i,j) == -1) then
            fluxm_x(i,j) = 0
        
        ! Case: interior domain (parallel only)
        else
            fluxm_x(i,j) = -Dm * (nm(i+1,j+1,2) - nm(i,j+1,2)) / g%dx(i)
        end if
        
        ! Right boundary
        i = g%bx+1
        
        ! Case: electrode
        if (g%type_x(i-1,j) == 2) then
            fluxm_x(i,j) = 0.25 * vi * nm(i,j+1,2)

        ! Case: vacuum
        else if (g%type_x(i-1,j) == 1) then
            fluxm_x(i,j) = 0
            
        ! Case: interior domain (parallel only)
        else
            fluxm_x(i,j) = -Dm * (nm(i+1,j+1,2) - nm(i,j+1,2)) / g%dx(i)
        end if
    end do
    
    ! Y-Dir Boundary
    if (g%ny > 1) then
        do i = 1, g%bx
            ! Left boundary
            j = 1
            
            ! Case: vacuum
            if (g%type_y(i,j) == -1) then
                fluxm_y(i,j) = 0
            
            ! Case: interior domain (parallel only)
            else
                fluxm_y(i,j) = -Dm * (nm(i+1,j+1,2) - nm(i+1,j,2)) / g%dy(j)
            end if
            
            ! Right boundary
            j = g%by+1
            
            ! Case: vacuum
            if (g%type_y(i,j-1) == 1) then
                fluxm_y(i,j) = 0
                
            ! Case: interior domain (parallel only)
            else
                fluxm_y(i,j) = -Dm * (nm(i+1,j+1,2) - nm(i+1,j,2)) / g%dy(j)
            end if
        end do
    end if
    end subroutine
    
    
    end module
