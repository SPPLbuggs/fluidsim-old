    module ptcl_lib
    use props
    use ptcl_props
    implicit none
    
    real(8), allocatable :: ki(:,:,:), ni(:,:,:), &
                            ke(:,:,:), ne(:,:,:), &
                            km(:,:,:), nm(:,:,:), &
                            kt(:,:,:), nt(:,:,:)
    real(8) :: err_prev = 1, tau_m
    
    ! variables
    public  :: ni, ne, nm, nt, tau_m
    private :: ki, ke, km, kt, err_prev
    
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
        
        ! Update particle densities
        do j = 2, g%by+1
            do i = 2, g%bx+1
                call continuity(g, i, j, stage, phi)
            end do
        end do
        
        call rk_step(g, ni, ki, stage, g%dt, err_ni, 4, n_zero)
        call rk_step(g, ne, ke, stage, g%dt, err_ne, 4, n_zero)
        call rk_step(g, nt, kt, stage, g%dt, err_nt, 4, n_zero/phi0/100.)
        call rk_step(g, nm, km, stage, g%dt, err_nm, 4, n_zero)
        
        if (stage == 5) then
            err_cur = (err_ni**2 + err_ne**2 + &
                       err_nm**2 + err_nt**2)**0.5
            scfac = 0.8 * err_cur**(-0.7 / 4.) &
                    * err_prev**( 0.4 / 4.)
            scfac = min(2.5d0, max(3d-1, scfac))
            g%dt = scfac * g%dt
            
            call MPI_Allreduce(MPI_In_Place, tau_m, 1, etype, &
                               MPI_Min, comm, ierr)
            if (g%ny > 1) then
                g%dt = min(max(g%dt, 1d-12), tau_m*sqrt(2.0))
            else
                g%dt = min(max(g%dt, 1d-12), tau_m*2.0)
            end if
            
            call MPI_Bcast(g%dt, 1, etype, 0, comm, ierr)
           
            if (g%dt <= 1.1d-12) then
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
              ke(g%bx+2, g%by+2, 5), ne(g%bx+2, g%by+2, 3), &
              kt(g%bx+2, g%by+2, 5), nt(g%bx+2, g%by+2, 3), &
              km(g%bx+2, g%by+2, 5), nm(g%bx+2, g%by+2, 3) )
    
    ki = 0
    ke = 0
    kt = 0
    km = 0
    ni = n_init
    ne = n_init
    nt = n_init / phi0 / 100.
    nm = n_init
    end subroutine

! *** Particle Continuity ***
    subroutine continuity(g, i, j, stage, phi)
    type(grid), intent(in) :: g
    integer, intent(in)    :: i, j, stage
    real(8), intent(in)    :: phi(:,:)
    real(8) :: Ex = 0, Ey = 0, Te, &
               dfluxi_dx = 0, dfluxi_dy = 0, fluxi_x(2), fluxi_y(2), &
               dfluxe_dx = 0, dfluxe_dy = 0, fluxe_x(2), fluxe_y(2), &
               dfluxt_dx = 0, dfluxt_dy = 0, fluxt_x(2), fluxt_y(2), &
               dfluxm_dx = 0, dfluxm_dy = 0, fluxm_x(2), fluxm_y(2), &
               term_sie, term_st1, term_st2, term_st3, term_sm, &
               k_ir, k_ex, k_sc, k_si, nu
    
    ! Calculate particle fluxes
    call calc_fluxi(g, i, j, phi, fluxi_x, fluxi_y)
    call calc_fluxe(g, i, j, phi, fluxe_x, fluxe_y, fluxt_x, fluxt_y)
    call calc_fluxm(g, i, j, fluxm_x, fluxm_y)
    
    ! X-dir flux gradients
    dfluxi_dx = (fluxi_x(2) - fluxi_x(1)) / g%dlx(i-1)
    dfluxe_dx = (fluxe_x(2) - fluxe_x(1)) / g%dlx(i-1)
    dfluxt_dx = (fluxt_x(2) - fluxt_x(1)) / g%dlx(i-1)
    dfluxm_dx = (fluxm_x(2) - fluxm_x(1)) / g%dlx(i-1)
    
    ! X-dir midpoint fluxes
    fluxe_x(1) = 0.5 * (fluxe_x(2) + fluxe_x(1))
    Ex = -((phi(i+1,j) - phi(i,j)) / g%dx(i) &
         +(phi(i,j) - phi(i-1,j)) / g%dx(i-1) &
         ) / 2.0
    
    if (g%ny > 1) then
        
        ! Y-dir flux gradients
        if (cyl) then
            dfluxi_dy = (fluxi_y(2) * (g%r(j+1) + g%r(j)) / 2.0   &
                         - fluxi_y(1) * (g%r(j) + g%r(j-1)) / 2.0) &
                         / g%dly(j-1) / g%r(j)
            dfluxe_dy = (fluxe_y(2) * (g%r(j+1) + g%r(j)) / 2.0   &
                         - fluxe_y(1) * (g%r(j) + g%r(j-1)) / 2.0) &
                         / g%dly(j-1) / g%r(j)
            dfluxt_dy = (fluxt_y(2) * (g%r(j+1) + g%r(j)) / 2.0   &
                         - fluxt_y(1) * (g%r(j) + g%r(j-1)) / 2.0) &
                         / g%dly(j-1) / g%r(j)
            dfluxm_dy = (fluxm_y(2) * (g%r(j+1) + g%r(j)) / 2.0   &
                         - fluxm_y(1) * (g%r(j) + g%r(j-1)) / 2.0) &
                         / g%dly(j-1) / g%r(j)
        else
            dfluxi_dy = (fluxi_y(2) - fluxi_y(1)) / g%dly(j-1)
            dfluxe_dy = (fluxe_y(2) - fluxe_y(1)) / g%dly(j-1)
            dfluxt_dy = (fluxt_y(2) - fluxt_y(1)) / g%dly(j-1)
            dfluxm_dy = (fluxm_y(2) - fluxm_y(1)) / g%dly(j-1)
        end if
        
        ! Y-dir midpoint fluxes
        fluxe_y(1) = 0.5 * (fluxe_y(2) + fluxe_y(1))
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
    term_st1 = - (fluxe_x(1) * Ex + fluxe_y(1) * Ey)

    ! -me/mg nu_e (Te - Tg)
    term_st2 = - nt(i,j,2) * nu * me/mi

    ! reactions
    term_st3 = - h_ir * k_ir * ninf * ne(i,j,2) &
               - h_si * k_si * nm(i,j,2) * ne(i,j,2) &
               - h_ex * k_ex * ninf * ne(i,j,2) &
               - h_sc * k_sc * nm(i,j,2) * ne(i,j,2)
    
    ! evaluate expression
    ki(i,j,stage) = -dfluxi_dx - dfluxi_dy + term_sie
    ke(i,j,stage) = -dfluxe_dx - dfluxe_dy + term_sie
    kt(i,j,stage) = -dfluxt_dx - dfluxt_dy + term_st1 + term_st2 + term_st3
    km(i,j,stage) = -dfluxm_dx - dfluxm_dy + term_sm
    
    if (stage == 5) then
        tau_m = min(tau_m, 1.0 / (get_mue(Te) * ne(i,j,2) + mui * ni(i,j,2)))
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

! *** Calculate Ion Density Flux ***
    subroutine calc_fluxi(g, i, j, phi, fluxi_x, fluxi_y)
    type(grid), intent(in) :: g
    integer, intent(in) :: i, j
    real(8), intent(in) :: phi(:,:)
    real(8), intent(out) :: fluxi_x(2), fluxi_y(2)
    real(8) :: a, Ex(2), Ey(2)
    
    fluxi_x = 0
    fluxi_y = 0
    
    ! X-dir fields:
    Ex(1) = -(phi(i,j) - phi(i-1,j)) / g%dx(i-1)
    Ex(2) = -(phi(i+1,j) - phi(i,j)) / g%dx(i)
    
    ! X-dir Fluxes:
    ! - center -
    if (g%type_x(i-1,j-1) == 0) then
        ! Flux at i - 1/2
        call get_flux(fluxi_x(1), Ex(1), g%dx(i-1), 1, mui, Di, &
                      ni(i-1,j,2), ni(i,j,2))
        !call get_fluxi(fluxi_x(1), Ex(1), mui, Di, g%dx(i-1),&
        !               ni(i-1:i,j,2))

        ! Flux at i + 1/2
        call get_flux(fluxi_x(2), Ex(2), g%dx(i), 1, mui, Di, &
                      ni(i,j,2), ni(i+1,j,2))
        !call get_fluxi(fluxi_x(2), Ex(2), mui, Di, g%dx(i),&
        !               ni(i:i+1,j,2))
    
    ! - left -
    else if (g%type_x(i-1,j-1) < 0) then
        ! Flux at i + 1/2
        call get_flux(fluxi_x(2), Ex(2), g%dx(i), 1, mui, Di, &
                      ni(i,j,2), ni(i+1,j,2))
        !call get_fluxi(fluxi_x(2), Ex(2), mui, Di, g%dx(i),&
        !               ni(i:i+1,j,2))
        
        ! - electrode -
        if (g%type_x(i-1,j-1) == -2) then
            if (Ex(1) < 0) then
                a = 1
            else
                a = 0
            end if
            
            ! Flux at i - 1/2
            fluxi_x(1) = a * mui * Ex(1) * ni(i,j,2) - 0.25 * vi * ni(i,j,2)

        ! - vacuum -
        else if (g%type_x(i-1,j-1) == -1) then
            ! Flux at i - 1/2
            fluxi_x(1) = 0
        end if
    
    ! - right -
    else if (g%type_x(i-1,j-1) > 0) then
        ! Flux at i - 1/2
        call get_flux(fluxi_x(1), Ex(1), g%dx(i-1), 1, mui, Di, &
                      ni(i-1,j,2), ni(i,j,2))
        !call get_fluxi(fluxi_x(1), Ex(1), mui, Di, g%dx(i-1),&
        !               ni(i-1:i,j,2))
        
        ! - electrode -
        if (g%type_x(i-1,j-1) == 2) then
            if (Ex(2) > 0) then
                a = 1
            else
                a = 0
            end if
            
            ! Flux at i + 1/2
            fluxi_x(2) = a * mui * Ex(2) * ni(i,j,2) + 0.25 * vi * ni(i,j,2)
            
        ! - vacuum -
        else if (g%type_x(i-1,j-1) == 1) then
            ! Flux at i + 1/2
            fluxi_x(2) = 0
        end if
    end if
    
    ! Y-dir Fluxes
    if (g%ny > 1) then
        Ey(1) = -(phi(i,j) - phi(i,j-1)) / g%dy(j-1)
        Ey(2) = -(phi(i,j+1) - phi(i,j)) / g%dy(j)
        
        ! - center -
        if (g%type_y(i-1,j-1) == 0) then
            ! Flux at j - 1/2
            call get_flux(fluxi_y(1), Ey(1), g%dy(j-1), 1, mui, Di, &
                          ni(i,j-1,2), ni(i,j,2))
            !call get_fluxi(fluxi_y(1), Ey(1), mui, Di, g%dy(j-1),&
            !           ni(i,j-1:j,2))

            ! Flux at j + 1/2
            call get_flux(fluxi_y(2), Ey(2), g%dy(j), 1, mui, Di, &
                          ni(i,j,2), ni(i,j+1,2))
            !call get_fluxi(fluxi_y(2), Ey(2), mui, Di, g%dy(j),&
            !           ni(i,j:j+1,2))
            
        
        ! - left -
        else if (g%type_y(i-1,j-1) < 0) then
            ! Flux at j + 1/2
            call get_flux(fluxi_y(2), Ey(2), g%dy(j), 1, mui, Di, &
                          ni(i,j,2), ni(i,j+1,2))
            !call get_fluxi(fluxi_y(2), Ey(2), mui, Di, g%dy(j),&
            !           ni(i,j:j+1,2))

            ! Flux at j - 1/2
            fluxi_y(1) = 0
        
        ! - right -
        else if (g%type_y(i-1,j-1) > 0) then
            ! Flux at j - 1/2
            call get_flux(fluxi_y(1), Ey(1), g%dy(j-1), 1, mui, Di, &
                          ni(i,j-1,2), ni(i,j,2))
            !call get_fluxi(fluxi_y(1), Ey(1), mui, Di, g%dy(j-1),&
            !           ni(i,j-1:j,2))
            
            ! Flux at j + 1/2
            if (right_wall) then
                if (Ey(2) > 0) then
                    a = 1
                else
                    a = 0
                end if
                
                fluxi_y(2) = a * mui * Ey(2) * ni(i,j,2) + 0.25 * vi * ni(i,j,2)
            else
                fluxi_y(2) = 0
            end if
        end if
    end if
    
    end subroutine

! *** Update Electron Flux ***
    subroutine calc_fluxe(g, i, j, phi, fluxe_x, fluxe_y, fluxt_x, fluxt_y)
    type(grid), intent(in) :: g
    integer, intent(in) :: i, j
    real(8), intent(in) :: phi(:,:)
    real(8), intent(out) :: fluxe_x(2), fluxe_y(2), fluxt_x(2), fluxt_y(2)
    real(8) :: a, Te(3), mue(2), mut(2), ve, Ex(2), Ey(2), fluxi

    fluxe_x = 0
    fluxe_y = 0
    fluxt_x = 0
    fluxt_y = 0
    
    ! X-dir fields:
    Ex(1) = -(phi(i,j) - phi(i-1,j)) / g%dx(i-1)
    Ex(2) = -(phi(i+1,j) - phi(i,j)) / g%dx(i)
    
    ! X-dir Fluxes:
    ! - center -
    if (g%type_x(i-1,j-1) == 0) then
        ! rates and coefficients
        Te(1) = get_Te(nt(i-1,j,2), ne(i-1,j,2))
        Te(2) = get_Te(nt(i,j,2),   ne(i,j,2))
        Te(3) = get_Te(nt(i+1,j,2), ne(i+1,j,2))
        
        mue(1) = 5d-1 * (get_mue(Te(1)) + get_mue(Te(2)))
        mue(2) = 5d-1 * (get_mue(Te(2)) + get_mue(Te(3)))
        
        mut = 5./3. * mue
        
        ! Flux at i - 1/2
        call get_fluxe(fluxe_x(1), Ex(1), mue(1), g%dx(i-1), &
                       ne(i-1:i,j,2), Te(1:2))
        call get_fluxt(fluxt_x(1), fluxe_x(1), mut(1), g%dx(i-1), &
                       nt(i-1:i,j,2), Te(1:2))

        ! Flux at i + 1/2
        call get_fluxe(fluxe_x(2), Ex(2), mue(2), g%dx(i), &
                       ne(i:i+1,j,2), Te(2:3))
        call get_fluxt(fluxt_x(2), fluxe_x(2), mut(2), g%dx(i), &
                       nt(i:i+1,j,2), Te(2:3))
        
    
    ! - left -
    else if (g%type_x(i-1,j-1) < 0) then
        ! rates and coefficients
        Te(2) = get_Te(nt(i,j,2),   ne(i,j,2))
        Te(3) = get_Te(nt(i+1,j,2), ne(i+1,j,2))
        
        mue(2) = 5d-1 * (get_mue(Te(2)) + get_mue(Te(3)))
        mut(2) = 5./3. * mue(2)

        ! Flux at i + 1/2
        call get_fluxe(fluxe_x(2), Ex(2), mue(2), g%dx(i), &
                       ne(i:i+1,j,2), Te(2:3))
        call get_fluxt(fluxt_x(2), fluxe_x(2), mut(2), g%dx(i), &
                       nt(i:i+1,j,2), Te(2:3))
        
        ! - electrode -
        if (g%type_x(i-1,j-1) == -2) then
            if (Ex(1) > 0) then
                a = 1
            else
                a = 0
            end if
            
            mue(1) = get_mue(Te(2))
            mut(1) = 5./3. * mue(1)
            ve = sqrt((16.0 * e * phi0 * Te(2)) / (3.0 * pi * me)) * t0 / x0
            
            ! Flux at i - 1/2
            fluxi = (1 - a) * mui * Ex(1) * ni(i,j,2) - 0.25 * vi * ni(i,j,2)
            
            fluxe_x(1) = - a * mue(1) * Ex(1) * ne(i,j,2) &
                         - 0.25 * ve * ne(i,j,2) &
                         - gam * fluxi
            
            fluxt_x(1) = - a * mut(1) * Ex(1) * nt(i,j,2) &
                         - 1.0/3.0 * ve * nt(i,j,2) &
                         - gam * Te(2) * fluxi

        ! - vacuum -
        else if (g%type_x(i-1,j-1) == -1) then
            ! Flux at i - 1/2
            fluxe_x(1) = 0
            fluxt_x(1) = 0
        end if
    
    ! - right -
    else if (g%type_x(i-1,j-1) > 0) then
        ! rates and coefficients
        Te(1) = get_Te(nt(i-1,j,2), ne(i-1,j,2))
        Te(2) = get_Te(nt(i,j,2),   ne(i,j,2))
        
        mue(1) = 5d-1 * (get_mue(Te(1)) + get_mue(Te(2)))
        mut(1) = 5./3. * mue(1)
        
        ! Flux at i - 1/2
        call get_fluxe(fluxe_x(1), Ex(1), mue(1), g%dx(i-1), &
                       ne(i-1:i,j,2), Te(1:2))
        call get_fluxt(fluxt_x(1), fluxe_x(1), mut(1), g%dx(i-1), &
                       nt(i-1:i,j,2), Te(1:2))
        
        ! - electrode -
        if (g%type_x(i-1,j-1) == 2) then
            if (-Ex(2) > 0) then
                a = 1
            else
                a = 0
            end if
            
            mue(2) = get_mue(Te(2))
            mut(2) = 5./3. * mue(2)
            ve = sqrt((16.0 * e * phi0 * Te(2)) / (3.0 * pi * me)) * t0 / x0
            
            ! Flux at i + 1/2
            fluxi = (1 - a) * mui * Ex(2) * ni(i,j,2) + 0.25 * vi * ni(i,j,2)
            
            fluxe_x(2) = - a * mue(2) * Ex(2) * ne(i,j,2) &
                         + 0.25 * ve * ne(i,j,2) &
                         - gam * fluxi

            fluxt_x(2) = - a * mut(2) * Ex(2) * nt(i,j,2) &
                         + 1.0/3.0 * ve * nt(i,j,2) &
                         - gam * Te(2) * fluxi

        ! - vacuum -
        else if (g%type_x(i-1,j-1) == 1) then
            ! Flux at i + 1/2
            fluxe_x(2) = 0
            fluxt_x(2) = 0
        end if
    end if
    
    ! Y-dir Fluxes
    if (g%ny > 1) then
        Ey(1) = -(phi(i,j) - phi(i,j-1)) / g%dy(j-1)
        Ey(2) = -(phi(i,j+1) - phi(i,j)) / g%dy(j)
        
        ! - center -
        if (g%type_y(i-1,j-1) == 0) then
            ! rates and coefficients
            Te(1) = get_Te(nt(i,j-1,2), ne(i,j-1,2))
            Te(2) = get_Te(nt(i,j,2),   ne(i,j,2))
            Te(3) = get_Te(nt(i,j+1,2), ne(i,j+1,2))
            
            mue(1) = 5d-1 * (get_mue(Te(1)) + get_mue(Te(2)))
            mue(2) = 5d-1 * (get_mue(Te(2)) + get_mue(Te(3)))
            
            mut = 5./3. * mue
            
            ! Flux at j - 1/2
            call get_fluxe(fluxe_y(1), Ey(1), mue(1), g%dy(j-1), &
                           ne(i,j-1:j,2), Te(1:2))
            call get_fluxt(fluxt_y(1), fluxe_y(1), mut(1), g%dy(j-1), &
                           nt(i,j-1:j,2), Te(1:2))

            ! Flux at j + 1/2
            call get_fluxe(fluxe_y(2), Ey(2), mue(2), g%dy(j), &
                           ne(i,j:j+1,2), Te(2:3))
            call get_fluxt(fluxt_y(2), fluxe_y(2), mut(2), g%dy(j), &
                           nt(i,j:j+1,2), Te(2:3))
            
        
        ! - left -
        else if (g%type_y(i-1,j-1) < 0) then
            ! rates and coefficients
            Te(2) = get_Te(nt(i,j,2),   ne(i,j,2))
            Te(3) = get_Te(nt(i,j+1,2), ne(i,j+1,2))
            
            mue(2) = 5d-1 * (get_mue(Te(2)) + get_mue(Te(3)))
            mut(2) = 5./3. * mue(2)

            ! Flux at j + 1/2
            call get_fluxe(fluxe_y(2), Ey(2), mue(2), g%dy(j), &
                           ne(i,j:j+1,2), Te(2:3))
            call get_fluxt(fluxt_y(2), fluxe_y(2), mut(2), g%dy(j), &
                           nt(i,j:j+1,2), Te(2:3))
            
            ! Flux at j - 1/2
            fluxe_y(1) = 0
            fluxt_y(1) = 0
        
        ! - right -
        else if (g%type_y(i-1,j-1) > 0) then
            ! rates and coefficients
            Te(1) = get_Te(nt(i,j-1,2), ne(i,j-1,2))
            Te(2) = get_Te(nt(i,j,2),   ne(i,j,2))
            
            mue(1) = 5d-1 * (get_mue(Te(1)) + get_mue(Te(2)))
            
            mut(1) = 5./3. * mue(1)
            
            ! Flux at j - 1/2
            call get_fluxe(fluxe_y(1), Ey(1), mue(1), g%dy(j-1), &
                           ne(i,j-1:j,2), Te(1:2))
            call get_fluxt(fluxt_y(1), fluxe_y(1), mut(1), g%dy(j-1), &
                           nt(i,j-1:j,2), Te(1:2))
            
            ! Flux at j + 1/2
            if (right_wall) then
                if (-Ey(2) > 0) then
                    a = 1
                else
                    a = 0
                end if
                
                mue(2) = get_mue(Te(2))
                mut(2) = 5./3. * mue(2)
                ve = sqrt((16.0 * e * phi0 * Te(2)) / (3.0 * pi * me)) * t0 / x0
                
                ! Flux at i + 1/2
                fluxe_y(2) = - a * mue(2) * Ey(2) * ne(i,j,2) &
                             + 0.25 * ve * ne(i,j,2)

                fluxt_y(2) = - a * mut(2) * Ey(2) * nt(i,j,2) &
                             + 1.0/3.0 * ve * nt(i,j,2)
            else
                fluxe_y(2) = 0
                fluxt_y(2) = 0
            end if
        end if
    end if
    end subroutine

! *** Calculate Metastable Flux ***
    subroutine calc_fluxm(g, i, j, fluxm_x, fluxm_y)
    type(grid), intent(in) :: g
    integer, intent(in) :: i, j
    real(8), intent(out) :: fluxm_x(2), fluxm_y(2)
    
    fluxm_x = 0
    fluxm_y = 0
    
    ! X-dir Fluxes:
    ! - center -
    if (g%type_x(i-1,j-1) == 0) then
        ! Flux at i - 1/2
        fluxm_x(1) = -Dm * (nm(i,j,2) - nm(i-1,j,2)) / g%dx(i-1)
        !call get_fluxm(fluxm_x(1), Dm, g%dx(i-1), nm(i-1:i,j,2))
        
        ! Flux at i + 1/2
        fluxm_x(2) = -Dm * (nm(i+1,j,2) - nm(i,j,2)) / g%dx(i)
        !call get_fluxm(fluxm_x(2), Dm, g%dx(i), nm(i:i+1,j,2))
    
    ! - left -
    else if (g%type_x(i-1,j-1) < 0) then
        ! Flux at i + 1/2
        fluxm_x(2) = -Dm * (nm(i+1,j,2) - nm(i,j,2)) / g%dx(i)
        !call get_fluxm(fluxm_x(2), Dm, g%dx(i), nm(i:i+1,j,2))
        
        ! - electrode -
        if (g%type_x(i-1,j-1) == -2) then
            
            ! Flux at i - 1/2
            fluxm_x(1) = - 0.25 * vi * nm(i,j,2)

        ! - vacuum -
        else if (g%type_x(i-1,j-1) == -1) then
            ! Flux at i - 1/2
            fluxm_x(1) = 0
        end if
    
    ! - right -
    else if (g%type_x(i-1,j-1) > 0) then
        ! Flux at i - 1/2
        fluxm_x(1) = -Dm * (nm(i,j,2) - nm(i-1,j,2)) / g%dx(i-1)
        !call get_fluxm(fluxm_x(1), Dm, g%dx(i-1), nm(i-1:i,j,2))
        
        ! - electrode -
        if (g%type_x(i-1,j-1) == 2) then
            
            ! Flux at i + 1/2
            fluxm_x(2) = 0.25 * vi * nm(i,j,2)
            
        ! - vacuum -
        else if (g%type_x(i-1,j-1) == 1) then
            ! Flux at i + 1/2
            fluxm_x(2) = 0
        end if
    end if
    
    ! Y-dir Fluxes
    if (g%ny > 1) then
        ! - center -
        if (g%type_y(i-1,j-1) == 0) then
            ! Flux at j - 1/2
            fluxm_y(1) = -Dm * (nm(i,j,2) - nm(i,j-1,2)) / g%dy(j-1)
            !call get_fluxm(fluxm_y(1), Dm, g%dy(j-1), nm(i,j-1:j,2))
            
            ! Flux at j + 1/2
            fluxm_y(2) = -Dm * (nm(i,j+1,2) - nm(i,j,2)) / g%dy(j)
            !call get_fluxm(fluxm_y(2), Dm, g%dy(j), nm(i,j:j+1,2))
        
        ! - left -
        else if (g%type_y(i-1,j-1) < 0) then
            ! Flux at j + 1/2
            fluxm_y(2) = -Dm * (nm(i,j+1,2) - nm(i,j,2)) / g%dy(j)
            !call get_fluxm(fluxm_y(2), Dm, g%dy(j), nm(i,j:j+1,2))

            ! Flux at j - 1/2
            fluxm_y(1) = 0
        
        ! - right -
        else if (g%type_y(i-1,j-1) > 0) then
            ! Flux at j - 1/2
            fluxm_y(1) = -Dm * (nm(i,j,2) - nm(i,j-1,2)) / g%dy(j-1)
            !call get_fluxm(fluxm_y(1), Dm, g%dy(j-1), nm(i,j-1:j,2))

            ! Flux at j + 1/2
            fluxm_y(2) = 0
        end if
    end if
    
    end subroutine
    end module
