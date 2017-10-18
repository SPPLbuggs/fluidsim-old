    module ptcl_props
    use props
    implicit none
    
    ! Particle density initialization
    real(8) :: n_init = 1e11 * x0**3
    
    ! electron properties
    real(8), parameter:: me = 9.10938188e-31

    ! positive ion properties
    real(8), parameter:: mi  = 1.67262178d-27 * 39.948, &
                         mui = 1.45d3 / p * 1d-4 * phi0 * t0 / x0**2, &
                         Ti  = Tg, &
                         vi  = Sqrt( (8d0 * kb*Ti) / (pi * mi) ) * t0 / x0, &
                         Di  = mui * (kb*Ti / e) / phi0

    ! metastable argon properties
    real(8), parameter:: mm   = mi, &
                         Dm   = 2.42d18 * x0 / ninf / 1d2 * t0, &
                         k_r  = 2d-7 / 1d6 / x0**3 * t0, &
                         k_mp = 6.2d-10 / 1.0d6 / x0**3 * t0, &
                         k_2q = 3d-15 / 1d6 / x0**3 * t0, &
                         k_3q = 1.1d-31 / 1.0d12 / x0**6 * t0

    ! reactions
    real(8), parameter:: H_ir  =  15.8d0 / phi0, &
                         H_ex  =  11.56d0 / phi0, &
                         H_si  =   4.14d0 / phi0, &
                         H_sc  = -11.56d0 / phi0, &
                         beta  =   1d-13 * t0 / x0**3, &
                         gam   =   0.075

    contains

! *** Calculate Particle Flux ***
!     (scharfetter-Gummel)
    subroutine get_flux(flux, E, dh, q, mu, D, n_left, n_right)
    real(8), intent(inout) :: flux
    integer, intent(in) :: q
    real(8), intent(in) :: E, dh, mu, D, n_left, n_right
    real(8) :: v, tol, arg
    
    tol = 1e-12
    v = q * mu * E   
    arg = v * dh / D
    
    if (abs(q*E) < tol) then
        flux = D * (n_left - n_right) / dh
    
    ! Positive exponentials blow up,
    !  so rewrite always as negative exp.
    !  this is the same analytical expression
    else if (arg > 0) then
        flux = v * (n_left - n_right * exp(max(-arg,-100d0))) &
            / (1d0 - exp(max(-arg, -100d0)))
    else
        flux = v * (n_right - n_left * exp(max(arg, -100d0))) &
            / (1d0  - exp( max(arg,-100d0)))
    end if
    end subroutine

! *** Ion flux ***
    subroutine get_fluxi(flux, E, mu, D, dh, n)
    real(8), intent(out) :: flux
    real(8), intent(in)  :: E, mu, D, dh, n(2)
    
    flux = 0.5 * (n(1) + n(2)) * ( mu * E - D * log(n(2)/n(1)) / dh )
    end subroutine
    
! *** Electron flux ***
    subroutine get_fluxe(flux, E, mu, dh, n, T)
    real(8), intent(out) :: flux
    real(8), intent(in)  :: E, mu, dh, n(2), T(2)
    
    flux = 0.5 * (n(1) + n(2)) * mu * ( -E &
              - 1.0 / 3.0 * (T(1) + T(2)) * log(n(2)/n(1)) / dh &
              - 2.0 / 3.0 * (T(2) - T(1)) / dh )
    end subroutine
    
! *** Electron energy flux ***
    subroutine get_fluxt(flux, Je, mu, dh, n, T)
    real(8), intent(out) :: flux
    real(8), intent(in)  :: Je, mu, dh, n(2), T(2)
    
    flux = 5.0/6.0 * Je * (T(1) + T(2)) &
         - 1.0/3.0 * mu * (n(2) + n(1)) * (T(2) - T(1)) / dh
    end subroutine

! *** Metastable flux ***
    subroutine get_fluxm(flux, D, dh, n)
    real(8), intent(out) :: flux
    real(8), intent(in)  :: D, dh, n(2)
    
    flux = -0.5 * (n(1) + n(2)) * D * log(n(2)/n(1)) / dh
    end subroutine
    
! *** Runge-Kutta Adaptive Timestepping ***
    subroutine rk_step(g, n, k, stage, dt, nerr_n, order, n_min)
    type(grid), intent(in) :: g
    integer, intent(in) :: stage, order
    real(8), intent(in) :: dt, k(:,:,:), n_min
    real(8), intent(inout) :: n(:,:,:), nerr_n
    real(8) :: err_n(g%bx+2, g%by+2), abs_tol = 1e-4, rel_tol = 1e-4
    integer :: i,j
    
    if (order == 1) then
    ! euler scheme
        do j = 2, g%by+1
            do i = 2, g%bx+1
                n(i,j,1) = n(i,j,3) + k(i,j,1) * dt
            end do
        end do
        
        call comm_real(g%bx, g%by, n(:,:,1))
        
    else if (order == 2) then
    ! 2nd order
        if (stage == 1) then
            do j = 2, g%by+1
                do i = 2, g%bx+1
                    n(i,j,1) = n(i,j,3) + k(i,j,1) * dt / 2.0
                    n(i,j,2) = n(i,j,3) + k(i,j,1) * dt
                end do
            end do
            
        else
            do j = 2, g%by+1
                do i = 2, g%bx+1
                    n(i,j,1) = n(i,j,1) + k(i,j,2) * dt / 2d0
                end do
            end do
        end if
        
        n(:,:,1:2) = max(n(:,:,1:2), n_min)
        
        call comm_real(g%bx, g%by, n(:,:,1))
        call comm_real(g%bx, g%by, n(:,:,2))
        
    else if (order == 4) then
        ! merson 4("5") adaptive time-stepping
        if (stage == 1) then
            do j = 2, g%by+1
                do i = 2, g%bx+1
                    n(i,j,2) = n(i,j,3) + k(i,j,1) * dt / 3.0
                end do
            end do
            
            n(:,:,2) = max(n(:,:,2), n_min)
            call comm_real(g%bx, g%by, n(:,:,2))
            
        else if (stage == 2) then
            do j = 2, g%by+1
                do i = 2, g%bx+1
                    n(i,j,2) = n(i,j,3) + dt * ( k(i,j,1) + k(i,j,2) ) / 6.0
                 end do
            end do
               
            n(:,:,2) = max(n(:,:,2), n_min)
            
            call comm_real(g%bx, g%by, n(:,:,2))
            
        else if (stage == 3) then
            do j = 2, g%by+1
                do i = 2, g%bx+1
                    n(i,j,2) = n(i,j,3) + dt * (k(i,j,1) + k(i,j,3) * 3.0) / 8.0
                end do
            end do
            
            n(:,:,2) = max(n(:,:,2), n_min)
            call comm_real(g%bx, g%by, n(:,:,2))
            
        else if (stage == 4) then
            do j = 2, g%by+1
                do i = 2, g%bx+1
                    n(i,j,2) = n(i,j,3) + dt * (k(i,j,1) &
                               - k(i,j,3) * 3.0 + k(i,j,4) * 4.0 ) / 2.0
                end do
            end do
            
            n(:,:,2) = max(n(:,:,2), n_min)
            call comm_real(g%bx, g%by, n(:,:,2))
        else
            do j = 2, g%by+1
                do i = 2, g%bx+1
                    n(i,j,1) = n(i,j,3) + dt * (k(i,j,1) &
                               + k(i,j,4) * 4.0 + k(i,j,5)) / 6.0
                end do
            end do
            
            n(:,:,1) = max(n(:,:,1), n_min)
            call comm_real(g%bx, g%by, n(:,:,1))
            
            err_n = 0
            
            err_n = abs(dt * (k(:,:,1) * 2.0 / 30.0 - k(:,:,3) * 3.0 / 10.0 &
                        + k(:,:,4) * 4.0 / 15.0 - k(:,:,5) / 30.0 ))
            
            nerr_n = maxval(err_n/(abs_tol+rel_tol*abs(n(:,:,3))))
            
            call MPI_Allreduce(MPI_In_Place, nerr_n, 1, etype, &
                               MPI_Max, comm, ierr)
        end if
    end if
    end subroutine
           
    function get_Te(nt,ne)
    real(8):: get_Te
    real(8), intent(in) :: nt, ne
    
    if ((nt >= 0d0) .and. (ne >= 0d0)) then
        get_Te = nt / ne
    else
        write(*,*) "Error, negative density. Stop."
        stop
        get_Te = 1d-8
    end if
    
    return
    end function
    
    function get_mue(T)
    real(8):: get_mue
    real(8), intent(in):: T
    real(8):: x, &
              a = 55.0126564455,     &
              b =  0.685594575551,   &
              c = -0.383637563328,   &
              d = -0.0340209762821,  &
              f =  0.0276748887423,  &
              g =  0.00242420301108, &
              h = -0.0012183881312,  &
              i = -5.6471677382e-05, &
              j =  2.11510269874e-05
    
    x = log(min(2d2, max(1d-2, T * phi0)))
    
    get_mue = exp(a + b*x + c*x**2. + d*x**3. + f*x**4. + g*x**5. &
                  + h*x**6. + i*x**7. + j*x**8.) * x0 / ninf * t0 * phi0
    return
    end function get_mue
    
    function get_k_ex(T)
    real(8):: get_k_ex
    real(8), intent(in):: T
    real(8):: x, &
              a = -50.3785102239,    &
              b =  19.0129183764,    &
              c = -11.7950315424,    &
              d =   7.41674013553,   &
              f =  -3.84148086698,   &
              g =   1.2962229976,    &
              h =  -0.259359346989,  &
              i =   0.0279182131315, &
              j =  -0.00124438710099
              
    x = log(min(2d2, max(5d-1, T * phi0)))
    
    get_k_ex = exp(a + b*x + c*x**2. + d*x**3. + f*x**4. + g*x**5. &
                   + h*x**6. + i*x**7. + j*x**8.) * t0 / x0**3
    return
    end function get_k_ex

    function get_k_ir(T)
    real(8):: get_k_ir
    real(8), intent(in):: T
    real(8):: x, &
              a = -56.2478139553,    &
              b =  26.6123052468,    &
              c = -15.9868576469,    &
              d =   7.97316507041,   &
              f =  -3.18287109994,   &
              g =   0.890666459851,  &
              h =  -0.156771317788,  &
              i =   0.0153819279555, &
              j =  -0.000638729430911
              
    x = log(min(2d2, max(5d-1, T * phi0)))

    get_k_ir = exp(a + b*x + c*x**2. + d*x**3. + f*x**4. + g*x**5. &
                   + h*x**6. + i*x**7. + j*x**8.) * t0 / x0**3
    return
    end function get_k_ir

    function get_nu(T)
    real(8):: get_nu
    real(8), intent(in):: T
    real(8):: x, &
              a = -32.3199262825,      &
              b =   1.69388044254,     &
              c =   0.0842378277404,   &
              d =  -0.164556371047,    &
              f =   0.00861307304011,  &
              g =   0.00855257716132,  &
              h =  -0.000983504760785, &
              i =  -0.000160952834008, &
              j =   2.37965210684e-05
              
    x = log(min(2d2, max(1d-2, T * phi0)))
    
    get_nu = exp(a + b*x + c*x**2. + d*x**3. + f*x**4. + g*x**5. &
                 + h*x**6. + i*x**7 + j*x**8.) / x0**3 * ninf * t0
    return
    end function get_nu

    function get_k_sc(T)
    real(8):: get_k_sc
    real(8), intent(in):: T
    real(8):: x, &
              a = -21.4827864151,   &
              b =   0.457356923276, &
              c =  -0.555439231606, &
              d =   1.27257798891,  &
              f =  -0.67840685073,  &
              g =   0.10591014464
              
    x = log(min(16d0, max(5d-1, T * phi0)))
    
    get_k_sc = exp(a + b*x + c*x**2. + d*x**3. + f*x**4. &
                   + g*x**5. ) / 1.0d6 / x0**3 * t0
    return
    end function get_k_sc

    function get_k_si(T)
    real(8):: get_k_si
    real(8), intent(in):: T
    real(8):: x, &
              a = -43.1347385848, &
              b =  43.9905424566, &
              c = -28.1169537586, &
              d =   8.28853856817, &
              f =  -0.931626144207
              
    x = log(min(16d0, max(5d-1, T * phi0)))
    
    get_k_si = exp(a + b*x + c*x**2. + d*x**3. + f*x**4.) &
                  / 1.0d6 / x0**3 * t0
    return
    end function get_k_si

    end module
