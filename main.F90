    program main
    use props
    use petsc_lib
    use lapl_lib
    use ptcl_lib
    use sfc_lib
    use circ_lib
    implicit none
    
    type(grid) :: g
    integer :: ts = 0, nx, ny, i, narg, rftype
    real(8) :: l, w, ew, vl, t_fin, t_pr, t_sv, t_sv0, &
               sim_start, sim_fin, time1, time2, Res
    character(80):: arg, path
    
    call PetscInitialize(petsc_null_character, ierr)
    comm = PETSC_COMM_WORLD
    
    call MPI_Comm_rank(comm, my_id, ierr)
    call MPI_Comm_size(comm, nproc, ierr)
    
    call cpu_time(sim_start)
    call cpu_time(time1)
    
    nx = 60
    ny = 60
    px = 1
    py = 1
    l  = 1e-2   / x0
    w  = 1.5e-2 / x0
    ew = 2e-2   / x0
    vl = 500    / phi0
    rftype = 0
    t_fin = 20
    t_pr = 2.77778e-3
    t_sv = 1e-3
    t_sv0 = 1e-3
    res = 1e6
    
    ! Check for -help
    narg = iargc()
    if (mod(narg,2) .ne. 0) then
        if (my_id == 0) then
            write(*,*)
            write(*,*) 'Usage:   mpiexec -n <nproc> ./main <options>'
            write(*,*) 'Options: -nx <nx>, -ny <ny>, -px <px>, -py <py>, -len <len>'
            write(*,*) '         -elen <elen>, -tfin <tfin>, -vl <vl>'
            write(*,*)
        end if
        call MPI_Finalize(ierr)
        stop
    end if
    
    ! Read input arguments
    do i = 1, narg/2
        call getarg(2 * (i - 1) + 1, arg)
        select case (arg)
            case ('-nx')
                call getarg(2 * (i - 1) + 2, arg)
                read(arg,*) nx
            case ('-ny')
                call getarg(2 * (i - 1) + 2, arg)
                read(arg,*) ny
            case ('-len')
                call getarg(2 * (i - 1) + 2, arg)
                read(arg,*) l
            case ('-elen')
                call getarg(2 * (i - 1) + 2, arg)
                read(arg,*) ew
            case ('-rf')
                call getarg(2 * (i - 1) + 2, arg)
                read(arg,*) rftype
            case ('-tfin')
                call getarg(2 * (i - 1) + 2, arg)
                read(arg,*) t_fin
            case ('-px')
                call getarg(2 * (i - 1) + 2, arg)
                read(arg,*) px
            case ('-py')
                call getarg(2 * (i - 1) + 2, arg)
                read(arg,*) py
            case ('-vl')
                call getarg(2 * (i - 1) + 2, arg)
                read(arg,*) vl
                vl = vl / phi0
            case ('-n0')
                call getarg(2 * (i - 1) + 2, arg)
                read(arg,*) n_init
                n_init = n_init * x0**3
            case ('-res')
                call getarg(2 * (i - 1) + 2, arg)
                read(arg,*) res
        end select
    end do
    
    if (px * py .ne. nproc) then
        if (my_id == 0) then
            write(*,*) 'Error: px * py must equal nproc. Stop.'
            write(*,*) '       Enter ./main -help for usage.'
            call MPI_Abort(comm,2,ierr)
        end if
    end if
    if (mod(nx,px) .ne. 0) then
        if (my_id == 0) then
            write(*,*) 'Error: px needs to divide nx. Stop.'
            write(*,*) '       Enter ./main -help for usage.'
            call MPI_Abort(comm,3,ierr)
        end if
    end if
    if (mod(ny,py) .ne. 0) then
        if (my_id == 0) then
            write(*,*) 'Error: py needs to divide ny. Stop.'
            write(*,*) '       Enter ./main -help for usage.'
            call MPI_Abort(comm,4,ierr)
        end if
    end if
    
    call MPI_Barrier(comm, ierr)
    
    write(path,4) int(res / 10**floor(log10(res))), floor(log10(res))
    4 format('output/res_',i0,'e',i0,'/')
    
    call system('mkdir '//trim(path))
    
    call g_init(g, nx, ny, px, py, l, w, ew, trim(path))
    call p_init(g)
    call lapl_init(g, vl, rftype)
    call sfc_init(g)
    call circ_init(vl, res)
    
    do
        ts = ts + 1
        g%t = g%t + g%dt
        if (g%t >= t_fin) exit
        
        ! Solve poisson system
        call lapl_solve(g, Vd, sig)
        
        ! Solve particle system
        call p_step(g, phi_pl)
        
        ! Solve surface charge
        if ((right_wall) .and. (g%ny > 1)) call sfc_step(g, phi_pl)
        
        ! Solve electrode potential
        call circ_step(g, phi_pl)
        
        if ((t_pr <= g%t) .and. (my_id == 0)) then
            call cpu_time(time2)
            write(*,*)
            write(*,11) float(ts), g%t, g%dt, tau_m
            write(*,12) Vd*phi0, Id * e / t0, (time2-time1)/10.0
            t_pr = t_pr + 2.7778e-3
            call cpu_time(time1)
        end if
    
        if (t_sv <= g%t) then
        !if (mod(ts,int(3)) == 0) then
            call savedat(trim(path)//'phi.dat', phi_pl * phi0)
            call savedat(trim(path)//'ni.dat', ni(:,:,1) / x0**3)
            call savedat(trim(path)//'ne.dat', ne(:,:,1) / x0**3)
            call savedat(trim(path)//'nt.dat', nt(:,:,1) * phi0 / x0**3)
            call savedat(trim(path)//'nm.dat', nm(:,:,1) / x0**3)
            
            call MPI_File_Open(comm, trim(path)//'time.dat', &
                MPI_MODE_WRonly + MPI_Mode_Append,  info, fh, ierr)
            if (my_id == 0) call MPI_File_Write(fh, g%t, 1, etype, stat, ierr)
            call MPI_File_Close(fh, ierr)
            
            call MPI_File_Open(comm, trim(path)//'vd.dat', &
                MPI_MODE_WRonly + MPI_Mode_Append,  info, fh, ierr)
            if (my_id == 0) call MPI_File_Write(fh, Vd, 1, etype, stat, ierr)
            call MPI_File_Close(fh, ierr)
            
            call MPI_File_Open(comm, trim(path)//'id.dat', &
                MPI_MODE_WRonly + MPI_Mode_Append,  info, fh, ierr)
            if (my_id == 0) call MPI_File_Write(fh, Id, 1, etype, stat, ierr)
            call MPI_File_Close(fh, ierr)
            
            t_sv  = t_sv + t_sv0
            t_sv0 = t_sv0 * 1.01
        end if
    end do
    
    if (my_id == 0) then
        call cpu_time(sim_fin)
        write(*,*)
        write(*,9) int(sim_fin - sim_start) / 3600, &
                    mod(int(sim_fin - sim_start)/60,60)
        write(*,*)
    end if
    
    call petsc_destroy
    call PetscFinalize(ierr)

11 format('Timestep:', es10.2, '  Time:', es10.2, '  dT:', es10.2, '  Tm:', es10.2)
12 format('  Vd: ', f8.2, '  Id: ', es10.2, '  time/100ns:', f7.2, ' hr')   
9  format('Simulation finished in ', i0, ' hr ', i0, ' min')

    end program
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
