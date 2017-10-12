    module props
#include <petsc/finclude/petscksp.h>
    use petscksp
    implicit none
    
    ! mpi variables
    integer :: comm, my_id, nproc, ierr, &
               rx, ry, px, py, north, south, east, west, &
               fh, etype, amode, info, stat(MPI_Status_Size), core_array, glob_array
    integer(kind=MPI_Offset_Kind) :: dispx, dispy
    
    type :: grid
        integer :: nx, ny, bx, by, offx, offy, nglob, nloc, dof
        integer, allocatable :: type_x(:,:), type_y(:,:), node(:,:,:)
        real(8) :: dt, t, l, w, ew
        real(8), allocatable :: dx(:), dlx(:), dy(:), dly(:), r(:)
    end type
    
    ! non-dimensional parameters
    real(8), parameter:: x0   = 1e-3,    &
                         phi0 = 1e3,     &
                         t0 = 1e-6

    ! fundamental constants
    real(8), parameter:: pi     = 4d0 * atan(1d0),        &
                         eps0   = 8.85418782e-12,         &
                         mu0    = 4 * pi * 1e-7,          &
                         c0     = 1d0 / sqrt(eps0 * mu0), &
                         c      = c0 / x0 * t0,           &
                         e      = 1.60217646e-19,         &
                         kb     = 1.3806503e-23

    ! case properties
    real(8), parameter:: Tg     = 350,                            & ! kelvin
                         p      = 3,                              & ! torr
                         ninf   = p * 101325d0 / 760d0 / kb / Tg * x0**3, &
                         n_init = 1e12 * x0**3, &
                         n_zero = 1e8 * x0**3
    
    real(8) :: phiL = 0, phiR = 0
    logical, parameter:: cyl = .True.

    contains
    
    ! *** Initialize Grid ***
    subroutine g_init(g, nx, ny, px, py, l, w, ew)
    type(grid), intent(inout) :: g
    real(8), intent(in) :: l, w, ew
    integer, intent(in) :: nx, ny, px, py
    integer :: i, j, d
    real(8) :: xtemp, ytemp
    real(8), allocatable :: x(:), y(:)
    
    g%nx = nx
    g%ny = ny
    
    ! Coordinates of processor (my_id = ry * px + rx)
    rx = mod(my_id, px)
    ry = my_id / px
    
    ! Determine neighbor processors
    north = (ry - 1) * px + rx
    if (ry - 1 < 0) north = MPI_Proc_Null
    south = (ry + 1) * px + rx
    if (ry + 1 >= py) south = MPI_Proc_Null
    west = ry * px + rx - 1
    if (rx - 1 < 0) west = MPI_Proc_Null
    east = ry * px + rx + 1
    if (rx + 1 >= px) east = MPI_Proc_Null
    
    ! Domain decomposition
    g%bx = nx / px
    g%by = ny / py
    g%offx = rx * g%bx
    g%offy = ry * g%by
    
    ! Define node types
    allocate(g%type_x(g%bx,g%by), g%type_y(g%bx,g%by))
    g%type_x = 0
    g%type_y = 0
    
    g%l   = l
    g%ew  = ew
    g%w   = w
    
    xtemp = 2.5 / float(g%nx+1)
    x = (/ ( tanh(-1.25 + xtemp * (g%offx + i - 1)), i = 1, g%bx+2) /)
    !x = (/ ( g%l / float(g%bx+1) * (i - 1), i = 1, g%bx+2) /)
    
    xtemp = x(1)
    call MPI_Bcast( xtemp, 1, MPI_Real8, 0, comm, ierr)
    x = x - xtemp
    
    xtemp = x(g%bx+2)
    call MPI_Bcast( xtemp, 1, MPI_Real8, nproc-1, comm, ierr)
    x = x / xtemp
    x = x * g%l
    
    g%dx = (/ ( x(i+1) - x(i), i = 1, g%bx+1) /)
    g%dlx = (/ (0.5 * (g%dx(i+1) + g%dx(i)), i = 1, g%bx) /)
    
    if (g%ny > 1) then
        ytemp = 1.25 / float(g%ny+1)
        y = (/ ( tanh(-1.25 + ytemp * (g%offy + j - 1)), j = 1, g%by+2) /)
        !y = (/ ( g%w / float(g%by+1) * (j - 1), j = 1, g%by+2) /)
    
        ytemp = y(1)
        call MPI_Bcast( ytemp, 1, MPI_Real8, 0, comm, ierr)
        y = y - ytemp
    
        ytemp = y(g%by+2)
        call MPI_Bcast( ytemp, 1, MPI_Real8, nproc-1, comm, ierr)
        y = y / ytemp
        y = y * g%w
    
        g%dy = (/ ( y(j+1) - y(j), j = 1, g%by+1) /)
        g%dly = (/ (0.5 * (g%dy(j+1) + g%dy(j)), j = 1, g%by) /)
    else
        allocate(g%dy(1))
        g%dy = g%dx(1)
        ytemp = g%offy * g%dy(1)
        y = (/ (ytemp + g%dy * (j - 1), j = 1, g%by) /)
    end if
    
    g%r = y
    
    do j = 1, g%by
        do i = 1, g%bx      
            if ((ry == 0) .and. (j == 1)) g%type_y(i,j) = -1
            if ((ry == py-1) .and. (j == g%by)) g%type_y(i,j) =  1
            
            if ((rx == 0) .and. (i == 1)) then
                if (y(j) .le. g%ew) then
                    g%type_x(i,j) = -2
                else
                    g%type_x(i,j) = -1
                end if
            end if
            
            if ((rx == px-1) .and. (i == g%bx)) then
                if (y(j) .le. g%ew) then
                    g%type_x(i,j) = 2
                else
                    g%type_x(i,j) = 1
                end if
            end if
        end do
    end do
    
    g%t = 0
    g%dt = 1e-6
    
    g%dof   = 1
    g%nloc  = g%bx * g%by * g%dof
    g%nglob = g%nx * g%ny * g%dof
    
    allocate( g%node(g%bx+2, g%by+2, g%dof) )
    g%node = 0
    do d = 1, g%dof
        do j = 2, g%by+1
            do i = 2, g%bx+1
                g%node(i,j,d) = g%nloc * my_id + (d - 1) &
                            + ((i - 2) + (j - 2) * g%bx) * g%dof
            end do
        end do
        
        call comm_int(g%bx, g%by, g%node(:,:,d))
    end do
    
    ! MPI-IO Variables
    amode = MPI_Mode_WRonly + MPI_Mode_Create + MPI_Mode_EXCL
    etype = MPI_Real8
    stat  = MPI_Status_Ignore
    info  = MPI_Info_Null
    dispx  = rx*g%bx*8
    dispy  = ry*g%by*8
    
    ! Save mesh to disk
    call MPI_File_Open(comm, 'output/meshx.dat', amode,  info, fh, ierr)
    if (ierr .ne. MPI_SUCCESS) then
        if (my_id == 0) call MPI_File_Delete('output/meshx.dat', info, ierr);
        call MPI_File_Open(comm, 'output/meshx.dat', amode,  info, fh, ierr)
    end if
    
    call MPI_File_Set_View(fh, dispx, etype, etype, 'native', info, ierr)
    if (ry == 0) call MPI_File_Write(fh, x, g%bx, etype, stat, ierr)
    call MPI_File_Close(fh, ierr)
    
    call MPI_File_Open(comm, 'output/meshy.dat', amode,  info, fh, ierr)
    if (ierr .ne. MPI_SUCCESS) then
        if (my_id == 0) call MPI_File_Delete('output/meshy.dat', info, ierr);
        call MPI_File_Open(comm, 'output/meshy.dat', amode,  info, fh, ierr)
    end if
    
    call MPI_File_Set_View(fh, dispy, etype, etype, 'native', info, ierr)
    if (rx == 0) call MPI_File_Write(fh, y, g%by, etype, stat, ierr)
    call MPI_File_Close(fh, ierr)
    
    call MPI_File_Open(comm, 'output/time.dat', amode,  info, fh, ierr)
    if (ierr .ne. MPI_SUCCESS) then
        if (my_id == 0) call MPI_File_Delete('output/time.dat', info, ierr);
        call MPI_File_Open(comm, 'output/time.dat', amode,  info, fh, ierr)
    end if
    call MPI_File_Close(fh, ierr)
    
    call MPI_File_Open(comm, 'output/phi.dat', amode,  info, fh, ierr)
    if (ierr .ne. MPI_SUCCESS) then
        if (my_id == 0) call MPI_File_Delete('output/phi.dat', info, ierr);
        call MPI_File_Open(comm, 'output/phi.dat', amode,  info, fh, ierr)
    end if
    call MPI_File_Close(fh, ierr)
    
    call MPI_File_Open(comm, 'output/ni.dat', amode,  info, fh, ierr)
    if (ierr .ne. MPI_SUCCESS) then
        if (my_id == 0) call MPI_File_Delete('output/ni.dat', info, ierr);
        call MPI_File_Open(comm, 'output/ni.dat', amode,  info, fh, ierr)
    end if
    call MPI_File_Close(fh, ierr)
    
    call MPI_File_Open(comm, 'output/ne.dat', amode,  info, fh, ierr)
    if (ierr .ne. MPI_SUCCESS) then
        if (my_id == 0) call MPI_File_Delete('output/ne.dat', info, ierr);
        call MPI_File_Open(comm, 'output/ne.dat', amode,  info, fh, ierr)
    end if
    call MPI_File_Close(fh, ierr)
    
    call MPI_File_Open(comm, 'output/nt.dat', amode,  info, fh, ierr)
    if (ierr .ne. MPI_SUCCESS) then
        if (my_id == 0) call MPI_File_Delete('output/nt.dat', info, ierr);
        call MPI_File_Open(comm, 'output/nt.dat', amode,  info, fh, ierr)
    end if
    call MPI_File_Close(fh, ierr)
    
    call MPI_File_Open(comm, 'output/nm.dat', amode,  info, fh, ierr)
    if (ierr .ne. MPI_SUCCESS) then
        if (my_id == 0) call MPI_File_Delete('output/nm.dat', info, ierr);
        call MPI_File_Open(comm, 'output/nm.dat', amode,  info, fh, ierr)
    end if
    call MPI_File_Close(fh, ierr)
    
    call MPI_Type_Create_Subarray(2, (/ g%bx+2, g%by+2 /), (/ g%bx, g%by /), &
        (/ 1, 1 /), MPI_Order_Fortran, etype, core_array, ierr)
    call MPI_Type_Commit(core_array, ierr)
    
    call MPI_Type_Create_Subarray(2, (/ g%nx, g%ny /), (/ g%bx, g%by /), &
        (/ rx * g%bx, ry * g%by /), MPI_Order_Fortran, etype, glob_array, ierr)
    call MPI_Type_Commit(glob_array, ierr)
    
    call MPI_Barrier(comm, ierr)
    
    end subroutine

! *** Communicate data across processors
    subroutine comm_real(bx, by, A)
    integer, intent(in) :: bx, by
    real(8), intent(inout) :: A(:,:)
    
    call MPI_Send(A(2:bx+1,2), bx, etype, north, 9, comm, ierr)
    call MPI_Send(A(2:bx+1,by+1), bx, etype, south, 9, comm, ierr)
    call MPI_Recv(A(2:bx+1,1), bx, etype, north, 9, comm, stat, ierr)
    call MPI_Recv(A(2:bx+1,by+2), bx, etype, south, 9, comm, stat, ierr)
    
    call MPI_Send(A(bx+1,2:by+1), by, etype, east, 9, comm, ierr)
    call MPI_Send(A(2,2:by+1), by, etype, west, 9, comm, ierr)
    call MPI_Recv(A(bx+2,2:by+1), by, etype, east, 9, comm, stat, ierr)
    call MPI_Recv(A(1,2:by+1), by, etype, west, 9, comm, stat, ierr)
    
    end subroutine
    
! *** Communicate data across processors
    subroutine comm_int(bx, by, A)
    integer, intent(in) :: bx, by
    integer, intent(inout) :: A(:,:)
    
    call MPI_Send(A(2:bx+1,2), bx, MPI_INT, north, 9, comm, ierr)
    call MPI_Send(A(2:bx+1,by+1), bx, MPI_INT, south, 9, comm, ierr)
    call MPI_Recv(A(2:bx+1,1), bx, MPI_INT, north, 9, comm, stat, ierr)
    call MPI_Recv(A(2:bx+1,by+2), bx, MPI_INT, south, 9, comm, stat, ierr)
    
    call MPI_Send(A(bx+1,2:by+1), by, MPI_INT, east, 9, comm, ierr)
    call MPI_Send(A(2,2:by+1), by, MPI_INT, west, 9, comm, ierr)
    call MPI_Recv(A(bx+2,2:by+1), by, MPI_INT, east, 9, comm, stat, ierr)
    call MPI_Recv(A(1,2:by+1), by, MPI_INT, west, 9, comm, stat, ierr)
    
    end subroutine

! *** Save Data ***
    subroutine savedat(path,dat)
    character(*), intent(in) :: path
    real(8), intent(in) :: dat(:,:)
    integer (kind = MPI_Offset_Kind) :: offset
    
    call MPI_File_Open(comm, path, MPI_MODE_RDWR + MPI_MODE_APPEND, info, fh, ierr)
    call MPI_File_Get_Position(fh, offset, ierr)
    call MPI_File_Set_View(fh, offset, etype, glob_array, 'native', info, ierr)
    call MPI_File_Write_All(fh, dat, 1, core_array, stat, ierr)
    call MPI_File_Close(fh, ierr)
    end subroutine
    
    end module
