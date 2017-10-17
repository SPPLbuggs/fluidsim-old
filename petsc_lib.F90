    module petsc_lib
#include <petsc/finclude/petscksp.h>
    use petscksp
    use props
    implicit none
    
    ! Petsc Variables
    Mat A
    Vec b, x
    KSP ksp
    PetscInt Istart, Iend, ii, jj, nn(1)
    
    abstract interface
        subroutine subIn(g, i, j, n, m, f, b_temp)
            use props
            type(grid), intent(in) :: g
            integer, intent(in) :: i, j, n, m
            real(8), intent(in) :: f(n,m)
            real(8), intent(out) :: b_temp(1)
        end subroutine
    end interface
    
    contains
    
! *** Create PETSc Objects ***
    subroutine petsc_create(g)
    type(grid), intent(in) :: g

    ! Petsc Objects A and b
    call MatCreate(comm, A, ierr)
    call MatSetSizes(A, g%nloc, g%nloc, g%nglob, g%nglob, ierr)
    call MatSetUp(A, ierr)
    call MatSetFromOptions(A, ierr)
    call MatSeqAIJSetPreallocation(A, 5, petsc_null_integer, ierr)
    call MatSetOption(A, mat_ignore_zero_entries, petsc_true, ierr)

    ! Find parallel partitioning range
    call MatGetOwnershipRange(A, Istart, Iend, ierr)

    ! Create parallel vectors
    call VecCreateMPI(comm, g%nloc, g%nglob, b, ierr)
    call VecSetFromOptions(b, ierr)
    call VecSetOption(b, vec_ignore_negative_indices, petsc_true, ierr)
    
    call VecCreateMPI(comm, g%nloc, g%nglob, x, ierr)
    call VecSetFromOptions(x, ierr)
    call VecSetOption(x, vec_ignore_negative_indices, petsc_true, ierr)
    
    ! Create Linear Solver
    call KSPCreate(comm, ksp, ierr)
    call KSPSetOperators(ksp, A, A, ierr)
    call KSPSetInitialGuessNonzero(ksp, PETSC_TRUE, ierr)
    call KSPSetTYpe(ksp, KSPIBCGS, ierr)
    call KSPSetFromOptions(ksp, ierr)
    
    end subroutine

! *** Assemble A and b ***
    subroutine assem_Ab(g, f, feval)
    type(grid), intent(in) :: g
    real(8), intent(inout) :: f(:,:)
    procedure(subIn) :: feval
    real(8):: b_temp(1), A_temp(5,1)
    integer :: i, j, cols(5), rows(1)
    
    ! Assemble A and b
    do j = 2, g%by+1
        do i = 2, g%bx+1
            cols = -1
            rows = g%node(i,j)
            
            call feval(g, i, j, g%bx+2, g%by+2, f, b_temp)
            call jacob(g, i, j, cols, f, feval, b_temp, A_temp)
            
            ii = 1
            jj = 5
            call VecSetValues(b, ii, rows, -b_temp, Insert_Values, ierr)
            call MatSetValues(A, ii, rows, jj, cols, A_temp, &
                              Insert_Values, ierr)
      end do
    end do
    
    call VecAssemblyBegin(b, ierr)
    call VecAssemblyEnd(b, ierr)

    call MatAssemblyBegin(A, Mat_Final_Assembly, ierr)
    call MatAssemblyEnd(A, Mat_Final_Assembly, ierr)
    
    end subroutine
    
! *** Assemble b ***
    subroutine assem_b(g, f, feval)
    type(grid), intent(in) :: g
    real(8), intent(inout) :: f(:,:)
    procedure(subIn) :: feval
    real(8) :: b_temp(1)
    integer :: i, j, rows(1)
    
    ! Assemble A and b
    do j = 2, g%by+1
        do i = 2, g%bx+1
            rows = g%node(i,j)
            
            call feval(g, i, j, g%bx+2, g%by+2, f, b_temp)
            
            ii = 1
            call VecSetValues(b, ii, rows, -b_temp, Insert_Values, ierr)
      end do
    end do
    
    call VecAssemblyBegin(b, ierr)
    call VecAssemblyEnd(b, ierr)
    end subroutine
    
! *** Numerical Jacobian ***
    subroutine jacob(g, i_loc, j_loc, cols, f, feval, b_temp, A_temp)
    type(grid), intent(in) :: g
    integer, intent(in):: i_loc, j_loc
    integer, intent(inout):: cols(:)
    procedure(subIn) :: feval
    real(8), intent(in):: b_temp(1)
    real(8), intent(inout):: f(:,:), A_temp(:,:)
    real(8) :: perturb, temp, b_pert(1)
    integer :: i,j,k, width, k_start, k_stop
    integer, dimension(5,2):: stencil

    ! initialize
    temp = 0
    width = 0
    perturb = 1e-4
    b_pert = 0
        
    k_start = 1
    k_stop  = 3
    
    if (g%ny > 1) k_stop = 5
    
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
        i = i_loc + stencil(k,1)
        j = j_loc + stencil(k,2)

        temp = f(i,j)
        f(i,j) = f(i,j) + perturb
            
        call feval(g, i_loc, j_loc, g%bx+2, g%by+2, f, b_pert)
        
        cols(k) = g%node(i,j)
        A_temp(k, :) = (b_pert - b_temp) / perturb
        
        f(i,j) = temp
    end do
    end subroutine
    
! *** Update Solution ***
    subroutine upd_soln(g, f)
    type(grid), intent(in) :: g
    real(8), intent(inout) :: f(:,:)
    integer :: i,j
    real(8) :: soln(1)
    
    do j = 2, g%by+1
        do i = 2, g%bx+1
            nn = g%node(i,j)
            call VecGetValues(x, 1, nn, soln, ierr)        
            f(i,j) = f(i,j) + soln(1)
        end do
    end do
    
    call comm_real(g%bx, g%by, f)
    end subroutine
    
! *** Destroy PETSc Objects ***
    subroutine petsc_destroy
    call KSPDestroy(ksp,ierr)
    call VecDestroy(b,ierr)
    call VecDestroy(x,ierr)
    call MatDestroy(A,ierr)
    end subroutine

! *** View Matrix and Vector ***
    subroutine view
    integer :: wait
    call MatView(A,PETSC_VIEWER_STDOUT_WORLD,ierr)
    call VecView(b,PETSC_VIEWER_STDOUT_WORLD,ierr)
    if (my_id == 0) read(*,*) wait
    call MPI_Barrier(comm, ierr)
    end subroutine
    
    end module
