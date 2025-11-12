module matrix_module
  ! This module contains procedures for matrix operations.
  implicit none

  ! Define a kind parameter for double precision (8-byte real/complex)
  INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(15, 307)
  
  private
  
  ! --- Generic Interface Declaration ---
  ! The name 'print_matrix' will resolve to either the real or complex specific routine
  ! based on the type of the argument passed.
  public :: print_matrix, norm_diff, generate_pencil
  INTERFACE print_matrix
    MODULE PROCEDURE print_real_matrix, print_complex_matrix
  END INTERFACE
  
contains

  ! ==========================================================
  ! Specific Subroutine for REAL (Double Precision) Matrices
  ! ==========================================================
  subroutine print_real_matrix(mat)
    implicit none

    REAL(KIND=DP), DIMENSION(:, :), INTENT(IN) :: mat

    INTEGER :: i, j
    INTEGER :: num_rows, num_cols

    num_rows = SIZE(mat, 1)
    num_cols = SIZE(mat, 2)

    ! Using a specific format (I0) to prevent spacing around the integers
    WRITE(*, '(A, I0, A, I0, A)') 'Real Matrix (', num_rows, 'x', num_cols, ')'
    
    ! Prints rows with F12.4 format for clarity and alignment
    DO i = 1, num_rows
      WRITE(*, '(*(F12.4))') (mat(i, j), j = 1, num_cols)
    END DO

  end subroutine print_real_matrix

  ! ==========================================================
  ! Specific Subroutine for COMPLEX (Double Precision) Matrices
  ! ==========================================================
  subroutine print_complex_matrix(mat)
    implicit none

    COMPLEX(KIND=DP), DIMENSION(:, :), INTENT(IN) :: mat

    INTEGER :: i, j
    INTEGER :: num_rows, num_cols

    num_rows = SIZE(mat, 1)
    num_cols = SIZE(mat, 2)

    ! Using a specific format (I0) to prevent spacing around the integers
    WRITE(*, '(A, I0, A, I0, A)') 'Complex Matrix (', num_rows, 'x', num_cols, ')'
    
    ! For complex numbers, the output format must be wide enough 
    ! to accommodate (real, imaginary). G format (General) handles this gracefully.
    ! G24.12 gives 24 characters of space, which is typically sufficient 
    ! for a (r, i) double-precision complex number.
    DO i = 1, num_rows
      WRITE(*, '(*(G24.12))') (mat(i, j), j = 1, num_cols)
    END DO

  end subroutine print_complex_matrix

  function norm_diff(A, B) result(norm)
        implicit none
        complex(8), intent(in) :: A(:), B(:)
        real(8) :: norm

        if (size(A) /= size(B)) then
            print *, "Error: sizes of A and B differ"
            stop
        end if

        norm = sqrt(sum(abs(A - B)**2))
  end function norm_diff

  ! ==========================================================
  ! Generate matrix pencil (A, B) with specified eigenvalues
  ! A = P * D * Q,  B = P * Q
  ! P, Q are random orthogonal matrices (via QR)
  ! ==========================================================
  subroutine generate_pencil(N, A, B, EIGS)
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none

    integer, intent(in) :: N
    real(real64), allocatable, intent(out) :: A(:,:), B(:,:)
    complex(real64), allocatable, intent(out) :: EIGS(:)

    real(real64), allocatable :: D(:,:), P(:,:), Q(:,:), tau(:), work(:)
    integer :: i, j, info, lwork

    ! --- allocate internal matrices ---
    allocate(A(N,N), B(N,N))
    allocate(D(N,N), P(N,N), Q(N,N))
    allocate(EIGS(N))
    allocate(tau(N))

    ! workspace for QR
    lwork = max(1,N)
    allocate(work(lwork))

    ! diagonal matrix with eigenvalues 1..N
    D = 0.0_real64
    do i = 1, N
        D(i,i) = real(i,kind=real64)
        EIGS(i) = cmplx(real(i,kind=real64),0.0_real64,kind=real64)
    end do

    ! random P,Q in [-0.5, 0.5]
    call random_number(P);  call random_number(Q)
    P = P - 0.5_real64
    Q = Q - 0.5_real64

    ! QR -> orthogonal P
    call dgeqrf(N,N,P,N,tau,work,lwork,info); if(info/=0) stop "dgeqrf(P) failed"
    call dorgqr(N,N,N,P,N,tau,work,lwork,info); if(info/=0) stop "dorgqr(P) failed"

    ! QR -> orthogonal Q
    call dgeqrf(N,N,Q,N,tau,work,lwork,info); if(info/=0) stop "dgeqrf(Q) failed"
    call dorgqr(N,N,N,Q,N,tau,work,lwork,info); if(info/=0) stop "dorgqr(Q) failed"


    ! A = P * D * Q
    A = matmul(P, matmul(D, Q))

    ! B = P * Q
    B = matmul(P, Q)

    ! free workspace
    deallocate(D, P, Q, tau, work)

  end subroutine generate_pencil

end module matrix_module