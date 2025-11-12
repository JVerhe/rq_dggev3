module matrix_module
  ! This module contains procedures for matrix operations.
  implicit none

  ! Define a kind parameter for double precision (8-byte real/complex)
  INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(15, 307)
  
  private
  
  ! --- Generic Interface Declaration ---
  ! The name 'print_matrix' will resolve to either the real or complex specific routine
  ! based on the type of the argument passed.
  public :: print_matrix, norm_diff
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

end module matrix_module