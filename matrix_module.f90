module matrix_module
  ! This module contains procedures for matrix operations.
  implicit none
  private
  public :: print_matrix

contains

  subroutine print_matrix(mat)
    implicit none

    real(kind=8), dimension(:, :), intent(in) :: mat

    integer :: i, j
    integer :: num_rows, num_cols

    num_rows = size(mat, 1)
    num_cols = size(mat, 2)

    print *, "(" // &
             char(ichar('0') + num_rows) // "x" // char(ichar('0') + num_cols) // ")"

    do i = 1, num_rows
      write(*, '(*(F12.4))') (mat(i, j), j = 1, num_cols)
    end do

  end subroutine print_matrix

end module matrix_module
