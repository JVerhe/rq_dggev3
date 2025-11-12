module quicksort

    implicit none
    private
    public :: quicksort_complex

contains

    recursive subroutine quicksort_complex(a, left, right)
        implicit none
        complex(kind=8), intent(inout) :: a(:)
        integer, intent(in) :: left, right
        integer :: i, j
        complex(kind=8) :: pivot, temp

        i = left
        j = right
        pivot = a((left + right) / 2)

        do
            do while (abs(a(i)) < abs(pivot))
                i = i + 1
            end do
            do while (abs(a(j)) > abs(pivot))
                j = j - 1
            end do

            if (i <= j) then
                temp = a(i)
                a(i) = a(j)
                a(j) = temp
                i = i + 1
                j = j - 1
            end if

            if (i > j) exit
        end do
        if (left < j) call quicksort_complex(a, left, j)
        if (i < right) call quicksort_complex(a, i, right)
    end subroutine quicksort_complex

end module quicksort