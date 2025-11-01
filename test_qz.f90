! qz_dggev3.f90
!
! Program to perform the QZ algorithm on a matrix pencil (A, B)
! using the LAPACK routine DGGEV3.
! The generalized eigenvalues lambda are computed as:
! lambda = (ALPHAR + i * ALPHAI) / BETA
!
program test_qz
    implicit none

    ! External declaration for the LAPACK routine
    external DGGEV3

    ! --- Parameters ---
    integer, parameter :: N = 3      ! Dimension of the square matrices
    integer, parameter :: LWORK_SIZE = 8*N ! Recommended minimum workspace size

    ! --- Matrix and Vector Declarations ---
    ! Note: Fortran stores arrays in column-major order.
    real(kind=8) :: A(N, N)     ! Coefficient matrix A
    real(kind=8) :: A2(N, N)     ! Coefficient matrix A
    real(kind=8) :: B(N, N)     ! Coefficient matrix B (the "pencil")
    real(kind=8) :: B2(N, N)     ! Coefficient matrix B (the "pencil")
    
    ! --- DGGEV3 Output Arrays for Generalized Eigenvalues ---
    ! ALPHAR: Real parts of alpha components (numerator of lambda)
    ! ALPHAI: Imaginary parts of alpha components (numerator of lambda)
    ! BETA: Denominator of lambda. BETA is always real and non-negative.
    real(kind=8) :: ALPHAR(N)
    real(kind=8) :: ALPHAR2(N)
    real(kind=8) :: ALPHAI(N)
    real(kind=8) :: ALPHAI2(N)
    real(kind=8) :: BETA(N)
    real(kind=8) :: BETA2(N)
    
    ! --- Workspace and Status Variables ---
    real(kind=8) :: WORK(LWORK_SIZE) ! Workspace array
    real(kind=8) :: WORK2(LWORK_SIZE) ! Workspace array
    real(kind=8) :: VL(N, N)         ! Left eigenvectors (not used here)
    real(kind=8) :: VL2(N, N)         ! Left eigenvectors (not used here)
    real(kind=8) :: VR(N, N)         ! Right eigenvectors (not used here)
    real(kind=8) :: VR2(N, N)         ! Right eigenvectors (not used here)
    integer :: LDA = N               ! Leading dimension of A
    integer :: LDB = N               ! Leading dimension of B
    integer :: LDVL = N              ! Leading dimension of VL
    integer :: LDVR = N              ! Leading dimension of VR
    integer :: INFO                  ! Status code
    integer :: INFO2                  ! Status code
    
    ! --- Loop index and temporary variables ---
    integer :: i, j
    ! Variables moved from inside the DO loop to resolve the "Unexpected data declaration" error.
    real(kind=8) :: lambda_real, lambda_imag 
    real(kind=8) :: lambda_real2, lambda_imag2 

    ! --- 1. Initialize Matrices A and B ---
    ! Example: A simple diagonal pencil (A-lambda*B).
    ! A is upper triangular with eigenvalues 1, 2, 3. B is identity.
    ! Expected generalized eigenvalues (lambda): 1.0, 2.0, 3.0
    
    A = reshape([1.0, 4.0, 5.0, & ! Column 1
                 1.0, 2.0, 6.0, & ! Column 2
                 1.0, 1.0, 3.0], & ! Column 3
                 shape=[N, N])

    B = reshape([1.0, 4.0, 3.0, & ! Column 1
                 1.0, 1.0, 2.0, & ! Column 2
                 1.0, 1.0, 1.0], & ! Column 3
                 shape=[N, N])

    A2 = reshape([1.0, 4.0, 5.0, & ! Column 1
                 1.0, 2.0, 6.0, & ! Column 2
                 1.0, 1.0, 3.0], & ! Column 3
                 shape=[N, N])
                 
    B2 = reshape([1.0, 4.0, 3.0, & ! Column 1
                 1.0, 1.0, 2.0, & ! Column 2
                 1.0, 1.0, 1.0], & ! Column 3
                 shape=[N, N])

    ! --- 2. Print Input Matrices ---
    write(*, '("--- Input Matrix Pencil (A, B) ---")')
    write(*, '("Matrix A:")')
    do i = 1, N
        write(*, '(3F10.4)') A(i, 1:N)
    end do

    write(*, '("Matrix B:")')
    do i = 1, N
        write(*, '(3F10.4)') B(i, 1:N)
    end do

    ! --- 3. Call DGGEV3 (Generalized Eigenvalue Problem) ---
    ! JOBVL='N', JOBVR='N': Do not compute eigenvectors, only eigenvalues.
    call DGGEV3( &
        'N', 'N', & ! JOBVL, JOBVR (No eigenvectors)
        N, &        ! N (Matrix size)
        A, LDA, &   ! A, LDA
        B, LDB, &   ! B, LDB
        ALPHAR, ALPHAI, BETA, & ! Outputs: ALPHAR, ALPHAI, BETA
        VL, LDVL, & ! VL, LDVL (ignored since JOBVL='N')
        VR, LDVR, & ! VR, LDVR (ignored since JOBVR='N')
        WORK, LWORK_SIZE, & ! WORK, LWORK
        INFO &
    )

    ! Make copy of input arguments: A B ALPHAR ALPHAI BETA VL VR WORK INFO
    call DGGEV3_RQ( &
        'N', 'N', & ! JOBVL, JOBVR (No eigenvectors)
        N, &        ! N (Matrix size)
        A2, LDA, &   ! A, LDA
        B2, LDB, &   ! B, LDB
        ALPHAR2, ALPHAI2, BETA2, & ! Outputs: ALPHAR, ALPHAI, BETA
        VL2, LDVL, & ! VL, LDVL (ignored since JOBVL='N')
        VR2, LDVR, & ! VR, LDVR (ignored since JOBVR='N')
        WORK2, LWORK_SIZE, & ! WORK, LWORK
        INFO2 &
    )

    ! --- 4. Check INFO and Print Results ---
    write(*, '(" ")')
    write(*, '("--- Results from DGGEV3 ---")')

    if (INFO == 0) then
        write(*, '("---------------------------------------------------------------")')
        write(*, '("Index |   Alpha Real |   Alpha Imag |    Beta    |   Lambda (Approx.) ")')
        write(*, '("---------------------------------------------------------------")')
        
        do i = 1, N
            ! Calculate the approximate lambda value for printing clarity
            
            if (BETA(i) .eq. 0.0_8) then
                ! Infinite eigenvalue (or near-infinite)
                write(*, '(I5, " | ", F12.6, " | ", F12.6, " | ", F10.6, " | INFINITE")') &
                     i, ALPHAR(i), ALPHAI(i), BETA(i)
            else
                ! Finite eigenvalue: lambda = (ALPHAR + i * ALPHAI) / BETA
                lambda_real = ALPHAR(i) / BETA(i)
                lambda_imag = ALPHAI(i) / BETA(i)
                
                write(*, '(I5, " | ", F12.6, " | ", F12.6, " | ", F10.6, " | ", F10.6, " + ", F10.6, "i")') &
                     i, ALPHAR(i), ALPHAI(i), BETA(i), lambda_real, lambda_imag
            end if
        end do
        
        write(*, '("---------------------------------------------------------------")')

        write(*, '(" ")')
        write(*, '("* Generalized Eigenvectors (VL, VR)")')
        
        i = 1
        do while (i <= N)
            
            write(*, '(" ")')
            write(*, '("Eigenvector Pair ", I0, ": Lambda = (", F6.3, " + ", F6.3, "i) / ", F6.3, ")")') &
                 i, ALPHAR(i), ALPHAI(i), BETA(i)

            if (ALPHAI(i) .eq. 0.0_8) then
                ! Real Eigenvalue: Eigenvector is in column i
                write(*, '("  Right Eigenvector (VR):")')
                do j = 1, N
                    write(*, '("    ", I2, ": ", F10.6, " + 0.000000i")') j, VR(j, i)
                end do
                write(*, '("  Left Eigenvector (VL):")')
                do j = 1, N
                    write(*, '("    ", I2, ": ", F10.6, " + 0.000000i")') j, VL(j, i)
                end do
                i = i + 1
            else
                ! Complex Eigenvalue Pair: Eigenvectors are in columns i and i+1
                ! First eigenvector (lambda) = VR(:, i) + i * VR(:, i+1)
                write(*, '("  Right Eigenvector (VR) - First of Conjugate Pair:")')
                do j = 1, N
                    write(*, '("    ", I2, ": ", F10.6, " + ", F10.6, "i")') VR(j, i), VR(j, i+1)
                end do
                ! Second eigenvector (lambda_conjugate) = VR(:, i) - i * VR(:, i+1)
                write(*, '("  Right Eigenvector (VR) - Second of Conjugate Pair:")')
                do j = 1, N
                    write(*, '("    ", I2, ": ", F10.6, " - ", F10.6, "i")') VR(j, i), VR(j, i+1)
                end do
                ! Repeat for Left Eigenvector (VL)
                write(*, '("  Left Eigenvector (VL) - First of Conjugate Pair:")')
                do j = 1, N
                    write(*, '("    ", I2, ": ", F10.6, " + ", F10.6, "i")') VL(j, i), VL(j, i+1)
                end do
                write(*, '("  Left Eigenvector (VL) - Second of Conjugate Pair:")')
                do j = 1, N
                    write(*, '("    ", I2, ": ", F10.6, " - ", F10.6, "i")') VL(j, i), VL(j, i+1)
                end do
                i = i + 2 ! Skip the next column as it was handled in this pair
            end if
        end do

    else if (INFO < 0) then
        write(*, '("Error: Argument ", I0, " in DGGEV3 had an illegal value.")') -INFO
    else
        write(*, '("Error: QZ iteration failed to converge. INFO = ", I0)') INFO
    end if

    write(*, '(" ")')
    write(*, '("--- Results from DGGEV3_RQ ---")')

    if (INFO2 == 0) then
        write(*, '("---------------------------------------------------------------")')
        write(*, '("Index |   Alpha Real |   Alpha Imag |    Beta    |   Lambda (Approx.) ")')
        write(*, '("---------------------------------------------------------------")')
        
        do i = 1, N
            ! Calculate the approximate lambda value for printing clarity
            
            if (BETA2(i) .eq. 0.0_8) then
                ! Infinite eigenvalue (or near-infinite)
                write(*, '(I5, " | ", F12.6, " | ", F12.6, " | ", F10.6, " | INFINITE")') &
                     i, ALPHAR2(i), ALPHAI2(i), BETA2(i)
            else
                ! Finite eigenvalue: lambda = (ALPHAR + i * ALPHAI) / BETA
                lambda_real2 = ALPHAR2(i) / BETA2(i)
                lambda_imag2 = ALPHAI2(i) / BETA2(i)
                
                write(*, '(I5, " | ", F12.6, " | ", F12.6, " | ", F10.6, " | ", F10.6, " + ", F10.6, "i")') &
                     i, ALPHAR2(i), ALPHAI2(i), BETA2(i), lambda_real2, lambda_imag2
            end if
        end do
        
        write(*, '("---------------------------------------------------------------")')

        write(*, '(" ")')
        write(*, '("* Generalized Eigenvectors (VL, VR)")')
        
        i = 1
        do while (i <= N)
            
            write(*, '(" ")')
            write(*, '("Eigenvector Pair ", I0, ": Lambda = (", F6.3, " + ", F6.3, "i) / ", F6.3, ")")') &
                 i, ALPHAR2(i), ALPHAI2(i), BETA2(i)

            if (ALPHAI2(i) .eq. 0.0_8) then
                ! Real Eigenvalue: Eigenvector is in column i
                write(*, '("  Right Eigenvector (VR):")')
                do j = 1, N
                    write(*, '("    ", I2, ": ", F10.6, " + 0.000000i")') j, VR2(j, i)
                end do
                write(*, '("  Left Eigenvector (VL):")')
                do j = 1, N
                    write(*, '("    ", I2, ": ", F10.6, " + 0.000000i")') j, VL2(j, i)
                end do
                i = i + 1
            else
                ! Complex Eigenvalue Pair: Eigenvectors are in columns i and i+1
                ! First eigenvector (lambda) = VR(:, i) + i * VR(:, i+1)
                write(*, '("  Right Eigenvector (VR) - First of Conjugate Pair:")')
                do j = 1, N
                    write(*, '("    ", I2, ": ", F10.6, " + ", F10.6, "i")') VR2(j, i), VR2(j, i+1)
                end do
                ! Second eigenvector (lambda_conjugate) = VR(:, i) - i * VR(:, i+1)
                write(*, '("  Right Eigenvector (VR) - Second of Conjugate Pair:")')
                do j = 1, N
                    write(*, '("    ", I2, ": ", F10.6, " - ", F10.6, "i")') VR2(j, i), VR2(j, i+1)
                end do
                ! Repeat for Left Eigenvector (VL)
                write(*, '("  Left Eigenvector (VL) - First of Conjugate Pair:")')
                do j = 1, N
                    write(*, '("    ", I2, ": ", F10.6, " + ", F10.6, "i")') VL2(j, i), VL2(j, i+1)
                end do
                write(*, '("  Left Eigenvector (VL) - Second of Conjugate Pair:")')
                do j = 1, N
                    write(*, '("    ", I2, ": ", F10.6, " - ", F10.6, "i")') VL2(j, i), VL2(j, i+1)
                end do
                i = i + 2 ! Skip the next column as it was handled in this pair
            end if
        end do

    else if (INFO2 < 0) then
        write(*, '("Error: Argument ", I0, " in DGGEV3 had an illegal value.")') -INFO2
    else
        write(*, '("Error: QZ iteration failed to converge. INFO = ", I0)') INFO2
    end if

end program test_qz