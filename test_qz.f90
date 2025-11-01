! qz_dggev3.f90
!
! Program to perform the QZ algorithm on a matrix pencil (A, B)
! using the LAPACK routine DGGEV3.
! The generalized eigenvalues lambda are computed as:
! lambda = (ALPHAR + i * ALPHAI) / BETA
!
program test_qz
    use matrix_module
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
    write(*, '("Matrix A:")')
    call print_matrix(A)

    write(*, '("Matrix B:")')
    call print_matrix(B)

    ! --- 3. Call DGGEV3 (Generalized Eigenvalue Problem) ---
    ! JOBVL='N', JOBVR='N': Do not compute eigenvectors, only eigenvalues.
    call DGGEV3( &
        'V', 'V', & ! JOBVL, JOBVR (No eigenvectors)
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
        'V', 'V', & ! JOBVL, JOBVR (No eigenvectors)
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
        write(*, '("Index |    ALPHAR    |    ALPHAI    |    BETA    |   LAMBDA ")')
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

        call print_matrix(A)
        call print_matrix(B)
        call print_matrix(VL)
        call print_matrix(VR)

    else if (INFO < 0) then
        write(*, '("Error: Argument ", I0, " in DGGEV3 had an illegal value.")') -INFO
    else
        write(*, '("Error: QZ iteration failed to converge. INFO = ", I0)') INFO
    end if

    write(*, '(" ")')
    write(*, '("--- Results from DGGEV3_RQ ---")')

    if (INFO2 == 0) then
        write(*, '("---------------------------------------------------------------")')
        write(*, '("Index |    ALPHAR    |    ALPHAI    |    BETA    |   LAMBDA ")')
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

        call print_matrix(A2)
        call print_matrix(B2)
        call print_matrix(VL2)
        call print_matrix(VR2)
        
    else if (INFO2 < 0) then
        write(*, '("Error: Argument ", I0, " in DGGEV3 had an illegal value.")') -INFO2
    else
        write(*, '("Error: QZ iteration failed to converge. INFO = ", I0)') INFO2
    end if

end program test_qz