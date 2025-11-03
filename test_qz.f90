program test_qz
    use matrix_module
    implicit none
    external DGGEV3

    ! --- Parameters ---
    integer, parameter :: N = 3      ! Dimension of the square matrices
    integer, parameter :: LWORK_SIZE = 8*N ! Recommended minimum workspace size

    real(kind=8) :: A(N, N)    
    real(kind=8) :: A2(N, N)   
    real(kind=8) :: B(N, N)    
    real(kind=8) :: B2(N, N)   
    
    real(kind=8) :: ALPHAR(N)
    real(kind=8) :: ALPHAR2(N)
    real(kind=8) :: ALPHAI(N)
    real(kind=8) :: ALPHAI2(N)
    real(kind=8) :: BETA(N)
    real(kind=8) :: BETA2(N)
    
    real(kind=8) :: WORK(LWORK_SIZE) 
    real(kind=8) :: WORK2(LWORK_SIZE)
    real(kind=8) :: VL(N, N)         
    real(kind=8) :: VL2(N, N)        
    real(kind=8) :: VR(N, N)         
    real(kind=8) :: VR2(N, N)        
    integer :: LDA = N               
    integer :: LDB = N               
    integer :: LDVL = N              
    integer :: LDVR = N              
    integer :: INFO                  
    integer :: INFO2                 
    
    integer :: i, j
    real(kind=8) :: lambda_real, lambda_imag 
    real(kind=8) :: lambda_real2, lambda_imag2 

    ! --- 1. Initialize Matrices A and B ---
    
    A = reshape([1.0, 1.0, 5.0, & 
                 1.0, 2.0, 6.0, & 
                 1.0, 1.0, 3.0], & 
                 shape=[N, N])

    B = reshape([2.0, 4.0, 3.0, & 
                 1.0, 1.0, 2.0, &
                 1.0, 1.0, 1.0], & 
                 shape=[N, N])

    A2 = reshape([1.0, 1.0, 5.0, &
                 1.0, 2.0, 6.0, &
                 1.0, 1.0, 3.0], &
                 shape=[N, N])
                 
    B2 = reshape([2.0, 4.0, 3.0, &
                 1.0, 1.0, 2.0, & 
                 1.0, 1.0, 1.0], & 
                 shape=[N, N])


    ! --- 2. Print Input Matrices ---
    write(*, '("Matrix A:")')
    call print_matrix(A)

    write(*, '("Matrix B:")')
    call print_matrix(B)

    ! --- 3. Call DGGEV3 (Generalized Eigenvalue Problem) ---
    call DGGEV3( &
        'V', 'V', & ! JOBVL, JOBVR 
        N, &        ! N (Matrix size)
        A, LDA, &   ! A, LDA
        B, LDB, &   ! B, LDB
        ALPHAR, ALPHAI, BETA, & ! Outputs: ALPHAR, ALPHAI, BETA
        VL, LDVL, & ! VL, LDVL 
        VR, LDVR, & ! VR, LDVR
        WORK, LWORK_SIZE, & ! WORK, LWORK
        INFO &
    )

    ! Call own implementation with copies of the input arguments
    call DGGEV3_RQ( &
        'V', 'V', & ! JOBVL, JOBVR 
        N, &        ! N (Matrix size)
        A2, LDA, &   ! A, LDA
        B2, LDB, &   ! B, LDB
        ALPHAR2, ALPHAI2, BETA2, & ! Outputs: ALPHAR, ALPHAI, BETA
        VL2, LDVL, & ! VL, LDVL
        VR2, LDVR, & ! VR, LDVR
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
        
    else if (INFO2 < 0) then
        write(*, '("Error: Argument ", I0, " in DGGEV3 had an illegal value.")') -INFO2
    else
        write(*, '("Error: QZ iteration failed to converge. INFO = ", I0)') INFO2
    end if

    print *, "A:"
    call print_matrix(A)
    print *, "B:"
    call print_matrix(B)
    print *, "VL:"
    call print_matrix(VL)
    print *, "VR:"
    call print_matrix(VR)

    write(*, '("---------------------------------------------------------------")')
    write(*, '(" ")')

    print *, "A2:"
    call print_matrix(A2)
    print *, "B2:"
    call print_matrix(B2)
    print *, "VL2:"
    call print_matrix(VL2)
    print *, "VR2:"
    call print_matrix(VR2)

end program test_qz