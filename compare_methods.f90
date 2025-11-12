program compare_methods
    use matrix_module
    use quicksort
    use, intrinsic :: ieee_arithmetic
    implicit none
    external DGGEV3

    integer :: N, LWORK_SIZE, LDA, LDB, LDVL, LDVR, INFO
    integer, parameter :: N_MAX = 10, METHODS = 2
    integer :: i, j, k

    real(kind=8), allocatable :: tau(:)

    real(kind=8), allocatable :: A(:,:), B(:,:), D(:,:), P(:,:), Q(:,:)
    real(kind=8), allocatable :: ALPHAR(:), ALPHAI(:), BETA(:), WORK(:)
    real(kind=8) :: lambda_real, lambda_imag, norm(METHODS)
    complex(kind=8), allocatable :: lambda(:)
    real(kind=8), allocatable :: VL(:,:), VR(:,:)
    complex(kind=8), allocatable :: EIGS(:)

    

    ! Initialize Random Seed
    integer :: n_size, seed_id = 123456789
    integer, allocatable :: seed_array(:)
    call RANDOM_SEED(SIZE=n_size)
    allocate(seed_array(n_size))
    seed_array = seed_id
    call RANDOM_SEED(PUT=seed_array)

    open(unit=10, file="output.txt", status="unknown")

    do N = 1, N_MAX

        ! --- Parameters ---
        LWORK_SIZE = 8*N ! Recommended minimum workspace size

        allocate( A(N,N), B(N,N), D(N,N), P(N,N), Q(N,N) )
        allocate( VL(N,N), VR(N,N) )
        allocate( ALPHAR(N), ALPHAI(N), BETA(N), EIGS(N) )
        allocate( lambda(N) )
        allocate( WORK(LWORK_SIZE) )    
        LDA = N             
        LDB = N               
        LDVL = N              
        LDVR = N              
        
        allocate( tau(N) )
        
        ! --- 1. Initialize Matrices A and B ---
        ! Initialize P and Q first
        ! A = PDQ, B = PQ, with P and Q orthogonal matrices
        D = 0.0_8
        do i = 1, N
            EIGS(i) = i
            D(i,i) = i
        end do
        CALL RANDOM_NUMBER(P)
        CALL RANDOM_NUMBER(Q)
        do i=1,n
            do j=1,n
                P(i,j) = P(i,j) - 0.5d0
                Q(i,j) = Q(i,j) - 0.5d0
            end do
        end do
        ! Use QR to generate orthogonal matrices
        CALL DGEQRF(N,N,P,N,tau,WORK,LWORK_SIZE,INFO)
        CALL DORGQR(N,N,N,P,N,tau,WORK,LWORK_SIZE,INFO)

        call DGEQRF(n, n, Q, n, tau, WORK, LWORK_SIZE, INFO)
        call DORGQR(n, n, n, Q, n, tau, WORK, LWORK_SIZE, INFO)

        do k=1,METHODS
            ! A = P*D
            do i=1,n
                do j=1,n
                    A(i,j) = sum(P(i,:)*D(:,j))
                end do
            end do
            A = matmul(A,Q)
            ! B = P*Q
            B = matmul(P,Q)

            if (K == 1) then
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
            else 
                call DGGEV3_RQ( &
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
            end if

            if (INFO == 0) then
                do i = 1, N
                    if (BETA(i) .eq. 0.0_8) then
                        ! Infinite eigenvalue (or near-infinite)
                        lambda(i) = ieee_value(1.0_8, IEEE_QUIET_NAN)
                    else
                        ! Finite eigenvalue: lambda = (ALPHAR + i * ALPHAI) / BETA
                        lambda_real = ALPHAR(i) / BETA(i)
                        lambda_imag = ALPHAI(i) / BETA(i)
                        lambda(i) = cmplx(lambda_real, lambda_imag, kind=8)
                    end if
                end do

                ! Sort eigenvalues
                call quicksort_complex(lambda,1,N)
            else if (INFO < 0) then
                write(*, '("Error: Argument ", I0, " in DGGEV3 had an illegal value.")') -INFO
            else
                write(*, '("Error: QZ iteration failed to converge. INFO = ", I0)') INFO
            end if

            norm(k) = norm_diff(EIGS, lambda)

            ! do i = 1, N
            !     print *, lambda(i)
            ! end do

            ! print *, "Norm of Error", norm

        end do

        ! N qr rq
        write(10,*) N, norm

        deallocate( A, B, D, P, Q, VL, VR, ALPHAR, ALPHAI, BETA, EIGS, WORK, tau, lambda )

    end do

    close(10)

end program compare_methods