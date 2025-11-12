program compare_methods
    use matrix_module
    use quicksort
    use, intrinsic :: ieee_arithmetic
    implicit none
    external DGGEV3

    integer :: N, LWORK_SIZE, LDA, LDB, LDVL, LDVR, INFO
    integer, parameter :: N_MAX = 2500, METHODS = 2
    integer :: i, j, k, l

    real(kind=8), allocatable :: tau(:)

    real(kind=8), allocatable :: A(:,:), B(:,:), D(:,:), P(:,:), Q(:,:)
    real(kind=8), allocatable :: ALPHAR(:), ALPHAI(:), BETA(:)
    real(kind=8) :: lambda_real, lambda_imag, norm(METHODS)
    complex(kind=8), allocatable :: lambda(:)
    real(kind=8), allocatable :: VL(:,:), VR(:,:)
    complex(kind=8), allocatable :: EIGS(:)

    real(kind=8), allocatable :: WORK_QR(:)    ! workspace for QR (DGEQRF/DORGQR)
    real(kind=8), allocatable :: WORK_GE(:)    ! workspace for DGGEV3
    real(kind=8), allocatable :: work_query(:)
    integer :: LWORK_QR, LWORK_GE

    ! Initialize Random Seed
    integer :: n_size, seed_id = 123456789
    integer, allocatable :: seed_array(:)
    call RANDOM_SEED(SIZE=n_size)
    allocate(seed_array(n_size))
    seed_array = seed_id
    call RANDOM_SEED(PUT=seed_array)

    open(unit=10, file="raw.txt", status="unknown")

    do l = 2, 15
        N = FLOOR(1.8**l)
        print *,N
        if (N .gt. N_MAX) then
            exit
        end if

        ! --- Parameters ---
        LWORK_GE = 8*N ! Recommended minimum workspace size
        LWORK_QR = N

        allocate( A(N,N), B(N,N), D(N,N), P(N,N), Q(N,N) )
        allocate( VL(N,N), VR(N,N) )
        allocate( ALPHAR(N), ALPHAI(N), BETA(N), EIGS(N) )
        allocate( lambda(N) )
        allocate( WORK_QR(LWORK_QR) ) 

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
        ! Transform range to [-0.5, 0.5]
        P = P - 0.5d0
        Q = Q - 0.5d0


        ! Use QR to generate orthogonal matrices
        CALL DGEQRF(N,N,P,N,tau,WORK_QR,LWORK_QR,INFO)
        if (INFO /= 0) then
            write(*,*) "DGEQRF(P) failed, INFO=", INFO; stop
        end if
        CALL DORGQR(N,N,N,P,N,tau,WORK_QR,LWORK_QR,INFO)
        if (INFO /= 0) then
            write(*,*) "DORGQR(P) failed, INFO=", INFO; stop
        end if

        call DGEQRF(n, n, Q, n, tau, WORK_QR, LWORK_QR, INFO)
        if (INFO /= 0) then
            write(*,*) "DGEQRF(Q) failed, INFO=", INFO; stop
        end if
        call DORGQR(n, n, n, Q, n, tau, WORK_QR, LWORK_QR, INFO)
        if (INFO /= 0) then
            write(*,*) "DORGQR(Q) failed, INFO=", INFO; stop
        end if

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
                allocate(work_query(1))
                work_query(1) = 0.0d0
                LWORK_GE = -1

                call DGGEV3('N','N', N, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, work_query, LWORK_GE, INFO)
                if (INFO /= 0) then
                    write(*,*) 'Workspace query for DGGEV3 returned INFO=', INFO
                    stop
                end if

                LWORK_GE = int(work_query(1))
                deallocate(work_query)
                allocate(WORK_GE(LWORK_GE))

                call DGGEV3( &
                    'V', 'V', & ! JOBVL, JOBVR 
                    N, &        ! N (Matrix size)
                    A, LDA, &   ! A, LDA
                    B, LDB, &   ! B, LDB
                    ALPHAR, ALPHAI, BETA, & ! Outputs: ALPHAR, ALPHAI, BETA
                    VL, LDVL, & ! VL, LDVL 
                    VR, LDVR, & ! VR, LDVR
                    WORK_GE, LWORK_GE, & ! WORK, LWORK
                    INFO &
                )
                if (INFO /= 0) then
                    write(*,*) 'DGGEV3 returned INFO=', INFO
                end if
            else 
                allocate(work_query(1))
                work_query(1) = 0.0d0
                LWORK_GE = -1

                call DGGEV3('N','N', N, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, work_query, LWORK_GE, INFO)
                if (INFO /= 0) then
                    write(*,*) 'Workspace query for DGGEV3 returned INFO=', INFO
                    stop
                end if

                LWORK_GE = int(work_query(1))
                deallocate(work_query)
                allocate(WORK_GE(LWORK_GE))

                call DGGEV3_RQ( &
                    'V', 'V', & ! JOBVL, JOBVR 
                    N, &        ! N (Matrix size)
                    A, LDA, &   ! A, LDA
                    B, LDB, &   ! B, LDB
                    ALPHAR, ALPHAI, BETA, & ! Outputs: ALPHAR, ALPHAI, BETA
                    VL, LDVL, & ! VL, LDVL 
                    VR, LDVR, & ! VR, LDVR
                    WORK_GE, LWORK_GE, & ! WORK, LWORK
                    INFO &
                )
            end if

            deallocate( WORK_GE )

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

            ! norm(1) = qr_error, norm(2) = rq_error in iteration N
            norm(k) = norm_diff(EIGS, lambda)
            
        end do

        ! N qr rq
        write(10,*) N, norm
        deallocate( A, B, D, P, Q, WORK_QR, VL, VR, ALPHAR, ALPHAI, BETA, EIGS, tau, lambda )
    end do

    close(10)

end program compare_methods