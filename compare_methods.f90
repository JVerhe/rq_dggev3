program compare_methods
    use matrix_module
    use quicksort
    use, intrinsic :: ieee_arithmetic
    implicit none

    integer :: N, LDA, LDB, LDVL, LDVR, INFO
    integer, parameter :: N_MAX = 700, METHODS = 2
    integer :: i, j, k, l

    character(len=30) :: timestamp
    character(len=50) :: fname
    integer :: year, month, day, hour, minute
    integer :: date_vals(8)

    real(kind=8), allocatable :: tau(:)

    real(kind=8), allocatable :: A(:,:), B(:,:)
    real(kind=8), allocatable :: ALPHAR(:), ALPHAI(:), BETA(:)
    real(kind=8) :: lambda_real, lambda_imag, norm(METHODS)
    complex(kind=8), allocatable :: lambda(:)
    real(kind=8), allocatable :: VL(:,:), VR(:,:)
    complex(kind=8), allocatable :: EIGS(:)

    real(kind=8), allocatable :: WORK_GE(:)    ! workspace for DGGEV3
    real(kind=8), allocatable :: work_query(:)
    integer :: LWORK_GE

    ! Initialize Random Seed
    integer :: n_size, seed_id = 123456789
    integer, allocatable :: seed_array(:)
    call RANDOM_SEED(SIZE=n_size)
    allocate(seed_array(n_size))
    seed_array = seed_id
    call RANDOM_SEED(PUT=seed_array)

    ! Create new file with current timestamp
    call date_and_time(values=date_vals)
    year   = date_vals(1)
    month  = date_vals(2)
    day    = date_vals(3)
    hour   = date_vals(5)
    minute = date_vals(6)
    write(timestamp,'(I4.4,"_",I2.2,"_",I2.2,"_",I2.2,"-",I2.2)') &
        year, month, day, hour, minute
    fname = 'results/' // trim(timestamp) // "_output.txt"
    open(unit=10, file=fname, status="unknown", action="write")

    do l = 2, 15
        N = FLOOR(1.8**l)
        
        if (N .gt. N_MAX) then
            print *, "Writing result to ", fname
            exit
        end if
        print *,N

        ! --- Parameters ---
        LWORK_GE = 8*N ! Recommended minimum workspace size

        allocate( A(N,N), B(N,N) )
        allocate( VL(N,N), VR(N,N) )
        allocate( ALPHAR(N), ALPHAI(N), BETA(N), EIGS(N) )
        allocate( lambda(N) )

        LDA = N             
        LDB = N               
        LDVL = N              
        LDVR = N              
        
        ! --- 1. Initialize Matrices A and B ---
        ! Initialize P and Q first
        ! A = PDQ, B = PQ, with P and Q orthogonal matrices
        do k=1,METHODS
            
            call generate_reg_pencil(N,A,B,EIGS)

            if (K == 1) then
                allocate(work_query(1))
                work_query(1) = 0.0d0
                LWORK_GE = -1

                call DGGEV3_QR('N','N', N, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, work_query, LWORK_GE, INFO)
                if (INFO /= 0) then
                    write(*,*) 'Workspace query for DGGEV3 returned INFO=', INFO
                    stop
                end if

                LWORK_GE = int(work_query(1))
                deallocate(work_query)
                allocate(WORK_GE(LWORK_GE))

                call DGGEV3_QR( &
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

                call DGGEV3_RQ('N','N', N, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, work_query, LWORK_GE, INFO)
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
        deallocate( A, B, VL, VR, ALPHAR, ALPHAI, BETA, EIGS, lambda )
    end do

    close(10)

end program compare_methods