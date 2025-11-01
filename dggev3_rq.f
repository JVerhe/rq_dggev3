*  Definition:
*  ===========
*
*       SUBROUTINE DGGEV3_rq( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR,
*      $                   ALPHAI, BETA, VL, LDVL, VR, LDVR, WORK, LWORK,
*      $                   INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          JOBVL, JOBVR
*       INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * ), ALPHAI( * ), ALPHAR( * ),
*      $                   B( LDB, * ), BETA( * ), VL( LDVL, * ),
*      $                   VR( LDVR, * ), WORK( * )
*       ..
*  =====================================================================

      SUBROUTINE dggev3_rq( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR,
     $                   ALPHAI, BETA, VL, LDVL, VR, LDVR, WORK, LWORK,
     $                   INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          JOBVL, JOBVR
      INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), ALPHAI( * ), ALPHAR( * ),
     $                   B( LDB, * ), BETA( * ), VL( LDVL, * ),
     $                   vr( ldvr, * ), work( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0d+0, one = 1.0d+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ILASCL, ILBSCL, ILV, ILVL, ILVR, LQUERY
      CHARACTER          CHTEMP
      INTEGER            ICOLS, IERR, IHI, IJOBVL, IJOBVR, ILEFT, ILO,
     $                   in, iright, irows, itau, iwrk, jc, jr, lwkopt,
     $                   lwkmin
      DOUBLE PRECISION   ANRM, ANRMTO, BIGNUM, BNRM, BNRMTO, EPS,
     $                   SMLNUM, TEMP
*     ..
*     .. Local Arrays ..
      LOGICAL            LDUMMA( 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           dgeqrf, dggbak, dggbal,
     $                   dgghd3, dlaqz0, dlacpy,
     $                   dlascl, dlaset, dorgqr,
     $                   dormqr, dtgevc, xerbla
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANGE
      EXTERNAL           lsame, dlamch, dlange
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, max, sqrt
*     ..
*     .. Executable Statements ..
*
*     Decode the input arguments
*
      IF( lsame( jobvl, 'N' ) ) THEN
         ijobvl = 1
         ilvl = .false.
      ELSE IF( lsame( jobvl, 'V' ) ) THEN
         ijobvl = 2
         ilvl = .true.
      ELSE
         ijobvl = -1
         ilvl = .false.
      END IF
*
      IF( lsame( jobvr, 'N' ) ) THEN
         ijobvr = 1
         ilvr = .false.
      ELSE IF( lsame( jobvr, 'V' ) ) THEN
         ijobvr = 2
         ilvr = .true.
      ELSE
         ijobvr = -1
         ilvr = .false.
      END IF
      ilv = ilvl .OR. ilvr
*
*     Test the input arguments
*
      info = 0
      lquery = ( lwork.EQ.-1 )
      lwkmin = max( 1, 8*n )
      IF( ijobvl.LE.0 ) THEN
         info = -1
      ELSE IF( ijobvr.LE.0 ) THEN
         info = -2
      ELSE IF( n.LT.0 ) THEN
         info = -3
      ELSE IF( lda.LT.max( 1, n ) ) THEN
         info = -5
      ELSE IF( ldb.LT.max( 1, n ) ) THEN
         info = -7
      ELSE IF( ldvl.LT.1 .OR. ( ilvl .AND. ldvl.LT.n ) ) THEN
         info = -12
      ELSE IF( ldvr.LT.1 .OR. ( ilvr .AND. ldvr.LT.n ) ) THEN
         info = -14
      ELSE IF( lwork.LT.lwkmin .AND. .NOT.lquery ) THEN
         info = -16
      END IF
*
*     Compute workspace
*
      IF( info.EQ.0 ) THEN
         CALL dgeqrf( n, n, b, ldb, work, work, -1, ierr )
         lwkopt = max( lwkmin, 3*n+int( work( 1 ) ) )
         CALL dormqr( 'L', 'T', n, n, n, b, ldb, work, a, lda,
     $                work, -1, ierr )
         lwkopt = max( lwkopt, 3*n+int( work( 1 ) ) )
         IF( ilvl ) THEN
            CALL dorgqr( n, n, n, vl, ldvl, work, work, -1, ierr )
            lwkopt = max( lwkopt, 3*n+int( work( 1 ) ) )
         END IF
         IF( ilv ) THEN
            CALL dgghd3( jobvl, jobvr, n, 1, n, a, lda, b, ldb, vl,
     $                   ldvl, vr, ldvr, work, -1, ierr )
            lwkopt = max( lwkopt, 3*n+int( work( 1 ) ) )
            CALL dlaqz0( 'S', jobvl, jobvr, n, 1, n, a, lda, b, ldb,
     $                   alphar, alphai, beta, vl, ldvl, vr, ldvr,
     $                   work, -1, 0, ierr )
            lwkopt = max( lwkopt, 2*n+int( work( 1 ) ) )
         ELSE
            CALL dgghd3( 'N', 'N', n, 1, n, a, lda, b, ldb, vl, ldvl,
     $                   vr, ldvr, work, -1, ierr )
            lwkopt = max( lwkopt, 3*n+int( work( 1 ) ) )
            CALL dlaqz0( 'E', jobvl, jobvr, n, 1, n, a, lda, b, ldb,
     $                   alphar, alphai, beta, vl, ldvl, vr, ldvr,
     $                   work, -1, 0, ierr )
            lwkopt = max( lwkopt, 2*n+int( work( 1 ) ) )
         END IF
         IF( n.EQ.0 ) THEN
            work( 1 ) = 1
         ELSE
            work( 1 ) = lwkopt
         END IF
      END IF
*
      IF( info.NE.0 ) THEN
         CALL xerbla( 'DGGEV3 ', -info )
         RETURN
      ELSE IF( lquery ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( n.EQ.0 )
     $   RETURN
*
*     Get machine constants
*
      eps = dlamch( 'P' )
      smlnum = dlamch( 'S' )
      bignum = one / smlnum
      smlnum = sqrt( smlnum ) / eps
      bignum = one / smlnum
*
*     Scale A if max element outside range [SMLNUM,BIGNUM]
*
      anrm = dlange( 'M', n, n, a, lda, work )
      ilascl = .false.
      IF( anrm.GT.zero .AND. anrm.LT.smlnum ) THEN
         anrmto = smlnum
         ilascl = .true.
      ELSE IF( anrm.GT.bignum ) THEN
         anrmto = bignum
         ilascl = .true.
      END IF
      IF( ilascl )
     $   CALL dlascl( 'G', 0, 0, anrm, anrmto, n, n, a, lda, ierr )
*
*     Scale B if max element outside range [SMLNUM,BIGNUM]
*
      bnrm = dlange( 'M', n, n, b, ldb, work )
      ilbscl = .false.
      IF( bnrm.GT.zero .AND. bnrm.LT.smlnum ) THEN
         bnrmto = smlnum
         ilbscl = .true.
      ELSE IF( bnrm.GT.bignum ) THEN
         bnrmto = bignum
         ilbscl = .true.
      END IF
      IF( ilbscl )
     $   CALL dlascl( 'G', 0, 0, bnrm, bnrmto, n, n, b, ldb, ierr )
*
*     Permute the matrices A, B to isolate eigenvalues if possible
*
      ileft = 1
      iright = n + 1
      iwrk = iright + n
      CALL dggbal( 'P', n, a, lda, b, ldb, ilo, ihi, work( ileft ),
     $             work( iright ), work( iwrk ), ierr )
*
*     Reduce B to triangular form (QR decomposition of B)
*
      irows = ihi + 1 - ilo
      IF( ilv ) THEN
         icols = n + 1 - ilo
      ELSE
         icols = irows
      END IF
      itau = iwrk
      iwrk = itau + irows
      CALL dgeqrf( irows, icols, b( ilo, ilo ), ldb, work( itau ),
     $             work( iwrk ), lwork+1-iwrk, ierr )
*
*     Apply the orthogonal transformation to matrix A
*
      CALL dormqr( 'L', 'T', irows, icols, irows, b( ilo, ilo ), ldb,
     $             work( itau ), a( ilo, ilo ), lda, work( iwrk ),
     $             lwork+1-iwrk, ierr )
*
*     Initialize VL
*
      IF( ilvl ) THEN
         CALL dlaset( 'Full', n, n, zero, one, vl, ldvl )
         IF( irows.GT.1 ) THEN
            CALL dlacpy( 'L', irows-1, irows-1, b( ilo+1, ilo ), ldb,
     $                   vl( ilo+1, ilo ), ldvl )
         END IF
         CALL dorgqr( irows, irows, irows, vl( ilo, ilo ), ldvl,
     $                work( itau ), work( iwrk ), lwork+1-iwrk, ierr )
      END IF
*
*     Initialize VR
*
      IF( ilvr )
     $   CALL dlaset( 'Full', n, n, zero, one, vr, ldvr )
*
*     Reduce to generalized Hessenberg form
*
      IF( ilv ) THEN
*
*        Eigenvectors requested -- work on whole matrix.
*
         CALL dgghd3( jobvl, jobvr, n, ilo, ihi, a, lda, b, ldb, vl,
     $                ldvl, vr, ldvr, work( iwrk ), lwork+1-iwrk, ierr )
      ELSE
         CALL dgghd3( 'N', 'N', irows, 1, irows, a( ilo, ilo ), lda,
     $                b( ilo, ilo ), ldb, vl, ldvl, vr, ldvr,
     $                work( iwrk ), lwork+1-iwrk, ierr )
      END IF
*
*     Perform QZ algorithm (Compute eigenvalues, and optionally, the
*     Schur forms and Schur vectors)
*
      iwrk = itau
      IF( ilv ) THEN
         chtemp = 'S'
      ELSE
         chtemp = 'E'
      END IF
      CALL dlaqz0( chtemp, jobvl, jobvr, n, ilo, ihi, a, lda, b, ldb,
     $             alphar, alphai, beta, vl, ldvl, vr, ldvr,
     $             work( iwrk ), lwork+1-iwrk, 0, ierr )
      IF( ierr.NE.0 ) THEN
         IF( ierr.GT.0 .AND. ierr.LE.n ) THEN
            info = ierr
         ELSE IF( ierr.GT.n .AND. ierr.LE.2*n ) THEN
            info = ierr - n
         ELSE
            info = n + 1
         END IF
         GO TO 110
      END IF
*
*     Compute Eigenvectors
*
      IF( ilv ) THEN
         IF( ilvl ) THEN
            IF( ilvr ) THEN
               chtemp = 'B'
            ELSE
               chtemp = 'L'
            END IF
         ELSE
            chtemp = 'R'
         END IF
         CALL dtgevc( chtemp, 'B', ldumma, n, a, lda, b, ldb, vl,
     $                ldvl,
     $                vr, ldvr, n, in, work( iwrk ), ierr )
         IF( ierr.NE.0 ) THEN
            info = n + 2
            GO TO 110
         END IF
*
*        Undo balancing on VL and VR and normalization
*
         IF( ilvl ) THEN
            CALL dggbak( 'P', 'L', n, ilo, ihi, work( ileft ),
     $                   work( iright ), n, vl, ldvl, ierr )
            DO 50 jc = 1, n
               IF( alphai( jc ).LT.zero )
     $            GO TO 50
               temp = zero
               IF( alphai( jc ).EQ.zero ) THEN
                  DO 10 jr = 1, n
                     temp = max( temp, abs( vl( jr, jc ) ) )
   10             CONTINUE
               ELSE
                  DO 20 jr = 1, n
                     temp = max( temp, abs( vl( jr, jc ) )+
     $                      abs( vl( jr, jc+1 ) ) )
   20             CONTINUE
               END IF
               IF( temp.LT.smlnum )
     $            GO TO 50
               temp = one / temp
               IF( alphai( jc ).EQ.zero ) THEN
                  DO 30 jr = 1, n
                     vl( jr, jc ) = vl( jr, jc )*temp
   30             CONTINUE
               ELSE
                  DO 40 jr = 1, n
                     vl( jr, jc ) = vl( jr, jc )*temp
                     vl( jr, jc+1 ) = vl( jr, jc+1 )*temp
   40             CONTINUE
               END IF
   50       CONTINUE
         END IF
         IF( ilvr ) THEN
            CALL dggbak( 'P', 'R', n, ilo, ihi, work( ileft ),
     $                   work( iright ), n, vr, ldvr, ierr )
            DO 100 jc = 1, n
               IF( alphai( jc ).LT.zero )
     $            GO TO 100
               temp = zero
               IF( alphai( jc ).EQ.zero ) THEN
                  DO 60 jr = 1, n
                     temp = max( temp, abs( vr( jr, jc ) ) )
   60             CONTINUE
               ELSE
                  DO 70 jr = 1, n
                     temp = max( temp, abs( vr( jr, jc ) )+
     $                      abs( vr( jr, jc+1 ) ) )
   70             CONTINUE
               END IF
               IF( temp.LT.smlnum )
     $            GO TO 100
               temp = one / temp
               IF( alphai( jc ).EQ.zero ) THEN
                  DO 80 jr = 1, n
                     vr( jr, jc ) = vr( jr, jc )*temp
   80             CONTINUE
               ELSE
                  DO 90 jr = 1, n
                     vr( jr, jc ) = vr( jr, jc )*temp
                     vr( jr, jc+1 ) = vr( jr, jc+1 )*temp
   90             CONTINUE
               END IF
  100       CONTINUE
         END IF
*
*        End of eigenvector calculation
*
      END IF
*
*     Undo scaling if necessary
*
  110 CONTINUE
*
      IF( ilascl ) THEN
         CALL dlascl( 'G', 0, 0, anrmto, anrm, n, 1, alphar, n,
     $                ierr )
         CALL dlascl( 'G', 0, 0, anrmto, anrm, n, 1, alphai, n,
     $                ierr )
      END IF
*
      IF( ilbscl ) THEN
         CALL dlascl( 'G', 0, 0, bnrmto, bnrm, n, 1, beta, n, ierr )
      END IF
*
      work( 1 ) = lwkopt
      RETURN
*
*     End of DGGEV3
*

      END