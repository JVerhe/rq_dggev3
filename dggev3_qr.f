*> \brief <b> DGGEV3 computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE matrices (blocked algorithm)</b>
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download DGGEV3 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dggev3.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dggev3.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dggev3.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE DGGEV3( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR,
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
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DGGEV3 computes for a pair of N-by-N real nonsymmetric matrices (A,B)
*> the generalized eigenvalues, and optionally, the left and/or right
*> generalized eigenvectors.
*>
*> A generalized eigenvalue for a pair of matrices (A,B) is a scalar
*> lambda or a ratio alpha/beta = lambda, such that A - lambda*B is
*> singular. It is usually represented as the pair (alpha,beta), as
*> there is a reasonable interpretation for beta=0, and even for both
*> being zero.
*>
*> The right eigenvector v(j) corresponding to the eigenvalue lambda(j)
*> of (A,B) satisfies
*>
*>                  A * v(j) = lambda(j) * B * v(j).
*>
*> The left eigenvector u(j) corresponding to the eigenvalue lambda(j)
*> of (A,B) satisfies
*>
*>                  u(j)**H * A  = lambda(j) * u(j)**H * B .
*>
*> where u(j)**H is the conjugate-transpose of u(j).
*>
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] JOBVL
*> \verbatim
*>          JOBVL is CHARACTER*1
*>          = 'N':  do not compute the left generalized eigenvectors;
*>          = 'V':  compute the left generalized eigenvectors.
*> \endverbatim
*>
*> \param[in] JOBVR
*> \verbatim
*>          JOBVR is CHARACTER*1
*>          = 'N':  do not compute the right generalized eigenvectors;
*>          = 'V':  compute the right generalized eigenvectors.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrices A, B, VL, and VR.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA, N)
*>          On entry, the matrix A in the pair (A,B).
*>          On exit, A has been overwritten.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is DOUBLE PRECISION array, dimension (LDB, N)
*>          On entry, the matrix B in the pair (A,B).
*>          On exit, B has been overwritten.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of B.  LDB >= max(1,N).
*> \endverbatim
*>
*> \param[out] ALPHAR
*> \verbatim
*>          ALPHAR is DOUBLE PRECISION array, dimension (N)
*> \endverbatim
*>
*> \param[out] ALPHAI
*> \verbatim
*>          ALPHAI is DOUBLE PRECISION array, dimension (N)
*> \endverbatim
*>
*> \param[out] BETA
*> \verbatim
*>          BETA is DOUBLE PRECISION array, dimension (N)
*>          On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will
*>          be the generalized eigenvalues.  If ALPHAI(j) is zero, then
*>          the j-th eigenvalue is real; if positive, then the j-th and
*>          (j+1)-st eigenvalues are a complex conjugate pair, with
*>          ALPHAI(j+1) negative.
*>
*>          Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j)
*>          may easily over- or underflow, and BETA(j) may even be zero.
*>          Thus, the user should avoid naively computing the ratio
*>          alpha/beta.  However, ALPHAR and ALPHAI will be always less
*>          than and usually comparable with norm(A) in magnitude, and
*>          BETA always less than and usually comparable with norm(B).
*> \endverbatim
*>
*> \param[out] VL
*> \verbatim
*>          VL is DOUBLE PRECISION array, dimension (LDVL,N)
*>          If JOBVL = 'V', the left eigenvectors u(j) are stored one
*>          after another in the columns of VL, in the same order as
*>          their eigenvalues. If the j-th eigenvalue is real, then
*>          u(j) = VL(:,j), the j-th column of VL. If the j-th and
*>          (j+1)-th eigenvalues form a complex conjugate pair, then
*>          u(j) = VL(:,j)+i*VL(:,j+1) and u(j+1) = VL(:,j)-i*VL(:,j+1).
*>          Each eigenvector is scaled so the largest component has
*>          abs(real part)+abs(imag. part)=1.
*>          Not referenced if JOBVL = 'N'.
*> \endverbatim
*>
*> \param[in] LDVL
*> \verbatim
*>          LDVL is INTEGER
*>          The leading dimension of the matrix VL. LDVL >= 1, and
*>          if JOBVL = 'V', LDVL >= N.
*> \endverbatim
*>
*> \param[out] VR
*> \verbatim
*>          VR is DOUBLE PRECISION array, dimension (LDVR,N)
*>          If JOBVR = 'V', the right eigenvectors v(j) are stored one
*>          after another in the columns of VR, in the same order as
*>          their eigenvalues. If the j-th eigenvalue is real, then
*>          v(j) = VR(:,j), the j-th column of VR. If the j-th and
*>          (j+1)-th eigenvalues form a complex conjugate pair, then
*>          v(j) = VR(:,j)+i*VR(:,j+1) and v(j+1) = VR(:,j)-i*VR(:,j+1).
*>          Each eigenvector is scaled so the largest component has
*>          abs(real part)+abs(imag. part)=1.
*>          Not referenced if JOBVR = 'N'.
*> \endverbatim
*>
*> \param[in] LDVR
*> \verbatim
*>          LDVR is INTEGER
*>          The leading dimension of the matrix VR. LDVR >= 1, and
*>          if JOBVR = 'V', LDVR >= N.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER.
*>          The dimension of the array WORK. LWORK >= MAX(1,8*N).
*>          For good performance, LWORK should generally be larger.
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the WORK array, and no error
*>          message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
*>          = 1,...,N:
*>                The QZ iteration failed.  No eigenvectors have been
*>                calculated, but ALPHAR(j), ALPHAI(j), and BETA(j)
*>                should be correct for j=INFO+1,...,N.
*>          > N:  =N+1: other than QZ iteration failed in DLAQZ0.
*>                =N+2: error return from DTGEVC.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup ggev3
*
*  =====================================================================

      SUBROUTINE dggev3_qr( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR,
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
      CALL dggbal( 'N', n, a, lda, b, ldb, ilo, ihi, work( ileft ),
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