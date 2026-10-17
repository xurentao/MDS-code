*> \brief \b DLAEXC_G swaps adjacent diagonal blocks of a real upper quasi-triangular matrix in Schur canonical form, by an orthogonal similarity transformation.
*
*  =========== DOCUMENTATION ===========
*
*  Definition:
*  ===========
*
*       SUBROUTINE DLAEXC_G( WANTQ, N, T, LDT, Q, LDQ, J1, N1, N2, WORK,
*                            INFO, TOL )
*
*       .. Scalar Arguments ..
*       LOGICAL            WANTQ
*       INTEGER            INFO, J1, LDQ, LDT, N, N1, N2
*       DOUBLE             TOL
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   Q( LDQ, * ), T( LDT, * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLAEXC_G is a test subroutine to swap adjacent diagonal blocks T11
*> and T22 of order 1 or 2 in an upper quasi-triangular matrix T by an
*> orthogonal similarity transformation.
*>
*> T must be in Schur canonical form, that is, block upper triangular
*> with 1-by-1 and 2-by-2 diagonal blocks; each 2-by-2 diagonal block
*> has its diagonal elements equal and its off-diagonal elements of
*> opposite sign.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] WANTQ
*> \verbatim
*>          WANTQ is LOGICAL
*>          = .TRUE. : accumulate the transformation in the matrix Q;
*>          = .FALSE.: do not accumulate the transformation.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix T. N >= 0.
*> \endverbatim
*>
*> \param[in,out] T
*> \verbatim
*>          T is DOUBLE PRECISION array, dimension (LDT,N)
*>          On entry, the upper quasi-triangular matrix T, in Schur
*>          canonical form.
*>          On exit, the updated matrix T, again in Schur canonical form.
*> \endverbatim
*>
*> \param[in] LDT
*> \verbatim
*>          LDT is INTEGER
*>          The leading dimension of the array T. LDT >= max(1,N).
*> \endverbatim
*>
*> \param[in,out] Q
*> \verbatim
*>          Q is DOUBLE PRECISION array, dimension (LDQ,N)
*>          On entry, if WANTQ is .TRUE., the orthogonal matrix Q.
*>          On exit, if WANTQ is .TRUE., the updated matrix Q.
*>          If WANTQ is .FALSE., Q is not referenced.
*> \endverbatim
*>
*> \param[in] LDQ
*> \verbatim
*>          LDQ is INTEGER
*>          The leading dimension of the array Q.
*>          LDQ >= 1; and if WANTQ is .TRUE., LDQ >= N.
*> \endverbatim
*>
*> \param[in] J1
*> \verbatim
*>          J1 is INTEGER
*>          The index of the first row of the first block T11.
*> \endverbatim
*>
*> \param[in] N1
*> \verbatim
*>          N1 is INTEGER
*>          The order of the first block T11. N1 = 0, 1 or 2.
*> \endverbatim
*>
*> \param[in] N2
*> \verbatim
*>          N2 is INTEGER
*>          The order of the second block T22. N2 = 0, 1 or 2.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (N)
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0: successful exit
*>          = 1: the transformed matrix T would be too far from Schur
*>               form. The transformation would still be applied to the 
*>               test matrix, but Q would not be changed.
*> \endverbatim
*>
*> \param[out] TOL
*> \verbatim
*>          TOL is DOUBLE
*>          Used to quantify the algorithmic error in swapping adjacent
*>          2-by-2 diagonal blocks. For numerical experiments only.
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
*> \ingroup doubleOTHERauxiliary
*
*  =====================================================================
      SUBROUTINE DLAEXC_G( WANTQ, N, T, LDT, Q, LDQ, J1, N1, N2, 
     $                     WORK, INFO, TOL )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      LOGICAL            WANTQ
      INTEGER            INFO, J1, LDQ, LDT, N, N1, N2
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   Q( LDQ, * ), T( LDT, * ), WORK( * ), TOL
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   TEN
      PARAMETER          ( TEN = 1.0D+1 )
      INTEGER            LDD, LDX
      PARAMETER          ( LDD = 4, LDX = 2 )
*     ..
*     .. Local Scalars ..
      INTEGER            IERR, J2, J3, J4, K, ND
      DOUBLE PRECISION   CS, DNORM, EPS, SCALE, SMLNUM, SN, T11, T22,
     $                   T33, TAU, TEMP, THRESH, WI1, WI2, WR1, WR2,
     $                   XNORM, MAX_COL1, MAX_COL2, M_SCALE, M_ZERO 
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   D( LDD, 4 ), U( 3 ), X( LDX, 2 )
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   S( 4 ), C( 4 )
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLANGE
      EXTERNAL           DLAMCH, DLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLACPY, DLANV2, DLARFG, DLARFX, DLARTG, DLASY2,
     $                   DROT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     ..
*     .. Executable Statements ..
*
      INFO = 0
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. N1.EQ.0 .OR. N2.EQ.0 )
     $   RETURN
      IF( J1+N1.GT.N )
     $   RETURN
*
      J2 = J1 + 1
      J3 = J1 + 2
      J4 = J1 + 3
*
      IF( N1.EQ.1 .AND. N2.EQ.1 ) THEN
*
*        Swap two 1-by-1 blocks.
*
         T11 = T( J1, J1 )
         T22 = T( J2, J2 )
*
*        Determine the transformation to perform the interchange.
*
         CALL DLARTG( T( J1, J2 ), T22-T11, CS, SN, TEMP )
*
*        Apply transformation to the matrix T.
*
         IF( J3.LE.N )
     $      CALL DROT( N-J1-1, T( J1, J3 ), LDT, T( J2, J3 ), LDT, CS,
     $                 SN )
         CALL DROT( J1-1, T( 1, J1 ), 1, T( 1, J2 ), 1, CS, SN )
*
         T( J1, J1 ) = T22
         T( J2, J2 ) = T11
*
         IF( WANTQ ) THEN
*
*           Accumulate transformation in the matrix Q.
*
            CALL DROT( N, Q( 1, J1 ), 1, Q( 1, J2 ), 1, CS, SN )
         END IF
*
      ELSE
*
*        Swapping involves at least one 2-by-2 block.
*
*        Copy the diagonal block of order N1+N2 to the local array D
*        and compute its norm.
*
         ND = N1 + N2
         CALL DLACPY( 'Full', ND, ND, T( J1, J1 ), LDT, D, LDD )
         DNORM = DLANGE( 'Max', ND, ND, D, LDD, WORK )
*
*        Compute machine-dependent threshold for test for accepting
*        swap.
*
         EPS = DLAMCH( 'P' )
         SMLNUM = DLAMCH( 'S' ) / EPS
         THRESH = MAX( TEN*EPS*DNORM, SMLNUM )
*
*        Solve T11*X - X*T22 = scale*T12 for X.
*
         CALL DLASY2( .FALSE., .FALSE., -1, N1, N2, D, LDD,
     $                D( N1+1, N1+1 ), LDD, D( 1, N1+1 ), LDD, SCALE, X,
     $                LDX, XNORM, IERR )
*
*        Swap the adjacent diagonal blocks.
*
         K = N1 + N1 + N2 - 3
         GO TO ( 10, 20, 30 )K
*
   10    CONTINUE
*
*        N1 = 1, N2 = 2: generate elementary reflector H so that:
*
*        ( scale, X11, X12 ) H = ( 0, 0, * )
*
         U( 1 ) = SCALE
         U( 2 ) = X( 1, 1 )
         U( 3 ) = X( 1, 2 )
         CALL DLARFG( 3, U( 3 ), U, 1, TAU )
         U( 3 ) = ONE
         T11 = T( J1, J1 )
*
*        Perform swap provisionally on diagonal block in D.
*
         CALL DLARFX( 'L', 3, 3, U, TAU, D, LDD, WORK )
         CALL DLARFX( 'R', 3, 3, U, TAU, D, LDD, WORK )
*
*        Test whether to reject swap.
*
         IF( MAX( ABS( D( 3, 1 ) ), ABS( D( 3, 2 ) ), ABS( D( 3,
     $       3 )-T11 ) ).GT.THRESH )GO TO 50
*
*        Accept swap: apply transformation to the entire matrix T.
*
         CALL DLARFX( 'L', 3, N-J1+1, U, TAU, T( J1, J1 ), LDT, WORK )
         CALL DLARFX( 'R', J2, 3, U, TAU, T( 1, J1 ), LDT, WORK )
*
         T( J3, J1 ) = ZERO
         T( J3, J2 ) = ZERO
         T( J3, J3 ) = T11
*
         IF( WANTQ ) THEN
*
*           Accumulate transformation in the matrix Q.
*
            CALL DLARFX( 'R', N, 3, U, TAU, Q( 1, J1 ), LDQ, WORK )
         END IF
         GO TO 40
*
   20    CONTINUE
*
*        N1 = 2, N2 = 1: generate elementary reflector H so that:
*
*        H (  -X11 ) = ( * )
*          (  -X21 ) = ( 0 )
*          ( scale ) = ( 0 )
*
         U( 1 ) = -X( 1, 1 )
         U( 2 ) = -X( 2, 1 )
         U( 3 ) = SCALE
         CALL DLARFG( 3, U( 1 ), U( 2 ), 1, TAU )
         U( 1 ) = ONE
         T33 = T( J3, J3 )
*
*        Perform swap provisionally on diagonal block in D.
*
         CALL DLARFX( 'L', 3, 3, U, TAU, D, LDD, WORK )
         CALL DLARFX( 'R', 3, 3, U, TAU, D, LDD, WORK )
*
*        Test whether to reject swap.
*
         IF( MAX( ABS( D( 2, 1 ) ), ABS( D( 3, 1 ) ), ABS( D( 1,
     $       1 )-T33 ) ).GT.THRESH )GO TO 50
*
*        Accept swap: apply transformation to the entire matrix T.
*
         CALL DLARFX( 'R', J3, 3, U, TAU, T( 1, J1 ), LDT, WORK )
         CALL DLARFX( 'L', 3, N-J1, U, TAU, T( J1, J2 ), LDT, WORK )
*
         T( J1, J1 ) = T33
         T( J2, J1 ) = ZERO
         T( J3, J1 ) = ZERO
*
         IF( WANTQ ) THEN
*
*           Accumulate transformation in the matrix Q.
*
            CALL DLARFX( 'R', N, 3, U, TAU, Q( 1, J1 ), LDQ, WORK )
         END IF
         GO TO 40
*
   30    CONTINUE
*
*        N1 = 2, N2 = 2: generate Givens rotations G(1), G(2), G(3) and G(4) so
*        that:
*
*        G(4) G(3) G(2) G(1) (  -X11  -X12 ) = (  *  * )
*                            (  -X21  -X22 )   (  0  * )
*                            ( scale    0  )   (  0  0 )
*                            (    0  scale )   (  0  0 )
*        when the largest absolute value in X is attained at X11 or X21
*
*        otherwise:
* 
*        G(4) G(3) G(2) G(1) (  -X11  -X12 ) = (  *  * )
*                            (  -X21  -X22 )   (  *  0 )
*                            ( scale    0  )   (  0  0 )
*                            (    0  scale )   (  0  0 )
*
         X = -X;
         MAX_COL1 = MAX( ABS( X( 1, 1 ) ), ABS( X( 2, 1 ) ) )
         MAX_COL2 = MAX( ABS( X( 1, 2 ) ), ABS( X( 2, 2 ) ) )
         M_SCALE = SCALE 
         M_ZERO = ZERO
*
         IF ( MAX_COL1 .GE. MAX_COL2 ) THEN
*
            CALL DLARTG( X( 1, 1 ), X( 2, 1 ), C( 1 ), S( 1 ),
     $                    X( 1, 1 ) )
            CALL DROT( 1, X( 1, 2 ), 1, X( 2, 2 ), 1, C( 1 ), S( 1 ) )
            CALL DLARTG( X( 1, 1 ), M_SCALE, C( 2 ), S( 2 ), 
     $                    X( 1, 1 ) )
            CALL DROT( 1, X( 1, 2 ), 1, M_ZERO, 1, C( 2 ), S( 2 ) )
            CALL DLARTG( X( 2, 2 ), M_SCALE, C( 3 ), S( 3 ), 
     $                    X( 2, 2 ) )
            CALL DLARTG( X( 2, 2 ), M_ZERO, C( 4 ), S( 4 ), X( 2, 2 ) )
*
         ELSE
*
            CALL DLARTG( X( 1, 2 ), X( 2, 2 ), C( 1 ), S( 1 ),
     $                    X( 1, 2 ) )
            CALL DROT( 1, X( 1, 1 ), 1, X( 2, 1 ), 1, C( 1 ), S( 1 ) )
            CALL DLARTG( X( 1, 2 ), M_SCALE, C( 2 ), S( 2 ), 
     $                    X( 1, 2 ) )
            CALL DROT( 1, X( 1, 1 ), 1, M_ZERO, 1, C( 2 ), S( 2 ) )
            CALL DLARTG( X( 2, 1 ), M_SCALE, C( 3 ), S( 3 ), 
     $                    X( 2, 1 ) )
            CALL DLARTG( X( 2, 1 ), M_ZERO, C( 4 ), S( 4 ), X( 2, 1 ) ) 
*
         END IF       
*
*        Perform swap provisionally on diagonal block in D.
*
         CALL DROT( 4, D( 1, 1 ), 4, D( 2, 1 ), 4, C( 1 ), S( 1 ) )
         CALL DROT( 4, D( 1, 1 ), 1, D( 1, 2 ), 1, C( 1 ), S( 1 ) )
*
         IF ( MAX_COL1 .GE. MAX_COL2 ) THEN
*
            CALL DROT( 4, D( 1, 1 ), 4, D( 3, 1 ), 4, C( 2 ), S( 2 ) )
            CALL DROT( 4, D( 1, 1 ), 1, D( 1, 3 ), 1, C( 2 ), S( 2 ) )
            CALL DROT( 4, D( 2, 1 ), 4, D( 4, 1 ), 4, C( 3 ), S( 3 ) )
            CALL DROT( 4, D( 1, 2 ), 1, D( 1, 4 ), 1, C( 3 ), S( 3 ) )
            CALL DROT( 4, D( 2, 1 ), 4, D( 3, 1 ), 4, C( 4 ), S( 4 ) )
            CALL DROT( 4, D( 1, 2 ), 1, D( 1, 3 ), 1, C( 4 ), S( 4 ) )
*
         ELSE
*
            CALL DROT( 4, D( 1, 1 ), 4, D( 4, 1 ), 4, C( 2 ), S( 2 ) )
            CALL DROT( 4, D( 1, 1 ), 1, D( 1, 4 ), 1, C( 2 ), S( 2 ) )
            CALL DROT( 4, D( 2, 1 ), 4, D( 3, 1 ), 4, C( 3 ), S( 3 ) )
            CALL DROT( 4, D( 1, 2 ), 1, D( 1, 3 ), 1, C( 3 ), S( 3 ) )
            CALL DROT( 4, D( 2, 1 ), 4, D( 4, 1 ), 4, C( 4 ), S( 4 ) )
            CALL DROT( 4, D( 1, 2 ), 1, D( 1, 4 ), 1, C( 4 ), S( 4 ) )
*
         END IF
*
*        Experimental error metric.
*
         TOL = MAX( ABS( D( 3, 1 ) ), ABS( D( 3, 2 ) ), 
     $              ABS( D( 4, 1 ) ), ABS( D( 4, 2 ) ) ) / DNORM
*
*        Test whether to reject swap.
*
         IF( MAX( ABS( D( 3, 1 ) ), ABS( D( 3, 2 ) ), ABS( D( 4, 1 ) ),
     $       ABS( D( 4, 2 ) ) ).GT.THRESH )GO TO 50
*
*
*        Accept swap: apply transformation to the entire matrix T.
*
         CALL DROT( N-J1+1, T( J1, J1 ), LDT, T( J2, J1 ), 
     $              LDT, C( 1 ), S( 1 ) )
         CALL DROT( J2, T( 1, J1 ), 1, T( 1, J2 ), 1,
     $              C( 1 ), S( 1 ) )
*
         IF (MAX_COL1 .GE. MAX_COL2) THEN
*
            CALL DROT( N-J1+1, T( J1, J1 ), LDT, T( J3, J1 ), 
     $                 LDT, C( 2 ), S( 2 ) )
            CALL DROT( J4, T( 1, J1 ), 1, T( 1, J3 ), 1,
     $                 C( 2 ), S( 2 ) )
            CALL DROT( N-J1+1, T( J2, J1 ), LDT, T( J4, J1 ), 
     $                 LDT, C( 3 ), S( 3 ) )
            CALL DROT( J4, T( 1, J2 ), 1, T( 1, J4 ), 1,
     $                 C( 3 ), S( 3 ) )
            CALL DROT( N-J1+1, T( J2, J1 ), LDT, T( J3, J1 ), 
     $                 LDT, C( 4 ), S( 4 ) )
            CALL DROT( J4, T( 1, J2 ), 1, T( 1, J3 ), 1,
     $                 C( 4 ), S( 4 ) )
         ELSE
*
            CALL DROT( N-J1+1, T( J1, J1 ), LDT, T( J4, J1 ), 
     $                 LDT, C( 2 ), S( 2 ) )
            CALL DROT( J4, T( 1, J1 ), 1, T( 1, J4 ), 1,
     $                 C( 2 ), S( 2 ) )
            CALL DROT( N-J1+1, T( J2, J1 ), LDT, T( J3, J1 ), 
     $                 LDT, C( 3 ), S( 3 ) )
            CALL DROT( J4, T( 1, J2 ), 1, T( 1, J3 ), 1,
     $                 C( 3 ), S( 3 ) )
            CALL DROT( N-J1+1, T( J2, J1 ), LDT, T( J4, J1 ), 
     $                 LDT, C( 4 ), S( 4 ) )
            CALL DROT( J4, T( 1, J2 ), 1, T( 1, J4 ), 1,
     $                 C( 4 ), S( 4 ) )
         END IF
*
         T( J3, J1 ) = ZERO
         T( J3, J2 ) = ZERO
         T( J4, J1 ) = ZERO
         T( J4, J2 ) = ZERO
*
         IF( WANTQ ) THEN
*
*           Accumulate transformation in the matrix Q.
*
             CALL DROT( N, Q( 1, J1 ), 1, Q( 1, J2 ), 1,
     $                  C( 1 ), S( 1 ) )
*
             IF ( MAX_COL1 .GE. MAX_COL2 ) THEN
*           
                 CALL DROT( N, Q( 1, J1 ), 1, Q( 1, J3 ), 1,
     $                      C( 2 ), S( 2 ) )
                 CALL DROT( N, Q( 1, J2 ), 1, Q( 1, J4 ), 1,
     $                  C( 3 ), S( 3 ) )
                 CALL DROT( N, Q( 1, J2 ), 1, Q( 1, J3 ), 1,
     $                  C( 4 ), S( 4 ) ) 
*
             ELSE
*
                 CALL DROT( N, Q( 1, J1 ), 1, Q( 1, J4 ), 1,
     $                      C( 2 ), S( 2 ) )
                 CALL DROT( N, Q( 1, J2 ), 1, Q( 1, J3 ), 1,
     $                  C( 3 ), S( 3 ) )
                 CALL DROT( N, Q( 1, J2 ), 1, Q( 1, J4 ), 1,
     $                  C( 4 ), S( 4 ) )
*
             END IF
*
         END IF
*
   40    CONTINUE
*
         IF( N2.EQ.2 ) THEN
*
*           Standardize new 2-by-2 block T11
*
            CALL DLANV2( T( J1, J1 ), T( J1, J2 ), T( J2, J1 ),
     $                   T( J2, J2 ), WR1, WI1, WR2, WI2, CS, SN )
            CALL DROT( N-J1-1, T( J1, J1+2 ), LDT, T( J2, J1+2 ), LDT,
     $                 CS, SN )
            CALL DROT( J1-1, T( 1, J1 ), 1, T( 1, J2 ), 1, CS, SN )
            IF( WANTQ )
     $         CALL DROT( N, Q( 1, J1 ), 1, Q( 1, J2 ), 1, CS, SN )
         END IF
*
         IF( N1.EQ.2 ) THEN
*
*           Standardize new 2-by-2 block T22
*
            J3 = J1 + N2
            J4 = J3 + 1
            CALL DLANV2( T( J3, J3 ), T( J3, J4 ), T( J4, J3 ),
     $                   T( J4, J4 ), WR1, WI1, WR2, WI2, CS, SN )
            IF( J3+2.LE.N )
     $         CALL DROT( N-J3-1, T( J3, J3+2 ), LDT, T( J4, J3+2 ),
     $                    LDT, CS, SN )
            CALL DROT( J3-1, T( 1, J3 ), 1, T( 1, J4 ), 1, CS, SN )
            IF( WANTQ )
     $         CALL DROT( N, Q( 1, J3 ), 1, Q( 1, J4 ), 1, CS, SN )
         END IF
*
      END IF
      RETURN
*
*     Exit with INFO = 1 if swap was rejected.
*
   50 CONTINUE
      INFO = 1
*
*     For numerical experiments only.
*
      CALL DROT( N-J1+1, T( J1, J1 ), LDT, T( J2, J1 ), 
     $                 LDT, C( 1 ), S( 1 ) )
      CALL DROT( J2, T( 1, J1 ), 1, T( 1, J2 ), 1,
     $                 C( 1 ), S( 1 ) )
*
      IF (MAX_COL1 .GE. MAX_COL2) THEN
*
         CALL DROT( N-J1+1, T( J1, J1 ), LDT, T( J3, J1 ), 
     $              LDT, C( 2 ), S( 2 ) )
         CALL DROT( J4, T( 1, J1 ), 1, T( 1, J3 ), 1,
     $              C( 2 ), S( 2 ) )
         CALL DROT( N-J1+1, T( J2, J1 ), LDT, T( J4, J1 ), 
     $              LDT, C( 3 ), S( 3 ) )
         CALL DROT( J4, T( 1, J2 ), 1, T( 1, J4 ), 1,
     $              C( 3 ), S( 3 ) )
         CALL DROT( N-J1+1, T( J2, J1 ), LDT, T( J3, J1 ), 
     $              LDT, C( 4 ), S( 4 ) )
         CALL DROT( J4, T( 1, J2 ), 1, T( 1, J3 ), 1,
     $              C( 4 ), S( 4 ) )
      ELSE
*
         CALL DROT( N-J1+1, T( J1, J1 ), LDT, T( J4, J1 ), 
     $              LDT, C( 2 ), S( 2 ) )
         CALL DROT( J4, T( 1, J1 ), 1, T( 1, J4 ), 1,
     $              C( 2 ), S( 2 ) )
         CALL DROT( N-J1+1, T( J2, J1 ), LDT, T( J3, J1 ), 
     $              LDT, C( 3 ), S( 3 ) )
         CALL DROT( J4, T( 1, J2 ), 1, T( 1, J3 ), 1,
     $              C( 3 ), S( 3 ) )
         CALL DROT( N-J1+1, T( J2, J1 ), LDT, T( J4, J1 ), 
     $              LDT, C( 4 ), S( 4 ) )
         CALL DROT( J4, T( 1, J2 ), 1, T( 1, J4 ), 1,
     $              C( 4 ), S( 4 ) )
      END IF
*
      RETURN
*
*     End of DLAEXC_G
*
      END