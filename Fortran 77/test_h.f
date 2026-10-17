* ========================================================
* Program：reorder the eigenvalue blocks of a 4-by-4 matrix in the Schur form by subroutine DLAEXC_H.
* ========================================================
      PROGRAM MATRIX_EXAMPLE
      INTEGER N, I, J
      PARAMETER ( N = 4 )
      DOUBLE PRECISION A( N, N ), A_WORK_H( N, N ), Q( N, N ),
     $                 WORK( 100 ), TOL
      CHARACTER*64 FNAME
      LOGICAL WANTQ
      DOUBLE PRECISION WR( N ), WI( N )
      DOUBLE PRECISION VL( N, N ), VR( N, N )  
      INTEGER LWORK, INFO
      CHARACTER JOBVL, JOBVR
*
*     Initialization.
*
      JOBVL = 'N'    
      JOBVR = 'N'   
      LWORK = 100  
      TOL = 0.0D0
      WR( N ) = 0.0D0
      WI( N ) = 0.0D0
      WORK( 100 ) = 0.0D0
*
*     Not wanted Q for our test.
*
      DO 10 I = 1, N
          DO 20 J = 1, N
              Q( I, J ) = 0.0D0
   20     CONTINUE
          Q( I, I ) = 1.0D0
   10 CONTINUE
      WANTQ = .FALSE.
*
*     Read and copy the test matrix.
*
      CALL GETARG( 1, FNAME )
*     If no argument provided, use default: 'test matrix 1.dat'.
      IF ( FNAME.EQ.' ' ) THEN
          FNAME = 'TEST MATRICES/TEST MATRIX 1.DAT'
      END IF
      CALL LOAD_TEST( FNAME, A )
      CALL COPY_MATRIX( A, A_WORK_H )
*
*     Print the test matrix before swapping.
*
      WRITE( *, '(A)' ) 'A before swapping: '
      CALL PRINT_A( A )
*
*     Print eigenvalues of the test matrix.
*
      CALL DGEEV( JOBVL, JOBVR, N, A, N, WR, WI,
     $            VL, N, VR, N, WORK, LWORK, INFO )
*
      IF ( INFO.EQ.0 ) THEN
          WRITE( *, '(A)' ) 'eigenvalues: '
          CALL PRINT_EIG( WR, WI )
      END IF
*
*     Call DLAEXC_H to swap adjacent blocks.
*
      CALL DLAEXC_H( WANTQ, 4, A_WORK_H, 4, Q, 4, 1, 2, 2, 
     $               WORK, INFO, TOL )
*
*     Check and print outcomes.
*
      IF ( INFO.EQ.0 ) THEN
*
*         If INFO = 0, print the test matrix after swapping
*         and its eigenvalues.
*
          WRITE( *, '(A)' ) 'A after swapping by DLAEXC_H: '
          CALL PRINT_A( A_WORK_H )
          CALL DGEEV( JOBVL, JOBVR, N, A_WORK_H, N, WR, WI,
     $                VL, N, VR, N, WORK, LWORK, INFO )
          IF ( INFO.EQ.0 ) THEN
              WRITE( *, '(A)' ) 'eigenvalues: '
              CALL PRINT_EIG( WR, WI )
          END IF
      ELSE
*
*         If INFO = 1, still print the test matrix after orthogonal 
*         similarity transformation. Then zero the lower left 2-by-2
*         block of the test matrix and print its eigenvalues.
*
          WRITE( *, '(A)' ) 'A fails to be reordered by DLAEXC_H: '
          CALL PRINT_A( A_WORK_H )
          A_WORK_H( 3, 1 ) = 0.0D0
          A_WORK_H( 3, 2 ) = 0.0D0
          A_WORK_H( 4, 1 ) = 0.0D0
          A_WORK_H( 4, 2 ) = 0.0D0
          CALL DGEEV( JOBVL, JOBVR, N, A_WORK_H, N, WR, WI,
     $                VL, N, VR, N, WORK, LWORK, INFO )    
          IF ( INFO.EQ.0 ) THEN
              WRITE( *, '(A)' ) 'eigenvalues: '
              CALL PRINT_EIG( WR, WI )
          ENDIF
      END IF
*
      WRITE( *, '(A)' ) 'tol: '
      WRITE( *, '(ES23.15)' ) TOL
      WRITE( *, * )
*
*     End of test_h.
*
      END
*
* ========================================================
* subroutine: load a 4-by-4 matrix from specified file.
* ========================================================
      SUBROUTINE LOAD_TEST( FNAME, MAT )
      CHARACTER*(*) FNAME
      DOUBLE PRECISION MAT( 4, 4 )
      INTEGER I, J, U
      PARAMETER ( U = 20 )
*
      OPEN( UNIT = U, FILE = FNAME, STATUS = 'OLD' )
      DO 30 I = 1, 4
          READ( U, * ) ( MAT( I, J ), J = 1, 4 )
   30 CONTINUE
      CLOSE( U )
      RETURN
      END
*
* ========================================================
* subroutine: print a 4-by-4 matrix.
* ========================================================
      SUBROUTINE PRINT_A( MAT )
      DOUBLE PRECISION MAT( 4, 4 )
      INTEGER I, J
*
      DO 40 I = 1, 4
          WRITE( *, '(4ES14.6)' ) ( MAT( I, J ), J = 1, 4 )
   40 CONTINUE
      WRITE( *, * )
      RETURN
      END
*
* ========================================================
* subroutine: copy a 4-by-4 matrix.
* ========================================================
      SUBROUTINE COPY_MATRIX( SRC, DEST )
      DOUBLE PRECISION SRC( 4, 4 ), DEST( 4, 4 )
      INTEGER I, J
*
      DO 50 I = 1, 4
          DO 60 J = 1, 4
              DEST( I, J ) = SRC( I, J )
   60     CONTINUE
   50 CONTINUE
      RETURN
      END
*
* ========================================================
* subroutine: print eigenvalues.
* ========================================================
      SUBROUTINE PRINT_EIG( WR, WI )
      INTEGER N
      PARAMETER ( N = 4 )
      DOUBLE PRECISION WR( N ), WI( N )   
*      
      DO 30 I = 1, N
          IF ( WI( I ).EQ.0.0D0 ) THEN
              WRITE( *, 40 ) I, WR( I )
   40         FORMAT( 'Lambda', I1,' = ', ES22.15, ' (real number)' )
              WRITE( *, * )
          ELSE
              WRITE( *, 50 ) I, WR( I ), WI( I )
   50         FORMAT( 'Lambda', I1,' = ', ES22.15,' + ',ES22.15, 'i' )
              WRITE( *, * )
          ENDIF
   30 CONTINUE
      END