SUBROUTINE PACKINGBANDMATRIX(N,A,KD,AB,INFO)

! brute force packing of a banded matrix
  
  IMPLICIT NONE
  INTEGER, INTENT(INOUT) :: INFO
  INTEGER, INTENT(IN)    :: N,KD
  COMPLEX*16, DIMENSION(N,N)    :: A
  COMPLEX*16, DIMENSION(KD+1,N) :: AB

  INTEGER i,j,N_,I_

  
  
  DO j=1,N
     IF(1.GE.(j-KD)) THEN
        I_=1
     ELSE
        I_ = j-KD
     END IF
     ! U STORAGE
     DO i=I_,j
        AB(KD+1+i-j,j) = A(i,j)
     END DO
  END DO

END SUBROUTINE PACKINGBANDMATRIX


SUBROUTINE LAPACK_FULLEIGENVALUESBAND(AB,Z,KD,N,W,INFO)
!eigenvalues/vectors of banded matrix ab
!AB, inout, packed banded matrix
!Z, out,eigenvectors
!KD out, calcuated eigenvectors
!N, in,matrix dimension
!W, out, eigenvalues
!INFO,inout, error flag

  !H is COMPLEX*16 array, dimension (N, N)
  !  69 *>          On entry, the Hermitian matrix A.  If UPLO = 'U', the
  !  70 *>          leading N-by-N upper triangular part of A contains the
  !  71 *>          upper triangular part of the matrix A.  If UPLO = 'L',
  !  72 *>          the leading N-by-N lower triangular part of A contains
  !  73 *>          the lower triangular part of the matrix A.
  !  74 *>          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
  !  75 *>          orthonormal eigenvectors of the matrix A.
  !  76 *>          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
  !  77 *>          or the upper triangle (if UPLO='U') of A, including the
  !  78 *>          diagonal, is destroyed.
  !
  ! The eigenvector H(:,r) corresponds to the eigenvalue W_SPACE(r)
  !
  IMPLICIT NONE
  INTEGER,                                INTENT(IN)    :: N,KD
  COMPLEX*16,       DIMENSION(KD+1,N), INTENT(INOUT)    :: AB
  COMPLEX*16,       DIMENSION(N,N),       INTENT(INOUT) :: Z
  DOUBLE PRECISION, DIMENSION(N),         INTENT(INOUT) :: W
  INTEGER,                                INTENT(OUT)   :: INFO


  !---SETTING  LAPACK VARIABLES: START ---------!
  CHARACTER         JOBZ, UPLO
  INTEGER           LDAB,LDZ
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: RWORK
  COMPLEX*16,       DIMENSION(:),   ALLOCATABLE :: WORK
  
  JOBZ = 'V'
  UPLO = 'U'  

  LDAB = KD+1
  LDZ  = N

  ALLOCATE(WORK(N))
  ALLOCATE(RWORK(3*N-2))

!  WRITE(*,*) JOBZ, UPLO, N, KD,  LDAB, LDZ,  INFO
  CALL ZHBEV(JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, WORK, RWORK, INFO )
  IF(INFO /= 0) WRITE(*,*) "# DIAG FAIL 1",INFO

  DEALLOCATE(WORK)
  DEALLOCATE(RWORK) 
 
END SUBROUTINE LAPACK_FULLEIGENVALUESBAND




SUBROUTINE LAPACK_FULLEIGENVALUES(H,N,W_SPACE,INFO)
!eigenvalues/vectors of matrix ab
!H, inout, packed banded matrix
! , out,eigenvectors
!N, in,matrix dimension
!W_space, out, eigenvalues
!INFO,inout, error flag

  !H is COMPLEX*16 array, dimension (N, N)
  !  69 *>          On entry, the Hermitian matrix A.  If UPLO = 'U', the
  !  70 *>          leading N-by-N upper triangular part of A contains the
  !  71 *>          upper triangular part of the matrix A.  If UPLO = 'L',
  !  72 *>          the leading N-by-N lower triangular part of A contains
  !  73 *>          the lower triangular part of the matrix A.
  !  74 *>          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
  !  75 *>          orthonormal eigenvectors of the matrix A.
  !  76 *>          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
  !  77 *>          or the upper triangle (if UPLO='U') of A, including the
  !  78 *>          diagonal, is destroyed.
  !
  ! The eigenvector H(:,r) corresponds to the eigenvalue W_SPACE(r)
  !
  IMPLICIT NONE
  INTEGER,                          INTENT(IN)    :: N
  COMPLEX*16,       DIMENSION(:,:), INTENT(INOUT) :: H
  DOUBLE PRECISION, DIMENSION(:),   INTENT(INOUT) :: W_SPACE
  INTEGER,                          INTENT(OUT)   :: INFO


  !---SETTING  LAPACK VARIABLES: START ---------!
  CHARACTER         JOBZ, UPLO
  INTEGER           LWORK,LDA
  INTEGER,          DIMENSION(:),   ALLOCATABLE :: IWORK,ISUPPZ
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: RWORK
  COMPLEX*16,       DIMENSION(:),   ALLOCATABLE :: WORK
  COMPLEX*16,       DIMENSION(:,:), ALLOCATABLE :: H_AUX_  
  JOBZ = 'V'
  UPLO = 'L'  

  !  write(*,*) N
  ALLOCATE(H_AUX_(N,N))
  H_AUX_ = 0.0
  LDA   =  N
  LWORK = 2*N

  ALLOCATE(WORK(2*N))
  ALLOCATE(RWORK(3*N-2))

  !---- use zheev to get the optimun value of LWORK
  WORK = -1
  CALL ZHEEV(JOBZ,UPLO,N,H_AUX_,LDA,W_SPACE,WORK,LWORK,RWORK,INFO)
  LWORK = INT(WORK(1))
  DEALLOCATE(WORK)
  ALLOCATE(WORK(LWORK))
  IF(INFO /= 0) WRITE(*,*) "# DIAG FAIL 0",INFO

  CALL ZHEEV(JOBZ,UPLO,N,H,LDA,W_SPACE, WORK, LWORK, RWORK,INFO)       
  IF(INFO /= 0) WRITE(*,*) "# DIAG FAIL 1",INFO

  DEALLOCATE(WORK)
  DEALLOCATE(RWORK) 
  DEALLOCATE(H_AUX_)
END SUBROUTINE LAPACK_FULLEIGENVALUES


SUBROUTINE LAPACK_SELECTEIGENVALUES(H,N,W_SPACE,L1,L2,Z,INFO)
!selected eigenvalues/vectors of hermitian matrix
!H, inout, packed banded matrix
! , out,eigenvectors
!N, in,matrix dimension
!W_space, out, eigenvalues
!L1 ordinal lowest eigenvalue
!L2 ordinal highest eigenvlaue
!Z : eigenvectors
!INFO,inout, error flag

  !USE FLOQUET
  IMPLICIT NONE
  INTEGER,                        INTENT(IN)    :: N,L1,L2
  COMPLEX*16, DIMENSION(:,:),     INTENT(INOUT) :: H
  COMPLEX*16, DIMENSION(:,:),     INTENT(OUT)   :: Z
  DOUBLE PRECISION, DIMENSION(:), INTENT(OUT)   :: W_SPACE
  INTEGER,                        INTENT(OUT)   :: INFO


  !---SETTING  LAPACK VARIABLES: START ---------!
  CHARACTER         JOBZ, UPLO, RANGE
  DOUBLE PRECISION ABSTOL, VL, VU
  INTEGER           LWORK,LDA,LRWORK, LIWORK, IL, IU, MP,LDZ,LWMAX
  INTEGER,          DIMENSION(:),   ALLOCATABLE :: IWORK,ISUPPZ
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: RWORK
  COMPLEX*16,       DIMENSION(:),   ALLOCATABLE :: WORK
  COMPLEX*16,       DIMENSION(:,:), ALLOCATABLE :: H_AUX_  
  JOBZ = 'V'
  UPLO = 'L'
  RANGE = 'A'
  LWMAX  = 2000
  ABSTOL = -1.0

  ALLOCATE(H_AUX_(N,N))
  H_AUX_ = 0.0
  LDA   =   N
  LWORK = 2*N
  LDZ   =   N
  ALLOCATE(WORK(LWMAX))
  ALLOCATE(RWORK(LWMAX))
  ALLOCATE(IWORK(LWMAX))
  ALLOCATE(ISUPPZ(N))
  LWORK = -1
  LRWORK = -1
  LIWORK = -1      
  VL = 0.749
  VU = 0.7503! 45.0
  IL = 1!L1!(2*N_Floquet+1)*(3) + 5*(int(0.5*(2*N_Floquet+1)) - SUBSPACES) + 1
  IU = 1!L2!IL + (2*SUBSPACES + 1)*5 - 1
  !  write(*,*) IL,IU,N
!!$  IF(JOBZ.EQ.'V')  THEN
!!$     ALLOCATE(Z(N,IU-IL + 1))
!!$  ELSE
!!$     ALLOCATE(Z(1,1))
!!$  END IF

  CALL ZHEEVR( JOBZ, RANGE, 'L', N, H_AUX_, LDA, VL, VU, IL, &
       &   IU, ABSTOL, MP, W_SPACE, Z, LDZ, ISUPPZ, WORK, LWORK, RWORK, &
       &   LRWORK, IWORK, LIWORK, INFO )
  IF(INFO /= 0) WRITE(*,*) "DIAG FAIL"
  LWORK  = WORK(1)
  LRWORK = RWORK(1)
  LIWORK = IWORK(1)
  DEALLOCATE(WORK,RWORK,IWORK)
  ALLOCATE(WORK(LWORK))
  ALLOCATE(RWORK(LRWORK))
  ALLOCATE(IWORK(LIWORK))
  !write(*,*)INFO,LWORK,LRWORK,LIWORK
  CALL ZHEEVR( JOBZ, RANGE, 'L', N, H, LDA, VL, VU, IL, &
       &   IU, ABSTOL, MP, W_SPACE, Z, LDZ, ISUPPZ, WORK, LWORK, RWORK, &
       &   LRWORK, IWORK, LIWORK, INFO )
  IF(INFO /= 0) WRITE(*,*) "DIAG FAIL"
  !  write(*,*) info

  DEALLOCATE(WORK,RWORK,IWORK,H_AUX_,ISUPPZ)


END SUBROUTINE LAPACK_SELECTEIGENVALUES

!!$SUBROUTINE MKLSPARSE_FULLEIGENVALUES(D,DV,VALUES,ROW_INDEX,COLUMN,E_L,E_R,E_FLOQUET,U_F,INFO)
!!$
!!$!CALCULATES THE ENERGY SPECTRUM OF THE MATRIX REPRESENTED BY VALUES, ROW_INDEX AND COLUMN
!!$! D (IN), MATRIX DIMENSION == NUMBER OF EIGENVALUES
!!$! DV (IN), NUMBER OF VALUES != 0
!!$! VALUES (IN) ARRAY OF VALUES
!!$! ROW_INDEX (IN), ARRAY OF INDICES
!!$! COLUMN (IN),    ARRAY OF COLUMN NUMBERS
!!$! E_L (IN),       LEFT BOUNDARY OF THE SEARCH INTERVAL
!!$! E_R (IN),       RIGHT BOUNDARY OF THE SEARCH INTERVAL
!!$! E_FLOQUET (OUT), ARRAY OF EIGENVALUES
!!$! INFO     (INOUT)  ERROR FLAG
!!$
!!$  USE FEAST
!!$  IMPLICIT NONE
!!$  INTEGER,                          INTENT(IN)    :: D,DV
!!$  COMPLEX*16,       DIMENSION(DV),  INTENT(INOUT) :: VALUES
!!$  INTEGER,          DIMENSION(DV),  INTENT(INOUT) :: COLUMN
!!$  INTEGER,          DIMENSION(D+1), INTENT(INOUT) :: ROW_INDEX
!!$  DOUBLE PRECISION,                 INTENT(IN)    :: E_L,E_R
!!$  DOUBLE PRECISION, DIMENSION(D),   INTENT(OUT)   :: E_FLOQUET
!!$  COMPLEX*16,       DIMENSION(D,D), INTENT(OUT)   :: U_F
!!$  INTEGER,                          INTENT(INOUT) :: INFO
!!$
!!$  CHARACTER*1 UPLO
!!$!  DOUBLE PRECISION :: Emin,Emax
!!$  
!!$
!!$  ! ----- 2. RUN FEASTINIT
!!$  
!!$  CALL feastinit(fpm)
!!$  
!!$  ! ----- 3. SOLVE THE STANDARD EIGENVALUE PROBLEM
!!$
!!$  M0  = D  ! number of eignevalues requested
!!$  
!!$  ALLOCATE(E(M0)) ! array of eigenvalues
!!$  ALLOCATE(RES(M0)) ! array of residuals
!!$  ALLOCATE(X(D,M0)) ! matrix with eigenvectors
!!$  E   = 0
!!$  RES = 0
!!$  X   = 0
!!$
!!$
!!$  info_FEAST = 0
!!$  Emin = E_L!-15.0 ! SEARCH INTERVAL: LOWER BOUND
!!$  Emax = E_R! 15.0 ! SEARCH INTERVAL: UPPER BOUND
!!$  UPLO = 'F'
!!$!  WRITE(*,*) d,SIZE(VALUES,1),SIZE(ROW_INDEX,1),SIZE(COLUMN,1),SIZE(FPM,1)
!!$  CALL zfeast_hcsrev(UPLO,D,VALUES,ROW_INDEX,COLUMN,fpm,epsout,loop, &
!!$       &   Emin,Emax,M0,E,X,M1,res,info_FEAST)
!!$  !PRINT  *,'FEAST OUTPUT INFO ',info
!!$  !WRITE(*,*) 'GUESSED NUMBER OF EIGENVALUES:', m0
!!$  !WRITE(*,*) 'NUMBER OF EIGENVALUES FOUND:', m1
!!$  !WRITE(*,*) 'EIGENVALUES:', E, D
!!$  !WRITE(*,*) 'EIGENVECTORS:', ABS(X)
!!$  !CALL WRITE_MATRIX(D_MULTIFLOQUET,D_MULTIFLOQUET,ABS(X))
!!$  E_FLOQUET = E
!!$  U_F       = X
!!$
!!$  DEALLOCATE(RES)
!!$  DEALLOCATE(X)
!!$  DEALLOCATE(E)
!!$  
!!$END SUBROUTINE MKLSPARSE_FULLEIGENVALUES
!!$
