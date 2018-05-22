SUBROUTINE MULTIMODEFLOQUETMATRIX(ATOM_,NM,NF,MODES_NUM,FIELD,INFO)
  !ID,size(modes_num,1),total_frequencies,MODES_NUM,FIELDS,INFO
  !  USE FLOQUET
  !ATOM_ type atom, -> dimension of the bare Hilbert space
  !NM -> number of modes
  !NF -> Number of Fields
  !MODES_NUM -> number of harmonics of each mode
  !FIELD -> Field couplings
  !INFO

  USE ARRAYS
  USE ATOMIC_PROPERTIES
  USE TYPES
  USE SUBINTERFACE_LAPACK

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NM,NF
  INTEGER,                  INTENT(INOUT) :: INFO
  INTEGER,   DIMENSION(NM), INTENT(IN)     :: MODES_NUM
  TYPE(MODE),DIMENSION(NF), INTENT(IN)    :: FIELD
  TYPE(ATOM),               INTENT(IN)    :: ATOM_                       

  INTEGER m,n,D,r,o
  COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: H_TEMP,H_STATIC,COUPLING,Z_M_COPY
  DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE :: E_DRESSED
  INTEGER index_l_lower,index_l_upper,index_r_lower,index_r_upper

  INTEGER, DIMENSION(NF) :: N_FLOQUET

  N_FLOQUET = 0
  DO n=2,NF
     N_FLOQUET(n)=FIELD(n)%N_Floquet
  END DO

  D        = ATOM_%D_BARE
  ALLOCATE(H_FLOQUET_COPY(D,D))

  H_FLOQUET_COPY = FIELD(1)%V  ! STATIC HAMILTONIAN


  DO n=2,NM  ! RUN OVER EACH FREQUENCY

     ! D : UPDATED AT THE ENDO OF THE LOOP. DIMENSION OF THE MULTIMODE FLOQUET MATRIX

     H_STATIC  = H_FLOQUET_COPY  
     DEALLOCATE(H_FLOQUET_COPY)

     ALLOCATE(IDENTITY(D,D))
     IDENTITY  = 0.0
     DO m= 1,D 
        IDENTITY(m,m) = 1.0
     END DO

     ALLOCATE(COUPLING(D,D))
     COUPLING  = 0.0

     DO r=1,(2*N_Floquet(n-1)+1)

        index_l_lower = ATOM_%D_BARE*(r - 1) + 1
        index_l_upper = ATOM_%D_BARE*(r - 1) + ATOM_%D_BARE
        index_r_lower = index_l_lower
        index_r_upper = index_l_upper
        COUPLING(index_l_lower:index_l_upper, index_r_lower:index_r_upper) = &
             &     FIELD(n)%V  ! COUPLING MATRIX OF MODE n
     END DO

     D = D*(2*N_FLOQUET(n)+1)
     ALLOCATE(H_FLOQUET(D,D))
     H_FLOQUET = 0.0


     DO m=-N_floquet(n),N_Floquet(n)

        index_l_lower = (m + N_Floquet(n)    )*SIZE(COUPLING,1) + 1
        index_l_upper = index_l_lower + SIZE(COUPLING,1) - 1
        index_r_lower =  index_l_lower
        index_r_upper =  index_l_upper
        H_FLOQUET(index_l_lower:index_l_upper, index_r_lower:index_r_upper) = &
             &  1.0*H_STATIC + 1.0*m*1.0*FIELD(n)%OMEGA*IDENTITY/A
        IF(m.LT.N_floquet(n)) THEN

           index_l_lower =  (m + N_Floquet(n) + 1)*SIZE(COUPLING,1) + 1
           index_l_upper =  index_l_lower + SIZE(COUPLING,1) - 1
           index_r_lower =  (m + N_Floquet(n)    )*SIZE(COUPLING,1) + 1
           index_r_upper =  index_r_lower + SIZE(COUPLING,1) - 1
           H_FLOQUET(index_l_lower:index_l_upper, index_r_lower:index_r_upper) = &
                &     0.5*COUPLING

           index_l_lower =  (m + N_Floquet(n)    )*SIZE(COUPLING,1) + 1
           index_l_upper =  index_r_lower + SIZE(COUPLING,1)  - 1        
           index_r_lower =  (m + N_Floquet(n) + 1)*SIZE(COUPLING,1) + 1
           index_r_upper =  index_l_lower + SIZE(COUPLING,1)  - 1         
           H_FLOQUET(index_l_lower:index_l_upper, index_r_lower:index_r_upper) = &
                &     0.5*TRANSPOSE(CONJG(COUPLING))
        END IF

     END DO

     DEALLOCATE(IDENTITY)
     DEALLOCATE(COUPLING)

     IF(n.LT.NF) THEN
        ALLOCATE(H_FLOQUET_COPY(D,D))
        H_FLOQUET_COPY = H_FLOQUET
        DEALLOCATE(H_FLOQUET)
     END IF
  END DO
END SUBROUTINE MULTIMODEFLOQUETMATRIX

