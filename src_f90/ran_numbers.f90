MODULE RAN_NUMBERS ! new wersion with the elemental subroutione srun1()
  !$ USE OMP_LIB
  IMPLICIT NONE
  PUBLIC ::  RAN_INIT, RAN1, SRAN1, GRAN1, IRAN1,MRAN1,NRAN1, RAN_SAVE,OMP_RAN_INIT,GAU_OMP,BYTECOUNT
  REAL ,PRIVATE,PARAMETER:: R12=3.4641016151377544 !sqrt(12)

  INTEGER(4), PRIVATE, PARAMETER :: HG=HUGE(1_4), HGM=-HG, HGNG=HGM-1
  PRIVATE :: RAN_HASH
  TYPE RAN_SAVE
     SEQUENCE
     INTEGER(4) :: I,J,K,N,M,R
  END TYPE RAN_SAVE
  TYPE (RAN_SAVE), ALLOCATABLE,TARGET,PUBLIC  :: OMPSEED(:)


  REAL(4),PRIVATE :: AMM, XMM,GXMM
  REAL(8),PRIVATE :: DAMM,DXAMM

CONTAINS

  SUBROUTINE RAN_INIT(S,X)
    TYPE (RAN_SAVE) , INTENT(INOUT) :: X
    INTEGER(4), INTENT(IN) ::S

    INTEGER(4) :: J,HGT,RANSEEDS(5)
    HGT=HG
    !The following lines check that kind value 4 is in fact a 32-bit integer with the usual properties
    !that we expect it to have (under negation and wrap-around addition).If all of these tests are
    !satis bed,then the routines that use this module are portable,even though they go eyond
    !Fortran 90 ’s integer model.
    IF (HG /= 2147483647) PRINT*, 'ran_init: arith assump 1 fails'
    IF (HGNG >= 0) PRINT*,'ran_init: arith assump 2 fails'
    IF (HGT+1 /= HGNG) PRINT*,'ran_init: arith assump 3 fails'
    IF (NOT(HG) >= 0) PRINT*,'ran_init: arith assump 4 fails'
    IF (NOT(HGNG) < 0) PRINT*,'ran_init: arith assump 5 fails'
    IF (HG+HGNG >= 0) PRINT*,'ran_init: arith assump 6 fails'
    IF (NOT(-1_4) < 0) PRINT*,'ran_init: arith assump 7 fails'
    IF (NOT(0_4) >= 0) PRINT*,'ran_init: arith assump 8 fails'
    IF (NOT(1_4) >= 0) PRINT*,'ran_init: arith assump 9 fails'

    AMM=NEAREST(1.0_4,-1.0_4)/HGNG
    DAMM=NEAREST(1.0_8,-1.0_8)/HGNG
    DXAMM=NEAREST(.5_8,-1.0_8)/HGNG
    XMM=NEAREST(.5_4,-1.0_4)/(REAL(HG,4))
    GXMM=XMM*R12
    IF (AMM*HGNG >= 1.0 .OR. AMM*HGNG <= 0.0) PRINT*,'ran_init: arth assump 10 fails'
    RANSEEDS(1)=S
    RANSEEDS(2:5)=(/10,20,30,40/)

    DO J=1,4 
       CALL RAN_HASH(RANSEEDS(J),RANSEEDS(J+1))
    END DO
    WHERE (RANSEEDS(1:3) < 0) RANSEEDS(1:3)=NOT(RANSEEDS(1:3)) !Enforce nonnegativity.
       WHERE (RANSEEDS(4:5) == 0) RANSEEDS(4:5)=1 !Enforce nonzero.

       X%I=RANSEEDS(1)
       X%J=RANSEEDS(2)
       X%K=RANSEEDS(3)
       X%M=RANSEEDS(4)
       X%N=RANSEEDS(5)
       X%R=X%N


     END SUBROUTINE RAN_INIT

     SUBROUTINE RAN_HASH(IL,IR)
       IMPLICIT NONE
       INTEGER(4), INTENT(INOUT) :: IL,IR
       !DES-like hashing of two 32-bit integers,using shifts,xor ’s,and adds to make the internal nonlinear function.
       INTEGER(4) :: IS,J
       DO J=1,4
          IS=IR
          IR=IEOR(IR,ISHFT(IR,5))+1422217823     ! THE VARIOUS CONSTANTS ARE CHOSEN TO GIVE GOOD BIT 
          IR=IEOR(IR,ISHFT(IR,-16))+1842055030   ! MIXING AND SHOULD NOT BE CHANGED.
          IR=IEOR(IR,ISHFT(IR,9))+80567781
          IR=IEOR(IL,IR)
          IL=IS
       END DO
     END SUBROUTINE RAN_HASH

     REAL(4) FUNCTION  RAN1(X) RESULT(R) 
       TYPE (RAN_SAVE) , INTENT(INOUT) :: X
       !REAL(4), INTENT(OUT) :: R

       X%R=X%I-X%K !Update Fibonacci generator,which has period p^2+p+1, p=2^31-69
       IF (X%R < 0) X%R=X%R+2147483579_4
       X%I=X%J
       X%J=X%K
       X%K=X%R
       X%N=IEOR(X%N,ISHFT(X%N,13)) !Update Marsaglia shift sequence.
       X%N=IEOR(X%N,ISHFT(X%N,-17))
       X%N=IEOR(X%N,ISHFT(X%N,5))!Once only per cycle,advance sequence y 1,shortening its period to 2 32?2

       if (X%N == 1) X%N=270369_4
       X%M=IEOR(X%M,ISHFT(X%M,5))    !Update Marsaglia shift sequence with period 2**32-1
       X%M=IEOR(X%M,ISHFT(X%M,-13))
       X%M=IEOR(X%M,ISHFT(X%M,6))
       X%R=IEOR(X%N,X%R)+X%M !Combine the generators.The a ove statement has wrap-around addition.
       R=AMM*MERGE(X%R,NOT(X%R), X%R<0 ) !Make the result positive definite (note that amm is negative).
     END FUNCTION RAN1


     INTEGER(4) FUNCTION  MRAN1(X,MAX) RESULT(R) ! generates integer random number from 0<=R<max
       INTEGER,INTENT(IN)::MAX 
       TYPE (RAN_SAVE) , INTENT(INOUT) :: X
       !REAL(4), INTENT(OUT) :: R

       X%R=X%I-X%K !Update Fibonacci generator,which has period p^2+p+1, p=2^31-69
       IF (X%R < 0) X%R=X%R+2147483579_4
       X%I=X%J
       X%J=X%K
       X%K=X%R
       X%N=IEOR(X%N,ISHFT(X%N,13)) !Update Marsaglia shift sequence.
       X%N=IEOR(X%N,ISHFT(X%N,-17))
       X%N=IEOR(X%N,ISHFT(X%N,5))!Once only per cycle,advance sequence y 1,shortening its period to 2 32?2

       if (X%N == 1) X%N=270369_4
       X%M=IEOR(X%M,ISHFT(X%M,5))    !Update Marsaglia shift sequence with period 2**32-1
       X%M=IEOR(X%M,ISHFT(X%M,-13))
       X%M=IEOR(X%M,ISHFT(X%M,6))
       X%R=IEOR(X%N,X%R)+X%M !Combine the generators.The a ove statement has wrap-around addition.
       R=FLOOR((DAMM*REAL(MAX,8))*REAL(MERGE(X%R,NOT(X%R),X%R<0),8),4)
       !R=MOD(ABS(X%R),MAX);
     END FUNCTION MRAN1

     INTEGER(4) FUNCTION  NRAN1(X,MAX) RESULT(R) ! generates integer random number from 1<=R<=max
       INTEGER,INTENT(IN)::MAX 
       TYPE (RAN_SAVE) , INTENT(INOUT) :: X
       X%R=X%I-X%K !Update Fibonacci generator,which has period p^2+p+1, p=2^31-69
       IF (X%R < 0) X%R=X%R+2147483579_4
       X%I=X%J
       X%J=X%K
       X%K=X%R
       X%N=IEOR(X%N,ISHFT(X%N,13)) !Update Marsaglia shift sequence.
       X%N=IEOR(X%N,ISHFT(X%N,-17))
       X%N=IEOR(X%N,ISHFT(X%N,5))!Once only per cycle,advance sequence y 1,shortening its period to 2 32?2
       IF (X%N == 1) X%N=270369_4
       X%M=IEOR(X%M,ISHFT(X%M,5))    !Update Marsaglia shift sequence with period 2**32-1
       X%M=IEOR(X%M,ISHFT(X%M,-13))
       X%M=IEOR(X%M,ISHFT(X%M,6))
       X%R=IEOR(X%N,X%R)+X%M !Combine the generators.The a ove statement has wrap-around addition.
       !R=CEILING((DAMM*REAL(MAX,8))*REAL(MERGE(X%R,NOT(X%R),X%R<0),8),4)
       R=CEILING((REAL(MAX,8))*(DXAMM*REAL(X%R,8)+.5_8),4)
     END FUNCTION NRAN1


     ELEMENTAL SUBROUTINE  SRAN1(R,X) ! 0<r<1
       REAL(4),INTENT(OUT) :: R 
       TYPE (RAN_SAVE) , INTENT(INOUT) :: X
       !REAL(4), INTENT(OUT) :: R

       X%R=X%I-X%K !Update Fibonacci generator,which has period p^2+p+1, p=2^31-69
       IF (X%R < 0) X%R=X%R+2147483579_4
       X%I=X%J
       X%J=X%K
       X%K=X%R
       X%N=IEOR(X%N,ISHFT(X%N,13)) !Update Marsaglia shift sequence.
       X%N=IEOR(X%N,ISHFT(X%N,-17))
       X%N=IEOR(X%N,ISHFT(X%N,5))!Once only per cycle,advance sequence y 1,shortening its period to 2 32?2

       if (X%N == 1) X%N=270369_4
       X%M=IEOR(X%M,ISHFT(X%M,5))    !Update Marsaglia shift sequence with period 2**32-1
       X%M=IEOR(X%M,ISHFT(X%M,-13))
       X%M=IEOR(X%M,ISHFT(X%M,6))
       X%R=IEOR(X%N,X%R)+X%M !Combine the generators.The a ove statement has wrap-around addition.
       R=AMM*MERGE(X%R,NOT(X%R), X%R<0 ) !Make the result positive definite (note that amm is negative).
     END SUBROUTINE SRAN1


     ELEMENTAL SUBROUTINE IRAN1(X) 
       TYPE (RAN_SAVE) , INTENT(INOUT) :: X
       !INTEGER(4), INTENT(OUT) :: IR

       X%R=X%I-X%K !Update Fibonacci generator,which has period p^2+p+1, p=2^31-69
       IF (X%R < 0) X%R=X%R+2147483579_4
       X%I=X%J
       X%J=X%K
       X%K=X%R
       X%N=IEOR(X%N,ISHFT(X%N,13)) !Update Marsaglia shift sequence.
       X%N=IEOR(X%N,ISHFT(X%N,-17))
       X%N=IEOR(X%N,ISHFT(X%N,5))!Once only per cycle,advance sequence y 1,shortening its period to 2 32?2

       if (X%N == 1) X%N=270369_4
       X%M=IEOR(X%M,ISHFT(X%M,5))    !Update Marsaglia shift sequence with period 2**32-1
       X%M=IEOR(X%M,ISHFT(X%M,-13))
       X%M=IEOR(X%M,ISHFT(X%M,6))
       X%R=IEOR(X%N,X%R)+X%M !Combine the generators.The a ove statement has wrap-around addition.
     END SUBROUTINE IRAN1

     REAL(4) FUNCTION GRAN1(X) RESULT(RN)
       TYPE (RAN_SAVE) , INTENT(INOUT) :: X
       X%R=X%I-X%K !Update Fibonacci generator,which has period p^2+p+1, p=2^31-69
       IF (X%R < 0) X%R=X%R+2147483579_4
       X%I=X%J
       X%J=X%K
       X%K=X%R
       X%N=IEOR(X%N,ISHFT(X%N,13)) !Update Marsaglia shift sequence.
       X%N=IEOR(X%N,ISHFT(X%N,-17))
       X%N=IEOR(X%N,ISHFT(X%N,5))!Once only per cycle,advance sequence y 1,shortening its period to 2 32?2

       if (X%N == 1) X%N=270369_4
       X%M=IEOR(X%M,ISHFT(X%M,5))    !Update Marsaglia shift sequence with period 2**32-1
       X%M=IEOR(X%M,ISHFT(X%M,-13))
       X%M=IEOR(X%M,ISHFT(X%M,6))
       X%R=IEOR(X%N,X%R)+X%M !Combine the generators.The a ove statement has wrap-around addition.
       RN=XMM*X%R                !Make  -0.5<=RN<=0.5
     END FUNCTION GRAN1

     REAL(4) FUNCTION GRAN_OMP(X) RESULT(RN)
       TYPE (RAN_SAVE) , INTENT(INOUT) :: X
       X%R=X%I-X%K !Update Fibonacci generator,which has period p^2+p+1, p=2^31-69
       IF (X%R < 0) X%R=X%R+2147483579_4
       X%I=X%J
       X%J=X%K
       X%K=X%R
       X%N=IEOR(X%N,ISHFT(X%N,13)) !Update Marsaglia shift sequence.
       X%N=IEOR(X%N,ISHFT(X%N,-17))
       X%N=IEOR(X%N,ISHFT(X%N,5))!Once only per cycle,advance sequence y 1,shortening its period to 2 32?2

       if (X%N == 1) X%N=270369_4
       X%M=IEOR(X%M,ISHFT(X%M,5))    !Update Marsaglia shift sequence with period 2**32-1
       X%M=IEOR(X%M,ISHFT(X%M,-13))
       X%M=IEOR(X%M,ISHFT(X%M,6))
       X%R=IEOR(X%N,X%R)+X%M !Combine the generators.The above statement has wrap-around addition.
       RN=XMM*X%R                !Make  -0.5<=RN<=0.5
     END FUNCTION GRAN_OMP

     REAL(4) FUNCTION GAU_OMP(X) RESULT(RN)
       TYPE (RAN_SAVE) , INTENT(INOUT) :: X
       X%R=X%I-X%K !Update Fibonacci generator,which has period p^2+p+1, p=2^31-69
       IF (X%R < 0) X%R=X%R+2147483579_4
       X%I=X%J
       X%J=X%K
       X%K=X%R
       X%N=IEOR(X%N,ISHFT(X%N,13)) !Update Marsaglia shift sequence.
       X%N=IEOR(X%N,ISHFT(X%N,-17))
       X%N=IEOR(X%N,ISHFT(X%N,5))!Once only per cycle,advance sequence y 1,shortening its period to 2 32?2

       if (X%N == 1) X%N=270369_4
       X%M=IEOR(X%M,ISHFT(X%M,5))    !Update Marsaglia shift sequence with period 2**32-1
       X%M=IEOR(X%M,ISHFT(X%M,-13))
       X%M=IEOR(X%M,ISHFT(X%M,6))
       X%R=IEOR(X%N,X%R)+X%M !Combine the generators.The a ove statement has wrap-around addition.
       RN=GXMM*X%R                !Make  <RN^2>=1
     END FUNCTION  GAU_OMP


     SUBROUTINE OMP_RAN_INIT(NTR)
       INTEGER, INTENT(OUT) :: NTR
       INTEGER :: I,J,TNR
       NTR=0
       TNR=0
       !$OMP  PARALLEL DEFAULT(SHARED) PRIVATE(TNR)
       !$OMP MASTER
       NTR=OMP_GET_NUM_THREADS()-1
       !$OMP END MASTER
       TNR=OMP_GET_THREAD_NUM()
       !PRINT*,'THREAD_NUM',TNR, ' OF ', NTR 
       !$OMP  END PARALLEL
       ALLOCATE(OMPSEED(0:NTR))
       DO I=0,NTR
          CALL SYSTEM_CLOCK (COUNT=J)
          J=MAX((J+16290047*I),16290047*(I+1))
          CALL RAN_INIT(J,OMPSEED(I))
       END DO
     END SUBROUTINE OMP_RAN_INIT
     !==================================================================================
     INTEGER(1) ELEMENTAL FUNCTION BYTECOUNT(IB) RESULT(N)
       INTEGER(1),INTENT(IN) ::IB
       INTEGER :: I
       N=1_1.AND.IB
       do i=1,7
          N=N+(1_1.AND.RSHFT(IB,I))
       end do
     END FUNCTION BYTECOUNT

   END MODULE RAN_NUMBERS







