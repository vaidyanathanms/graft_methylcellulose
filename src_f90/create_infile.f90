!--- To generate pair coeff files for CG MC--------------------------
!--- Input required temperature -------------------------------------
!--- Version: Nov-29-2017--------------------------------------------
!********************************************************************

MODULE INPPARAMS

  IMPLICIT NONE
  
  INTEGER, PARAMETER :: N = 1000
  INTEGER, PARAMETER :: natomtypes = 8
  INTEGER, PARAMETER :: grafts = 2
  REAL, PARAMETER :: epsgraft = 1.0

  REAL, DIMENSION(1:natomtypes) :: epsval, sigma, rcut
  REAL, PARAMETER :: bondcoeff = 500
  REAL, PARAMETER :: bondeq = 1.0
  REAL, PARAMETER :: anglcoeff =15
  REAL, PARAMETER :: thetaeq = 165
  REAL, PARAMETER :: dihdcoeff = 2.0
  INTEGER, PARAMETER :: neq = 1
  INTEGER, PARAMETER :: deq = 1

  REAL :: siggraft, rcutgraft
  REAL :: tempval

END MODULE INPPARAMS

!--------------------------------------------------------------------

PROGRAM CREATE_PAIRCOEFF_FILE

  USE INPPARAMS
  IMPLICIT NONE
  
  CALL FINDTEMP()
  CALL FIND_GENERATE_AND_WRITE_PAIRCOEFF()
  
  PRINT *, "Calculations complete .."

END PROGRAM CREATE_PAIRCOEFF_FILE

!--------------------------------------------------------------------

SUBROUTINE FINDTEMP()

  USE INPPARAMS
  IMPLICIT NONE
  
  INTEGER :: narg
  CHARACTER(LEN=10) :: tval

  narg = IARGC()
  IF(narg .NE. 1) STOP "Unknown number of extensions"

  CALL GETARG(narg,tval)

  READ(tval,"(F5.2)") tempval

  PRINT *, "Input temperature : ", tempval

END SUBROUTINE FINDTEMP

!--------------------------------------------------------------------

SUBROUTINE FIND_GENERATE_AND_WRITE_PAIRCOEFF()

  USE INPPARAMS
  IMPLICIT NONE
  
  INTEGER :: i,ierr,u, flag
  CHARACTER(LEN=10) :: tval
  REAL :: dumval, aij, bij

  flag = 0
  OPEN(unit = 10,file="cgparams.txt",action="read",status="old"&
       &,iostat=ierr)

  IF(ierr .NE. 0) STOP "cgparams.txt not found"

  sigma = -1; rcut = -1; epsval = -1

  DO

     READ(10,*) dumval

     IF(dumval == tempval) THEN

        DO i = 1,natomtypes
           
           READ(10,*) u, aij, bij, sigma(u), rcut(u)
           
           epsval(u) = aij*(N**(-bij))
           
        END DO

        flag = 1
        EXIT
        
     END IF

  END DO

  IF(flag .NE. 1) STOP "Reqd temperature not found"

  CLOSE(10)

  OPEN(unit=20,file="in.cgpair_MC",action="write",status="replace"&
       &,iostat=ierr)

  IF(ierr /= 0) STOP "Could not open in.cgpair_MC"

  WRITE(20,'(A60,1X,F12.4)') "# Pair Coefficients for CGMC simulation&
       &s at",tempval

  WRITE(20,*)
  WRITE(20,'(A)') " # Pair_Coeffs" 
  DO i = 1,natomtypes

     WRITE(20,"(A,1X,2(I0,1X),3(F14.6,1X))") "pair_coeff",i,i&
          &,epsval(i),sigma(i),rcut(i)

  END DO



  siggraft = 0.0; rcutgraft = 0.0
  DO i = 1,natomtypes
     siggraft  = siggraft + sigma(i)
     rcutgraft = rcutgraft + rcut(i)
  END DO
  
  IF(grafts == 1) THEN
   
     siggraft  = siggraft/REAL(natomtypes)
     rcutgraft = rcutgraft/REAL(natomtypes)

  ELSEIF(grafts == 2) THEN

     siggraft  = (0.67)*siggraft/REAL(natomtypes) 
     rcutgraft = (0.67)*rcutgraft/REAL(natomtypes)
     
  END IF

  WRITE(20,"(A,1X,2(I0,1X),3(F14.6,1X))") "pair_coeff"&
       &,natomtypes+1, natomtypes+1,epsgraft,siggraft,rcutgraft


  WRITE(20,*)     
  WRITE(20,'(A)') " # Bond_Coeffs" 
  WRITE(20,'(A,1X,2(F12.8,1X))') "bond_coeff * ",bondcoeff,bondeq

  WRITE(20,*)     
  WRITE(20,'(A)') " # Angle_Coeffs" 
  WRITE(20,'(A,1X,2(F12.8,1X))') "angle_coeff * ",anglcoeff,thetaeq

  WRITE(20,*)     
  WRITE(20,'(A)') " # Dihedral_Coeffs" 
  WRITE(20,'(A,1X,F12.8,1X,2(I0,1X))') "dihedral_coeff * ",dihdcoeff&
       &,deq,neq
  

END SUBROUTINE FIND_GENERATE_AND_WRITE_PAIRCOEFF

!--------------------------------------------------------------------
