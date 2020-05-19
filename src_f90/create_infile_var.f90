!--- To generate pair coeff files for CG MC--------------------------
!--- Input required temperature -------------------------------------
!--- Version: Nov-29-2017--------------------------------------------
!********************************************************************

MODULE INPPARAMS

  IMPLICIT NONE
  
  INTEGER, PARAMETER :: N = py_nmons
  INTEGER, PARAMETER :: natomtypes = 8
  REAL, PARAMETER :: graftdens = py_graftdens
  INTEGER, PARAMETER :: graft_opt = py_graftopt
  INTEGER, PARAMETER :: key_eps_polygraft = py_key_epspg
  REAL, PARAMETER :: eps_poly_graft = py_eps_pg
  REAL, PARAMETER :: eps_graft_graft = py_epsval

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
  REAL :: sig_pg, rc_pg

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


  IF(graftdens .NE. 0) THEN
     siggraft = 0.0; rcutgraft = 0.0
     DO i = 1,natomtypes
        siggraft  = siggraft + sigma(i)
        rcutgraft = rcutgraft + rcut(i)
     END DO
     
     IF(graft_opt == 1) THEN
        
        siggraft  = siggraft/REAL(natomtypes)
        rcutgraft = rcutgraft/REAL(natomtypes)
        
     ELSEIF(graft_opt == 2) THEN
        
        siggraft  = (0.67)*siggraft/REAL(natomtypes) 
        rcutgraft = (0.67)*rcutgraft/REAL(natomtypes)
        
     END IF

     WRITE(20,"(A,1X,2(I0,1X),3(F14.6,1X))") "pair_coeff"&
          &,natomtypes+1, natomtypes+1,eps_graft_graft,siggraft&
          &,rcutgraft

  END IF

  
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
  
  CLOSE(20)

  ! Prof Ilja's suggestion. If not just add a comment line in the
  ! pairfile so that when it is read from LAMMPS infile, it will be
  ! ignored. Add it after LB rules are applied 

  OPEN(unit=22,file="in.cgpair_polygraft_MC",action="write",status="re&
       &place",iostat=ierr)

  IF(ierr /= 0) STOP "Could not open in.cgpair_polygraft_MC"

  IF(key_eps_polygraft == 1 .AND. graftdens .NE. 0) THEN

     DO i = 1, natomtypes

        sig_pg = sqrt(siggraft*sigma(i))
        rc_pg  = sqrt(rcutgraft*rcut(i))
        WRITE(22,"(A,1X,2(I0,1X),3(F14.6,1X))") "pair_coeff"&
             &,i, natomtypes+1,eps_poly_graft,sig_pg,rc_pg

     END DO

  ELSE
     
     WRITE(22, "(A80)") "## Dummy file: No data needed from this file"

  END IF

END SUBROUTINE FIND_GENERATE_AND_WRITE_PAIRCOEFF

!--------------------------------------------------------------------
