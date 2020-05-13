!-----------Use this in conjunction with SMD of LAMMPS---------------
!---------------To extract files for US after SMD--------------------
!---------------Version 2: Apr-23-2018-------------------------------
!---------------Param File: extract_params.f90-----------------------
!********************************************************************

PROGRAM EXTRACTCONF

  USE EXTRACT_PARAMS
  IMPLICIT NONE

  PRINT *, "Static analysis to generate configuration .."
  PRINT *, "Starting OMP Threads .."

!$OMP PARALLEL
  nproc = OMP_GET_NUM_THREADS()
  PRINT *, "Number of threads are: ", nproc
!$OMP END PARALLEL

  CALL READ_ANA_IP_FILE()
  CALL READ_DATAFILE()
  CALL ANALYZE_TRAJECTORYFILE()
  CALL OUTPUT_ALL()

END PROGRAM EXTRACTCONF

!--------------------------------------------------------------------

SUBROUTINE READ_ANA_IP_FILE()

  USE EXTRACT_PARAMS

  IMPLICIT NONE
  
  INTEGER :: nargs,ierr,logflag,AllocateStatus,i,j
  CHARACTER(256) :: dumchar

  CALL DEFAULTVALUES()

  nargs = IARGC()
  IF(nargs .NE. 1) STOP "Input incorrect"

  logflag = 0

  CALL GETARG(nargs,ana_fname)

  OPEN(unit = anaread,file=trim(ana_fname),action="read",status="old"&
       &,iostat=ierr)
  
  IF(ierr /= 0) THEN

     PRINT *, trim(ana_fname), "not found"
     STOP

  END IF

  DO

     READ(anaread,*,iostat=ierr) dumchar

     IF(ierr .LT. 0) EXIT

     IF(dumchar == 'datafile') THEN
        
        READ(anaread,*,iostat=ierr) data_fname

     ELSEIF(dumchar == 'trajectory_file') THEN

        READ(anaread,*,iostat=ierr) traj_fname

     ELSEIF(dumchar == 'nframes') THEN

        READ(anaread,*,iostat=ierr) nframes
     
     ELSEIF(dumchar == 'group1') THEN

        READ(anaread,*,iostat=ierr) g1_init, g1_fin
        g1_nmons = g1_fin - g1_init + 1

        PRINT *, "Number of monomers in group 1: ", g1_nmons

        ALLOCATE(grp1_array(g1_nmons),stat=AllocateStatus)
        IF(AllocateStatus /= 0) STOP "Did not allocate grp1_array"

     ELSEIF(dumchar == 'group2') THEN

        READ(anaread,*,iostat=ierr) g2_init, g2_fin
        g2_nmons = g2_fin - g2_init + 1

        PRINT *, "Number of monomers in group 1: ", g2_nmons

        ALLOCATE(grp2_array(g2_nmons),stat=AllocateStatus)
        IF(AllocateStatus /= 0) STOP "Did not allocate grp2_array"

     ELSEIF(dumchar == 'multiconfig') THEN

        READ(anaread,*,iostat=ierr) savedist, tol, pullaxis,comflag
        multiposflag = 1

     ELSEIF(dumchar == 'singleconfig') THEN

        READ(anaread,*,iostat=ierr) savedist, tol, pullaxis,comflag&
             &,periodicity
        multiposflag = 0; netgenflag = 0

     ELSEIF(dumchar == 'log_file') THEN

        READ(anaread,*,iostat=ierr) log_fname
        logflag  = 1

     ELSE
        
        PRINT *, "unknown keyword: ", trim(dumchar)
        STOP

     END IF

  END DO

  IF(pullaxis == -1 .OR. pullaxis .GT. 3) THEN
     PRINT *, "ERROR: Unidentified pull axis", pullaxis
     STOP
  END IF

  IF(logflag == 0) log_fname = "logextract."//trim(adjustl(traj_fname))
  OPEN(unit = logout,file=trim(log_fname),action="write",status="repla&
       &ce",iostat=ierr)

  IF(comflag) THEN
     dum_fname = "com_"//trim(adjustl(traj_fname))
     OPEN(unit = comout,file=trim(dum_fname),action="write",status="re&
          &place",iostat=ierr)
  END IF


  PRINT *, "Analysis input file read finished .."

END SUBROUTINE READ_ANA_IP_FILE

!--------------------------------------------------------------------

SUBROUTINE DEFAULTVALUES()

  USE EXTRACT_PARAMS
  IMPLICIT NONE

  ! Frame, Molecules and Processor Details
  nframes = 0; skipfr = 0; freqfr = 1; nfrcntr = 0

  ! Flags

  multiposflag = 0; comflag = 0; netgenflag = 0

  ! Initialize variables

  g1_init = 0; g1_fin = 0; g2_init =0; g2_fin = 0; g1_nmons = 0
  g2_nmons =0; savedist = 0; tol = 0  
  pullaxis = -1
  periodicity = 0 ! if NO pbc then periodicity = 0, else = 1

END SUBROUTINE DEFAULTVALUES

!--------------------------------------------------------------------

SUBROUTINE READ_DATAFILE()

  USE EXTRACT_PARAMS

  IMPLICIT NONE

  INTEGER :: i,j,k,ierr,u,AllocateStatus,imax
  INTEGER :: flag, cntr, nwords
  INTEGER :: aid,molid,atype,ix,iy,iz
  REAL    :: charge,rx,ry,rz
  REAL    :: xlo,xhi,ylo,yhi,zlo,zhi
  CHARACTER(256) :: rline,dumchar

  CALL COMPUTE_INIT_NLINES(imax)

  OPEN(unit=inpread,file = trim(data_fname),action =&
       & 'read', status='old',iostat=ierr) 
  
  IF(ierr .NE. 0) STOP "Data file not found"

  WRITE(logout,*) "Datafile used is :", trim(adjustl(data_fname))

  ntotatoms = 0;ntotbonds=0;ntotangls=0;ntotdihds=0;ntotimprs=0
  atomflag =0;velflag = 0;bondflag=0;anglflag=0;dihdflag=0;imprflag=0

  READ(inpread,*)
  READ(inpread,*)

  DO i = 1,imax-2 
       
     READ(inpread,*) u, dumchar
     
        IF(dumchar == "atoms") THEN
           ntotatoms = u
        ELSEIF(dumchar == "bonds") THEN
           ntotbonds = u
        ELSEIF(dumchar == "angles") THEN
           ntotangls = u
        ELSEIF(dumchar == "dihedrals") THEN
           ntotdihds = u
        ELSEIF(dumchar == "atom" .OR. dumchar == "atomtypes") THEN
           ntotatomtypes = u
        ELSEIF(dumchar == "bond" .OR. dumchar == "bondtypes") THEN
           ntotbondtypes = u
        ELSEIF(dumchar == "angle" .OR. dumchar == "atomtypes") THEN
           ntotangltypes = u
        ELSEIF(dumchar == "dihedral" .OR. dumchar == "dihedraltypes") THEN
           ntotdihdtypes = u
        ELSEIF(dumchar == "improper" .OR. dumchar == "impropertypes") THEN
           ntotimprtypes = u
        ELSEIF(dumchar == "Masses") THEN
           
           ALLOCATE(masses(ntotatomtypes,1),stat = AllocateStatus)
           IF(AllocateStatus/=0) STOP "did not allocate masses"
           
           DO j = 1,ntotatomtypes
              
              READ(inpread,*) u, masses(u,1)
              
           END DO
           
        END IF
        
  END DO

  READ(inpread,*)
  READ(inpread,*) xlo, xhi
  READ(inpread,*) ylo, yhi
  READ(inpread,*) zlo, zhi
  
  box_xl = xhi - xlo
  box_yl = yhi - ylo
  box_zl = zhi - zlo

  PRINT *, "x-box  ", "y-box  ", "z-box  "
  PRINT *, box_xl, box_yl, box_zl

  PRINT *, "STATISTICS"
  PRINT *, "Number of atoms/atomtypes: " , ntotatoms,ntotatomtypes
  PRINT *, "Number of bonds/bondtypes: " , ntotbonds,ntotbondtypes
  PRINT *, "Number of angles/angletypes: " , ntotangls,ntotangltypes
  PRINT *, "Number of diheds/dihedtypes: " , ntotdihds,ntotdihdtypes
  flag = 0; cntr = 0

  CALL ALLOCATE_ARRAYS()

  DO 

     READ(inpread,*,iostat=ierr) dumchar

     IF(ierr .LT. 0) EXIT

     !READ DATA HERE FOR CHARGES AND MOLID
     !READ EVERYTHING AND OVERWRITE LATER
     IF(trim(dumchar) == "Atoms") THEN
             
        atomflag = 1
        print *, "Reading ", trim(dumchar), " info"

        DO j = 1,ntotatoms

           READ(inpread,*) aid,molid,atype,rx,ry,rz

           rx = rx - xlo
           ry = ry - ylo
           rz = rz - zlo

           aidvals(aid,1)     = aid
           aidvals(aid,2)     = molid
           aidvals(aid,3)     = atype
           rxyz_lmp(aid,1)    = rx
           rxyz_lmp(aid,2)    = ry
           rxyz_lmp(aid,3)    = rz

        END DO

     END IF

     IF(trim(dumchar) == "Masses") THEN

        ALLOCATE(masses(ntotatomtypes,1),stat = AllocateStatus)
        IF(AllocateStatus/=0) STOP "did not allocate masses"
         
        DO j = 1,ntotatomtypes

           READ(inpread,*) u, masses(u,1)

        END DO

     END IF

     IF(trim(dumchar) == "Velocities") THEN
             
        velflag = 1
        print *, "Reading ", trim(dumchar), " info"

        DO j = 1,ntotatoms

           READ(inpread,*) vel_xyz(j,1),vel_xyz(j,2),vel_xyz(j,3)&
                &,vel_xyz(j,4)

        END DO


     END IF

     IF(trim(dumchar) == "Bonds") THEN
             
        bondflag = 1
        print *, "Reading ", trim(dumchar), " info"

        DO j = 1,ntotbonds

           READ(inpread,*) bond_lmp(j,1),bond_lmp(j,2),bond_lmp(j,3)&
                &,bond_lmp(j,4)

        END DO

     END IF

     IF(trim(dumchar) == "Angles") THEN
             
        anglflag = 1
        print *, "Reading ", trim(dumchar), " info"

        DO j = 1,ntotangls

           READ(inpread,*) angl_lmp(j,1),angl_lmp(j,2),angl_lmp(j,3)&
                &,angl_lmp(j,4),angl_lmp(j,5)

        END DO

     END IF

     IF(trim(dumchar) == "Dihedrals") THEN
             
        dihdflag = 1
        print *, "Reading", trim(dumchar), "info"

        DO j = 1,ntotdihds

           READ(inpread,*) dihd_lmp(j,1),dihd_lmp(j,2),dihd_lmp(j,3)&
                &,dihd_lmp(j,4),dihd_lmp(j,5),dihd_lmp(j,6)

        END DO

     END IF
  
     IF(trim(dumchar) == "Impropers") THEN
             
        imprflag = 1
        print *, "Reading", trim(dumchar), "info"

        DO j = 1,ntotimprs

           READ(inpread,*) impr_lmp(j,1),impr_lmp(j,2),impr_lmp(j,3)&
                &,impr_lmp(j,4),impr_lmp(j,5),impr_lmp(j,6)

        END DO

     END IF

  END DO
  
  PRINT *, "Fileread finish .."


END SUBROUTINE READ_DATAFILE

!--------------------------------------------------------------------

SUBROUTINE COMPUTE_INIT_NLINES(imax)

  USE EXTRACT_PARAMS

  IMPLICIT NONE

  INTEGER, INTENT(OUT) :: imax
  INTEGER :: init, pos, ipos,u,nwords,lcnt,ierr
  CHARACTER(LEN=120) :: charline

  OPEN(unit=inpread,file = trim(data_fname),action =&
       & 'read', status='old',iostat=ierr) 
  
  IF(ierr .NE. 0) STOP "Data file not found"
  
  lcnt = 0

  READ(inpread,*)

  DO 

     READ(inpread,'(A)',iostat=ierr) charline     

     lcnt = lcnt + 1
     pos = 1
     nwords = 0

     DO

        ipos = VERIFY(charline(pos:),' ')
        IF(ipos == 0) EXIT
        nwords = nwords + 1
        pos = pos + ipos - 1
        ipos = SCAN(charline(pos:),' ')
        IF(ipos == 0) EXIT
        pos = pos + ipos - 1
        
     END DO

     IF(nwords .GE. 4) THEN

        imax = lcnt - 1
        EXIT
        
     END IF

  END DO

  CLOSE(inpread)

END SUBROUTINE COMPUTE_INIT_NLINES

!--------------------------------------------------------------------

SUBROUTINE ANALYZE_TRAJECTORYFILE()

  USE EXTRACT_PARAMS

  IMPLICIT NONE

  INTEGER :: i,j,aid,ierr,atchk,atype,jumpfr
  REAL :: xlo,xhi,ylo,yhi,zlo,zhi
  CHARACTER(LEN = 8) :: dumchar

  OPEN(unit = 15,file =trim(traj_fname),action="read",status="old"&
       &,iostat=ierr)

  IF(ierr /= 0) STOP "trajectory file not found"

  PRINT *, "trajectory file used is :",trim(adjustl(traj_fname))
  WRITE(logout,*) "trajectory file used is :"&
       &,trim(adjustl(traj_fname))


  CALL SORT_GROUPS()

  DO i = 1,skipfr

     DO j = 1,ntotatoms+9

        READ(15,*) 

     END DO

     IF(mod(i,100) == 0) PRINT *, "Skipped ", i, "frames"

  END DO

  i = 1

  DO WHILE(i .LE. nframes)

     IF(netgenflag == 1 .AND. multiposflag == 0) THEN
        PRINT *, "Success at ", i-1, "th frame"
        EXIT
     END IF

     nfrcntr = nfrcntr + 1
     IF(mod(i,100) == 0) PRINT *, "Processing ", i,"th frame"

     READ(15,*)
     READ(15,*) timestep

     READ(15,*) 
     READ(15,*) atchk
!!$     IF(atchk /= ntotatoms) STOP "Number of atoms do not match"

     READ(15,*) 
     READ(15,*) xlo, xhi
     READ(15,*) ylo, yhi
     READ(15,*) zlo, zhi

     READ(15,*)

     box_xl = xhi - xlo
     box_yl = yhi - ylo
     box_zl = zhi - zlo
     
     boxx_arr(i)  = box_xl
     boxy_arr(i)  = box_yl
     boxz_arr(i)  = box_zl

     DO j = 1,atchk

        READ(15,*) aid,atype,rxyz_lmp(aid,1),rxyz_lmp(aid,2)&
             &,rxyz_lmp(aid,3)

        IF(atype .NE. aidvals(aid,3)) THEN

           PRINT *, "Incorrect atom ids"
           PRINT *, i,j,aid,atype,aidvals(aid,3)
           STOP

        END IF

     END DO

     DO j = 1,atchk

        rxyz_lmp(j,1) = rxyz_lmp(j,1) - xlo
        rxyz_lmp(j,2) = rxyz_lmp(j,2) - ylo
        rxyz_lmp(j,3) = rxyz_lmp(j,3) - zlo
        
     END DO

     CALL COMPUTE_COM_DIST(i)

     DO jumpfr = 1,freqfr

        READ(15,*)
        READ(15,*)        
        READ(15,*)
 
        READ(15,*) atchk

        DO j = 1,atchk+5

           READ(15,*) 

        END DO

     END DO
     
     i = i + 1

  END DO

  CLOSE(15)

END SUBROUTINE ANALYZE_TRAJECTORYFILE

!--------------------------------------------------------------------

SUBROUTINE SORT_GROUPS()

  USE EXTRACT_PARAMS
  IMPLICIT NONE

  INTEGER :: i

  DO i = 1, g1_nmons
     
     grp1_array(i) = g1_init + i - 1

  END DO

  IF(grp1_array(g1_nmons) .NE. g1_fin) STOP "Final monomer in g1 did  &
       &not match"

  DO i = 1, g2_nmons
     
     grp2_array(i) = g2_init + i - 1

  END DO

  IF(grp2_array(g2_nmons) .NE. g2_fin) STOP "Final monomer in g1 did  &
       &not match"
    

END SUBROUTINE SORT_GROUPS

!--------------------------------------------------------------------

SUBROUTINE COMPUTE_COM_DIST(iframe)

  USE EXTRACT_PARAMS
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: iframe
  REAL :: xdist, ydist, zdist
  REAL, DIMENSION(2,3) :: comgrp
  REAL, DIMENSION(2) :: totmas_grp
  INTEGER :: aid, atyp,i,j,genflag
  REAL :: dval,absdval,remdval

  comgrp = 0.0; totmas_grp = 0.0; genflag = 0

  ! Compute COM of groups

  DO i = 1,g1_nmons

     aid = grp1_array(i); atyp = aidvals(aid,3) 
     comgrp(1,1) = comgrp(1,1) + rxyz_lmp(aid,1)*masses(atyp,1)
     comgrp(1,2) = comgrp(1,2) + rxyz_lmp(aid,2)*masses(atyp,1)
     comgrp(1,3) = comgrp(1,3) + rxyz_lmp(aid,3)*masses(atyp,1)
     totmas_grp(1) = totmas_grp(1) + masses(atyp,1)
     
  END DO

  DO i = 1,g2_nmons

     aid = grp2_array(i); atyp = aidvals(aid,3) 
     comgrp(2,1) = comgrp(2,1) + rxyz_lmp(aid,1)*masses(atyp,1)
     comgrp(2,2) = comgrp(2,2) + rxyz_lmp(aid,2)*masses(atyp,1)
     comgrp(2,3) = comgrp(2,3) + rxyz_lmp(aid,3)*masses(atyp,1)
     totmas_grp(2) = totmas_grp(2) + masses(atyp,1)
     
  END DO

  DO i = 1,2

     DO j = 1,3

        comgrp(i,j) = comgrp(i,j)/totmas_grp(i)

     END DO

  END DO


  ! Compute Distance

  xdist = comgrp(2,1) - comgrp(1,1)
  ydist = comgrp(2,2) - comgrp(1,2)
  zdist = comgrp(2,3) - comgrp(1,3)


  IF(periodicity == 1) THEN

     !An extra floor command is not necessary - ANINT will take care
     xdist = xdist - box_xl*ANINT(xdist/box_xl)
     ydist = ydist - box_yl*ANINT(ydist/box_yl)
     zdist = zdist - box_zl*ANINT(zdist/box_zl)

  END IF

  ! Whether to extract and write

  IF(pullaxis == 1) THEN
     dval = xdist
  ELSEIF(pullaxis == 2) THEN
     dval = ydist
  ELSEIF(pullaxis == 3) THEN
     dval = zdist
  ELSEIF(pullaxis == 0) THEN
     dval = sqrt(xdist**2 + ydist**2  + zdist**2)
  END IF

  IF(multiposflag) THEN

     absdval = abs(dval)
     remdval = mod(absdval,savedist)

     IF(remdval .LE. tol)THEN

        CALL EXTRACT_AND_WRITE()
        genflag = 1

     END IF


  ELSE

     absdval = abs(dval)
     remdval = absdval - savedist
     
     IF(abs(remdval) .LE. tol) THEN
        
        CALL EXTRACT_AND_WRITE()
        genflag = 1
        netgenflag = 1

     END IF

  END IF


  WRITE(logout,'(2(I0,1X),2(F14.8,1X))') iframe,genflag,absdval&
       &,remdval
  
  IF(comflag) THEN
     
     WRITE(comout,'(I0,1X,7(F14.8,1X))') iframe,comgrp(1,1),comgrp(1&
          &,2),comgrp(1,3),comgrp(2,1),comgrp(2,2),comgrp(2,3),absdval

  END IF

END SUBROUTINE COMPUTE_COM_DIST

!--------------------------------------------------------------------

SUBROUTINE MAKEDIR(dum_dval)

  USE EXTRACT_PARAMS
  IMPLICIT NONE

  REAL, INTENT(IN) :: dum_dval
  CHARACTER (LEN = 6) :: dchar
  CHARACTER (LEN = 30) :: dirname

  WRITE(dchar,'(F6.2)') dum_dval
  dirname = 'mkdir -p '//trim(adjustl('US_dir_'//trim(adjustl(dchar))))

  CALL SYSTEM(dirname)

END SUBROUTINE MAKEDIR

!--------------------------------------------------------------------

SUBROUTINE EXTRACT_AND_WRITE()

  USE EXTRACT_PARAMS
  IMPLICIT NONE

  INTEGER :: ierr,i,aid
  INTEGER, PARAMETER :: nullval = 0
  REAL*8 :: csum
  CHARACTER(LEN=256) :: lmpdatafile
  CHARACTER (LEN = 6) :: dchar

  lmpdatafile = trim(adjustl('init_datafile'))

  OPEN(unit = outfile,file=trim(lmpdatafile),action="write"&
       &,status="replace",iostat=ierr)

  IF(ierr /= 0) STOP "LAMMPS data file not found"


2000 FORMAT(2(F15.8,1X),2X,A)

  WRITE(outfile,*) "Input Configuration for Methylcellulose"
  WRITE(outfile,*)

  WRITE(outfile,"(I0,1X,A5)") ntotatoms, "atoms"
  WRITE(outfile,"(I0,1X,A5)") ntotbonds, "bonds"
  WRITE(outfile,"(I0,1X,A6)") ntotangls, "angles"
  WRITE(outfile,"(I0,1X,A9)") ntotdihds, "dihedrals"
  WRITE(outfile,"(I0,1X,A9)") nullval,"impropers"

  WRITE(outfile,"(I0,1X,A10)") ntotatomtypes,"atom types"
  WRITE(outfile,"(I0,1X,A10)") ntotbondtypes,"bond types"
  WRITE(outfile,"(I0,1X,A11)") ntotangltypes,"angle types"
  WRITE(outfile,"(I0,1X,A14)") ntotdihdtypes,"dihedral types"
  WRITE(outfile,"(I0,1X,A14)") nullval,"improper types"

  
  WRITE (outfile,'(2(F15.8,1X),2X,A)') 0.00, box_xl, "xlo xhi"
  WRITE (outfile,'(2(F15.8,1X),2X,A)') 0.00, box_yl, "ylo yhi"
  WRITE (outfile,'(2(F15.8,1X),2X,A)') 0.00, box_zl, "zlo zhi"

  WRITE(outfile,*)

  WRITE(outfile,*) " Masses"

  WRITE(outfile,*)

  
  DO i = 1, ntotatomtypes

     WRITE(outfile,"(I0,1X,F14.9)") i,masses(i,1)

  END DO

  WRITE(outfile,*)

  WRITE(outfile,*) " Atoms"

  WRITE(outfile,*)

  csum = 0.0

  DO i = 1,ntotatoms
 
     aid = aidvals(i,1)
     WRITE(outfile,'(3(I0,1X),4(F15.9,1X))') aidvals(aid,1)&
          &,aidvals(aid,2),aidvals(aid,3),rxyz_lmp(aid,1)&
          &,rxyz_lmp(aid,2),rxyz_lmp(aid,3)
     
  END DO
     
  IF(ntotbonds /= 0) THEN

     WRITE(outfile,*)
     WRITE(outfile,*) " Bonds"
     WRITE(outfile,*)
     
     DO i = 1,ntotbonds
        
        WRITE(outfile,'(4(I0,1X))') bond_lmp(i,1),bond_lmp(i,2)&
             &,bond_lmp(i,3),bond_lmp(i,4)
        
     END DO
     
  END IF

  IF(ntotangls /= 0) THEN

     WRITE(outfile,*)
     WRITE(outfile,*) " Angles"
     WRITE(outfile,*)

     DO i = 1,ntotangls
     
        WRITE(outfile,'(5(I0,1X))') angl_lmp(i,1),angl_lmp(i,2)&
             &,angl_lmp(i,3),angl_lmp(i,4),angl_lmp(i,5)
     
     END DO

  END IF


  IF(ntotdihds /= 0) THEN

     WRITE(outfile,*)
     WRITE(outfile,*) " Dihedrals"
     WRITE(outfile,*)
     

     DO i = 1,ntotdihds
        
        WRITE(outfile,'(6(I0,1X))') dihd_lmp(i,1),dihd_lmp(i,2)&
             &,dihd_lmp(i,3),dihd_lmp(i,4),dihd_lmp(i,5),dihd_lmp(i,6)
        
     END DO

  END IF

  CLOSE(outfile)

END SUBROUTINE EXTRACT_AND_WRITE

!--------------------------------------------------------------------

SUBROUTINE OUTPUT_ALL()

  USE EXTRACT_PARAMS
  IMPLICIT NONE

  INTEGER :: ierr

  OPEN(unit = 25,file="US_outstat.dat",action="write",status&
       &="replace",iostat = ierr)

  IF(ierr /= 0) STOP "Could not open US_outstat.dat"
  
  IF(pullaxis == 0) WRITE(25,*) "Random Pull axis"
  IF(pullaxis == 1) WRITE(25,*) "Pull axis: X-axis"
  IF(pullaxis == 2) WRITE(25,*) "Pull axis: Y-axis"
  IF(pullaxis == 3) WRITE(25,*) "Pull axis: Z-axis"

  IF(multiposflag) THEN
     WRITE(25,*) "Saving interval: ", savedist
  ELSE
     WRITE(25,*) "Saved configuration at ", savedist
  END IF

  WRITE(25,*) "Tolerance: ", tol

  IF(multiposflag) THEN
     WRITE(25,*) "Output configuration style: Multiple Configurations"
  ELSE
     WRITE(25,*) "Output configuration style: Single Configuration"
  END IF
  
  IF(multiposflag == 0) THEN
     IF(netgenflag == 0) THEN
        WRITE(25,*) "Could not find a suitable configuration"
     ELSE
        WRITE(25,*) "Written initial configuration"
     END IF
  END IF

  CLOSE(25)
  
END SUBROUTINE OUTPUT_ALL

!--------------------------------------------------------------------

SUBROUTINE ALLOCATE_ARRAYS()

  USE EXTRACT_PARAMS
  IMPLICIT NONE

  INTEGER :: AllocateStatus

! Allocate LAMMPS Structure

  ALLOCATE(aidvals(ntotatoms,3),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate aidvals"
  ALLOCATE(rxyz_lmp(ntotatoms,3),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate rxyz_lmp"
  ALLOCATE(charge_lmp(ntotatoms,1),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate charge_lmp"
  DEALLOCATE(charge_lmp)
  ALLOCATE(vel_xyz(ntotatoms,4),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate vel_xyz"

  IF(ntotbonds /= 0) THEN
     ALLOCATE(bond_lmp(ntotbonds,4),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate bond_lmp"
  ELSE
     PRINT *, "Warning: No bonds - Not correct for bonded systems"
     ALLOCATE(bond_lmp(1,1),stat = AllocateStatus)
     DEALLOCATE(bond_lmp)
  END IF
  
  IF(ntotangls /= 0) THEN
     ALLOCATE(angl_lmp(ntotangls,5),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate angl_lmp"
  ELSE
     ALLOCATE(angl_lmp(1,1),stat = AllocateStatus)
     DEALLOCATE(angl_lmp)
  END IF
     
  IF(ntotdihds /= 0) THEN
     ALLOCATE(dihd_lmp(ntotdihds,6),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate dihd_lmp"
  ELSE
     ALLOCATE(dihd_lmp(1,1),stat = AllocateStatus)
     DEALLOCATE(dihd_lmp)
  END IF
  
  IF(ntotimprs /= 0) THEN
     ALLOCATE(impr_lmp(ntotimprs,6),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate zlmp"
  ELSE
     ALLOCATE(impr_lmp(1,1),stat = AllocateStatus)
     DEALLOCATE(impr_lmp)
  END IF


! Allocate Box details

  ALLOCATE(boxx_arr(nframes),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate boxx_arr"
  ALLOCATE(boxy_arr(nframes),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate boxy_arr"
  ALLOCATE(boxz_arr(nframes),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate boxz_arr"

  PRINT *, "Successfully allocated memory"

END SUBROUTINE ALLOCATE_ARRAYS

!--------------------------------------------------------------------

SUBROUTINE DEALLOCATE_ARRAYS()

  USE EXTRACT_PARAMS

  IMPLICIT NONE

  DEALLOCATE(aidvals)
  DEALLOCATE(rxyz_lmp)
  DEALLOCATE(vel_xyz)

  IF(ntotbonds /= 0) DEALLOCATE(bond_lmp)
  IF(ntotangls /= 0) DEALLOCATE(angl_lmp)
  IF(ntotdihds /= 0) DEALLOCATE(dihd_lmp)
  IF(ntotimprs /= 0) DEALLOCATE(impr_lmp)

END SUBROUTINE DEALLOCATE_ARRAYS

!--------------------------------------------------------------------
