!********************************************************************
!-------------------------Main Program File--------------------------
!--------------------------Ver:Jun-28-2019---------------------------
!-------For analyzing static properties of Solvated MethylCellulose--
!--------------------Parameter file - allparams.f90------------------
!-----------Definition of subroutines in the parameter file----------
!********************************************************************

PROGRAM MAIN

  USE PARAMS_SOLVATEDMC

  IMPLICIT NONE

  PRINT *, "Static analysis of PS-PEO system with ions .."
  PRINT *, "Starting OMP Threads .."

  PRINT *, "MAKE SURE THE INPUT DATAFILE HAS THE FOLLOWING ORDER FOR  &
       & ATOM IDs : BACKBONES (1-2-3..N) - GRAFT OF (1-2-3..N) - SOLVE&
       &NT"

!$OMP PARALLEL
  nproc = OMP_GET_NUM_THREADS()
!$OMP END PARALLEL

  PRINT *, "Number of threads are: ", nproc

  CALL READ_ANA_IP_FILE()
  CALL READ_DATAFILE()
  CALL ANALYZE_TRAJECTORYFILE()
  CALL ALLOUTPUTS()
  CALL DEALLOCATE_ARRAYS()

  PRINT *, "All Calculations Completed Succesfully :)"

END PROGRAM MAIN

!--------------------------------------------------------------------

SUBROUTINE READ_ANA_IP_FILE()

  USE PARAMS_SOLVATEDMC

  IMPLICIT NONE
  
  INTEGER :: nargs,ierr,logflag,AllocateStatus,i
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

     ELSEIF(dumchar == 'skipfr') THEN

        READ(anaread,*,iostat=ierr) skipfr
        
     ELSEIF(dumchar == 'nchains') THEN

        READ(anaread,*,iostat=ierr) nchains

     ELSEIF(dumchar == 'solv_type') THEN

        READ(anaread,*,iostat=ierr) num_at_type_solv
        ALLOCATE(solv_attypes(num_at_type_solv),stat =&
             & AllocateStatus)
        IF(AllocateStatus/=0) STOP "did not allocate solv_attypes"

        READ(anaread,*,iostat=ierr) (solv_attypes(i),i=1&
             &,num_at_type_solv)

     ELSEIF(dumchar == 'chain_at_types') THEN

        !This includes backbone + graft
        READ(anaread,*,iostat=ierr) num_at_type_perchain
        ALLOCATE(chain_attypes(num_at_type_perchain),stat =&
             & AllocateStatus)
        IF(AllocateStatus/=0) STOP "did not allocate chain_attypes"

        READ(anaread,*,iostat=ierr) (chain_attypes(i),i=1&
             &,num_at_type_perchain)
        
     ELSEIF(dumchar == 'graft_type') THEN

        READ(anaread,*,iostat=ierr) graft_type
        
     ELSEIF(dumchar == 'nmonsbackbone') THEN
        !Only backbone number of monomers
        READ(anaread,*,iostat=ierr) nmons_backbone

     ELSEIF(dumchar == 'compute_rdf') THEN

        rdfcalc = 1
        READ(anaread,*,iostat=ierr) rdffreq, rmaxbin, rdomcut,npairs
        
        ALLOCATE(pairs_rdf(npairs,3),stat = AllocateStatus)
        IF(AllocateStatus/=0) STOP "did not allocate pairs_rdf"
      
        DO i = 1,npairs

           READ(anaread,*,iostat=ierr) pairs_rdf(i,1), pairs_rdf(i,2)

        END DO

     ELSEIF(dumchar == 'compute_rg') THEN

        rgcalc = 1
        READ(anaread,*,iostat=ierr) rgfreq,rgall

     ELSEIF(dumchar == 'compute_com') THEN

        comcalc = 1
        READ(anaread,*,iostat=ierr) comfreq


     ELSEIF(dumchar == 'compute_sysrg') THEN

        rgsys = 1
        READ(anaread,*,iostat=ierr) rg_s_freq

     ELSEIF(dumchar == 'compute_oxydist') THEN

        oxycalc = 1
        READ(anaread,*,iostat=ierr) oxytype, oxyfreq, oxycut

     ELSEIF(dumchar == 'compute_resauto') THEN

        rescalc = 1
        READ(anaread,*,iostat=ierr) oxyrestype, oxyrescut

     ELSEIF(dumchar == 'log_file') THEN

        READ(anaread,*,iostat=ierr) log_fname
        logflag  = 1

     ELSEIF(dumchar == 'compute_perslen') THEN

        !skip_tailmon just tells to cutoff that part of the chain to
        !avoid end effects. If it is 0, then no part of the chain is
        !cutoff

        READ(anaread,*,iostat=ierr) perscalc, skip_tailmon

        IF(nmons_backbone == 0) STOP "No: of monomers in MC unspecifie&
             &d"

        IF(skip_tailmon .GE. INT(0.5*nmons_backbone) .OR.&
             & skip_tailmon .LT. 0) STOP "Unphysical end monomer cutof&
             &f"
        nanglmain = nmons_backbone  - 2*(skip_tailmon+1)
        PRINT *, "Number of angles utilized for computing persistence &
             &length and the number of tail monomers (on each side) sk&
             &ipped",nanglmain,skip_tailmon

        WRITE(logout,*) "Number of angles utilized for computing persi&
             &stence length and the number of tail monomers skipped on&
             & each side",nanglmain,skip_tailmon

     ELSEIF(dumchar == 'segment_re') THEN

        READ(anaread,*,iostat=ierr) resegfreq
        resegcalc = 1

     ELSEIF(dumchar == 'compute_eigvals') THEN

        READ(anaread,*,iostat=ierr) eigcalc, indeig, eigMC
        IF(eigcalc .OR. indeig .OR. eigMC) eigcalc  = 1

     ELSEIF(dumchar == 'compute_density') THEN

        denscalc = 1
        READ(anaread,*,iostat=ierr) densfreq,ndentypes,dens_axis&
             &,maxden_bin

        ALLOCATE(dentyp_arr(ndentypes),stat = AllocateStatus)
        IF(AllocateStatus/=0) STOP "did not allocate dentyp_arr"
      
        READ(anaread,*,iostat=ierr) (dentyp_arr(i),i=1,ndentypes)

     ELSE
        
        PRINT *, "unknown keyword", trim(dumchar)
        STOP

     END IF

  END DO

  IF(eigMC == 1 .AND. nmons_backbone == 0) STOP "Undefined number of mons"
  
  IF(logflag == 0) log_fname = "log."//trim(adjustl(traj_fname))
  OPEN(unit = logout,file=trim(log_fname),action="write",status="repla&
       &ce",iostat=ierr)

  PRINT *, "Analysis input file read finished .."

END SUBROUTINE READ_ANA_IP_FILE

!--------------------------------------------------------------------

SUBROUTINE DEFAULTVALUES()

  USE PARAMS_SOLVATEDMC
  IMPLICIT NONE

  ! Frame, Molecules and Processor Details
  nframes = 0; skipfr = 0; nsolvent = 0; nchains = 0
  nproc = 0; nmons_backbone = 0
  num_at_type_perchain = 0; num_at_type_solv = 0


  ! Initialize Flags
  rgall = 0; rgcalc = 0;rdfcalc = 0; oxycalc = 0; rescalc = 0
  eigcalc = 0; rgsys = 0; denscalc = 0; indeig = 0; eigMC = 0
  perscalc = 0; comcalc = 0; resegcalc = 0
  
  ! Initialize distributions and frequencies
  rdffreq = 1; rgfreq = 1; oxyfreq = 1;  rg_s_freq=1
  densfreq = 1; dens_axis = 0; ndentypes = 0;comfreq = 1
  resegfreq = 1;

  ! Initialzie Extra Structural Quantities
  rdomcut = 0; rclus_cut = 0; rmaxbin = 0; rbinval = 0
  oxycut = 0; oxytype = 0
  oxyrescut = 0; oxyrestype = 0

  !Averages
  rvolavg = 0

  !Structural quantities
  skip_tailmon = 0

END SUBROUTINE DEFAULTVALUES

!--------------------------------------------------------------------

SUBROUTINE READ_DATAFILE()

  USE PARAMS_SOLVATEDMC

  IMPLICIT NONE

  INTEGER :: i,j,k,ierr,u,AllocateStatus
  INTEGER :: flag, cntr, nwords, imax
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

  DO i = 1,imax-2 !Change here according to convenience
       
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

           READ(inpread,*) aid,molid,atype,charge,rx,ry,rz

           rx = rx - xlo
           ry = ry - ylo
           rz = rz - zlo

           aidvals(aid,1)     = aid
           aidvals(aid,2)     = molid
           aidvals(aid,3)     = atype
           charge_lmp(aid,1)  = charge
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

  USE PARAMS_SOLVATEDMC

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

  USE PARAMS_SOLVATEDMC

  IMPLICIT NONE

  INTEGER :: i,j,aid,ierr,atchk,atype
  REAL :: xlo,xhi,ylo,yhi,zlo,zhi

  OPEN(unit = 15,file =trim(traj_fname),action="read",status="old"&
       &,iostat=ierr)

  IF(ierr /= 0) STOP "trajectory file not found"

  WRITE(logout,*) "trajectory file used is :",&
       & trim(adjustl(traj_fname))


  CALL STRUCT_INIT()
  CALL OPEN_STRUCT_OUTPUT_FILES()

  DO i = 1,skipfr

     DO j = 1,ntotatoms+9

        READ(15,*) 

     END DO

     IF(mod(i,100) == 0) PRINT *, "Skipped ", i, "frames"

  END DO

  DO i = 1,nframes

     IF(mod(i,50) == 0) PRINT *, "Processing ", i,"th frame"

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

     CALL STRUCT_MAIN(i)

  END DO

  CLOSE(15)

END SUBROUTINE ANALYZE_TRAJECTORYFILE

!--------------------------------------------------------------------

SUBROUTINE OPEN_STRUCT_OUTPUT_FILES()

  USE PARAMS_SOLVATEDMC

  IMPLICIT NONE

  IF(rgcalc) THEN
     
     dum_fname = "rgavgall_"//trim(adjustl(traj_fname))
     OPEN(unit = rgavgwrite,file =trim(dum_fname),action="write"&
          &,status="replace")
     
     dum_fname = "rgMConly_"//trim(adjustl(traj_fname))
     OPEN(unit = rergwrite,file =trim(dum_fname),action="write"&
          &,status="replace")

     IF(rgall) THEN
        dum_fname = "rgall_"//trim(adjustl(traj_fname))
        OPEN(unit = rgwrite,file =trim(dum_fname),action="write"&
             &,status="replace")
     END IF

  END IF

  IF(rgsys) THEN
     dum_fname = "rgsys_"//trim(adjustl(traj_fname))
     OPEN(unit = rgswrite,file =trim(dum_fname),action="write",status&
          &="replace")
  END IF

  IF(indeig) THEN
     dum_fname = "indeig_"//trim(adjustl(traj_fname))
     OPEN(unit = eigwrite,file =trim(dum_fname),action="write",status&
          &="replace")
  END IF

  IF(eigMC) THEN
     dum_fname = "eigMCavg_"//trim(adjustl(traj_fname))
     OPEN(unit = eigMCwrite,file =trim(dum_fname),action="write"&
          &,status="replace")
  END IF

  IF(comcalc) THEN
     dum_fname = "compos_"//trim(adjustl(traj_fname))
     OPEN(unit = comwrite,file =trim(dum_fname),action="write"&
          &,status="replace")

     IF(nchains == 2) THEN
        dum_fname = "distcom_"//trim(adjustl(traj_fname))
        OPEN(unit = distcomwrite,file =trim(dum_fname),action="write"&
             &,status="replace")
     END IF

  END IF

END SUBROUTINE OPEN_STRUCT_OUTPUT_FILES

!--------------------------------------------------------------------

SUBROUTINE STRUCT_INIT()

  USE PARAMS_SOLVATEDMC
  IMPLICIT NONE

  INTEGER :: i,j,t1,t2
  INTEGER :: AllocateStatus

  CALL COUNT_ATOMS_PER_CHAIN()
  CALL COUNT_SOLVENT()
  CALL LOGWRITE()

  IF(rdfcalc) THEN

     rdfarray = 0.0
     rbinval = rdomcut/REAL(rmaxbin)

     DO i = 1, npairs

        t1 = 0; t2 = 0

        DO j = 1,ntotatoms

           IF(aidvals(j,3) == pairs_rdf(i,1)) t1 = t1+1
           IF(aidvals(j,3) == pairs_rdf(i,2)) t2 = t2+1
           
        END DO

        pairs_rdf(i,3) = t1*t2

     END DO

  END IF


  IF(oxycalc) THEN

     histoxyarray = 0.0

     PRINT *, "Number of oxygen types: ", nsolvent

     ALLOCATE(oxyarray(nsolvent),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate oxyarray"
     j = 0

     DO i = 1,ntotatoms

        j = j + 1
        IF(aidvals(i,3) == oxytype) oxyarray(j) = aidvals(i,1)
   
     END DO

     IF(j .NE. nsolvent) THEN
        
        PRINT *, "Oxyarray not counted correctly", j, nsolvent
        STOP

     END IF

  END IF

  IF(rescalc) THEN

     ALLOCATE(autocf(nsolvent,nframes),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate autocf"
     ALLOCATE(oxyresarray(nsolvent),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate oxyresarray"

     autocf = 0.0; oxyresarray = 0.0

     j = 0

     DO i = 1,ntotatoms
        
        j = j + 1
        IF(aidvals(i,3) == oxyrestype) oxyresarray(j) = aidvals(i,1)
        
     END DO
        
     IF(j .NE. nsolvent) THEN
        
        PRINT *, "Oxyresarray not counted correctly", j, nsolvent
        STOP
        
     END IF

  END IF

  IF(indeig .OR. eigcalc) THEN

     ALLOCATE(eigarray(nframes,3),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate eigarray"
     ALLOCATE(timestep_arr(nframes),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate timestep_arr"     
     eigarray = 0.0
     timestep_arr = 0

  END IF

  IF(perscalc) THEN
     
     ALLOCATE(avgtheta_main(0:nanglmain),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate avgtheta_main"
     ALLOCATE(avgdot_main(nmons_backbone-1),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate avgdot_main"

     avgtheta_main = 0.0;avgdot_main = 0.0

  END IF

  IF(denscalc) THEN

     densarray = 0.0; normdens = 0.0; denbinavg = 0.0

  END IF


  IF(resegcalc) THEN
     ALLOCATE(resegarr(nmons_backbone), stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate resegarr"
     resegarr = 0.0
  END IF


END SUBROUTINE STRUCT_INIT

!--------------------------------------------------------------------

SUBROUTINE STRUCT_MAIN(tval)

  USE PARAMS_SOLVATEDMC
  IMPLICIT NONE

  INTEGER, INTENT(IN):: tval

  IF(rgcalc .AND. mod(tval-1,rgfreq)==0) CALL COMPUTE_RADGYR(tval)
  IF(rgcalc .AND. mod(tval-1,rgfreq)==0) CALL COMPUTE_RADMC(tval)
  IF(rdfcalc .AND. mod(tval-1,rdffreq)==0) CALL COMPUTE_RDF(tval)
  IF(oxycalc .AND. mod(tval-1,oxyfreq)==0) CALL COMPUTE_OXYCLUS(tval)
  IF(rescalc) CALL COMPUTE_RESTIME(tval)
  IF(rescalc .AND. tval == nframes) CALL COMPUTE_RESTIMESPECTRA()
  IF(perscalc) CALL COMPUTE_PERSISTENCELEN(tval)
  IF(perscalc) CALL COMPUTE_PERSLENMETHOD2(tval)
  IF(indeig .OR. eigcalc) CALL COMPUTE_EIGVALS(tval)
  IF(eigMC) CALL COMPUTE_EIGVALSOFMCCHAIN(tval)
  IF(comcalc) CALL COMPUTE_COM(tval)
  IF(rgsys .AND. mod(tval-1,rg_s_freq)==0) CALL COMPUTE_RGSYS(tval)
  IF(denscalc .AND. mod(tval-1,densfreq)==0) CALL COMPUTE_DENS(tval)
  IF(resegcalc) CALL COMPUTE_RESEGDIST(tval)

END SUBROUTINE STRUCT_MAIN

!--------------------------------------------------------------------

SUBROUTINE COUNT_SOLVENT()

  USE PARAMS_SOLVATEDMC
  IMPLICIT NONE

  INTEGER :: i,j,atype

  nsolvent = 0

  DO i = 1,ntotatoms

     DO j = 1,num_at_type_solv

        IF(aidvals(i,3) == solv_attypes(j)) THEN

           nsolvent = nsolvent + 1

        END IF

     END DO
     
  END DO

END SUBROUTINE COUNT_SOLVENT

!--------------------------------------------------------------------

SUBROUTINE COUNT_ATOMS_PER_CHAIN()

  USE PARAMS_SOLVATEDMC
  IMPLICIT NONE

  INTEGER :: i,j,molid,atype
  
  DO i = 1,ntotatoms

     DO j = 1,num_at_type_perchain

        IF(aidvals(i,3) == chain_attypes(j)) THEN

           molid = aidvals(i,2)
           num_atoms_perchain(molid) = num_atoms_perchain(molid)+1

        END IF

     END DO
     
  END DO

END SUBROUTINE COUNT_ATOMS_PER_CHAIN

!--------------------------------------------------------------------

SUBROUTINE LOGWRITE()

  USE PARAMS_SOLVATEDMC
  IMPLICIT NONE

  WRITE(logout,*) "------------------------------------------------"
  WRITE(logout,*) "num_atoms_perchain details: NCHAINS = ", nchains
  WRITE(logout,*) num_atoms_perchain(:)
  WRITE(logout,*) "------------------------------------------------"

  WRITE(logout,*) "Number of solvents: ", nsolvent
  WRITE(logout,*) "Number of backbone atoms/chain: ", nmons_backbone

END SUBROUTINE LOGWRITE

!--------------------------------------------------------------------

SUBROUTINE COMPUTE_RADGYR(iframe)
  !Inclusive of graft
  USE PARAMS_SOLVATEDMC
  IMPLICIT NONE

  INTEGER :: i,j,molid,atype,tot_chain_atoms
  REAL, DIMENSION(1:nchains) :: rgxx, rgyy, rgzz, rgsq
  REAL, DIMENSION(1:nchains) :: rxcm, rycm, rzcm, totmass
  REAL :: rgsqavg, rgxxavg, rgyyavg, rgzzavg, resqavg
  INTEGER, INTENT(IN) :: iframe

  rgxx = 0.0; rgyy =0.0; rgzz = 0.0; rgsq = 0.0
  totmass = 0.0
  rgsqavg = 0.0; rgxxavg = 0.0; rgyyavg = 0.0; rgzzavg = 0.0
  rxcm = 0.0; rycm = 0.0; rzcm = 0.0
  
  tot_chain_atoms = SUM(num_atoms_perchain) !Relevant when solvent is
  !present, or else this is equal to the total number of particles in
  !the system. 

  IF(iframe == 1) THEN
     WRITE(logout,*) "Total chain atoms: ",tot_chain_atoms
     WRITE(logout,*) "Masses per type: "
     WRITE(logout,*) masses
  END IF

  DO i = 1,tot_chain_atoms

     molid = aidvals(i,2)
     atype = aidvals(i,3)
     totmass(molid) = totmass(molid) + masses(atype,1)

     rxcm(molid) = rxcm(molid)+ rxyz_lmp(i,1)*masses(atype,1)
     rycm(molid) = rycm(molid)+ rxyz_lmp(i,2)*masses(atype,1)
     rzcm(molid) = rzcm(molid)+ rxyz_lmp(i,3)*masses(atype,1)

  END DO

  DO i = 1,nchains

     rxcm(i) = rxcm(i)/totmass(i)
     rycm(i) = rycm(i)/totmass(i)
     rzcm(i) = rzcm(i)/totmass(i)

  END DO

  
  IF(iframe == 1) THEN

     OPEN(unit = 98,file ="totmasschk.txt",action="write",status="repl&
          &ace")

     DO i = 1,nchains

        WRITE(98,'(I0,1X,4(F14.9,1X))') i, totmass(i),rxcm(i),&
             & rycm(i), rzcm(i)

     END DO

     CLOSE(98)

     OPEN(unit = 98,file ="molidchk.txt",action="write",status="repl&
          &ace")

     DO i = 1,ntotatoms

        WRITE(98,'(I0,1X,I0)') i, aidvals(i,2)

     END DO

     CLOSE(98)

  END IF
  
  DO i = 1,tot_chain_atoms

     molid = aidvals(i,2)
     atype = aidvals(i,3)

     rgxx(molid) = rgxx(molid) + masses(atype,1)*((rxyz_lmp(i,1)&
          &-rxcm(molid))**2)
     rgyy(molid) = rgyy(molid) + masses(atype,1)*((rxyz_lmp(i,2)&
          &-rycm(molid))**2)
     rgzz(molid) = rgzz(molid) + masses(atype,1)*((rxyz_lmp(i,3)&
          &-rzcm(molid))**2)

     rgsq(molid) = rgsq(molid) + masses(atype,1)*((rxyz_lmp(i,1)&
          &-rxcm(molid))**2 + (rxyz_lmp(i,2)-rycm(molid))**2 +&
          & (rxyz_lmp(i,3)-rzcm(molid))**2)

  END DO

  DO i = 1,nchains

     rgsq(i) = rgsq(i)/totmass(i)
     rgxx(i) = rgxx(i)/totmass(i)
     rgyy(i) = rgyy(i)/totmass(i)
     rgzz(i) = rgzz(i)/totmass(i)

  END DO


  DO i = 1,nchains

     rgsqavg = rgsqavg + rgsq(i)
     rgxxavg = rgxxavg + rgxx(i)
     rgyyavg = rgyyavg + rgyy(i)
     rgzzavg = rgzzavg + rgzz(i)
     
  END DO
  
  rgsqavg = rgsqavg/REAL(nchains)
  rgxxavg = rgxxavg/REAL(nchains)
  rgyyavg = rgyyavg/REAL(nchains)
  rgzzavg = rgzzavg/REAL(nchains)
  

  IF(iframe == 1) WRITE(rgavgwrite,'(A)') 'Timestep  sqrt(rgxxavg)&
       & sqrt(rgyyavg) sqrt(rgzzavg) sqrt(rgsqavg)'

  WRITE(rgavgwrite,'(I0,1X,4(F14.6,1X))') timestep, sqrt(rgxxavg),&
       & sqrt(rgyyavg), sqrt(rgzzavg), sqrt(rgsqavg)

  IF(rgall) THEN
     
     IF(iframe == 1) WRITE(rgwrite,'(A)') 'chain #, rgxx(i), rgyy(i), &
          &rgzz(i), sqrt(rgsq(i))'

     WRITE(rgwrite,'(2(I0,1X))') timestep, nchains

     DO i = 1,nchains
     
        WRITE(rgwrite,'(I0,1X,4(F14.6,1X))') i,rgxx(i),rgyy(i),&
             & rgzz(i),sqrt(rgsq(i))

     END DO

  END IF
     
END SUBROUTINE COMPUTE_RADGYR

!--------------------------------------------------------------------

SUBROUTINE COMPUTE_RESEGDIST(tval)

  USE PARAMS_SOLVATEDMC
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: tval

END SUBROUTINE COMPUTE_RESEGDIST

!--------------------------------------------------------------------

SUBROUTINE COMPUTE_RADMC(iframe)
  !Only the backbone
  USE PARAMS_SOLVATEDMC
  IMPLICIT NONE

  INTEGER :: i,j,k,molid,atype,initmon, finmon
  REAL, DIMENSION(1:nchains) :: rgxx, rgyy, rgzz, rgsq
  REAL, DIMENSION(1:nchains) :: rxcm, rycm, rzcm, totmass
  REAL :: rgsqavg, rgxxavg, rgyyavg, rgzzavg, resqavg
  REAL :: bsqavg, rx,ry,rz
  INTEGER, INTENT(IN) :: iframe

  rgxx = 0.0; rgyy =0.0; rgzz = 0.0; rgsq = 0.0; resqavg = 0.0
  totmass = 0.0; bsqavg = 0.0
  rgsqavg = 0.0; rgxxavg = 0.0; rgyyavg = 0.0; rgzzavg = 0.0
  rxcm = 0.0; rycm = 0.0; rzcm = 0.0
  
  IF(nmons_backbone == 0) STOP "Cannot compute Rg of MC chain"

  IF(iframe == 1) PRINT *, "Atoms in MC: ", nmons_backbone

  initmon = 0

!Resq, COM
  DO i = 1,nchains

     initmon = (i-1)*nmons_backbone + 1
     finmon  = i*nmons_backbone
     
     IF(aidvals(initmon,2) .NE. aidvals(finmon,2)) THEN

        PRINT *, "Not computing the same chain for Re**2",initmon&
             &,finmon
        STOP

     END IF
        
     resqavg = resqavg + ((rxyz_lmp(initmon,1)-rxyz_lmp(finmon,1))**2&
          & + (rxyz_lmp(initmon,2)-rxyz_lmp(finmon,2))**2 +&
          & (rxyz_lmp(initmon,3)-rxyz_lmp(finmon,3))**2)
     
     DO j = 1,nmons_backbone

        k = (i-1)*nmons_backbone + j
        molid = aidvals(k,2)
        atype = aidvals(k,3)
        totmass(molid) = totmass(molid) + masses(atype,1)

        rxcm(molid) = rxcm(molid)+ rxyz_lmp(k,1)*masses(atype,1)
        rycm(molid) = rycm(molid)+ rxyz_lmp(k,2)*masses(atype,1)
        rzcm(molid) = rzcm(molid)+ rxyz_lmp(k,3)*masses(atype,1)

     END DO

  END DO

  resqavg = resqavg/nchains

  DO i = 1,nchains

     rxcm(i) = rxcm(i)/totmass(i)
     rycm(i) = rycm(i)/totmass(i)
     rzcm(i) = rzcm(i)/totmass(i)

  END DO


  bsqavg = 0.0
! Rgsq, bsq
  DO i = 1,nchains

     initmon = (i-1)*nmons_backbone

     DO j = 1, nmons_backbone !separate into two loops so that the
        !grafts can be separated out while computing rg

        k = initmon + j

        molid = aidvals(k,2)
        atype = aidvals(k,3)

        rgxx(molid) = rgxx(molid) + masses(atype,1)*((rxyz_lmp(k,1)&
             &-rxcm(molid))**2)
        rgyy(molid) = rgyy(molid) + masses(atype,1)*((rxyz_lmp(k,2)&
             &-rycm(molid))**2)
        rgzz(molid) = rgzz(molid) + masses(atype,1)*((rxyz_lmp(k,3)&
             &-rzcm(molid))**2)

        rgsq(molid) = rgsq(molid) + masses(atype,1)*((rxyz_lmp(k,1)&
             &-rxcm(molid))**2 + (rxyz_lmp(k,2)-rycm(molid))**2 +&
             & (rxyz_lmp(k,3)-rzcm(molid))**2)



        IF(j .LT. nmons_backbone) THEN

           rx = rxyz_lmp(k,1) - rxyz_lmp(k+1,1)
           ry = rxyz_lmp(k,2) - rxyz_lmp(k+1,2)
           rz = rxyz_lmp(k,3) - rxyz_lmp(k+1,3)
           bsqavg = bsqavg + rx**2 + ry**2 + rz**2

        END IF
        
     END DO
     
  END DO

  bsqavg = bsqavg/REAL(nmons_backbone-1)

  DO i = 1,nchains

     rgsq(i) = rgsq(i)/totmass(i)
     rgxx(i) = rgxx(i)/totmass(i)
     rgyy(i) = rgyy(i)/totmass(i)
     rgzz(i) = rgzz(i)/totmass(i)
     
  END DO


  DO i = 1,nchains

     rgsqavg = rgsqavg + rgsq(i)
     rgxxavg = rgxxavg + rgxx(i)
     rgyyavg = rgyyavg + rgyy(i)
     rgzzavg = rgzzavg + rgzz(i)
     
  END DO
  
  rgsqavg = rgsqavg/REAL(nchains)
  rgxxavg = rgxxavg/REAL(nchains)
  rgyyavg = rgyyavg/REAL(nchains)
  rgzzavg = rgzzavg/REAL(nchains)

  
  IF(iframe == 1) WRITE(rergwrite,'(A)') 'Timestep  sqrt(rgxxavg)&
       & sqrt(rgyyavg) sqrt(rgzzavg) sqrt(rgsqavg)&
       & sqrt(resqavg) sqrt(bsqavg)'
  
  WRITE(rergwrite,'(I0,1X,6(F14.6,1X))') timestep, sqrt(rgxxavg),&
       & sqrt(rgyyavg), sqrt(rgzzavg), sqrt(rgsqavg),&
       & sqrt(resqavg), sqrt(bsqavg)


END SUBROUTINE COMPUTE_RADMC

!--------------------------------------------------------------------

SUBROUTINE COMPUTE_RGSYS(iframe)

  USE PARAMS_SOLVATEDMC
  IMPLICIT NONE

  INTEGER :: i,j,molid,atype,tot_chain_atoms
  REAL :: rgxx, rgyy, rgzz, rgsq
  REAL :: rxcm, rycm, rzcm, totmass
  INTEGER, INTENT(IN) :: iframe

  rgxx = 0.0; rgyy =0.0; rgzz = 0.0; rgsq = 0.0
  totmass = 0.0
  rxcm = 0.0; rycm =0.0; rzcm = 0.0 
  
  tot_chain_atoms = SUM(num_atoms_perchain)

  DO i = 1,tot_chain_atoms

     molid = aidvals(i,2)
     atype = aidvals(i,3)
     totmass = totmass + masses(atype,1)

     rxcm = rxcm + rxyz_lmp(i,1)*masses(atype,1)
     rycm = rycm + rxyz_lmp(i,2)*masses(atype,1)
     rzcm = rzcm + rxyz_lmp(i,3)*masses(atype,1)

  END DO

  rxcm = rxcm/totmass
  rycm = rycm/totmass
  rzcm = rzcm/totmass

  DO i = 1,tot_chain_atoms

     molid = aidvals(i,2)
     atype = aidvals(i,3)
     rgxx = rgxx + masses(atype,1)*((rxyz_lmp(i,1)-rxcm)**2)
     rgyy = rgyy + masses(atype,1)*((rxyz_lmp(i,2)-rycm)**2)
     rgzz = rgzz + masses(atype,1)*((rxyz_lmp(i,3)-rzcm)**2)

     rgsq = rgsq + masses(atype,1)*((rxyz_lmp(i,1)-rxcm)**2 +&
          & (rxyz_lmp(i,2)-rycm)**2 + (rxyz_lmp(i,3)-rzcm)**2)

  END DO

  rgsq = rgsq/totmass
  rgxx = rgxx/totmass
  rgyy = rgyy/totmass
  rgzz = rgzz/totmass


  IF(iframe == 1) WRITE(rgswrite,'(A)') 'Timestep, rgxx, rgyy, rgzz, r&
       &gsq'
  WRITE(rgswrite,'(I0,1X,4(F14.6,1X))') timestep,rgxx,rgyy,rgzz,rgsq
     
END SUBROUTINE COMPUTE_RGSYS

!--------------------------------------------------------------------

SUBROUTINE COMPUTE_RDF(iframe)

  USE PARAMS_SOLVATEDMC
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: iframe
  INTEGER :: i,j,a1type,a2type,ibin,a1id,a2id,paircnt
  REAL :: rxval,ryval,rzval,rval
  INTEGER :: a1ref,a2ref, AllocateStatus
  INTEGER,ALLOCATABLE, DIMENSION(:,:) :: dumrdfarray

  rvolval = box_xl*box_yl*box_zl
  rvolavg = rvolavg + rvolval 

  ALLOCATE(dumrdfarray(0:rmaxbin-1,npairs),stat=AllocateStatus)
  IF(AllocateStatus/=0) STOP "dumrdfarray not allocated"
  dumrdfarray = 0

  DO paircnt = 1,npairs

     a1ref = pairs_rdf(paircnt,1); a2ref = pairs_rdf(paircnt,2)
     
     DO i = 1,ntotatoms
        
        a1id   = aidvals(i,1)     
        a1type = aidvals(i,3)
        
        DO j = 1,ntotatoms

           a2id   = aidvals(j,1)        
           a2type = aidvals(j,3)

           IF((a1type == a1ref .AND. a2type == a2ref) .AND. (a1id&
                & .NE. a2id)) THEN        

              rxval = rxyz_lmp(a1id,1) - rxyz_lmp(a2id,1) 
              ryval = rxyz_lmp(a1id,2) - rxyz_lmp(a2id,2) 
              rzval = rxyz_lmp(a1id,3) - rxyz_lmp(a2id,3) 
              
              rxval = rxval - box_xl*ANINT(rxval/box_xl)
              ryval = ryval - box_yl*ANINT(ryval/box_yl)
              rzval = rzval - box_zl*ANINT(rzval/box_zl)
              
              rval = sqrt(rxval**2 + ryval**2 + rzval**2)
              ibin = FLOOR(rval/rbinval)
              
              IF(ibin .LT. rmaxbin) THEN
              
                 dumrdfarray(ibin,paircnt) = dumrdfarray(ibin&
                      &,paircnt) + 1
              
              END IF

           END IF

        END DO

     END DO

  END DO

!$OMP DO PRIVATE(i,j)
  DO j = 1,npairs
     
     DO i = 0,rmaxbin-1

        rdfarray(i,j) = rdfarray(i,j) + REAL(dumrdfarray(i,j))&
             &*rvolval/(REAL(pairs_rdf(j,3)))
        
     END DO
     
  END DO
!$OMP END DO

  DEALLOCATE(dumrdfarray)

END SUBROUTINE COMPUTE_RDF

!--------------------------------------------------------------------

SUBROUTINE COMPUTE_OXYCLUS(iframe)

  USE PARAMS_SOLVATEDMC
  IMPLICIT NONE
  
  INTEGER, INTENT(IN) :: iframe
  INTEGER :: AllocateStatus
  INTEGER :: i,j,a1id,molid,atcnt,cnt,a2id,chcnt,initmon
  REAL :: rxval, ryval, rzval, rval
  INTEGER, ALLOCATABLE, DIMENSION(:) :: inst_avg_num_chains

  ALLOCATE(inst_avg_num_chains(nsolvent),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate inst_avg_num_chains"

  inst_avg_num_chains = 0; atcnt = 0; cnt = 0
  box_xl = boxx_arr(iframe)
  box_yl = boxy_arr(iframe)
  box_zl = boxz_arr(iframe)

  DO i = 1,nsolvent
     
     a1id = oxyarray(i)
     initmon = 0

     DO chcnt = 1,nchains
        
        atcnt = 1
        
        DO WHILE(atcnt .LE. num_atoms_perchain(chcnt))
           
           j = initmon + atcnt
           
           rxval = rxyz_lmp(a1id,1) - rxyz_lmp(j,1)
           ryval = rxyz_lmp(a1id,2) - rxyz_lmp(j,2)
           rzval = rxyz_lmp(a1id,3) - rxyz_lmp(j,3)
           
           rxval = rxval - box_xl*ANINT(rxval/box_xl)
           ryval = ryval - box_yl*ANINT(ryval/box_yl)
           rzval = rzval - box_zl*ANINT(rzval/box_zl)
           
           rval = sqrt(rxval**2 + ryval**2 + rzval**2)
           
           IF(rval .LE. oxycut) THEN
              
              a2id  = aidvals(j,1)
              molid = aidvals(j,2)
              
              IF(molid .NE. chcnt) THEN
                 PRINT *, "Something wrong in assigning chain num"
                 PRINT *, j,chcnt,a2id,molid
                 STOP
              END IF
              
              inst_avg_num_chains(i) = inst_avg_num_chains(i) + 1
              atcnt = num_atoms_perchain(chcnt) + 1
              
           ELSE
              
              atcnt = atcnt + 1
              
           END IF
           
        END DO

        initmon = initmon + num_atoms_perchain(chcnt) - 1
        
     END DO
     
  END DO
  
!Add it to histoxyarray  

!$OMP PARALLEL DO PRIVATE(i,cnt)
  DO i = 1,nsolvent

     cnt = inst_avg_num_chains(i)
     histoxyarray(cnt) = histoxyarray(cnt)+ 1

  END DO
!$OMP END PARALLEL DO


END SUBROUTINE COMPUTE_OXYCLUS

!--------------------------------------------------------------------

SUBROUTINE COMPUTE_RESTIME(iframe)

  USE PARAMS_SOLVATEDMC
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: iframe
  INTEGER :: i,j,tval,a1id,a2id,ierr,ifin,chcnt,initmon
  REAL :: rxval,ryval,rzval,rval

  j = 1

  DO i = 1,nsolvent !populate autocorrelation fn array

     a1id = oxyresarray(i)
     initmon = 0

     DO WHILE(chcnt .LE. nchains)
        
        DO WHILE(j .LE. num_atoms_perchain(chcnt))
           
           a2id = initmon + j
           
           rxval = rxyz_lmp(a1id,1) - rxyz_lmp(a2id,1) 
           ryval = rxyz_lmp(a1id,2) - rxyz_lmp(a2id,2) 
           rzval = rxyz_lmp(a1id,3) - rxyz_lmp(a2id,3) 
           
           rxval = rxval - box_xl*ANINT(rxval/box_xl)
           ryval = ryval - box_yl*ANINT(ryval/box_yl)
           rzval = rzval - box_zl*ANINT(rzval/box_zl)
           
           rval = sqrt(rxval**2 + ryval**2 + rzval**2)
           
           IF(rval .LT. oxyrescut) THEN
              
              autocf(i,iframe) = 1
              j = num_atoms_perchain(chcnt)+1
              chcnt = chcnt + 1
              
           ELSE
              
              j = j +1
              
           END IF
           
        END DO

        initmon = initmon + num_atoms_perchain(chcnt) - 1
        chcnt = chcnt + 1
     
     END DO

  END DO
  

END SUBROUTINE COMPUTE_RESTIME

!--------------------------------------------------------------------

SUBROUTINE COMPUTE_RESTIMESPECTRA()

  USE PARAMS_SOLVATEDMC
  IMPLICIT NONE

  INTEGER :: tinc, ifin, tim, i, j
  
!$OMP PARALLEL DO PRIVATE(tinc,ifin,tim,i,j)

  DO tinc = 0, nframes-1 !compute spectral product

     ifin = nframes - tinc
     
     DO i = 1,ifin
        
        tim = i + tinc
      
        DO j = 1,nsolvent

           tplot_cf(tinc) = tplot_cf(tinc) + REAL(autocf(j,tim)&
                &*autocf(j,i))
           
        END DO

     END DO

     tplot_cf(tinc) = tplot_cf(tinc)/REAL(ifin*nsolvent)
     
  END DO
!$OMP END PARALLEL DO

END SUBROUTINE COMPUTE_RESTIMESPECTRA

!--------------------------------------------------------------------

SUBROUTINE COMPUTE_EIGVALS(iframe)

  !For computing eigenvalue of chains INCLUDING grafts
  USE PARAMS_SOLVATEDMC
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: iframe
  INTEGER :: i,j,a1id,u,v,AllocateStatus,initmon,natchain
  INTEGER, PARAMETER :: N_LA = 3
  INTEGER, PARAMETER :: NMAX_LA = 3
  INTEGER, PARAMETER :: LWMAX = 1000
  INTEGER :: LWORK
  INTEGER :: INFO,LWKOPT
  REAL :: SMN(N_LA,3),W(N_LA),WORK(LWMAX)
  REAL :: SMN_ALL(N_LA,3)
  REAL :: rxcm, rycm, rzcm
  REAL,ALLOCATABLE,DIMENSION(:,:) :: rshift_lmp 

  !Note 1: Ordered in such a way that first chains, then the graft in
  !each chain. THIS IS VERY IMPORTANT.

  !IF(indeig)  WRITE(eigwrite,'(2(I0,1X))') timestep, nchains

  ! Create SMN (rij matrix)
  SMN_ALL = 0.0

  initmon = 0
  DO i = 1,nchains

     natchain = num_atoms_perchain(i)
     ALLOCATE(rshift_lmp(natchain,3),stat=AllocateStatus)
     IF(AllocateStatus/= 0) STOP "Did not allocate rshift_lmp"

     rshift_lmp = 0.0; SMN = 0.0
     rxcm = 0.0; rycm = 0.0; rzcm = 0.0
     CALL COMPUTE_CENTROID(i,rxcm,rycm,rzcm,num_atoms_perchain(i))
     
     DO j = 1, num_atoms_perchain(i)
        
        a1id = initmon + j
        rshift_lmp(j,1) = rxyz_lmp(a1id,1)-rxcm
        rshift_lmp(j,2) = rxyz_lmp(a1id,2)-rycm
        rshift_lmp(j,3) = rxyz_lmp(a1id,3)-rzcm
   
        DO u = 1,3

           DO v = 1,3

              SMN(u,v) = SMN(u,v)+rshift_lmp(j,u)*rshift_lmp(j,v)

           END DO
           
        END DO

     END DO


     DEALLOCATE(rshift_lmp)

     SMN = SMN/REAL(num_atoms_perchain(i))
     
     !Individual Eigenvalues
     IF(indeig /= 0) THEN
        
        ! Create Eigen Value Matrix for each chain
        ! Query Workspace and then compute lambda
        LWORK = -1
        CALL DSYEV('Vectors','Upper',N_LA,SMN,N_LA,W,WORK,LWORK,INFO)
        LWORK = MIN(LWMAX,INT(WORK(1)))
        CALL DSYEV('N','U',N_LA,SMN,N_LA,W,WORK,LWORK,INFO)
        IF(INFO .NE. 0) STOP 'Physically irrelevant Eigenvalues'
        WRITE(eigwrite,'(2(I0,1X),3(F14.6,1X))') i, timestep, W(1), W(2), W(3)
        
     END IF

     SMN_ALL = SMN + SMN_ALL

     initmon = initmon + num_atoms_perchain(i)

  END DO

  SMN_ALL = SMN_ALL/REAL(nchains)

  !Overall Eigenvalues
  LWORK = -1 !Query workspace and then compute lambda
  CALL DSYEV('Vectors','Upper',N_LA,SMN_ALL,N_LA,W,WORK,LWORK,INFO)
  LWORK = MIN(LWMAX,INT(WORK(1))) 
  CALL DSYEV('N','U',N_LA,SMN_ALL,N_LA,W,WORK,LWORK,INFO)
  IF(INFO .NE. 0) STOP 'Physically irrelevant Eigenvalues'      

  DO i = 1,3

     eigarray(iframe,i) = W(i)
     
  END DO
  
  timestep_arr(iframe) = timestep

END SUBROUTINE COMPUTE_EIGVALS

!--------------------------------------------------------------------

SUBROUTINE COMPUTE_EIGVALSOFMCCHAIN(iframe)

  USE PARAMS_SOLVATEDMC
  IMPLICIT NONE

  !For computing eigenvalues of chain EXCLUDING grafts
  INTEGER, INTENT(IN) :: iframe
  INTEGER :: i,j,a1id,u,v,initmon
  INTEGER, PARAMETER :: N_LA = 3
  INTEGER, PARAMETER :: NMAX_LA = 3
  INTEGER, PARAMETER :: LWMAX = 1000
  INTEGER :: LWORK
  INTEGER :: INFO,LWKOPT
  REAL :: SMN(N_LA,3),W(N_LA),WORK(LWMAX)
  REAL :: SMN_ALL(N_LA,3)
  REAL :: rxcm, rycm, rzcm
  REAL,DIMENSION(1:nmons_backbone,3) :: rshift_lmp 


  !Note 1: Ordered in such a way that first chains, then the graft in
  !each chain. THIS IS VERY IMPORTANT.

  !Note 2: Overall eigenvalues do not make physical sense and hence
  !not computed for MC chains taken together.

  ! Create SMN (rij matrix)
  SMN_ALL = 0.0
  initmon = 0

  IF (iframe == 1) THEN
     WRITE(eigMCwrite,'(5(A,2X))') "Timestep","ChainID","LamXsq","LamY&
          &sq","LamZsq"
  END IF

  DO i = 1,nchains

     rshift_lmp = 0.0; SMN = 0.0
     rxcm = 0.0; rycm = 0.0; rzcm = 0.0
     CALL COMPUTE_CENTROID(i,rxcm,rycm,rzcm,nmons_backbone)
          
     DO j = 1, nmons_backbone
        
        a1id = initmon + j
        rshift_lmp(j,1) = rxyz_lmp(a1id,1)-rxcm
        rshift_lmp(j,2) = rxyz_lmp(a1id,2)-rycm
        rshift_lmp(j,3) = rxyz_lmp(a1id,3)-rzcm
   
        DO u = 1,3

           DO v = 1,3

              SMN(u,v) = SMN(u,v)+rshift_lmp(j,u)*rshift_lmp(j,v)

           END DO
           
        END DO

     END DO

     SMN = SMN/REAL(nmons_backbone)


     !Compute eigenvalues ONLY for MC chains -- individually
     LWORK = -1; W = 0
     CALL DSYEV('Vectors','Upper',N_LA,SMN,N_LA,W,WORK,LWORK,INFO)
     LWORK = MIN(LWMAX,INT(WORK(1)))
     CALL DSYEV('N','U',N_LA,SMN,N_LA,W,WORK,LWORK,INFO)
     IF(INFO .NE. 0) STOP 'Physically irrelevant Eigenvalues'

     WRITE(eigMCwrite,'(2(I0,1X),3(F14.6,1X))') timestep,i,W(1),W(2)&
          &,W(3)

     initmon = initmon + nmons_backbone
     !Ordered in such a way that first chains, then the graft in each
     !chain

  END DO
          
END SUBROUTINE COMPUTE_EIGVALSOFMCCHAIN

!--------------------------------------------------------------------

SUBROUTINE COMPUTE_DENS(iframe)

  USE PARAMS_SOLVATEDMC
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: iframe
  INTEGER :: i,j, typeval, flag
  REAL :: binval, boxdir, rval, rbin
  REAL, DIMENSION(0:maxden_bin-1,ndentypes) :: inst_array
  
  inst_array = 0

  IF(dens_axis == 1) THEN
     boxdir = box_xl
  ELSEIF(dens_axis == 2) THEN
     boxdir = box_yl
  ELSEIF(dens_axis == 3) THEN
     boxdir = box_zl
  ELSE
     STOP "Unknown box direction"
  END IF

  rbin = boxdir/REAL(maxden_bin)
  normdens = normdens + (boxdir/(box_xl*box_yl*box_zl*rbin))
  denbinavg = denbinavg + rbin

!$OMP PARALLEL DO PRIVATE(i,j,typeval,flag,binval,rval)

  DO i = 1,ntotatoms

     DO j = 1,ndentypes

        flag = -1
        IF(aidvals(i,3) == dentyp_arr(j)) THEN
           
           typeval = j
           flag = 1
           EXIT
           
        END IF

     END DO

     IF(flag == 1) THEN
        IF(dens_axis == 1) THEN
           rval = rxyz_lmp(i,1) - box_xl*ANINT(rxyz_lmp(i,1)/box_xl)
           binval = FLOOR(rval/rbin)
        ELSEIF(dens_axis == 2) THEN
           rval = rxyz_lmp(i,2) - box_yl*ANINT(rxyz_lmp(i,2)/box_yl)
           binval = FLOOR(rval/rbin)
        ELSEIF(dens_axis == 3) THEN
           rval = rxyz_lmp(i,3) - box_zl*ANINT(rxyz_lmp(i,3)/box_zl)
           binval = FLOOR(rval/rbin)
        END IF

        IF(binval .LE. maxden_bin) inst_array(binval,typeval) =&
             & inst_array(binval,typeval) + 1

     END IF

  END DO

!$OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE(i,j)
  DO i = 0, maxden_bin - 1

     DO j = 1, ndentypes

        densarray(i,j) = densarray(i,j) + inst_array(i,j)

     END DO

  END DO
  
!$OMP END PARALLEL DO

END SUBROUTINE COMPUTE_DENS

!--------------------------------------------------------------------

SUBROUTINE COMPUTE_COM(tval)

  USE PARAMS_SOLVATEDMC
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: tval
  INTEGER :: i,j,atype,a1id,moltype
  REAL    :: mval
  REAL,DIMENSION(1:nchains) :: rxcm,rycm,rzcm,totmass
  REAL :: xdist, ydist, zdist, rdist

  rxcm = 0.0; rycm = 0.0; rzcm = 0.0
  totmass = 0.0
  
  DO i = 1,ntotatoms

     DO j = 1,nchains

        IF(aidvals(i,2) == j) THEN

           a1id    = aidvals(i,1)
           moltype = aidvals(i,2)
           atype   = aidvals(i,3)
           mval    = masses(atype,1)
           rxcm(moltype) = rxcm(moltype) + rxyz_lmp(a1id,1)*mval
           rycm(moltype) = rycm(moltype) + rxyz_lmp(a1id,2)*mval
           rzcm(moltype) = rzcm(moltype) + rxyz_lmp(a1id,3)*mval
           totmass(moltype) = totmass(moltype) + mval

        END IF
     
     END DO
 
  END DO

  IF(tval == 1) THEN

     WRITE(comwrite,'(8(A,2X))') "timestep","chainID","xCOM","yCOM"&
          &,"zCOM","box_xl","box_yl","box_zl"

     IF(nchains == 2) WRITE(distcomwrite,'(5(A,2X))') "timestep","xdis&
          &t","yidst","zdist","rdist"

  END IF
  
  DO i = 1,nchains
     
     rxcm(i) = rxcm(i)/totmass(i)
     rycm(i) = rycm(i)/totmass(i)
     rzcm(i) = rzcm(i)/totmass(i)
     
     WRITE(comwrite,'(2(I0,1X),6(F16.8))') timestep,i,rxcm(i),rycm(i)&
          &,rzcm(i),box_xl,box_yl,box_zl
     
  END DO

  IF(nchains == 2) THEN

     xdist = rxcm(1) - rxcm(2)
     ydist = rycm(1) - rycm(2)
     zdist = rzcm(1) - rzcm(2)

     xdist = xdist - box_xl*ANINT(xdist/box_xl)
     ydist = ydist - box_yl*ANINT(ydist/box_yl)
     zdist = zdist - box_zl*ANINT(zdist/box_zl)
     
     rdist = sqrt(xdist**2 + ydist**2 + zdist**2)
     
     WRITE(distcomwrite,'(I0,1X,4(F16.8))') timestep,xdist,ydist&
          &,zdist,rdist

  END IF


END SUBROUTINE COMPUTE_COM

!--------------------------------------------------------------------

SUBROUTINE COMPUTE_CENTROID(chval,rxcm,rycm,rzcm,nmonscomp)

  USE PARAMS_SOLVATEDMC
  IMPLICIT NONE

  INTEGER :: i,a1id,molid,atcntr
  INTEGER,INTENT(IN) :: chval,nmonscomp
  REAL,INTENT(OUT) :: rxcm,rycm,rzcm
  
  rxcm = 0.0; rycm = 0.0; rzcm = 0.0; atcntr = 0

  DO i = 1,ntotatoms

     IF(nmonscomp == nmons_backbone) THEN

        IF(aidvals(i,2) == chval .AND. aidvals(i,3) .NE. graft_type)&
             & THEN

           a1id = aidvals(i,1)
           rxcm = rxcm + rxyz_lmp(a1id,1)
           rycm = rycm + rxyz_lmp(a1id,2)
           rzcm = rzcm + rxyz_lmp(a1id,3)
           atcntr = atcntr + 1
           
        END IF

     ELSE

        IF(aidvals(i,2) == chval) THEN

           a1id = aidvals(i,1)
           rxcm = rxcm + rxyz_lmp(a1id,1)
           rycm = rycm + rxyz_lmp(a1id,2)
           rzcm = rzcm + rxyz_lmp(a1id,3)
           atcntr = atcntr + 1
           
        END IF

     END IF

  END DO




  IF(atcntr .NE. nmonscomp) THEN

     PRINT *, "Unequal numb of monomers for centroid calcaulation"
     PRINT *, "Input monomer cnt, Loop count"
     PRINT *, nmonscomp, atcntr
     STOP

  END IF

  rxcm = rxcm/REAL(nmonscomp)
  rycm = rycm/REAL(nmonscomp)
  rzcm = rzcm/REAL(nmonscomp)

END SUBROUTINE COMPUTE_CENTROID

!--------------------------------------------------------------------

SUBROUTINE COMPUTE_PERSISTENCELEN(iframe)

  USE PARAMS_SOLVATEDMC
  IMPLICIT NONE
  
  INTEGER, INTENT(IN) :: iframe
  INTEGER :: i,j,k,ierr
  INTEGER :: chcnt, binc, bfin, bvec1, bvec2, anavg
  REAL :: xdot, ydot, zdot
  REAL :: modvec1, modvec2
  REAL :: costheta
  REAL, DIMENSION(0:nanglmain) :: coschain, cosavg
  INTEGER, DIMENSION(0:nanglmain) :: norigins
  REAL, DIMENSION(1:nmons_backbone-1,nchains) :: xbvec, ybvec, zbvec
  
  IF(nmons_backbone == 0) STOP "Number of MC monomers required"

! Generate bond vectors

  DO i = 1, nchains

     DO j = 1, nmons_backbone - 1
        
        k = (i-1)*nchains + j

        xbvec(j,i) = rxyz_lmp(k+1,1) - rxyz_lmp(k,1)
        ybvec(j,i) = rxyz_lmp(k+1,2) - rxyz_lmp(k,2)
        zbvec(j,i) = rxyz_lmp(k+1,3) - rxyz_lmp(k,3)

     END DO

  END DO

! Persistence length for the MC chain
  cosavg = 0.0

  DO chcnt = 1,nchains !Chain loop (for averaging)

     coschain = 0.0; norigins = 0  
     ! default skip_tailmon = 0.
     ! binc goes from 0 to nbonds if skip_tailmon = 0 (default).
     DO binc = 0,nanglmain !shifted time origin

        bfin = nanglmain+1-binc!total bond_vecs available
        !this includes the bond vector dotted with itself

        DO bvec1 = 1+skip_tailmon, bfin+skip_tailmon!bond_vec 1
           !bvec1 corresponds to the possible range of array indices
           !for xbvec, ybvec and zbvec woth given bfin and skipmon
           bvec2 = bvec1 + binc !bond_vec 2
           modvec1 = sqrt(xbvec(bvec1,chcnt)**2 + ybvec(bvec1,chcnt)&
                &**2 + zbvec(bvec1,chcnt)**2)

           modvec2 = sqrt(xbvec(bvec2,chcnt)**2 + ybvec(bvec2,chcnt)&
                &**2 + zbvec(bvec2,chcnt)**2)
                      
           xdot = xbvec(bvec1,chcnt)*xbvec(bvec2,chcnt)
           ydot = ybvec(bvec1,chcnt)*ybvec(bvec2,chcnt)
           zdot = zbvec(bvec1,chcnt)*zbvec(bvec2,chcnt)
           
           !The definition of cos theta is such that the vectors are
           !NOT aligned with their tails together. In normal math
           !terms this will be equal to 180-theta
           costheta = (xdot + ydot + zdot)/(modvec1*modvec2)
           coschain(binc) = coschain(binc) + costheta
           norigins(binc) = norigins(binc) + 1

!!$           if (bvec1 == 1+skip_tailmon) then
!!$              print *, binc, bfin, bvec1, bvec2
!!$              pause
!!$           end if
!!$
!!$           if (binc == 0) then
!!$              print *, bvec1, bvec2
!!$              print *, costheta,modvec1, modvec2
!!$              pause
!!$           end if
           
        END DO

     END DO

     DO anavg = 0,nanglmain

        cosavg(anavg) = cosavg(anavg) + coschain(anavg)&
             &/REAL(norigins(anavg))
          
     END DO
   
  END DO

  DO i = 0,nanglmain
     
     avgtheta_main(i) = avgtheta_main(i) + cosavg(i)/(REAL(nchains))
     
  END DO

  IF(iframe == nframes) THEN

     dum_fname = "mainpersistautocf_"//trim(traj_fname)
     OPEN(unit =dumwrite,file = trim(dum_fname), status="replace",&
          & action ="write",iostat=ierr)
     
     IF(ierr /= 0) STOP "persistence time file not found"

     WRITE(dumwrite,'(3(A,5X))') "i-j","theta_(i,i+j)","theta_(i,i+j)/&
          &theta_(i,i)" 
     DO i = 0, nanglmain
        
        WRITE(dumwrite,"(I0,1X,2(F14.6,1X))") i, avgtheta_main(i)&
             &/REAL(nframes), avgtheta_main(i)/avgtheta_main(0)

     END DO

     CLOSE(dumwrite)

  END IF

END SUBROUTINE COMPUTE_PERSISTENCELEN

!--------------------------------------------------------------------

SUBROUTINE COMPUTE_PERSLENMETHOD2(iframe)

  USE PARAMS_SOLVATEDMC
  IMPLICIT NONE
  
  INTEGER, INTENT(IN) :: iframe
  INTEGER :: i,j,k,a1id,a2id,a3id,ierr,init_mon
  REAL :: xvec1, rexvec, yvec1, reyvec, zvec1, rezvec
  REAL :: xdot, ydot, zdot, cdot
  REAL :: modvec1, modvec2
  REAL :: netdot
  REAL, DIMENSION(1:nmons_backbone-1) :: dotchain, dotavg
  INTEGER, DIMENSION(1:nmons_backbone-1) :: norigins

! Persistence length for the MC chain
  IF(nmons_backbone == 0) STOP "Number of MC monomers required"

  dotchain = 0.0; dotavg = 0.0; norigins = 0
  init_mon = 0

  DO j = 1,nchains
  
     dotchain = 0

     a1id = init_mon + j
     a2id = init_mon + nmons_backbone - 1

     rexvec = rxyz_lmp(a2id,1) - rxyz_lmp(a1id,1)
     reyvec = rxyz_lmp(a2id,2) - rxyz_lmp(a1id,2)
     rezvec = rxyz_lmp(a2id,3) - rxyz_lmp(a1id,3)

!!$     rexvec = rexvec - box_xl*ANINT(rexvec/box_xl)
!!$     reyvec = reyvec - box_yl*ANINT(reyvec/box_yl)
!!$     rezvec = rezvec - box_zl*ANINT(rezvec/box_zl)

     DO i = 1,nmons_backbone-1
        
        a1id = init_mon + i
        
        DO k = i+1,nmons_backbone
           
           a2id = init_mon + k

           xvec1 = rxyz_lmp(a2id,1) - rxyz_lmp(a1id,1)
           yvec1 = rxyz_lmp(a2id,2) - rxyz_lmp(a1id,2)
           zvec1 = rxyz_lmp(a2id,3) - rxyz_lmp(a1id,3)

!!$           xvec1 = xvec1 - box_xl*ANINT(xvec1/box_xl)
!!$           yvec1 = yvec1 - box_yl*ANINT(yvec1/box_yl)
!!$           zvec1 = zvec1 - box_zl*ANINT(zvec1/box_zl)

           modvec1 = sqrt(xvec1**2 + yvec1**2 + zvec1**2)
           
           xdot = xvec1*rexvec
           ydot = yvec1*reyvec
           zdot = zvec1*rezvec

           netdot = (xdot + ydot + zdot)/(modvec1)
           dotchain(k-i) = dotchain(k-i) + netdot
           norigins(k-i) = norigins(k-i) + 1

        END DO

     END DO

     DO k = 1,nmons_backbone-1

        dotavg(k) = dotavg(k) + dotchain(k)/REAL(norigins(k))
          
     END DO
     
     init_mon = init_mon + num_atoms_perchain(j) - 1
     
  END DO
    

  DO i = 1,nmons_backbone-1
     
     avgdot_main(i) = avgdot_main(i) + dotavg(i)/(REAL(nchains))
     
  END DO

  IF(iframe == nframes) THEN

     dum_fname = "persist2autocf_"//trim(traj_fname)
     OPEN(unit =dumwrite,file = trim(dum_fname), status="replace",&
          & action ="write",iostat=ierr)
     
     IF(ierr /= 0) STOP "persistence time file not found"

     DO i = 1, nmons_backbone-1
        
        WRITE(dumwrite,"(I0,1X,2(F14.6,1X))") i, avgdot_main(i)&
             &/REAL(nframes)

     END DO

     CLOSE(dumwrite)

  END IF

END SUBROUTINE COMPUTE_PERSLENMETHOD2

!--------------------------------------------------------------------

SUBROUTINE ALLOUTPUTS()

  USE PARAMS_SOLVATEDMC
  IMPLICIT NONE

  IF(rdfcalc) THEN
     PRINT *, "Writing RDFs .."
     CALL OUTPUT_ALLRDF()
  END IF

  IF(oxycalc) THEN
     PRINT *, "Writing Chain Coordination Data .."
     CALL OUTPUT_OXYCLUS()
  END IF

  IF(rescalc) THEN
     PRINT *, "Writing Residence Time Data .."
     CALL OUTPUT_RESTIME()
  END IF

  IF(eigcalc) THEN
     PRINT *, "Writing Eigenvalue Data .."
     CALL OUTPUT_EIGTIME()
  END IF

  IF(denscalc) THEN
     PRINT *, "Writing Density Data .."
     CALL OUTPUT_DENS()
  END IF

END SUBROUTINE ALLOUTPUTS

!--------------------------------------------------------------------

SUBROUTINE OUTPUT_ALLRDF()

  USE PARAMS_SOLVATEDMC
  IMPLICIT NONE

  INTEGER :: i,j,ierr
  REAL, PARAMETER :: vconst = 4.0*pival/3.0
  REAL :: rlower,rupper,nideal,rdffrnorm

  IF(rdffreq == 1) rdffrnorm = nframes
  IF(rdffreq .NE. 1) rdffrnorm = INT(nframes/rdffreq)
  rvolavg = rvolavg/REAL(rdffrnorm)
  PRINT *, "Average volume of box", rvolavg

  dum_fname = "rdf_"//trim(adjustl(traj_fname))
  OPEN(unit = dumwrite,file =trim(dum_fname),action="write"&
       &,status="replace",iostat=ierr)

  IF(ierr /= 0) THEN
     PRINT *, "Could not open", trim(dum_fname)
  END IF

  WRITE(dumwrite,'(A,2X)',advance="no") "r"

  DO j = 1,npairs

     WRITE(dumwrite,'(2(I0,1X))',advance="no") pairs_rdf(j,1)&
          &,pairs_rdf(j,2)

  END DO

  WRITE(dumwrite,*)

  DO i = 0,rmaxbin-1

     rlower = real(i)*rbinval
     rupper = rlower + rbinval
     nideal = vconst*(rupper**3 - rlower**3)

     WRITE(dumwrite,'(F16.5,2X)',advance="no") 0.5*rbinval*(REAL(2*i&
          &+1))

     DO j = 1,npairs

        WRITE(dumwrite,'(F16.9,1X)',advance="no")rdfarray(i,j)&
             &/(rdffrnorm*nideal)

     END DO

     WRITE(dumwrite,*)
     
  END DO

  CLOSE(dumwrite)

END SUBROUTINE OUTPUT_ALLRDF

!--------------------------------------------------------------------

SUBROUTINE OUTPUT_OXYCLUS()

  USE PARAMS_SOLVATEDMC
  IMPLICIT NONE

  INTEGER :: i,j,ierr,oxynorm

  IF(oxyfreq == 1) oxynorm = nframes
  IF(oxyfreq .NE. 1) oxynorm = INT(nframes/oxyfreq)

  dum_fname = "histchainoxy_"//trim(adjustl(traj_fname))
  OPEN(unit = dumwrite,file =trim(dum_fname),action="write"&
       &,status="replace",iostat=ierr)
  IF(ierr /= 0) THEN
     PRINT *, "Could not open", trim(dum_fname)
  END IF

  DO i = 0, nchains

     WRITE(dumwrite,*) i, histoxyarray(i)/REAL(nsolvent*oxynorm)

  END DO

  CLOSE(dumwrite)


END SUBROUTINE OUTPUT_OXYCLUS

!--------------------------------------------------------------------

SUBROUTINE OUTPUT_RESTIME()

  USE PARAMS_SOLVATEDMC
  IMPLICIT NONE

  INTEGER :: i,ierr

  dum_fname = "resautocf_"//trim(traj_fname)
  OPEN(unit =dumwrite,file = trim(dum_fname), status="replace",&
       & action ="write",iostat=ierr)

  IF(ierr /= 0) STOP "residence time file not found"

  DO i = 0, nframes-1

     WRITE(dumwrite,"(I0,1X,F14.6)") i, tplot_cf(i)

  END DO

  CLOSE(dumwrite)

END SUBROUTINE OUTPUT_RESTIME

!--------------------------------------------------------------------


SUBROUTINE OUTPUT_EIGTIME()

  USE PARAMS_SOLVATEDMC
  IMPLICIT NONE

  INTEGER :: i,ierr,j

  dum_fname = "avgeigval_"//trim(traj_fname)
  OPEN(unit =dumwrite,file = trim(dum_fname), status="replace",&
       & action ="write",iostat=ierr)

  IF(ierr /= 0) STOP "avg eigenvalue file not found"

  DO i = 1, nframes

     WRITE(dumwrite,"(I0,1X,3(F14.6,1X))") timestep_arr(i),&
          & eigarray(i,1),eigarray(i,2), eigarray(i,3)

  END DO

  CLOSE(dumwrite)

END SUBROUTINE OUTPUT_EIGTIME

!--------------------------------------------------------------------

SUBROUTINE OUTPUT_DENS()

  USE PARAMS_SOLVATEDMC
  IMPLICIT NONE

  INTEGER :: i,j,ierr
  REAL :: rlower,rupper,nideal,densfrnorm
  REAL, DIMENSION(1:ndentypes) :: frac

  IF(densfreq == 1) densfrnorm = nframes
  IF(densfreq .NE. 1) densfrnorm = INT(nframes/densfreq)

  normdens = normdens/REAL(densfrnorm)
  denbinavg = denbinavg/REAL(densfrnorm)
  frac = 0
  
  DO i = 1,ntotatoms
     
     DO j = 1,ndentypes

        IF(aidvals(i,3) == dentyp_arr(j)) THEN
           
           frac(j) = frac(j) + 1

        END IF

     END DO

  END DO

  frac = frac/REAL(ntotatoms)

  dum_fname = "dens_"//trim(adjustl(traj_fname))
  OPEN(unit = dumwrite,file =trim(dum_fname),action="write"&
       &,status="replace",iostat=ierr)

  IF(ierr /= 0) THEN
     PRINT *, "Could not open", trim(dum_fname)
  END IF

  WRITE(dumwrite,'(A,2X)',advance="no") "r"

  DO j = 1,ndentypes

     WRITE(dumwrite,'(I0,1X)',advance="no") dentyp_arr(j)

  END DO

  WRITE(dumwrite,*)

  DO i = 0,maxden_bin-1

     WRITE(dumwrite,'(F16.5,2X)',advance="no") 0.5*denbinavg*(REAL(2&
          &*i+1))

     DO j = 1,npairs

        WRITE(dumwrite,'(F16.9,1X)',advance="no") frac(j)*densarray(i&
             &,j)*normdens

     END DO

     WRITE(dumwrite,*)
     
  END DO

  CLOSE(dumwrite)

END SUBROUTINE OUTPUT_DENS

!--------------------------------------------------------------------

SUBROUTINE ALLOCATE_ARRAYS()

  USE PARAMS_SOLVATEDMC
  IMPLICIT NONE

  INTEGER :: AllocateStatus

! Allocate LAMMPS Structure

  ALLOCATE(aidvals(ntotatoms,3),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate aidvals"
  ALLOCATE(rxyz_lmp(ntotatoms,3),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate rxyz_lmp"
  ALLOCATE(charge_lmp(ntotatoms,1),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate charge_lmp"
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


  ALLOCATE(num_atoms_perchain(nchains),stat=AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate num_atoms_perchain"

  IF(rdfcalc) THEN
     ALLOCATE(rdfarray(0:rmaxbin-1,npairs),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate rdfarray"
  ELSE
     ALLOCATE(rdfarray(1,1),stat = AllocateStatus)
     DEALLOCATE(rdfarray)
  END IF
  
  IF(oxycalc .NE. 1) THEN
     ALLOCATE(oxyarray(1),stat = AllocateStatus)
     DEALLOCATE(oxyarray)
     ALLOCATE(histoxyarray(1),stat = AllocateStatus)
     DEALLOCATE(histoxyarray)
  END IF

  IF(rescalc .NE. 1) THEN
     ALLOCATE(oxyresarray(1),stat = AllocateStatus)
     DEALLOCATE(oxyresarray)
     ALLOCATE(autocf(1,1),stat = AllocateStatus)
     DEALLOCATE(autocf)
     ALLOCATE(tplot_cf(1),stat = AllocateStatus)
     DEALLOCATE(tplot_cf)
  END IF


  IF(resegcalc .NE. 1) THEN
     ALLOCATE(resegarr(1),stat = AllocateStatus)
     DEALLOCATE(resegarr)
  END IF
     
     

  IF(perscalc .NE. 1) THEN
     ALLOCATE(avgtheta_main(1), stat=AllocateStatus)
     DEALLOCATE(avgtheta_main)
     ALLOCATE(avgdot_main(1), stat=AllocateStatus)
     DEALLOCATE(avgdot_main)
  END IF
     
  IF(eigcalc .NE. 1) THEN
     ALLOCATE(eigarray(1,1),stat = AllocateStatus)
     DEALLOCATE(eigarray)
     ALLOCATE(timestep_arr(1),stat = AllocateStatus)
     DEALLOCATE(timestep_arr)
  END IF

  IF(denscalc) THEN
     ALLOCATE(densarray(0:maxden_bin-1,ndentypes),stat =&
          & AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate densarray"
  ELSE
     ALLOCATE(densarray(1,1),stat = AllocateStatus)
     DEALLOCATE(densarray)
  END IF
     

  PRINT *, "Successfully allocated memory"

END SUBROUTINE ALLOCATE_ARRAYS

!--------------------------------------------------------------------

SUBROUTINE DEALLOCATE_ARRAYS()

  USE PARAMS_SOLVATEDMC

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

