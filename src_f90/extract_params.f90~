!---------------To extract files for US after SMD--------------------
!---------------Version 1: Mar-29-2018-------------------------------
!---------------Main File: extract_conf.f90--------------------------
!********************************************************************

MODULE EXTRACT_PARAMS

  USE OMP_LIB
  IMPLICIT NONE

  ! Required Input Variables

  INTEGER :: initdist
  INTEGER :: nframes, skipfr, freqfr, nfrcntr
  INTEGER :: nproc

  ! Logical Flags
  
  INTEGER :: multiposflag, comflag, netgenflag

  ! Variables for extracting files

  INTEGER :: g1_init, g1_fin, g2_init, g2_fin, g1_nmons, g2_nmons
  INTEGER :: pullaxis
  REAL :: savedist, tol

  ! File names and unit Numbers
  
  CHARACTER(LEN = 256) :: ana_fname,data_fname,traj_fname,log_fname
  CHARACTER(LEN = 256) :: dum_fname
  INTEGER, PARAMETER :: anaread = 2,logout = 3, dumwrite = 200
  INTEGER, PARAMETER :: outfile = 250, inpread = 100,comout = 210

  !Math Constants

  REAL*8, PARAMETER :: pival  = 3.14159265359
  REAL*8, PARAMETER :: pi2val = 2.0*pival

  !Global analysis variables and arrays

  INTEGER :: atomflag, velflag, bondflag, anglflag, dihdflag,imprflag
  INTEGER :: ntotatoms, ntotbonds, ntotangls,ntotdihds,ntotimprs
  INTEGER :: ntotatomtypes,ntotbondtypes,ntotangltypes,ntotdihdtypes&
       &,ntotimprtypes

  !Lammps trajectory file read details

  REAL :: box_xl,box_yl,box_zl, boxval
  INTEGER*8 :: timestep

  !Structural Average Variables

  REAL :: re2ave, re4ave, rg2ave, rg4ave, b2ave
  
  !Required Arrays - LAMMPS

  REAL*8, ALLOCATABLE, DIMENSION(:,:) :: rxyz_lmp, vel_xyz, charge_lmp&
       &,masses
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: bond_lmp, angl_lmp,&
       & dihd_lmp, impr_lmp,aidvals
  CHARACTER,ALLOCATABLE,DIMENSION(:,:) :: keywords
  REAL,ALLOCATABLE,DIMENSION(:):: boxx_arr, boxy_arr,boxz_arr

  !Required Arrays - Structural Quantities

  INTEGER, ALLOCATABLE,DIMENSION(:):: grp1_array, grp2_array
  
END MODULE EXTRACT_PARAMS
!--------------------------------------------------------------------
