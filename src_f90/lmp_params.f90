!-----To generate input parameters for LAMMPS CG-MethylCellulose-----
!----- Main File: lammps_inp.f90-------------------------------------
!------Version: Jan-25-2018------------------------------------------
!--------------------------------------------------------------------

MODULE PARAMS

  USE RAN_NUMBERS

  IMPLICIT NONE

  ! Parameter data for creating the data file
  ! All the parameters can be reset

  INTEGER :: N = 0 
  INTEGER :: M = 0
  REAL :: DS_MC = 1.8
  REAL :: perc_graft_per_chain = 0.0
  REAL :: mass_gr2bb = 1.0
  REAL :: polywtperc = 1.0
  INTEGER :: ngraft_mons = 0
  INTEGER :: ngraft_per_bb = 0
  INTEGER :: totpart = 0
  INTEGER :: bondcalc = 0
  INTEGER :: nsolvent = 0
  INTEGER :: total_poly = 0
  REAL    :: init_com_dist = 0.0
  REAL    :: com_tol = 1.0

  ! Box details-- All parameters can be reset or else use this

  REAL :: insidebox = 100.0 !start of all polymers
  REAL :: density = 5000.0/(600.0**3.0)
  REAL :: boxl_x, boxl_y, boxl_z
  REAL :: volbox

  ! Math Constants

  REAL, PARAMETER :: math_pi = 3.14159265359

  ! MC constants

  REAL,PARAMETER :: DS_fac = 2.1
  REAL, PARAMETER :: mean_tol = 0.01

  ! Chain Constants

  REAL, PARAMETER :: r0init  = 0.97
  REAL, PARAMETER :: r0sq3   = r0init/sqrt(3.0)
  REAL, PARAMETER :: rmaxsq  = r0init*r0init

  ! Flags for creating the data file

  !grafts=-1 => only atoms
  !grafts=0 => no grafts
  !grafts=1 => MC with subsituted monomers
  !grafts=2 => MC with graft
  !grafts=3 => Semiflexible chain with graft (not MC)

  INTEGER :: grafts = 3
  INTEGER :: numatomtypes = 2
  INTEGER :: numbondtypes = 1
  INTEGER :: numangltypes = 1
  INTEGER :: numdihdtypes = 1

! Global Arrays involved in creating data file
  
  REAL,ALLOCATABLE,DIMENSION(:,:) :: rxyz, uxyz
  INTEGER,ALLOCATABLE,DIMENSION(:,:) :: im_xyz
  INTEGER,ALLOCATABLE,DIMENSION(:,:) :: atype
  INTEGER,ALLOCATABLE,DIMENSION(:) :: cntatomtype
  REAL,    ALLOCATABLE, DIMENSION(:,:) :: rgrafts,rsolvent
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: im_xyzgraft,agraftid&
       &,graftmonvals
  REAL,    ALLOCATABLE, DIMENSION(:) :: mass_arr

  ! Character Arrays for creating the data file name and fileptr
  
  CHARACTER (LEN = 256):: datafile, dumchar,ana_fname, log_fname
  INTEGER, PARAMETER :: outfile = 4,logfile = 8,anaread=6

  ! Random number details

  TYPE (RAN_SAVE) :: X
  INTEGER(4) :: S
  integer :: init_rannum = -1
  
END MODULE PARAMS
