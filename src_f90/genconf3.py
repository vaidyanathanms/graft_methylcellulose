# To generate initial configurations for semiflexible chains
#!C: Version: V_Feb_14_2019
# New version: V_Mar_21_2020 see #!C for major change in algorithm
# New version looks for restart files rather than dump files

import numpy
import os
import shutil
import subprocess
import sys
import glob
import re

from subprocess import call
from my_python_functions import find_recent_file
from my_python_functions import cpy_lammps_files
from my_python_functions import cpy_main_files
from my_python_functions import generate_backup
from my_python_functions import create_infile
from my_python_functions import check_integer
from resubmit import start_from_beginning
from resubmit import run_long_equil_cycle
from resubmit import run_production_cycle

#------------------Input Options---------------------------------


# 0 - No grafts
# 1 - Make PEG as MC substitute (implicit solv)
# 2 - Grow side chains for MC (implicit solv)
# 3 - Grow side chains for semiflexible chains (explicit solv)
# 3.1 - Fredrickson's case (explicit/implicit solv)

graftopt = 2 # see options above
wlcchec  = 1 # only valid when graftopt = 3, simple WLC model
mindump  = 20000000 # end point for equilibrium
maxdump  = 80000000 # maximum dump time
anglstyl = 'harmonic' # cosine or harmonic
config   = 5 # Trial number
if graftopt != 3:
    wlccheck = 0

#--- Prof. Ilja's suggestion: Use genconf_ilja.py---------------
key_eps_pg = 0 # 
epsval_pg  = 0 #

#-----------------Input Arrays------------------------------------

#Interaction details

epsarr_gg   = [0.8,1.0,1.2]#,1.0,1.2]  # graft-graft
epsarr_sg   = [1.0,1.0,1.2]#,1.0,1.0]  # solvent-graft 
epsarr_ps   = [1.0,1.0,1.2]#,1.0,1.0]  # polymer-solvent

#epsarr_gg   = [1.0]#,1.2]#,1.2]  # graft-graft
#epsarr_sg   = [1.0]#,1.0]#,1.0]  # solvent-graft
#epsarr_ps   = [1.0]#,1.0]#,1.0]  # polymer-solvent
rcut        = pow(2,1/6)
#
#Chain and initial configuration details
nchains     = 1 # Number of backbone chains
nmons       = [1000]#,500,1000,2000]#,1000,2000] # Number of backbone monomers
initcom     = 20.0 # Only for 2 chain systems - d_COM as fn(t)
graftMW     = 25  # number of graft monomers per graft
polywtperc  = 1.0 # tot_poly wt% - explicit generic 
polydens    = 0.5 # overall density of polymers

#Graft percentage/chain
#graftarr    = [0.01,0.05,0.10,0.15,0.20,0.25,0.30] 
#graftarr    = [0.20,0.25,0.30]
graftarr    = [0.2]
#graftarr     = [0.0]
#graftarr    = [0.0,0.25,0.30]#,0.03,0.05,0.10,0.15,0.25,0.30]
#For methylcellulose only
DS_MC   = '1.80' # String Value - only for methylcellulose
temp    = '50.0' # String Value - for computing parameters for MC

#dump and other filenames
dumpname = 'dump_stage_*' # include the *
restname = 'archive_*'    # include the *

#---------------Topology Info ------------------------------------

# All constants are for generic semiflexible system.
# For MC it is in built
k_b = 500.0  # Stretching constant
r_b = 1.0    # Equilibrium distance
k_t = 10.0   # Bending constant
teq = 0.0  # Equilibrium angle
k_phi = 0.0  # Dihedral constant

blist = [k_b, r_b]
alist = [k_t, teq]

#----------------Directory Settings-------------------------------

maindir     = os.getcwd() # Also files where f90 files are present
scratchdir  = '/scratch.global/vaidya/' # Directory for calculations
#srclmp - LAMMPS input files directory
srclmp      = '/home/dorfmank/vsethura/allfiles/files_graft/src_lmp'
lmp_execdir = '/home/dorfmank/vsethura/mylammps/src' #LAMMPS src

if not os.path.isdir(scratchdir):
    print(scratchdir, "does not exist")
    sys.exit()

if not os.path.isdir(srclmp):
    print(srclmp, "does not exist")
    sys.exit()

if not os.path.isdir(lmp_execdir):
    print(lmp_execdir, "does not exist")
    sys.exit()


#-----------------File List---------------------------------------

geninp_list = ('lmp_params.f90','lammps_inp2.f90','ran_numbers.f90',
              'params_input_var.dat')

#--------------Main Analysis--------------------------------------

for bblen in range(len(nmons)): #Backbone length loop

    print( "Backbone length: ", nmons[bblen])
    workdir_main = scratchdir + 'grafts_all' 
    if not os.path.isdir(workdir_main):
        os.mkdir(workdir_main)

    if graftopt == 0:
        workdir_main = workdir_main + '/no_grafts' 
        print ( "No graft system")
    elif graftopt == 1:
        workdir_main = workdir_main + '/substitute_grafts_MC' 
        print ( "Substituted graft system")
    elif graftopt == 2:
        workdir_main = workdir_main + '/grafts_MC' 
        print ( "Side chain graft system")
    elif graftopt == 3:
        if anglstyl == 'cosine':
            workdir_main = workdir_main + '/semiflex_grafts_cosstyle' 
        else:
            workdir_main = workdir_main + '/semiflex_grafts' 

        print ( "Semiflexible graft system")
    elif graftopt == 3.1:
        workdir_main = workdir_main + '/flex_bb' 
        print ( "Flexible backbone system")

    if not os.path.isdir(workdir_main):
        os.mkdir(workdir_main)

    #Make number of monomers in backbone directory
    workdir_bb_main = workdir_main + '/n_bb_' + str(nmons[bblen])
    if not os.path.isdir(workdir_bb_main):
        os.mkdir(workdir_bb_main)

    #Make temperature directory
    if int(graftopt) != 3:
        workdir_temp = workdir_bb_main + '/temp_' + temp
        if not os.path.isdir(workdir_temp):
            os.mkdir(workdir_temp)
    else:
        workdir_temp = workdir_bb_main

    #Make trial configuration directory
    workdir_config = workdir_temp + '/config_' + str(config)
    if not os.path.isdir(workdir_config):
        os.mkdir(workdir_config)

    #Make number of chains directory
    workdir_chain = workdir_config + '/nchains_' + str(nchains)
    if not os.path.isdir(workdir_chain):
        os.mkdir(workdir_chain)

    #Make wt% directory
    workdir_bb  = workdir_chain + '/backbonewtperc_'+ str(polydens)
    if not os.path.isdir(workdir_bb):
        os.mkdir(workdir_bb)

    #Make graftMW directory
    workdir_graft = workdir_bb + '/n_graft_' + str(graftMW)
    if not os.path.isdir(workdir_graft):
        os.mkdir(workdir_graft)

    for glen in range(len(graftarr)): #Graft loop

        print( "Graft Percentage: ", graftarr[glen])            

        workdir1 = workdir_graft + '/graftperc_' + str(graftarr[glen])

        #Continue iff at least one graft is present
        num_graft_molecules = int(graftarr[glen]*nmons[bblen])

        print("Number of graft molecules: ",num_graft_molecules)
#        if num_graft_molecules < 1 and wlccheck == 0:
#            print("number of graft less than 1 for sigma:", graftarr[glen])
#            continue

        if not os.path.isdir(workdir1):
            os.mkdir(workdir1)


        for eps in range(len(epsarr_gg)): #eps_pg loop
            
            if polywtperc != 1.0:
                data_fname = "input_" + str(epsarr_ps[eps]) + "_" + \
                             str(epsarr_sg[eps]) + ".dat"
            else:
                data_fname = "input_" + str(epsarr_gg[eps]) + ".dat"
                
            print( "EpsValue_polymer_graft: ", epsarr_gg[eps])

            workdir3 = workdir1 + '/epsval_' + str(epsarr_gg[eps])

            if not os.path.isdir(workdir3):
                os.mkdir(workdir3)

            os.chdir(workdir3)
            destdir = os.getcwd()
            
            print ("Currently in: ", destdir)


            list_of_files = glob.glob('*')


            if list_of_files == []:
                begin_fresh = 1 #No files here to start with
                print("No files found in: ", destdir)
                print("Beginning new simulation")
                start_from_beginning(maindir,srclmp,lmp_execdir,destdir,\
                                     data_fname,temp,polywtperc,rcut,blist,\
                                     alist,k_phi,wlccheck,anglstyl,nchains,\
                                     polydens,graftMW,initcom,nmons[bblen],\
                                     graftarr[glen],epsarr_gg[eps],graftopt,\
                                     geninp_list,begin_fresh,key_eps_pg,\
                                     epsval_pg)
                
                
            else:

                restflag = 1 # restart file found = 1                
                latest_restfyl = find_recent_file(destdir,restname)
#!C:            trajflag = 1 # trajectory found=1
#!C:            flagstr  = -1 # is a string = 1
#!C:            latest_trajfyl = find_recent_file(destdir,
#!C:                                               dumpname)

                
                #!C: dumpfile should be of type dump_stage_*
                #restart file should be of type archive_*

                # check whether archive files are present
                if latest_restfyl == "nil": #!C: latest_trajfyl == "nil"
                    restflag = -1 #!C: trajflag = -1
                else:
                    delimited_vals = re.split("\W+|_",latest_restfyl)
                    #!C re.split("\W+|_",latest_trajfyl)
                    timeval = delimited_vals[len(delimited_vals)-1] #!C: [len(delimited_vals)-2]
                    #!C: check whether the dumpfile is dump_stage1.*
                    #!C: flagstr = check_integer(timeval)
                    #!C: if flagstr == 1:
                    #!C:    print("Did not find matching dumpfile, restarting")
                    #!C:    trajflag = -1
                    if int(timeval) < 1000000: #!C: elif
                        print("Timestep less than 1000000, restarting")
                        restflag = -1 #!C: trajflag = -1
                    else:
                        print("Latest timestep: ", timeval)


                if restflag == -1: #!C: trajflag == -1

                    #Files already here needs to be deleted

                    begin_fresh = -1 
                    
#!C:                    if flagstr == -1:
#!C:                        print("No trajectory file in")
#!C:                    else:
#!C:                        print("Stopped before 1 million steps")
                    print("No restart files/ stopped before 1 million steps")
                    print("Restarting simulation...")

                    start_from_beginning(maindir,srclmp,lmp_execdir,destdir,\
                                         data_fname,temp,polywtperc,rcut,blist,\
                                         alist,k_phi,wlccheck,anglstyl,nchains,\
                                         polydens,graftMW,initcom,nmons[bblen],\
                                         graftarr[glen],epsarr_gg[eps],graftopt,\
                                         geninp_list,begin_fresh,key_eps_pg,\
                                         epsval_pg)
                    
                elif int(timeval) > maxdump:
                    print("Maximum timestep reached")
                    continue

                elif int(timeval) > mindump:
                    print("Running production cycle")
                    run_production_cycle(maindir,destdir,srclmp,nmons[bblen],\
                                         graftarr[glen],graftopt)

                else:
                    print("Running long equil+prod cycles for: ",destdir)
                    run_long_equil_cycle(maindir,destdir,srclmp,nmons[bblen],\
                                         graftarr[glen],graftopt)


                

                
