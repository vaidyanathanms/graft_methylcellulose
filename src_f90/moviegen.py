# To generate initial configurations for semiflexible chains
# Version: V_Feb_14_2019

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
from resubmit import run_movie_cycle

#------------------Input Options---------------------------------


# 0 - No grafts
# 1 - Make PEG as MC substitute (implicit solv)
# 2 - Grow side chains for MC (implicit solv)
# 3 - Grow side chains for semiflexible chains (explicit solv)
# 3.1 - Fredrickson's case (explicit/implicit solv)

graftopt = 2 # see options above
wlccheck = 1 # only valid when graftopt = 3, simple WLC model
mindump  = 20000000 # end point for equilibrium
anglstyl = 'harmonic' # cosine or harmonic
config   = 5 # Trial number
if graftopt != 3:
    wlccheck = 0

#--- Prof. Ilja's suggestion: Use genconf_ilja.py---------------
key_eps_pg = 0 # 
epsval_pg  = 0 #

#-----------------Input Arrays------------------------------------

#Interaction details

#epsarr_gg   = [1.2]#,1.0,1.2]  # graft-graft
#epsarr_sg   = [1.0]#,1.0,1.0]  # solvent-graft 
#epsarr_ps   = [1.0,1.0,1.0]  # polymer-solvent

epsarr_gg   = [0.8,1.0,1.2]  # graft-graft
epsarr_sg   = [0.8,1.0,1.0]  # solvent-graft
epsarr_ps   = [0.8,1.0,1.0]  # polymer-solvent
rcut        = pow(2,1/6)

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
#graftarr    = [0.0]
graftarr     = [0.0,0.01,0.03,0.05,0.10,0.15,0.20,0.25,0.30]
#graftarr    = [0.0,0.25,0.30]#,0.03,0.05,0.10,0.15,0.25,0.30]
#For methylcellulose only
DS_MC   = '1.80' # String Value - only for methylcellulose
temp    = '50.0' # String Value - for computing parameters for MC

#dump and other filenames
dumpname = 'dump_stage_*' #include the *

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
        print(workdir_bb_main," does not exist")
        continue

    #Make temperature directory
    if int(graftopt) != 3:
        workdir_temp = workdir_bb_main + '/temp_' + temp
        if not os.path.isdir(workdir_temp):
            print(workdir_temp, " does not exist")
            continue
    else:
        workdir_temp = workdir_bb_main

    #Make trial configuration directory
    workdir_config = workdir_temp + '/config_' + str(config)
    if not os.path.isdir(workdir_config):
        print(workdir_config, " does not exist")
        continue

    #Make number of chains directory
    workdir_chain = workdir_config + '/nchains_' + str(nchains)
    if not os.path.isdir(workdir_chain):
        print(workdir_chain, " does not exist")
        continue


    #Make wt% directory
    workdir_bb  = workdir_chain + '/backbonewtperc_'+ str(polydens)
    if not os.path.isdir(workdir_bb):
        print(workdir_bb, " does not exist")
        continue

    #Make graftMW directory
    workdir_graft = workdir_bb + '/n_graft_' + str(graftMW)
    if not os.path.isdir(workdir_graft):
        print(workdir_graft, " does not exist")
        continue

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
            print(workdir1, " does not exist")
            continue


        for eps in range(len(epsarr_gg)): #eps_pg loop
            
            if polywtperc != 1.0:
                data_fname = "input_" + str(epsarr_ps[eps]) + "_" + \
                             str(epsarr_sg[eps]) + ".dat"
            else:
                data_fname = "input_" + str(epsarr_gg[eps]) + ".dat"
                
            print( "EpsValue_polymer_graft: ", epsarr_gg[eps])

            workdir3 = workdir1 + '/epsval_' + str(epsarr_gg[eps])

            if not os.path.isdir(workdir3):
                print(workdir3, " does not exist")
                continue

            os.chdir(workdir3)
            destdir = os.getcwd()
            
            print ("Currently in: ", destdir)


            list_of_files = glob.glob('*')


            if list_of_files == []:
                print("No files here: Check path")
                continue
                
            else:
                
                trajflag = 1 #trajectory found=1
                flagstr  = -1 #is a string = 1
                latest_trajfyl = find_recent_file(destdir,
                                                  dumpname)
                #dumpfile should be of type dump_stage_*
                if latest_trajfyl == "nil":
                    trajflag = -1
                else:
                    delimited_vals = re.split("\W+|_",latest_trajfyl)
                    timeval = delimited_vals[len(delimited_vals)-2]
                    #check whether the dumpfile is dump_stage1.*
                    flagstr = check_integer(timeval)
                    if flagstr == 1:
                        print("Did not find matching dumpfile, restarting")
                        trajflag = -1
                    elif int(timeval) < 1000000:
                        print("Timestep less than 1000000, restarting")
                        trajflag = -1
                    else:
                        print("Latest timestep: ", timeval)


                if trajflag == -1:

                    print("Not equilibrated: Run more")
                    
                elif int(timeval) > mindump:
                    print("Compiling movie files for: ",destdir)
                    fig_dir = destdir + '/fig_files'
                    if not os.path.isdir(fig_dir):
                        os.mkdir(fig_dir)
                    run_movie_cycle(maindir,destdir,srclmp,nmons[bblen],\
                                   graftarr[glen],graftopt)

                else:
                    print("yet to run production cycle")

                

                
