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

#------------------Input Options---------------------------------

# 0 - No grafts
# 1 - Make PEG as MC substitute (implicit solv)
# 2 - Grow side chains for MC (implicit solv)
# 3 - Grow side chains for semiflexible chains (explicit solv)
# 3.1 - Fredrickson's case (explicit/implicit solv)

graftopt = 2 # see options above
wlccheck = 1 # only valid when graftopt = 3, simple WLC model
restart  = 0 # 0-From beginning, 1-long_equil_prod, 2-long_prod
mindump  = 20000000 # only applicable for restart cases
anglstyl = 'harmonic' #cosine or harmonic
config   = 2 # Trial number

if graftopt != 3:
    wlccheck = 0

#-----------------Input Arrays------------------------------------

#Interaction details

#epsarr_pg   = [1.2]  # polymer-graft
#epsarr_ps   = [1.0]  # polymer-solvent
#epsarr_sg   = [1.0]  # solvent-graft


epsarr_pg   = [0.8,1.0,1.2]  # polymer-graft
epsarr_ps   = [1.0,1.0,1.2]  # polymer-solvent
epsarr_sg   = [1.0,1.0,1.2]  # solvent-graft
rcut        = pow(2,1/6)

#Chain and initial configuration details
nchains     = 1 # Number of backbone chains
nmons       = [1000]#,500,1000,2000]#,1000,2000] # Number of backbone monomers
initcom     = 20.0 # Only for 2 chain systems - d_COM as fn(t)
graftMW     = 25  # number of graft monomers per graft
polywtperc  = 1.0 # tot_poly wt% - explicit generic 
polydens    = 0.5 # overall density of polymers

#Graft percentage/chain
graftarr    = [0.01,0.05,0.10,0.15,0.20,0.25,0.30] 

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


        for eps in range(len(epsarr_pg)): #eps_pg loop
            
            if polywtperc != 1.0:
                data_fname = "input_" + str(epsarr_ps[eps]) + "_" + \
                             str(epsarr_sg[eps]) + ".dat"
            else:
                data_fname = "input_" + str(epsarr_pg[eps]) + ".dat"
                
            print( "EpsValue_polymer_graft: ", epsarr_pg[eps])

            workdir3 = workdir1 + '/epsval_' + str(epsarr_pg[eps])

            if not os.path.isdir(workdir3):
                os.mkdir(workdir3)

            os.chdir(workdir3)
            destdir = os.getcwd()


            if restart == 0: 
                
                for fylcnt in range(len(geninp_list)):
                    cpy_main_files(maindir, destdir,geninp_list[fylcnt])

                # Manipulate Files
                os.chdir(destdir)
                print( "Manipulating input parameter file..")
                
                launch_fyl = 'params_input.dat'
                fr = open('params_input_var.dat','r')
                fw = open(launch_fyl,'w')
                fid = fr.read().replace("py_inputname",str(data_fname)).\
                      replace("py_graftstyle",str(int(graftopt))).\
                      replace("py_nchains_bb",str(nchains)).\
                      replace("py_nmons_bb",str(nmons[bblen])).\
                      replace("py_graftperc",str(graftarr[glen])).\
                      replace("py_polydens",str(polydens)).\
                      replace("py_graftmons",str(graftMW)).\
                      replace("py_comdist",str(initcom)).\
                      replace("py_polyperc",str(polywtperc))
                fw.write(fid)
                fw.close()
                fr.close()

                # Generate Required Files for LAMMPS
                
                print( "Running FORTRAN script for generating Datafile")
                subprocess.call(["ifort","-r8","-qopenmp","-mkl",\
                                 "ran_numbers.f90","lmp_params.f90",\
                                 "lammps_inp2.f90","-o","lmpinp.o"])
                

                subprocess.call(["./lmpinp.o","params_input.dat"])                    

                if not os.path.exists(data_fname):
                    print(data_fname, "not found - check input code")
                    continue

                cpy_lammps_files(maindir,srclmp,lmp_execdir,destdir,graftopt,
                                 data_fname,nmons[bblen],graftarr[glen],anglstyl)
                create_infile(maindir,srclmp,destdir,nmons[bblen],
                              epsarr_pg[eps],graftopt,temp,graftarr[glen],
                              polywtperc,rcut,blist,alist,k_phi,wlccheck,anglstyl)
                generate_backup(maindir,srclmp,lmp_execdir,destdir)

            elif restart == 1:
                
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
                    elif int(timeval) < 100000:
                        print("Timestep less than 1000000, restarting")
                        trajflag = -1
                    else:
                        print("Latest timestep: ", timeval)


                if trajflag == -1:

                    print("Could not find a trajectory file in "
                          ,destdir,"\n")
                    print("Restarting simulation...")

                    list_of_files = glob.glob('*')
                    if list_of_files != []:
                        for f in list_of_files:
                            if not os.path.isdir(f):
                                os.remove(f)

                    for fylcnt in range(len(geninp_list)):
                        cpy_main_files(maindir, destdir,geninp_list[fylcnt])

                    # Manipulate Files
                    os.chdir(destdir)
                    print( "Manipulating input parameter file..")
                
                    launch_fyl = 'params_input.dat'
                    fr = open('params_input_var.dat','r')
                    fw = open(launch_fyl,'w')
                    fid = fr.read().replace("py_inputname",str(data_fname)).\
                          replace("py_graftstyle",str(int(graftopt))).\
                          replace("py_nchains_bb",str(nchains)).\
                          replace("py_nmons_bb",str(nmons[bblen])).\
                          replace("py_graftperc",str(graftarr[glen])).\
                          replace("py_polydens",str(polydens)).\
                          replace("py_graftmons",str(graftMW)).\
                          replace("py_comdist", str(initcom)).\
                          replace("py_polyperc",str(polywtperc))
                    fw.write(fid)
                    fw.close()
                    fr.close()
                    
                    # Generate Required Files for LAMMPS
                    
                    print( "Running FORTRAN script for generating Datafile")
                    subprocess.call(["ifort","-r8","-qopenmp","-mkl",\
                                     "ran_numbers.f90","lmp_params.f90",\
                                     "lammps_inp2.f90","-o","lmpinp.o"])
                    
                    subprocess.call(["./lmpinp.o","params_input.dat"])

                    
                    if not os.path.exists(data_fname):
                        print(data_fname, "not found - check input code")
                        continue

                    cpy_lammps_files(maindir,srclmp,lmp_execdir,destdir,graftopt,
                                     data_fname,nmons[bblen],graftarr[glen],anglstyl)
                    create_infile(maindir,srclmp,destdir,nmons[bblen],
                                  epsarr_pg[eps],graftopt,temp,graftarr[glen],
                                  polywtperc,rcut,blist,alist,k_phi,wlccheck,anglstyl)
                    generate_backup(maindir,srclmp,lmp_execdir,destdir)

                    
                elif int(timeval) > mindump:
                        print("System already at production stage \n")
                        continue

                else:
                    print("Copying data for", destdir)

                    if int(graftopt) == 3:
                        srcfyl = srclmp + '/in.equil2'
                    else:
                        srcfyl = srclmp + '/in.mc_equil2'

                    desfyl = destdir + '/in.equil2'
                    shutil.copy2(srcfyl, desfyl)

                    if int(graftopt) == 3:
                        srcfyl = srclmp + '/in.equil3'
                    else:
                        srcfyl = srclmp + '/in.mc_equil3'

                    desfyl = destdir + '/in.equil3'
                    shutil.copy2(srcfyl, desfyl)

                    if int(graftopt) == 3:
                        srcfyl = srclmp + '/in.prod'
                    else:
                        srcfyl = srclmp + '/in.mc_prod'
                    desfyl = destdir + '/in.prod'
                    shutil.copy2(srcfyl, desfyl)

                    restartset = 'restart*'
                    restartlist = glob.glob(restartset)
                    if restartlist != []:
                        srcfyl = srclmp + '/jobmain_long_var.sh'
                        desfyl = destdir + '/jobmain_long_var.sh'
                        shutil.copy2(srcfyl, desfyl)
                        mainfyl = 'jobmain_long_var.sh'
                        print("found restart files in", destdir)
                    else:
                        srcfyl = srclmp + '/jobmain_var.sh'
                        desfyl = destdir + '/jobmain_var.sh'
                        shutil.copy2(srcfyl, desfyl)
                        mainfyl = 'jobmain_var.sh'
                        print("Did not find restart files in", destdir)

                    launch_fyl = 'jobmain_long.sh'
                    fr = open(mainfyl,'r')
                    fw = open(launch_fyl,'w')
                    jobstr = "joblong_"+ str(nmons[bblen]) + "_" \
                             + str(graftarr[glen])
                    fid = fr.read().replace("py_jobname",jobstr)
                    fw.write(fid)
                    fw.close()
                    fr.close()
                    
                    srcfyl = '/home/dorfmank/vsethura/mylammps/src/lmp_mesabi'
                    desfyl = destdir + '/lmp_mesabi'
                    shutil.copy2(srcfyl, desfyl)
                    print( "Submitting Jobs..")
                    
                    subprocess.call(["qsub","jobmain_long.sh"])                
                    print("\n")
                    os.chdir(maindir)
                
            else:


                trajflag = 1 #trajectory found=1
                flagstr  = -1 #is a string = 1
                latest_trajfyl = find_recent_file(destdir,dumpname)

                #dumpfile should be of type dump_stage_*
                if latest_trajfyl == "nil":
                    print("Here")
                    trajflag = -1
                else:
                    delimited_vals = re.split("\W+|_",latest_trajfyl)
                    timeval = delimited_vals[len(delimited_vals)-2]
                    #check whether the dumpfile is dump_stage1.*
                    flagstr = check_integer(timeval)
                    if flagstr == 1:
                        print("Did not find matching dumpfile, do rerun")
                        trajflag = -1
                    elif int(timeval) < 20000000:
                        print("Timestep less than 20000000, do rerun")
                        trajflag = -1
                    else:
                        print("Latest timestep: ", timeval)


                if trajflag == -1:
                    print ("Here")
                    continue

                if int(graftopt) == 3:
                    srcfyl = srclmp + '/in.equil2'
                else:
                    srcfyl = srclmp + '/in.mc_equil2'

                desfyl = destdir + '/in.equil2'
                shutil.copy2(srcfyl, desfyl)

                if int(graftopt) == 3:
                    srcfyl = srclmp + '/in.prod'
                else:
                    srcfyl = srclmp + '/in.mc_prod'
                desfyl = destdir + '/in.prod'
                shutil.copy2(srcfyl, desfyl)

                restartset = 'restart*'
                restartlist = glob.glob(restartset)
                if restartlist != []:
                    srcfyl = srclmp + '/jobmain_long_prod_var.sh'
                    desfyl = destdir + '/jobmain_long_prod_var.sh'
                    shutil.copy2(srcfyl, desfyl)
                    mainfyl = 'jobmain_long_prod_var.sh'
                    print("found restart files in", destdir)
                else:
                    print("Did not find restart files in", destdir)
                    continue

                launch_fyl = 'jobmain_prod_long.sh'
                fr = open(mainfyl,'r')
                fw = open(launch_fyl,'w')
                jobstr = "joblong_"+ str(nmons[bblen]) + "_" \
                         + str(graftarr[glen])
                fid = fr.read().replace("py_jobname",jobstr)
                fw.write(fid)
                fw.close()
                fr.close()

                srcfyl = '/home/dorfmank/vsethura/mylammps/src/lmp_mesabi'
                desfyl = destdir + '/lmp_mesabi'
                shutil.copy2(srcfyl, desfyl)
                print( "Submitting Jobs..")
                
                subprocess.call(["qsub","jobmain_prod_long.sh"])
                
                print("\n")
                os.chdir(maindir)
                

                
