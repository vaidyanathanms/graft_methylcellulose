# To perform SMD/US simulations
# Version: V2_May_14_2019
# Added umbrella sampling capabilities

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
from my_python_functions import create_datafyle
from my_python_functions import gen_extract_inp
from my_python_functions import gen_launch_fyl
from my_python_functions import cleanup_files
from my_python_functions import create_umbcolfile
from my_python_functions import remove_chunks

#------------------Input Options---------------------------------

# 0 - No grafts
# 1 - Make PEG as MC substitute (implicit solv)
# 2 - Grow side chains for MC (implicit solv)
# 3 - Grow side chains for semiflexible chains (explicit solv)
# 3.1 - Fredrickson's case (explicit/implicit solv)

graftopt = 2 # see options above
ana_flag = 'umbrella' #options -> smd, umbrella

#-----------------Input Arrays------------------------------------

epsarr_pg   = [0.8]#,1.0,1.2]  # polymer-graft
epsarr_ps   = [1.0]#,1.0,1.0]  # polymer-solvent
epsarr_sg   = [1.0]#,1.0,1.0]  # solvent-graft
rcut        = pow(2,1/6)

nchains     = 2 # Number of backbone chains
nmons       = [1000] # Number of backbone monomers

graftMW     = 25  # number of graft monomers per graft
polywtperc  = 1.0 # tot_poly wt% - explicit generic 
polydens    = 0.5 # overall density of polymers

#Graft percentage/chain
graftarr    = [0.03,0.08,0.16,0.20] 

DS_MC   = '1.80' # String Value - only for methylcellulose
temp    = '50.0' # String Value - for computing parameters for MC

#----------------File List----------------------------------------

#for SMD
#in.cgpair_MC will be recreated in the SMD directory
fyls_main = ['in.smd','colfile.inp','jobsmd_var.sh']
smd_fylename = 'smdinp.data'
lmp_fylename = 'lmp_mesabi'
res_prefix   = 'archive_restart_prod.*'

#for umbrella sampling
umblmp_fyls  = ['in.umbrella','umbcolfile_var','jobumb_var.sh']
umbmain_fyls = ['extract_conf.f90','extract_params.f90',
                'extract_inp_var.txt'] #This should go to SMD dir
US_dirprefix = 'umbmain'
umb_centers  = [4.5,5.5,6.5,8.0,10.0,11.0,12.0,14.0,16.0,18.0,\
                19.0,20.0,21.0,22.0,23.0]
smdtraj_pref = 'smdtraj_*'
tolval       = 0.3
axval        = 0 # no preference
comflag      = 1
forcecon     = 20.0
pbcval       = 1 # pbc ON

#----------------Directory Settings-------------------------------

srclmpdir   = '/home/dorfmank/vsethura/allfiles/files_graft/src_lmp'
lmp_execdir = '/home/dorfmank/vsethura/mylammps/src'
maindir     = os.getcwd()
scratchdir  = '/scratch.global/vaidya/'

if not os.path.isdir(scratchdir):
    print(scratchdir, "does not exist")
    sys.exit()

#--------------Main Analysis--------------------------------------

for bblen in range(len(nmons)): #Backbone length loop

    print( "Backbone length: ", nmons[bblen])
    workdir_main = scratchdir + 'grafts_all' 
    if not os.path.isdir(workdir_main):
        print(workdir_main," does not exist")
        sys.exit()

    if graftopt == 0:
        workdir_main = workdir_main + '/no_grafts' 
    elif graftopt == 1:
        workdir_main = workdir_main + '/substitute_grafts_MC' 
    elif graftopt == 2:
        workdir_main = workdir_main + '/grafts_MC' 
    elif graftopt == 3:
        workdir_main = workdir_main + '/semiflex_grafts' 
    elif graftopt == 3.1:
        workdir_main = workdir_main + '/flex_bb' 
        
    if not os.path.isdir(workdir_main):
        print(workdir_main," does not exist")
        continue
    
    workdir_bb_main = workdir_main + '/n_bb_' + str(nmons[bblen])
    if not os.path.isdir(workdir_bb_main):
        print(workdir_bb_main," does not exist")
        continue

    if int(graftopt) != 3:
        workdir_temp = workdir_bb_main + '/temp_' + temp
        if not os.path.isdir(workdir_temp):
            print(workdir_temp," does not exist")
            continue
    else:
        workdir_temp = workdir_bb_main


    workdir_chain = workdir_temp + '/nchains_' + str(nchains)
    if not os.path.isdir(workdir_chain):
        print(workdir_chain, "not found")
        continue
        
    workdir_bb  = workdir_chain + '/backbonewtperc_' + \
                  str(polydens)
    
    if not os.path.isdir(workdir_bb):
        print(workdir_bb, "not found")
        continue

    workdir_graft = workdir_bb + '/n_graft_' + str(graftMW)
        
    if not os.path.isdir(workdir_graft):
        print(workdir_graft, " does not exist")
        continue


    for glen in range(len(graftarr)): #Graft loop

        print( "Graft Percentage: ", graftarr[glen])            

        workdir1 = workdir_graft + '/graftperc_' + str(graftarr[glen])
        if not os.path.isdir(workdir1):
            print(workdir1, " does not exist")
            continue

        #Continue iff at least one graft is present

        num_graft_molecules = int(graftarr[glen]*nmons[bblen])
        print("Number of graft molecules: ",num_graft_molecules)
        if num_graft_molecules < 1:
            print("number of graft less than 1 for sigma:", graftarr[glen])
            continue

        for eps in range(len(epsarr_pg)): #eps_pg loop
            
            workdir3 = workdir1 + '/epsval_' + str(epsarr_pg[eps])

            if not os.path.isdir(workdir3):
                print(workdir3, " does not exist")
                continue

            if polywtperc != 1.0:
                data_fname = "input_" + str(epsarr_ps[eps]) + "_" + \
                             str(epsarr_sg[eps]) + ".dat"
            else:
                data_fname = "input_" + str(epsarr_pg[eps]) + ".dat"
                
            print( "EpsValue_polymer_graft: ", epsarr_pg[eps])


            if ana_flag == 'smd':
                
                print ("Generating SMD initial conditions ..")
                #------Create SMD work directory-----------------------

                traj_workdir = workdir3
                smd_workdir  = traj_workdir + '/smd_main'
                if not os.path.isdir(smd_workdir):
                    os.mkdir(smd_workdir)

                #------All copying/manipulations-----------------------

                print("Finding latest archived restart file..")
                refyle = find_recent_file(traj_workdir,res_prefix)

                if refyle == "nil":
                    print("no restart file found in ", traj_workdir)
                    continue
                else:
                    os.chdir(smd_workdir)
                    print("Creating datafile..")
                    cpy_main_files(traj_workdir,smd_workdir,refyle)
                    lmp_execfile = lmp_execdir + '/' + lmp_fylename
                    create_datafyle(lmp_execfile,refyle,smd_fylename)

                for fylcnt in range(len(fyls_main)):
                    cpy_main_files(srclmpdir, smd_workdir,fyls_main[fylcnt])

                cpy_main_files(lmp_execdir,smd_workdir,lmp_fylename)

                create_infile(maindir,srclmpdir,smd_workdir,nmons[bblen],
                              epsarr_pg[eps],graftopt,temp)



                gen_launch_fyl(smd_workdir,'jobsmd_var.sh','jobsmd.sh',
                               graftarr[glen],epsarr_pg[eps])
                subprocess.call(['qsub','jobsmd.sh'])
                cleanup_files(umbdir,"f90")
                cleanup_files(umbdir,"txt")
                os.chdir(maindir)

            elif ana_flag == 'umbrella':

                print ("Generating umbrella sampling init config ..")
                
                #-----Check for SMD directory ----------------------#
                
                traj_workdir = workdir3
                smd_workdir  = traj_workdir + '/smd_main'
                if not os.path.isdir(smd_workdir):
                    print(smd_workdir, "not found")
                    continue

                ref_data = smd_workdir + "/" + smd_fylename
            
                if not os.path.exists(ref_data):
                    print("SMD reference file not found: ", smd_workdir)
                    continue

                for ucen in range(len(umb_centers)):
                
                    umbdir = traj_workdir + '/' + US_dirprefix + \
                             '_' + str(umb_centers[ucen])
            
                    print("US directory under construction", umbdir)
                    if not os.path.isdir(umbdir):
                        os.mkdir(umbdir)
                    else:
                        remove_chunks(umbdir, '*')

                    os.chdir(umbdir)
            
                    #umbmain_fyls should go into SMD directory
                    for fcnt in range(len(umbmain_fyls)):
                        cpy_main_files(maindir,smd_workdir,umbmain_fyls[fcnt])

                    for f2cnt in range(len(umblmp_fyls)):
                        cpy_main_files(srclmpdir,umbdir,umblmp_fyls[f2cnt])

                    cpy_main_files(lmp_execdir,umbdir,lmp_fylename)
                    
                    #--Generating input files-----------------------------------

                    fextout = gen_extract_inp('extract_inp_var.txt',\
                                              smdtraj_pref,smd_workdir,
                                              umb_centers[ucen],\
                                              tolval,axval,comflag,pbcval)

                    if fextout == "nil":
                        print("No trajectory files exist for ", umbcenters[ucen])


                    subprocess.call(["ifort","-mkl","-qopenmp","-r8",\
                                     "extract_params.f90","extract_conf.f90",\
                                     "-o","extract.o"])
                    subprocess.call(["./extract.o",fextout])


                    checkinitdata = smd_workdir + '/init_datafile'
                    desfyl = umbdir + '/init_datafile'
                
                    if os.path.exists(checkinitdata):

                        print(checkinitdata, "prepared for US simulations.")
                        shutil.copy2(checkinitdata,desfyl)
                    
                        #Delete & keep a copy under different name for
                        #future copying if necessary

                        prev_US_fyl = smd_workdir + '/previnit_datafile'
                        shutil.copy2(checkinitdata,prev_US_fyl)
                        os.remove(checkinitdata)
                        
                    else: #if US could not find a file near that place -
                        #use SMD initialized file or latest US file
                    
                        print(checkinitdata, "not prepared for",
                              umb_centers[ucen])

                        prev_US_fyl = smd_workdir + '/previnit_datafile'

                        if os.path.exists(prev_US_fyl):

                            print("Copying previous US file")
                            shutil.copy2(prev_US_fyl,desfyl)

                        else:
                            
                            print("Copying SMD file")
                            ref_data = smd_workdir + '/' + smd_fylename
                            shutil.copy2(ref_data,desfyl)
                        
                    create_infile(maindir,srclmpdir,umbdir,nmons[bblen],
                                  epsarr_pg[eps],graftopt,temp)

                    gen_launch_fyl(umbdir,'jobumb_var.sh','jobumb.sh',
                                   graftarr[glen],epsarr_pg[eps])

                    create_umbcolfile(umbdir,'umbcolfile_var',
                                      umb_centers[ucen],forcecon)
                    print("Submitting US jobs for",umb_centers[ucen])
                    subprocess.call(["qsub", "jobumb.sh"])

                    cleanup_files(umbdir,"f90")
                    cleanup_files(umbdir,"txt")
                    os.chdir(maindir)

            else:

                print("Command not recongized", ana_flag)
                sys.exit()
                
