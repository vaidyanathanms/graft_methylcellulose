# To analyze structural properties for semiflexible chains
# Version_2: V_Mar_26_2019

import numpy
import os
import shutil
import subprocess
import sys
import glob
import datetime

#---------------Import Functions-----------------------------------

from subprocess import call
from my_python_functions import cpy_main_files 
from my_python_functions import find_datafyl
from my_python_functions import find_recent_traj_file
from my_python_functions import write_to_file
from my_python_functions import find_datafyl_mc
from my_python_functions import mult_launch_fyl
from my_python_functions import write_many_files

#------------------Input Options----------------------------------

# 0 - No grafts
# 1 - Make PEG as MC substitute 
# 2 - Grow side chains for MC
# 3 - Grow side chains for semiflexible chains
# 3.1 - Flexible back bone

graftopt = 2
nframes  = 1000
skframes = 0
nsolv    = 0
allana   = 1 #1 - analysis of all dumpfiles, 0 - only latest
prefix_traj = 'dump_stage_*.lammpstrj'
wlccheck = 1 #only valid when graftopt = 3, simple WLC model
anglstyl = 'harmonic' #cosine or harmonic
config  = 2

#-----------------Input Arrays------------------------------------

epsarr_pg   = [0.8]#,1.0,1.2]
epsarr_ps   = [1.0]#,1.0,1.0]
epsarr_sg   = [1.0]#,1.0,1.0]

nchains     = 2 # Number of backbone chains
nmons       = [1000]#,500,1000,2000]#,500,1000,2000] # Number of backbone monomers
graftMW     = 25  # number of graft monomers per graft
polywtperc  = 1.0 # total polymer weight percentage
polydens    = 0.1

#Graft percentage/chain
graftarr    = [0.01]#,0.05,0.1,0.15,0.2,0.25,0.3]
#graftarr    = [0.01,0.03,0.05,0.08,0.12,0.16,0.2,0.24,0.28,0.32] 
#graftarr     = [0.00]

DS_MC   = '1.80' # String Value
temp    = '50.0' # String Value 

#-----------------File List---------------------------------------

if wlccheck == 1 and graftopt == 3:
    inpananame = 'anainp_wlc_var.txt'
else:
    inpananame = 'anainp_mc_var.txt'

print("input analysis file", inpananame)

fyl_list    = ['main.f90','ran_numbers.f90','params.f90',
               inpananame,'jobana.sh','jobana_var.sh']
    
#----------------Directory Settings-------------------------------

maindir     = os.getcwd()
scratchdir  = '/scratch.global/vaidya/'

if not os.path.isdir(scratchdir):
    print(scratchdir, "does not exist")
    sys.exit()

#----------------Outfile Settings---------------------------------

cur_year = datetime.datetime.now().strftime("%y")
cur_mon  = datetime.datetime.now().strftime("%m")
cur_date = datetime.datetime.now().strftime("%d")
outfyl_ana = 'analysis_data_' + str(cur_mon) + '_' + str(nchains) \
             + '_' + str(cur_date) + '_' + str(cur_year) + '.dat'

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
        if anglstyl == 'cosine':
            workdir_main = workdir_main + '/semiflex_grafts_cosstyle' 
        else:
            workdir_main = workdir_main + '/semiflex_grafts' 
    elif graftopt == 3.1:
        workdir_main = workdir_main + '/flex_bb' 
    else:
        print("unknown graft option")
        
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


    workdir_config = workdir_temp + '/config_' + str(config)
    if not os.path.isdir(workdir_config):
        print(workdir_config," does not exist")
        continue


    workdir_chain = workdir_config + '/nchains_' + str(nchains)
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


    if bblen == 0: #Create or append log files
        analys_fyl = workdir_temp+'/'+outfyl_ana
        if os.path.exists(analys_fyl):
            fana_out = open(analys_fyl,'a')
        else:
            fana_out   = open(analys_fyl,'w')
            fana_out.write('%s\t%s\t%s\t%s\t%s\n' % ('mw_bb','mw_graft',\
                                                     'graft_%','eps_pg',\
                                                     'dumpfile_name'))
    for glen in range(len(graftarr)): #Graft loop

        print( "Graft Percentage: ", graftarr[glen])            

        workdir1 = workdir_graft + '/graftperc_' + str(graftarr[glen])
        if not os.path.isdir(workdir1):
            print(workdir1, " does not exist")
            continue

        #Continue iff at least one graft is present

        num_graft_molecules = int(graftarr[glen]*nmons[bblen])
        print("Number of graft molecules: ",num_graft_molecules)
        if num_graft_molecules < 1 and wlccheck == 0:
            print("number of graft less than 1 for sigma:", graftarr[glen])
            continue

        for eps in range(len(epsarr_pg)): #eps_pg loop
            
            if polywtperc != 1.0:
                data_fname = "input_" + str(epsarr_ps[eps]) + "_" + \
                             str(epsarr_sg[eps]) + ".dat"
            else:
                data_fname = "input_" + str(epsarr_pg[eps]) + ".dat"
                
            print( "EpsValue_polymer_graft: ", epsarr_pg[eps])

            workdir3 = workdir1 + '/epsval_' + str(epsarr_pg[eps])

            if not os.path.isdir(workdir3):
                print(workdir3, " does not exist")
                continue

            #------All copying/manipulations-------------------------------------
            os.chdir(workdir3)
            destdir = os.getcwd()

            print( "Copying Files..")
            
            for fylcnt in range(len(fyl_list)):
                cpy_main_files(maindir, destdir, fyl_list[fylcnt])

            print("Creating/finding Datafile..")
            if int(graftopt) == 3:
                [datafyle,flag] = find_datafyl(destdir,
                                               epsarr_ps[eps],epsarr_sg[eps],epsarr_pg[eps],wlccheck)
            else:
                [datafyle,flag] = find_datafyl_mc(destdir,
                                                  epsarr_pg[eps])

            if flag == -1:
                print("datafile not found in ", destdir)
                continue

            if allana == 0: #Only analyzing the latest trajectory

                print("Finding latest trajectory file ..")
                latest_traj_fyle = find_recent_traj_file(destdir,
                                                     prefix_traj)
                if latest_traj_fyle == "nil":
                    print("trajectory file not found in", destdir)
                    continue
                      

                print("Manipulating Files..")
                write_to_file(inpananame,datafyle,latest_traj_fyle,
                              nframes,skframes,nchains,nmons[bblen])
                    
                fana_out.write('%d\t%d\t%f\t%f\t%s\n' % (nmons[bblen],graftMW,\
                                                         graftarr[glen],\
                                                         epsarr_pg[eps],\
                                                         latest_traj_fyle))

                print("Compiling files and submitting..")
                subprocess.call(["ifort","-r8","-qopenmp","-mkl","-check",
                                 "-traceback","params.f90","main.f90",
                                 "-o","ana.o"])

                subprocess.call(["qsub","jobana.sh"])


            else: #analyzing all trajectory files

                trajnames = destdir + '/' + prefix_traj
                list_of_files = glob.glob(trajnames)
                list_cntr = 0

                modfyl = glob.glob('*.o')
                for f in modfyl:
                    os.remove(f)

                for fname in range(len(list_of_files)):

                    list_cntr = list_cntr + 1
                    traj_name_full = list_of_files[fname]
                    traj_name = traj_name_full.split('/')
                    traj_name = traj_name[len(traj_name) - 1]
                    print("Compiling files and submitting for:",traj_name)

                    ana_name  = write_many_files(inpananame,
                                                 datafyle,traj_name,
                                                 nframes,skframes,nchains,
                                                 nmons[bblen],list_cntr)
                    
                        
                    if ana_name == "nil":
                        print("Error: Did not generate analysis input")
                        continue

                    objname = "ana_" + str(list_cntr) + ".o"
                    subprocess.call(["ifort","-r8","-qopenmp","-mkl","-check",
                                     "-traceback","params.f90","main.f90",
                                     "-o",objname])

                    jobfile = "jobana_" + str(list_cntr) + ".sh"
                    mult_launch_fyl(destdir,'jobana_var.sh',ana_name,objname,
                                    graftarr[glen], epsarr_pg[eps],jobfile)

                    subprocess.call(["qsub",jobfile])

                os.remove('jobana_var.sh')

fana_out.close()


            

