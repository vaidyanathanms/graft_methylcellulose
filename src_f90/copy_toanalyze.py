# To copy the resultant files for further analysis

import numpy
import os
import shutil
import subprocess
import sys
import glob

#---------------Import Functions-----------------------------------

from subprocess import call
from my_python_functions import cpy_main_files 
from my_python_functions import find_datafyl
from my_python_functions import find_recent_traj_file
from my_python_functions import write_to_file
from my_python_functions import my_cpy_generic
from my_python_functions import find_recent_file

#------------------Input Options----------------------------------

# 0 - No grafts
# 1 - Make PEG as MC substitute 
# 2 - Grow side chains for MC
# 3 - Grow side chains for semiflexible chains
# 3.1 - Flexible back bone

graftopt = 2
wlccheck = 1 #valid only for graftopt = 3
anglstyl = 'harmonic' #cosine or harmonic
config   = 1 #configuration to copy
figflag  = 0 #for copying figure_files
fylflag  = 1 #for copying analysis_files

#-----------------Input Arrays------------------------------------

#epsarr_gg   = [0.8,1.0,1.2]
#epsarr_ps   = [1.0,1.0,1.0]
#epsarr_sg   = [1.0,1.0,1.0]

epsarr_gg   = [1.0]
epsarr_ps   = [1.0]
epsarr_sg   = [1.0]

nchains     = 2 # Number of backbone chains
nmons       = [1000] # Number of backbone monomers
graftMW     = 25  # number of graft monomers per graft
polywtperc  = 1.0 # total polymer weight percentage
polydens    = 0.1 # density of the system

#Graft percentage/chain
#graftarr    = [0.0]
graftarr    = [0.15]#,0.01,0.05,0.1,0.15,0.2,0.25,0.3]#,0.24,0.28,0.32] 
#graftarr    = [0.20,0.25,0.30]
# polydens = totpart/box_volume - explict generic
# polydens = wt percentage "MC" - implicit MC

DS_MC   = '1.80' # String Value
temp    = '50.0' # String Value 

#-----------------File List---------------------------------------

#Give prefix for files to be copied followed by *
fyl_list    = ['avgeival_*','compos_*','eigMCavg_*','indeig_*',\
               'mainpersistautocf_*','rgall_*','rgavgall_*',\
               'rgMConly_*','rgsys_*','distcom_*','job*']
figname     = 'dump_stage_*'

#----------------Directory Settings-------------------------------

maindir     = os.getcwd()
scratchdir  = '/scratch.global/vaidya/'

if not os.path.isdir(scratchdir):
    print(scratchdir, "does not exist")
    sys.exit()

#--------------Main Analysis------------------------------------

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
        
    if not os.path.isdir(workdir_main):
        print(workdir_main," does not exist")
        continue

    #Backbone directory
    workdir_bb_main = workdir_main + '/n_bb_' + str(nmons[bblen])
    if not os.path.isdir(workdir_bb_main):
        print(workdir_bb_main," does not exist (backbone dir)")

    #Temperature directory
    if int(graftopt) != 3:
        workdir_temp = workdir_bb_main + '/temp_' + temp
        if not os.path.isdir(workdir_temp):
            print(workdir_temp," does not exist (temp dir)")
            continue
    else:
        workdir_temp = workdir_bb_main

    #Config directory
    workdir_config = workdir_temp + '/config_' + str(config)
    if not os.path.isdir(workdir_config):
        print(workdir_config, "not found (chains dir)")
        continue
    

    #Chains directory
    workdir_chain = workdir_config + '/nchains_' + str(nchains)
    if not os.path.isdir(workdir_chain):
        print(workdir_chain, "not found (chains dir)")
        continue

    #wt% directory        
    workdir_bb  = workdir_chain + '/backbonewtperc_'+str(polydens)
    if not os.path.isdir(workdir_bb):
        print(workdir_bb, "not found (wt% dir)")
        continue
    
    #Grafts MW directory
    workdir_graft = workdir_bb + '/n_graft_' + str(graftMW)
    if not os.path.isdir(workdir_graft):
        print(workdir_graft, " does not exist (graftMW dir)")
        continue
    
    #Make analysis files output directory

    if fylflag == 1:
        out_dir = 'out_'+'bbMW_' + str(nmons[bblen]) + \
                  '_ngMW_' + str(graftMW) + '_rho_'+\
                  str(polydens) + '_nch_' + str(nchains)
        anafyl_dir = workdir_graft + '/' + out_dir 
        if not os.path.isdir(anafyl_dir):
            os.mkdir(anafyl_dir)


    for glen in range(len(graftarr)): #Graft loop

        print( "Graft Percentage: ", graftarr[glen])

        workdir1 = workdir_graft + '/graftperc_' + str(graftarr[glen])

        #Continue iff at least one graft is present
        num_graft_molecules = int(graftarr[glen]*nmons[bblen])
        print("Number of graft molecules: ",num_graft_molecules)
        if num_graft_molecules < 1 and wlccheck == 0:
            print("number of graft less than 1 for sigma:", graftarr[glen])
            continue

        if not os.path.isdir(workdir1):
            print(workdir1, " does not exist (grafts per bb)")
            continue

        for eps in range(len(epsarr_gg)): #eps_pg loop
            
            workdir3 = workdir1 + '/epsval_' + str(epsarr_gg[eps])

            if not os.path.isdir(workdir3):
                print(workdir3, " does not exist")
                continue

            #Make output directory for figure files if needed
            if figflag == 1:
                outfig_name = 'snapshots_'+'config_' + str(config)\
                              + '_bbMW_' + str(nmons[bblen]) + \
                              '_ngMW_' + str(graftMW) + '_rho_'+\
                              str(polydens) + '_nch_' + str(nchains) \
                              + '_sig_' + str(graftarr[glen]) + \
                              '_eps_' + str(epsarr_gg[eps])
 
                outfig_dir  = workdir_temp + '/allsnapshots'
                if not os.path.isdir(outfig_dir):
                    os.mkdir(outfig_dir)
                

            #------All copying/manipulations--------------------------
            if fylflag == 1:

                workdir4 = workdir3 + '/all_output_data'
                if not os.path.isdir(workdir4):
                    print(workdir4, " does not exist")
                    continue

                print("Copying output files from", workdir4)

                for fylcnt in range(len(fyl_list)):

                    os.chdir(workdir4)
                    destdir = os.getcwd()
                    list_of_files = glob.glob(fyl_list[fylcnt])

                    #search in all_output_data
                    if list_of_files == []:

                        os.chdir(workdir3)
                        destdir = workdir3
                        list_of_files = glob.glob(fyl_list[fylcnt])

                        #search in previous directory
                        if list_of_files == []:
                            print("Did not find files of type", \
                                  fyl_list[fylcnt])
                            continue

                    print( "Copying files of type: ", fyl_list[fylcnt])

                    for filenum in range(len(list_of_files)):

                        fylname = list_of_files[filenum]
                        if fyl_list[fylcnt] == 'job*':
                            anafylname = fylname
                        else:
                            outfyl  = fylname.split("_")[0]
                            anafylname = outfyl + '_' + str(epsarr_gg[eps]) + '_' \
                                         + str(graftarr[glen]) + '_' + \
                                         str(nmons[bblen])  + '_' + \
                                         str(filenum+1) + '.dat'

                        my_cpy_generic(destdir,anafyl_dir,fylname,anafylname)


            #Copy all figures if needed
            if figflag == 1:

                specfig_dir = outfig_dir + '/' + outfig_name
                if not os.path.isdir(specfig_dir):
                    os.mkdir(specfig_dir)

                os.chdir(workdir3)
                
                #check if figure directory is present
                fig_maindir = workdir3 + '/fig_files'
                if not os.path.isdir(fig_maindir):
                    print(fig_maindir, " does not exist")
                    continue
                else:
                    os.chdir(fig_maindir)

                list_of_files = glob.glob(figname)

                if list_of_files == []:
                    print("No figure files with ", figname, "in ", fig_maindir) 

                for filenum in range(len(list_of_files)):

                    fylname = list_of_files[filenum]
                    my_cpy_generic(fig_maindir,outfig_dir,fylname,fylname)
                
