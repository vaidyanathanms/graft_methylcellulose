# WHAM Analysis (Recalculate with error bars)
# Use with either COLVARS output or from previously stored output traj 

import numpy
import os
import shutil
import subprocess
import sys
import glob
import os.path
from subprocess import call
import re
import random

#--------------------------------Inputs-----------------------------
skipfr  = 10000
kspr = 20.0
correl_time = 20
US_dirprefix = 'umbmain'
umb_centers  = ['3.0','4.0','5.0','6.0','7.0','9.0','11.0','13.0',\
                '15.0','18.0','21.0','24.0'] #must be str array

#---------------------------Path Details---------------------------

bb_MW    = 1000
temp     = '50.0'
nchains  = 2
bb_wtper = '0.5'
gr_MW    = 25
graftarr = [0.03,0.05,0.08,0.12,0.16,0.2] 

#----------------------------WHAM Analysis--------------------------

perform_wham = 1 #1-perform, 0-not perform
minhist = 3.0
maxhist = 21.0
maxbins = 100
maxtol  = 0.01
temperature = 1.0
padding = 0
outfile = "whamout.txt"
maxiter = 1000
stripterm = 'new' #new/none. If new - change before the harmonic cntr

#------------------------------Directories--------------------------

f90src   = '/home/dorfmank/vsethura/allfiles/files_graft/src_f90'
scratchdir = '/scratch.global/vaidya/'
whamexec = '/home/dorfmank/vsethura/wham/wham/wham'

if frombackup == 0:
    maindir = os.getcwd()
else:
    maindir = f90src

#-----------------------------Files-----------------------------------------

whampref = 'whaminp_' #Keep * here
whammain = 'whaminp.txt'
colfile = 'out.colvars.traj'
colfile_name = 'umbsampling'

if frombackup == 0:
    # Keep the prefix without number
    if stripterm:
        prefix  = 'US_dir_' + stripterm 
    else:
        prefix  = 'US_dir_new' 
else:
    prefix = 'USbackup_colvar'

whamcopy = maindir + '/wham'

if not os.path.isfile(whamcopy):
    shutil.copy2(whamexec,whamcopy)

#---------------------------Main Loop Starts--------------------------------

curr_dir = os.getcwd()

print("Starting WHAM Analysis")        


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


            wham_dir = workdir3
            whamcopy = wham_dir + '/wham' #WHAM executable
            whamhead = whammain #files containing WHAM paths
            whaminp  = whampref + '*' #individual WHAM files
            whamout  = outfile #output WHAM file

            rvalarr = []
            if not os.path.isfile(whamcopy):
                shutil.copy2(whamexec,whamcopy)
            else:
                os.remove(whamcopy)  #remove and copy again
                shutil.copy2(whamexec,whamcopy)

            for rvals in range(len(umb_centers)):

                cval = umb_centers[rvals]
                wham_dir2 = workdir3 + '/' + US_dirprefix + '_' + \
                            umb_centers[rvals]

                if not os.path.isdir(wham_dir2):
                    print(wham_dir2, "not found")
                
                os.chdir(wham_dir2)


                list_of_files = glob.glob(whaminp) 


        for fylval in range(len(list_of_files)):
            
            
            splice  = list_of_files[fylval].split('_')
            rval = splice[len(splice)-1]


        fmain = open(whamhead,'w')
        for fylval in sorted(rvalarr):

            if stripterm == 'none':
                fylname = wham_dir2 + '/' + whampref + str(fylval)
            else:
                fylname = wham_dir2 + '/' + whamnew_pref + str(fylval)


            if os.path.isfile(fylname):
                
                fmain.write(fylname + '\t' +str(fylval) +'\t'+ \
                                str(kspr)+'\t' + str(correl_time)+'\n')

        fmain.close()
        print(os.getcwd())
        if perform_wham:

            ranseed = int(10000*random.random())   
            fsub = open('out.txt','w')
            subprocess.call(['./wham',str(minhist),str(maxhist),str(maxbins)\
                             ,str(maxtol),str(temperature),str(padding),\
                             str(whammain),str(whamout),str(maxiter),str(ranseed)]\
                            ,stdout = fsub)

            fsub.close()


