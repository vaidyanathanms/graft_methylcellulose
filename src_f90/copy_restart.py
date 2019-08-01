#Code to copy all restart or any type of files to a common directory
import numpy
import os
import shutil
import subprocess
import sys
import glob

from subprocess import call

graftopt = 2 # 1- Make CG as graft, 2 -Grow side chains

graftarr = [0.05, 0.1,0.15,0.2,0.25, 0.3, 0.5, 0.7]
epsarr   = [0.2,1.0,1.5]  #, 2.5]

nchains = 1
nmons   = 1000
DS_MC   = '1.80' # String Value
temp    = '50.0' # String Value 

maindir = os.getcwd()

if graftopt == 1:
    archivedir = maindir + '/' + 'archive_restart_replaceCG'
else:
    archivedir = maindir + '/' + 'archive_restart_growside'

if not os.path.isdir(archivedir):
    print("Generating Archive Directory")
    os.mkdir(archivedir)

scratchdir = '/scratch.global/vaidya/'
   
for glen in range(len(graftarr)):

    print( "Graft Percentage: ", graftarr[glen])

    workdir1 = scratchdir + 'cg_grafts'

    if graftopt == 1:

        workdir1 = workdir1 + '/switchmon'
    
    elif graftopt == 2:

        workdir1 = workdir1 + '/growsidechains'

    if not os.path.isdir(workdir1):
        print(workdir1, "not found")
        os.mkdir(workdir1)

    workdir2 = workdir1 + '/graft_' + str(graftarr[glen])

    if not os.path.isdir(workdir2):
        print(workdir2, "not found")
        os.mkdir(workdir2)

    for eps in range(len(epsarr)):

        print( "EpsValue: ", epsarr[eps])
        workdir3 = workdir2 + '/epsval_' + str(epsarr[eps])

        if not os.path.isdir(workdir3):
            print(workdir3, "not found")
            os.mkdir(workdir3)

        os.chdir(workdir3)
        destdir = os.getcwd()
    
        print("Copying Files")
            
        # * Heart of the Code

        restart_names = destdir + '/archival*'
        list_of_files = glob.glob(restart_names) 
        
        print("Dataname for LAMMPS: ", restart_names)
        print("Generating Restart Dir")

        restart_dir = destdir + '/restartfiles_sig_' + \
            str(graftarr[glen]) + '_eps_' + str(epsarr[eps])

        if not os.path.exists(restart_dir):
                
            os.mkdir(restart_dir)


        if list_of_files:
	
            for fyl in list_of_files:
                
                subprocess.call(["mv", fyl, restart_dir])
                
        else:
		
            print(restart_names, "not found in", destdir)

        if restart_dir:
                
            print("Zipping Directories ..")
            tar_restart = restart_dir + ".tar.gz"
            print(tar_restart)
            subprocess.call(["tar", "-cvzf", tar_restart, restart_dir])
            desttar = archivedir + '/'
            subprocess.call(["mv", tar_restart, desttar])


