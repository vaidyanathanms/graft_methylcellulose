#Select latest files of a given type and remove the coordinate from
#file name. For ease in MATLAB postprocessing
import numpy
import os
import shutil
import subprocess
import sys
import glob
import os.path
from subprocess import call

graftopt = 1 # 1- Make CG as graft, 2 -Grow side chains

graftarr = [0.05, 0.1,0.15,0.2,0.25,0.3]#, 0.5, 0.7]
epsarr   = [0.2,1.0,1.5]  #, 2.5]

nchains = 1
nmons   = 1000
DS_MC   = '1.80' # String Value
temp    = '50.0' # String Value 

maindir = os.getcwd()
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
    
        # * Heart of the Code
        
        resdir = destdir

        if not os.path.exists(resdir):
            print(resdir, "not found")
            sys.exit()

        superdir = workdir1 + '/finconf_all'

        if not os.path.exists(superdir):
            os.mkdir(superdir)


        f_name = resdir + '/restart*'
        list_of_files = glob.glob(f_name) 

        if list_of_files:

            latest_fyl = max(list_of_files, key=os.path.getctime)
            finf_name = superdir + '/finconf_' + str(graftarr[glen]) \
                + '_eps_' + str(epsarr[eps]) + '.txt'

            if os.path.isfile(finf_name):
                print(finf_name)
                os.remove(finf_name)
                subprocess.call(["mpirun","-np","24","./lmp_mesabi",
                                 "-r", latest_fyl, finf_name])
            else:
                subprocess.call(["mpirun","-np","24","./lmp_mesabi",
                                 "-r", latest_fyl, finf_name])


