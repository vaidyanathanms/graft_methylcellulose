# Functions for restarting simulations for genconf3.py
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

#------------------Functions---------------------------------

def start_from_beginning(maindir,srclmp,lmp_execdir,destdir,\
                         data_fname,temp,polywtperc,rcut,blist,\
                         alist,k_phi,wlccheck,anglstyl,nchains,\
                         polydens,graftMW,initcom,nmonval,\
                         graftval,epsval,graftopt,geninp_list,\
                         begin_fresh,key_epsval,epsval_pg):
    
    if begin_fresh == -1:
        list_of_files = glob.glob('*')
        if list_of_files != []:
            for f in list_of_files:
                if not os.path.isdir(f):
                    os.remove(f)
                
    for fylcnt in range(len(geninp_list)):
        cpy_main_files(maindir,destdir,geninp_list[fylcnt])

    # Manipulate Files
    os.chdir(destdir)
    print( "Manipulating input parameter file..")
                
    launch_fyl = 'params_input.dat'
    fr = open('params_input_var.dat','r')
    fw = open(launch_fyl,'w')
    fid = fr.read().replace("py_inputname",str(data_fname)).\
          replace("py_graftstyle",str(int(graftopt))).\
          replace("py_nchains_bb",str(nchains)).\
          replace("py_nmons_bb",str(nmonval)).\
          replace("py_graftperc",str(graftval)).\
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

    datafind = 1
    if not os.path.exists(data_fname):
        print(data_fname, "not found - check input code")
        datafind = -1

    if datafind == 1:
        cpy_lammps_files(maindir,srclmp,lmp_execdir,destdir,graftopt,\
                         data_fname,nmonval,graftval,anglstyl)
        create_infile(maindir,srclmp,destdir,nmonval,epsval,graftopt,\
                      temp,graftval,polywtperc,rcut,blist,alist,k_phi,\
                      wlccheck,anglstyl,graftval,key_epsval,epsval_pg)
        generate_backup(maindir,srclmp,lmp_execdir,destdir)


def run_long_equil_cycle(maindir,destdir,srclmp,nmonval,graftval,graftopt):

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


    #in case the polygraft file has been deleted
    #works only for those cases where Prof. IS suggestion is NOT valid
    srcfyl = srclmp + '/in.cgpair_polygraft_MC'
    desfyl = destdir + '/in.cgpair_polygraft_MC'
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
    jobstr = "joblong_"+ str(nmonval) + "_" + str(graftval)
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


def run_production_cycle(maindir,destdir,srclmp,nmonval,graftval,graftopt):

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

    #in case the polygraft file has been deleted
    #works only for those cases where Prof. IS suggestion is NOT valid
    srcfyl = srclmp + '/in.cgpair_polygraft_MC'
    desfyl = destdir + '/in.cgpair_polygraft_MC'
    shutil.copy2(srcfyl, desfyl)

    restartset = 'restart*'
    restartlist = glob.glob(restartset)
    refind = 1
    if restartlist != []:
        print("found restart files in", destdir)
        srcfyl = srclmp + '/jobmain_long_prod_var.sh'
        desfyl = destdir + '/jobmain_long_prod_var.sh'
        shutil.copy2(srcfyl, desfyl)
        mainfyl = 'jobmain_long_prod_var.sh'
    else:
        print("Did not find restart files in", destdir)
        refind = -1        

    if refind == 1:
        launch_fyl = 'jobmain_prod_long.sh'
        fr = open(mainfyl,'r')
        fw = open(launch_fyl,'w')
        jobstr = "joblong_"+ str(nmonval) + "_" + str(graftval)
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



def run_movie_cycle(maindir,destdir,srclmp,nmonval,graftval,graftopt):

    os.chdir(destdir)
    restartset1 = 'restart1'
    restartlist1 = glob.glob(restartset1)
    restartset2 = 'restart2'
    restartlist2 = glob.glob(restartset2)
    
    refind = 1
    if restartlist1 != []:
        print("found restart1 files in", destdir)
        srcfyl = srclmp + '/in.makemovie1'
        desfyl = destdir + '/in.makemovie'
        shutil.copy2(srcfyl, desfyl)
    elif restartlist2 != []:
        print("found restart2 files in", destdir)
        srcfyl = srclmp + '/in.makemovie2'
        desfyl = destdir + '/in.makemovie'
        shutil.copy2(srcfyl, desfyl)
    else:
        print("Did not find restart files in", destdir)
        refind = -1        

    #in case the polygraft file has been deleted
    #works only for those cases where Prof. IS suggestion is NOT valid
    srcfyl = srclmp + '/in.cgpair_polygraft_MC'
    desfyl = destdir + '/in.cgpair_polygraft_MC'
    shutil.copy2(srcfyl, desfyl)

    if refind == 1:
        srcfyl = srclmp + '/jobmain_movie_var.sh'
        desfyl = destdir + '/jobmain_movie_var.sh'
        shutil.copy2(srcfyl, desfyl)

        mainfyl = 'jobmain_movie_var.sh'
        launch_fyl = 'jobmain_movie.sh'
        fr = open(mainfyl,'r')
        fw = open(launch_fyl,'w')
        jobstr = "joblong_"+ str(nmonval) + "_" + str(graftval)
        fid = fr.read().replace("py_jobname",jobstr)
        fw.write(fid)
        fw.close()
        fr.close()
    
        srcfyl = '/home/dorfmank/vsethura/mylammps/src/lmp_mesabi'
        desfyl = destdir + '/lmp_mesabi'
        shutil.copy2(srcfyl, desfyl)
        print( "Submitting Jobs..")
        
        subprocess.call(["qsub","jobmain_movie.sh"])
        
        print("\n")
        os.chdir(maindir)

def run_autocf(maindir,destdir,srclmp,nmonval,graftval,graftopt):

    os.chdir(destdir)
    restartset1 = 'restart1'
    restartlist1 = glob.glob(restartset1)
    restartset2 = 'restart2'
    restartlist2 = glob.glob(restartset2)
    
    refind = 1
    if restartlist1 != []:
        print("found restart1 files in", destdir)
        srcfyl = srclmp + '/in.runacf1'
        desfyl = destdir + '/in.runacf'
        shutil.copy2(srcfyl, desfyl)
    elif restartlist2 != []:
        print("found restart2 files in", destdir)
        srcfyl = srclmp + '/in.runacf2'
        desfyl = destdir + '/in.runacf'
        shutil.copy2(srcfyl, desfyl)
    else:
        print("Did not find restart files in", destdir)
        refind = -1        

    if refind == 1:
        srcfyl = srclmp + '/jobmain_acf_var.sh'
        desfyl = destdir + '/jobmain_acf_var.sh'
        shutil.copy2(srcfyl, desfyl)

        mainfyl = 'jobmain_acf_var.sh'
        launch_fyl = 'jobmain_acf.sh'
        fr = open(mainfyl,'r')
        fw = open(launch_fyl,'w')
        jobstr = "joblong_"+ str(nmonval) + "_" + str(graftval)
        fid = fr.read().replace("py_jobname",jobstr)
        fw.write(fid)
        fw.close()
        fr.close()
    
        srcfyl = '/home/dorfmank/vsethura/mylammps/src/lmp_mesabi'
        desfyl = destdir + '/lmp_mesabi'
        shutil.copy2(srcfyl, desfyl)
        print( "Submitting Jobs..")
        
        subprocess.call(["qsub","jobmain_acf.sh"])
        
        print("\n")
        os.chdir(maindir)
