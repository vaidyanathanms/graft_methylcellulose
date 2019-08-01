# Generic Function Definitions
# Version_1: V_Mar_26_2019

import numpy
import os
import shutil
import subprocess
import sys
import glob
import re

def cpy_main_files(maindir,destdir,fylname):

    srcfyl = maindir + '/' + fylname
    desfyl = destdir + '/' + fylname
    shutil.copy2(srcfyl, desfyl)

def find_datafyl(destdir,eps_ps,eps_sg,eps_pg,wlccheck):

    srclmp      = '/home/dorfmank/vsethura/allfiles/files_graft/src_lmp'
    lmp_execdir = '/home/dorfmank/vsethura/mylammps/src'

    if wlccheck == 0:
        dataname = "input_" + str(eps_ps) + "_" + str(eps_sg) + ".dat"
        datapath = destdir + '/' + dataname
    else:
        dataname = "input_" + str(eps_pg) + ".dat"
        datapath = destdir + '/' + dataname
        

    restartfiles = destdir + '/archive_restart*'
    outflag = 1

    if os.path.exists(datapath):
        
        print( "Dataname for LAMMPS found: ", dataname)
        outstr = dataname

    elif glob.glob(restartfiles):
        
        lmpsrcfyl = lmp_execdir + '/lmp_mesabi'
        if not os.path.exists(lmpsrcfyl):
            print("Could not determine path to LAMMPS src\n")
            outstr = "Could not determine path to LAMMPS src\n"
            outflag = -1
        else:
            print( "Creating datafile from restart files: ")
            restartlist = glob.glob(restartfiles)
            latest_fyl = max(restartlist,key=os.path.getctime)
            print(latest_fyl)
            subprocess.call(["./lmp_mesabi","-r",latest_fyl,dataname])
            outstr  = dataname
            print( "Successfully created datafile ")
    else:
            
        print( dataname, "not found in ", datapath)
        outstr  = dataname + "not found in" +  datapath
        outflag = -1
        
    return dataname, outflag


def find_recent_traj_file(destdir,trajkeyword):

    trajnames = destdir + '/' + trajkeyword
    list_of_files = glob.glob(trajnames)
    if list_of_files != []:
        traj_str = max(list_of_files, key=os.path.getctime)
        traj_arr = traj_str.split("/")
        print( "trajval", traj_arr[len(traj_arr)-1])
        return traj_arr[len(traj_arr)-1]
    else:
        return "nil"

def write_to_file(inpfyle_name,dataname,trajname,nframes,skframes
                  ,nchains,nmons):

    logout = "log_" + trajname
    launch_fyl = 'anainp.txt'
    fr  = open(inpfyle_name,'r')
    fw  = open(launch_fyl,'w')
    fid = fr.read().replace("py_datafyl",dataname).\
          replace("py_trajfyl",trajname).\
          replace("py_numframes", str(nframes)).\
          replace("py_skipframes", str(skframes)).\
          replace("py_ntotchains",str(nchains)).\
          replace("py_backbonemons",str(nmons)).\
          replace("py_logname",logout)
    fw.write(fid)
    fw.close()
    fr.close()

def write_many_files(inpfyle_name,dataname,trajname,nframes,skframes
                     ,nchains,nmons,cntrnum):

    logout = "log_" + trajname
    launch_fyl = 'anainp_' + str(cntrnum) + '.txt'
    fr  = open(inpfyle_name,'r')
    fw  = open(launch_fyl,'w')
    fid = fr.read().replace("py_datafyl",dataname).\
          replace("py_trajfyl",trajname).\
          replace("py_numframes", str(nframes)).\
          replace("py_skipframes", str(skframes)).\
          replace("py_ntotchains",str(nchains)).\
          replace("py_backbonemons",str(nmons)).\
          replace("py_logname",logout)
    fw.write(fid)
    fw.close()
    fr.close()

    if os.stat(launch_fyl).st_size == 0:
        return "nil"
    else:
        return launch_fyl

def my_cpy_generic(srcdir,destdir,inpfylname,destfylname):
    src_fyl  = srcdir  + '/' + inpfylname
    dest_fyl = destdir + '/' + destfylname
    shutil.copy2(src_fyl, dest_fyl)

def find_recent_file(destdir,keyword): #A replica of find_recent_traj_file

    fylnames = destdir + '/' + keyword
    list_of_files = glob.glob(fylnames)
    if list_of_files != []:
        fyl_str = max(list_of_files, key=os.path.getctime)
        fyl_arr = fyl_str.split("/")
        print( "File Name: ", fyl_arr[len(fyl_arr)-1])
        return fyl_arr[len(fyl_arr)-1]
    else:
        return "nil"

def find_datafyl_mc(destdir,eps_pg):

    srclmp      = '/home/dorfmank/vsethura/allfiles/files_graft/src_lmp'
    lmp_execdir = '/home/dorfmank/vsethura/mylammps/src'

    dataname = "input_" + str(eps_pg) + ".dat"
    datapath = destdir + '/' + dataname
    
    restartfiles = destdir + '/archive_restart*'
    outflag = 1

    if os.path.exists(datapath):
        
        print( "Dataname for LAMMPS found: ", dataname)
        outstr = dataname

    elif glob.glob(restartfiles):
        
        srcfyl = lmp_execdir + '/lmp_mesabi'

        if not os.path.exists(srcfyl):
            print("Could not determine path to LAMMPS src\n")
            outstr = "Could not determine path to LAMMPS src\n"
            outflag = -1
        else:
            print( "Creating datafile from restart files: ")
            latest_fyl = max(restartfiles,key=os.path.getctime)
            subprocess.call(["./lmp_mesabi","-r",latest_fyl,dataname])
            outstr  = dataname
    
    else:
            
        print( dataname, "not found in ", datapath)
        outstr  = dataname + "not found in" +  datapath
        outflag = -1
        
    return dataname, outflag


def cpy_lammps_files(maindir,srclmp,lmp_execdir,destdir,graftopt,
                     data_fname,nmonval,graftval):

    os.chdir(destdir)
    print( "Copying LAMMPS files..")                
    
    srcfyl = srclmp + '/jobmain_var.sh'
    desfyl = destdir + '/jobmain_var.sh'
    shutil.copy2(srcfyl, desfyl)

    if int(graftopt) == 3:
    
        srcfyl = srclmp + '/in.prod'
        desfyl = destdir + '/in.prod'
        shutil.copy2(srcfyl, desfyl)

        if nmonval < 1000:
            srcfyl = srclmp + '/in.equil_var'
            desfyl = destdir + '/in.equil_var'
        else:
            srcfyl = srclmp + '/in.equil_gt_1000_var'
            desfyl = destdir + '/in.equil_var'
    else:

        srcfyl = srclmp + '/in.mc_prod'
        desfyl = destdir + '/in.prod'
        shutil.copy2(srcfyl, desfyl)

        srcfyl = srclmp + '/in.mc_equil_var'
        desfyl = destdir + '/in.mc_equil_var'
                
                    
    shutil.copy2(srcfyl, desfyl)

    srcfyl = lmp_execdir + '/lmp_mesabi'
    desfyl = destdir + '/lmp_mesabi'
    shutil.copy2(srcfyl, desfyl)
        

    print("Copying in.equil and jobmain.sh files..")

    equil_fyl = 'in.equil'

    if int(graftopt) == 3:
        fr = open('in.equil_var','r')
    else:
        fr = open('in.mc_equil_var','r')
    
    fw = open(equil_fyl,'w')
    fid = fr.read().replace("py_data",data_fname)
    fw.write(fid)
    fw.close()
    fr.close()
                
    launch_fyl = 'jobmain.sh'
    fr = open('jobmain_var.sh','r')
    fw = open(launch_fyl,'w')
    jobstr = "job_"+ str(nmonval) + "_" + str(graftval)
    fid = fr.read().replace("py_jobname",jobstr)
    fw.write(fid)
    fw.close()
    fr.close()

def create_infile(maindir,srclmp,destdir,nmonval,epsval,
                  graftopt,temp,graft_perc,polywtperc,rcut,
                  blist,alist,k_phi,wlcflag):

    os.chdir(destdir)
        
    if int(graftopt) == 1 or int(graftopt) == 2:
                    
        print("Copying files to create the pair coefficients..")
        srcfyl = maindir + '/create_infile_var.f90'
        desfyl = destdir + '/create_infile_var.f90'
        shutil.copy2(srcfyl, desfyl)
        
        srcfyl = maindir + '/cgparams.txt'
        desfyl = destdir + '/cgparams.txt'
        shutil.copy2(srcfyl, desfyl)
        
        launch_fyl = 'create_infile.f90'
        fr = open('create_infile_var.f90','r')
        fw = open(launch_fyl,'w')
        fid = fr.read().replace("py_nmons",str(nmonval)).\
              replace("py_epsval",str(epsval)).\
              replace("py_graftopt",str(int(graftopt)))
        fw.write(fid)
        fw.close()
        fr.close()
        
        print( "Creating PairCoeff File ..")
        
        subprocess.call(["ifort","-r8","-qopenmp","create_infile.f90",
                         "-o","infile.o"])
        subprocess.call(["./infile.o", temp])
        

    elif int(graftopt) == 3:
        
        
        pair_fyl = 'in.pair'
        fw = open(pair_fyl,'w')
        fw.write("#-----Interaction Coeffiecients-- \n")
        fw.write("#PairCoeff 1-Backbone, 2-Graft, 3-Solvent\n")
        
        fw.write('%s\t %d\t %d\t %g\t %g\t %g\n'\
                 %("pair_coeff",1,1,epsval,1.0,rcut))

        if graft_perc != 0.0: #graft present
            fw.write('%s\t %d\t %d\t %g\t %g\t %g\n'\
                     %("pair_coeff",2,2,1.0,1.0,rcut))
            fw.write('%s\t %d\t %d\t %g\t %g\t %g\n'\
                     %("pair_coeff",1,2,epsarr_pg[eps],1.0,rcut))

        if polywtperc != 1.0 and graft_perc != 0.0: #solv,graft
            fw.write('%s\t %d\t %d\t %g\t %g\t %g\n'\
                     %("pair_coeff",3,3,1.0,1.0,rcut))
            fw.write('%s\t %d\t %d\t %g\t %g\t %g\n'\
                     %("pair_coeff",1,3,epsarr_ps[eps],1.0,rcut))
            fw.write('%s\t %d\t %d\t %g\t %g\t %g\n'\
                     %("pair_coeff",2,3,epsarr_sg[eps],1.0,rcut))
        elif polywtperc != 1.0: #only solvent
            fw.write('%s\t %d\t %d\t %g\t %g\t %g\n'\
                     %("pair_coeff",2,2,1.0,1.0,rcut))
            fw.write('%s\t %d\t %d\t %g\t %g\t %g\n'\
                     %("pair_coeff",1,2,epsarr_ps[eps],1.0,rcut))
            
            
        fw.write('#BondCoeffs \n')
        fw.write('%s\t %s\t %g\t %g\n' %("bond_coeff","*",blist[0],blist[1]))
        fw.write('#AngleCoeffs \n')
        fw.write('%s\t %s\t %g\t %g\n' %("angle_coeff","*",alist[0],alist[1]))
        

        fw.write('#DihedralCoeffs \n')
        fw.write('%s\t %s\t %d\t %d\t %d\t'\
                 %("dihedral_coeff","*",k_phi,1,1))
        
        fw.close()
        
        
def generate_backup(maindir,srclmp,lmp_execdir,destdir):
        
    os.chdir(destdir)
    backup_init_dir = destdir + '/init_files'
    if not os.path.isdir(backup_init_dir):
        os.mkdir(backup_init_dir)
        
            
    files = glob.glob(destdir +'/*.txt')
    for f in files:
        shutil.copy2(f,backup_init_dir)
        os.remove(f)
    files = glob.glob(destdir +'/*.f90')
    for f in files:
        shutil.copy2(f,backup_init_dir)
        os.remove(f)
    files = glob.glob(destdir +'/*var*')
    for f in files:
        shutil.copy2(f,backup_init_dir)
        os.remove(f)
    files = glob.glob(destdir +'/*.mod')
    for f in files:
        shutil.copy2(f,backup_init_dir)
        os.remove(f)
    files = glob.glob(destdir +'/*.o')
    for f in files:
        shutil.copy2(f,backup_init_dir)
        os.remove(f)

    print( "Submitting Jobs..")

    subprocess.call(["qsub","jobmain.sh"])

    os.chdir(maindir)

def create_datafyle(lmp_execfile,restartfyle,outfylename):

    subprocess.call([lmp_execfile,'-r',restartfyle,outfylename])
    

def gen_extract_inp(inpfyle,trajpref,smddir,savedist,
                    tolerance,axis,comflag,pbcflag):

    os.chdir(smddir)
    smdfyles = smddir + '/' + trajpref
    traj_fyles = glob.glob(smdfyles)

    if traj_fyles == []:
        outfyle = "nil"
    else:
        latest_fyle = max(traj_fyles,key=os.path.getctime)
        extract_latestfyle = latest_fyle.split('/')
        latesttraj = extract_latestfyle[len(extract_latestfyle)-1]
        outfyle = 'extract_inp.txt'
        fr = open(inpfyle,'r')
        fw = open(outfyle,'w')

        fid = fr.read().replace("py_smdtrajfyle",latesttraj).\
              replace("py_savedist",str(savedist)).\
              replace("py_tol", str(tolerance)).\
              replace("py_axis",str(axis)).\
              replace("py_comflag", str(comflag)).\
              replace("py_pbc", str(pbcflag))
        fw.write(fid)
        fw.close()
        fr.close()

    return outfyle

def gen_launch_fyl(dirname,inpfile,outfile,graftval,epsval):
    
    os.chdir(dirname)
    headstr = dirname.split('/')
    headstr = headstr[len(headstr)-1]
    fr = open(inpfile,'r')
    fw = open(outfile,'w')
    jobstr = headstr+ "_" + str(graftval) + "_" + str(epsval)
    fid = fr.read().replace("py_jobname",jobstr)
    fw.write(fid)
    fw.close()
    fr.close()
    os.remove(inpfile)


def mult_launch_fyl(dirname,inpfile,anafile,objfile,graftval,epsval,
                    submit_file):
    os.chdir(dirname)
    fr = open(inpfile,'r')
    fw = open(submit_file,'w')
    jobstr =  "ana_" + str(graftval) + "_" + str(epsval)
    obj_name = "./" + objfile
    fid = fr.read().replace("py_jobname",jobstr).\
          replace("py_objfile",obj_name).\
          replace("py_anafile",anafile)
    fw.write(fid)
    fw.close()
    fr.close()



def cleanup_files(dirname, extval):

    os.chdir(dirname)
    
    backupdir = dirname + '/files_' + extval
    if not os.path.isdir(backupdir):
        os.mkdir(backupdir)

    filestring = './*.' + extval
    allfiles = glob.glob(filestring)

    if allfiles != []:

        for fscnt in range(len(allfiles)):
        
            desfyl = backupdir + '/' + allfiles[fscnt]
            shutil.move(allfiles[fscnt],desfyl)


def create_umbcolfile(dirumb,umbinpfyle,centerval,fconst):
    
    os.chdir(dirumb)
    outfyle = 'umbcolfile'
    fr = open(umbinpfyle,'r')
    fw = open(outfyle,'w')
    
    fid = fr.read().replace("py_cval",str(centerval)).\
          replace("py_fcon",str(fconst))
    fw.write(fid)
    fw.close()
    fr.close()


def remove_chunks(dirname, fylstr):
    fyles = dirname + '/' + fylstr
    fyllist = glob.glob(fyles)
    
    for rcnt in range(len(fyllist)):
        if not os.path.isdir(fyllist[rcnt]):
            os.remove(fyllist[rcnt])


def check_integer(inpval):

    try:
        number = int(inpval)
        isstr = -1
    except ValueError:
        isstr = 1

    return isstr
