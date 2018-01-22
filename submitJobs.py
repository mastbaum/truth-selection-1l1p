import argparse
import os.path
import sys
import tarfile
import re
import subprocess
import glob
import itertools

parser = argparse.ArgumentParser(description='Submit beam data jobs.')

# arguments -- pretty much all of these pass an argument to RunSelection.cxx
parser.add_argument("-f", "--first-file", dest="firstfile",
                    required=True,
                    help="First run.")
parser.add_argument("-l", "--last-file", dest="lastfile",
                    required=True,
                    help="Last run.")
parser.add_argument("-i", "--input-file-list", dest="f_list_name",
                    required=True,
                    help="Path to input file list.")
parser.add_argument("-o", "--output-path", dest="outputpath",
                    default="/pnfs/uboone/scratch/users/%s/truth_selection_out"%os.environ['USER'],
                    help="Path where to copy final output. (default=/pnfs/uboone/scratch/users/%s/bnb_redecay_to_gsimple)"%os.environ['USER'])
parser.add_argument("-d", "--debug",action='store_true',
                    help="Will not delete submission files in the end. Useful for debugging and will only print the submission command on screen.")
parser.add_argument("-n", "--n_files_per_run", type=int, default=1)
parser.add_argument("-T", "--n_trials",type=int, default=1)
parser.add_argument("-r", "--run_no", type=int, default=0)

parser.add_argument("-c", "--root_config_file", default=None)
parser.add_argument("-t", "--track_energy_distortion", type=float, nargs="+")
parser.add_argument("-s", "--shower_energy_distortion", type=float, nargs="+")
parser.add_argument("--track_energy_distortion_by_percent", action="store_true")
parser.add_argument("--shower_energy_distortion_by_percent", action="store_true")
parser.add_argument("--drop_np", action="store_true")
parser.add_argument("--drop_ntrk", action="store_true")
parser.add_argument("--all_combinations", action="store_true")

args = parser.parse_args()


# coherence checks on the passed in parameters 
track_energy_distortion = args.track_energy_distortion
shower_energy_distortion = args.shower_energy_distortion
if args.all_combinations:
   s_len = len(shower_energy_distortion) 
   t_len = len(track_energy_distortion)
   track_energy_distortion = list(itertools.chain.from_iterable(itertools.repeat(x, s_len) for x in track_energy_distortion))
   shower_energy_distortion *= t_len

if not track_energy_distortion is None:
    assert(len(track_energy_distortion) == len(shower_energy_distortion))
else:
    assert(shower_energy_distortion is None)

n_selections = 1 if track_energy_distortion is None else len(track_energy_distortion)

# tar up the input files

#now create jobfiles_*.tar that is shipped with the job
#this includes the executable
#tar -cf jobfiles.tar --transform '!^[^/]*/!!' file1 file2 file3
tarfilename="jobfiles_%i.tar.bz2"%os.getpid()
outtar = tarfile.open(tarfilename, mode='w:bz2')
outtar.add(args.f_list_name, arcname=args.f_list_name)
outtar.add("truth_selection",arcname="truth_selection")
if os.path.isfile("config.root"):
    outtar.add("config.root", arcname="config.root")
outtar.close()

# pass required arguments to RunSelection.cxx
truth_selection_args = ""
truth_selection_args = "-o "
for d in range(n_selections):
   truth_selection_args += "out%i.root " % d
truth_selection_args += " -n %i " % n_selections
if track_energy_distortion is not None:
    truth_selection_args += " --T_energy_distortion "
    for dist in track_energy_distortion:
        truth_selection_args += "%f " % dist
    if args.track_energy_distortion_by_percent:
        truth_selection_args += " --T_edist_by_percent"
if shower_energy_distortion is not None:
    truth_selection_args += " --S_energy_distortion "
    for dist in shower_energy_distortion:
        truth_selection_args += "%f " % dist
    if args.shower_energy_distortion_by_percent:
        truth_selection_args += " --S_edist_by_percent"
if args.root_config_file is not None:
    truth_selection_args += " --root_config_file %s" % args.root_config_file

truth_selection_args += " -d "
for d in range(n_selections):
    truth_selection_args += "%i " % d
if args.drop_np:
    truth_selection_args += " --drop_np"
if args.drop_ntrk:
    truth_selection_args += " --drop_ntrk"
truth_selection_args += " -t %i" % args.n_trials

# the shell script that will be run on the grid
ofstr='''
#!/bin/bash

export _RUN_NUMBER=$((PROCESS+%(firstrun)s))

#export THESFILE=$(sed -n "$((_RUN_NUMBER+1))q;d" %(f_list_name)s)
export MINFILENO=$((%(min_file_no)s+1+$PROCESS*%(files_per_job)s))
export THESEFILES=$(sed -n "$MINFILENO,$(($MINFILENO+%(files_per_job)s-1))p" %(f_list_name)s)

echo "RunArgs: " >>${_RUN_NUMBER}.out 2>&1
echo "%(args)s" >>${_RUN_NUMBER}.out 2>&1

echo "Running $0 on "$HOSTNAME >>${_RUN_NUMBER}.out 2>&1
echo "Cluster: " ${CLUSTER} >>${_RUN_NUMBER}.out 2>&1
echo "Process: " ${PROCESS} >>${_RUN_NUMBER}.out 2>&1

echo " Sourcing everything...." >>${_RUN_NUMBER}.out 2>&1

source /grid/fermiapp/products/uboone/setup_uboone.sh >>${_RUN_NUMBER}.out 2>&1
source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setup >>${_RUN_NUMBER}.out 2>&1
source /cvmfs/fermilab.opensciencegrid.org/products/larsoft/setup >>${_RUN_NUMBER}.out 2>&1

#setup larsoftobj v1_13_00 -q e10:prof && echo "Setup Larsoftobj" >>${_RUN_NUMBER}.out 2>&1
#setup uboonecode v06_26_01_07 -q e10:prof && echo "Setup Uboonecode" >>${_RUN_NUMBER}.out 2>&1

setup gallery v1_05_03 -q prof:e14:nu
setup uboonecode v06_55_00 -q prof:e14

mv ${_RUN_NUMBER}.out ${_CONDOR_SCRATCH_DIR}/.
cd ${_CONDOR_SCRATCH_DIR}

cp $INPUT_TAR_FILE . 
tar -jvxf `basename ${INPUT_TAR_FILE}` >>${_RUN_NUMBER}.out 2>&1

echo $THISFILE >>${_RUN_NUMBER}.out 2>&1

export FILES=($THESEFILES)
echo " Copying File...." >>${_RUN_NUMBER}.out 2>&1

export FILELIST=""
for file in "${FILES[@]}"
do 
    echo "$file ${PWD}/$(basename $file)" >> filelist
    FILELIST="$FILELIST $(basename $file)"
done
ifdh cp -f filelist

echo "What's in here? " >>${_RUN_NUMBER}.out 2>&1
ls >>${_RUN_NUMBER}.out 2>&1

echo " Run time! " >>${_RUN_NUMBER}.out 2>&1
echo "With args: %(args)s" >>${_RUN_NUMBER}.out 2>&1
echo "With Input Files: ${FILELIST}" >>${_RUN_NUMBER}.out 2>&1

./truth_selection $FILELIST %(args)s >>${_RUN_NUMBER}.out 2>&1

mkdir ${_RUN_NUMBER}

cp ${_RUN_NUMBER}.out ${_RUN_NUMBER}/.
for i in {0..%(max_i_selection)s} 
do
  cp out${i}.root ${_RUN_NUMBER}/.
done

ifdh mkdir %(outputdir)s/${_RUN_NUMBER}
ifdh cp -r ${_RUN_NUMBER} %(outputdir)s/${_RUN_NUMBER}/

'''%{'firstrun':args.run_no,'min_file_no':args.firstfile,'outputdir':args.outputpath, 'args': truth_selection_args, 'f_list_name': args.f_list_name,'max_i_selection':n_selections-1,'files_per_job':args.n_files_per_run}

runjobfname="runjob_%i.sh"%os.getpid()
of=open(runjobfname,'w')
of.write(ofstr)
of.close()

n_files = int(args.lastfile)-int(args.firstfile)+1
n_jobs = n_files / args.n_files_per_run
if n_files % args.n_files_per_run != 0:
    n_jobs += 1

# the jobsub command for actually running stuff on grid
cmd="jobsub_submit --memory=1000MB --disk=10GB --group=uboone -N %i --tar_file_name=dropbox://%s file://%s"%(n_jobs,os.path.abspath(tarfilename),os.path.abspath(runjobfname))

if (not args.debug):
    print "Running submit cmd:"
    print cmd
    os.system(cmd)
else:
    print "Would have ran:"
    print cmd

#Delete temp files unless debugging
if (not args.debug):
    os.remove(tarfilename)
    os.remove(runjobfname)



