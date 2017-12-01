import argparse
import os.path
import sys
import tarfile
import re
import subprocess
import glob
import itertools

parser = argparse.ArgumentParser(description='Submit beam data jobs.')

parser.add_argument("-f", "--first-run", dest="firstrun",
                    required=True,
                    help="First run.")
parser.add_argument("-l", "--last-run", dest="lastrun",
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

parser.add_argument("-t", "--track_energy_distortion", type=float, nargs="+")
parser.add_argument("-s", "--shower_energy_distortion", type=float, nargs="+")
parser.add_argument("--track_energy_distortion_by_percent", action="store_true")
parser.add_argument("--shower_energy_distortion_by_percent", action="store_true")
parser.add_argument("--accept_np", action="store_true")
parser.add_argument("--accept_ntrk", action="store_true")
parser.add_argument("--all_combinations", action="store_true")

args = parser.parse_args()

track_energy_distortion = args.track_energy_distortion
shower_energy_distortion = args.shower_energy_distortion
if args.all_combinations:
   s_len = len(shower_energy_distortion) 
   t_len = len(track_energy_distortion)
   track_energy_distortion = list(itertools.chain.from_iterable(itertools.repeat(x, s_len) for x in track_energy_distortion))
   shower_energy_distortion *= t_len

assert(len(track_energy_distortion) == len(shower_energy_distortion))

n_selections = 1 if track_energy_distortion is None else len(track_energy_distortion)

#now create jobfiles_*.tar that is shipped with the job
#this includes the executable
#tar -cf jobfiles.tar --transform '!^[^/]*/!!' file1 file2 file3
tarfilename="jobfiles_%i.tar.bz2"%os.getpid()
outtar = tarfile.open(tarfilename, mode='w:bz2')
outtar.add(args.f_list_name, arcname=args.f_list_name)
outtar.add("truth_selection_multiple_selections",arcname="truth_selection")
outtar.add("data/dedx_pdfs.root", arcname="dedx_pdfs.root")
outtar.close()

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
truth_selection_args += " -d "
for d in range(n_selections):
    truth_selection_args += "%i " % d
if args.accept_np:
    truth_selection_args += " --accept_np"
if args.accept_ntrk:
    truth_selection_args += " --accept_ntrk"

ofstr='''
#!/bin/bash

export _RUN_NUMBER=$((PROCESS+%(firstrun)s))

export THISFILE=$(sed "$((_RUN_NUMBER+1))q;d" %(f_list_name)s)

echo "RunArgs: " >>${_RUN_NUMBER}.out 2>&1
echo "%(args)s" >>${_RUN_NUMBER}.out 2>&1

echo "Running $0 on "$HOSTNAME >>${_RUN_NUMBER}.out 2>&1
echo "Cluster: " ${CLUSTER} >>${_RUN_NUMBER}.out 2>&1
echo "Process: " ${PROCESS} >>${_RUN_NUMBER}.out 2>&1

echo " Sourcing everything...." >>${_RUN_NUMBER}.out 2>&1

#source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setup >>${_RUN_NUMBER}.out 2>&1
#source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone.sh >>${_RUN_NUMBER}.out 2>&1
#source /cvmfs/fermilab.opensciencegrid.org/products/larsoft/setup >>${_RUN_NUMBER}.out 2>&1

source /grid/fermiapp/products/uboone/setup_uboone.sh >>${_RUN_NUMBER}.out 2>&1
source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setup >>${_RUN_NUMBER}.out 2>&1
source /cvmfs/fermilab.opensciencegrid.org/products/larsoft/setup >>${_RUN_NUMBER}.out 2>&1

setup larsoftobj v1_13_00 -q e10:prof && echo "Setup Larsoftobj" >>${_RUN_NUMBER}.out 2>&1
setup uboonecode v06_26_01_07 -q e10:prof && echo "Setup Uboonecode" >>${_RUN_NUMBER}.out 2>&1

mv ${_RUN_NUMBER}.out ${_CONDOR_SCRATCH_DIR}/.
cd ${_CONDOR_SCRATCH_DIR}

cp $INPUT_TAR_FILE . 
tar -jvxf `basename ${INPUT_TAR_FILE}` >>${_RUN_NUMBER}.out 2>&1

echo $THISFILE >>${_RUN_NUMBER}.out 2>&1

echo " Copying File...." >>${_RUN_NUMBER}.out 2>&1
ifdh cp $THISFILE ${PWD}/test.root >>${_RUN_NUMBER}.out 2>&1

echo "What's in here? "
ls >>${_RUN_NUMBER}.out 2>&1

echo " Run time! " >>${_RUN_NUMBER}.out 2>&1
echo "With args: %(args)s" >>${_RUN_NUMBER}.out 2>&1

./truth_selection test.root %(args)s >>${_RUN_NUMBER}.out 2>&1

mkdir ${_RUN_NUMBER}

cp ${_RUN_NUMBER}.out ${_RUN_NUMBER}/.
for i in {0..%(max_i_selection)s} 
do
  cp out${i}.root ${_RUN_NUMBER}/.
done

ifdh mkdir %(outputdir)s/
ifdh mkdir %(outputdir)s/${_RUN_NUMBER}
ifdh cp -r ${_RUN_NUMBER} %(outputdir)s/${_RUN_NUMBER}/

'''%{'firstrun':args.firstrun,'outputdir':args.outputpath, 'args': truth_selection_args, 'f_list_name': args.f_list_name,'max_i_selection':n_selections-1}

runjobfname="runjob_%i.sh"%os.getpid()
of=open(runjobfname,'w')
of.write(ofstr)
of.close()

n_jobs = int(args.lastrun)-int(args.firstrun)+1

cmd="jobsub_submit --memory=1000MB --group=uboone -N %i --tar_file_name=dropbox://%s file://%s"%(int(args.lastrun)-int(args.firstrun)+1,os.path.abspath(tarfilename),os.path.abspath(runjobfname))

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



