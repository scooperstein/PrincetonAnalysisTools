import os
import sys

cwd = os.getcwd()
ifile = open(sys.argv[1])

nJobs = 0
for line in ifile:
    submit_text = ''' universe = vanilla
Executable = condor_runscript_%i.sh
Should_Transfer_Files = YES
Output = stdout_%i
Error  = stderr_%i
Log    = log_%i
Notification = never
WhenToTransferOutput=On_Exit
Queue  1
    ''' % (nJobs,nJobs,nJobs,nJobs)
    runscript_text = '''export ORIG_DIR=$PWD
cd %s
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh` \n''' % cwd
    runscript_text += line
    submitfile = open("submit_%i.txt" % nJobs, "w")
    runscriptfile = open("condor_runscript_%i.sh" % nJobs, "w")
    submitfile.write(submit_text)
    runscriptfile.write(runscript_text)
    submitfile.close()
    runscriptfile.close()
    nJobs += 1
print "submitting %i jobs..." % nJobs
for i in range(nJobs):
    os.system("condor_submit submit_%i.txt" %i)
