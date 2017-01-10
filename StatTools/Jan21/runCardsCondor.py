import os
import sys

cwd = os.getcwd()
ifile = open(sys.argv[1])
jobname = ""
if (len(sys.argv)>2):
    jobname = sys.argv[2]
print sys.argv[1],sys.argv[2]
print jobname
nJobs = 0
for line in ifile:
    submit_text = ''' universe = vanilla
Executable = condor_runscript_%s_%i.sh
Should_Transfer_Files = YES
Output = stdout_%s_%i
Error  = stderr_%s_%i
Log    = log_%s_%i
Notification = never
WhenToTransferOutput=On_Exit
Queue  1
    ''' % (jobname,nJobs,jobname,nJobs,jobname,nJobs,jobname,nJobs)
    runscript_text = '''export ORIG_DIR=$PWD
cd %s
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh` \n''' % cwd
    runscript_text += line
    if (jobname == ""):
        submitfile = open("submit_%i.txt" % nJobs, "w")
        runscriptfile = open("condor_runscript_%i.sh" % nJobs, "w")
    else:
        submitfile = open("submit_%s_%i.txt" % (jobname,nJobs), "w")
        runscriptfile = open("condor_runscript_%s_%i.sh" % (jobname,nJobs), "w")
    submitfile.write(submit_text)
    runscriptfile.write(runscript_text)
    submitfile.close()
    runscriptfile.close()
    nJobs += 1
print "submitting %i jobs..." % nJobs
for i in range(nJobs):
    os.system("condor_submit submit_%s_%i.txt" %(jobname,i))
