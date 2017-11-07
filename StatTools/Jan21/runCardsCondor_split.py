import os
import sys

samples = ["ZH_hbb","WH_hbb","s_Top","TT","Wj0b","Wj1b","Wj2b","VVHF","VVLF","Zj0b","Zj1b","Zj2b","data_obs"]

cwd = os.getcwd()
ifile = open(sys.argv[1])
jobname = ""
doEWK = False

if (len(sys.argv)>2):
    jobname = sys.argv[2]
if (len(sys.argv)>3):
    samples_to_run = sys.argv[3].split(',')
if (len(sys.argv)>4):
    doEWK = bool(sys.argv[4])

if doEWK:
    samples = ["s_Top","TT","Wjets","VV","Zjets","data_obs","EWKWJets","IntEWKWJets","QCD_data"]

if len(samples_to_run) > 0:
    ## run on sub-set of analysis samples
    samples = []
    for sample in samples_to_run:
        samples.append(sample)
print "samples = ",samples
    
print sys.argv[1],sys.argv[2]
print jobname
nJobs = 0
for line in ifile:
    print line
    for sample in samples:
        cmd = line.replace("SAMP",sample)
        submit_text = '''universe = vanilla
Executable     = condor_runscript_%s_%i_%s.sh
Should_Transfer_Files     = YES
Output     = stdout_%s_%i_%s
Error      = stderr_%s_%i_%s
Log        = log_%s_%i_%s
Notification     = never
transfer_input_files = ../../splitSamples.py,../../systematics_ewk_Wmn.txt,../../systematics_ewk_Wen.txt,../../systematics_Wmn_demo.txt
WhenToTransferOutput=On_Exit
Queue      1
        ''' % (jobname,nJobs,sample,jobname,nJobs,sample,jobname,nJobs,sample,jobname,nJobs,sample)
        runscript_text = '''export ORIG_DIR=$PWD
#cd     %s
cd /cvmfs/cms.cern.ch/slc6_amd64_gcc493/cms/cmssw-patch/CMSSW_7_6_3_patch2/src
source     /cvmfs/cms.cern.ch/cmsset_default.sh
eval     `scramv1 runtime -sh`
cd $ORIG_DIR
 \n''' % cwd
        runscript_text += cmd
        submitfile = open("submit_%s_%i_%s.txt" % (jobname,nJobs,sample), "w")
        runscriptfile = open("condor_runscript_%s_%i_%s.sh" % (jobname,nJobs,sample), "w")
        submitfile.write(submit_text)
        runscriptfile.write(runscript_text)
        submitfile.close()
        runscriptfile.close()
    nJobs += 1
print "submitting %i jobs times %i samples" % (nJobs,len(samples))
for i in range(nJobs):
    for sample in samples:
        #print "condor_submit submit_%s_%i_%s.txt" %(jobname,i,sample)
        os.system("condor_submit submit_%s_%i_%s.txt" %(jobname,i,sample))
