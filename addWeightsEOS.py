import os
import sys

idir = sys.argv[1]

submitfiles = []
for subdir, dirs, files in os.walk(idir):
    for filename in files:
        if filename.find(".root")==-1 or filename.find("sum_")==-1: continue
        sample = filename[4:filename.find(".root")]
        print sample
        filepath = idir + "/" + filename
        filepath = filepath.replace("/eos/uscms/","root://cmseos.fnal.gov//")
        jobtext = "export ORIG_DIR=$PWD\n"
        jobtext += "cd /uscms_data/d3/sbc01/HbbAnalysis13TeV/CMSSW_7_6_3_patch2/src/\n"
        jobtext += "source /cvmfs/cms.cern.ch/cmsset_default.sh\n"
        jobtext += "eval `scramv1 runtime -sh`\n"
        jobtext += "cd - \n"
        jobtext += "cd /uscms_data/d3/sbc01/HbbAnalysis13TeV/CMSSW_7_6_3_patch2/src/PrincetonAnalysisTools/scripts/bTagWeighting\n"
        jobtext += "python bTagSF.py %s $ORIG_DIR/%s\n" % (filepath,filename.replace(".root","_2.root"))
        jobtext += "cd - \n"
        jobtext += "cd /uscms_data/d3/sbc01/HbbAnalysis13TeV/CMSSW_7_6_3_patch2/src/PrincetonAnalysisTools/scripts/\n"
        jobtext += "python addVPtCorrFactors_split.py $ORIG_DIR/%s $ORIG_DIR/%s\n" % (filename.replace(".root","_2.root"),filename.replace(".root","_3.root"))
        jobtext += "cd - \n"
        jobtext += "rm %s\n" % filename.replace(".root","_2.root")
        jobtext += "cd /uscms_data/d3/sbc01/HbbAnalysis13TeV/CMSSW_7_6_3_patch2/src/PrincetonAnalysisTools/scripts/VHSignalCorrections\n"
        jobtext += "python addVHSignalCorrections.py  $ORIG_DIR/%s $ORIG_DIR/%s\n" % (filename.replace(".root","_3.root"),filename.replace(".root","_weighted.root"))
        jobtext += "cd - \n"
        jobtext += "rm %s\n" % filename.replace(".root","_3.root")
        #jobtext += "rm %s %s" % (filename.replace(".root","_2.root"),filename.replace(".root","_3.root"))
        
        runscript = open(idir + "/" +sample+"_weights.sh","w")
        runscript.write(jobtext)
        submittext =  "universe = vanilla\n"
        submittext += "Executable = %s/%s_weights.sh\n" % (idir,sample)
        submittext += "initialdir = %s \n" %  idir
        submittext += "Should_Transfer_Files = YES\n"
        submittext += "Output = %s_weights.stdout\n" % sample
        submittext += "Error  = %s_weights.stderr\n" % sample
        submittext += "Log    = %s_weights.log\n"    % sample
        submittext += "Notification = never\n"
        submittext += "WhenToTransferOutput=On_Exit\n"
        submittext += "Queue  1\n"
        submitFile = open(idir + "/" + sample+"_weights.submit","w")
        submitFile.write(submittext)
        submitFile.close()
        submitfiles.append(idir + "/" + sample+"_weights.submit")

for submitfile in submitfiles:
    print "condor_submit ",submitfile
    os.system("condor_submit %s" %submitfile)
