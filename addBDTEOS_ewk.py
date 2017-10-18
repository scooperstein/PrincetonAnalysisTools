import os
import sys

idir = sys.argv[1]
curr_dir = os.getcwd()
os.system("mkdir -p addBDTEOS_ewk")
curr_dir = curr_dir + "/addBDTEOS_ewk"

submitfiles = []
for subdir, dirs, files in os.walk(idir):
    #print subdir, dirs, files
    i = 0 
    for filename in files:
        if subdir.find("haddjobs")!=-1: continue
        if filename.find(".root") == -1: continue
        #if filename.find("_3.root")==-1 or filename.find("sum_")==-1: continue
        #if filename.find("_weighted2.root")==-1 or filename.find("sum_")==-1: continue
        label = "Aug30/"
        sample = subdir[subdir.find(label)+len(label):]
        #sample = filename[4:filename.find("_weighted2.root")]
        #sample = filename[4:filename.find("_3.root")]
        #print sample
        if subdir.find("ZZ_fil")==-1: continue
        #print subdir
        filepath = idir + "/" + sample + "/" + filename
        filepath = filepath.replace("/eos/uscms/","root://cmseos.fnal.gov//")
        jobtext = "export ORIG_DIR=$PWD\n"
        jobtext += "cd /uscms_data/d3/sbc01/HbbAnalysis13TeV/CMSSW_7_6_3_patch2/src/\n"
        jobtext += "source /cvmfs/cms.cern.ch/cmsset_default.sh\n"
        jobtext += "eval `scramv1 runtime -sh`\n"
        jobtext += "cd - \n"
        jobtext += "cd /uscms_data/d3/sbc01/HbbAnalysis13TeV/CMSSW_7_6_3_patch2/src/PrincetonAnalysisTools/scripts\n"
        #jobtext += "python bTagSF_2.py %s $ORIG_DIR/%s\n" % (filepath,filename.replace("_3.root","_weighted2.root"))
        #jobtext += 'root -l .x addBDTToTree.C\(\"%s\" \"${ORIG_DIR}/%s\"\) \n' % (filepath,filename.replace("_weighted2.root","_weighted3.root"))
        ofilename = filename.replace("output","output2")
        jobtext += './addBDTToTree_ewk %s $ORIG_DIR/%s \n' % (filepath,ofilename)
        jobtext += "cd - \n"
        jobtext += "xrdcp -f $ORIG_DIR/%s root://cmseos.fnal.gov/%s/%s \n" % (ofilename,subdir.replace("/eos/uscms",""),ofilename)
        jobtext += "rm %s \n" % ofilename
        #jobtext += "cd /uscms_data/d3/sbc01/HbbAnalysis13TeV/CMSSW_7_6_3_patch2/src/PrincetonAnalysisTools/scripts/\n"
        #jobtext += "python addVPtCorrFactors_split.py $ORIG_DIR/%s $ORIG_DIR/%s\n" % (filename.replace(".root","_2.root"),filename.replace(".root","_3.root"))
        #jobtext += "cd - \n"
        #jobtext += "rm %s\n" % filename.replace(".root","_2.root")
        #jobtext += "cd /uscms_data/d3/sbc01/HbbAnalysis13TeV/CMSSW_7_6_3_patch2/src/PrincetonAnalysisTools/scripts/VHSignalCorrections\n"
        #jobtext += "python addVHSignalCorrections.py  $ORIG_DIR/%s $ORIG_DIR/%s\n" % (filename.replace(".root","_3.root"),filename.replace(".root","_weighted.root"))
        #jobtext += "cd - \n"
        #jobtext += "rm %s\n" % filename.replace(".root","_3.root")
        #jobtext += "rm %s %s" % (filename.replace(".root","_2.root"),filename.replace(".root","_3.root"))
        
        runscript = open(curr_dir + "/" +sample+"_addbdt_"+str(i)+".sh","w")
        runscript.write(jobtext)
        runscript.close()
        submittext =  "universe = vanilla\n"
        submittext += "Executable = %s/%s_addbdt_%i.sh\n" % (curr_dir,sample,i)
        submittext += "initialdir = %s \n" %  curr_dir
        submittext += "Should_Transfer_Files = YES\n"
        submittext += "Output = %s_addbdt_%i.stdout\n" % (sample,i)
        submittext += "Error  = %s_addbdt_%i.stderr\n" % (sample,i)
        submittext += "Log    = %s_addbdt_%i.log\n"    % (sample,i)
        submittext += "Notification = never\n"
        submittext += "WhenToTransferOutput=On_Exit\n"
        submittext += "Queue  1\n"
        submitFile = open(curr_dir + "/" + sample+"_addbdt_"+str(i)+".submit","w")
        submitFile.write(submittext)
        submitFile.close()
        submitfiles.append(curr_dir + "/" + sample+"_addbdt_"+str(i)+".submit")
        i = i + 1

for submitfile in submitfiles:
    print "condor_submit ",submitfile
    os.system("condor_submit %s" %submitfile)
