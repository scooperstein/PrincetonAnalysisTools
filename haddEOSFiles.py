
def findAllRootFiles(value, site,doFilter=False,doRecursive=True):
    samplepaths = []
    #if value.find("Run")!=-1 or value.find("QCD")!=-1 or value.find("output_allmc.root")!=-1: return samplepaths
    #if value.find("Run")==-1: return samplepaths
    if doFilter and (value.find("TT_powheg")!=-1 or value.find("H125")!=-1 or value.find("Run")!=-1): return samplepaths
    #if value.find("TT_powheg")==-1: return samplepaths
    #print value
    if value.find("/store") is 0:
        import subprocess
        siteIP = "root://cmseos.fnal.gov"
        #siteIP = "root://131.225.207.127:1094"
        if (site == "CERN"): 
            #siteIP = "root://188.184.38.46:1094"
            ##siteIP = "root://cmseos.cern.ch"
            onlyFiles = subprocess.check_output(["eos", "ls", value]).split('\n')
        elif (site == "FNAL"):
            #onlyFiles = subprocess.check_output(["/cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw/CMSSW_7_4_14/external/slc6_amd64_gcc491/bin/xrdfs", siteIP, "ls", value]).split('\n')
            onlyFiles = subprocess.check_output(["xrdfs", siteIP, "ls", value]).split('\n')
        #print "onlyFiles = ",onlyFiles
        for filepath in onlyFiles:
            if (filepath == ""): continue
            #filepath = "root://xrootd-cms.infn.it/" + filepath
            #filepath = "root://cmsxrootd.fnal.gov/" + filepath
            if filepath.find(".root") != -1:
                ##if filepath.find("V27_EWK_SR_Oct17/output")!=-1: continue
                filepath = siteIP + "/" + filepath
                samplepaths.append(filepath)
            elif filepath.find("haddjobs")==-1 and filepath.find("sum_")==-1 and filepath.find("/log/")==-1 and filepath.find(".log")==-1 and filepath.find(".submit")==-1 and filepath.find(".stderr")==-1 and filepath.find(".stdout")==-1 and filepath.find(".sh")==-1:
                if doRecursive:
                    samplepaths.extend(findAllRootFiles(filepath,site,doFilter,doRecursive))
                else:
                    filepath = siteIP + "/" + filepath
                    samplepaths.append(filepath)
    else:
        from os import listdir
        from os.path import isfile, join
        print "here"
        onlyfiles = [ f for f in listdir(str(value)) if isfile(join(str(value),f)) ]
        for rootfile in onlyfiles:
            print rootfile
            if rootfile.find(".root") != -1:
                samplepaths.append(str(value)+"/"+str(rootfile))
            else:
                if doRecursive:
                    samplepaths.extend(findAllRootFiles(str(value)+"/"+str(rootfile)),site,doFilter,doRecursive)
    return samplepaths

import sys
import os
##region = sys.argv[1]
##kind = sys.argv[2]

useCondor = True

#ipath = "/store/user/sbc01/VHbbAnalysisNtuples/V25_Wlnu_%s_Feb13" % region
#ipath = "/eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V25_Wlnu_%s_March28" % region
#ipath = "/store/user/sbc01/VHbbAnalysisNtuples/V25_Wlnu_%s_March28" % region
##ipath = "/store/user/sbc01/VHbbAnalysisNtuples/V27_EWK_SR_Oct17/"
ipath = sys.argv[1]
site = "FNAL"
if len(sys.argv)>2:
    site = sys.argv[2]
if site != "FNAL" and site != "CERN":
    print "unknown site: %s, quitting..." % site

if site == "FNAL":
    ipath_short = ipath.replace("/eos/uscms/store","/store")
elif site == "CERN":
    ipath_short = ipath.replace("/eos/cms/store","/store")

samplepaths = findAllRootFiles(ipath_short,site,False,False)
#print samplepaths
outstring = ""
os.mkdir(ipath+"/haddjobs")
tmp = ipath.split('/')
submission_dir = tmp[len(tmp)-2] + "_haddjobsubmission"
print "job submission files will be created in: ",submission_dir
os.mkdir(submission_dir)
submitfiles = []
for samplepath in samplepaths:
    if samplepath.find("sum_")!=-1: continue
    print samplepath
    #sample = samplepath[samplepath.find("March28/")+8:]
    #sample = samplepath[samplepath.find("Oct17/")+6:]
    sample = samplepath[samplepath.rfind("/")+1:]
    print sample
    ##if (sample.find("WJets_madgraph")==-1 and sample.find("WJets-HT")==-1): continue
    #if (sample.find("Wplus")==-1 and sample.find("Wminus")==-1): continue
    samplefiles = findAllRootFiles("%s/%s" % (ipath_short,sample),"FNAL")
    outstring += "hadd sum_%s.root " % sample
    #jobtext = "cd /uscms_data/d3/sbc01/HbbAnalysis13TeV/CMSSW_7_6_3_patch2/src/\n" 
    jobtext = "cd /cvmfs/cms.cern.ch/slc6_amd64_gcc493/cms/cmssw-patch/CMSSW_7_6_3_patch2/src/\n" 
    jobtext += "source /cvmfs/cms.cern.ch/cmsset_default.sh\n"
    jobtext += "eval `scramv1 runtime -sh`\n"
    jobtext += "cd - \n"
    jobtext += "hadd sum_%s.root " % sample
    for samplefile in samplefiles:
        outstring += samplefile + " "
        jobtext += samplefile + " "
    jobtext += "\n"
    outstring += "\n"
    xrdcp_string = ""
    if site == "FNAL":
        xrdcp_string = "xrdcp sum_%s.root root://cmseos.fnal.gov:/%s/haddjobs/sum_%s.root\n" % (sample,ipath_short,sample)
    elif site == "CERN":
        xrdcp_string = "xrdcp sum_%s.root root://eoscms.cern.ch:/%s/haddjobs/sum_%s.root\n" % (sample,ipath_short,sample)
    jobtext += xrdcp_string
    jobtext += "rm sum_%s.root\n" % sample
    runscript = open(submission_dir + "/" +sample+".sh","w")
    runscript.write(jobtext)
    submittext =  "universe = vanilla\n"
    submittext += "Executable = %s/%s.sh\n" % (submission_dir,sample)
    submittext += "initialdir = %s \n" %  submission_dir
    submittext += "Should_Transfer_Files = YES\n"
    submittext += "Output = %s.stdout\n" % sample
    submittext += "Error  = %s.stderr\n" % sample
    submittext += "Log    = %s.log\n"    % sample
    submittext += "Notification = never\n"
    submittext += "WhenToTransferOutput=On_Exit\n"
    submittext += "Queue  1\n"
    submitFile = open(submission_dir + "/" +sample+".submit","w")
    submitFile.write(submittext)
    submitFile.close()
    runscript.close()
    submitfiles.append(submission_dir + "/" +sample+".submit")
    
for submitfile in submitfiles:
    if useCondor:
        print "condor_submit ",submitfile
        os.system("condor_submit %s" %submitfile)
    else:
        print 'bsub -R "pool>30000" -q 1nh -J job1 < %s' % submitfile.replace(".submit",".sh")
        os.system('bsub -R "pool>30000" -q 1nh -J job1 < %s' % submitfile.replace(".submit",".sh"))
#print outstring


#if (kind == "ttpowheg"):
#    rootfiles = findAllRootFiles("%s/TT_powheg" % ipath,"FNAL")
#elif (kind == "signal"):
#    rootfiles = findAllRootFiles("%s/WplusH125_powheg" %ipath,"FNAL")
#    rootfiles.extend(findAllRootFiles("%s/WminusH125_powheg" %ipath,"FNAL"))
#    rootfiles.extend(findAllRootFiles("%s/WH125"%ipath,"FNAL"))
#    rootfiles.extend(findAllRootFiles("%s/ZH125_powheg"%ipath,"FNAL"))
#elif (kind == "data"):
#    rootfiles = findAllRootFiles("%s/Run2016BToG_Ele" % ipath,"FNAL")
#    rootfiles.extend(findAllRootFiles("%s/Run2016BToG_Mu" %ipath,"FNAL"))
#elif (kind == "mudata"):
#    rootfiles = findAllRootFiles("%s/Run2016BToG_Mu" % ipath,"FNAL")
#elif (kind == "eldata"):
#    rootfiles = findAllRootFiles("%s/Run2016BToG_Ele" % ipath,"FNAL")
#elif (kind == "mc"):
#    rootfiles = findAllRootFiles("%s/" %ipath,"FNAL",True)
##outstring = "hadd -f output_%s.root " % kind
##outstring = "hadd -f output_allmc.root "
##outstring = "hadd -f root://cmseos.fnal.gov//store/user/sbc01/VHbbAnalysisNtuples/V25_Wlnu_SR_Feb9/output_allmc.root "
#for rootfile in rootfiles:
#    outstring += rootfile + " "
#print outstring
