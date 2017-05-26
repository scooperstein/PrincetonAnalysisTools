
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
            siteIP = "root://188.184.38.46:1094"
        onlyFiles = subprocess.check_output(["/cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw/CMSSW_7_4_14/external/slc6_amd64_gcc491/bin/xrdfs", siteIP, "ls", value]).split('\n')
        #print "onlyFiles = ",onlyFiles
        for filepath in onlyFiles:
            if (filepath == ""): continue
            #filepath = "root://xrootd-cms.infn.it/" + filepath
            #filepath = "root://cmsxrootd.fnal.gov/" + filepath
            if filepath.find(".root") != -1:
                if filepath.find("V25_EWK_SR_May3/output")!=-1: continue
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
region = sys.argv[1]
kind = sys.argv[2]

#ipath = "/store/user/sbc01/VHbbAnalysisNtuples/V25_Wlnu_%s_Feb13" % region
#ipath = "/eos/uscms/store/user/sbc01/VHbbAnalysisNtuples/V25_Wlnu_%s_March28" % region
#ipath = "/store/user/sbc01/VHbbAnalysisNtuples/V25_Wlnu_%s_March28" % region
ipath = "/store/user/sbc01/VHbbAnalysisNtuples/V25_EWK_SR_May3"

samplepaths = findAllRootFiles(ipath,"FNAL",False,False)
#print samplepaths
outstring = ""
os.mkdir("/eos/uscms"+ipath+"/haddjobs")
submitfiles = []
for samplepath in samplepaths:
    if samplepath.find("sum_")!=-1: continue
    print samplepath
    #sample = samplepath[samplepath.find("March28/")+8:]
    sample = samplepath[samplepath.find("May3/")+5:]
    samplefiles = findAllRootFiles("%s/%s" % (ipath,sample),"FNAL")
    outstring += "hadd sum_%s.root " % sample
    jobtext = "cd /uscms_data/d3/sbc01/HbbAnalysis13TeV/CMSSW_7_6_3_patch2/src/\n" 
    jobtext += "source /cvmfs/cms.cern.ch/cmsset_default.sh\n"
    jobtext += "eval `scramv1 runtime -sh`\n"
    jobtext += "cd - \n"
    jobtext += "hadd sum_%s.root " % sample
    for samplefile in samplefiles:
        outstring += samplefile + " "
        jobtext += samplefile + " "
    outstring += "\n"
    runscript = open("/eos/uscms"+ipath+"/haddjobs/"+sample+".sh","w")
    runscript.write(jobtext)
    submittext =  "universe = vanilla\n"
    submittext += "Executable = /eos/uscms%s/haddjobs/%s.sh\n" % (ipath,sample)
    submittext += "initialdir = /eos/uscms%s/haddjobs \n" %  ipath
    submittext += "Should_Transfer_Files = YES\n"
    submittext += "Output = %s.stdout\n" % sample
    submittext += "Error  = %s.stderr\n" % sample
    submittext += "Log    = %s.log\n"    % sample
    submittext += "Notification = never\n"
    submittext += "WhenToTransferOutput=On_Exit\n"
    submittext += "Queue  1\n"
    submitFile = open("/eos/uscms"+ipath+"/haddjobs/"+sample+".submit","w")
    submitFile.write(submittext)
    submitFile.close()
    runscript.close()
    submitfiles.append("/eos/uscms"+ipath+"/haddjobs/"+sample+".submit")
    
for submitfile in submitfiles:
    print "condor_submit ",submitfile
    os.system("condor_submit %s" %submitfile)
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
