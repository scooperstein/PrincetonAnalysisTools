import ROOT
import sys
import ReadInput
from optparse import OptionParser
import os

parser = OptionParser()
parser.add_option("-c", "--config", dest="configFile", default="ewk_config.txt",
                  help="specify config file for this analysis"
)
parser.add_option("-b", "--batch", dest="runBatch", default=0, type=int, 
                  help="run analysis jobs on condor (default=0)"
)
parser.add_option("-n", "--jobName", dest="jobName", default="condor_jobs",
                  help="Specify label for condor jobs. Only to be used when running batch jobs"
)
parser.add_option("-f", "--nFilesPerJob", dest="nFilesPerJob", default=10, type=int, help="Number of input files per batch job")
parser.add_option("-s","--sample", dest="sample", default="", type=str, help="Run on only a specific sample (can be comma-separated list)")
parser.add_option("-d","--doData", dest="doData", default=-1, type=int, help="If -1 run all samples, if 0 run only MC, if 1 run only data")
(options, args) = parser.parse_args()

ROOT.gSystem.Load("AnalysisDict.so")

# reads samples, existing branches and new branches
if (options.sample != ""):
    samplesToSubmit = options.sample.split(',')
else:
    samplesToSubmit = [] # if empty run on all samples
am=ReadInput.ReadTextFile(options.configFile, "cfg", samplesToSubmit,"",options.runBatch)
am.debug=2

doData = options.doData


if (options.runBatch == False):
    print "Running locally over all samples"

    if(am.debug>100):
        am.PrintBranches()

    # loop over all the samples
    # FIXME - need to add the possibility of doing a small portion of files
    am.Loop()

else:
    print "Running analysis jobs on condor. Note that this will only work from LPC machines"
    jobName = options.jobName

    #if os.path.exists(jobName):
    #    print "Attempting to create jobs under jobName %s, which already exists!" % jobName
    #    sys.exit(0)
    os.system("mkdir -p %s" % jobName)
    os.system("mkdir -p /eos/uscms/store/user/sbc01/SkimmedHeppyNtuples/%s" % jobName)

    submitFiles = []

    #for sampleName in am.ListSampleNames():
    nFilesPerJob = options.nFilesPerJob 
    for sample in am.samples:
        if (options.sample != ""):
            #if (sample.sampleName != options.sample): continue
            if (sample.sampleName not in samplesToSubmit ):
                print "sample: ",sample.sampleName," not in list, skipping..."
                continue
        if (options.doData == 0 and sample.sampleName.find("Run")!=-1):
            print "skipping data sample: ",sample.sampleName
            continue
        if (options.doData == 1 and sample.sampleName.find("Run")==-1):
            print "skipping MC sample: ",sample.sampleName
            continue
        sampleName = sample.sampleName
        print sampleName
        os.system("mkdir -p %s/%s" % (jobName,sampleName))
        os.system("mkdir -p /eos/uscms/store/user/sbc01/SkimmedHeppyNtuples/%s/%s" % (jobName,sampleName))
        nProcJobs = 0
        nFiles = len(sample.files)
        nJobs = nFiles / nFilesPerJob + 1
        #for filename in sample.files:
        for i in range(nJobs):
            filesToRun = ""
            for j in range(nFilesPerJob):
                index = i*nFilesPerJob + j
                if (index >= nFiles): continue
                filesToRun += "%s," % sample.files[index]
            filesToRun = filesToRun[:-1] # remove trailing ','
            if (filesToRun == ""): continue

            nProcJobs += 1
            fname = "%s/%s/job%i.submit" % (jobName, sampleName,nProcJobs)
            submitFile = open(fname, "w")
            content =  "universe = vanilla\n"
            content += "Executable = %s/condor_runscript.sh\n" % jobName
            content += "Arguments = %s %s %s output_%s_%i.root\n" % (options.configFile, sampleName, filesToRun,sampleName, nProcJobs)
            #content += "Arguments = %s %s %s %i\n" % (options.configFile, sampleName, filename, nProcJobs)
            #content += "Requirements   =  OpSys == 'LINUX' && (Arch =='INTEL' || Arch =='x86_64')\n"
            content += "initialdir = %s/%s\n" % (jobName,sampleName)
            #content += "Should_Transfer_Files = YES\n"
            content += "Output = %i.stdout\n" % nProcJobs
            content += "Error  = %i.stderr\n" % nProcJobs
            content += "Log    = %i.log\n"    % nProcJobs
            content += "Notification = never\n"
            #content += "WhenToTransferOutput=On_Exit\n"
            #content += "transfer_input_files = ../../%s,../../cfg/samples.txt,../../cfg/earlybranches.txt,../../cfg/existingbranches.txt,../../cfg/newbranches.txt,../../cfg/bdtsettings.txt,../../cfg/reg1_settings.txt,../../cfg/reg2_settings.txt,../../cfg/settings.txt,../../aux/TMVARegression_BDTG_ttbar_Nov23.weights.xml,../../aux/TMVA_13TeV_Dec14_3000_5_H125Sig_0b1b2bWjetsTTbarBkg_Mjj_BDT.weights.xml,../../aux/MuonIso_Z_RunCD_Reco74X_Dec1.json,../../aux/SingleMuonTrigger_Z_RunCD_Reco74X_Dec1.json,../../aux/MuonID_Z_RunCD_Reco74X_Dec1.json,../../aux/CutBasedID_TightWP.json,../../aux/CutBasedID_LooseWP.json,../../RunSample.py,../../../AnalysisDict.so,../../cfg/systematics.txt,../../cfg/scalefactors.txt\n" % options.configFile
            content += "Queue  1\n"
            #print content
            submitFile.write(content)
            submitFile.close()
            submitFiles.append(fname)

    condor_runscript_text = '''

        #!/bin/bash
        ##
        ## Script to run Analysis Manager jobs on Condor from the LPC
        ## Author: Stephane Cooperstein
        ##
        ## Argument 1: Analysis config file
        ## Argument 2: Sample name to run on
        ##

        export ORIG_DIR=$PWD
        # Set up environment
        echo "setting up the environment"
        cd %s/..
        source /cvmfs/cms.cern.ch/cmsset_default.sh
        eval `scramv1 runtime -sh`
        pwd
        ls
        #source env.sh
        source %s/../env.sh
        
        echo "successfully set up the enviroment"
  
        cd $ORIG_DIR  
    
        echo "copying over necessary files"
        cp -r %s/cfg .
        cp -r %s/aux .
        cp %s/RunSkim.py .
        cp %s/../AnalysisDict.so .
        cp -r %s/*.txt . 

        echo "running RunSkim.py"
        echo $ORIG_DIR/$4
        python RunSkim.py $1 $2 $3 $ORIG_DIR/$4
        echo "done skimming, now copying output to EOS"
        xrdcp -f $ORIG_DIR/$4 root://cmseos.fnal.gov//store/user/sbc01/SkimmedHeppyNtuples/%s/$2
        rm $ORIG_DIR/$4
        echo "all done!" ''' % (os.getcwd(), os.getcwd(), os.getcwd(), os.getcwd(), os.getcwd(), os.getcwd(), os.getcwd(), jobName )

    runscript = open("%s/condor_runscript.sh" % (jobName) , "w")
    runscript.write(condor_runscript_text)
    runscript.close()

    ## Send job to condor
    print "Submit files created, sending jobs to Condor..."
    for file in submitFiles:
        print("condor_submit %s" % file)
        os.system("condor_submit %s" % file)

        

