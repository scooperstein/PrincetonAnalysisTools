import ROOT
import sys
import ReadInput
from optparse import OptionParser
import os

parser = OptionParser()
parser.add_option("-c", "--config", dest="configFile", default="vhbb_config.txt",
                  help="specify config file for this analysis"
)
parser.add_option("-b", "--batch", dest="runBatch", default=0, type=int, 
                  help="run analysis jobs on condor (default=0)"
)
parser.add_option("-n", "--jobName", dest="jobName", default="condor_jobs",
                  help="Specify label for condor jobs. Only to be used when running batch jobs"
)
parser.add_option("-f", "--nFilesPerJob", dest="nFilesPerJob", default=10, type=int, help="Number of input files per batch job")
(options, args) = parser.parse_args()

ROOT.gSystem.Load("AnalysisDict.so")

# reads samples, existing branches and new branches
samplesToRun = [] # if empty run on all samples
am=ReadInput.ReadTextFile(options.configFile, "cfg", samplesToRun,"",options.runBatch)
am.debug=2

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

    submitFiles = []

    #for sampleName in am.ListSampleNames():
    nFilesPerJob = options.nFilesPerJob 
    for sample in am.samples:
        sampleName = sample.sampleName
        print sampleName
        os.system("mkdir -p %s/%s" % (jobName,sampleName))
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
            content += "Should_Transfer_Files = YES\n"
            content += "Output = %i.stdout\n" % nProcJobs
            content += "Error  = %i.stderr\n" % nProcJobs
            content += "Log    = %i.log\n"    % nProcJobs
            content += "Notification = never\n"
            content += "WhenToTransferOutput=On_Exit\n"
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
        cd VHbbAnalysis
        #source /cvmfs/cms.cern.ch/slc6_amd64_gcc481/lcg/root/5.34.22-cms/bin/thisroot.sh
        #source /cvmfs/cms.cern.ch/slc6_amd64_gcc481/lcg/root/5.34.18-cms/bin/thisroot.sh
        echo "successfully set up the enviroment"

        #echo "moving text files to their respective directories"
        #mkdir -p cfg
        #mv samples.txt earlybranches.txt existingbranches.txt newbranches.txt bdtsettings.txt reg1_settings.txt reg2_settings.txt settings.txt systematics.txt scalefactors.txt cfg
        #mkdir -p aux
        #mv TMVARegression_BDTG_ttbar_Nov23.weights.xml  TMVA_13TeV_Dec14_3000_5_H125Sig_0b1b2bWjetsTTbarBkg_Mjj_BDT.weights.xml MuonIso_Z_RunCD_Reco74X_Dec1.json SingleMuonTrigger_Z_RunCD_Reco74X_Dec1.json MuonID_Z_RunCD_Reco74X_Dec1.json CutBasedID_TightWP.json CutBasedID_LooseWP.json aux


        echo "running RunSample.py"
        echo $ORIG_DIR/$4
        python RunSample.py $1 $2 $3 $ORIG_DIR/$4
        echo "all done!" ''' % (os.getcwd(), os.getcwd() )

    runscript = open("%s/condor_runscript.sh" % (jobName) , "w")
    runscript.write(condor_runscript_text)
    runscript.close()

    ## Send job to condor
    print "Submit files created, sending jobs to Condor..."
    for file in submitFiles:
        print("condor_submit %s" % file)
        os.system("condor_submit %s" % file)

        

