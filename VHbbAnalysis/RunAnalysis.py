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
parser.add_option("-o", "--outputDir",    dest="outputDir",    default="", type=str, help="Output directory for the jobs")
parser.add_option("-s","--sample", dest="sample", default="", type=str, help="Run on only a specific sample (can be comma-separated list)")
parser.add_option("-d","--doData", dest="doData", default=-1, type=int, help="If -1 run all samples, if 0 run only MC, if 1 run only data")
parser.add_option("--site","--site", dest="site", default="FNAL", type=str, help="If running on lxplus include option --site CERN, otherwise assumes you are running on FNAL")
parser.add_option("--useSGE","--useSGE", dest="useSGE", default=0, type=int, help="If 0 (default) use condor, if 1 use SGE job submission")
parser.add_option("--doSkim","--doSkim", dest="doSkim", default=0, type=int, help="If 0 (default) run analysis jobs, if 1 run skimming jobs")
parser.add_option("--submitJobs","--submitJobs", dest="submitJobs", default=1, type=int, help="If 1 (default) submit jobs to batch queue, if 0 then only create submission files")
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
site = options.site
useSGE = options.useSGE
doSkim = options.doSkim
submitJobs = options.submitJobs

if (options.runBatch == False):
    print "Running locally over all samples"

    if(am.debug>100):
        am.PrintBranches()

    # loop over all the samples
    # FIXME - need to add the possibility of doing a small portion of files
    am.Loop()

else:
    print "Running analysis jobs on batch queue. Site: %s" % site
    if useSGE:
        print "Using SGE batch system..."
    else:
        print "Using Condor batch system..."
    jobName = options.jobName

    if site not in ("FNAL", "CERN", "DESY"):
        print "unknown site: ",site,"exiting..."
        sys.exit(1)

    if site == "FNAL":  ## let's do a little autodetect for desy and cern..
        import socket
        hostname = socket.gethostname()
        if hostname.endswith('.desy.de'):
            print 'detected site: DESY'
            site = 'DESY'
        if hostname.endswith('.cern.ch'):
            print 'detected site: CERN'
            site = 'CERN'

    if site == "FNAL":
        output_dir = "/eos/uscms/store/user/sbc01/VHbbAnalysisNtuples" # parent directory in user's EOS space to store all output files
    elif site == "CERN":
        output_dir = "/eos/cms/store/user/scoopers/VHbbAnalysisNtuples"
    elif site == "DESY":
        output_dir = "/nfs/dust/cms/user/%s/VHbbAnalysisNtuples" % os.getlogin()
    if doSkim:
        output_dir = output_dir.replace("VHbbAnalysisNtuples","SkimmedAnalysisNtuples")

    if options.outputDir:
        print "changing output_dir to",options.outputDir
        output_dir = options.outputDir


    #if os.path.exists(jobName):
    #    print "Attempting to create jobs under jobName %s, which already exists!" % jobName
    #    sys.exit(0)
    os.system("mkdir -p %s" % jobName)
    os.system("mkdir -p %s/%s" % (output_dir,jobName))

    submitFiles = []

    #for sampleName in am.ListSampleNames():
    inputs_to_transfer = [
        options.configFile,
        'RunSkim.py',
        'RunSample.py',
        '../AnalysisDict.so',
        '../env.sh',
        'cfg',
        'aux',
        '../python/ReadInput.py',
        '../AnalysisDict_rdict.pcm',
        '../HelperClasses',
        '../AnalysisManager.h',
        '../plugins',
    ]
    # copy inputs to jobdir
    if site == "DESY":
        os.system('cp -r %s %s' % (' '.join(inputs_to_transfer), jobName))
        # prepend paths with '../' for the condor jobs
        inputs_to_transfer = list('../'+p if not p.startswith('../') else p for p in inputs_to_transfer)
        inputs_to_transfer.remove('../python/ReadInput.py')
        inputs_to_transfer.append('../ReadInput.py')
    else:
        inputs_to_transfer = list('../../'+p for p in inputs_to_transfer)

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
        os.system("mkdir -p %s/%s/%s" % (output_dir,jobName,sampleName))
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
            if not useSGE:
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
                content += "transfer_input_files = %s\n" % ','.join(inputs_to_transfer)
                content += "Queue  1\n"
            else:
                fname = "%s/%s/job%iSubmit.sh" % (jobName, sampleName,nProcJobs)
                submitFile = open(fname, "w")
                content = "source %s/%s/condor_runscript.sh %s %s %s output_%s_%i.root\n" % (os.getcwd(),jobName,options.configFile, sampleName, filesToRun,sampleName, nProcJobs)
            #print content
            submitFile.write(content)
            submitFile.close()
            submitFiles.append(fname)

    xrdcp_string = "xrdcp -f $ORIG_DIR/$4 %s/%s/$2" %(output_dir.replace("/eos/uscms","root://cmseos.fnal.gov/"),jobName)
    if site == "CERN":
        xrdcp_string = "xrdcp -f $ORIG_DIR/$4 root://eoscms.cern.ch:/%s/%s/$2" %(output_dir,jobName)
    if site == "DESY":
        local_path = "%s/%s" % (output_dir,jobName)
        xrdcp_string = "mkdir -p %s/$2; cp -vf $4 %s/$2" % (local_path, local_path)
    copy_string = "" # not sure how to automatically transfer files to SGE job nodes, for now do it manually
    runfile = "RunSample.py"
    if doSkim:
        runfile = "RunSkim.py"
    if useSGE:
        copy_string = ''' cp -r %s/cfg .
        cp -r %s/aux .
        cp %s/%s .
        cp %s/../AnalysisDict.so .
        cp -r %s/*.txt . ''' % (os.getcwd(),os.getcwd(),os.getcwd(),runfile,os.getcwd(),os.getcwd())

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
        #cd %s/..
        cd /cvmfs/cms.cern.ch/slc6_amd64_gcc493/cms/cmssw-patch/CMSSW_7_6_3_patch2/src
        source /cvmfs/cms.cern.ch/cmsset_default.sh
        eval `scramv1 runtime -sh`
        #source env.sh
        #source %s/../env.sh
        #source env.sh

        echo "successfully set up the enviroment"

        cd $ORIG_DIR
        pwd
        ls

        #echo "copying over necessary files"
        #cp -r %s/cfg .
        #cp -r %s/aux .
        #cp %s/RunSample.py .
        #cp %s/../AnalysisDict.so .
        #cp -r %s/*.txt .
        %s

        echo "running %s"
        echo $ORIG_DIR/$4
        python %s $1 $2 $3 $ORIG_DIR/$4
        echo "done running, now copying output to EOS"
        ###xrdcp -f $ORIG_DIR/$4 root://cmseos.fnal.gov//store/user/sbc01/VHbbAnalysisNtuples/%s/$2
        %s
        rm $ORIG_DIR/$4
        echo "all done!" ''' % (os.getcwd(), os.getcwd(), os.getcwd(), os.getcwd(), os.getcwd(), os.getcwd(), os.getcwd(), copy_string, runfile, runfile, jobName, xrdcp_string )

    condor_runscript_text_desy = '''

        export PATH=/afs/desy.de/common/passwd:/usr/local/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/cvmfs/grid.cern.ch/emi3ui-latest/bin:/cvmfs/grid.cern.ch/emi3ui-latest/sbin:/cvmfs/grid.cern.ch/emi3ui-latest/usr/bin:/cvmfs/grid.cern.ch/emi3ui-latest/usr/sbin:$PATH
        echo "echo PATH:"
        echo $PATH
        echo "arguments: " $1 $2 $3 $4
        echo "username and group"
        id -n -u
        id -n -g

        echo "creating tempdir and copy"
        tmp_dir=$(mktemp -d)
        cp -r %s $tmp_dir

        echo "setting up the environment"
        cd /cvmfs/cms.cern.ch/slc6_amd64_gcc493/cms/cmssw-patch/CMSSW_7_6_3_patch2/src
        source /cvmfs/cms.cern.ch/cmsset_default.sh
        eval `scramv1 runtime -sh`
        echo "echo PATH:"
        echo $PATH

        echo "changing to tempdir"
        cd $tmp_dir
        pwd
        ls

        echo "running %s"
        echo $4
        python %s $1 $2 $3 $4
        echo "done running, now copying output to EOS"

        echo "copying output (%s)"
        %s

        echo "delete tmp dir"
        cd $TMP
        rm -r $tmp_dir

        echo "all done!"
    ''' % (' '.join(inputs_to_transfer), runfile, runfile, xrdcp_string, xrdcp_string)

    runscript = open("%s/condor_runscript.sh" % (jobName) , "w")
    runscript.write(condor_runscript_text_desy if site=='DESY' else condor_runscript_text)
    runscript.close()
    os.system("chmod 755 %s/condor_runscript.sh" % (jobName))

    if submitJobs:
        if not useSGE:
            ## Send job to condor
            print "Submit files created, sending jobs to Condor..."
            for file in submitFiles:
                print("condor_submit %s" % file)
                os.system("condor_submit %s" % file)
        else:
            ## Send jobs with SGE
            print "Submit files created, sending jobs to SGE batch..."
            for file in submitFiles:
                print 'bsub -R "pool>30000" -q 1nh -J job1 < %s' %file
                os.system('bsub -R "pool>30000" -q 1nh -J job1 < %s' %file)
    else:
        print "Submit files created, but will not submit jobs to batch queue"



