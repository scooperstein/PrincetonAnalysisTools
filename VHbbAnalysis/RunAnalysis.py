import ROOT
import sys
import ReadInput
from optparse import OptionParser
import os

parser = OptionParser()
parser.add_option("-c", "--config", dest="configFile", default="vhbb_config.txt",
                  help="specify config file for this analysis"
)
parser.add_option("-b", "--batch", dest="runBatch", default=False,
                  help="run analysis jobs on condor (default=False)"
)
parser.add_option("-n", "--jobName", dest="jobName", default="condor_jobs",
                  help="Specify label for condor jobs. Only to be used when running batch jobs"
)
(options, args) = parser.parse_args()

ROOT.gSystem.Load("AnalysisDict.so")

# reads samples, existing branches and new branches
samplesToRun = [] # if empty run on all samples
am=ReadInput.ReadTextFile(options.configFile, "cfg", samplesToRun)
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
    for sample in am.samples:
        sampleName = sample.sampleName
        print sampleName
        os.system("mkdir -p %s/%s" % (jobName,sampleName))
        nFiles = 0
        for filename in sample.files:
            nFiles += 1
            fname = "%s/%s/job%i.submit" % (jobName, sampleName,nFiles)
            submitFile = open(fname, "w")
            content =  "universe = vanilla\n"
            content += "Executable = condor_runscript.sh\n"
            content += "Arguments = %s %s %s %i\n" % (options.configFile, sampleName, filename, nFiles)
            #content += "Requirements   =  OpSys == 'LINUX' && (Arch =='INTEL' || Arch =='x86_64')\n"
            content += "initialdir = %s/%s\n" % (jobName,sampleName)
            content += "Should_Transfer_Files = YES\n"
            content += "Output = %i.stdout\n" % nFiles
            content += "Error  = %i.stderr\n" % nFiles
            content += "Log    = %i.log\n"    % nFiles
            content += "Notification = never\n"
            content += "WhenToTransferOutput=On_Exit\n"
            content += "transfer_input_files = ../../%s,../../cfg/samples.txt,../../cfg/earlybranches.txt,../../cfg/existingbranches.txt,../../cfg/newbranches.txt,../../cfg/bdtsettings.txt,../../cfg/reg1_settings.txt,../../cfg/reg2_settings.txt,../../cfg/settings.txt,../../aux/new-weights-23Jan.xml,../../aux/TMVA_8TeV_H125Sig_LFHFWjetsNewTTbarVVBkg_newCuts4_BDT.weights.xml,../../RunSample.py,../../../AnalysisDict.so\n" % options.configFile
            content += "Queue = 1\n"
            print content
            submitFile.write(content)
            submitFile.close()
            submitFiles.append(fname)

    # Send job to condor
    print "Submit files created, sending jobs to Condor..."
    for file in submitFiles:
        print("condor_submit %s" % file)
        os.system("condor_submit %s" % file)

        

