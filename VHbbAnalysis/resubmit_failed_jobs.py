## Resubmit analysis jobs which failed to successfully transfer output
##
## Author: Stephane Cooperstein
##

import subprocess
import sys
import os
from ROOT import TFile, TTree

if (len(sys.argv) != 2):
    print "give one argument, the job output directory"
    sys.exit(1)

path = sys.argv[1]

filesToResubmit = []

for subdir, dirs, files in os.walk(path):
    for file in files:
        if (".submit" in file):
            sample = subdir.split('/')[1]
            rootfilename = "output_%s_%s" % (sample, file.replace(".submit",".root").replace("job","") ) 
            #print "rootfilename = ",rootfilename
            if (rootfilename not in files):
                # root output does not exist
                filesToResubmit.append(os.path.join(subdir, file) )
                continue
            rootfile = TFile("%s/%s" % (subdir, rootfilename), "r")
            otree = rootfile.Get("tree")
            try:
                # make sure the proper output ntuple exists in the output root file
                otree.GetEntries()
            except AttributeError:
                "%s/%s exists but the output tree is invalid" % (subdir, rootfilename)
                filesToResubmit.append(os.path.join(subdir, file) )
            rootfile.Close()
                
print "resubmitting %i failed jobs" % len(filesToResubmit)
for failed_file in filesToResubmit:
    cmd = ["condor_submit", failed_file]
    #subprocess.Popen(cmd)
    print "condor_submit %s" % failed_file
            
