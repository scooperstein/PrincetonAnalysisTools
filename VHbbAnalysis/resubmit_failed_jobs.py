## Resubmit analysis jobs which failed to successfully transfer output
##
## Author: Stephane Cooperstein
##

import subprocess
import sys
import os
import ROOT
import argparse
import time

parser = argparse.ArgumentParser(description='Resubmit jobs with bad and/or missing output.')
parser.add_argument('-d', '--dir',    type=str, default="", help='Directory with would-be output to check.')
parser.add_argument('-m', '--missing', default=True,  help="Resubmit when output ROOT files are missing.  (Default:  True)")
parser.add_argument('-e', '--empty',   default=False, help="Resubmit when output ROOT files are empty.  (Default:  False)")
parser.add_argument('-c', '--check',   default=False, help="Just check if output files are present and/or valid.  (Default:  False)")
args = parser.parse_args()



if args.dir=="":
    print "Tell me which job output directory to check!! -d DIRECTORY"
    sys.exit(1)


filesToResubmit = []

for subdir, dirs, files in os.walk(args.dir):
    for file in files:
        if (".submit" in file):
            sample = subdir.split('/')[1]
            jobNum=file.replace(".submit",".root").replace("job","")
            rootfilename = "output_%s_%s" % (sample, jobNum ) 
            if (rootfilename not in files):
                if args.missing:
                    print "root output does not exist",rootfilename
                    filesToResubmit.append(os.path.join(subdir, file) )
                continue

            if args.empty:
                try:
                    rootfile = ROOT.TFile("%s/%s" % (subdir, rootfilename), "r")
                    otree = rootfile.Get("tree")
                    # make sure the proper output ntuple exists in the output root file
                    otree.GetEntries()
                    rootfile.Close()
                except AttributeError:
                    filesToResubmit.append(os.path.join(subdir, file) )
                

for failed_file in filesToResubmit:
    cmd = ["condor_submit", failed_file]
    if not args.check:
        try:
            subprocess.Popen(cmd)
            print "condor_submit %s" % failed_file
            time.sleep(0.20)
        except:
            print "What, what?!"
            raw_input()
    else:
        print "missing: ",failed_file
