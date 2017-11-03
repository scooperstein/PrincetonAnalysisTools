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
parser.add_argument('-o', '--odir',    type=str, default="", help='Directory with would-be output to check, if none specified assume it is the same as the submission.')
parser.add_argument('-d', '--dir',    type=str, default="", help='Directory with submit files to check.')
parser.add_argument('-m', '--missing', default=True,  help="Resubmit when output ROOT files are missing.  (Default:  True)")
parser.add_argument('-e', '--empty',   default=False, help="Resubmit when output ROOT files are empty.  (Default:  False)")
parser.add_argument('-c', '--check',   default=False, help="Just check if output files are present and/or valid.  (Default:  False)")
args = parser.parse_args()



if args.dir=="":
    print "Tell me which job output directory to check!! -d DIRECTORY"
    sys.exit(1)

if args.odir=="":
    args.odir = args.dir

filesToResubmit = []

rootFiles = []
for subdir, dirs, files in os.walk(args.odir):
    for file in files:
        if (".root" in file):
             tmp = file.split('/')
             filename = tmp[len(tmp)-1]
             rootFiles.append(filename)

for subdir, dirs, files in os.walk(args.dir):
    for file in files:
        #if subdir.find("QCD")!=-1 or subdir.find("TT_DiLep")!=-1 or subdir.find("TT_SingleLep")!=-1: continue
        #if subdir.find("QCD")!=-1: continue
        if (".submit" in file):
            dirpaths = subdir.split('/')
            sample = dirpaths[len(dirpaths)-1]
            jobNum=file.replace(".submit",".root").replace("job","")
            #jobNum = file[file.find("output_%s" % sample)+7 + len(sample):file.find(".root")]
            rootfilename = "output_%s_%s" % (sample, jobNum ) 
            if (rootfilename not in rootFiles):
                if args.missing:
                    print "root output does not exist",rootfilename
                    filesToResubmit.append(os.path.join(subdir, file) )
                continue

            if args.empty:
                try:
                    #rootfile = ROOT.TFile("%s/%s" % (subdir, rootfilename), "r")
                    rootfile = ROOT.TFile("%s/%s/%s" % (args.odir,sample,rootfilename),"r")
                    otree = rootfile.Get("tree")
                    # make sure the proper output ntuple exists in the output root file
                    nentries = otree.GetEntries()
                    #print nentries
                    #if (nentries == 0):
                    #    print rootfilename + " has 0 entries!"
                    #    filesToResubmit.append(os.path.join(subdir, file) )
                    rootfile.Close()
                except AttributeError:
                    print "caught"
                    filesToResubmit.append(os.path.join(subdir, file) )
                
print "resubmitting %i failed jobs" % len(filesToResubmit)
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
