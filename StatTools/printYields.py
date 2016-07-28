import ROOT
#import sys
import argparse

parser = argparse.ArgumentParser("Create datacards from histograms")
parser.add_argument('-c', '--channelName', type=str, default="Chan", help="The label for the channel in the datacard")
parser.add_argument('-i', '--ihistfile', type=str, default="hists.root", help="The input root file with the histograms")
parser.add_argument('-s', '--systematics', type=str, default="", help="The systematics config. file (default no systematics)")
parser.add_argument('-b', '--binstats', type=str, default="", help="Text file listing all the individual bin. stat. uncertainties to include")
parser.add_argument('-r', '--rateParams', type=str, default="", help="Comma-separated list of samples for which to include freely floating scale factors") 
parser.add_argument('-o', '--outputfilename', type=str, default="", help="Name of the output datacard")
parser.add_argument('-vv','--doVV', type=bool, default=False, help="If true do VV analysis (default False)")
args = parser.parse_args()


rates = "rate                                                        " # yields in datacard format
dc_string = "" # write datacard

#ifile = ROOT.TFile("hists.root")

#cats = ["WenLowPt","WenHighPt","WmnLowPt","WmnMidPt","WmnHighPt"]
#cats = ["WenLowPt_0","WenLowPt_1","WenHighPt_0","WenHighPt_1","WmnLowPt_0","WmnLowPt_1","WmnMidPt_0","WmnMidPt_1","WmnHighPt_0","WmnHighPt_1"]
#cat_labels = ["ch1_Wenu","ch1_Wenu3","ch2_Wmunu","ch2_Wmunu2","ch2_Wmunu3"]
#cat_labels = ["ch1_WenuL","ch1_WenuH","ch1_Wenu3L","ch1_Wenu3H","ch2_WmunuL","ch2_WmunuH","ch2_Wmunu2L","ch2_Wmunu2H","ch2_Wmunu3L","ch2_Wmunu3H"]
#samples = ["ZH","WH","s_Top","Zj1b","TT","Zj0b","Wj0b","Wj1b","Wj2b","Zj2b"]
#samples = ["ZH","WH","s_Top","Zj1b","TT","Zj0b","Wj0b","Wj1b","Wj2b","Zj2b"]
#samples = ["ZH","WH","s_Top","TT","Wj0b","Wj1b","Wj2b","VVHF","VVLF","QCD","Zj0b","Zj1b","Zj2b"]
#samples = ["ZH","WH","s_Top","TT","Wj0b","Wj1b","Wj2b","VVHF","VVLF","Zj0b","Zj1b","Zj2b"]
samples = ["ZH","WH","s_Top","TT","Wj0b","Wj1b","Wj2b","VVHF","VVLF","Zj0b","Zj1b","Zj2b","QCD"]
#samples = ["ZH","WH","Bkg"]
#samples = ["ZH","WH","s_Top","TT","Wj0b","Wj1b","Wj2b","VVHF","VVLF"]
#samples = ["ZH","WH","s_Top","TT","Wj0b","Wj1b","Wj2b"]
#samples = ["WH","TT","s_Top"]
#samples = ["WH","TT"]
#cats = ["WmnLowPt","WmnMidPt", "WmnHighPt"]
#cat_labels = ["ch2_Wmunu", "ch2_Wmunu2","ch2_Wmunu3"]
#cats = ["WenLowPt", "WenHighPt"]
#cat_labels = ["ch1_Wenu","ch1_Wenu3"]
#if (len(sys.argv) > 2):
#    s = sys.argv[2].strip(',')
#    cats = s.split(',')
#    print cats
#    cat_labels = list(cats)

if args.doVV:
    samples = ["VVHF","VVLF","ZH","WH","s_Top","TT","Wj0b","Wj1b","Wj2b","Zj0b","Zj1b","Zj2b","QCD"]

# It probably makes sense to just run this script once per channel and
# then combine the datacards with combineDatacard.py, but the code is set
# up for multiple channels at a time if ever necessary.
cats = []
cat_labels = []
channelNames = args.channelName.split(',')
print channelNames
print args.channelName
for channelName in channelNames:
    cats.append(channelName)
    cat_labels.append(channelName)

print "cat_labels = "
print cat_labels
fake_obs_line = "observation  "

cats_to_remove = [] # categories with zero yields for all samples, remove from datacard
for i in range(len(cats)):
    print "Calculating Yields for Category: %s: %s" % (cats[i],cat_labels[i])
    nBkgTot = 0.0;
    nSigTot = 0.0;
    ifile = ROOT.TFile("hists_%s.root" % cats[i])
    #ifile = ROOT.TFile(args.ihistfile, "r")
    zeroYield = True # don't write to datacard if all samples in cat are zero
    cat_rates = ""
    print cat_labels[i]
    nData = ifile.Get("BDT_%s_data_obs" % cat_labels[i]).Integral()
    cat_fake_obs_line = "   %.4f" % nData
    for sample in samples: 
        print sample
        #nyield = ifile.Get("%s/%s" % (cat,sample)).Integral()
        nyield = ifile.Get("BDT_%s_%s" % (cat_labels[i],sample)).Integral()
        if (nyield > 0): zeroYield = False
        print "%s = %.4f" % (sample, nyield )
        if (sample != "WH" and sample != "ZH"): nBkgTot += nyield
        else: nSigTot += nyield
        # for readable aligment
        if (nyield < 10):
            cat_rates += "%.4f      " % nyield
        elif (nyield < 100):
            cat_rates += "%.4f     " % nyield
        elif (nyield < 1000):
            cat_rates += "%.4f    " % nyield
        else:
            cat_rates += "%.4f   " % nyield
    if not zeroYield: 
        rates += cat_rates
        fake_obs_line += cat_fake_obs_line
        #fake_obs_line += "   %.4f" % nBkgTot
    else: cats_to_remove.append(cats[i])
    print "Total Signal Yield: %.4f" % nSigTot
    print "Total Background Yield: %.4f" % nBkgTot

for cat in cats_to_remove:
    cat_labels.remove(cat_labels[cats.index(cat)])
    cats.remove(cat)

print rates

systematics = ""

scale_factors = args.rateParams.split(',')
for scale_factor in scale_factors:
    for sample in samples:
        if (scale_factor.find(sample) != -1):
            systematics += "SF_%s  rateParam  %s %s  1\n" % (scale_factor,cats[0], sample) 

nSys = 0

if (args.systematics != ""):
    sys_file = open(args.systematics,"r")
    for line in sys_file:
        if (line[0] == '#'): continue
        line = line.strip()
        #print line
        params = line.split(' ')
        # remove extra spaces, there's probably a smarter way to do this
        paramsToKeep = []
        for param in params:
            if (param != ''):
                paramsToKeep.append(param)
        params = paramsToKeep
        name = params[0]
        sysType = params[1]
        sysSamples = params[3].split(',')
        vals = []
        if (params[2].find(',') == -1):
            for sample in sysSamples:
                vals.append(params[2])
        else:
            # separate values for different processes
            vals = params[2].split(',')
        sysLine = name
        
        for i in range(48 - len(name)):
            sysLine += " "
        sysLine += sysType
        for i in range(12 - len(sysType)):
            sysLine += " "
        nSys += 1
        for i in range(len(cats)):
            j = 0
            for sample in samples:
                if sample in sysSamples:
                    sysLine += "%s        " % vals[j]
                    j += 1
                else:
                    sysLine += "-           "  
        sysLine += "\n"
        systematics += sysLine

if (args.binstats != ""):
    for cat in cat_labels:
        try:
            #binStats_file = open(args.binstats, "r")
            binStats_file = open("binStats_%s.txt" % cat,"r")
        except FileOpenError:
            print "Couldn't open binstats file %s" % args.binstats
            continue
        for line in binStats_file:
            if (line[0] == '#'): continue
            line = line.strip()
            sysLine = line
            for i in range(60 - len(line)):
                sysLine += " "
            sysLine += "shape      "
            matchesSample = False
            for cat_i in cat_labels:
                for sample in samples:
                    if (cat_i == cat and sysLine.find(sample) != -1):
                        sysLine += "1.0        "
                        matchesSample = True
                    else:
                        sysLine += "-           "
            sysLine += "\n"
            if matchesSample:
                nSys += 1
                systematics += sysLine

dc_string += "imax %i number of bins\n" % len(cats)
dc_string += "jmax %i number of processes minus 1\n" % (len(samples) - 1)
dc_string += "kmax %i number of nuisance parameters\n" % nSys
dc_string += "----------------------------------------------------------------------------------------------------------------------------------\n"
for i in range(len(cats)):
    dc_string += "shapes *           %s    hists_%s.root BDT_$CHANNEL_$PROCESS BDT_$CHANNEL_$PROCESS_$SYSTEMATIC\n" % (cat_labels[i],cat_labels[i])
    #dc_string += "shapes *           %s    %s BDT_$CHANNEL_$PROCESS BDT_$CHANNEL_$PROCESS_$SYSTEMATIC\n" % (cat_labels[i],args.ihistfile)
dc_string += "----------------------------------------------------------------------------------------------------------------------------------\n"
dc_string += "bin          "
for label in cat_labels:
    dc_string += "%s    " % label
dc_string += "\n"
#dc_string += "observation  "
#for i in range(len(cats)):
#    dc_string += "0.0000      "
dc_string += "%s\n" % fake_obs_line
dc_string += "\n"
dc_string += "bin                                                         "
for label in cat_labels:
    for sample in samples:
        dc_string += label
        for i in range(12 - len(label)):
            dc_string += " "
dc_string += "\n"
dc_string += "process                                                     "
for label in cat_labels:
    for sample in samples:
        dc_string += sample
        for i in range(12 - len(sample)):
            dc_string += " "
dc_string += "\n"
dc_string += "process                                                     "
for label in cat_labels:
    for i in range(len(samples)):
        dc_string += str(i-1)
        if (i < 10): dc_string += "           "
        else: dc_string += "          "
dc_string += "\n"
dc_string += rates
# dummy systematic
dc_string += "\n"
#dc_string += "lumi    lnN      "
#for i in range(len(cat_labels)):
#    for j in range(len(samples)):
#        dc_string += "1.0001    "
dc_string += "\n"

dc_string += systematics

if (args.outputfilename != ""):
    ofile = open(args.outputfilename ,"write")
else: ofile = open("vhbb_%s_13TeV.txt" % args.channelName ,"write")
ofile.write(dc_string)
ofile.close()
print dc_string
ifile.Close()
