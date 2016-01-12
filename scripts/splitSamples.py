import ROOT
import sys

## Create histogram with histograms of a given BDT shape, one histogram for each sample. Used in creation of datacards.
##
## Author: Stephane Cooperstein
##

if (len(sys.argv) != 4):
    print "Takes two arguments:"
    print "python splitSamples.py [input_root_file] [catName] [presel]"
    sys.exit(1)

ROOT.gROOT.SetBatch(True)

ifile = ROOT.TFile(sys.argv[1], "r")
catName = sys.argv[2]
presel = sys.argv[3]
tree = ifile.Get("tree")
#bdtname = "BDT_wMass_Dec14_3000_5"
#bdtname = "BDT_wMass_Dec4"
bdtname = "CMS_vhbb_BDT_Wln_13TeV"
#nBins = 1000
nBins = 20

sampleMap = {} # map sampleNames to list of sampleIndex's

#sampleMap["data_obs"] = [16,17,10,11,20,21,50,51,52,2200,2201,2202,2300,2301,2302,3500,3501,3502,4400,4401,4402,4500,4501,4502,4600,4601,4602,4700,4701,4702,24,25,26,27,28,29,30,31] # dummy filler until we unblind
sampleMap["data_obs"] = [50,51,52]
sampleMap["TT"] = [50,51,52]
sampleMap["s_Top"] = [16,17,10,11,20,21]
sampleMap["WH"] = [-12501]
sampleMap["ZH"] = [-12502]
sampleMap["Wj0b"] = [2200,4400,4500,4600,4700]
sampleMap["Wj1b"] = [2201,4401,4501,4601,4701]
sampleMap["Wj2b"] = [2202,4402,4502,4602,4702]
sampleMap["VVHF"] = [3501,3502]
sampleMap["VVLF"] = [3500]
sampleMap["QCD"]  = [24,25,26,27,28,29,30,31]
sampleMap["Zj0b"] = [2300]
sampleMap["Zj1b"] = [2301]
sampleMap["Zj2b"] = [2302]

ofile = ROOT.TFile("hists_%s.root" % catName, "RECREATE")
for sample in sampleMap:
    cutString = presel + "&&("
    for index in sampleMap[sample]:
            cutString += "(sampleIndex==%i)||" % index
    cutString = cutString[0:len(cutString)-2]
    cutString += ")"
    #print cutString
    hBDT = ROOT.TH1F(sample,sample,nBins,-1,1)
    #tree.Draw("%s>>%s" % (bdtname, sample),"(%s)*weight" % cutString)
    tree.Draw("%s>>%s" % (bdtname, sample),"(%s)*weight*(2.2/1.28)" % cutString) # temp hack to avoid rerunning just to change lumi
    ofile.cd()
    hBDT.Write()

ofile.Close()
    
