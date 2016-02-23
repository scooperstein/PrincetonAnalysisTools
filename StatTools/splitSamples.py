import ROOT
import sys
from math import sqrt
import numpy
import argparse

## Create histogram with histograms of a given BDT shape, one histogram for each sample. Used in creation of datacards.
##
## Author: Stephane Cooperstein
##

parser = argparse.ArgumentParser("Create histograms for datacards")
parser.add_argument('-i', '--inputfile', type=str, default="", help="The input root file (ntuple)")
parser.add_argument('-c', '--catName', type=str, default="Cat", help="The category label for datacard")
parser.add_argument('-p', '--preselection', type=str, default="", help="The category selection cuts on top of ntuple")
parser.add_argument('-s', '--systematics', type=str, default="", help="The systematics config file")
parser.add_argument('-w', '--weights', type=str, default="", help="Comma-separated list of weights to apply on top of nominal 'weight'")
args = parser.parse_args()
print args

#if (len(sys.argv) != 4 and len(sys.argv) != 5):
#    print "Takes three or four arguments:"
#    print "python splitSamples.py [input_root_file] [catName] [presel]"
#    print "python splitSamples.py [input_root_file] [catName] [presel] [sys.txt] (for shape systematics besides bin-by-bin)"
#    sys.exit(1)

ROOT.gROOT.SetBatch(True)

ifile = ROOT.TFile(args.inputfile, "r")
catName = args.catName
presel = args.preselection

systematics = {}
if (args.systematics != ""):
    # prepare shape systematics
    sysfile = open(args.systematics)
    for line in sysfile:
        if (line[0] == '#'): continue
        line = line.strip()
        params = line.split(' ')
        paramsToKeep = []
        for param in params:
            if (param != ''):
                paramsToKeep.append(param)
        print paramsToKeep
        if (paramsToKeep[1] != "shape"): continue
        sysSamples = paramsToKeep[3].split(',')
        # map systematic name to the corresponding weight in the ntuple and possibly an extension if it has a different BDT/Mjj output than the nominal
        # In that case, the BDT name used is [bdtname]_[paramsToKeep[5]]UpDown
        if (len(paramsToKeep) > 5):
            systematics[paramsToKeep[0]] = (paramsToKeep[4], paramsToKeep[5], sysSamples)
        else:
            systematics[paramsToKeep[0]] = (paramsToKeep[4], "", sysSamples)
    
print systematics

weights = [""]
if (args.weights != ""):
    set_of_weights = args.weights.split(',')
    weights.append(set_of_weights)
print weights
weight_string = "weight"
for weight in set_of_weights:
    weight_string += "*" + weight
print weight_string

tree = ifile.Get("tree")
#bdtname = "BDT_wMass_Dec14_3000_5"
#bdtname = "BDT_wMass_Dec4"
bdtname = "CMS_vhbb_BDT_Wln_13TeV"
#nBins = 1000
nBins = 20
tolerance = 0.5 # dB/B tolerance for bin-by-bin stat. uncertainties

sampleMap = {} # map sampleNames to list of sampleIndex's
sampleMapAltModel = {} # alternate MC samples for model shape systematics

#sampleMap["data_obs"] = [16,17,10,11,20,21,50,51,52,2200,2201,2202,2300,2301,2302,3500,3501,3502,4400,4401,4402,4500,4501,4502,4600,4601,4602,4700,4701,4702,24,25,26,27,28,29,30,31] # dummy filler until we unblind
sampleMap["data_obs"] = [50,51,52]
#sampleMap["TT"] = [50,51,52]
sampleMap["TT"] = [12]
sampleMapAltModel["TT"] = [13]
sampleMap["s_Top"] = [16,17,10,11,20,21]
sampleMap["WH"] = [-12501]
sampleMap["ZH"] = [-12502]
sampleMap["Wj0b"] = [2200,4400,4500,4600,4700]
#sampleMapAltModel["Wj0b"] = [2200]
sampleMap["Wj1b"] = [2201,4401,4501,4601,4701]
#sampleMapAltModel["Wj1b"] = [2201]
sampleMap["Wj2b"] = [2202,4402,4502,4602,4702]
#sampleMapAltModel["Wj2b"] = [2202]
sampleMap["VVHF"] = [3501,3502]
sampleMap["VVLF"] = [3500]
sampleMap["QCD"]  = [24,25,26,27,28,29,30,31]
#sampleMap["Zj0b"] = [2300]
#sampleMap["Zj1b"] = [2301]
#sampleMap["Zj2b"] = [2302]


def makeCutString(sample, sMap):
    cutString = presel + "&&("
    for index in sMap[sample]:
            cutString += "(sampleIndex==%i)||" % index
    cutString = cutString[0:len(cutString)-2]
    cutString += ")"
    return cutString

# first rebin the histogram so that the first and last bins are not empty and have less than 35% stat. uncertainty
hBkg = ROOT.TH1F("hBkg","hBkg",nBins,-1,1)
hSig = ROOT.TH1F("hSig","hSig",nBins,-1,1)
#tree.Draw("%s>>hBkg" % bdtname,"((%s)&&sampleIndex>0)*weight*(2.2/1.28)" % presel)
tree.Draw("%s>>hBkg" % bdtname,"((%s)&&sampleIndex>0&&sampleIndex!=13)*bTagWeightNew*%s" % (presel,weight_string))
tree.Draw("%s>>hSig" % bdtname,"((%s)&&sampleIndex==-12501)*bTagWeightNew*%s" % (presel,weight_string))
binBoundaries = numpy.zeros(nBins+1,dtype=float)
binBoundaries[0] = -1.0
binBoundaries[nBins] = 1.0
foundLowBinEdge = False
foundHighBinEdge = False
for ibin in range(2,nBins):
    if not foundLowBinEdge:
        B_low = hBkg.Integral(1,ibin)
        #S_low = hSig.Integral(1,ibin)
        #print "edge = ",hBkg.GetBinLowEdge(ibin+1)
        #print "B_low = ",B_low
        #if (B_low > 0):
        #    print "1./sqrt(B_low) = ",1./sqrt(B_low)
        if (B_low > 0 and 1./sqrt(B_low) < 0.35):
            binBoundaries[1] = hBkg.GetBinLowEdge(ibin+1)
            foundLowBinEdge = True
    if not foundHighBinEdge:
        B_high = hBkg.Integral(nBins-ibin+1, nBins) 
        #S_high = hSig.Integral(nBins-ibin+1, nBins) 
        if (B_high > 0 and 1./sqrt(B_high) < 0.35):
            binBoundaries[nBins-1] = hBkg.GetBinLowEdge(nBins-ibin+1)
            foundHighBinEdge = True
# split the middle bins equidistantly
for i in range(2,nBins-1):
    binBoundaries[i] = binBoundaries[1] + (i-1)*((binBoundaries[nBins-1] - binBoundaries[1])/(nBins-2)) 

hBkg = hBkg.Rebin(nBins, "", binBoundaries)
print binBoundaries

ofile = ROOT.TFile("hists_%s.root" % catName, "RECREATE")
otextfile = open("binStats_%s.txt" % catName, "w")
hBkg.Write()
for sample in sampleMap:
    #cutString = presel + "&&("
    #for index in sampleMap[sample]:
    #        cutString += "(sampleIndex==%i)||" % index
    #cutString = cutString[0:len(cutString)-2]
    #cutString += ")"
    #print cutString
    cutString = makeCutString(sample, sampleMap)
    print cutString
    #hBDT = ROOT.TH1F(sample,sample,nBins,-1,1)
    hBDT = ROOT.TH1F(sample,sample,nBins,binBoundaries)
    #tree.Draw("%s>>%s" % (bdtname, sample),"(%s)*weight" % cutString)
    #tree.Draw("%s>>%s" % (bdtname, sample),"(%s)*weight*(2.2/1.28)" % cutString) # temp hack to avoid rerunning just to change lumi
    tree.Draw("%s>>%s" % (bdtname, sample),"((%s)&&Pass_nominal)*bTagWeightNew*%s" % (cutString,weight_string)) # temp hack to avoid rerunning just to change lumi
    #hBDT = hBDT.Rebin(nBins, "", binBoundaries)
    # Add bin-by-bin stat. uncertainties
    if (sample not in ["WH","ZH","data_obs"]):
        for ibin in range(1, hBDT.GetNbinsX()):
            B = hBDT.GetBinContent(ibin)
            B_err = hBDT.GetBinError(ibin)
            #hBDT.GetXaxis().SetRange(ibin,ibin)
            #B_eff = hBDT.GetEffectiveEntries() # effective statistical number of entries in bin considering weights
            if ( B > 0 and ( ( B >=1 and (B_err/sqrt(B)) > tolerance) or (B < 1 and B_err/B > tolerance) ) ):
                #hBinStat = ROOT.TH1F("CMS_vhbb_stat%s_%s_bin%i_13TeV" % (sample,catName,ibin),"CMS_vhbb_stat%s_%s_bin%i_13TeV" % (sample,catName,ibin),nBins,-1,1)
                hBinStatUp = hBDT.Clone()
                hBinStatDown = hBDT.Clone()
                hBinStatUp.SetName("%s_CMS_vhbb_stat%s_%s_bin%i_13TeVUp" % (sample,sample,catName,ibin))
                hBinStatDown.SetName("%s_CMS_vhbb_stat%s_%s_bin%i_13TeVDown" % (sample,sample,catName,ibin))
                hBinStatUp.SetBinContent(ibin, B + sqrt(B))
                hBinStatDown.SetBinContent(ibin, max(B - sqrt(B),0.000001))
                otextfile.write("CMS_vhbb_stat%s_%s_bin%i_13TeV\n" % (sample,catName,ibin))
                ofile.cd()
                hBinStatUp.Write()
                hBinStatDown.Write()
    # other shape systematics
    for syst in systematics:
        sysWeight, sysName, sysSamples = systematics[syst]
        if sample not in sysSamples: continue
        hBDTSystUp = ROOT.TH1F("%s_%sUp" % (sample,syst), "%s_%sUp" % (sample,syst),nBins,binBoundaries)
        hBDTSystDown = ROOT.TH1F("%s_%sDown" % (sample,syst), "%s_%sDown" % (sample,syst),nBins,binBoundaries)
        sysBDTNameUp = bdtname
        sysBDTNameDown = bdtname
        passSys = "Pass_nominal"
        if (sysName != ""):
            sysBDTNameUp += "_%sUp" % sysName
            sysBDTNameDown += "_%sDown" % sysName
            pasSys = "Pass_%s" % sysName
        if (sysWeight != "1.0"):
            #tree.Draw("%s>>%s_%sUp" % (sysBDTNameUp, sample, syst),"(%s)*weight*(2.2/1.28)*(%sUp)" % (cutString,sysWeight)) # temp hack to avoid rerunning just to change lumi 
            tree.Draw("%s>>%s_%sUp" % (sysBDTNameUp, sample, syst),"((%s)&&%s)*%s*(%sUp)" % (cutString,passSys,weight_string,sysWeight)) # temp hack to avoid rerunning just to change lumi 
            #tree.Draw("%s>>%s_%sDown" % (sysBDTNameDown, sample, syst),"(%s)*weight*(2.2/1.28)*(%sDown)" % (cutString,sysWeight)) # temp hack to avoid rerunning just to change lumi 
            tree.Draw("%s>>%s_%sDown" % (sysBDTNameDown, sample, syst),"((%s)&&%s)*%s*(%sDown)" % (cutString,passSys,weight_string,sysWeight)) # temp hack to avoid rerunning just to change lumi 
        else:
            #tree.Draw("%s>>%s_%sUp" % (sysBDTNameUp, sample, syst),"(%s)*weight*(2.2/1.28)" % (cutString)) # temp hack to avoid rerunning just to change lumi 
            tree.Draw("%s>>%s_%sUp" % (sysBDTNameUp, sample, syst),"((%s)&&%s)*%s*bTagWeightNew" % (cutString,passSys,weight_string)) # temp hack to avoid rerunning just to change lumi 
            #tree.Draw("%s>>%s_%sDown" % (sysBDTNameDown, sample, syst),"(%s)*weight*(2.2/1.28)" % (cutString)) # temp hack to avoid rerunning just to change lumi 
            tree.Draw("%s>>%s_%sDown" % (sysBDTNameDown, sample, syst),"((%s)&&%s)*%s*bTagWeightNew" % (cutString,passSys,weight_string)) # temp hack to avoid rerunning just to change lumi 
            if (syst.find("Model") != -1):
                # shape from different sample (amc@NLO vs. POWHEG, for example)
                cutStringAltModel = makeCutString(sample,sampleMapAltModel)
                hBDTAltModel = ROOT.TH1F("%s_%sAltModel" % (sample,syst), "%s_%sAltModel" % (sample,syst),nBins,binBoundaries)
                tree.Draw("%s>>%s_%sAltModel" % (sysBDTNameUp, sample, syst),"((%s)&&%s)*%s*bTagWeightNew" % (cutStringAltModel,passSys,weight_string))
                # hnom + (hnom - halt) = 2*hnom - halt
                hBDTSystUp.Scale(2.0)
                hBDTSystUp.Add(hBDTAltModel, -1.0)
                # hnom - (hnom - halt) = halt
                hBDTSystDown = hBDTAltModel
                hBDTSystDown.SetName("%s_%sDown" % (sample, syst))
            
        ofile.cd()
        hBDTSystUp.Write()
        hBDTSystDown.Write()
    ofile.cd()
    hBDT.Write()

ofile.Close()
otextfile.close()
   
