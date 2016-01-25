import ROOT
import sys
from math import sqrt
import numpy

## Create histogram with histograms of a given BDT shape, one histogram for each sample. Used in creation of datacards.
##
## Author: Stephane Cooperstein
##

if (len(sys.argv) != 4 and len(sys.argv) != 5):
    print "Takes three or four arguments:"
    print "python splitSamples.py [input_root_file] [catName] [presel]"
    print "python splitSamples.py [input_root_file] [catName] [presel] [sys.txt] (for shape systematics besides bin-by-bin)"
    sys.exit(1)

ROOT.gROOT.SetBatch(True)

ifile = ROOT.TFile(sys.argv[1], "r")
catName = sys.argv[2]
presel = sys.argv[3]

systematics = {}
if (len(sys.argv) == 5):
    # prepare shape systematics
    sysfile = open(sys.argv[4])
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

tree = ifile.Get("tree")
#bdtname = "BDT_wMass_Dec14_3000_5"
#bdtname = "BDT_wMass_Dec4"
bdtname = "CMS_vhbb_BDT_Wln_13TeV"
#nBins = 1000
nBins = 20
tolerance = 0.25 # dB/B tolerance for bin-by-bin stat. uncertainties

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

# first rebin the histogram so that the first and last bins are not empty and have less than 35% stat. uncertainty
hBkg = ROOT.TH1F("hBkg","hBkg",nBins,-1,1)
#tree.Draw("%s>>hBkg" % bdtname,"((%s)&&sampleIndex>0)*weight*(2.2/1.28)" % presel)
tree.Draw("%s>>hBkg" % bdtname,"((%s)&&sampleIndex>0)*weight*1000" % presel)
binBoundaries = numpy.zeros(nBins+1,dtype=float)
binBoundaries[0] = -1.0
binBoundaries[nBins] = 1.0
foundLowBinEdge = False
foundHighBinEdge = False
for ibin in range(2,nBins):
    if not foundLowBinEdge:
        B_low = hBkg.Integral(1,ibin)
        #print "edge = ",hBkg.GetBinLowEdge(ibin+1)
        #print "B_low = ",B_low
        #if (B_low > 0):
        #    print "1./sqrt(B_low) = ",1./sqrt(B_low)
        if (B_low > 0 and 1./sqrt(B_low) < 0.35):
            binBoundaries[1] = hBkg.GetBinLowEdge(ibin+1)
            foundLowBinEdge = True
    if not foundHighBinEdge:
        B_high = hBkg.Integral(nBins-ibin+1, nBins) 
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
    cutString = presel + "&&("
    for index in sampleMap[sample]:
            cutString += "(sampleIndex==%i)||" % index
    cutString = cutString[0:len(cutString)-2]
    cutString += ")"
    #print cutString
    hBDT = ROOT.TH1F(sample,sample,nBins,-1,1)
    #tree.Draw("%s>>%s" % (bdtname, sample),"(%s)*weight" % cutString)
    #tree.Draw("%s>>%s" % (bdtname, sample),"(%s)*weight*(2.2/1.28)" % cutString) # temp hack to avoid rerunning just to change lumi
    tree.Draw("%s>>%s" % (bdtname, sample),"((%s)&&Pass_nominal)*weight*1000" % cutString) # temp hack to avoid rerunning just to change lumi
    hBDT = hBDT.Rebin(nBins, "", binBoundaries)
    # Add bin-by-bin stat. uncertainties
    if (sample not in ["WH","ZH","data_obs"]):
        for ibin in range(1, hBDT.GetNbinsX()):
            B = hBDT.GetBinContent(ibin)
            if ( B > 0 and (1./sqrt(B)) > tolerance ):
                #hBinStat = ROOT.TH1F("CMS_vhbb_stat%s_%s_bin%i_13TeV" % (sample,catName,ibin),"CMS_vhbb_stat%s_%s_bin%i_13TeV" % (sample,catName,ibin),nBins,-1,1)
                hBinStatUp = hBDT.Clone()
                hBinStatDown = hBDT.Clone()
                hBinStatUp.SetName("%s_CMS_vhbb_stat%s_%s_bin%i_13TeVUp" % (sample,sample,catName,ibin))
                hBinStatDown.SetName("%s_CMS_vhbb_stat%s_%s_bin%i_13TeVDown" % (sample,sample,catName,ibin))
                hBinStatUp.SetBinContent(ibin, B + sqrt(B))
                hBinStatDown.SetBinContent(ibin, B - sqrt(B))
                otextfile.write("CMS_vhbb_stat%s_%s_bin%i_13TeV\n" % (sample,catName,ibin))
                ofile.cd()
                hBinStatUp.Write()
                hBinStatDown.Write()
    # other shape systematics
    for syst in systematics:
        sysWeight, sysName, sysSamples = systematics[syst]
        if sample not in sysSamples: continue
        hBDTSystUp = ROOT.TH1F("%s_%sUp" % (sample,syst), "%s_%sUp" % (sample,syst),nBins,-1,1)
        hBDTSystDown = ROOT.TH1F("%s_%sDown" % (sample,syst), "%s_%sDown" % (sample,syst),nBins,-1,1)
        sysBDTNameUp = bdtname
        sysBDTNameDown = bdtname
        passSys = "Pass_nominal"
        if (sysName != ""):
            sysBDTNameUp += "_%sUp" % sysName
            sysBDTNameDown += "_%sDown" % sysName
            pasSys = "Pass_%s" % sysName 
        if (sysWeight != "1.0"):
            #tree.Draw("%s>>%s_%sUp" % (sysBDTNameUp, sample, syst),"(%s)*weight*(2.2/1.28)*(%sUp)" % (cutString,sysWeight)) # temp hack to avoid rerunning just to change lumi 
            tree.Draw("%s>>%s_%sUp" % (sysBDTNameUp, sample, syst),"((%s)&&%s)*weight*1000*(%sUp)" % (cutString,passSys,sysWeight)) # temp hack to avoid rerunning just to change lumi 
        else:
            #tree.Draw("%s>>%s_%sUp" % (sysBDTNameUp, sample, syst),"(%s)*weight*(2.2/1.28)" % (cutString)) # temp hack to avoid rerunning just to change lumi 
            tree.Draw("%s>>%s_%sUp" % (sysBDTNameUp, sample, syst),"((%s)&&%s)*weight*1000" % (cutString,passSys)) # temp hack to avoid rerunning just to change lumi 
        if (sysWeight != "1.0"):
            #tree.Draw("%s>>%s_%sDown" % (sysBDTNameDown, sample, syst),"(%s)*weight*(2.2/1.28)*(%sDown)" % (cutString,sysWeight)) # temp hack to avoid rerunning just to change lumi 
            tree.Draw("%s>>%s_%sDown" % (sysBDTNameDown, sample, syst),"((%s)&&%s)*weight*1000*(%sDown)" % (cutString,passSys,sysWeight)) # temp hack to avoid rerunning just to change lumi 
        else:
            #tree.Draw("%s>>%s_%sDown" % (sysBDTNameDown, sample, syst),"(%s)*weight*(2.2/1.28)" % (cutString)) # temp hack to avoid rerunning just to change lumi 
            tree.Draw("%s>>%s_%sDown" % (sysBDTNameDown, sample, syst),"((%s)&&%s)*weight*1000" % (cutString,passSys)) # temp hack to avoid rerunning just to change lumi 
        ofile.cd()
        hBDTSystUp.Write()
        hBDTSystDown.Write()
    ofile.cd()
    hBDT.Write()

ofile.Close()
otextfile.close()
   
