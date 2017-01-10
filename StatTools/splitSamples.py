import ROOT
import sys
from math import sqrt
from math import pow
import numpy
import argparse
import re

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
parser.add_argument('-d', '--doData', type=int, default=0, help="If non-zero run on real data, if 0 compute data as sum of background pdf's (default 0)")
parser.add_argument('-t', '--dataTree', type=str, default="", help="If doData is true, you specify here the ntuple with the real data events")
parser.add_argument('-tt', '--ttbarTree', type=str, default="", help="Specify separate file with ttbar powheg so it only takes forever when running tt")
parser.add_argument('-wj', '--wjetsTree', type=str, default="", help="Specify separate file with W+jets MC")
parser.add_argument('-v', '--varname', type=str, default="CMS_vhbb_BDT_Wln_13TeV", help="The name of the variable shape which goes into the histograms and will be fitted")
parser.add_argument('-xl','--xlow', type=float, default=-1.0, help="Lowest bin edge of fitted distribution")
parser.add_argument('-xh','--xhigh', type=float, default=1.0, help="Highest bin edge of fitted distribution")
parser.add_argument('-bh','--binhigh',type=float,default=1.0, help="manually set low bin edge of most sensitive bin")
parser.add_argument('-bl','--binlow',type=float,default=-1.0, help="manually set low bin edge of least sensitive bin")
parser.add_argument('-o', '--ofilename',type=str,default="", help="output histogram file name (default hists_[cat].root)")
parser.add_argument('-b', '--binstats', type=str, default="", help="Text file listing all the individual bin. stat. uncertainties to include")
parser.add_argument('-n', '--nbins', type=int, default=20, help="number of bins in datacard shapes")
parser.add_argument('-r', '--doRebin',type=int, default=1, help="If not 0 then rebin histogram so that bins are sufficiently populated in background (default=1)")
parser.add_argument('-tol','--tolerance',type=float,default=0.5,help="Tolerance threshold on dB/sqrt(B) for inclusion of bin-by-bin stat. shape uncertainties (default=0.50)")
parser.add_argument('-bb','--binBoundaries',type=str,default="",help="Manually specify the bin edges for the histogram that goes into the histogram")
parser.add_argument('-vv','--doVV', type=bool, default=False, help="If true do VV analysis (default False)")
args = parser.parse_args()
print args

#if (len(sys.argv) != 4 and len(sys.argv) != 5):
#    print "Takes three or four arguments:"
#    print "python splitSamples.py [input_root_file] [catName] [presel]"
#    print "python splitSamples.py [input_root_file] [catName] [presel] [sys.txt] (for shape systematics besides bin-by-bin)"
#    sys.exit(1)

ROOT.gROOT.SetBatch(True)

print "going to open this file now: ",args.inputfile
ifile = ROOT.TFile.Open(args.inputfile, "r")
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

set_of_weights = []
if (args.weights != ""):
    set_of_weights = args.weights.split(',')
print set_of_weights
#weight_string = "weight"
#weight_string = "weight*(bTagWeightEF/bTagWeight)*(1+isWenu*(-1+(SF_HLT_Ele23_WPLoose[lepInd]/EffHLT_Ele27_WPLoose_Eta2p1[lepInd])))"
#weight_string = "(1./CS_SF)*weight"
weight_string = "(1./CS_SF)*(bTagWeightICHEP/bTagWeight)*weight"
#weight_string = "(1./CS_SF)*weight*(1+isWenu*(-1+(SF_HLT_Ele23_WPLoose[lepInd]/EffHLT_Ele27_WPLoose_Eta2p1[lepInd])))"
#weight_string = "weight*(1+isWenu*(-1+(1./SF_egammaEffi_tracker[lepInd])))*(1+isWenu*(-1+(SF_HLT_Ele23_WPLoose[lepInd]/EffHLT_Ele27_WPLoose_Eta2p1[lepInd])))"
#weight_string = "weight*(1./SF_egammaEffi_tracker[lepInd])*(1+isWenu*(-1+(SF_HLT_Ele23_WPLoose[lepInd]/EffHLT_Ele27_WPLoose_Eta2p1[lepInd])))"
#weight_string = "weight*(puWeight/weight_PU)*(1./CS_SF)"
#weight_string = "weight*(1./CS_SF)*(puWeight/weight_PU)*(1 + (sampleIndex>=24&&sampleIndex<=31)*(-0.63)) * ( 1 + (sampleIndex>=24&&sampleIndex<=31&&isWmunu)*(-1 + 1./(0.0673 *SF_SMuTrig_Block1[lepInd] + 0.9327 * SF_SMuTrig_Block2[lepInd]))  ) * ( 1 + (sampleIndex>=24&&sampleIndex<=31&&isWenu)*(-1 + 1./(SF_HLT_Ele23_WPLoose[lepInd]))  )"
#weight_string = "weight*(1./CS_SF)*(puWeight/weight_PU)*(1 + (sampleIndex>=24&&sampleIndex<=31)*(-0.63)) * ( 1 + (sampleIndex>=24&&sampleIndex<=31&&isWmunu)*(-1 + 1./(0.0673 *SF_SMuTrig_Block1[lepInd] + 0.9327 * SF_SMuTrig_Block2[lepInd]))  ) * ( 1 + (sampleIndex>=24&&sampleIndex<=31&&isWenu)*(-1 + 1./(SF_HLT_Ele23_WPLoose[lepInd]))  )"
#weight_string = "weight*(CS_SF_new/CS_SF)"
#weight_string = "weight*((sampleIndex!=50&&sampleIndex!=51&&sampleIndex!=52) + (sampleIndex==50||sampleIndex==51||sampleIndex==52)*(245.79/831.76))" ## FIXME: change it back!!!!!
#weight_string = "weight*((sampleIndex!=50&&sampleIndex!=51&&sampleIndex!=52) + (sampleIndex==50||sampleIndex==51||sampleIndex==52)*(245.79/831.76)*1.956)" ## FIXME: change it back!!!!!
for weight in set_of_weights:
    weight_string += "*" + weight
print weight_string

tree_mc = ifile.Get("tree")
#bdtname = "BDT_wMass_Dec14_3000_5"
#bdtname = "BDT_wMass_Dec4"
#bdtname = "CMS_vhbb_BDT_Wln_13TeV"
bdtname = args.varname # not necessarily the bdt shape, can fit whatever shape is specified
#nBins = 1000
nBins = args.nbins  # number of bins in final histogram after rebinning
nBinsFine = 1000 # number of candidate bin edges to begin with
#tolerance = 1.5 # dB/B tolerance for bin-by-bin stat. uncertainties
#tolerance = 0.5 # dB/B tolerance for bin-by-bin stat. uncertainties
#tolerance = 0.75 # dB/B tolerance for bin-by-bin stat. uncertainties
tolerance = args.tolerance

sampleMap = {} # map sampleNames to list of sampleIndex's
sampleMapAltModel = {} # alternate MC samples for model shape systematics

#sampleMap["data_obs"] = [16,17,10,11,20,21,50,51,52,2200,2201,2202,2300,2301,2302,3500,3501,3502,4400,4401,4402,4500,4501,4502,4600,4601,4602,4700,4701,4702,24,25,26,27,28,29,30,31] # dummy filler until we unblind
#sampleMap["data_obs"] = [50,51,52]
sampleMap["data_obs"] = [0]
sampleMapAltModel["TT"] = [50,51,52]
sampleMap["TT"] = [120]
sampleMap["s_Top"] = [16,17,20,21]
#sampleMap["WH"] = [-12501]
sampleMap["WH"] = [-12500,-12501]
sampleMapAltModel["WH"] = [-125010, -125011]
sampleMap["ZH"] = [-12502]
sampleMap["Wj0b"] = [2200,4100,4200,4300,4400,4500,4600,4700,4800,4900]
sampleMapAltModel["Wj0b"] = [6000]
sampleMap["Wj1b"] = [2201,4101,4201,4301,4401,4501,4601,4701,4801,4901]
sampleMapAltModel["Wj1b"] = [6001]
sampleMap["Wj2b"] = [2202,4102,4202,4302,4402,4502,4602,4702,4802,4902]
sampleMapAltModel["Wj2b"] = [6002]
#sampleMap["VVHF"] = [35000,36000,37000,35001,36001,37001,35002,36002,37002]
#sampleMap["VVLF"] = [35000,36000,37000,35001,36001,37001,35002,36002,37002]
sampleMap["VVHF"] = [3500,3600,3700,3501,3601,3701,3502,3602,3702]
sampleMap["VVLF"] = [3500,3600,3700,3501,3601,3701,3502,3602,3702]
#sampleMap["VVHF"] = [3501,3502,3601,3602,3701,3702]
#sampleMap["VVLF"] = [3500,3600,3700]
sampleMap["QCD"]  = [24,25,26,27,28,29,30,31]
sampleMap["Zj0b"] = [2300,6100,6200,6300,6400]
sampleMap["Zj1b"] = [2301,6101,6201,6301,6401]
sampleMap["Zj2b"] = [2302,6102,6202,6302,6402]
#sampleMapAltModel["Zj0b"] = [23000]
#sampleMapAltModel["Zj1b"] = [23001]
#sampleMapAltModel["Zj2b"] = [23002]

sampleNameMap = {}
sampleNameMap["TT"] = ["TT_powheg"]
sampleNameMap["s_Top"] = ["TToLeptons_s","TToLeptons_t","T_tW","Tbar_tW"]
#sampleNameMap["WH"] = ["WH125_powheg"]
sampleNameMap["WH"] = ["WplusH125_powheg","WminusH125_powheg"]
sampleNameMap["ZH"] = ["ZH125_powheg"]
sampleNameMap["Wj0b"] = ["WJets_madgraph","WJets-HT800To1200","WJets-HT1200To2500","WJets-HT2500ToInf","WJets-HT100To200","WJets-HT200To400","WJets-HT400To600","WJets-HT600To800","WBJets","WJets_BGenFilter"]
sampleNameMap["Wj1b"] = ["WJets_madgraph","WJets-HT800To1200","WJets-HT1200To2500","WJets-HT2500ToInf","WJets-HT100To200","WJets-HT200To400","WJets-HT400To600","WJets-HT600To800","WBJets","WJets_BGenFilter"]
sampleNameMap["Wj2b"] = ["WJets_madgraph","WJets-HT800To1200","WJets-HT1200To2500","WJets-HT2500ToInf","WJets-HT100To200","WJets-HT200To400","WJets-HT400To600","WJets-HT600To800","WBJets","WJets_BGenFilter"]
#sampleNameMap["VVHF"] = ["WZ","WW","ZZ"]
#sampleNameMap["VVLF"] = ["WZ","WW","ZZ"]
sampleNameMap["VVHF"] = ["WZ_fil","WW_fil","ZZ_fil"]
sampleNameMap["VVLF"] = ["WZ_fil","WW_fil","ZZ_fil"]
sampleNameMap["QCD"]  = []
sampleNameMap["Zj0b"] = ["DYToLL_madgraph","DYToLL_HT100to200","DYToLL_HT200to400","DYToLL_HT400to600","DYToLL_HT600toInf"]
sampleNameMap["Zj1b"] = ["DYToLL_madgraph","DYToLL_HT100to200","DYToLL_HT200to400","DYToLL_HT400to600","DYToLL_HT600toInf"]
sampleNameMap["Zj2b"] = ["DYToLL_madgraph","DYToLL_HT100to200","DYToLL_HT200to400","DYToLL_HT400to600","DYToLL_HT600toInf"]

allSampInd = [] # list of all indices for all backgrounds
sigSamps = ["WH","ZH"]
if args.doVV:
    sigSamps = ["VVHF","VVLF"]
print "sigSamps = "
print sigSamps
for sample in sampleMap:
    if (sample in sigSamps or sample == "data_obs"): continue
    allSampInd.extend(sampleMap[sample])
sampleMap["Bkg"] = allSampInd
print "all Bkg sampleInd = "
print allSampInd
if not args.doData:
    sampleMap["data_obs"] = allSampInd # fake data as sum of all background MC

def makeCutString(sample, sMap):
    cutString = presel + "&&("
    if (sample == "VVHF"):
        cutString = presel + "&&Sum$(abs(GenWZQuark_pdgId)==5)>=2&&("
    if (sample == "VVLF"):
        cutString = presel + "&&Sum$(abs(GenWZQuark_pdgId)==5)<2&&("
    for index in sMap[sample]:
            #if (index == 4900 or index == 4901 or index == 4902): continue # HACK HACK FIXME
            cutString += "(sampleIndex==%i)||" % index
    cutString = cutString[0:len(cutString)-2]
    cutString += ")"
    return cutString

print "Bkg string = ",makeCutString("Bkg", sampleMap)
#print "Wj0b string = ",makeCutString("Wj0b", sampleMap)
print "Thank you for your attention!"

# first rebin the histogram so that the first and last bins are not empty and have less than 35% stat. uncertainty
hBkg = ROOT.TH1F("hBkg","hBkg",nBinsFine,args.xlow,args.xhigh)
hSig = ROOT.TH1F("hSig","hSig",nBinsFine,args.xlow,args.xhigh)
#tree.Draw("%s>>hBkg" % bdtname,"((%s)&&sampleIndex>0)*weight*(2.2/1.28)" % presel)
bkgCutString = makeCutString("Bkg", sampleMap)
print bkgCutString
sigCutString = makeCutString("WH", sampleMap)
if args.doVV:
    sigCutString = makeCutString("VVHF", sampleMap)
tree_mc.Draw("%s>>hBkg" % bdtname,"((%s)&&(%s))*%s" % (presel,bkgCutString,weight_string))
if args.doVV:
    tree_mc.Draw("%s>>hSig" % bdtname,"((%s)&&(%s))*%s" % (presel,sigCutString,weight_string))

ifile_sig = ROOT.TFile.Open(args.inputfile.replace("output_mc","output_signal"))
tree_sig = ifile_sig.Get("tree")
ifile.cd()
if args.doVV:
    hBkgWH = ROOT.TH1F("hBkgWH","hBkgWH",nBinsFine,args.xlow,args.xhigh)
    tree_sig.Draw("%s>>hBkgWH" % bdtname,"((%s)&&(%s))*%s" % (presel,bkgCutString,weight_string))
    print hBkg.Integral()
    print "adding VH contribution to background.."
    hBkg.Add(hBkgWH) 
    print hBkg.Integral()
else: 
    tree_sig.Draw("%s>>hSig" % bdtname,"((%s)&&(%s))*%s" % (presel,sigCutString,weight_string))
ifile_sig.Close()

if (args.ttbarTree != ""):
    print "got here"
    ifile_tt = ROOT.TFile.Open(args.ttbarTree,"r")
    hBkg2 = ROOT.TH1F("hBkg2","hBkg2",nBinsFine,args.xlow,args.xhigh)
    tree_tt = ifile_tt.Get("tree")
    # make sure we don't weight actual data by puWeight, SF's, etc.
    #print tree_tt.GetEntries("sampleIndex==120")
    tree_tt.Draw("%s>>hBkg2" % bdtname,"((%s)&&(%s))*%s" % (presel,bkgCutString,weight_string))
    print hBkg.Integral()
    print hBkg2.Integral()
    hBkg.Add(hBkg2)
    print hBkg.Integral()
    ifile_tt.Close() 

print "total background, total sig:",hSig.Integral(),", ",hBkg.Integral()

binBoundaries = numpy.zeros(nBins+1,dtype=float)
binBoundaries[0] = args.xlow
binBoundaries[nBins] = args.xhigh
foundLowBinEdge = False
foundHighBinEdge = False
if (args.binlow != -1):
    foundLowBinEdge = True
    binBoundaries[1] = args.binlow
if (args.binhigh != 1):
    foundHighBinEdge = True
    binBoundaries[nBins-1] = args.binhigh
B_err2_low = 0.
B_err2_high = 0.
for ibin in range(1,nBinsFine):
    if not foundLowBinEdge:
        B_low = hBkg.Integral(1,ibin)
        B_err2_low += pow(hBkg.GetBinError(ibin),2) 
        #S_low = hSig.Integral(1,ibin)
        #print "edge = ",hBkg.GetBinLowEdge(ibin+1)
        #print "B_low = ",B_low
        #if (B_low > 0):
        #    print "1./sqrt(B_low) = ",1./sqrt(B_low)
        #print B_low,ibin,hBkg.Integral(),args.xlow,args.xhigh
        #print "B_low,sqrt(B_err2_low)",B_low,sqrt(B_err2_low)
        #if (B_low > 0 and 1./sqrt(B_low) < 0.35):
        if (B_low > 0 and sqrt(B_err2_low)/B_low < 0.35):
            binBoundaries[1] = hBkg.GetBinLowEdge(ibin+1)
            foundLowBinEdge = True
            #print "found low bin!",ibin
    if not foundHighBinEdge:
        B_high = hBkg.Integral(nBinsFine-ibin+1, nBinsFine) 
        B_err2_high += pow(hBkg.GetBinError(nBinsFine-ibin+1),2)
        S_high = hSig.Integral(nBinsFine-ibin+1, nBinsFine) 
        #S_high = hSig.Integral(nBinsFine-ibin+1, nBinsFine) 
        #print "B_high,sqrt(B_err2_high)",B_high,sqrt(B_err2_high)
        #if (B_high > 0 and 1./sqrt(B_high) < 0.35):
        #if (B_high > 0 and sqrt(B_err2_high)/B_high < 0.35 and S_high >= 1.0 and abs(hBkg.GetBinLowEdge(nBinsFine-ibin+1)-0.608)<0.001):
        #if (B_high > 0 and sqrt(B_err2_high)/B_high < 0.35 and S_high >= 1.0):
        #if (B_high > 0 and sqrt(B_err2_high)/B_high < 0.35 and S_high >= 6.0):
        if (B_high > 0 and sqrt(B_err2_high)/B_high < 0.35):
            binBoundaries[nBins-1] = hBkg.GetBinLowEdge(nBinsFine-ibin+1)
            foundHighBinEdge = True
            #print "found high bin!",(nBinsFine-ibin+1)
# split the middle bins equidistantly
## HACK FIXME
#binBoundaries[nBins-2] = 0.376
#for i in range(2,nBins-2):
#    binBoundaries[i] = binBoundaries[1] + (i-1)*((binBoundaries[nBins-2] - binBoundaries[1])/(nBins-3)) 
for i in range(2,nBins-1):
    binBoundaries[i] = binBoundaries[1] + (i-1)*((binBoundaries[nBins-1] - binBoundaries[1])/(nBins-2)) 

#if True:
if not args.doRebin:
    for i in range(1,nBins):
        binBoundaries[i] = binBoundaries[0] + (i)*((binBoundaries[nBins] - binBoundaries[0])/(nBins))
if (args.binBoundaries != ""):
    boundaryString = args.binBoundaries
    print "parsing boundaries from ",boundaryString
    boundaries = boundaryString.split(',')
    nBins = len(boundaries)-1
    binBoundaries = numpy.zeros(nBins+1,dtype=float)
    for i in range(len(boundaries)):
        binBoundaries[i] = float(boundaries[i].replace('m','-'))
    print binBoundaries
    
hBkg = hBkg.Rebin(nBins, "", binBoundaries)
print binBoundaries
if args.ofilename == "":
    ofile = ROOT.TFile.Open("hists_%s.root" % catName, "RECREATE")
else:
    ofile = ROOT.TFile.Open(args.ofilename, "RECREATE")
if args.binstats == "":
    otextfile = open("binStats_%s.txt" % catName, "w")
else:
    otextfile = open(args.binstats,"w")
tree = ROOT.TTree("tree","tree")
#hBkg.Write()
for sample in sampleMap:
    #if (sample!="s_Top" and sample!="VVHF" and sample!="VVLF"): continue
    #if (sample != "Wj0b" and sample!="Wj1b" and sample!="Wj2b"): continue
    #if (sample != "Bkg" and sample != "data_obs" and sample != "WH" and sample != "ZH"): continue
    #if (sample != "WH"): continue
    if (sample == "TT" and args.ttbarTree != ""): 
        ifile_tt = ROOT.TFile.Open(args.ttbarTree,"r")
        tree_tt = ifile_tt.Get("tree")
        tree = tree_tt
        #ofile.cd()
        #ifile_tt.Close()
    elif (sample == "WH" or sample == "ZH"):
        ifile_sig = ROOT.TFile.Open(args.inputfile.replace("output_mc","output_signal"))
        tree_sig = ifile_sig.Get("tree")
        tree = tree_sig
    #elif (sample == "VVHF" or sample == "VVLF"):
    #    ifile_vv = ROOT.TFile.Open(args.inputfile.replace("output_mc","output_vv"))
    #    tree_vv = ifile_vv.Get("tree")
    #    tree = tree_vv
    elif (args.wjetsTree != "" and (sample == "Wj0b" or sample == "Wj1b" or sample=="Wj2b")):
        ifile_wjets = ROOT.TFile.Open(args.wjetsTree)
        tree_wjets = ifile_wjets.Get("tree")
        tree = tree_wjets
    else:
        tree = tree_mc  
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
    #tree.Draw("%s>>%s" % (bdtname, sample),"(%s)*weight*(2.2/1.28)" % cutString)  
    if (sample == "data_obs" and args.doData):
        ifile_data = ROOT.TFile.Open(args.dataTree,"r")
        tree_data = ifile_data.Get("tree")
        ofile.cd()
        # make sure we don't weight actual data by puWeight, SF's, etc.
        tree_data.Draw("%s>>%s" % (bdtname, sample),"((%s)&&Pass_nominal)" % (cutString)) 
        ifile_data.Close()
    elif (sample == "data_obs"):
        # fake data which is sum of all MC
        #hBkg.Write(sample)
        ofile.cd()
        hBkg.Write("BDT_%s_%s" % (catName,sample))
         
    #    continue
    #elif (sample == "TT" and args.ttbarTree != ""):
    #    print "got here"
    #    ifile_tt = ROOT.TFile.Open(args.ttbarTree,"r")
    #    tree_tt = ifile_tt.Get("tree")
    #    ofile.cd()
    #    # make sure we don't weight actual data by puWeight, SF's, etc.
    #    tree_tt.Draw("%s>>%s" % (bdtname, sample),"((%s)&&Pass_nominal)*%s" % (cutString,weight_string)) 
    #    ifile_tt.Close() 
    else:
        tree.Draw("%s>>%s" % (bdtname, sample),"((%s)&&Pass_nominal)*%s" % (cutString,weight_string)) 
    print "tree.Draw(\"%s>>%s\",\"((%s)&&Pass_nominal)*%s\"" % (bdtname,sample,cutString,weight_string) 
    #hBDT = hBDT.Rebin(nBins, "", binBoundaries)
    # Add bin-by-bin stat. uncertainties
    if (sample not in ["data_obs"] and (tolerance <= 0.5 or sample not in ["QCD","VVLF","VVHF","WH","ZH"])): # assuming for SR the tolerance does not go above 0.50
    #if (sample not in ["QCD","VVLF","VVHF","WH","ZH","data_obs"]): # can exclude these when running on background to reduce the number of nuisances
        for ibin in range(1, hBDT.GetNbinsX()):
            B = hBDT.GetBinContent(ibin)
            B_err = hBDT.GetBinError(ibin)
            #print "ibin = %i" % ibin
            #print "B = %f" % B
            #print "B_err = %f" % B_err
            #print "GetEntries = %i" % hBDT.GetEntries()
            #print "GetSumOfWeights = %f" % hBDT.GetSumOfWeights()
            #if (B > 0):
                #print "B_err/sqrt(B) = %f" % (B_err/sqrt(B))
                #print "B_err/B = %f" % (B_err/B)
            #B_eff = hBDT.GetEffectiveEntries() # effective statistical number of entries in bin considering weights
            if ( B > 0 and ( ( B >=1 and (B_err/sqrt(B)) > tolerance) or (B < 1 and B_err/B > tolerance) ) ):
                print "qualified as a bin stat. shape uncertainty! Bin %i, sample %s" % (ibin,sample)
                if (B >=1): print B_err/sqrt(B)
                else: print B_err/B
                print "Bin content: %f +/- %f" % (B,B_err)
                #hBinStat = ROOT.TH1F("CMS_vhbb_stat%s_%s_bin%i_13TeV" % (sample,catName,ibin),"CMS_vhbb_stat%s_%s_bin%i_13TeV" % (sample,catName,ibin),nBins,-1,1)
                hBinStatUp = hBDT.Clone()
                hBinStatDown = hBDT.Clone()
                hBinStatUp.SetName("BDT_%s_%s_CMS_vhbb_stat%s_%s_bin%i_13TeVUp" % (catName,sample,sample,catName,ibin))
                hBinStatDown.SetName("BDT_%s_%s_CMS_vhbb_stat%s_%s_bin%i_13TeVDown" % (catName,sample,sample,catName,ibin))
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
        #if (syst.find("LHE")==-1): continue
        ofile.cd()
        hBDTSystUp = ROOT.TH1F("BDT_%s_%s_%sUp" % (catName,sample,syst), "%s_%sUp" % (sample,syst),nBins,binBoundaries)
        hBDTSystDown = ROOT.TH1F("BDT_%s_%s_%sDown" % (catName,sample,syst), "%s_%sDown" % (sample,syst),nBins,binBoundaries)
        sysBDTNameUp = bdtname
        sysBDTNameDown = bdtname
        passSys = "Pass_nominal"
        if (sysName != ""):
            if (bdtname.find('[') != -1):
                # drawing an array variable
                sysBDTNameUp = bdtname[:bdtname.find('[')] + "_" + sysName + "Up" + bdtname[bdtname.find('['):] 
            else:
                sysBDTNameUp += "_%sUp" % sysName
                sysBDTNameDown += "_%sDown" % sysName
            pasSys = "Pass_%s" % sysName
        if (sysWeight != "1.0" and sysWeight.find(".root") == -1):
            if (sysWeight.find(',') == -1):
                if (sysWeight.find("bTagWeight") != -1):
                    bTagWeightNom = sysWeight[:sysWeight.find('_')]
                    print "bTagWeightNom = "+bTagWeightNom
                    tree.Draw("%s>>BDT_%s_%s_%sUp" % (sysBDTNameUp, catName,sample, syst),"((%s)&&%s)*(1./%s)*%s*(%sUp)" % (cutString,passSys,bTagWeightNom,weight_string,sysWeight))   
                    tree.Draw("%s>>BDT_%s_%s_%sDown" % (sysBDTNameDown, catName, sample, syst),"((%s)&&%s)*(1./%s)*%s*(%sDown)" % (cutString,passSys,bTagWeightNom,weight_string,sysWeight))   
                #elif (sysWeight.find("weight_PU") != -1):
                #    tree.Draw("%s>>BDT_%s_%s_%sUp" % (sysBDTNameUp, catName,sample, syst),"((%s)&&%s)*(1./weight_PU)*%s*(%sUp)" % (cutString,passSys,weight_string,sysWeight))   
                #    tree.Draw("%s>>BDT_%s_%s_%sDown" % (sysBDTNameDown, catName, sample, syst),"((%s)&&%s)*(1./weight_PU)*%s*(%sDown)" % (cutString,passSys,weight_string,sysWeight))   
                #elif (sysWeight.find("weight_PU") != -1): ## FIXMEEE
                #     if (sample != "Wj0b" and sample != "Wj1b" and sample != "Wj2b"):
                #         tree.Draw("%s>>BDT_%s_%s_%sUp" % (sysBDTNameUp, catName,sample, syst),"((%s)&&%s)*(1./weight_PU)*%s*(%sUp)" % (cutString,passSys,weight_string,sysWeight))   
                #         tree.Draw("%s>>BDT_%s_%s_%sDown" % (sysBDTNameDown, catName, sample, syst),"((%s)&&%s)*(1./weight_PU)*%s*(%sDown)" % (cutString,passSys,weight_string,sysWeight))
                #     else:
                #         tree.Draw("%s>>BDT_%s_%s_%sUp" % (sysBDTNameUp, catName,sample, syst),"((%s)&&%s)*%s*(%sUp)" % (cutString,passSys,weight_string,sysWeight))   
                #         tree.Draw("%s>>BDT_%s_%s_%sDown" % (sysBDTNameDown, catName, sample, syst),"((%s)&&%s)*%s*(%sDown)" % (cutString,passSys,weight_string,sysWeight))
                #elif (sysWeight.find("puWeight") != -1):
                #    tree.Draw("%s>>BDT_%s_%s_%sUp" % (sysBDTNameUp, catName,sample, syst),"((%s)&&%s)*(1./puWeight)*%s*(%sUp)" % (cutString,passSys,weight_string,sysWeight))   
                #    tree.Draw("%s>>BDT_%s_%s_%sDown" % (sysBDTNameDown, catName, sample, syst),"((%s)&&%s)*(1./puWeight)*%s*(%sDown)" % (cutString,passSys,weight_string,sysWeight))   
                else:
                    tree.Draw("%s>>BDT_%s_%s_%sUp" % (sysBDTNameUp, catName,sample, syst),"((%s)&&%s)*%s*(%sUp)" % (cutString,passSys,weight_string,sysWeight))   
                    tree.Draw("%s>>BDT_%s_%s_%sDown" % (sysBDTNameDown, catName, sample, syst),"((%s)&&%s)*%s*(%sDown)" % (cutString,passSys,weight_string,sysWeight))  
            else:
                # separate branches for up/down variation
                sysWeightUp = sysWeight.split(',')[0]
                sysWeightDown = sysWeight.split(',')[1]
                # special case for LHE scale and pdf variations, we should probably find a way to make this cleaner
                if (syst.find("LHE_weights_scale")!=-1 or syst.find("LHE_weights_pdf")!=-1):
                    #if (sample != "TT"):
                    #    sysWeightUp += "*LHE_weights_scale_normwgt[%s]" % re.findall('\d+', sysWeightUp)[0]
                    #    sysWeightDown += "*LHE_weights_scale_normwgt[%s]" % re.findall('\d+', sysWeightDown)[0]
                    #else:
                    #    # have to do something special here since we didn't hadd everything on the ttbar jobs so that it doesn't take 100 million years to run the jobs
                    #    ifile_ttpowhegcounts = ROOT.TFile.Open(args.inputfile.replace("output_mc","ttpowheg_counts"))
                    #    CountWeightedLHEWeightScale_TT_powheg = ifile_ttpowhegcounts.Get("CountWeightedLHEWeightScale_TT_powheg")
                    #    normweightUp = CountWeightedLHEWeightScale_TT_powheg.GetBinContent(CountWeightedLHEWeightScale_TT_powheg.FindBin(int(re.findall('\d+', sysWeightUp)[0])))
                    #    normweightDown = CountWeightedLHEWeightScale_TT_powheg.GetBinContent(CountWeightedLHEWeightScale_TT_powheg.FindBin(int(re.findall('\d+', sysWeightDown)[0])))
                    #    sysWeightUp += "*(nProcEvents/%f)" % normweightUp
                    #    sysWeightDown += "*(nProcEvents/%f)" % normweightDown
                    #    ifile_ttpowhegcounts.Close()
                    #    ofile.cd()
                    print "doing LHE_weights_scale for sample: ",sample
                    #if (sample == "TT"):
                    #    # have to do something special here since we didn't hadd everything on the ttbar jobs so that it doesn't take 100 million years to run the jobs
                    #    ifile_counts = ROOT.TFile.Open(args.inputfile.replace("output_mc","ttpowheg_counts"))
                    #elif (sample == "WH" or sample == "ZH"):
                    #    ifile_counts = ROOT.TFile.Open(args.inputfile.replace("output_mc","output_signal"))
                    #else:
                    #    ifile_counts = ROOT.TFile.Open(args.inputfile.replace("output_mc","%s/output_%s_1.root" %(sample,sample)))
                    #    #ifile_counts = ROOT.TFile.Open(args.inputfile)
                    #print "opened file fine"
                    #CountWeightedLHEWeightScale = ROOT.TH1F("CountWeightedLHEWeightScale_%s" % sample,"CountWeightedLHEWeightScale_%s" % sample,6,-0.5,5.5)
                    #CountWeightedLHEWeightPdf = ROOT.TH1F("CountWeightedLHEWeightPdf_%s" % sample,"CountWeightedLHEWeightPdf_%s" % sample,103,-0.5,102.5)
                    print "initiated histogram"
                    print "looping through..."
                    print sampleNameMap[sample]
                    for i in range(len(sampleNameMap[sample])):
                        sampleName = sampleNameMap[sample][i]
                        print sampleName
                        #if (sample == "TT"):
                        #    # have to do something special here since we didn't hadd everything on the ttbar jobs so that it doesn't take 100 million years to run the jobs
                        #    ifile_counts = ROOT.TFile.Open(args.inputfile.replace("output_mc","ttpowheg_counts"))
                        #elif (sample == "WH" or sample == "ZH"):
                        #    ifile_counts = ROOT.TFile.Open(args.inputfile.replace("output_mc","output_signal"))
                        #else:
                        #    print args.inputfile.replace("output_mc","%s/output_%s_1.root" %(sampleName,sampleName))
                        #    ifile_counts = ROOT.TFile.Open(args.inputfile.replace("output_mc","%s/output_%s_1.root" %(sampleName,sampleName)))
                        #    #ifile_counts = ROOT.TFile.Open(args.inputfile)
                        ifilename_counts = ""
                        if (sample == "TT"):
                            #print args.inputfile.replace("output_mc.root","TT_Powheg/%s/output_%s_1.root" %(sampleName,sampleName))
                            #ifilename_counts = args.inputfile.replace("output_mc.root","TT_Powheg/%s/output_%s_1.root" %(sampleName,sampleName))
                            print args.inputfile.replace("output_mc.root","ttpowheg_counts.root")
                            ifilename_counts = args.inputfile.replace("output_mc.root","ttpowheg_counts.root")
                        elif (sample == "WH" or sample == "ZH"):
                             print args.inputfile.replace("output_mc.root","Signal/%s/output_%s_1.root" %(sampleName,sampleName))
                             ifilename_counts = args.inputfile.replace("output_mc.root","Signal/%s/output_%s_1.root" %(sampleName,sampleName))
                        elif (args.wjetsTree != "" and (sample=="Wj0b" or sample=="Wj1b" or sample=="Wj2b")):
                            print args.wjetsTree.replace("output_wjets.root","%s/output_%s_1.root" %(sampleName,sampleName))
                            ifilename_counts = args.wjetsTree.replace("output_wjets.root","%s/output_%s_1.root" %(sampleName,sampleName))
                        else:
                            print args.inputfile.replace("output_mc.root","%s/output_%s_1.root" %(sampleName,sampleName))
                            ifilename_counts = args.inputfile.replace("output_mc.root","%s/output_%s_1.root" %(sampleName,sampleName))
                        ifile_counts = ROOT.TFile.Open(ifilename_counts)
                        #ifile_counts = ROOT.TFile.Open(args.inputfile.replace("output_mc","%s/output_%s_1.root" %(sampleName,sampleName)))
                        #ifile_counts = ROOT.TFile.Open(args.inputfile)
                        if (syst.find("LHE_weights_scale")!=-1):
                            CountWeightedLHEWeightScale = ifile_counts.Get("CountWeightedLHEWeightScale_%s" % sampleName)
                        else:
                            CountWeightedLHEWeightScale = ifile_counts.Get("CountWeightedLHEWeightPdf_%s" % sampleName)
                        print "got tmp histogram"
                        #CountWeightedLHEWeightPdf_tmp = ifile_counts.Get("CountWeightedLHEWeightPdf_%s" % sampleName)
                        #CountWeightedLHEWeightScale.Add(CountWeightedLHEWeightScale_tmp)
                        #print "add tmp histogram"
                        #CountWeightedLHEWeightPdf.Add(CountWeightedLHEWeightPdf_tmp)
                        print CountWeightedLHEWeightScale.GetName()
                        print CountWeightedLHEWeightScale.GetBinContent(1)
                        print CountWeightedLHEWeightScale.FindBin(int(re.findall('\d+', sysWeightUp)[0]))
                        normweightUp = CountWeightedLHEWeightScale.GetBinContent(CountWeightedLHEWeightScale.FindBin(int(re.findall('\d+', sysWeightUp)[0])))
                        normweightDown = CountWeightedLHEWeightScale.GetBinContent(CountWeightedLHEWeightScale.FindBin(int(re.findall('\d+', sysWeightDown)[0])))
                        print "calculated norm weights: ",normweightUp,", ",normweightDown
                        sysWeightUp += "*(1 + (sampleIndex==%i)*(-1+(nProcEvents/%f)))" % (int(sampleMap[sample][i]),normweightUp)
                        sysWeightDown += "*(1 + (sampleIndex==%i)*(-1+(nProcEvents/%f)))" % (int(sampleMap[sample][i]),normweightDown)
                        CountWeightedLHEWeightScale.Reset()
                        ifile_counts.Close()
                    #normweightUp = CountWeightedLHEWeightScale.GetBinContent(CountWeightedLHEWeightScale.FindBin(int(re.findall('\d+', sysWeightUp)[0])))
                    #normweightDown = CountWeightedLHEWeightScale.GetBinContent(CountWeightedLHEWeightScale.FindBin(int(re.findall('\d+', sysWeightDown)[0])))
                    #normweightAv = (normweightUp + normweightDown) / 2.
                    #sysWeightUp += "*(%f/%f)" % (normweightAv,normweightUp)
                    #sysWeightDown += "*(%f/%f)" % (normweightAv,normweightUp)
                    #sysWeightUp += "*(nProcEvents/%f)" % normweightUp
                    #sysWeightDown += "*(nProcEvents/%f)" % normweightDown
                    #ifile_counts.Close()
                    ofile.cd()
                    print "all done"

                #if (syst.find("LHE_weights_pdf")!=-1):
                #    if (sample != "TT"):
                #        sysWeightUp += "*LHE_weights_pdf_normwgt[%s]" % re.findall('\d+', sysWeightUp)[0] 
                #        sysWeightDown += "*LHE_weights_pdf_normwgt[%s]" % re.findall('\d+', sysWeightDown)[0]
                #    else:
                #        # have to do something special here since we didn't hadd everything on the ttbar jobs so that it doesn't take 100 million years to run the jobs
                #        ifile_ttpowhegcounts = ROOT.TFile.Open(args.inputfile.replace("output_mc","ttpowheg_counts"))
                #        CountWeightedLHEWeightPdf_TT_powheg = ifile_ttpowhegcounts.Get("CountWeightedLHEWeightPdf_TT_powheg")
                #        normweightUp = CountWeightedLHEWeightPdf_TT_powheg.GetBinContent(CountWeightedLHEWeightPdf_TT_powheg.FindBin(int(re.findall('\d+', sysWeightUp)[0])))
                #        normweightDown = CountWeightedLHEWeightPdf_TT_powheg.GetBinContent(CountWeightedLHEWeightPdf_TT_powheg.FindBin(int(re.findall('\d+', sysWeightDown)[0])))
                #        sysWeightUp += "*(nProcEvents/%f)" % normweightUp
                #        sysWeightDown += "*(nProcEvents/%f)" % normweightDown
                #        ifile_ttpowhegcounts.Close()
                #        ofile.cd()
                #print sysWeightUp,sysWeightDown 
                      
                #print "tree.Draw(\"%s>>BDT_%s_%s_%sUp\",\"((%s)&&%s)*%s*(%s)\")" % (sysBDTNameUp, catName,sample, syst,cutString,passSys,weight_string,sysWeightUp)
                tree.Draw("%s>>BDT_%s_%s_%sUp" % (sysBDTNameUp, catName,sample, syst),"((%s)&&%s)*%s*(%s)" % (cutString,passSys,weight_string,sysWeightUp)) 
                #print hBDTSystUp.Integral() 
                #print "tree.Draw(\"%s>>BDT_%s_%s_%sUp\",\"((%s)&&%s)*%s*(%s)\")" % (sysBDTNameDown, catName,sample, syst,cutString,passSys,weight_string,sysWeightDown)
                tree.Draw("%s>>BDT_%s_%s_%sDown" % (sysBDTNameDown, catName, sample, syst),"((%s)&&%s)*%s*(%s)" % (cutString,passSys,weight_string,sysWeightDown))  
                #print hBDTSystDown.Integral()
                
                if (syst.find("LHE_weights_pdf")!=-1): 
                    # new prescription for LHE PDF nuisances, treat each variation like an alternative generator shape
                    #if (hBDTSystDown.Integral() != 0):
                    #    hBDTSystDown.Scale(hBDT.Integral()/hBDTSystDown.Integral())
                    # hnom + (hnom - halt) = 2*hnom - halt
                    hBDTSystUp = hBDT.Clone("BDT_%s_%s_%sUp" %(catName, sample, syst))
                    hBDTSystUp.Scale(2.0)
                    hBDTSystUp.Add(hBDTSystDown, -1.0)
                    # hnom - (hnom - halt) = halt
                    for ibin in range(hBDTSystUp.GetNbinsX()):
                        B_Up = hBDTSystUp.GetBinContent(ibin)
                        B_Down = hBDTSystDown.GetBinContent(ibin)
                        if (B_Up < 0):
                            B_Up_Err = hBDTSystUp.GetBinError(ibin)
                            hBDTSystUp.SetBinContent(ibin, 0.0)
                            hBDTSystUp.SetBinError(ibin, B_Up + B_Up_Err)
                        if (B_Down < 0):
                            B_Down_Err = hBDTSystDown.GetBinError(ibin)
                            hBDTSystDown.SetBinContent(ibin, 0.0)
                            hBDTSystDown.SetBinError(ibin, B_Down + B_Down_Err)
        elif (sysWeight != "1.0"):
            print "made it here"
            # read in up/down weights from a root file
            weightfilename,weighthistname = sysWeight.split(':')
            weightfile = ROOT.TFile.Open(weightfilename)
            weighthist = weightfile.Get(weighthistname)
            ofile.cd()
            hBDTSystUp = hBDT + (hBDT * weighthist)
            hBDTSystDown = hBDT - (hBDT * weighthist) 
            # name gets set to the name of hBDT when you do as above, need to change it back
            hBDTSystUp.SetName("BDT_%s_%s_%sUp" % (catName,sample,syst))
            hBDTSystDown.SetName("BDT_%s_%s_%sDown" % (catName,sample,syst))
            print hBDTSystUp.GetName(),hBDTSystDown.GetName()
            weightfile.Close()
        else:
            #tree.Draw("%s>>%s_%sUp" % (sysBDTNameUp, sample, syst),"(%s)*weight*(2.2/1.28)" % (cutString))   
            #tree.Draw("%s>>%s_%sUp" % (sysBDTNameUp, sample, syst),"((%s)&&%s)*%s" % (cutString,passSys,weight_string))   
            tree.Draw("%s>>BDT_%s_%s_%sUp" % (sysBDTNameUp, catName,sample, syst),"((%s)&&%s)*%s" % (cutString,passSys,weight_string))   
            #tree.Draw("%s>>%s_%sDown" % (sysBDTNameDown, sample, syst),"(%s)*weight*(2.2/1.28)" % (cutString))   
            #tree.Draw("%s>>%s_%sDown" % (sysBDTNameDown, sample, syst),"((%s)&&%s)*%s" % (cutString,passSys,weight_string))   
            tree.Draw("%s>>BDT_%s_%s_%sDown" % (sysBDTNameDown, catName,sample, syst),"((%s)&&%s)*%s" % (cutString,passSys,weight_string))   
            if (syst.find("Model") != -1):
                # shape from different sample (amc@NLO vs. POWHEG, for example)
                cutStringAltModel = makeCutString(sample,sampleMapAltModel)
                hBDTAltModel = ROOT.TH1F("%s_%sAltModel" % (sample,syst), "%s_%sAltModel" % (sample,syst),nBins,binBoundaries)
                tree.Draw("%s>>%s_%sAltModel" % (sysBDTNameUp, sample, syst),"((%s)&&%s)*%s" % (cutStringAltModel,passSys,weight_string))
                if (hBDTAltModel.Integral() != 0):
                    hBDTAltModel.Scale(hBDT.Integral()/hBDTAltModel.Integral())
                # hnom + (hnom - halt) = 2*hnom - halt
                hBDTSystUp.Scale(2.0)
                hBDTSystUp.Add(hBDTAltModel, -1.0)
                # hnom - (hnom - halt) = halt
                hBDTSystDown = hBDTAltModel
                hBDTSystDown.SetName("BDT_%s_%s_%sDown" % (catName,sample, syst))
                for ibin in range(hBDTSystUp.GetNbinsX()):
                    B_Up = hBDTSystUp.GetBinContent(ibin)
                    B_Down = hBDTSystDown.GetBinContent(ibin)
                    if (B_Up < 0):
                        B_Up_Err = hBDTSystUp.GetBinError(ibin)
                        hBDTSystUp.SetBinContent(ibin, 0.0)
                        hBDTSystUp.SetBinError(ibin, B_Up + B_Up_Err)
                    if (B_Down < 0):
                        B_Down_Err = hBDTSystDown.GetBinError(ibin)
                        hBDTSystDown.SetBinContent(ibin, 0.0)
                        hBDTSystDown.SetBinError(ibin, B_Down + B_Down_Err)
        # Combine cannot handle shapes with negative norms
        if (hBDT.Integral() < 0):
            hBDT.Scale(-0.00001) ## if it was just 0.0 we would lose the shape
        if (hBDTSystUp.Integral() < 0):
            hBDTSystUp.Scale(-0.00001) ## if it was just 0.0 we would lose the shape
        if (hBDTSystDown.Integral() < 0):
            hBDTSystDown.Scale(-0.00001) ## if it was just 0.0 we would lose the shape

        ofile.cd()
        hBDTSystUp.Write()
        hBDTSystDown.Write()
        print syst
        print sysWeight
        print hBDT.Integral()
        print hBDTSystUp.Integral(),hBDTSystDown.Integral()
        hBDTSystUp.Delete()
        hBDTSystDown.Delete()
    ofile.cd()
    if (sample != "data_obs"):
        print True
    else:
        print False
    print args.doData
    print sample
    if (sample != "data_obs" or args.doData):
        print "doing this"
        hBDT.Write("BDT_%s_%s" % (catName,sample))
        print hBDT.Integral()
    else: print hBkg.Integral()
    print bdtname
    print "((%s)&&Pass_nominal)*%s" % (cutString,weight_string)
    if (sample == "TT" and args.ttbarTree != ""): 
        ifile_tt.Close()

ifile.Close()
ofile.Close()
otextfile.close()
   
