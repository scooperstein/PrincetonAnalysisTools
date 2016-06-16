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
parser.add_argument('-d', '--doData', type=bool, default=False, help="If true run on real data, if false compute data as sum of background pdf's")
parser.add_argument('-t', '--dataTree', type=str, default="", help="If doData is true, you specify here the ntuple with the real data events")
parser.add_argument('-tt', '--ttbarTree', type=str, default="", help="Specify separate file with ttbar powheg so it only takes forever when running tt")
parser.add_argument('-v', '--varname', type=str, default="CMS_vhbb_BDT_Wln_13TeV", help="The name of the variable shape which goes into the histograms and will be fitted")
parser.add_argument('-xl','--xlow', type=float, default=-1.0, help="Lowest bin edge of fitted distribution")
parser.add_argument('-xh','--xhigh', type=float, default=1.0, help="Highest bin edge of fitted distribution")
parser.add_argument('-bh','--binhigh',type=float,default=1.0, help="manually set low bin edge of most sensitive bin")
parser.add_argument('-bl','--binlow',type=float,default=-1.0, help="manually set low bin edge of least sensitive bin")
parser.add_argument('-o', '--ofilename',type=str,default="", help="output histogram file name (default hists_[cat].root)")
parser.add_argument('-b', '--binstats', type=str, default="", help="Text file listing all the individual bin. stat. uncertainties to include")
parser.add_argument('-n', '--nbins', type=int, default=20, help="number of bins in datacard shapes")
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
nBinsFine = 100 # number of candidate bin edges to begin with
tolerance = 1.5 # dB/B tolerance for bin-by-bin stat. uncertainties
#tolerance = 0.5 # dB/B tolerance for bin-by-bin stat. uncertainties

sampleMap = {} # map sampleNames to list of sampleIndex's
sampleMapAltModel = {} # alternate MC samples for model shape systematics

#sampleMap["data_obs"] = [16,17,10,11,20,21,50,51,52,2200,2201,2202,2300,2301,2302,3500,3501,3502,4400,4401,4402,4500,4501,4502,4600,4601,4602,4700,4701,4702,24,25,26,27,28,29,30,31] # dummy filler until we unblind
#sampleMap["data_obs"] = [50,51,52]
sampleMap["data_obs"] = [0]
sampleMapAltModel["TT"] = [50,51,52]
sampleMap["TT"] = [120]
sampleMap["s_Top"] = [16,17,20,21]
sampleMap["WH"] = [-12501]
sampleMapAltModel["WH"] = [-125010, -125011]
sampleMap["ZH"] = [-12502]
sampleMap["Wj0b"] = [2200,4400,4500,4600,4700,4800,4900]
sampleMapAltModel["Wj0b"] = [6000]
sampleMap["Wj1b"] = [2201,4401,4501,4601,4701,4801,4901]
sampleMapAltModel["Wj1b"] = [6001]
sampleMap["Wj2b"] = [2202,4402,4502,4602,4702,4802,4902]
sampleMapAltModel["Wj2b"] = [6002]
sampleMap["VVHF"] = [3501,3502,3601,3602,3701,3702]
sampleMap["VVLF"] = [3500,3600,3700]
sampleMap["QCD"]  = [24,25,26,27,28,29,30,31]
sampleMap["Zj0b"] = [2300,6100,6200,6300,6400]
sampleMap["Zj1b"] = [2301,6101,6201,6301,6401]
sampleMap["Zj2b"] = [2302,6102,6202,6302,6402]
#sampleMapAltModel["Zj0b"] = [23000]
#sampleMapAltModel["Zj1b"] = [23001]
#sampleMapAltModel["Zj2b"] = [23002]

allSampInd = [] # list of all indices for all backgrounds
for sample in sampleMap:
    if (sample == "WH" or sample == "ZH" or sample == "data_obs"): continue
    allSampInd.extend(sampleMap[sample])
sampleMap["Bkg"] = allSampInd
if not args.doData:
    sampleMap["data_obs"] = allSampInd # fake data as sum of all background MC

def makeCutString(sample, sMap):
    cutString = presel + "&&("
    for index in sMap[sample]:
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
tree_mc.Draw("%s>>hBkg" % bdtname,"((%s)&&(%s))*bTagWeight*%s" % (presel,bkgCutString,weight_string))
tree_mc.Draw("%s>>hSig" % bdtname,"((%s)&&sampleIndex==-12501)*bTagWeight*%s" % (presel,weight_string))


if (args.ttbarTree != ""):
    print "got here"
    ifile_tt = ROOT.TFile(args.ttbarTree,"r")
    hBkg2 = ROOT.TH1F("hBkg2","hBkg2",nBinsFine,args.xlow,args.xhigh)
    tree_tt = ifile_tt.Get("tree")
    # make sure we don't weight actual data by puWeight, SF's, etc.
    #print tree_tt.GetEntries("sampleIndex==120")
    tree_tt.Draw("%s>>hBkg2" % bdtname,"((%s)&&(%s))*bTagWeight*%s" % (presel,bkgCutString,weight_string))
    print hBkg.Integral()
    print hBkg2.Integral()
    hBkg.Add(hBkg2)
    print hBkg.Integral()
    ifile_tt.Close() 

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
for ibin in range(2,nBinsFine):
    if not foundLowBinEdge:
        B_low = hBkg.Integral(1,ibin)
        #S_low = hSig.Integral(1,ibin)
        #print "edge = ",hBkg.GetBinLowEdge(ibin+1)
        #print "B_low = ",B_low
        #if (B_low > 0):
        #    print "1./sqrt(B_low) = ",1./sqrt(B_low)
        #print B_low,ibin,hBkg.Integral(),args.xlow,args.xhigh
        if (B_low > 0 and 1./sqrt(B_low) < 0.35):
            binBoundaries[1] = hBkg.GetBinLowEdge(ibin+1)
            foundLowBinEdge = True
    if not foundHighBinEdge:
        B_high = hBkg.Integral(nBinsFine-ibin+1, nBinsFine) 
        #S_high = hSig.Integral(nBinsFine-ibin+1, nBinsFine) 
        if (B_high > 0 and 1./sqrt(B_high) < 0.35):
            binBoundaries[nBins-1] = hBkg.GetBinLowEdge(nBinsFine-ibin+1)
            foundHighBinEdge = True
# split the middle bins equidistantly
for i in range(2,nBins-1):
    binBoundaries[i] = binBoundaries[1] + (i-1)*((binBoundaries[nBins-1] - binBoundaries[1])/(nBins-2)) 

### FIXME get rid of rebinning with this, should put it back in, it's just a test
#for i in range(1,nBins):
#    binBoundaries[i] = binBoundaries[0] + (i)*((binBoundaries[nBins] - binBoundaries[0])/(nBins))
hBkg = hBkg.Rebin(nBins, "", binBoundaries)
print binBoundaries
if args.ofilename == "":
    ofile = ROOT.TFile("hists_%s.root" % catName, "RECREATE")
else:
    ofile = ROOT.TFile(args.ofilename, "RECREATE")
if args.binstats == "":
    otextfile = open("binStats_%s.txt" % catName, "w")
else:
    otextfile = open(args.binstats,"w")
tree = ROOT.TTree("tree","tree")
#hBkg.Write()
for sample in sampleMap:
    if (sample == "Bkg"): continue
    if (sample == "TT" and args.ttbarTree != ""): 
        ifile_tt = ROOT.TFile(args.ttbarTree,"r")
        tree_tt = ifile_tt.Get("tree")
        tree = tree_tt
        #ofile.cd()
        #ifile_tt.Close()
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
        ifile_data = ROOT.TFile(args.dataTree,"r")
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
    #    ifile_tt = ROOT.TFile(args.ttbarTree,"r")
    #    tree_tt = ifile_tt.Get("tree")
    #    ofile.cd()
    #    # make sure we don't weight actual data by puWeight, SF's, etc.
    #    tree_tt.Draw("%s>>%s" % (bdtname, sample),"((%s)&&Pass_nominal)*bTagWeight*%s" % (cutString,weight_string)) 
    #    ifile_tt.Close() 
    else:
        tree.Draw("%s>>%s" % (bdtname, sample),"((%s)&&Pass_nominal)*bTagWeight*%s" % (cutString,weight_string)) 
    print "tree.Draw(\"%s>>%s\",\"((%s)&&Pass_nominal)*bTagWeight*%s\"" % (bdtname,sample,cutString,weight_string) 
    #hBDT = hBDT.Rebin(nBins, "", binBoundaries)
    # Add bin-by-bin stat. uncertainties
    if (sample not in ["WH","ZH","data_obs"]):
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
            #hBDT.GetXaxis().SetRange(ibin,ibin)
            #B_eff = hBDT.GetEffectiveEntries() # effective statistical number of entries in bin considering weights
            if ( B > 0 and ( ( B >=1 and (B_err/sqrt(B)) > tolerance) or (B < 1 and B_err/B > tolerance) ) ):
                print "qualified as a bin stat. shape uncertainty! Bin %i, sample %s" % (ibin,sample)
                if (B >=1): print B_err/sqrt(B)
                else: print B_err/B
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
        hBDTSystUp = ROOT.TH1F("BDT_%s_%s_%sUp" % (catName,sample,syst), "%s_%sUp" % (sample,syst),nBins,binBoundaries)
        hBDTSystDown = ROOT.TH1F("BDT_%s_%s_%sDown" % (catName,sample,syst), "%s_%sDown" % (sample,syst),nBins,binBoundaries)
        sysBDTNameUp = bdtname
        sysBDTNameDown = bdtname
        passSys = "Pass_nominal"
        if (sysName != ""):
            sysBDTNameUp += "_%sUp" % sysName
            sysBDTNameDown += "_%sDown" % sysName
            pasSys = "Pass_%s" % sysName
        if (sysWeight != "1.0"):
            #tree.Draw("%s>>%s_%sUp" % (sysBDTNameUp, sample, syst),"(%s)*weight*(2.2/1.28)*(%sUp)" % (cutString,sysWeight))   
            tree.Draw("%s>>BDT_%s_%s_%sUp" % (sysBDTNameUp, catName,sample, syst),"((%s)&&%s)*%s*(%sUp)" % (cutString,passSys,weight_string,sysWeight))   
            #tree.Draw("%s>>%s_%sDown" % (sysBDTNameDown, sample, syst),"(%s)*weight*(2.2/1.28)*(%sDown)" % (cutString,sysWeight))   
            tree.Draw("%s>>BDT_%s_%s_%sDown" % (sysBDTNameDown, catName, sample, syst),"((%s)&&%s)*%s*(%sDown)" % (cutString,passSys,weight_string,sysWeight))   
        else:
            #tree.Draw("%s>>%s_%sUp" % (sysBDTNameUp, sample, syst),"(%s)*weight*(2.2/1.28)" % (cutString))   
            #tree.Draw("%s>>%s_%sUp" % (sysBDTNameUp, sample, syst),"((%s)&&%s)*%s*bTagWeight" % (cutString,passSys,weight_string))   
            tree.Draw("%s>>BDT_%s_%s_%sUp" % (sysBDTNameUp, catName,sample, syst),"((%s)&&%s)*%s*bTagWeight" % (cutString,passSys,weight_string))   
            #tree.Draw("%s>>%s_%sDown" % (sysBDTNameDown, sample, syst),"(%s)*weight*(2.2/1.28)" % (cutString))   
            #tree.Draw("%s>>%s_%sDown" % (sysBDTNameDown, sample, syst),"((%s)&&%s)*%s*bTagWeight" % (cutString,passSys,weight_string))   
            tree.Draw("%s>>BDT_%s_%s_%sDown" % (sysBDTNameDown, catName,sample, syst),"((%s)&&%s)*%s*bTagWeight" % (cutString,passSys,weight_string))   
            if (syst.find("Model") != -1):
                # shape from different sample (amc@NLO vs. POWHEG, for example)
                cutStringAltModel = makeCutString(sample,sampleMapAltModel)
                hBDTAltModel = ROOT.TH1F("%s_%sAltModel" % (sample,syst), "%s_%sAltModel" % (sample,syst),nBins,binBoundaries)
                tree.Draw("%s>>%s_%sAltModel" % (sysBDTNameUp, sample, syst),"((%s)&&%s)*%s*bTagWeight" % (cutStringAltModel,passSys,weight_string))
                # hnom + (hnom - halt) = 2*hnom - halt
                hBDTSystUp.Scale(2.0)
                hBDTSystUp.Add(hBDTAltModel, -1.0)
                # hnom - (hnom - halt) = halt
                hBDTSystDown = hBDTAltModel
                hBDTSystDown.SetName("BDT_%s_%s_%sDown" % (catName,sample, syst))
            
        ofile.cd()
        hBDTSystUp.Write()
        hBDTSystDown.Write()
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
    print "((%s)&&Pass_nominal)*bTagWeight*%s" % (cutString,weight_string)
    if (sample == "TT" and args.ttbarTree != ""): 
        ifile_tt.Close()

ifile.Close()
ofile.Close()
otextfile.close()
   
