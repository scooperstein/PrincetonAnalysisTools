## code form https://github.com/vhbb/cmssw/blob/vhbbHeppy80X/VHbbAnalysis/Heppy/python/btagSF.py

import ROOT
import os
import sys
import numpy as np

debug_btagSF = False

# load the BTagCalibrationStandalone.cc macro from https://twiki.cern.ch/twiki/bin/view/CMS/BTagCalibration
#csvpath = os.environ['CMSSW_BASE']+"/src/VHbbAnalysis/Heppy/data/csv/"
ROOT.gSystem.Load("./BTagCalibrationStandalone.so")

# CSVv2
#calib_csv = ROOT.BTagCalibration("csvv2", "./CSVv2_ichep.csv")
calib_csv = ROOT.BTagCalibration("csvv2", "./ttH_BTV_CSVv2_13TeV_2016EF_7p2_2016_09_23.csv")

# cMVAv2
#calib_cmva = ROOT.BTagCalibration("cmvav2", "./cMVAv2_ichep.csv")
calib_cmva = ROOT.BTagCalibration("cmvav2", "./ttH_BTV_cMVAv2_13TeV_2016EF_7p2_2016_09_23.csv")

# map between algo/flavour and measurement type
sf_type_map = {
    "CSV" : {
        "file" : calib_csv,
        "bc" : "comb",
        "l" : "incl",
        },
    "CMVAV2" : {
        "file" : calib_cmva,
        "bc" : "ttbar",
        "l" : "incl",
        },
    }

# map of calibrators. E.g. btag_calibrators["CSVM_nominal_bc"], btag_calibrators["CSVM_up_l"], ...
btag_calibrators = {}
#for algo in ["CSV", "CMVAV2"]:
#    for wp in [ [0, "L"],[1, "M"], [2,"T"] ]:
#        for syst in ["central", "up", "down"]:
#            for fl in ["bc", "l"]:
#                print "[btagSF]: Loading calibrator for algo:", algo, ", WP:", wp[1], ", systematic:", syst, ", flavour:", fl
#                btag_calibrators[algo+wp[1]+"_"+syst+"_"+fl] = ROOT.BTagCalibrationReader(sf_type_map[algo]["file"], wp[0], sf_type_map[algo][fl], syst)

for algo in ["CSV", "CMVAV2"]:
    for syst in ["central", "up_jes", "down_jes", "up_lf", "down_lf", "up_hf", "down_hf", "up_hfstats1", "down_hfstats1", "up_hfstats2", "down_hfstats2", "up_lfstats1", "down_lfstats1", "up_lfstats2", "down_lfstats2", "up_cferr1", "down_cferr1", "up_cferr2", "down_cferr2"]:
        print "[btagSF]: Loading calibrator for algo:", algo, "systematic:", syst
        btag_calibrators[algo+"_iterative_"+syst] = ROOT.BTagCalibrationReader(sf_type_map[algo]["file"], 3 , "iterativefit", syst)

# depending on flavour, only a sample of systematics matter
def applies( flavour, syst ):
    if flavour==5 and syst not in ["central", "up_jes", "down_jes",  "up_lf", "down_lf",  "up_hfstats1", "down_hfstats1", "up_hfstats2", "down_hfstats2"]:
        return False
    elif flavour==4 and syst not in ["central", "up_cferr1", "down_cferr1", "up_cferr2", "down_cferr2" ]:
        return False
    elif flavour==0 and syst not in ["central", "up_jes", "down_jes", "up_hf", "down_hf",  "up_lfstats1", "down_lfstats1", "up_lfstats2", "down_lfstats2" ]:
        return False

    return True


# function that reads the SF
def get_SF(pt=30., eta=0.0, fl=5, val=0.0, syst="central", algo="CSV", wp="M", shape_corr=False, btag_calibrators=btag_calibrators):

    #print "get SF for jet with pt %f, eta %f, flavor %i, systematic %s" % (pt, eta, fl,syst)
    #print "val %f, algo %s, wp %s, shape_corr %i" % (val,algo,wp,int(shape_corr))
    # no SF for pT<20 GeV or pt>1000 or abs(eta)>2.4
    if abs(eta)>2.4 or pt>1000. or pt<20.:
        return 1.0

    # the .csv files use the convention: b=0, c=1, l=2. Convert into hadronFlavour convention: b=5, c=4, f=0
    fl_index = min(-fl+5,2)
    #print "fl_index = %i" % fl_index
    # no fl=1 in .csv for CMVAv2 (a bug???)
    if not shape_corr and "CMVAV2" in algo and fl==4:
        fl_index = 0

    if shape_corr:
        #print "passed shape_corr check"
        if applies(fl,syst):
            #print "passed applies(fl,syst) check"
            sf = btag_calibrators[algo+"_iterative_"+syst].eval(fl_index ,eta, pt, val)
            #print "going to return %f" % sf
            #print sf
            return sf
        else:
            sf = btag_calibrators[algo+"_iterative_central"].eval(fl_index ,eta, pt, val)
            #print sf
            #print "going to return %f" % sf
            return sf 

    # pt ranges for bc SF: needed to avoid out_of_range exceptions
    pt_range_high_bc = 670.-1e-02 if "CSV" in algo else 320.-1e-02
    pt_range_low_bc = 30.+1e-02

    # b or c jets
    if fl>=4:
        # use end_of_range values for pt in [20,30] or pt in [670,1000], with double error
        out_of_range = False
        if pt>pt_range_high_bc or pt<pt_range_low_bc:
            out_of_range = True        
        pt = min(pt, pt_range_high_bc)
        pt = max(pt, pt_range_low_bc)
        sf = btag_calibrators[algo+wp+"_"+syst+"_bc"].eval(fl_index ,eta, pt)
        # double the error for pt out-of-range
        if out_of_range and syst in ["up","down"]:
            sf = max(2*sf - btag_calibrators[algo+wp+"_central_bc"].eval(fl_index ,eta, pt), 0.)
        #print sf
        return sf
    # light jets
    else:
        sf = btag_calibrators[algo+wp+"_"+syst+"_l"].eval( fl_index ,eta, pt)
        #print sf
        return  sf

def get_event_SF(jets=[], syst="central", algo="CSV", btag_calibrators=btag_calibrators):
    weight = 1.0
    #print "called get_event_SF for syst %s" % syst
    for jet in jets:
        weight *= get_SF(pt=jet.pt, eta=jet.eta, fl=jet.hadronFlavour, val=jet.csv, syst=syst, algo=algo, wp="", shape_corr=True, btag_calibrators=btag_calibrators)
        #print "Jet with pt %f, eta %f, hadronFlavour %i gets weight %f" % (jet.pt, jet.eta, jet.hadronFlavour, (get_SF(pt=jet.pt, eta=jet.eta, fl=jet.hadronFlavour, val=jet.csv, syst=syst, algo=algo, wp="", shape_corr=True, btag_calibrators=btag_calibrators)))
    return weight                             

if debug_btagSF:
    #print "POG WP:"
    #for algo in ["CSV", "CMVAV2"]:
    #    for wp in [ "L", "M", "T" ]:
    #        print algo+wp+":"
    #        for syst in ["central", "up", "down"]:
    #            print "\t"+syst+":"
    #            for pt in [19.,25.,31.,330., 680.]:
    #                print ("\t\tB(pt=%.0f, eta=0.0): %.3f" % (pt, get_SF(pt=pt, eta=0.0, fl=5, val=0.0, syst=syst, algo=algo, wp=wp, shape_corr=False)))
    #                print ("\t\tC(pt=%.0f, eta=0.0): %.3f" % (pt, get_SF(pt=pt, eta=0.0, fl=4, val=0.0, syst=syst, algo=algo, wp=wp, shape_corr=False)))
    #                print ("\t\tL(pt=%.0f, eta=0.0): %.3f" % (pt, get_SF(pt=pt, eta=0.0, fl=0, val=0.0, syst=syst, algo=algo, wp=wp, shape_corr=False)))

    print "Iterative:"
    for algo in ["CSV"]:
    #for algo in ["CSV", "CMVAV2"]:
        print algo+":"
        for syst in ["central", "up_jes", "down_jes", "up_lf", "down_lf", "up_hf", "down_hf", "up_hfstats1", "down_hfstats1", "up_hfstats2", "down_hfstats2", "up_lfstats1", "down_lfstats1", "up_lfstats2", "down_lfstats2", "up_cferr1", "down_cferr1", "up_cferr2", "down_cferr2"]:
            print "\t"+syst+":"
            for pt in [50.]:
                print ("\t\tB(pt=%.0f, eta=0.0): %.3f" % (pt, get_SF(pt=pt, eta=0.0, fl=5, val=0.89, syst=syst, algo=algo, wp="", shape_corr=True)))
                print ("\t\tC(pt=%.0f, eta=0.0): %.3f" % (pt, get_SF(pt=pt, eta=0.0, fl=4, val=0.89, syst=syst, algo=algo, wp="", shape_corr=True)))
                print ("\t\tL(pt=%.0f, eta=0.0): %.3f" % (pt, get_SF(pt=pt, eta=0.0, fl=0, val=0.89, syst=syst, algo=algo, wp="", shape_corr=True)))

################################################################
# A dummy class of a jet object

class Jet :
    def __init__(self, pt, eta, fl, csv) :
        self.pt = pt
        self.eta = eta
        self.hadronFlavour = fl
        self.csv = csv

    #def hadronFlavour(self):
    #    return self.mcFlavour
    #def btag(self, name):
    #    return self.btagCSV
    #def pt(self):
    #    return self.pt
    #def eta(self):
    #    return self.eta


################################################################


ifile = ROOT.TFile(sys.argv[1], "r")
tree = ifile.Get("tree")
#settings = ifile.Get("settings")
#otree = tree.CopyTree(0)
#ofile = ROOT.TFile("output_wBTagWeightsNew.root", "RECREATE")
ofile = ROOT.TFile(sys.argv[2], "RECREATE")
ofile.cd()

Jet_pt = np.zeros(100, dtype=float)
Jet_eta = np.zeros(100, dtype=float)
Jet_hadronFlavour = np.zeros(100, dtype=float)
Jet_btagCSV = np.zeros(100, dtype=float)
hJetInd1 = np.zeros(1, dtype=int)
hJetInd2 = np.zeros(1, dtype=int)

#tree.SetBranchAddress("Jet_pt",Jet_pt)
#tree.SetBranchAddress("Jet_eta",Jet_eta)
#tree.SetBranchAddress("Jet_mcFlavour", Jet_hadronFlavour)
#tree.SetBranchAddress("Jet_btagCSV", Jet_btagCSV)
#tree.SetBranchAddress("hJetInd1", hJetInd1)
#tree.SetBranchAddress("hJetInd2", hJetInd2)

# since the SetBranchAddress call is not working
sysRefMap = {}

def MakeSysRefMap():
    sysRefMap["JESUp"] = tree.btagWeightCSV_up_jes
    sysRefMap["JESDown"] = tree.btagWeightCSV_down_jes
    sysRefMap["LFUp"] = tree.btagWeightCSV_up_lf
    sysRefMap["LFDown"] = tree.btagWeightCSV_down_lf
    sysRefMap["HFUp"] = tree.btagWeightCSV_up_hf
    sysRefMap["HFDown"] = tree.btagWeightCSV_down_hf
    sysRefMap["HFStats1Up"] = tree.btagWeightCSV_up_hfstats1
    sysRefMap["HFStats1Down"] = tree.btagWeightCSV_down_hfstats1
    sysRefMap["HFStats2Up"] = tree.btagWeightCSV_up_hfstats2
    sysRefMap["HFStats2Down"] = tree.btagWeightCSV_down_hfstats2
    sysRefMap["LFStats1Up"] = tree.btagWeightCSV_up_lfstats1
    sysRefMap["LFStats1Down"] = tree.btagWeightCSV_down_lfstats1
    sysRefMap["LFStats2Up"] = tree.btagWeightCSV_up_lfstats2
    sysRefMap["LFStats2Down"] = tree.btagWeightCSV_down_lfstats2
    sysRefMap["cErr1Up"] = tree.btagWeightCSV_up_cferr1
    sysRefMap["cErr1Down"] = tree.btagWeightCSV_down_cferr1
    sysRefMap["cErr2Up"] = tree.btagWeightCSV_up_cferr2
    sysRefMap["cErr2Down"] = tree.btagWeightCSV_down_cferr2

otree = tree.CloneTree(0)
#osettings = settings.CloneTree()

bTagWeights = {}
#bTagWeights["bTagWeight"] = np.zeros(1, dtype=float)
#otree.Branch("bTagWeight", bTagWeights["bTagWeight"], "bTagWeight/D")
bTagWeights["bTagWeightEF"] = np.zeros(1, dtype=float)
otree.Branch("bTagWeightEF", bTagWeights["bTagWeightEF"], "bTagWeightEF/D")
bTagWeights["bTagWeightBToG"] = np.zeros(1, dtype=float)
otree.Branch("bTagWeightBToG", bTagWeights["bTagWeightBToG"], "bTagWeightBToG/D")

# ICHEP weights already calculated in Heppy
bTagWeights["bTagWeight"] = np.zeros(1, dtype=float)
#tree.SetBranchAddress("bTagWeight",bTagWeights["bTagWeight"])

for syst in ["JES", "LF", "HF", "LFStats1", "LFStats2", "HFStats1", "HFStats2", "cErr1", "cErr2"]:
    for sdir in ["Up", "Down"]:
        bTagWeights["bTagWeightEF_"+syst+sdir] = np.zeros(1, dtype=float)
        bTagWeights["bTagWeightBToG_"+syst+sdir] = np.zeros(1, dtype=float)
        #bTagWeights["bTagWeight_"+syst+sdir] = np.zeros(1, dtype=float)
        otree.Branch("bTagWeightEF_"+syst+sdir, bTagWeights["bTagWeightEF_"+syst+sdir], "bTagWeightEF_"+syst+sdir+"/D")
        otree.Branch("bTagWeightBToG_"+syst+sdir, bTagWeights["bTagWeightBToG_"+syst+sdir], "bTagWeightBToG_"+syst+sdir+"/D")
        #tree.SetBranchAddress("bTagWeight_"+syst+sdir, bTagWeights["bTagWeight_"+syst+sdir])

## hack to add puWeight = 1.0 for data events
#puWeight = np.zeros(1)
#puWeight[0] = 1.0
#puWeight = array.array( 'f', [1.0] )
#puWeightUp = array.array( 'f', [1.0] )
#nentries = 100
nentries = tree.GetEntries()

# map from the systematic names we prefer to keep the same as before
# to the new names
sysMap = {}
sysMap["JESUp"] = "up_jes"
sysMap["JESDown"] = "down_jes"
sysMap["LFUp"] = "up_lf"
sysMap["LFDown"] = "down_lf"
sysMap["HFUp"] = "up_hf"
sysMap["HFDown"] = "down_hf"
sysMap["HFStats1Up"] = "up_hfstats1"
sysMap["HFStats1Down"] = "down_hfstats1"
sysMap["HFStats2Up"] = "up_hfstats2"
sysMap["HFStats2Down"] = "down_hfstats2"
sysMap["LFStats1Up"] = "up_lfstats1"
sysMap["LFStats1Down"] = "down_lfstats1"
sysMap["LFStats2Up"] = "up_lfstats2"
sysMap["LFStats2Down"] = "down_lfstats2"
sysMap["cErr1Up"] = "up_cferr1"
sysMap["cErr1Down"] = "down_cferr1"
sysMap["cErr2Up"] = "up_cferr2"
sysMap["cErr2Down"] = "down_cferr2"

print "nentries = ",nentries
for entry in range(nentries):
    if (entry%10000 == 0): print "processing entry: %i" % entry
    #print "processing entry: %i" % entry
    tree.GetEntry(entry)
    MakeSysRefMap()
    #if (tree.sampleIndex == 0):
    #    # data event
    #    #bTagWeights["bTagWeight"][0] = 1.0
    #    bTagWeights["bTagWeightEF"][0] = 1.0
    #    otree.Fill()
    #    continue
    #print Jet_pt, hJetInd1
    #print Jet_pt[hJetInd1], Jet_eta[hJetInd1], Jet_hadronFlavour[hJetInd1], Jet_btagCSV[hJetInd1], bweightcalc.btag
    #jet1 = Jet(Jet_pt[hJetInd1], Jet_eta[hJetInd1], Jet_hadronFlavour[hJetInd1], Jet_btagCSV[hJetInd1], bweightcalc.btag)
    #jet2 = Jet(Jet_pt[hJetInd2], Jet_eta[hJetInd2], Jet_hadronFlavour[hJetInd2], Jet_btagCSV[hJetInd2], bweightcalc.btag)
    
    jets = []
    for i in range(tree.nJet):
        if (tree.Jet_pt[i] > 25 and abs(tree.Jet_eta[i]) < 2.4): 
            jet = Jet(tree.Jet_pt[i], tree.Jet_eta[i], tree.Jet_hadronFlavour[i], tree.Jet_btagCSV[i])
            jets.append(jet)

    bTagWeights["bTagWeightEF"][0] = get_event_SF( jets, "central", "CSV", btag_calibrators)
    bTagWeights["bTagWeightBToG"][0] = (12.9/22.0)*tree.btagWeightCSV + (9.1/22.0)*bTagWeights["bTagWeightEF"][0] 
    for syst in ["JES", "LF", "HF", "LFStats1", "LFStats2", "HFStats1", "HFStats2", "cErr1", "cErr2"]:
        for sdir in ["Up", "Down"]:
            bTagWeights["bTagWeightEF_"+syst+sdir][0] = get_event_SF( jets, sysMap[syst+sdir], "CSV", btag_calibrators)
            bTagWeights["bTagWeightBToG_"+syst+sdir][0] = (12.9/22.0)*sysRefMap[syst+sdir] + (9.1/22.0)*bTagWeights["bTagWeightEF_"+syst+sdir][0]
            #bTagWeights["bTagWeightEF_"+syst+sdir][0] = bweightcalc.calcEventWeight(
            #    jets, kind="final", systematic=syst+sdir )
    #print bTagWeights
    #for bTagWeight in bTagWeights:
    #    print "%s: %f" % (bTagWeight,bTagWeights[bTagWeight])
    otree.Fill()

ofile.cd()
otree.Write()
#osettings.Write()
ofile.Close()
ifile.Close()
