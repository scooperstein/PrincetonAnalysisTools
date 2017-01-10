import os

import ROOT
import numpy as np
import sys
import array

class BTagWeightCalculator:
    """
    Calculates the jet and event correction factor as a weight based on the b-tagger shape-dependent data/mc 
    corrections.

    Currently, the recipe is only described in https://twiki.cern.ch/twiki/bin/viewauth/CMS/TTbarHbbRun2ReferenceAnalysis#Applying_CSV_weights

    In short, jet-by-jet correction factors as a function of pt, eta and CSV have been derived.
    This code accesses the flavour of MC jets and gets the correct weight histogram
    corresponding to the pt, eta and flavour of the jet.
    From there, the per-jet weight is just accessed according to the value of the discriminator.
    """
    def __init__(self, fn_hf, fn_lf) :
        self.pdfs = {}

        #bin edges of the heavy-flavour histograms
        #pt>=20 && pt<30 -> bin=0
        #pt>=30 && pt<40 -> bin=1
        #etc
        self.pt_bins_hf = np.array([20, 30, 40, 60, 100])
        self.eta_bins_hf = np.array([0, 2.41])

        #bin edges of the light-flavour histograms 
        self.pt_bins_lf = np.array([20, 30, 40, 60])
        self.eta_bins_lf = np.array([0, 0.8, 1.6, 2.41])

        #name of the default b-tagger
        self.btag = "pfCombinedInclusiveSecondaryVertexV2BJetTags"
        self.init(fn_hf, fn_lf)

        self.systematics_for_b = ["JESUp", "JESDown", "LFUp", "LFDown",
                                  "HFStats1Up", "HFStats1Down", "HFStats2Up", "HFStats2Down"]
        self.systematics_for_c = ["cErr1Up", "cErr1Down", "cErr2Up", "cErr2Down"]
        self.systematics_for_l = ["JESUp", "JESDown", "HFUp", "HFDown",
                                  "LFStats1Up", "LFStats1Down", "LFStats2Up", "LFStats2Down"]

    def getBin(self, bvec, val):
        return int(bvec.searchsorted(val, side="right")) - 1

    def init(self, fn_hf, fn_lf):
        """
        fn_hf (string) - path to the heavy flavour weight file
        fn_lf (string) - path to the light flavour weight file
        """
        print "[BTagWeightCalculator]: Initializing from files", fn_hf, fn_lf

        #print "hf"
        self.pdfs["hf"] = self.getHistosFromFile(fn_hf)
        #print "lf"
        self.pdfs["lf"] = self.getHistosFromFile(fn_lf)

        return True

    def getHistosFromFile(self, fn):
        """
        Initialized the lookup table for b-tag weight histograms based on jet
        pt, eta and flavour.
        The format of the weight file is similar to:
         KEY: TH1D     csv_ratio_Pt0_Eta0_final;1
         KEY: TH1D     csv_ratio_Pt0_Eta0_final_JESUp;1
         KEY: TH1D     csv_ratio_Pt0_Eta0_final_JESDown;1
         KEY: TH1D     csv_ratio_Pt0_Eta0_final_LFUp;1
         KEY: TH1D     csv_ratio_Pt0_Eta0_final_LFDown;1
         KEY: TH1D     csv_ratio_Pt0_Eta0_final_Stats1Up;1
         KEY: TH1D     csv_ratio_Pt0_Eta0_final_Stats1Down;1
         KEY: TH1D     csv_ratio_Pt0_Eta0_final_Stats2Up;1
         KEY: TH1D     csv_ratio_Pt0_Eta0_final_Stats2Down;1
         KEY: TH1D     c_csv_ratio_Pt0_Eta0_final;2
         KEY: TH1D     c_csv_ratio_Pt0_Eta0_final;1
         KEY: TH1D     c_csv_ratio_Pt0_Eta0_final_cErr1Up;1
         KEY: TH1D     c_csv_ratio_Pt0_Eta0_final_cErr1Down;1
         KEY: TH1D     c_csv_ratio_Pt0_Eta0_final_cErr2Up;1
         KEY: TH1D     c_csv_ratio_Pt0_Eta0_final_cErr2Down;1
        """
        ret = {}
        tf = ROOT.TFile(fn)
        if not tf or tf.IsZombie():
            raise FileNotFoundError("Could not open file {0}".format(fn))
        ROOT.gROOT.cd()
        for k in tf.GetListOfKeys():
            kn = k.GetName()
            if not (kn.startswith("csv_ratio") or kn.startswith("c_csv_ratio") ):
                continue
            spl = kn.split("_")
            is_c = 1 if kn.startswith("c_csv_ratio") else 0

            if spl[2+is_c] == "all":
                ptbin = -1
                etabin = -1
                kind = "all"
                syst = "nominal"
            else:
                ptbin = int(spl[2+is_c][2:])
                etabin = int(spl[3+is_c][3:])
                kind = spl[4+is_c]
                if len(spl)==(6+is_c):
                    syst = spl[5+is_c]
                else:
                    syst = "nominal"
            ret[(is_c, ptbin, etabin, kind, syst)] = k.ReadObj().Clone()
        return ret

    def calcJetWeight(self, jet, kind, systematic):
        """
        Calculates the per-jet correction factor.
        jet: either an object with the attributes pt, eta, mcFlavour, self.btag
             or a Heppy Jet
        kind: string specifying the name of the corrections. Usually "final".
        systematic: the correction systematic, e.g. "nominal", "JESUp", etc
     """
        #if jet is a simple class with attributes
        if isinstance(getattr(jet, "pt"), float):
            pt   = getattr(jet, "pt")
            aeta = abs(getattr(jet, "eta"))
            fl   = abs(getattr(jet, "hadronFlavour"))
            csv  = getattr(jet, "btag")
        #if jet is a heppy Jet object
        else:
            #print "could not get jet", e
            pt   = jet.pt()
            aeta = abs(jet.eta())
            fl   = abs(jet.hadronFlavour())
            csv  = jet.btag(self.btag)
        return self.calcJetWeightImpl(pt, aeta, fl, csv, kind, systematic)

    def calcJetWeightImpl(self, pt, aeta, fl, csv, kind, systematic):

        is_b = (fl == 5)
        is_c = (fl == 4)
        is_l = not (is_b or is_c)

        if systematic != "nominal":
            if (is_b and systematic not in self.systematics_for_b) or (is_c and systematic not in self.systematics_for_c) or (is_l and systematic not in self.systematics_for_l):
                systematic = "nominal"

        if "Stats" in systematic:
            systematic = systematic[2:]

        if is_b or is_c:
            ptbin = self.getBin(self.pt_bins_hf, pt)
            etabin = self.getBin(self.eta_bins_hf, aeta)
        else:
            ptbin = self.getBin(self.pt_bins_lf, pt)
            etabin = self.getBin(self.eta_bins_lf, aeta)

        if ptbin < 0 or etabin < 0:
            #print "pt or eta bin outside range", pt, aeta, ptbin, etabin
            return 1.0

        k = (is_c, ptbin, etabin, kind, systematic)
        hdict = self.pdfs["lf"]
        if is_b or is_c:
            hdict = self.pdfs["hf"]
        h = hdict.get(k, None)
        if not h:
            #print "no histogram", k
            return 1.0

        if csv > 1:
            csv = 1
            
        csvbin = 1
        csvbin = h.FindBin(csv)
        #This is to fix csv=-10 not being accounted for in CSV SF input hists
        if csvbin <= 0:
            csvbin = 1
        if csvbin > h.GetNbinsX():
            csvbin = h.GetNbinsX()

        w = h.GetBinContent(csvbin)
        return w

    def calcEventWeight(self, jets, kind, systematic):
        """
        The per-event weight is just a product of per-jet weights.
        """
        weights = np.array(
            [self.calcJetWeight(jet, kind, systematic)
            for jet in jets]
        )

        wtot = np.prod(weights)
        return wtot
      

################################################################
# A dummy class of a jet object

class Jet :
    def __init__(self, pt, eta, fl, csv, csvname) :
        self.pt = pt
        self.eta = eta
        self.hadronFlavour = fl
        setattr(self, csvname, csv)

    #def hadronFlavour(self):
    #    return self.mcFlavour
    #def btag(self, name):
    #    return self.btagCSV
    #def pt(self):
    #    return self.pt
    #def eta(self):
    #    return self.eta


################################################################
#Set up offline b-weight calculation

#csvpath = os.environ['CMSSW_BASE']+"/src/VHbbAnalysis/Heppy/data/csv"
csvpath = "./"
bweightcalc = BTagWeightCalculator(
    #csvpath + "/csv_rwt_fit_hf_2015_12_14.root",
    #csvpath + "/csv_rwt_fit_lf_2015_12_14.root"
    csvpath + "/csv_rwt_fit_hf_76x_2016_02_08.root",
    csvpath + "/csv_rwt_fit_lf_76x_2016_02_08.root"
)
#bweightcalc.btag = "btagCSV"
bweightcalc.btag = "btag"

ifile = ROOT.TFile(sys.argv[1], "r")
tree = ifile.Get("tree")
settings = ifile.Get("settings")
#otree = tree.CopyTree(0)
#ofile = ROOT.TFile("output_wBTagWeightsNew.root", "RECREATE")
ofile = ROOT.TFile(sys.argv[2], "RECREATE")
ofile.cd()

#Jet_pt = np.zeros(100, dtype=float)
#Jet_eta = np.zeros(100, dtype=float)
#Jet_hadronFlavour = np.zeros(100, dtype=float)
#Jet_btagCSV = np.zeros(100, dtype=float)
#hJetInd1 = np.zeros(1, dtype=int)
#hJetInd2 = np.zeros(1, dtype=int)

#tree.SetBranchAddress("Jet_pt",Jet_pt)
#tree.SetBranchAddress("Jet_eta",Jet_eta)
#tree.SetBranchAddress("Jet_mcFlavour", Jet_hadronFlavour)
#tree.SetBranchAddress("Jet_btagCSV", Jet_btagCSV)
#tree.SetBranchAddress("hJetInd1", hJetInd1)
#tree.SetBranchAddress("hJetInd2", hJetInd2)

otree = tree.CloneTree(0)
osettings = settings.CloneTree()

bTagWeights = {}
bTagWeights["bTagWeight"] = np.zeros(1, dtype=float)
otree.Branch("bTagWeight", bTagWeights["bTagWeight"], "bTagWeight/D")
#bTagWeights["bTagWeightNew"] = np.zeros(1, dtype=float)
#otree.Branch("bTagWeightNew", bTagWeights["bTagWeightNew"], "bTagWeightNew/D")
for syst in ["JES", "LF", "HF", "LFStats1", "LFStats2", "HFStats1", "HFStats2", "cErr1", "cErr2"]:
    for sdir in ["Up", "Down"]:
        bTagWeights["bTagWeightNew_"+syst+sdir] = np.zeros(1, dtype=float)
        #otree.Branch("bTagWeightNew_"+syst+sdir, bTagWeights["bTagWeightNew_"+syst+sdir], "bTagWeightNew_"+syst+sdir+"/D")

## hack to add puWeight = 1.0 for data events
#puWeight = np.zeros(1)
#puWeight[0] = 1.0
#puWeight = array.array( 'f', [1.0] )
#puWeightUp = array.array( 'f', [1.0] )
#puWeightDown = array.array( 'f', [1.0] )
#otree.SetBranchAddress("puWeight",puWeight)
#otree.Branch("puWeight",puWeight,"puWeight/F")
#otree.Branch("puWeightUp",puWeight,"puWeight/F")
#otree.Branch("puWeightDown",puWeight,"puWeight/F")

#nentries = 1000
nentries = tree.GetEntries()
print "nentries = ",nentries
for entry in range(nentries):
    if (entry%10000 == 0): print "processing entry: %i" % entry
    tree.GetEntry(entry)
    if (tree.sampleIndex == 0):
        # data event
        bTagWeights["bTagWeight"][0] = 1.0
        #bTagWeights["bTagWeightNew"][0] = 1.0
        otree.Fill()
        continue
    #print Jet_pt, hJetInd1
    #print Jet_pt[hJetInd1], Jet_eta[hJetInd1], Jet_hadronFlavour[hJetInd1], Jet_btagCSV[hJetInd1], bweightcalc.btag
    #jet1 = Jet(Jet_pt[hJetInd1], Jet_eta[hJetInd1], Jet_hadronFlavour[hJetInd1], Jet_btagCSV[hJetInd1], bweightcalc.btag)
    #jet2 = Jet(Jet_pt[hJetInd2], Jet_eta[hJetInd2], Jet_hadronFlavour[hJetInd2], Jet_btagCSV[hJetInd2], bweightcalc.btag)
    
    jets = []
    for i in range(tree.nJet):
        if (tree.Jet_pt[i] > 30 and abs(tree.Jet_eta[i]) < 2.4): 
            jet = Jet(tree.Jet_pt[i], tree.Jet_eta[i], tree.Jet_hadronFlavour[i], tree.Jet_btagCSV[i], bweightcalc.btag)
            jets.append(jet)

    bTagWeights["bTagWeightNew"][0] = bweightcalc.calcEventWeight(
    jets, kind="final", systematic="nominal",
    )
    for syst in ["JES", "LF", "HF", "LFStats1", "LFStats2", "HFStats1", "HFStats2", "cErr1", "cErr2"]:
        for sdir in ["Up", "Down"]:
            bTagWeights["bTagWeightNew_"+syst+sdir][0] = bweightcalc.calcEventWeight(
                jets, kind="final", systematic=syst+sdir )
    #print bTagWeights
    otree.Fill()

ofile.cd()
otree.Write()
osettings.Write()
ofile.Close()
ifile.Close()

## EXAMPLE (1): per-jet nominal weight
#jet = Jet(50., 1.2, 5, 0.89, bweightcalc.btag)
#
#jet_weight_nominal = bweightcalc.calcJetWeight(
#    jet, kind="final", systematic="nominal",
#    )
##print "Nominal jet weight: ", jet_weight_nominal
#
## EXAMPLE (2): per-jet systematic up/down weight
#for syst in ["JES", "LF", "HF", "LFStats1", "LFStats2", "HFStats1", "HFStats2", "cErr1", "cErr2"]:
#    for sdir in ["Up", "Down"]:
#        jet_weight_shift = bweightcalc.calcJetWeight(
#            jet, kind="final", systematic=syst+sdir
#            )
#        print syst, sdir, ": ", jet_weight_shift
#
#
## EXAMPLE (3): the nominal event weight 
#jet1 = Jet(50., -1.2, 5, 0.99, bweightcalc.btag)
#jet2 = Jet(30., 1.8, 4, 0.2, bweightcalc.btag)
#jet3 = Jet(100., 2.2, 0, 0.1, bweightcalc.btag)
#jet4 = Jet(20., 0.5,-5, 0.6, bweightcalc.btag)
#jets = [jet1,jet2,jet3,jet4]
#
#event_weight_nominal = bweightcalc.calcEventWeight(
#    jets, kind="final", systematic="nominal",
#    )
#print "Nominal event weight: ", event_weight_nominal
#
#
## EXAMPLE (4): the systematic up/down event weight 
#for syst in ["JES", "LF", "HF", "LFStats1", "LFStats2", "HFStats1", "HFStats2", "cErr1", "cErr2"]:
#    for sdir in ["Up", "Down"]:
#        event_weight_shift = bweightcalc.calcEventWeight(
#            jets, kind="final", systematic=syst+sdir
#            )
#        print syst, sdir, ": ", event_weight_shift
