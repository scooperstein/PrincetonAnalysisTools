import ROOT
import sys
import numpy as np

ROOT.gROOT.SetBatch(True)

filename=sys.argv[1]
tfile=ROOT.TFile(filename, "r")

treename="tree"
tree=tfile.Get(treename)

ofile = ROOT.TFile(sys.argv[2], "RECREATE")
otree = tree.CloneTree(0)

varsmap = {}
def AssignVarsToMap():
    varsmap["jet1csv"]        = tree.hJets_btagCSV_0
    varsmap["jet2csv"]        = tree.hJets_btagCSV_1
    varsmap["jet1pt"]         = tree.hJets_pt_0
    varsmap["jet2pt"]         = tree.hJets_pt_1
    varsmap["ptjj"]           = tree.H_pt
    varsmap["ptW"]            = tree.V_pt
    varsmap["met"]            = tree.met_pt
    varsmap["HVdPhi"]         = tree.HVdPhi
    varsmap["lepmetdphi"]     = tree.lepMetDPhi
    #varsmap["naddlep"]        = tree.nAddLeptons
    varsmap["naddjet"]        = tree.nAddJets252p9_puid
    varsmap["topmass"]        = tree.Top1_mass_fromLepton_regPT_w4MET
    varsmap["isWenu"]         = tree.isWenu
    varsmap["isWmunu"]        = tree.isWmunu
    varsmap["type"]           = tree.sampleIndex
    varsmap["Mjj"]            = tree.H_mass
    varsmap["ptlep"]          = tree.selLeptons_pt_0
    varsmap["etalep"]         = tree.selLeptons_eta_0
    varsmap["relisolep"]      = tree.selLeptons_relIso_0
    varsmap["drjj"]         = tree.HJ1_HJ2_dR
    varsmap["hwptbal"]        = tree.jjWPtBalance
    varsmap["ptaddjet"]      = tree.AddJets252p9_puid_leadJet_pt
    varsmap["mtW"]            = tree.V_mt
    varsmap["nsa5jet"]        = tree.softActivityVH_njets5

cutsets={}
cutsets["presel"] = {}
cutsets["presel"]["jet1csv"]       =   [0.935, 0.935]
cutsets["presel"]["jet2csv"]       =   [0.46, 0.46]
cutsets["presel"]["jet1pt"]        =   [25, 25]
cutsets["presel"]["jet2pt"]        =   [25, 25]
cutsets["presel"]["ptjj"]          =   [100, 100]
cutsets["presel"]["ptW"]           =   [100, 100]
cutsets["presel"]["met"]           =   [0, 0]
cutsets["presel"]["HVdPhi"]        =   [2.5, 2.5]
cutsets["presel"]["lepmetdphi"]    =   [2.0, 2.0]
#cutsets["presel"]["naddlep"]       =   [1, 1]
cutsets["presel"]["naddjet"]       =   [2, 2]
cutsets["presel"]["topmass"] =   [100, 100]
cutsets["presel"]["ptlep"]         =   [25, 25]
cutsets["presel"]["etalep"]       =   [2.5,2.5]
cutsets["presel"]["relisolep"]    =   [0.15, 0.12]
cutsets["presel"]["drjj"]         =   [3.5, 3.5]
#cutsets["presel"]["hwptbal"]      =   [3,3]
cutsets["presel"]["ptaddjet"]      =   [20,20]
cutsets["presel"]["mtW"]          = [200,200]
cutsets["presel"]["nsa5jet"]      = [10, 10]

cutsets["Analysis3Em3"]={}
cutsets["Analysis3Em3"]["etalep"]=[ 2.5,  2.5]
cutsets["Analysis3Em3"]["jet2pt"]=[ 25.0000,  25.0000]
cutsets["Analysis3Em3"]["lepmetdphi"]=[ 1.90,  1.90]
cutsets["Analysis3Em3"]["jet2csv"]=[ 0.46,  0.46]
cutsets["Analysis3Em3"]["jet1pt"]=[ 25.0000,  25.0000]
cutsets["Analysis3Em3"]["mtW"]=[ 200.0000,  200]
cutsets["Analysis3Em3"]["HVdPhi"]=[ 2.85,  2.85]
cutsets["Analysis3Em3"]["met"]=[ 0,  0]
cutsets["Analysis3Em3"]["topmass"]=[ 100,  100]
cutsets["Analysis3Em3"]["relisolep"]=[ 0.15,  0.12]
cutsets["Analysis3Em3"]["ptW"]=[ 115,  115.0000]
cutsets["Analysis3Em3"]["jet1csv"]=[ 0.935,  0.935]
cutsets["Analysis3Em3"]["ptaddjet"]=[ 20,  20]
cutsets["Analysis3Em3"]["naddjet"]=[ 2.0000,  2.0000]
cutsets["Analysis3Em3"]["ptjj"]=[ 110,  110]
cutsets["Analysis3Em3"]["ptlep"]=[ 25,  25]
cutsets["Analysis3Em3"]["drjj"]=[ 3.5,  3.5]
cutsets["Analysis3Em3"]["nsa5jet"]=[ 10.0000,  10.0000]

cutsets["Analysis3Em3p5"]={}
cutsets["Analysis3Em3p5"]["etalep"]=[ 2.5,  2.5]
cutsets["Analysis3Em3p5"]["jet2pt"]=[ 25.0000,  25.0000]
cutsets["Analysis3Em3p5"]["lepmetdphi"]=[ 1.85,  1.85]
cutsets["Analysis3Em3p5"]["jet2csv"]=[ 0.52,  0.52]
cutsets["Analysis3Em3p5"]["jet1pt"]=[ 25.0000,  25.0000]
cutsets["Analysis3Em3p5"]["mtW"]=[ 200.0000,  200]
cutsets["Analysis3Em3p5"]["HVdPhi"]=[ 2.91,  2.91]
cutsets["Analysis3Em3p5"]["met"]=[ 0,  0]
cutsets["Analysis3Em3p5"]["topmass"]=[ 100.0000,  100]
cutsets["Analysis3Em3p5"]["relisolep"]=[ 0.1500,  0.12]
cutsets["Analysis3Em3p5"]["ptW"]=[ 117,  117]
cutsets["Analysis3Em3p5"]["jet1csv"]=[ 0.935,  0.935]
cutsets["Analysis3Em3p5"]["ptaddjet"]=[ 20.0000,  20.0000]
cutsets["Analysis3Em3p5"]["naddjet"]=[ 2.0000,  2.0000]
cutsets["Analysis3Em3p5"]["ptjj"]=[ 115,  115]
cutsets["Analysis3Em3p5"]["ptlep"]=[ 25,  25.0000]
cutsets["Analysis3Em3p5"]["drjj"]=[ 3.5,  3.5]
cutsets["Analysis3Em3p5"]["nsa5jet"]=[ 10.0000,  10.0000]

cutsets["Analysis4Em3"]={}
cutsets["Analysis4Em3"]["etalep"]=[ 2.5,  2.5]
cutsets["Analysis4Em3"]["jet2pt"]=[ 25.0000,  25.0000]
cutsets["Analysis4Em3"]["lepmetdphi"]=[ 1.8500,  1.8500]
cutsets["Analysis4Em3"]["jet2csv"]=[ 0.54,  0.54]
cutsets["Analysis4Em3"]["jet1pt"]=[ 25.0000,  25.0000]
cutsets["Analysis4Em3"]["mtW"]=[ 200.0000,  200]
cutsets["Analysis4Em3"]["HVdPhi"]=[ 2.93,  2.93]
cutsets["Analysis4Em3"]["met"]=[ 0,  0]
cutsets["Analysis4Em3"]["topmass"]=[ 100.0000,  100]
cutsets["Analysis4Em3"]["relisolep"]=[ 0.1500,  0.12]
cutsets["Analysis4Em3"]["ptW"]=[ 120,  120]
cutsets["Analysis4Em3"]["jet1csv"]=[ 0.935,  0.935]
cutsets["Analysis4Em3"]["ptaddjet"]=[ 20.0000,  20.0000]
cutsets["Analysis4Em3"]["naddjet"]=[ 2.0000,  2.0000]
cutsets["Analysis4Em3"]["ptjj"]=[ 117,  117]
cutsets["Analysis4Em3"]["ptlep"]=[ 25,  25.0000]
cutsets["Analysis4Em3"]["drjj"]=[ 2.9,  2.9]
cutsets["Analysis4Em3"]["nsa5jet"]=[ 10,  10]

cutsets["Analysis5Em3"]={}
cutsets["Analysis5Em3"]["etalep"]=[ 2.5,  2.5]
cutsets["Analysis5Em3"]["jet2pt"]=[ 25.0000,  25.0000]
cutsets["Analysis5Em3"]["lepmetdphi"]=[ 1.8500,  1.8500]
cutsets["Analysis5Em3"]["jet2csv"]=[ 0.56,  0.56]
cutsets["Analysis5Em3"]["jet1pt"]=[ 25.0000,  25.0000]
cutsets["Analysis5Em3"]["mtW"]=[ 200.0000,  200]
cutsets["Analysis5Em3"]["HVdPhi"]=[ 2.93,  2.93]
cutsets["Analysis5Em3"]["met"]=[ 0,  0.0000]
cutsets["Analysis5Em3"]["topmass"]=[ 100,  100]
cutsets["Analysis5Em3"]["relisolep"]=[ 0.1500,  0.12]
cutsets["Analysis5Em3"]["ptW"]=[ 125,  125]
cutsets["Analysis5Em3"]["jet1csv"]=[ 0.935,  0.935]
cutsets["Analysis5Em3"]["ptaddjet"]=[ 20.0000,  20.0000]
cutsets["Analysis5Em3"]["naddjet"]=[ 2.0000,  2.0000]
cutsets["Analysis5Em3"]["ptjj"]=[ 125,  125]
cutsets["Analysis5Em3"]["ptlep"]=[ 25.0000,  25.0000]
cutsets["Analysis5Em3"]["drjj"]=[ 2.9,  2.9]
cutsets["Analysis5Em3"]["nsa5jet"]=[ 10.0000,  10.0000]

cutsets["Analysis6Em3"]={}
cutsets["Analysis6Em3"]["etalep"]=[ 2.5,  2.5]
cutsets["Analysis6Em3"]["jet2pt"]=[ 25.0000,  25.0000]
cutsets["Analysis6Em3"]["lepmetdphi"]=[ 1.8500,  1.85]
cutsets["Analysis6Em3"]["jet2csv"]=[ 0.60,  0.60]
cutsets["Analysis6Em3"]["jet1pt"]=[ 25.0000,  25.0000]
cutsets["Analysis6Em3"]["mtW"]=[ 200.0000,  200]
cutsets["Analysis6Em3"]["HVdPhi"]=[ 2.94,  2.94]
cutsets["Analysis6Em3"]["met"]=[ 0,  0]
cutsets["Analysis6Em3"]["topmass"]=[ 132,  132]
cutsets["Analysis6Em3"]["relisolep"]=[ 0.1500,  0.12]
cutsets["Analysis6Em3"]["ptW"]=[ 135,  135]
cutsets["Analysis6Em3"]["jet1csv"]=[ 0.935,  0.935]
cutsets["Analysis6Em3"]["ptaddjet"]=[ 20.0000,  20.0000]
cutsets["Analysis6Em3"]["naddjet"]=[ 2.0000,  2.0000]
cutsets["Analysis6Em3"]["ptjj"]=[ 138,  138]
cutsets["Analysis6Em3"]["ptlep"]=[ 25.0000,  25]
cutsets["Analysis6Em3"]["drjj"]=[ 2.3,  2.3]
cutsets["Analysis6Em3"]["nsa5jet"]=[ 9.0000,  9.0000]

cutsets["Analysis8Em3"]={}
cutsets["Analysis8Em3"]["etalep"]=[ 2.5,  2.5]
cutsets["Analysis8Em3"]["jet2pt"]=[ 25.0000,  25.0000]
cutsets["Analysis8Em3"]["lepmetdphi"]=[ 1.8500,  1.85]
cutsets["Analysis8Em3"]["jet2csv"]=[ 0.62,  0.62]
cutsets["Analysis8Em3"]["jet1pt"]=[ 25.0000,  25.0000]
cutsets["Analysis8Em3"]["mtW"]=[ 200.0000,  200]
cutsets["Analysis8Em3"]["HVdPhi"]=[ 2.95,  2.95]
cutsets["Analysis8Em3"]["met"]=[ 0,  0]
cutsets["Analysis8Em3"]["topmass"]=[ 190,  190]
cutsets["Analysis8Em3"]["relisolep"]=[ 0.15,  0.12]
cutsets["Analysis8Em3"]["ptW"]=[ 139,  139]
cutsets["Analysis8Em3"]["jet1csv"]=[ 0.935,  0.935]
cutsets["Analysis8Em3"]["ptaddjet"]=[ 20,  20]
cutsets["Analysis8Em3"]["naddjet"]=[ 2.0000,  2.0000]
cutsets["Analysis8Em3"]["ptjj"]=[ 141,  141]
cutsets["Analysis8Em3"]["ptlep"]=[ 25,  25.0000]
cutsets["Analysis8Em3"]["drjj"]=[ 1.52,  1.52]
cutsets["Analysis8Em3"]["nsa5jet"]=[ 4.0000,  4.0000]

cutsets["Analysis1Em2"]={}
cutsets["Analysis1Em2"]["etalep"]=[ 2.5,  2.5]
cutsets["Analysis1Em2"]["jet2pt"]=[ 25.0000,  25.0000]
cutsets["Analysis1Em2"]["lepmetdphi"]=[ 1.8500,  1.85]
cutsets["Analysis1Em2"]["jet2csv"]=[ 0.63,  0.63]
cutsets["Analysis1Em2"]["jet1pt"]=[ 25.0000,  25.0000]
cutsets["Analysis1Em2"]["mtW"]=[ 200.0000,  200]
cutsets["Analysis1Em2"]["HVdPhi"]=[ 2.96,  2.96]
cutsets["Analysis1Em2"]["met"]=[ 0,  0.0000]
cutsets["Analysis1Em2"]["topmass"]=[ 190,  190.0000]
cutsets["Analysis1Em2"]["relisolep"]=[ 0.15,  0.12]
cutsets["Analysis1Em2"]["ptW"]=[ 140,  140]
cutsets["Analysis1Em2"]["jet1csv"]=[ 0.935,  0.935]
cutsets["Analysis1Em2"]["ptaddjet"]=[ 20.0000,  20.0000]
cutsets["Analysis1Em2"]["naddjet"]=[ 2.0000,  2.0000]
cutsets["Analysis1Em2"]["ptjj"]=[ 145,  145]
cutsets["Analysis1Em2"]["ptlep"]=[ 25.0000,  25.0000]
cutsets["Analysis1Em2"]["drjj"]=[ 1.52,  1.52]
cutsets["Analysis1Em2"]["nsa5jet"]=[ 3.0000,  3.0000]

cutsets["Analysis1p5Em2"]={}
cutsets["Analysis1p5Em2"]["etalep"]=[ 2.5,  2.5]
cutsets["Analysis1p5Em2"]["jet2pt"]=[ 25.0000,  25.0000]
cutsets["Analysis1p5Em2"]["lepmetdphi"]=[ 1.85,  1.85]
cutsets["Analysis1p5Em2"]["jet2csv"]=[ 0.66,  0.66]
cutsets["Analysis1p5Em2"]["jet1pt"]=[ 25.0000,  25.0000]
cutsets["Analysis1p5Em2"]["mtW"]=[ 200.0000,  200]
cutsets["Analysis1p5Em2"]["HVdPhi"]=[ 2.97,  2.97]
cutsets["Analysis1p5Em2"]["met"]=[ 0,  0]
cutsets["Analysis1p5Em2"]["topmass"]=[ 205,  205]
cutsets["Analysis1p5Em2"]["relisolep"]=[ 0.15,  0.12]
cutsets["Analysis1p5Em2"]["ptW"]=[ 148,  148]
cutsets["Analysis1p5Em2"]["jet1csv"]=[ 0.935,  0.935]
cutsets["Analysis1p5Em2"]["ptaddjet"]=[ 20.0000,  20.0000]
cutsets["Analysis1p5Em2"]["naddjet"]=[ 2.0000,  2.0000]
cutsets["Analysis1p5Em2"]["ptjj"]=[ 145,  145.0000]
cutsets["Analysis1p5Em2"]["ptlep"]=[ 25,  25.0000]
cutsets["Analysis1p5Em2"]["drjj"]=[ 1.49,  1.49]
cutsets["Analysis1p5Em2"]["nsa5jet"]=[ 3.0000,  3.0000]

cutsets["Analysis2Em2"]={}
cutsets["Analysis2Em2"]["etalep"]=[ 2.5,  2.5]
cutsets["Analysis2Em2"]["jet2pt"]=[ 25,  25.0000]
cutsets["Analysis2Em2"]["lepmetdphi"]=[ 1.85,  1.85]
cutsets["Analysis2Em2"]["jet2csv"]=[ 0.66,  0.66]
cutsets["Analysis2Em2"]["jet1pt"]=[ 25,  25.0000]
cutsets["Analysis2Em2"]["mtW"]=[ 200,  200]
cutsets["Analysis2Em2"]["HVdPhi"]=[ 2.99,  2.99]
cutsets["Analysis2Em2"]["met"]=[ 0,  0]
cutsets["Analysis2Em2"]["topmass"]=[ 210,  210]
cutsets["Analysis2Em2"]["relisolep"]=[ 0.15,  0.12]
cutsets["Analysis2Em2"]["ptW"]=[ 152,  152]
cutsets["Analysis2Em2"]["jet1csv"]=[ 0.935,  0.935]
cutsets["Analysis2Em2"]["ptaddjet"]=[ 20,  20.0000]
cutsets["Analysis2Em2"]["naddjet"]=[ 1.0000,  1.0000]
cutsets["Analysis2Em2"]["ptjj"]=[ 146,  146]
cutsets["Analysis2Em2"]["ptlep"]=[ 25,  25.0000]
cutsets["Analysis2Em2"]["drjj"]=[ 1.49,  1.49]
cutsets["Analysis2Em2"]["nsa5jet"]=[ 3.0000,  3.0000]

cutsets["Analysis3Em2"]={}
cutsets["Analysis3Em2"]["etalep"]=[ 2.5,  2.5]
cutsets["Analysis3Em2"]["jet2pt"]=[ 25,  25]
cutsets["Analysis3Em2"]["lepmetdphi"]=[ 1.85,  1.85]
cutsets["Analysis3Em2"]["jet2csv"]=[ 0.66,  0.66]
cutsets["Analysis3Em2"]["jet1pt"]=[ 25,  25.0000]
cutsets["Analysis3Em2"]["mtW"]=[ 200,  200]
cutsets["Analysis3Em2"]["HVdPhi"]=[ 3.01,  3.01]
cutsets["Analysis3Em2"]["met"]=[ 0,  0]
cutsets["Analysis3Em2"]["topmass"]=[ 210,  210]
cutsets["Analysis3Em2"]["relisolep"]=[ 0.15,  0.12]
cutsets["Analysis3Em2"]["ptW"]=[ 155,  155]
cutsets["Analysis3Em2"]["jet1csv"]=[ 0.935,  0.935]
cutsets["Analysis3Em2"]["ptaddjet"]=[ 20.0000,  20.0000]
cutsets["Analysis3Em2"]["naddjet"]=[ 1.0000,  1.0000]
cutsets["Analysis3Em2"]["ptjj"]=[ 150,  150]
cutsets["Analysis3Em2"]["ptlep"]=[ 25,  25.0000]
cutsets["Analysis3Em2"]["drjj"]=[ 1.49,  1.49]
cutsets["Analysis3Em2"]["nsa5jet"]=[ 2.0000,  2.0000]

cutsets["Analysis4Em2"]={}
cutsets["Analysis4Em2"]["etalep"]=[ 2.5,  2.5]
cutsets["Analysis4Em2"]["jet2pt"]=[ 25,  25]
cutsets["Analysis4Em2"]["lepmetdphi"]=[ 1.85,  1.85]
cutsets["Analysis4Em2"]["jet2csv"]=[ 0.66,  0.66]
cutsets["Analysis4Em2"]["jet1pt"]=[ 25,  25.0000]
cutsets["Analysis4Em2"]["mtW"]=[ 200,  200]
cutsets["Analysis4Em2"]["HVdPhi"]=[ 3.01,  3.01]
cutsets["Analysis4Em2"]["met"]=[ 0,  0]
cutsets["Analysis4Em2"]["topmass"]=[ 210,  210.0000]
cutsets["Analysis4Em2"]["relisolep"]=[ 0.15,  0.12]
cutsets["Analysis4Em2"]["ptW"]=[ 155,  155.0000]
cutsets["Analysis4Em2"]["jet1csv"]=[ 0.935,  0.935]
cutsets["Analysis4Em2"]["ptaddjet"]=[ 20,  20.0000]
cutsets["Analysis4Em2"]["naddjet"]=[ 1.0000,  1.0000]
cutsets["Analysis4Em2"]["ptjj"]=[ 150.0000,  150.0000]
cutsets["Analysis4Em2"]["ptlep"]=[ 25,  25.0000]
cutsets["Analysis4Em2"]["drjj"]=[ 1.33,  1.33]
cutsets["Analysis4Em2"]["nsa5jet"]=[ 2.0000,  2.0000]

cutsets["Analysis5Em2"]={}
cutsets["Analysis5Em2"]["etalep"]=[ 2.5,  2.5]
cutsets["Analysis5Em2"]["jet2pt"]=[ 25,  2.5]
cutsets["Analysis5Em2"]["lepmetdphi"]=[ 1.85,  1.85]
cutsets["Analysis5Em2"]["jet2csv"]=[ 0.66,  0.66]
cutsets["Analysis5Em2"]["jet1pt"]=[ 25,  25]
cutsets["Analysis5Em2"]["mtW"]=[ 200,  200]
cutsets["Analysis5Em2"]["HVdPhi"]=[ 3.02,  3.02]
cutsets["Analysis5Em2"]["met"]=[ 0,  0]
cutsets["Analysis5Em2"]["topmass"]=[ 230,  230]
cutsets["Analysis5Em2"]["relisolep"]=[ 0.15,  0.12]
cutsets["Analysis5Em2"]["ptW"]=[ 162,  162]
cutsets["Analysis5Em2"]["jet1csv"]=[ 0.35,  0.935]
cutsets["Analysis5Em2"]["ptaddjet"]=[ 20,  20.0000]
cutsets["Analysis5Em2"]["naddjet"]=[ 1.0000,  1.0000]
cutsets["Analysis5Em2"]["ptjj"]=[ 150,  150]
cutsets["Analysis5Em2"]["ptlep"]=[ 25,  25]
cutsets["Analysis5Em2"]["drjj"]=[ 1.33,  1.33]
cutsets["Analysis5Em2"]["nsa5jet"]=[ 2.0000,  2.0000]

nBins=40
varstocut={}  # 0 keep less than, 1 keept greater than
varstocut["jet1csv"]     = [1,nBins,cutsets["presel"]["jet1csv"][0], 1.0, nBins, cutsets["presel"]["jet1csv"][1], 1.0 ]
varstocut["jet2csv"]     = [1,nBins,cutsets["presel"]["jet2csv"][0], 1.0, nBins, cutsets["presel"]["jet2csv"][1], 1.0 ]
varstocut["ptjj"]        = [1,nBins,cutsets["presel"]["ptjj"][0],300,nBins,cutsets["presel"]["ptjj"][1], 300]
varstocut["ptW"]         = [1,nBins,cutsets["presel"]["ptW"][0],300,nBins,cutsets["presel"]["ptW"][1], 300]
varstocut["met"]         = [1,nBins,cutsets["presel"]["met"][0],100,nBins,cutsets["presel"]["met"][1], 100]
varstocut["HVdPhi"]      = [1,nBins,cutsets["presel"]["HVdPhi"][0], 3.14, nBins,cutsets["presel"]["HVdPhi"][1],3.15]
varstocut["lepmetdphi"]  = [0,nBins,0.,cutsets["presel"]["lepmetdphi"][0],nBins,0.,cutsets["presel"]["lepmetdphi"][1]]
#varstocut["naddlep"]     = [0,5,0,5,5,0,5]
varstocut["naddjet"]     = [0,5,0,5,5,0,5]
varstocut["topmass"] = [1,nBins,cutsets["presel"]["topmass"][0],500,nBins,cutsets["presel"]["topmass"][1],500]
varstocut["ptlep"]       = [1,nBins,cutsets["presel"]["ptlep"][0], 100, nBins, cutsets["presel"]["ptlep"][1], 100 ]
varstocut["jet1pt"]       = [1,15,cutsets["presel"]["jet1pt"][0], 100, 15, cutsets["presel"]["jet1pt"][1], 100 ]
varstocut["jet2pt"]       = [1,15,cutsets["presel"]["jet2pt"][0], 100, 15, cutsets["presel"]["jet2pt"][1], 100 ]
varstocut["etalep"]     = [0,nBins,0.,2.5,nBins,0.,2.5]
varstocut["relisolep"]  = [0,nBins,0.,cutsets["presel"]["relisolep"][0],nBins,0.,cutsets["presel"]["relisolep"][1] ]
varstocut["drjj"]     =  [0,nBins,0.,cutsets["presel"]["drjj"][0],nBins,0.,cutsets["presel"]["drjj"][1] ]
#varstocut["hwptbal"]    = [0,nBins,0.,cutsets["presel"]["hwptbal"][0],nBins,0.,cutsets["presel"]["hwptbal"][1] ]
varstocut["ptaddjet"]   = [1,nBins,cutsets["presel"]["ptaddjet"][0],300,nBins,cutsets["presel"]["ptaddjet"][1],300 ]
varstocut["mtW"]        = [0,nBins,0.,cutsets["presel"]["mtW"][0],nBins,0.,cutsets["presel"]["mtW"][1] ]
varstocut["nsa5jet"]    = [0,10,0,10,10,0,10]

absCutVars=["HVdPhi","lepmetdphi","drjj","etalep"]

passFlags = {}
# assign a branch which keeps track of whether an event passes each given working point
for wp in cutsets:
    passFlags[wp] = np.zeros(1, dtype=float)
    otree.Branch("Pass_%s" % wp, passFlags[wp], "Pass_%s/D" % wp)

nentries = tree.GetEntries()
#nentries = 1000
print "Total entries: %i" % nentries
for ientry in range(nentries):
    if (ientry % 10000 == 0): print "processing entry: %i" % ientry
    tree.GetEntry(ientry)
    AssignVarsToMap()
    for wp in cutsets:
        passWP = True
        for var in cutsets[wp]:
            if var in absCutVars: varValue=abs(varsmap[var])
            else: varValue = varsmap[var]
            if (varstocut[var][0] == 1):
                # cut is greater than
                if varValue < cutsets[wp][var][0]: 
                    passWP=False
                    #print "failed %s" % var
            else:
                # cut is less than
                if varValue > cutsets[wp][var][0]:
                    passWP=False
                    #print "failed %s" % var
        if (passWP): passFlags[wp][0] = 1.0;
        else: passFlags[wp][0] = 0.;
    otree.Fill()

ofile.cd()
otree.Write()
ofile.Close()        

