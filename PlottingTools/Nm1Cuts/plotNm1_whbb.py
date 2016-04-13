import ROOT
import sys,os
import numpy,array
import math

ROOT.gROOT.SetBatch(True)

filename=sys.argv[1]
tfile=ROOT.TFile(filename)

treename="tree"
tree=tfile.Get(treename)

outDir="variablePlots_WHbb_April11_5Em2"
suffix="analysis"

if not os.path.exists(outDir):
    os.makedirs(outDir)

flatTree=True

varsmap={}
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


yieldsTextFile=open(outDir+"/yields.txt","w+")


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

cutsets["Analysis"] = {}
cutsets["Analysis"]["jet1csv"]       =   [0.935, 0.935]
cutsets["Analysis"]["jet2csv"]       =   [0.46, 0.46]
cutsets["Analysis"]["jet1pt"]        =   [25, 25]
cutsets["Analysis"]["jet2pt"]        =   [25, 25]
cutsets["Analysis"]["ptjj"]          =   [100, 100]
cutsets["Analysis"]["ptW"]           =   [100, 100]
cutsets["Analysis"]["met"]           =   [0, 0]
cutsets["Analysis"]["HVdPhi"]        =   [2.5, 2.5]
cutsets["Analysis"]["lepmetdphi"]    =   [2.0, 2.0]
#cutsets["Analysis"]["naddlep"]       =   [1, 1]
cutsets["Analysis"]["naddjet"]       =   [2, 2]
cutsets["Analysis"]["topmass"] =   [100, 100]
cutsets["Analysis"]["ptlep"]         =   [25, 25]
cutsets["Analysis"]["etalep"]       =   [2.5,2.5]
cutsets["Analysis"]["relisolep"]    =   [0.15, 0.12]
cutsets["Analysis"]["drjj"]         =   [3.5, 3.5]
#cutsets["Analysis"]["hwptbal"]      =   [3,3]
cutsets["Analysis"]["ptaddjet"]      =   [20,20]
cutsets["Analysis"]["mtW"]          = [200,200]
cutsets["Analysis"]["nsa5jet"]      = [10, 10]

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
cutsets["Analysis5Em2"]["jet1csv"]=[ 0.935,  0.935]
cutsets["Analysis5Em2"]["ptaddjet"]=[ 20,  20.0000]
cutsets["Analysis5Em2"]["naddjet"]=[ 1.0000,  1.0000]
cutsets["Analysis5Em2"]["ptjj"]=[ 150,  150]
cutsets["Analysis5Em2"]["ptlep"]=[ 25,  25]
cutsets["Analysis5Em2"]["drjj"]=[ 1.33,  1.33]
cutsets["Analysis5Em2"]["nsa5jet"]=[ 2.0000,  2.0000]

## [min, max] possible cut values
#boundaries["threshold"] = {}
#boundaries["threshold"]["jet1csv"]       =   [0.9, 1.0]
#boundaries["threshold"]["jet2csv"]       =   [0.65, 1.0]
#boundaries["threshold"]["jet1pt"]        =   [25, 25]
#boundaries["threshold"]["jet2pt"]        =   [25, 25]
#boundaries["threshold"]["ptjj"]          =   [55, 1000]
#boundaries["threshold"]["ptW"]           =   [60, 1000]
#boundaries["threshold"]["met"]           =   [0, 200]
#boundaries["threshold"]["HVdPhi"]        =   [2.8,3.15]
#boundaries["threshold"]["lepmetdphi"]    =   [0.0, 2.0]
#boundaries["threshold"]["naddlep"]       =   [1, 2]
#boundaries["threshold"]["naddjet"]       =   [1, 2]
#boundaries["threshold"]["topmass"] =   [0, 1000]
#boundaries["threshold"]["ptlep"]         =   [25, 25]

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

nCats=2

def deltaPhi(phi1,phi2):
    dphi=math.fabs(phi1-phi2)
    
    if dphi>math.pi:
        dphi=2*math.pi - dphi

    return dphi

def deltaR(eta1,phi1,eta2,phi2):
    dphi=deltaPhi(phi1,phi2)
    deta=math.fabs(eta1-eta2)
    dR=math.sqrt(dphi*dphi + deta*deta)
    return dR




absCutVars=["HVdPhi","lepmetdphi","drjj","etalep"]
def passVarCut(var,index,wp):
    if flatTree:
        varValue=varsmap[var]
    else:
        varValue=varsmap[var][index]

    if var in absCutVars:
        varValue=math.fabs(varValue)
 
    catInd=category(nCats,index)
    if varstocut[var][0]==0:
        if varValue >= cutsets[wp][var][catInd]:
            return False
    
    if varstocut[var][0]==1:
        if varValue <= cutsets[wp][var][catInd]:
            return False
    
    return True


def passNminus1(vartoplot,index,wp):
    for var in varstocut:
        if vartoplot != var:
            if passVarCut(var,index,wp) is False:
                return False

    return True


def passAnyNminus1(index,wp):
    cutsFailed=[]

    for var in varstocut:
        if passVarCut(var,index,wp) is False:
            cutsFailed.append(var)
        if len(cutsFailed)>1:
            break

    return cutsFailed


def category(nCats,index):
    isWmunu = varsmap["isWmunu"]
    isWenu  = varsmap["isWenu"]
    
    if nCats is 2:
        if (isWmunu): return 0
        elif (isWenu): return 1

    print "nCats is ",nCats,"what is happening?  returning -1"
    return -1

catMap={}
catMap[0]="isWmunu"
catMap[1]="isWenu"




hists={}
kinplots={}

tree.SetBranchStatus("*", 1)
      
nev = tree.GetEntries()
print "loop over ",nev

#mctypes =["WenuH","W","DY","TT","QCD"]
#mccolors=[632,633,634,417,402]

#mctypes =["DY","TTReal","TTFake","GJet"]
mctypes = ["VH","Bkg"]
mccolors=[880,416]
#mctypes =["DY","TTFake"]
#mccolors=[633,417]
bkgType=mctypes[1]
sigType=mctypes[0]


requireHLT=False
etaRestriction=False
matchL1=False
doSoB=True
doCutScan=False # only works if doSoB is True
makeTree=False

#eleHoE_PassNm1        = numpy.zeros(1, dtype=int)
#b_eleHoE_PassNm1=tree.Branch( 'eleHoE_PassNm1', eleHoE_PassNm1, 'eleHoE_PassNm1/I' )

#wps=["preselglobal","Veto","Loose","Medium"]
wps=["Analysis4Em2"]

boundaryWP = "Analysis4Em2" # don't let any of the cuts get looser than this wp
potentialCuts={}  #wp,cat,var

numIterations=5

# set boundaries so that cuts can never be looser than initial working point
boundaries = {}
for var in varstocut:
    boundaries[var] = {}
    for cat in catMap: 
       boundaries[var][cat] = cutsets[boundaryWP][var][cat] 

for wp in wps:
    print wp
    for iteration in range(numIterations):
        if iteration==0:
            originalwp=wp
        else:
            wp=wpPlus1
        wpPlus1=originalwp+"iteration"+str(iteration)
        print wp,wpPlus1
        potentialCuts[wpPlus1]={}
        #if makeTree:
        #    newfile = ROOT.TFile("merged_"+wp+"_"+filename.split("/")[-1],"recreate")
        #    newtree = tree.CloneTree(0)
        
        for var in varstocut:
            potentialCuts[wpPlus1][var]={}
            for cat in catMap:
                potentialCuts[wpPlus1][var][catMap[cat]]=[]
        
        for cat in catMap:
            for mctype in mctypes:
                kinplots[mctype+"_type_"+catMap[cat]+"_"+wp]=ROOT.TH1F(mctype+"_type_"+catMap[cat]+"_"+wp,";"+mctype+" MC type "+catMap[cat]+"_"+wp,50,0,50)
                kinplots[mctype+"_Mjj_"+catMap[cat]+"_"+wp]=ROOT.TH1F(mctype+"_Mjj_"+catMap[cat]+"_"+wp,";"+mctype+" di-jet mass "+catMap[cat]+"_"+wp,50,0,500)
        
        for var in varstocut:
            print var
            for cat in catMap:
                for mctype in mctypes:
                    #print "constructing th1f ",mctype+"_"+catMap[cat]+"_"+var+"_nm1_"+wp
                    hists[mctype+"_"+catMap[cat]+"_"+var+"_nm1_"+wp]=ROOT.TH1F(mctype+"_"+catMap[cat]+"_"+var+"_nm1_"+wp,mctype+"_"+catMap[cat]+"_"+var+"_nm1_"+wp,varstocut[var][1+3*cat],varstocut[var][2+3*cat],varstocut[var][3+3*cat])
                    #print hists[mctype+"_"+catMap[cat]+"_"+var+"_nm1_"+wp].Integral()
         
        for iev in range(0,nev):
        #for iev in range(0,300):
            if iev%100000==0:
                print "event",iev
            tree.GetEntry(iev)
            if flatTree:
                AssignVarsToMap()
            #if tree.nhltele == 0:
            #    continue

            hltindex=0
            hltptpassing=0
            

            if varsmap["type"]==-12501:
                mctype="VH"
            elif varsmap["type"]>0 and varsmap["type"]!=12:
                mctype="Bkg"
            else:
                continue
                print "sampleIndex",varsmap["type"]
                print "I don't want to live in a world of unknowns!!"
            # add by hand a mass window
            if (varsmap["Mjj"] < 90 or varsmap["Mjj"] > 150):
                continue
            if (tree.Vtype!=2 and tree.Vtype!=3):
                continue

            # preselection            
            #if varsmap["pt"]<20 or tree.passConversionVeto!=1 or math.fabs(tree.dz)>1.0 or (math.fabs(varsmap["eta"])<1.566 and math.fabs(varsmap["eta"])>1.4442):
            #    continue
                
            cutsFailed=passAnyNminus1(-1,wp)
           
            ## changed to 2 so I can see which cuts fail both iso cuts 
            #if len(cutsFailed)>2:
            #    if makeTree:
            #        b_eleHoE_PassNm1.Fill()
            #    continue
           
            #if hltindex >= tree.nhltele:
            #    print hltindex,tree.nhltele
            #    continue
 
            #if makeTree:
            #    if len(cutsFailed)==0:
            #        tree.eleHoE_PassNm1[hltindex]=1
            #            
            #    elif len(cutsFailed)==1:
            #        if cutsFailed[0]=="h_e":
            #            tree.eleHoE_PassNm1[hltindex]=1
            #        #elif cutsFailed[0]=="hcaliso_pt":
            #        #    tree.hltEleHCalIsoOPt_PassNm1[hltindex]=1
            #        #    tree.hltEleBothIsoOPt_PassNm2[hltindex]=1
            #        #elif cutsFailed[0]=="ecaliso_pt":
            #        #    tree.hltEleECalIsoOPt_PassNm1[iel]=1
            #        #    tree.hltEleBothIsoOPt_PassNm2[hltindex]=1
            #    
            #    elif len(cutsFailed)==2:
            #        #if (cutsFailed[0]=="hcaliso_pt" and cutsFailed[1]=="ecaliso_pt") or (cutsFailed[1]=="hcaliso_pt" and cutsFailed[0]=="ecaliso_pt"):
            #        #    tree.hltEleBothIsoOPt_PassNm2[hltindex]=1
            #        #else:
            #        continue         
            #    b_eleHoE_PassNm1.Fill()
                
            if len(cutsFailed)==0:
                if flatTree:
                    for var in varstocut:
                        if var in absCutVars:
                            hists[mctype+"_"+catMap[category(nCats,hltindex)]+"_"+var+"_nm1_"+wp].Fill(math.fabs(varsmap[var]),tree.weight*tree.weight_PU*tree.bTagWeight*tree.selLeptons_SF_IsoTight[tree.lepInd]*(tree.selLeptons_SF_IdCutTight[tree.lepInd]*tree.isWmunu + tree.selLeptons_SF_IdMVATight[tree.lepInd]*tree.isWenu)*tree.CS_SF)
                        else:
                            hists[mctype+"_"+catMap[category(nCats,hltindex)]+"_"+var+"_nm1_"+wp].Fill(varsmap[var],tree.weight*tree.weight_PU*tree.bTagWeight*tree.selLeptons_SF_IsoTight[tree.lepInd]*(tree.selLeptons_SF_IdCutTight[tree.lepInd]*tree.isWmunu + tree.selLeptons_SF_IdMVATight[tree.lepInd]*tree.isWenu)*tree.CS_SF)
                    kinplots[mctype+"_type_"+catMap[category(nCats,hltindex)]+"_"+wp].Fill(varsmap['type'],tree.weight*tree.weight_PU*tree.bTagWeight*tree.selLeptons_SF_IsoTight[tree.lepInd]*(tree.selLeptons_SF_IdCutTight[tree.lepInd]*tree.isWmunu + tree.selLeptons_SF_IdMVATight[tree.lepInd]*tree.isWenu)*tree.CS_SF)
                    kinplots[mctype+"_Mjj_"+catMap[category(nCats,hltindex)]+"_"+wp].Fill(varsmap['Mjj'],tree.weight*tree.weight_PU*tree.bTagWeight*tree.selLeptons_SF_IsoTight[tree.lepInd]*(tree.selLeptons_SF_IdCutTight[tree.lepInd]*tree.isWmunu + tree.selLeptons_SF_IdMVATight[tree.lepInd]*tree.isWenu)*tree.CS_SF)
                #else:
                #    for var in varstocut:
                #        hists[mctype+"_"+catMap[category(nCats,hltindex)]+"_"+var+"_nm1_"+wp].Fill(varsmap[var][hltindex],tree.weight*tree.weight_PU*tree.bTagWeight*tree.selLeptons_SF_IsoTight[tree.lepInd]*(tree.selLeptons_SF_IdCutTight[tree.lepInd]*tree.isWmunu + tree.selLeptons_SF_IdMVATight[tree.lepInd]*tree.isWenu))
                #    kinplots[mctype+"_eta_"+catMap[category(nCats,hltindex)]+"_"+wp].Fill(tree.hltEleEta[hltindex],tree.weight*tree.weight_PU*tree.bTagWeight*tree.selLeptons_SF_IsoTight[tree.lepInd]*(tree.selLeptons_SF_IdCutTight[tree.lepInd]*tree.isWmunu + tree.selLeptons_SF_IdMVATight[tree.lepInd]*tree.isWenu))
                #    kinplots[mctype+"_pt_"+catMap[category(nCats,hltindex)]+"_"+wp].Fill(tree.hltElePt[hltindex],tree.weight*tree.weight_PU*tree.bTagWeight*tree.selLeptons_SF_IsoTight[tree.lepInd]*(tree.selLeptons_SF_IdCutTight[tree.lepInd]*tree.isWmunu + tree.selLeptons_SF_IdMVATight[tree.lepInd]*tree.isWenu))
                #    kinplots[mctype+"_rho_"+catMap[category(nCats,hltindex)]+"_"+wp].Fill(tree.hltElePt[hltindex],tree.weight*tree.weight_PU*tree.bTagWeight*tree.selLeptons_SF_IsoTight[tree.lepInd]*(tree.selLeptons_SF_IdCutTight[tree.lepInd]*tree.isWmunu + tree.selLeptons_SF_IdMVATight[tree.lepInd]*tree.isWenu))
        
            if len(cutsFailed)==1:
                if flatTree:
                    if (cutsFailed[0] in absCutVars):
                        hists[mctype+"_"+catMap[category(nCats,hltindex)]+"_"+cutsFailed[0]+"_nm1_"+wp].Fill(math.fabs(varsmap[cutsFailed[0]]),tree.weight*tree.weight_PU*tree.bTagWeight*tree.selLeptons_SF_IsoTight[tree.lepInd]*(tree.selLeptons_SF_IdCutTight[tree.lepInd]*tree.isWmunu + tree.selLeptons_SF_IdMVATight[tree.lepInd]*tree.isWenu)*tree.CS_SF)
                    else:
                        hists[mctype+"_"+catMap[category(nCats,hltindex)]+"_"+cutsFailed[0]+"_nm1_"+wp].Fill(varsmap[cutsFailed[0]],tree.weight*tree.weight_PU*tree.bTagWeight*tree.selLeptons_SF_IsoTight[tree.lepInd]*(tree.selLeptons_SF_IdCutTight[tree.lepInd]*tree.isWmunu + tree.selLeptons_SF_IdMVATight[tree.lepInd]*tree.isWenu)*tree.CS_SF)
                else:
                    if (cutsFailed[0] in absCutVars):
                        hists[mctype+"_"+catMap[category(nCats,hltindex)]+"_"+cutsFailed[0]+"_nm1_"+wp].Fill(math.fabs(varsmap[cutsFailed[0]][hltindex]),tree.weight*tree.weight_PU*tree.bTagWeight*tree.selLeptons_SF_IsoTight[tree.lepInd]*(tree.selLeptons_SF_IdCutTight[tree.lepInd]*tree.isWmunu + tree.selLeptons_SF_IdMVATight[tree.lepInd]*tree.isWenu)*tree.CS_SF)
                    else:
                        hists[mctype+"_"+catMap[category(nCats,hltindex)]+"_"+cutsFailed[0]+"_nm1_"+wp].Fill(varsmap[cutsFailed[0]][hltindex],tree.weight*tree.weight_PU*tree.bTagWeight*tree.selLeptons_SF_IsoTight[tree.lepInd]*(tree.selLeptons_SF_IdCutTight[tree.lepInd]*tree.isWmunu + tree.selLeptons_SF_IdMVATight[tree.lepInd]*tree.isWenu)*tree.CS_SF)
                    
            #if makeTree:
            #    npass=0
            #    for iel in range(0,tree.nhltele):
            #        if newtree.hltEle_Pass[iel]==1:
            #            npass=npass+1
            #            if npass>1:
            #                print "WTF mate 3"
            #    newtree.Fill()
                    
                    
        
        
        for var in varstocut:
            for icat in catMap:
                can=ROOT.TCanvas(catMap[icat]+"_"+var+"_nm1_"+wp,"",700,500)
                if doSoB:
                    hists["sob_"+catMap[icat]+"_"+var+"_nm1_"+wp]=hists[sigType+"_"+catMap[icat]+"_"+var+"_nm1_"+wp].Clone()
                    bkgtemp=hists[bkgType+"_"+catMap[icat]+"_"+var+"_nm1_"+wp].Clone()
                    bkgtemp.Smooth(10)
                    if doCutScan:
                        hist_SoSqBScan = hists["sob_"+catMap[icat]+"_"+var+"_nm1_"+wp].Clone()
                        hist_SoSqBScan.Reset() # keep binning but clear contents
                        numberBins = hists["sob_"+catMap[icat]+"_"+var+"_nm1_"+wp].GetNbinsX()
                        if (bkgtemp.GetNbinsX() != numberBins):
                            print "Number of bins in the signal histogram (%i) is different than the number of bins in the background histogram (%i) for variable %s" % (numberBins,bkgtemp.GetNbinsX(),var)
                            sys.exit(1)
                        for ibin in range(numberBins):
                            S_tot = 0.
                            B_tot = 0.
                            if (varstocut[var][0] == 0):
                                # cut is less than 
                                S_tot = hists["sob_"+catMap[icat]+"_"+var+"_nm1_"+wp].Integral(0,ibin+1) 
                                B_tot = bkgtemp.Integral(1,ibin+1)
                            elif (varstocut[var][0] == 1):
                                # cut is greater than
                                S_tot = hists["sob_"+catMap[icat]+"_"+var+"_nm1_"+wp].Integral(ibin+1,numberBins+1) 
                                B_tot = bkgtemp.Integral(ibin+1,numberBins)
                            if (B_tot != 0):
                                hist_SoSqBScan.SetBinContent(ibin+1, (S_tot / math.sqrt(B_tot)))
                                #print S_tot, B_tot, (S_tot / math.sqrt(B_tot))
                            else:
                                hist_SoSqBScan.SetBinContent(ibin+1, 0.)
                        hists["sob_"+catMap[icat]+"_"+var+"_nm1_"+wp] = hist_SoSqBScan.Clone()
                        print "Integral = ", hists["sob_"+catMap[icat]+"_"+var+"_nm1_"+wp].Integral()
                    else:
                        hists["sob_"+catMap[icat]+"_"+var+"_nm1_"+wp].Divide(bkgtemp)
                        numberBins = hists["sob_"+catMap[icat]+"_"+var+"_nm1_"+wp].GetNbinsX()
                        hists["sob_"+catMap[icat]+"_"+var+"_nm1_"+wp].SetBinContent(numberBins, hists["sob_"+catMap[icat]+"_"+var+"_nm1_"+wp].GetBinContent(numberBins) + hists["sob_"+catMap[icat]+"_"+var+"_nm1_"+wp].GetBinContent(numberBins+1) )
                        hists["sob_"+catMap[icat]+"_"+var+"_nm1_"+wp].SetBinContent(1, hists["sob_"+catMap[icat]+"_"+var+"_nm1_"+wp].GetBinContent(0) + hists["sob_"+catMap[icat]+"_"+var+"_nm1_"+wp].GetBinContent(1) ) 
               
                itype=0 
                maxval=0
                for mctype in mctypes:
                    if hists[mctype+"_"+catMap[icat]+"_"+var+"_nm1_"+wp].Integral() > 0:
                        hists[mctype+"_"+catMap[icat]+"_"+var+"_nm1_"+wp].Scale(1/ hists[mctype+"_"+catMap[icat]+"_"+var+"_nm1_"+wp].Integral())
                    hists[mctype+"_"+catMap[icat]+"_"+var+"_nm1_"+wp].SetLineWidth(2)
                    hists[mctype+"_"+catMap[icat]+"_"+var+"_nm1_"+wp].SetLineColor(mccolors[itype])
                    itype=itype+1
                    maxval=max(maxval,hists[mctype+"_"+catMap[icat]+"_"+var+"_nm1_"+wp].GetMaximum())
        
                maxval=maxval*1.2
        
                first=True 
                leg=ROOT.TLegend(0.7,0.6,0.9,0.9)
                leg.SetHeader(catMap[icat])
                for mctype in mctypes:
                    if first is True:    
                        hists[mctype+"_"+catMap[icat]+"_"+var+"_nm1_"+wp].GetMaximum(maxval)
                        hists[mctype+"_"+catMap[icat]+"_"+var+"_nm1_"+wp].SetTitle(";"+var+" "+catMap[icat]+" "+wp)
                        hists[mctype+"_"+catMap[icat]+"_"+var+"_nm1_"+wp].Draw()
                        first=False
                    else:
                        hists[mctype+"_"+catMap[icat]+"_"+var+"_nm1_"+wp].Draw("sames")
                    leg.AddEntry(hists[mctype+"_"+catMap[icat]+"_"+var+"_nm1_"+wp],mctype,"l")
   
                leg.Draw("same") 
                line=ROOT.TLine(cutsets[wp][var][icat],0,cutsets[wp][var][icat], maxval/1.3)
                line.SetLineColor(4)    
                line.SetLineWidth(2)
                line.Draw("same")
         
                can.Update() 
                can.SaveAs(outDir+"/s_b_"+catMap[icat]+"_"+var+"_"+wp+"_"+suffix+".png")
                can.SaveAs(outDir+"/s_b_"+catMap[icat]+"_"+var+"_"+wp+"_"+suffix+".C")
                can.Write()
                can.Clear()
        
                
                if doSoB:
                    hists["sob_"+catMap[icat]+"_"+var+"_nm1_"+wp].SetLineWidth(2)
                    hists["sob_"+catMap[icat]+"_"+var+"_nm1_"+wp].Draw()
                    #hists["sob_"+catMap[icat]+"_"+var+"_nm1_"+wp].SetMinimum(0.1)
                    temp=hists["sob_"+catMap[icat]+"_"+var+"_nm1_"+wp].Clone()
                    line=ROOT.TLine(cutsets[wp][var][icat],0,cutsets[wp][var][icat], hists["sob_"+catMap[icat]+"_"+var+"_nm1_"+wp].GetMaximum()/1.3)
                    line.SetLineColor(4)    
                    line.SetLineWidth(2)
                    line.Draw("same")
        
                    hLineVal=[0.05,0.05]
                    #[300,15]
                    #[7.5,1.7]
                    #[900,400] 
                    #[250,145]
                    #[50,35]
                    #[7.5,7.5]
                    line2=ROOT.TLine(temp.GetBinLowEdge(1),hLineVal[icat],temp.GetBinLowEdge(temp.GetNbinsX()+1),hLineVal[icat])
                    line2.SetLineColor(633)    
                    line2.SetLineWidth(2)
                    line2.Draw("same")
                    
                    can.SetLogy(1)
                    
                    can.Update() 
                    can.SaveAs(outDir+"/soverb_"+catMap[icat]+"_"+var+"_"+wp+"_"+suffix+".png")
                    can.SaveAs(outDir+"/soverb_"+catMap[icat]+"_"+var+"_"+wp+"_"+suffix+".C")
                    #if (var=="ptlep" or var=="jet1pt" or var=="jet2pt"):
                    #    # manually keep the jet and lepton pt cuts fixed for now
                    #    potentialCuts[wpPlus1][var][catMap[icat]].append(25.)
                    #elif (var.find("nadd")>-1):
                    #    # add protection for high S/B targets so naddjet/lep cuts stay fixed
                    #    potentialCuts[wpPlus1][var][catMap[icat]].append(1)
                    #elif (varstocut[var][0] == 0):
                    if (varstocut[var][0] == 0):
                        # cut is less than
                        for iBin in range(hists["sob_"+catMap[icat]+"_"+var+"_nm1_"+wp].GetNbinsX(), 0, -1):
                            lowEdge=hists["sob_"+catMap[icat]+"_"+var+"_nm1_"+wp].GetBinLowEdge(iBin) + hists["sob_"+catMap[icat]+"_"+var+"_nm1_"+wp].GetBinWidth(iBin)
                            binVal=hists["sob_"+catMap[icat]+"_"+var+"_nm1_"+wp].GetBinContent(iBin)
                            #if binVal > hLineVal[icat] and lowEdge >= boundaries["threshold"][var][0] and lowEdge <= boundaries["threshold"][var][1]:
                            if binVal > hLineVal[icat] and lowEdge <= boundaries[var][icat]:
                                #potentialCuts[wpPlus1][var][catMap[icat]].append(lowEdge)
                                # instead only move half the distance so it's more stable
                                if (var.find("nadd")==-1 and var.find("nsa")==-1):
                                    potentialCuts[wpPlus1][var][catMap[icat]].append( (lowEdge + cutsets[wp][var][icat])/2.)
                                else:
                                    potentialCuts[wpPlus1][var][catMap[icat]].append(int((lowEdge + cutsets[wp][var][icat])/2. + 0.99))
                            elif (iBin==1 and len(potentialCuts[wpPlus1][var][catMap[icat]])==0):
                                # no bins have sufficient S/B for target, leave the cut at the boundary
                                potentialCuts[wpPlus1][var][catMap[icat]].append(boundaries[var][icat])
                    elif (varstocut[var][0] == 1):
                        # cut is greater than
                        for iBin in range(1, hists["sob_"+catMap[icat]+"_"+var+"_nm1_"+wp].GetNbinsX()+1):
                            lowEdge=hists["sob_"+catMap[icat]+"_"+var+"_nm1_"+wp].GetBinLowEdge(iBin)
                            binVal=hists["sob_"+catMap[icat]+"_"+var+"_nm1_"+wp].GetBinContent(iBin)
                            #if binVal > hLineVal[icat] and lowEdge >= boundaries["threshold"][var][0] and lowEdge <= boundaries["threshold"][var][1]:
                            if binVal > hLineVal[icat] and lowEdge >= boundaries[var][icat]:
                                #potentialCuts[wpPlus1][var][catMap[icat]].append(lowEdge)
                                # instead only move half the distance so it's more stable
                                if (var.find("nadd")==-1 and var.find("nsa")==-1):
                                    potentialCuts[wpPlus1][var][catMap[icat]].append( (lowEdge + cutsets[wp][var][icat])/2.)
                                else:
                                    potentialCuts[wpPlus1][var][catMap[icat]].append(int((lowEdge + cutsets[wp][var][icat])/2. + 0.99))
                            elif (iBin==hists["sob_"+catMap[icat]+"_"+var+"_nm1_"+wp].GetNbinsX() and len(potentialCuts[wpPlus1][var][catMap[icat]])==0):
                                # no bins have sufficient S/B for target, leave the cut at the boundary
                                potentialCuts[wpPlus1][var][catMap[icat]].append(boundaries[var][icat])
        can.SetLogy(0)

#        for wp in wps:
        if doSoB:
            for icat in catMap:
                for var in varstocut:
                    print wp,catMap[icat],var,
                    cutbin=int((cutsets[wp][var][icat]-hists[sigType+"_"+catMap[icat]+"_"+var+"_nm1_"+wp].GetBinLowEdge(1))/hists[sigType+"_"+catMap[icat]+"_"+var+"_nm1_"+wp].GetBinWidth(1))
                    finalbin=hists[sigType+"_"+catMap[icat]+"_"+var+"_nm1_"+wp].GetNbinsX()
        
                    print hists[sigType+"_"+catMap[icat]+"_"+var+"_nm1_"+wp].Integral(cutbin,finalbin),
                    print hists[bkgType+"_"+catMap[icat]+"_"+var+"_nm1_"+wp].Integral(cutbin,finalbin)
        
        for mctype in mctypes:
            print mctype,
        print
#for     wp in wps:
        for cat in catMap:
            yieldsTextFile.write(wp+"\t"+catMap[cat]+"\n")
            for mctype in mctypes:
                print kinplots[mctype+"_type_"+catMap[cat]+"_"+wp].Integral(),

            
            for ikin in ["Mjj", "type"]:
                if doSoB:
                    kinplots["soverb_"+ikin+"_"+catMap[cat]+"_"+wp]=kinplots[sigType+"_"+ikin+"_"+catMap[cat]+"_"+wp].Clone()
                    kinplots["soverb_"+ikin+"_"+catMap[cat]+"_"+wp].Divide(kinplots[bkgType+"_"+ikin+"_"+catMap[cat]+"_"+wp])
                    kinplots["soverb_"+ikin+"_"+catMap[cat]+"_"+wp].SetMinimum(0.01)
                
                itype=0 
                maxval=0
                leg=ROOT.TLegend(0.7,0.6,0.9,0.9)
                leg.SetHeader(catMap[cat]+" "+wp)
                for mctype in mctypes:
                    if kinplots[mctype+"_"+ikin+"_"+catMap[cat]+"_"+wp].Integral()>0:
                        if ikin=="Mjj":
                            yieldsTextFile.write(mctype+"\t"+str(kinplots[mctype+"_"+ikin+"_"+catMap[cat]+"_"+wp].Integral())+"\n")

                        kinplots[mctype+"_"+ikin+"_"+catMap[cat]+"_"+wp].Scale(1/kinplots[mctype+"_"+ikin+"_"+catMap[cat]+"_"+wp].Integral())
                    kinplots[mctype+"_"+ikin+"_"+catMap[cat]+"_"+wp].SetLineWidth(2)
                    kinplots[mctype+"_"+ikin+"_"+catMap[cat]+"_"+wp].SetLineColor(mccolors[itype])
                    itype=itype+1
                    maxval=max(maxval,kinplots[mctype+"_"+ikin+"_"+catMap[cat]+"_"+wp].GetMaximum())
                    leg.AddEntry(kinplots[mctype+"_"+ikin+"_"+catMap[cat]+"_"+wp],mctype,"l")
   
                leg.Draw("same") 
                line=ROOT.TLine(cutsets[wp][var][cat],0,cutsets[wp][var][cat], maxval/1.3)
        
                maxval=maxval*1.2
        
                first=True 
                for mctype in mctypes:

                    if first is True:    
                        kinplots[mctype+"_"+ikin+"_"+catMap[cat]+"_"+wp].SetMaximum(maxval)
                        kinplots[mctype+"_"+ikin+"_"+catMap[cat]+"_"+wp].Draw()
                        first=False
                    else:
                        kinplots[mctype+"_"+ikin+"_"+catMap[cat]+"_"+wp].Draw("sames")
        

                can.Update()
                can.SaveAs(outDir+"/s_b_"+ikin+"_"+catMap[cat]+"_"+wp+"_"+suffix+".png")
                can.SaveAs(outDir+"/s_b_"+ikin+"_"+catMap[cat]+"_"+wp+"_"+suffix+".C")
                can.Clear()

                if doSoB:
                    kinplots["soverb_"+ikin+"_"+catMap[cat]+"_"+wp].SetLineWidth(2)
                    kinplots["soverb_"+ikin+"_"+catMap[cat]+"_"+wp].Draw()
                    can.Update()
                    can.SaveAs(outDir+"/soverb_"+ikin+"_"+catMap[cat]+"_"+wp+"_"+suffix+".png")
                    can.SaveAs(outDir+"/soverb_"+ikin+"_"+catMap[cat]+"_"+wp+"_"+suffix+".C")
                    can.Clear()
            yieldsTextFile.write("\n")
        if makeTree:
            newfile.Write()
            newfile.Close()

        
        cutsets[wpPlus1]={}
        for var in potentialCuts[wpPlus1]:
            iCat=0
            if (varstocut[var][0] == 1):
                cutsets[wpPlus1][var]=[0,0]
            else: cutsets[wpPlus1][var]=[9999,9999]
            for cat in potentialCuts[wpPlus1][var]:
                iShow=0
                for i in potentialCuts[wpPlus1][var][cat]:
                    cutsets[wpPlus1][var][iCat]=i
                    iShow=iShow+1
                    if iShow>0:
                        break
                iCat=iCat+1

        for wp in potentialCuts:
            print "cutsets[\""+wp+"\"]={}"
            for var in potentialCuts[wp]:
                iCat=0
                for cat in potentialCuts[wp][var]:
                    iShow=0
                    if iCat==0:
                        print "cutsets[\""+wp+"\"][\""+var+"\"]=[",
                    val=9999
                    for i in potentialCuts[wp][var][cat]:
                        val=i
                        iShow=iShow+1
                        if iShow>0:
                            break
            
                    if iCat==0:
                        print "{0:.4f},".format(val)+" ",
                    elif iCat==1:
                        print "{0:.4f}]".format(val)
                    iCat=iCat+1
            print

yieldsTextFile.close()
wps=potentialCuts.keys()
wps.sort()

for wp in potentialCuts:
    print wp
    for var in potentialCuts[wp]:
        print var,
        for cat in potentialCuts[wp][var]:
            iShow=0
            print cat,   
            for i in potentialCuts[wp][var][cat]:
                print "{0:.4f}".format(i),
                iShow=iShow+1
                if iShow>6:
                    break
            print
    print

for wp in potentialCuts:
    print "cutsets[\""+wp+"\"]={}"
    for var in potentialCuts[wp]:
        iCat=0
        for cat in potentialCuts[wp][var]:
            iShow=0
            if iCat==0:
                print "cutsets[\""+wp+"\"][\""+var+"\"]=[",
            val=9999
            for i in potentialCuts[wp][var][cat]:
                val=i
                iShow=iShow+1
                if iShow>0:
                    break
            
            if iCat==0:
                print "{0:.4f},".format(val)+" ",
            elif iCat==1:
                print "{0:.4f}]".format(val)
            iCat=iCat+1
    print


outfile=ROOT.TFile("cutsets.root","recreate")
cuttree=ROOT.TTree("cuttree","each entry is a set of cuts")
setname = ROOT.std.string()


nCat=numpy.zeros(1, dtype=int)
cutsetvars={}
cuttree.Branch("nCat",nCat,"nCat/I") 
cuttree.Branch('cutsetname',setname)


for var in varstocut:
    cutsetvars[var]=numpy.zeros(10, dtype=float)
    cuttree.Branch(var,cutsetvars[var],var+"[nCat]/D")


for cutset in cutsets:
    setname.replace(0, ROOT.std.string.npos, cutset)
    for var in cutsets[cutset]:
        nCat[0]=int(len(cutsets[cutset][var]))
        #print "nCat",nCat
        count=0
        for cut in cutsets[cutset][var]:
            #print cutset,var,cut
            cutsetvars[var][count]=cut
            count=count+1
    cuttree.Fill()

cuttree.Write()
outfile.Save()
outfile.Close()


