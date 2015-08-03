import ROOT
import sys,os
import numpy,array
import math

filename=sys.argv[1]
tfile=ROOT.TFile(filename)

treename="electronTree"
tree=tfile.Get(treename)


kinPlots=["eta","pt","rho"]

cutSets={}

#cutSets["VetoRami"]={}
#cutSets["VetoRami"]["sieie"]        =   [0.012 , 0.0339] 
#cutSets["VetoRami"]["dEta"]         =   [0.0126, 0.0109] 
#cutSets["VetoRami"]["dPhi"]         =   [0.107 , 0.219 ] 
#cutSets["VetoRami"]["h_e"]          =   [0.186 , 0.0962] 
#cutSets["VetoRami"]["combiso_pt"]   =   [0.161 , 0.193 ] 
#cutSets["VetoRami"]["1_e_1_p"]      =   [0.239 , 0.141 ] 
#cutSets["VetoRami"]["d0"]           =   [0.0631, 0.279 ] 
#cutSets["VetoRami"]["dz"]           =   [0.613 , 0.947 ] 
#cutSets["VetoRami"]["missingHits"]  =   [2.1   , 3.1   ] 
#    
#cutSets["LooseRami"]={}
#cutSets["LooseRami"]["sieie"]        =   [0.0105 ,0.0318] 
#cutSets["LooseRami"]["dEta"]         =   [0.0098 ,0.0095] 
#cutSets["LooseRami"]["dPhi"]         =   [0.0929 ,0.181 ] 
#cutSets["LooseRami"]["h_e"]          =   [0.0765 ,0.0824] 
#cutSets["LooseRami"]["combiso_pt"]   =   [0.118  ,0.118 ] 
#cutSets["LooseRami"]["1_e_1_p"]      =   [0.184  ,0.125 ] 
#cutSets["LooseRami"]["d0"]           =   [0.0227 ,0.242 ] 
#cutSets["LooseRami"]["dz"]           =   [0.379  ,0.921 ] 
#cutSets["LooseRami"]["missingHits"]  =   [2.1    ,1.1   ] 
#    
#cutSets["MediumRami"]={}
#cutSets["MediumRami"]["sieie"]        =   [0.0101 ,0.0287 ] 
#cutSets["MediumRami"]["dEta"]         =   [0.00945,0.00773] 
#cutSets["MediumRami"]["dPhi"]         =   [0.0296 ,0.148  ] 
#cutSets["MediumRami"]["h_e"]          =   [0.0372 ,0.0546 ] 
#cutSets["MediumRami"]["combiso_pt"]   =   [0.0987 ,0.0902 ] 
#cutSets["MediumRami"]["1_e_1_p"]      =   [0.118  ,0.104  ] 
#cutSets["MediumRami"]["d0"]           =   [0.0151 ,0.0535 ] 
#cutSets["MediumRami"]["dz"]           =   [0.238  ,0.572  ] 
#cutSets["MediumRami"]["missingHits"]  =   [2.1    ,1.1    ] 
#    
#cutSets["TightRami"]={}
#cutSets["TightRami"]["sieie"]        =   [0.0101,0.0287 ] 
#cutSets["TightRami"]["dEta"]         =   [0.0095,0.00762] 
#cutSets["TightRami"]["dPhi"]         =   [0.0291,0.0439 ] 
#cutSets["TightRami"]["h_e"]          =   [0.0372,0.0544 ] 
#cutSets["TightRami"]["combiso_pt"]   =   [0.0468,0.0759 ] 
#cutSets["TightRami"]["1_e_1_p"]      =   [0.0174,0.01   ] 
#cutSets["TightRami"]["d0"]           =   [0.0144,0.0377 ] 
#cutSets["TightRami"]["dz"]           =   [0.323 ,0.571  ] 
#cutSets["TightRami"]["missingHits"]  =   [2.1   ,1.1    ] 
    



cutSets["Veto"]={}
cutSets["Veto"]["sieie"]        =   [0.0118,0.033]
cutSets["Veto"]["1_e_1_p"]      =   [9999 ,9999 ] 
cutSets["Veto"]["dEta"]         =   [0.0155,0.013]
cutSets["Veto"]["dPhi"]         =   [9999 ,9999 ] 
cutSets["Veto"]["combiso_pt"]   =   [0.20 ,0.20 ] 
cutSets["Veto"]["missingHits"]  =   [2    ,3    ] 
cutSets["Veto"]["h_e"]          =   [9999 ,0.100] 
cutSets["Veto"]["d0"]           =   [0.12 ,0.30 ] 
cutSets["Veto"]["dz"]           =   [9999 ,9999 ]

cutSets["Loose"]={}
cutSets["Loose"]["sieie"]=[ 0.0106,  0.0305]
cutSets["Loose"]["1_e_1_p"]=[9999,  9999]
cutSets["Loose"]["dEta"]=[ 0.0102,  0.0100]
cutSets["Loose"]["missingHits"]=[ 2.0000,  2.0000]
cutSets["Loose"]["combiso_pt"]=[ 0.1300,  0.1300]
cutSets["Loose"]["h_e"]=[ 0.0992,  0.0803]
cutSets["Loose"]["dz"]=[ 0.0600,  9999]  
cutSets["Loose"]["dPhi"]=[ 0.0880,  9999]
cutSets["Loose"]["d0"]=[ 0.0300,  0.1240]

cutSets["Medium"]={}
cutSets["Medium"]["sieie"]=[ 0.0102,  0.0290]
cutSets["Medium"]["1_e_1_p"]=[ 0.0295,  0.0355]
cutSets["Medium"]["dEta"]=[ 0.0094,  0.0088]
cutSets["Medium"]["missingHits"]=[ 2.0000,  2.0000]
cutSets["Medium"]["combiso_pt"]=[ 0.0900,  0.0900]
cutSets["Medium"]["h_e"]=[ 0.0530,  0.0593]
cutSets["Medium"]["dz"]=[ 0.0320,  0.1750]
cutSets["Medium"]["dPhi"]=[ 0.0310,  0.0510]
cutSets["Medium"]["d0"]=[ 0.0150,  0.0720]

cutSets["Tight"]={}
cutSets["Tight"]["sieie"]       =   [0.0101,0.0285]
cutSets["Tight"]["1_e_1_p"]     =   [0.013,0.017]  
cutSets["Tight"]["dEta"]        =   [0.0065,0.0090]
cutSets["Tight"]["dPhi"]        =   [0.018,0.026]  
cutSets["Tight"]["combiso_pt"]  =   [0.08 ,0.080]  
cutSets["Tight"]["missingHits"] =   [2    ,2    ]  
cutSets["Tight"]["h_e"]         =   [0.040,0.050]  
cutSets["Tight"]["d0"]          =   [0.010,0.020]  
cutSets["Tight"]["dz"]          =   [0.015,0.13 ]  



varToStr={}
varToStr["pt"]         = "pt"
varToStr["eta"]        = "etaSC"
varToStr["rho"]        = "rho"
varToStr["sieie"]      = "full5x5_sigmaIetaIeta"
varToStr["1_e_1_p"]    = "ooEmooP"
varToStr["dEta"]       = "abs(dEtaIn)"
varToStr["dPhi"]       = "abs(dPhiIn)"
varToStr["combiso_pt"] = "relIsoWithEA"
varToStr["missingHits"]= "expectedMissingInnerHits"
varToStr["h_e"]        = "hOverE"
varToStr["d0"]         = "abs(d0)"
varToStr["dz"]         = "abs(dz)"
varToStr["type"]       = "sampleType"


catMap=["Barrel","Endcaps"]

catToStr={}
catToStr["Barrel"]="abs(etaSC)<1.4442"
catToStr["Endcaps"]="abs(etaSC)>1.566"

sampleMap={}
sampleMap["DY"]="sampleType==10&&isTrueEle==1"
sampleMap["TTFake"]="sampleType==20&&(isTrueEle==0||isTrueEle==3)"
#sampleMap["GJet"]="sampleType==30&&(isTrueEle==0||isTrueEle==3)"


preselection="pt>20&&passConversionVeto==1&&abs(dz)<1.&&(abs(etaSC)<1.4442||abs(etaSC)>1.566)"

colors={}
colors["Veto"]=600
colors["Loose"]=416
colors["Medium"]=632
colors["Tight"]=1
colors["VetoBKG"]=602
colors["LooseBKG"]=418
colors["MediumBKG"]=634
colors["TightBKG"]=1
colors["VetoRami"]=600
colors["LooseRami"]=416
colors["MediumRami"]=632
colors["TightRami"]=1
colors["VetoRamiBKG"]=602
colors["LooseRamiBKG"]=418
colors["MediumRamiBKG"]=634
colors["TightRamiBKG"]=1

cutStrings={}
for sample in sampleMap:
    for cutSet in cutSets:
        for var in cutSets[cutSet]:
            iCat=0
            for cut in cutSets[cutSet][var]:
                strKey=cutSet+"_"+catMap[iCat]+"_"+sample
                if not cutStrings.has_key(strKey):
                    cutStrings[strKey]=preselection+"&&"+catToStr[catMap[iCat]]
                    cutStrings[strKey]=cutStrings[strKey]+"&&"+sampleMap[sample]
                cutStrings[strKey]=cutStrings[strKey]+"&&"+varToStr[var]+"<"+str(cutSets[cutSet][var][iCat])
                iCat=iCat+1

for strName in cutStrings:
    print strName
    print cutStrings[strName]

can=ROOT.TCanvas("can","",700,700)

hists={}
leg=ROOT.TLegend(0.40,0.22,0.60,0.60)
for cat in catMap:
    for sample in sampleMap:
        samplePre=preselection+"&&"+sampleMap[sample]+"&&"+catToStr[cat]
        hists[sample+"_eta_"+cat+"_presel"]=ROOT.TH1F(sample+"_eta_"+cat+"_presel",";"+sample+" #eta "+cat+"_presel",60,-2.566,2.566)
        hists[sample+"_pt_"+cat+"_presel"]=ROOT.TH1F(sample+"_pt_"+cat+"_presel",";"+sample+" P_{T} "+cat+"_presel",50,20,160)
        hists[sample+"_rho_"+cat+"_presel"]=ROOT.TH1F(sample+"_rho_"+cat+"_presel",";"+sample+" #rho "+cat+"_presel",80,0,60)
        tree.Draw("etaSC>>"+sample+"_eta_"+cat+"_presel",samplePre)
        tree.Draw("pt>>"+sample+"_pt_"+cat+"_presel",samplePre)
        tree.Draw("rho>>"+sample+"_rho_"+cat+"_presel",samplePre)
        for cutSet in cutSets:
            strKey=cutSet+"_"+cat+"_"+sample
            colorKey=cutSet
            if sample!="DY":
                colorKey=cutSet+"BKG"
            print "making",strKey,"plots",
            hists[sample+"_eta_"+cat+"_"+cutSet]=ROOT.TH1F(sample+"_eta_"+cat+"_"+cutSet,";"+sample+" #eta "+cat+"_"+cutSet,60,-2.566,2.566)
            hists[sample+"_pt_"+cat+"_"+cutSet]=ROOT.TH1F(sample+"_pt_"+cat+"_"+cutSet,";"+sample+" P_{T} "+cat+"_"+cutSet,50,20,160)
            hists[sample+"_rho_"+cat+"_"+cutSet]=ROOT.TH1F(sample+"_rho_"+cat+"_"+cutSet,";"+sample+" #rho "+cat+"_"+cutSet,80,0,60)
            tree.Draw("etaSC>>"+sample+"_eta_"+cat+"_"+cutSet,cutStrings[strKey])
            tree.Draw("pt>>"+sample+"_pt_"+cat+"_"+cutSet,cutStrings[strKey])
            tree.Draw("rho>>"+sample+"_rho_"+cat+"_"+cutSet,cutStrings[strKey])
            print hists[sample+"_eta_"+cat+"_"+cutSet].Integral()
            hists[sample+"_eta_"+cat+"_"+cutSet+"_eff"]=ROOT.TGraphAsymmErrors( hists[sample+"_eta_"+cat+"_"+cutSet],hists[sample+"_eta_"+cat+"_presel"])
            hists[sample+"_pt_"+cat+"_"+cutSet+"_eff"]=ROOT.TGraphAsymmErrors( hists[sample+"_pt_"+cat+"_"+cutSet],hists[sample+"_pt_"+cat+"_presel"])
            hists[sample+"_rho_"+cat+"_"+cutSet+"_eff"]=ROOT.TGraphAsymmErrors( hists[sample+"_rho_"+cat+"_"+cutSet],hists[sample+"_rho_"+cat+"_presel"])
            hists[sample+"_eta_"+cat+"_"+cutSet+"_eff"].SetMarkerColor(colors[colorKey])
            hists[sample+"_pt_"+cat+"_"+cutSet+"_eff"].SetMarkerColor(colors[colorKey])
            hists[sample+"_rho_"+cat+"_"+cutSet+"_eff"].SetMarkerColor(colors[colorKey])
            hists[sample+"_eta_"+cat+"_"+cutSet+"_eff"].SetLineColor(colors[colorKey])
            hists[sample+"_pt_"+cat+"_"+cutSet+"_eff"].SetLineColor(colors[colorKey])
            hists[sample+"_rho_"+cat+"_"+cutSet+"_eff"].SetLineColor(colors[colorKey])
            hists[sample+"_eta_"+cat+"_"+cutSet+"_eff"].SetMarkerStyle(21)
            hists[sample+"_pt_"+cat+"_"+cutSet+"_eff"].SetMarkerStyle(21)
            hists[sample+"_rho_"+cat+"_"+cutSet+"_eff"].SetMarkerStyle(21)
            if cat=="Barrel":
                leg.AddEntry(hists[sample+"_eta_"+cat+"_"+cutSet+"_eff"],colorKey,"pl")

mg_eta=ROOT.TMultiGraph()
mg_pt_eb=ROOT.TMultiGraph()
mg_pt_ee=ROOT.TMultiGraph()
mg_rho_eb=ROOT.TMultiGraph()
mg_rho_ee=ROOT.TMultiGraph()
can.cd()
iCat=0
for cat in catMap:
    for sample in sampleMap:
        for cutSet in cutSets:
            mg_eta.Add(hists[sample+"_eta_"+cat+"_"+cutSet+"_eff"])
            if iCat==0:
                mg_pt_eb.Add(hists[sample+"_pt_"+cat+"_"+cutSet+"_eff"])
                mg_rho_eb.Add(hists[sample+"_rho_"+cat+"_"+cutSet+"_eff"])
            else:
                mg_pt_ee.Add(hists[sample+"_pt_"+cat+"_"+cutSet+"_eff"])
                mg_rho_ee.Add(hists[sample+"_rho_"+cat+"_"+cutSet+"_eff"])
    iCat=iCat+1

mg_eta.Draw("AP");
mg_eta.SetTitle(";#eta_{SC};Efficiency")
leg.Draw("same")
can.Update()
can.SaveAs("eta_WP_comparison.png")

mg_pt_eb.Draw("AP");
mg_pt_eb.SetTitle("Barrel Electrons;P_{T} GeV;Efficiency")
leg.Draw("same")
can.Update()
can.SaveAs("pt_eb_WP_comparison.png")

mg_pt_ee.Draw("AP");
mg_pt_ee.SetTitle("Endcap Electrons;P_{T} GeV;Efficiency")
leg.Draw("same")
can.Update()
can.SaveAs("pt_ee_WP_comparison.png")

mg_rho_eb.Draw("AP");
mg_rho_eb.SetTitle("Barrel Electrons;#rho;Efficiency")
#leg.Draw("same")
can.Update()
can.SaveAs("rho_eb_WP_comparison.png")

mg_rho_ee.Draw("AP");
mg_rho_ee.SetTitle("Endcap Electrons;#rho;Efficiency")
#leg.Draw("same")
can.Update()
can.SaveAs("rho_ee_WP_comparison.png")



