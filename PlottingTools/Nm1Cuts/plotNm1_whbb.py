import ROOT
import sys,os
import numpy,array
import math

ROOT.gROOT.SetBatch(True)

filename=sys.argv[1]
tfile=ROOT.TFile(filename)

treename="tree"
tree=tfile.Get(treename)

outDir="variablePlots_WHbb_5E3"
suffix="analysis"

if not os.path.exists(outDir):
    os.makedirs(outDir)

flatTree=True

varsmap={}
def AssignVarsToMap():
    varsmap["jet1pt"]         = tree.hJets_pt_0
    varsmap["jet2pt"]         = tree.hJets_pt_1
    varsmap["jet1eta"]        = tree.hJets_eta_0
    varsmap["jet2eta"]        = tree.hJets_eta_1
    varsmap["jet1csv"]        = tree.hJets_btagCSV_0
    varsmap["jet2csv"]        = tree.hJets_btagCSV_1
    varsmap["ptjj"]           = tree.H_pt
    varsmap["ptW"]            = tree.V_pt
    varsmap["met"]            = tree.met_pt
    varsmap["HVdPhi"]         = tree.HVdPhi
    varsmap["lepID"]          = tree.selLeptons_tightId_0
    varsmap["leppt"]          = tree.selLeptons_pt_0
    varsmap["lepeta"]         = tree.selLeptons_eta_0
    varsmap["lepmetdphi"]     = tree.lepMetDPhi
    varsmap["naddlep"]        = tree.nAddLeptons
    varsmap["naddjet"]        = tree.nAddJets252p9_puid
    varsmap["topmass_nomet"]  = tree.Top1_mass_fromLepton
    varsmap["isWenu"]         = tree.isWenu
    varsmap["isWmunu"]        = tree.isWmunu
    varsmap["type"]           = tree.sampleIndex

yieldsTextFile=open(outDir+"/yields.txt","w+")


cutsets={}

cutsets["Analysis"] = {}
cutsets["Analysis"]["jet1pt"]        =   [30, 30]
cutsets["Analysis"]["jet2pt"]        =   [30, 30]
cutsets["Analysis"]["jet1eta"]       =   [2.5, 2.5]
cutsets["Analysis"]["jet2eta"]       =   [2.5, 2.5]
cutsets["Analysis"]["jet1csv"]       =   [0.9, 0.9]
cutsets["Analysis"]["jet2csv"]       =   [0.5, 0.5]
cutsets["Analysis"]["ptjj"]          =   [100, 100]
cutsets["Analysis"]["ptW"]           =   [100, 100]
cutsets["Analysis"]["met"]           =   [45, 45]
cutsets["Analysis"]["HVdPhi"]        =   [2.5, 2.5]
cutsets["Analysis"]["lepID"]         =   [1, 1]
cutsets["Analysis"]["leppt"]         =   [30, 30]
cutsets["Analysis"]["lepeta"]        =   [2.4, 2.5]
cutsets["Analysis"]["lepmetdphi"]    =   [1.571, 1.571]
cutsets["Analysis"]["naddlep"]       =   [1, 1]
cutsets["Analysis"]["naddjet"]       =   [1, 1]
cutsets["Analysis"]["topmass_nomet"] =   [0, 0]


cutsets["presel"] = {}
cutsets["presel"]["jet1pt"]        =   [15, 15]
cutsets["presel"]["jet2pt"]        =   [15, 15]
cutsets["presel"]["jet1eta"]       =   [2.5, 2.5]
cutsets["presel"]["jet2eta"]       =   [2.5, 2.5]
cutsets["presel"]["jet1csv"]       =   [0.5, 0.5]
cutsets["presel"]["jet2csv"]       =   [0.0, 0.0]
cutsets["presel"]["ptjj"]          =   [50, 50]
cutsets["presel"]["ptW"]           =   [50, 50]
cutsets["presel"]["met"]           =   [30, 30]
cutsets["presel"]["HVdPhi"]        =   [0.0, 0.0]
cutsets["presel"]["lepID"]         =   [0, 0]
cutsets["presel"]["leppt"]         =   [22, 30]
cutsets["presel"]["lepeta"]        =   [2.4, 2.5]
cutsets["presel"]["lepmetdphi"]    =   [3.14, 3.14]
cutsets["presel"]["naddlep"]       =   [4, 4]
cutsets["presel"]["naddjet"]       =   [10, 10]
cutsets["presel"]["topmass_nomet"] =   [0, 0]

nBins=40
varstocut={}  # 0 keep less than, 1 keept greater than
varstocut["jet1pt"]     = [1,nBins,cutsets["presel"]["jet1pt"][0], 100, nBins, cutsets["presel"]["jet1pt"][1], 100 ]
varstocut["jet2pt"]     = [1,nBins,cutsets["presel"]["jet2pt"][0], 100, nBins, cutsets["presel"]["jet2pt"][1], 100 ]
varstocut["jet1eta"]     = [0,nBins,0., cutsets["presel"]["jet1eta"][0], nBins, 0., cutsets["presel"]["jet1eta"][1] ]
varstocut["jet2eta"]     = [0,nBins,0., cutsets["presel"]["jet2eta"][0], nBins, 0., cutsets["presel"]["jet2eta"][1] ]
varstocut["jet1csv"]     = [1,nBins,cutsets["presel"]["jet1csv"][0], 1.0, nBins, cutsets["presel"]["jet1csv"][1], 1.0 ]
varstocut["jet2csv"]     = [1,nBins,cutsets["presel"]["jet2csv"][0], 1.0, nBins, cutsets["presel"]["jet2csv"][1], 1.0 ]
varstocut["ptjj"]        = [1,nBins,cutsets["presel"]["ptjj"][0],200,nBins,cutsets["presel"]["ptjj"][1], 200]
varstocut["ptW"]         = [1,nBins,cutsets["presel"]["ptW"][0],200,nBins,cutsets["presel"]["ptW"][1], 200]
varstocut["met"]         = [1,nBins,cutsets["presel"]["met"][0],100,nBins,cutsets["presel"]["met"][1], 100]
varstocut["HVdPhi"]      = [1,nBins,cutsets["presel"]["HVdPhi"][0], 3.14, nBins,cutsets["presel"]["HVdPhi"][1],3.14]
varstocut["lepID"]       = [1,4,0,4,4,0,4]
varstocut["leppt"]       = [1,nBins,cutsets["presel"]["leppt"][0], 100,nBins,cutsets["presel"]["leppt"][1],100]
varstocut["lepeta"]      = [0,nBins,0.,cutsets["presel"]["lepeta"][0],nBins,0.,cutsets["presel"]["lepeta"][1]]
varstocut["lepmetdphi"]  = [0,nBins,0.,cutsets["presel"]["lepmetdphi"][0],nBins,0.,cutsets["presel"]["lepmetdphi"][1]]
varstocut["naddlep"]     = [0,10,0,10,10,0,10]
varstocut["naddjet"]     = [0,10,0,10,10,0,10]
varstocut["topmass_nomet"] = [1,nBins,cutsets["presel"]["topmass_nomet"][0],500,nBins,cutsets["presel"]["topmass_nomet"][1],500]



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




absCutVars=["jet1eta","jet2eta"]
def passVarCut(var,index,wp):
    if flatTree:
        varValue=varsmap[var]
    else:
        varValue=varsmap[var][index]

    if var in absCutVars:
        varValue=math.fabs(varValue)
 
    catInd=category(nCats,index)
    if varstocut[var][0]==0:
        if varValue > cutsets[wp][var][catInd]:
            return False
    
    if varstocut[var][0]==1:
        if varValue < cutsets[wp][var][catInd]:
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
makeTree=False

#eleHoE_PassNm1        = numpy.zeros(1, dtype=int)
#b_eleHoE_PassNm1=tree.Branch( 'eleHoE_PassNm1', eleHoE_PassNm1, 'eleHoE_PassNm1/I' )

#wps=["preselglobal","Veto","Loose","Medium"]
wps=["Analysis"]

potentialCuts={}  #wp,cat,var

numIterations=5


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
                kinplots[mctype+"_jet1eta_"+catMap[cat]+"_"+wp]=ROOT.TH1F(mctype+"_jet1eta_"+catMap[cat]+"_"+wp,";"+mctype+" #eta "+catMap[cat]+"_"+wp,60,-2.8,2.8)
                kinplots[mctype+"_jet2eta_"+catMap[cat]+"_"+wp]=ROOT.TH1F(mctype+"_jet2eta_"+catMap[cat]+"_"+wp,";"+mctype+" #eta "+catMap[cat]+"_"+wp,60,-2.8,2.8)
                kinplots[mctype+"_jet1pt_"+catMap[cat]+"_"+wp]=ROOT.TH1F(mctype+"_jet1pt_"+catMap[cat]+"_"+wp,";"+mctype+" P_{T} "+catMap[cat]+"_"+wp,50,10,160)
                kinplots[mctype+"_jet2pt_"+catMap[cat]+"_"+wp]=ROOT.TH1F(mctype+"_jet2pt_"+catMap[cat]+"_"+wp,";"+mctype+" P_{T} "+catMap[cat]+"_"+wp,50,10,160)
                kinplots[mctype+"_jet1csv_"+catMap[cat]+"_"+wp]=ROOT.TH1F(mctype+"_jet1csv_"+catMap[cat]+"_"+wp,";"+mctype+" P_{T} "+catMap[cat]+"_"+wp,50,0,1.0)
                kinplots[mctype+"_jet2csv_"+catMap[cat]+"_"+wp]=ROOT.TH1F(mctype+"_jet2csv_"+catMap[cat]+"_"+wp,";"+mctype+" P_{T} "+catMap[cat]+"_"+wp,50,0,1.0)
        
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
            

            if varsmap["type"]<0:
                mctype="VH"
            elif varsmap["type"]>0:
                mctype="Bkg"
            else:
                continue
                print "sampleIndex",varsmap["type"]
                print "I don't want to live in a world of unknowns!!"


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
                        hists[mctype+"_"+catMap[category(nCats,hltindex)]+"_"+var+"_nm1_"+wp].Fill(varsmap[var],tree.weight)
                    kinplots[mctype+"_jet1eta_"+catMap[category(nCats,hltindex)]+"_"+wp].Fill(varsmap['jet1eta'],tree.weight)
                    kinplots[mctype+"_jet2eta_"+catMap[category(nCats,hltindex)]+"_"+wp].Fill(varsmap['jet2eta'],tree.weight)
                    kinplots[mctype+"_jet1pt_"+catMap[category(nCats,hltindex)]+"_"+wp].Fill(varsmap['jet1pt'],tree.weight)
                    kinplots[mctype+"_jet2pt_"+catMap[category(nCats,hltindex)]+"_"+wp].Fill(varsmap['jet2pt'],tree.weight)
                    kinplots[mctype+"_jet1csv_"+catMap[category(nCats,hltindex)]+"_"+wp].Fill(varsmap['jet1csv'],tree.weight)
                    kinplots[mctype+"_jet2csv_"+catMap[category(nCats,hltindex)]+"_"+wp].Fill(varsmap['jet2csv'],tree.weight)
                #else:
                #    for var in varstocut:
                #        hists[mctype+"_"+catMap[category(nCats,hltindex)]+"_"+var+"_nm1_"+wp].Fill(varsmap[var][hltindex],tree.weight)
                #    kinplots[mctype+"_eta_"+catMap[category(nCats,hltindex)]+"_"+wp].Fill(tree.hltEleEta[hltindex],tree.weight)
                #    kinplots[mctype+"_pt_"+catMap[category(nCats,hltindex)]+"_"+wp].Fill(tree.hltElePt[hltindex],tree.weight)
                #    kinplots[mctype+"_rho_"+catMap[category(nCats,hltindex)]+"_"+wp].Fill(tree.hltElePt[hltindex],tree.weight)
        
            if len(cutsFailed)==1:
                if flatTree:
                    hists[mctype+"_"+catMap[category(nCats,hltindex)]+"_"+cutsFailed[0]+"_nm1_"+wp].Fill(varsmap[cutsFailed[0]],tree.weight)
                else:
                    hists[mctype+"_"+catMap[category(nCats,hltindex)]+"_"+cutsFailed[0]+"_nm1_"+wp].Fill(varsmap[cutsFailed[0]][hltindex],tree.weight)
                    
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
                    hists["sob_"+catMap[icat]+"_"+var+"_nm1_"+wp].Divide(bkgtemp)
               
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
                        hists[mctype+"_"+catMap[icat]+"_"+var+"_nm1_"+wp].SetMaximum(maxval)
                        hists[mctype+"_"+catMap[icat]+"_"+var+"_nm1_"+wp].SetTitle(";"+var+" "+catMap[icat]+" "+wp)
                        hists[mctype+"_"+catMap[icat]+"_"+var+"_nm1_"+wp].Draw()
                        first=False
                    else:
                        hists[mctype+"_"+catMap[icat]+"_"+var+"_nm1_"+wp].Draw("sames")
                    leg.AddEntry(hists[mctype+"_"+catMap[icat]+"_"+var+"_nm1_"+wp],mctype,"l")
   
                leg.Draw("same") 
                line=ROOT.TLine(cutsets[wp][var][icat],0,cutsets[wp][var][icat], maxval/1.3)
                line.SetLineColor(1)    
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
                    hists["sob_"+catMap[icat]+"_"+var+"_nm1_"+wp].SetMinimum(0.1)
                    temp=hists["sob_"+catMap[icat]+"_"+var+"_nm1_"+wp].Clone()
                    line=ROOT.TLine(cutsets[wp][var][icat],0,cutsets[wp][var][icat], hists["sob_"+catMap[icat]+"_"+var+"_nm1_"+wp].GetMaximum()/1.3)
                    line.SetLineColor(1)    
                    line.SetLineWidth(2)
                    line.Draw("same")
        
                    hLineVal=[0.005,0.005]
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
                    if (varstocut[var][0] == 0):
                        # cut is less than
                        for iBin in range(hists["sob_"+catMap[icat]+"_"+var+"_nm1_"+wp].GetNbinsX()-1, 0, -1):
                            lowEdge=hists["sob_"+catMap[icat]+"_"+var+"_nm1_"+wp].GetBinLowEdge(iBin) + hists["sob_"+catMap[icat]+"_"+var+"_nm1_"+wp].GetBinWidth(iBin)
                            binVal=hists["sob_"+catMap[icat]+"_"+var+"_nm1_"+wp].GetBinContent(iBin)
                            if binVal > hLineVal[icat]:
                                potentialCuts[wpPlus1][var][catMap[icat]].append(lowEdge)
                    elif (varstocut[var][0] == 1):
                        # cut is greater than
                        for iBin in range(1, hists["sob_"+catMap[icat]+"_"+var+"_nm1_"+wp].GetNbinsX()):
                            lowEdge=hists["sob_"+catMap[icat]+"_"+var+"_nm1_"+wp].GetBinLowEdge(iBin)
                            binVal=hists["sob_"+catMap[icat]+"_"+var+"_nm1_"+wp].GetBinContent(iBin)
                            if binVal > hLineVal[icat]:
                                potentialCuts[wpPlus1][var][catMap[icat]].append(lowEdge)
                            
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
                print kinplots[mctype+"_jet1eta_"+catMap[cat]+"_"+wp].Integral(),

            
            for ikin in ["jet1pt","jet2pt","jet1eta","jet2eta","jet1csv","jet2csv"]:
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
                        if ikin=="jet1eta":
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


