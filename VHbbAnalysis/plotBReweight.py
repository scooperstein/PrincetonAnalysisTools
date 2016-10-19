import ROOT
import sys,os

ROOT.gROOT.SetBatch(ROOT.kTRUE)

filename=sys.argv[1]
tfile=ROOT.TFile.Open(filename)
tree=tfile.Get("tree")

pathName=filename.split(".r")[0][filename.rfind('/')+1:]
#pathName=filename.split(".")[0].replace("/","_")

print pathName
if not os.path.exists(pathName):
    os.mkdir(pathName)

tree.SetBranchStatus("*", 0)
tree.SetBranchStatus("*hJets_btagCSV_*", 1)
tree.SetBranchStatus("*eight*", 1)
tree.SetBranchStatus("*isW*", 1)
tree.SetBranchStatus("*Sample*", 1)
tree.SetBranchStatus("*selLeptons_relIso_0*", 1)

can=ROOT.TCanvas("can","",700,700)

ddsystNames=["bTagWeight","bTagWeight_HFDown","bTagWeight_HFUp","bTagWeight_HFStats1Down","bTagWeight_HFStats1Up","bTagWeight_HFStats2Down","bTagWeight_HFStats2Up","bTagWeight_JESDown","bTagWeight_JESUp","bTagWeight_LFDown","bTagWeight_LFUp","bTagWeight_LFStats1Down","bTagWeight_LFStats1Up","bTagWeight_LFStats2Down","bTagWeight_LFStats2Up","bTagWeight_cErr1Down","bTagWeight_cErr1Up","bTagWeight_cErr2Down","bTagWeight_cErr2Up"]
upOrDown=["Down","Up"]
systNames=["HF","HFStats1","HFStats2","JES","LF","LFStats1","LFStats2","cErr1","cErr2"]


for CS in range(1,4):
    CSVHists={}
    CSVHists["Nominal"]=ROOT.TH1F("Nominal","Nominal CSV",20,0,1)
    tree.Draw("hJets_btagCSV_1>>Nominal","(controlSample=="+str(CS)+"&&(hJets_btagCSV_0>0.935||controlSample==3))*weight")
    #tree.Draw("hJets_btagCSV_1>>Nominal","(controlSample=="+str(CS)+")*weight")
    for pm in upOrDown:
        for systName in systNames:
            fullName="bTagWeight_"+systName+pm
            CSVHists[fullName]=ROOT.TH1F(fullName,fullName,20,0,1)
            tree.Draw("hJets_btagCSV_1>>"+fullName,"(controlSample=="+str(CS)+"&&(hJets_btagCSV_0>0.935||controlSample==3))*weight/bTagWeight*"+fullName)
            print "Filled",fullName
    
    leg=ROOT.TLegend(0.6,0.5,0.9,0.9)
    colors=[633,417,921,601,402,617,434,801,881]
    
    CSVHists["Nominal"].SetLineWidth(4)
    CSVHists["Nominal"].SetMaximum(CSVHists["Nominal"].GetMaximum()*2)
    for pm in upOrDown:
        leg1=ROOT.TLegend(0.2,0.75,0.5,0.95)
        leg2=ROOT.TLegend(0.6,0.75,0.9,0.95)
        iHist=0
        CSVHists["Nominal"].Draw("hist")
        CSVHists["Nominal"].SetTitle(";CSV 2;Weighted MC")
        leg1.AddEntry(CSVHists["Nominal"],"Nominal","l")
        iHist=iHist+1
        for systName in systNames:
            histName="bTagWeight_"+systName+pm
            print iHist, 
            CSVHists[histName].SetLineWidth(2)
            CSVHists[histName].SetLineColor(colors[iHist%len(colors)])
            CSVHists[histName].Draw("histsame")
            if iHist<len(systNames)/2+1:
                leg1.AddEntry(CSVHists[histName],histName,"l")
            else:
                leg2.AddEntry(CSVHists[histName],histName,"l")
            iHist=iHist+1
        
        leg1.Draw("same")
        leg2.Draw("same")
        can.Update()
        can.SaveAs(pathName+"_"+pm+"_CS"+str(CS)+".png")
    
        iHist=0
        for systName in systNames:
            histName="bTagWeight_"+systName+pm
            CSVHists[histName].Divide(CSVHists["Nominal"])
            CSVHists[histName].SetTitle(";CSV 2;Syst / Nominal")
            if iHist==0:
                CSVHists[histName].Draw("hist")
                CSVHists[histName].SetMinimum(0.5)
                CSVHists[histName].SetMaximum(1.5)
            else:
                CSVHists[histName].Draw("histsame")
            iHist=iHist+1
        
        leg1.Draw("same")
        leg2.Draw("same")
        can.Update()
        can.SaveAs(pathName+"_"+pm+"_relative_CS"+str(CS)+".png")
    

    
