import ROOT
import sys
import math

filename=sys.argv[1]
tfile=ROOT.TFile(filename)
tree=tfile.Get("newtree")

tree.SetBranchStatus("*",0)
tree.SetBranchStatus("run",1)
tree.SetBranchStatus("LS",1)
tree.SetBranchStatus("hasBrilData",1)
tree.SetBranchStatus("hasCMSData",1)

runsToCheck=[246908,246919,246920,246923,246926,246930,246936,246951,246960,247068,247070,247073,247078,247079,247081,247252,247253,247262,247267,247377,247381]
missingRuns=[246908,246919,246920,246923,246926,246930,246936,246951,246960,247068,247070,247073,247078,247079,247081,247252,247253,247262,247267,247377,247381]

onlyBril=[]
onlyCMS=[]
bothSets=[]

runLSMax={}

nentries=tree.GetEntries()
for ient in range(nentries):
    tree.GetEntry(ient)
    if tree.run in missingRuns:
        missingRuns.remove(tree.run)
    if not runLSMax.has_key(tree.run):
        runLSMax[tree.run]=tree.LS
    elif runLSMax[tree.run]<tree.LS:
        runLSMax[tree.run]=tree.LS

    if tree.hasBrilData and tree.hasCMSData:
        bothSets.append(tree.run)
    elif tree.hasBrilData:
        onlyBril.append(tree.run)
    elif tree.hasCMSData:
        onlyCMS.append(tree.run)

print "missingRuns runs: ",missingRuns
print "Not both CMS and Bril: ",list(set(runsToCheck)-set(bothSets))

    
tree.SetBranchStatus("*",1)
histpix={}
histbest={}
hist2ndbest={}
histlayers={}

bestHF={}
bestBCM1f={}

#filling twice ????  FIXME!!!
for ient in range(nentries):
    tree.GetEntry(ient)
    if tree.run in bothSets:
        if not histpix.has_key(tree.run):
            histpix[tree.run]=ROOT.TH1F(str(tree.run),";Luminosity Section;PCC/BestLumi",runLSMax[tree.run],0,runLSMax[tree.run])
            histpix[tree.run].GetXaxis().SetTitleSize(0.07)
            histpix[tree.run].GetXaxis().SetTitleOffset(0.6)
            histpix[tree.run].GetYaxis().SetTitleSize(0.07)
            histpix[tree.run].GetYaxis().SetTitleOffset(0.3)
            for layer in range(0,5):
                layerkey=str(tree.run)+"_layer"+str(layer+1)
                histlayers[layerkey]=ROOT.TH1F(str(tree.run)+"_layer"+str(layer+1),";Luminosity Section;PCC/BestLumi",runLSMax[tree.run],0,runLSMax[tree.run])
                histlayers[layerkey].GetXaxis().SetTitleSize(0.10)
                histlayers[layerkey].GetXaxis().SetTitleOffset(0.3)
                histlayers[layerkey].GetYaxis().SetTitleSize(0.12)
                histlayers[layerkey].GetYaxis().SetTitleOffset(0.3)

        histpix[tree.run].Fill(tree.LS,tree.pixel_xsec)
        for layer in range(0,5):
            histlayers[str(tree.run)+"_layer"+str(layer+1)].Fill(tree.LS,tree.pixel_xsec_layers[layer])

    if tree.hasBrilData:
        key=tree.run
        if not histbest.has_key(key):
            histbest[key]=ROOT.TH1F(str(tree.run)+"best",";Luminosity Section;Integrated Luminosity(ub^{-1})",runLSMax[tree.run],0,runLSMax[tree.run])
            hist2ndbest[key]=ROOT.TH1F(str(tree.run)+"2ndbest",";Luminosity Section;Integrated Luminosity(ub^{-1})",runLSMax[tree.run],0,runLSMax[tree.run])
            histbest[key].GetXaxis().SetTitleSize(0.07)
            histbest[key].GetXaxis().SetTitleOffset(0.6)
            histbest[key].GetYaxis().SetTitleSize(0.07)
            histbest[key].GetYaxis().SetTitleOffset(0.3)
            hist2ndbest[key].GetXaxis().SetTitleSize(0.07)
            hist2ndbest[key].GetXaxis().SetTitleOffset(0.6)
            hist2ndbest[key].GetYaxis().SetTitleSize(0.07)
            hist2ndbest[key].GetYaxis().SetTitleOffset(0.3)
        if tree.BestLumi>0:
            histbest[key].Fill(tree.LS,tree.BestLumi)
        if tree.HFLumi>0 and tree.BCMFLumi>0 and tree.BestLumi>0:
            if tree.run not in bestHF.keys():
                bestHF[tree.run]=0
                bestBCM1f[tree.run]=0
            if math.fabs(tree.HFLumi-tree.BestLumi) < math.fabs(tree.BCMFLumi-tree.BestLumi):
                bestHF[tree.run]=bestHF[tree.run]+1
                hist2ndbest[key].Fill(tree.LS,tree.BCMFLumi)
            else:
                bestBCM1f[tree.run]=bestBCM1f[tree.run]+1
                hist2ndbest[key].Fill(tree.LS,tree.HFLumi)
           



tcan=ROOT.TCanvas("tcan","",1400,800)
padlumis =ROOT.TPad("padlumis", "",0.0,0.0,0.5,1.0)
padlayers=ROOT.TPad("padlayers","",0.5,0.0,1.0,1.0)

padlumis.Draw()        
padlayers.Draw()        

print bestHF
print bestBCM1f


for run in runsToCheck:
    padlumis.Divide(1,3)
    print run
    if run in histpix.keys():
        padlumis.cd(1)
        label=ROOT.TText(0,histpix[run].GetMaximum()*1.05,"   Pixel Cross Section - Run="+str(tree.run))
        label.SetTextSize(.1)
        histpix[run].SetMaximum(histpix[run].GetMaximum()*1.2)
        histpix[run].Draw("hist")
        label.Draw("same")
        padlumis.Update()

    if run in histbest.keys():
        padlumis.cd(2)
        histbest[run].Draw("hist")
        hist2ndbest[run].SetLineColor(633)
        hist2ndbest[run].Draw("histsame")
        leg=ROOT.TLegend(0.1,0.1,0.7,0.4)
        if bestHF.has_key(run) and  bestHF[run]>bestBCM1f[run]:
            bestfrac=float(bestHF[run])/(bestHF[run]+bestBCM1f[run])
            print "hf",bestfrac
            leg.AddEntry(histbest[run],"Best Lumi HF: "+"{0:.1f}".format(bestfrac*100)+"%","l")
            leg.AddEntry(hist2ndbest[run],"2nd Best Lumi BCM1f: "+"{0:.1f}".format((1-bestfrac)*100)+"%","l")
        if bestHF.has_key(run) and  bestHF[run]<bestBCM1f[run]:
            bestfrac=1-float(bestHF[run])/(bestHF[run]+bestBCM1f[run])
            print "bcm1f",bestfrac
            leg.AddEntry(histbest[run],"Best Lumi BCM1f: "+"{0:.1f}".format(bestfrac*100)+"%","l")
            leg.AddEntry(hist2ndbest[run],"2nd Best Lumi HF: "+"{0:.1f}".format((1-bestfrac)*100)+"%","l")
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        leg.Draw("same")
        padlumis.Update()

        padlumis.cd(3)
        myclone=histbest[run].Clone("temp")
        myclone.Divide(hist2ndbest[run])
        myclone.SetTitle(";Luminosity Section;Best/2ndBest")
        label2=ROOT.TText(0,myclone.GetMaximum()*1.05,"   Best/2ndBest Lumi - Run="+str(run))
        label2.SetTextSize(.1)
        myclone.SetMaximum(myclone.GetMaximum()*1.2)
        myclone.Draw("hist")
        label2.Draw("sames")
        padlumis.Update()
        
    
    if run in bothSets:
        padlayers.Divide(1,5)
        layertexts={}
        for layer in range(5):
            padlayers.cd(layer+1)
            key=str(run)+"_layer"+str(layer+1)
            layertexts[layer]=ROOT.TText(0,histlayers[key].GetMaximum()*1.05,"    Pixel Cross Section - Layer="+str(layer+1))
            layertexts[layer].SetTextSize(.1)
            histlayers[key].SetMaximum(histlayers[key].GetMaximum()*1.2)
            histlayers[key].Draw("hist")
            layertexts[layer].Draw("same")
        padlayers.Update()    
    
    tcan.Update()
    tcan.SaveAs(str(run)+".png")
    

    raw_input()
    #tcan.Clear()
    padlumis.Clear()
    padlayers.Clear()
    tcan.Update()

