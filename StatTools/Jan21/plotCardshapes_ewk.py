import ROOT
import sys

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True)

channel = sys.argv[1]
ifile = ROOT.TFile("hists_%s.root" % channel)
data = ifile.Get("BDT_%s_data_obs" % channel)

doBlind = False
#doBlind = True
if (doBlind):
    for i in range(data.GetNbinsX()-7,data.GetNbinsX()+1):
    #for i in range(data.GetNbinsX()-4,data.GetNbinsX()+1):
        data.SetBinContent(i,0.0)
        data.SetBinError(i,0.0)
data.SetMarkerStyle(20)

tt = ifile.Get("BDT_%s_TT" % channel)
stop = ifile.Get("BDT_%s_s_Top" % channel)
zjets = ifile.Get("BDT_%s_Zjets" % channel)
wjets = ifile.Get("BDT_%s_Wjets" % channel)
qcd = ifile.Get("BDT_%s_QCD_data" % channel)
vv = ifile.Get("BDT_%s_VV" % channel)
ewkwjets = ifile.Get("BDT_%s_EWKWJets" % channel)
intewkwjets = ifile.Get("BDT_%s_IntEWKWJets" % channel)

for i in range(1,qcd.GetNbinsX()+1):
    qcd.SetBinError(i,0.0)

rbin = 2
tt.Rebin(rbin)
stop.Rebin(rbin)
zjets.Rebin(rbin)
wjets.Rebin(rbin)
qcd.Rebin(rbin)
vv.Rebin(rbin)
ewkwjets.Rebin(rbin)
intewkwjets.Rebin(rbin)
data.Rebin(rbin)

stop.SetLineColor(432)
stop.SetFillColor(432)
tt.SetLineColor(600)
tt.SetFillColor(600)
zjets.SetLineColor(394)
zjets.SetFillColor(394)
wjets.SetLineColor(418)
wjets.SetFillColor(418)
qcd.SetLineColor(616)
qcd.SetFillColor(616)
vv.SetLineColor(922)
vv.SetFillColor(922)
ewkwjets.SetLineColor(ROOT.kRed)
ewkwjets.SetFillColor(ROOT.kRed)
intewkwjets.SetLineColor(416)
intewkwjets.SetFillColor(416)

#vv.SetLineWidth(3)
#vvlf.SetLineWidth(3)
ewkwjets.SetLineWidth(3)

#tt.Scale(SF_TT)
#wjets.Scale(SF_Wjets)
#ewkwjets.Scale(30.0)
#vv.Scale(0.6017)
#vvlf.Scale(0.6017)

stack = ROOT.THStack("stack","stack")
stack.Add(intewkwjets)
stack.Add(zjets)
stack.Add(vv)
stack.Add(qcd)
stack.Add(stop)
stack.Add(tt)
stack.Add(wjets)
stack.Add(ewkwjets)

print "ttnom, ewkwjetsnom: ", tt.GetBinContent(35),ewkwjets.GetBinContent(35)
print "zjets,vv,stop,tt,wjets,qcd,ewkwjets",zjets.GetBinContent(35),vv.GetBinContent(35),stop.GetBinContent(35),tt.GetBinContent(35),wjets.GetBinContent(35),qcd.GetBinContent(35),ewkwjets.GetBinContent(35)
#vvlf.Scale(10.)
#vv.Scale(10.)
#bkg.SetMaximum(1.3*max(bkg.GetMaximum(),ewkwjets.GetMaximum()))
#bkg.SetMinimum(0.)
#stack.SetMaximum(1.3*stack.GetMaximum())
#stack.SetMinimum(0.)

#bkg.SetTitle(channel)
#bkg.GetXaxis().SetTitle("Mjj [GeV]")

#stack.SetMaximum(4500)
#stack.SetMaximum(20000)
canv = ROOT.TCanvas("canv","canv")
canv.SetBottomMargin(0.3)
xmin = 0
xmax = 2.5
if (len(sys.argv)>2):
    xmin = float(sys.argv[2])
    xmax = float(sys.argv[3])
#xmin = 200
#xmax = 3000
frame = ROOT.TH1F("frame","",1,xmin,xmax)
frame.SetStats(0)
frame.GetXaxis().SetLabelSize(0)
frame.GetXaxis().SetTitleOffset(0.91);
frame.GetYaxis().SetTitle("Events")
#frame.GetXaxis().SetTitle("atanhBDT")
frame.GetYaxis().SetLabelSize(0.04)
#frame.GetYaxis().SetRangeUser(10E-03,16000)
frame.GetYaxis().SetRangeUser(1.0,stack.GetStack().Last().GetMaximum()*1.3)
#frame.GetYaxis().SetRangeUser(10E-03,stack.GetStack().Last().GetMaximum()*100)
#frame.GetYaxis().SetRangeUser(10E-03,200)
frame.Draw()

stack.Draw("hist same")
#print "Chi2/NDOF = ",data.Chi2Test( stack.GetHistogram() , "UWCHI2/NDF")
#print "Total bkg: ",stack.GetHistogram().Integral()
stack.Draw("ep same")
stack.SetTitle(channel)
#stack.GetXaxis().SetTitle("Mjj [GeV]")
data.Draw("ep same")
#print "Total data: ",data.Integral()
#bkg.Draw("hist")
#bkg.Draw("ep same")
#stackvv.Draw("hist same")
#stackvv.Draw("ep same")
#ewkwjets.Draw("hist same")
#ewkwjets.Draw("ep same")

leg = ROOT.TLegend(0.65,0.75,0.85,0.95)
#leg = ROOT.TLegend(0.7,0.75,0.9,0.95)
#leg = ROOT.TLegend(0.65,0.3,0.85,0.5)
#leg = ROOT.TLegend(0.7,0.7,0.9,0.9)
#leg = ROOT.TLegend(0.4,0.6,0.7,0.9)
#leg = ROOT.TLegend(0.7,0.6,0.9,0.9)
leg.AddEntry(ewkwjets,"EWK W+2-jet")
leg.AddEntry(wjets,"Wjets")
leg.AddEntry(tt,"TT")
leg.AddEntry(stop,"Single Top")
leg.AddEntry(qcd,"QCD (from data)")
leg.AddEntry(vv,"VV")
leg.AddEntry(zjets,"Zjets")
leg.AddEntry(intewkwjets,"Interference")
leg.AddEntry(data,"Data")
#leg.AddEntry(bkg,"Bkg")
leg.Draw("same")

pad2 = ROOT.TPad("pad2", "pad2", 0., 0., 1., 1.)
pad2.SetTopMargin(0.73)
#	pad2.SetBottomMargin(0.)
#	pad2.SetRightMargin(0.)
pad2.SetFillColor(0)
pad2.SetFillStyle(0)
pad2.Draw()
pad2.cd()
frame2 = ROOT.TH1F("frame2","",1,xmin,xmax)
frame2.SetMinimum(-0.2)	
frame2.SetMaximum(0.2) 
frame2.GetYaxis().SetLabelSize(0.02)
frame2.GetXaxis().SetLabelSize(0.04)
frame2.GetYaxis().SetTitleSize(0.04)
frame2.GetXaxis().SetTitle("Transformed BDT Score")
#frame2.GetXaxis().SetTitle("Sub-leading Jet CMVA")
if (len(sys.argv)>4):
    frame2.GetXaxis().SetTitle(sys.argv[4])
frame2.SetStats(0)
frame2.GetYaxis().SetTitle("Data/MC - 1")	
frame2.Draw()

doJEC = 1
if (len(sys.argv)>5):
   doJEC = int(sys.argv[5])

data2=data.Clone("data2")
mc_total=stack.GetStack().Last().Clone("mc_total")
mc2_total = ewkwjets.Clone("mc2_total")
mc2_total.Divide(mc_total)
data2.Add(mc_total,-1)
data2.Divide(mc_total)

mc_total_uncUp=stack.GetStack().Last().Clone("mc_total")
mc_total_uncDown=stack.GetStack().Last().Clone("mc_total")
for i in range(0,mc_total_uncUp.GetNbinsX()): 
	e=0.
	if mc_total.GetBinContent(i+1) != 0:
		e = mc_total.GetBinError(i+1)/mc_total.GetBinContent(i+1)
	mc_total_uncUp.SetBinContent(i+1,e)
	mc_total_uncDown.SetBinContent(i+1,-e)
mc_total_uncUp.SetLineColor(ROOT.kBlack)
mc_total_uncUp.SetLineWidth(1)
mc_total_uncUp.SetFillColor(ROOT.kBlack)
mc_total_uncUp.SetFillStyle(3004)
mc_total_uncDown.SetLineColor(ROOT.kBlack)
mc_total_uncDown.SetLineWidth(1)
mc_total_uncDown.SetFillColor(ROOT.kBlack)
mc_total_uncDown.SetFillStyle(3004)

#mc2_total.Draw("hist same")
mc_total_uncUp.Draw("HIST SAME")
mc_total_uncDown.Draw("HIST SAME")


line = ROOT.TLine(xmin,0,xmax,0)
line.SetLineStyle(3)
line.Draw("same")

#data2.Draw("PEsame")
#print "Data yield: ",data.Integral()
#print "MC yield: ",stack.GetHistogram().Integral()

# now add systematic bands
mc_total_systs = []
#colors = [ROOT.kRed,ROOT.kBlue,ROOT.kGreen,ROOT.kViolet] 
if doJEC:
    colors = [ROOT.kRed,ROOT.kBlue,ROOT.kViolet] 
else:
    colors = [ROOT.kBlue,ROOT.kViolet] 
#for syst in ["CMS_ewkwjet_scale_j_13TeVUp","CMS_ewkwjet_scale_j_13TeVDown","CMS_ewkwjet_qglUp","CMS_ewkwjet_qglDown","CMS_ewkwjet_LHE_weights_scale_muR_EWKWJetsUp","CMS_ewkwjet_LHE_weights_scale_muR_EWKWJetsDown","CMS_ewkwjet_LHE_weights_scale_muF_EWKWJetsUp","CMS_ewkwjet_LHE_weights_scale_muF_EWKWJetsDown"]:
for syst in ["CMS_ewkwjet_scale_j_13TeVUp","CMS_ewkwjet_scale_j_13TeVDown","CMS_ewkwjet_qglUp","CMS_ewkwjet_qglDown","CMS_ewkwjet_LHE_weights_scale_muF_EWKWJetsUp","CMS_ewkwjet_LHE_weights_scale_muF_EWKWJetsDown"]:
        if not doJEC and syst.find("scale_j")!=-1: continue
        tt = ifile.Get("BDT_%s_TT_%s" % (channel,syst.replace("EWKWJets","TT"))).Clone("tt_"+syst)
        if syst.find("scale_mu") == -1:
            stop = ifile.Get("BDT_%s_s_Top_%s" % (channel,syst)).Clone("stop_"+syst)
        else:
            stop = ifile.Get("BDT_%s_s_Top" % (channel))
        zjets = ifile.Get("BDT_%s_Zjets_%s" % (channel,syst.replace("EWKWJets","Zjets"))).Clone("zjets_"+syst)
        wjets = ifile.Get("BDT_%s_Wjets_%s" % (channel,syst.replace("EWKWJets","Wjets"))).Clone("wjets_"+syst)
        qcd = ifile.Get("BDT_%s_QCD_data" % (channel))
        for i in range(1,qcd.GetNbinsX()+1):
            qcd.SetBinError(i,0.0)
        vv = ifile.Get("BDT_%s_VV_%s" % (channel,syst.replace("EWKWJets","VV"))).Clone("vv_"+syst)
        ewkwjets = ifile.Get("BDT_%s_EWKWJets_%s" % (channel,syst)).Clone("ewkwjets_"+syst)
        intewkwjets = ifile.Get("BDT_%s_IntEWKWJets" % (channel)).Clone("intewkwjets_"+syst)

        tt.Rebin(rbin)
        stop.Rebin(rbin)
        zjets.Rebin(rbin)
        wjets.Rebin(rbin)
        qcd.Rebin(rbin)
        vv.Rebin(rbin)
        ewkwjets.Rebin(rbin)
        intewkwjets.Rebin(rbin)
        
        stack_syst = ROOT.THStack("stack_"+syst,"stack_"+syst)
        stack_syst.Add(zjets)
        stack_syst.Add(intewkwjets)
        stack_syst.Add(vv)
        stack_syst.Add(stop)
        stack_syst.Add(tt)
        stack_syst.Add(wjets)
        stack_syst.Add(qcd)
        stack_syst.Add(ewkwjets)
        mc_total_syst = stack_syst.GetStack().Last().Clone("stack_"+syst)
        print syst
        print "tot_syst,tot_nom,tt_syst,ewkwjets_syst",stack_syst.GetStack().Last().GetBinContent(35),mc_total.GetBinContent(35),tt.GetBinContent(35),ewkwjets.GetBinContent(35)
        print "zjets,vv,stop,tt,wjets,qcd,ewkwjets",zjets.GetBinContent(35),vv.GetBinContent(35),stop.GetBinContent(35),tt.GetBinContent(35),wjets.GetBinContent(35),qcd.GetBinContent(35),ewkwjets.GetBinContent(35)
        mc_total_syst.Add(mc_total,-1)
        mc_total_syst.Divide(mc_total)
        mc_total_systs.append(mc_total_syst)
        for i in range(1,mc_total_syst.GetNbinsX()+1):
            print mc_total_syst.GetBinContent(i)

for i in range(len(mc_total_systs)):
    hist = mc_total_systs[i]
    print i/2
    hist.SetLineColor(colors[i/2])
    hist.Draw("hist same")

data2.Draw("PEsame")
        
#leg2 = ROOT.TLegend(0.2,0.11,0.5,0.15)
leg2 = ROOT.TLegend(0.2,0.32,0.5,0.42)
#leg2 = ROOT.TLegend(0.11,0.77,0.41,0.87)
if doJEC:
    leg2.AddEntry(mc_total_systs[0],"JES Up/Down")
    leg2.AddEntry(mc_total_systs[2],"QGL Corr. Up/Down")
    #leg2.AddEntry(mc_total_systs[4],"QCD Scale #mu_{R} Up/Down")
    #leg2.AddEntry(mc_total_systs[6],"QCD Scale #mu_{F} Up/Down")
    leg2.AddEntry(mc_total_systs[4],"QCD Scale #mu_{F} Up/Down")
else:
    leg2.AddEntry(mc_total_systs[0],"QGL Corr. Up/Down")
    #leg2.AddEntry(mc_total_systs[4],"QCD Scale #mu_{R} Up/Down")
    #leg2.AddEntry(mc_total_systs[6],"QCD Scale #mu_{F} Up/Down")
    leg2.AddEntry(mc_total_systs[2],"QCD Scale #mu_{F} Up/Down")
leg2.Draw("same")


canv.SetLogy(True)
#raw_input()
canv.SaveAs("%s_prefit.pdf" % channel)
canv.SaveAs("%s_prefit.png" % channel)
canv.SaveAs("%s_prefit.C" % channel)
ifile.Close()


