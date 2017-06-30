import ROOT
import sys

ROOT.gStyle.SetOptStat(0)

#SF_TT = 0.99
#SF_Wj0b = 1.17
#SF_Wj1b = 1.65
#SF_Wj2b = 1.65
#SF_Wj1b = 2.38
#SF_Wj2b = 1.07

SF_TT = 1.0
SF_Wj0b = 1.0
SF_Wj1b = 1.0
SF_Wj2b = 1.0

ifile = ROOT.TFile(sys.argv[1])
channel = sys.argv[2]
datafile = ROOT.TFile(sys.argv[3])
#data = ifile.Get("shapes_fit_s/%s/data_obs" % channel)
tt = ifile.Get("shapes_fit_s/%s/TT" % channel)
stop = ifile.Get("shapes_fit_s/%s/s_Top" % channel)
wj0b = ifile.Get("shapes_fit_s/%s/Wj0b" % channel)
wj1b = ifile.Get("shapes_fit_s/%s/Wj1b" % channel)
wj2b = ifile.Get("shapes_fit_s/%s/Wj2b" % channel)
zj0b = ifile.Get("shapes_fit_s/%s/Zj0b" % channel)
zj1b = ifile.Get("shapes_fit_s/%s/Zj1b" % channel)
zj2b = ifile.Get("shapes_fit_s/%s/Zj2b" % channel)
wh = ifile.Get("shapes_fit_s/%s/WH_hbb" % channel)
vvhf = ifile.Get("shapes_fit_s/%s/VVHF" % channel)
vvlf = ifile.Get("shapes_fit_s/%s/VVLF" % channel)
#qcd = ifile.Get("shapes_fit_s/%s/QCD" % channel)
data_hist = datafile.Get("BDT_%s_data_obs" % channel)

doBlind = True
if (doBlind):
    for i in range(data_hist.GetNbinsX()-4,data_hist.GetNbinsX()+1):
        data_hist.SetBinContent(i,0.0)
        data_hist.SetBinError(i,0.0)

#data = ROOT.TGraphErrors()
#for ibin in range(1,data_hist.GetNbinsX()+1):
#   data.SetPoint(ibin, ibin-0.5, data_hist.GetBinContent(ibin))
#   data.SetPointError(ibin,0.5, data_hist.GetBinError(ibin))

nbins=tt.GetNbinsX()
data = ROOT.TH1F("data_hist","data_hist",nbins,tt.GetBinLowEdge(1),tt.GetBinCenter(nbins)+tt.GetXaxis().GetBinWidth(nbins)/2)
for ibin in range(1,data_hist.GetNbinsX()+1):
   data.SetBinContent(ibin, data_hist.GetBinContent(ibin))
   data.SetBinError(ibin,data_hist.GetBinError(ibin))

data.SetMarkerStyle(20)
stop.SetLineColor(432)
stop.SetFillColor(432)
tt.SetLineColor(600)
tt.SetFillColor(600)
wj0b.SetLineColor(418)
wj0b.SetFillColor(418)
wj1b.SetLineColor(410)
wj1b.SetFillColor(410)
wj2b.SetLineColor(416)
wj2b.SetFillColor(416)
zj0b.SetLineColor(394)
zj0b.SetFillColor(394)
zj1b.SetLineColor(398)
zj1b.SetFillColor(398)
zj2b.SetLineColor(400)
zj2b.SetFillColor(400)
wh.SetLineColor(808)
wh.SetFillColor(808)
#qcd.SetLineColor(616)
#qcd.SetFillColor(616)
#vvhf.SetLineColor(922)
#vvhf.SetFillColor(922)
vvhf.SetLineColor(632)
#vvhf.SetFillColor(632)
vvlf.SetLineColor(921)
vvlf.SetFillColor(921)

vvhf.SetLineWidth(3)
#vvlf.SetLineWidth(3)
#wh.SetLineWidth(3)

#tt.Scale(SF_TT)
#wj0b.Scale(SF_Wj0b)
#wj1b.Scale(SF_Wj1b)
#wj2b.Scale(SF_Wj2b)
#wh.Scale(30.0)

stack = ROOT.THStack("stack","stack")
stack.Add(stop)
stack.Add(zj0b)
stack.Add(zj1b)
stack.Add(zj2b)
stack.Add(wj0b)
stack.Add(wj1b)
stack.Add(wj2b)
stack.Add(tt)
stack.Add(wh)
#stack.Add(vvlf)
#stack.Add(vvhf)
#stack.Add(qcd)

#vvlf.Scale(10.)
#vvhf.Scale(10.)
stackvv = ROOT.THStack("stackvv","stackvv")
stackvv.Add(vvlf)
stackvv.Add(vvhf)
#bkg.SetMaximum(1.3*max(bkg.GetMaximum(),wh.GetMaximum()))
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
xmax = 15
frame = ROOT.TH1F("frame","",1,xmin,xmax)
frame.SetStats(0)
frame.GetXaxis().SetLabelSize(0)
frame.GetXaxis().SetTitleOffset(0.91);
frame.GetYaxis().SetTitle("Events")
#frame.GetXaxis().SetTitle("atanhBDT")
frame.GetYaxis().SetLabelSize(0.04)
#frame.GetYaxis().SetRangeUser(10E-03,stack.GetStack().Last().GetMaximum()*1.3)
frame.GetYaxis().SetRangeUser(10E-03,stack.GetStack().Last().GetMaximum()*100)
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
vvhf.Draw("hist same")
vvhf.Draw("ep same")

leg = ROOT.TLegend(0.7,0.7,0.9,0.9)
#leg = ROOT.TLegend(0.4,0.6,0.7,0.9)
#leg = ROOT.TLegend(0.7,0.6,0.9,0.9)
leg.AddEntry(stop,"Single Top")
leg.AddEntry(wj0b,"Wj0b")
leg.AddEntry(wj1b,"Wj1b")
leg.AddEntry(wj2b,"Wj2b")
leg.AddEntry(tt,"TT")
leg.AddEntry(vvlf,"VVLF")
leg.AddEntry(vvhf,"VZ(bb)")
#leg.AddEntry(vvlf,"VVLF X 10")
#leg.AddEntry(vvhf,"VVHF X 10")
leg.AddEntry(wh,"WH")
leg.AddEntry(zj0b,"Zj0b")
leg.AddEntry(zj1b,"Zj1b")
leg.AddEntry(zj2b,"Zj2b")
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
frame2.SetMinimum(-0.5)	
frame2.SetMaximum(0.5) 
frame2.GetYaxis().SetLabelSize(0.02)
frame2.GetXaxis().SetLabelSize(0.04)
frame2.GetYaxis().SetTitleSize(0.04)
frame2.GetXaxis().SetTitle("BDT Bin #")
#frame2.GetXaxis().SetTitle("Sub-leading Jet CMVA Bin #")
frame2.SetStats(0)
frame2.GetYaxis().SetTitle("Data/MC - 1")	
frame2.Draw()

data2=data.Clone("data2")
mc_total=stack.GetStack().Last().Clone("mc_total")
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

mc_total_uncUp.Draw("HIST SAME")
mc_total_uncDown.Draw("HIST SAME")


line = ROOT.TLine(xmin,0,xmax,0)
line.SetLineStyle(3)
line.Draw("same")

data2.Draw("PEsame")
#print "Data yield: ",data.Integral()
#print "MC yield: ",stack.GetHistogram().Integral()

canv.SetLogy(True)
#raw_input()
canv.SaveAs("%s_postfit_vv.pdf" % channel)
canv.SaveAs("%s_postfit_vv.png" % channel)
canv.SaveAs("%s_postfit_vv.C" % channel)
ifile.Close()


