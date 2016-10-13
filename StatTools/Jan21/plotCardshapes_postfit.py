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
data = ifile.Get("shapes_fit_b/%s/data_obs" % channel)
tt = ifile.Get("shapes_fit_b/%s/TT" % channel)
stop = ifile.Get("shapes_fit_b/%s/s_Top" % channel)
wj0b = ifile.Get("shapes_fit_b/%s/Wj0b" % channel)
wj1b = ifile.Get("shapes_fit_b/%s/Wj1b" % channel)
wj2b = ifile.Get("shapes_fit_b/%s/Wj2b" % channel)
wh = ifile.Get("shapes_fit_b/%s/WH" % channel)
#qcd = ifile.Get("shapes_fit_b/%s/QCD" % channel)
data_hist = datafile.Get("BDT_%s_data_obs" % channel)

data = ROOT.TGraphErrors()
for ibin in range(1,data_hist.GetNbinsX()+1):
   data.SetPoint(ibin, ibin-0.5, data_hist.GetBinContent(ibin))
   data.SetPointError(ibin,0.5, data_hist.GetBinError(ibin))

stop.SetLineColor(432)
tt.SetLineColor(600)
wj0b.SetLineColor(418)
wj1b.SetLineColor(410)
wj2b.SetLineColor(416)
wh.SetLineColor(808)
#qcd.SetLineColor(616)

#tt.Scale(SF_TT)
#wj0b.Scale(SF_Wj0b)
#wj1b.Scale(SF_Wj1b)
#wj2b.Scale(SF_Wj2b)
wh.Scale(10.0)

stack = ROOT.THStack("stack","stack")
stack.Add(stop)
stack.Add(wj0b)
stack.Add(wj1b)
stack.Add(wj2b)
stack.Add(tt)
#stack.Add(qcd)

#bkg.SetMaximum(1.3*max(bkg.GetMaximum(),wh.GetMaximum()))
#bkg.SetMinimum(0.)
stack.SetMaximum(1.3*stack.GetMaximum())
stack.SetMinimum(0.)

#bkg.SetTitle(channel)
#bkg.GetXaxis().SetTitle("Mjj [GeV]")

canv = ROOT.TCanvas("canv","canv")
stack.Draw("hist")
#print "Chi2/NDOF = ",data.Chi2Test( stack.GetHistogram() , "UWCHI2/NDF")
#print "Total bkg: ",stack.GetHistogram().Integral()
#stack.Draw("ep same")
stack.SetTitle(channel)
#stack.GetXaxis().SetTitle("Mjj [GeV]")
data.Draw("ep same")
#print "Total data: ",data.Integral()
#bkg.Draw("hist")
#bkg.Draw("ep same")
wh.Draw("hist same")
wh.Draw("ep same")

#stack.GetXaxis().SetRangeUser(0,4)

leg = ROOT.TLegend(0.55,0.55,0.9,0.9)
leg.AddEntry(stop,"Single Top")
leg.AddEntry(wj0b,"Wj0b")
leg.AddEntry(wj1b,"Wj1b")
leg.AddEntry(wj2b,"Wj2b")
leg.AddEntry(tt,"TT")
leg.AddEntry(wh,"WH X 10")
leg.AddEntry(data,"Data","epl")
#leg.AddEntry(qcd,"QCD")
leg.Draw("same")
canv.SaveAs("%s_postfit.pdf" % channel)
canv.SaveAs("%s_postfit.png" % channel)

#print "Data yield: ",data.Integral()
#print "MC yield: ",stack.GetHistogram().Integral()

#raw_input()
ifile.Close()


