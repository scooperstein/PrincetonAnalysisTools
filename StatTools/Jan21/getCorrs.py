import ROOT
import sys

#_input = 'mlfit.root' 
_input = sys.argv[1]

ROOT.gStyle.SetOptStat(0)

f = ROOT.TFile(_input, 'r')
r = ROOT.gDirectory.Get('fit_s')
#r.floatParsFinal().Print("s")

#For corr

#parf = r.Print()
corr = r.correlationHist()

nX = corr.GetNbinsX()
nY = corr.GetNbinsY()

for i in range(1,nX+1):
    xlabel = corr.GetXaxis().GetBinLabel(i)
    if xlabel.find("SF_") == -1: continue
    for j in range(1,nY+1):
        #if j <= i: continue
        ylabel = corr.GetYaxis().GetBinLabel(j)
        val = corr.GetBinContent(i,j)
        if abs(val) < 0.05: continue
        print "%s vs. %s: %f" % (xlabel,ylabel,val)
    


#c = ROOT.TCanvas('c','c')
#c.cd()
#corr.Draw("textcolz")
#corr.Draw("colz")
#c.SaveAs("corr.pdf")
#c.SaveAs("corr.png")
#c.SaveAs("corr.root")
