from ROOT import *
import sys

fitfile = TFile(sys.argv[1])

cov = fitfile.Get("covariance_fit_b")

xaxis = cov.GetXaxis()
yaxis = cov.GetYaxis()

#xbins = [179,180,181,182]
xbins = [368,369,370,371,372]
ybins = [3,4,5,6,7]

for xbin in xbins: 
    xlabel = xaxis.GetBinLabel(xbin)
    if (xlabel.find("SF")!=-1):
        xlabel = xlabel[xlabel.find("SF_")+3:]
        xlabel = xlabel[:xlabel.find("_")]
    for ybin in ybins:
        ylabel = yaxis.GetBinLabel(ybin)
        if (ylabel.find("SF")!=-1):
            ylabel = ylabel[ylabel.find("SF_")+3:]
            ylabel = ylabel[:ylabel.find("_")]
        print "%s vs. %s: %f +/- %f" % (xlabel, ylabel, cov.GetBinContent(xbin,ybin),cov.GetBinError(xbin,ybin))
