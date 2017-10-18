import ROOT
import sys

idir = sys.argv[1]

from os import listdir
from os.path import isfile, join

onlyfiles = [ f for f in listdir(idir) if isfile(join(idir,f)) ]
for rootfile in onlyfiles:
    if rootfile.find(".root") != -1 and rootfile.find("weighted")==-1:
        sample = rootfile[rootfile.find("_")+1:rootfile.find(".root")]
        rootfile = idir + rootfile
        ifile = ROOT.TFile(rootfile,"r")
        tree = ifile.Get("tree")
        histw = ROOT.TH1F("histw","histw",56,200,3000)
        histuw = ROOT.TH1F("histuw","histuw",56,200,3000)
        tree.Draw("Mjj>>histw","((Vtype==2||Vtype==3)&&selLeptons_relIso_0<0.06&&Pass_nominal)*weight") 
        tree.Draw("Mjj>>histuw","((Vtype==2||Vtype==3)&&selLeptons_relIso_0<0.06&&Pass_nominal)*weight*(1./corrQGL)")
        fac = 0.
        if (histw.Integral() > 0):
            fac = histuw.Integral() / histw.Integral() 
   
        s = 'sampleMap["%s"] = %f' % (sample, fac)
        print s
        ifile.Close()
