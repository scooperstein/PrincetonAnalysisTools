import sys
import ROOT
import numpy

wpfile = ROOT.TFile("Wp_nloEWK_weight_unnormalized.root")
hist_wp = wpfile.Get("SignalWeight_nloEWK")
wmfile = ROOT.TFile("Wm_nloEWK_weight_unnormalized.root")
hist_wm = wmfile.Get("SignalWeight_nloEWK")
zllfile = ROOT.TFile("Zll_nloEWK_weight_unnormalized.root")
hist_zll = zllfile.Get("SignalWeight_nloEWK")

hist_wp.Rebin(4)
hist_wp.Scale(1./4.)
hist_wm.Rebin(4)
hist_wm.Scale(1./4.)
hist_zll.Rebin(4)
hist_zll.Scale(1./4.)

def getVHCorrFactor(V_pt, typ):
    if typ == 0: hist = hist_wp
    if typ == 1: hist = hist_wm
    if typ == 2: hist = hist_zll

    ibin = hist.GetXaxis().FindBin(V_pt)
    if (ibin < 1): return hist.GetBinContent(1)
    if (ibin > hist.GetNbinsX()): return hist.GetBinContent(hist.GetNbinsX())
    return hist.GetBinContent(ibin)

ifile = ROOT.TFile.Open(sys.argv[1])
ofile = ROOT.TFile(sys.argv[2],"RECREATE")

tree = ifile.Get("tree")
nentries = tree.GetEntries()
#nentries = 100

otree = tree.CloneTree(0)
VHCorrFactor = numpy.zeros(1,dtype=float)
otree.Branch("VHCorrFactor",VHCorrFactor,"VHCorrFactor/D")

print "Processing %i entries" % nentries
for i in range(nentries):
    if (i%10000 == 0): print "Processing entry: %i" % i
    tree.GetEntry(i)
    si = tree.sampleIndex
    #V_pt = tree.V_pt
    if tree.nGenVbosons!= 1: 
        VHCorrFactor[0] = 1.0
    elif si == -12500:
        V_pt = tree.GenVbosons_pt[0]
        VHCorrFactor[0] = getVHCorrFactor(V_pt,0) + (0.172020708/0.159 - 1)
    elif si == -12501:
        V_pt = tree.GenVbosons_pt[0]
        VHCorrFactor[0] = getVHCorrFactor(V_pt,1) +  (0.108992575/0.101 - 1)
    elif si == -12502:
        V_pt = tree.GenVbosons_pt[0]
        VHCorrFactor[0] = getVHCorrFactor(V_pt,2) + (0.05438/0.053 - 1)
    else:
        VHCorrFactor[0] = 1.0
    otree.Fill()
    
ofile.cd()
otree.Write()

for key in ifile.GetListOfKeys():
     if key.GetName() == 'tree':
         continue
     obj = key.ReadObj()
     obj.Write(key.GetName())
