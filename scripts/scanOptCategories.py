import ROOT
from math import sqrt,pow
import numpy
import sys

#ifile = ROOT.TFile("TMVA_13TeV_Dec4_H125Sig_0b1b2bWjetsTTbarBkg_Mjj.root")
#ifile = ROOT.TFile("/uscms_data/d3/sbc01/HbbAnalysis13TeV/PrincetonAnalysisTools/PlottingTools/Nm1Cuts/Dec4_McForOpt/output_allsamples_wBDTS.root")
ifile = ROOT.TFile(sys.argv[1])
tree = ifile.Get("tree")

nCat = int(sys.argv[2])

ofile = ROOT.TFile("ofile.root","RECREATE")
otree = ROOT.TTree("otree","otree")
ncat = numpy.zeros(1, dtype=int)
ncat[0] = nCat
S_arr = numpy.zeros(nCat, dtype=float)
B_arr = numpy.zeros(nCat, dtype=float)
Sig   = numpy.zeros(1, dtype=float)
bound_arr = numpy.zeros(nCat+1, dtype=float)
otree.Branch("nCat",ncat,"nCat/I")
otree.Branch("S",S_arr,"S[nCat]/D")
otree.Branch("B",B_arr,"B[nCat]/D")
otree.Branch("Sig",Sig,"Sig/D")
otree.Branch("Boundaries",bound_arr,"Boundaries[%i]/D" % (nCat + 1))

nBinEdges = 30
doWPs = False
if (len(sys.argv) > 3):
    if (sys.argv[3] ==  "True"): doWPs = True
    print sys.argv[3]
    print doWPs

nbins = 1000
if (doWPs):
    nbins = 12

xlow = -1.
xhigh = 1.
if (doWPs):
    xlow = 1
    xhigh = 13

hs = ROOT.TH1F("hs","hs",nbins,xlow,xhigh)
hbkg = ROOT.TH1F("hbkg","hbkg",nbins,xlow,xhigh)
hSig = ROOT.TH1F("hSig","hSig",nbins, xlow, xhigh)
hSig2 = ROOT.TH2F("hSig","hSig",nbins,xlow,xhigh,nbins,xlow,xhigh)

if not doWPs:
    tree.Draw("BDT_wMass_Dec4>>hs","(sampleIndex<0)*weight")
    tree.Draw("BDT_wMass_Dec4>>hbkg","(sampleIndex>0)*weight")

else:
    hs_wp.SetBinContent(1,1.32)
    hs_wp.SetBinContent(2,2.432)
    hs_wp.SetBinContent(3,0.627)
    hs_wp.SetBinContent(4,1.097)
    hs_wp.SetBinContent(5,1.286)
    hs_wp.SetBinContent(6,1.212)
    hs_wp.SetBinContent(7,0.966)
    hs_wp.SetBinContent(8,0.397)
    hs_wp.SetBinContent(9,0.231)
    hs_wp.SetBinContent(10,0.198)
    hs_wp.SetBinContent(11,0.114)
    hs_wp.SetBinContent(12,0.136)
    hs_wp.SetBinContent(13,0.214)

    hbkg_wp.SetBinContent(1,523.6)
    hbkg_wp.SetBinContent(2,1090.86)
    hbkg_wp.SetBinContent(3,172.6)
    hbkg_wp.SetBinContent(4,233.87)
    hbkg_wp.SetBinContent(5,227.57)
    hbkg_wp.SetBinContent(6,140.8)
    hbkg_wp.SetBinContent(7,71.72)
    hbkg_wp.SetBinContent(8,18.37)
    hbkg_wp.SetBinContent(9,6.488)
    hbkg_wp.SetBinContent(10,4.035)
    hbkg_wp.SetBinContent(11,1.23)
    hbkg_wp.SetBinContent(12,1.904)
    hbkg_wp.SetBinContent(13,1.853)

boundaries = []
S_tot = hs.Integral()
# determine category boundaries based on incremental signal efficiency steps
for i in range(1,nbins+1):
    S_frac = hs.Integral(1,i-1)/S_tot;
    #print i, S_frac
    if (doWPs):
        boundaries.append(hs.GetBinLowEdge(i))
    elif (S_frac >= ((len(boundaries)+1)*1.0)/nBinEdges):
        boundaries.append(hs.GetBinLowEdge(i))
boundaries.append(xhigh)
print "boundary bin edges are..."
print boundaries

boundaries_array = numpy.zeros(len(boundaries),dtype=float)
for i in range(len(boundaries)):
    boundaries_array[i] = boundaries[i]
if not doWPs:
    hs = hs.Rebin(len(boundaries)-1,"",boundaries_array)
    hbkg = hbkg.Rebin(len(boundaries)-1,"",boundaries_array)
    hSig = hSig.Rebin(len(boundaries)-1,"",boundaries_array)
nbinsNew = hs.GetNbinsX()
print "nbinsNew = ",nbinsNew

# list of tuples of all possible category boundaries
bsets = []
if (nCat == 2):
    for bound in boundaries:
        bsets.append([boundaries[0],bound,boundaries[len(boundaries)-1]])
elif (nCat == 3):
    for bound1 in boundaries:
        for bound2 in boundaries:
            if (bound2 <= bound1): continue
            bsets.append([boundaries[0],bound1,bound2,boundaries[len(boundaries)-1]])
elif (nCat == 4):
    for bound1 in boundaries:
        for bound2 in boundaries:
            if (bound2 <= bound1): continue
            for bound3 in boundaries:
                if (bound3 <= bound2): continue
                bsets.append([boundaries[0],bound1,bound2,bound3,boundaries[len(boundaries)-1]])

elif (nCat == 5):
    for bound1 in boundaries:
        for bound2 in boundaries:
            if (bound2 <= bound1): continue
            for bound3 in boundaries:
                if (bound3 <= bound2): continue
                for bound4 in boundaries:
                    if (bound4 <= bound3): continue
                    bsets.append([boundaries[0],bound1,bound2,bound3,bound4,boundaries[len(boundaries)-1]])
elif (nCat == 6):
    for bound1 in boundaries:
        for bound2 in boundaries:
            if (bound2 <= bound1): continue
            for bound3 in boundaries:
                if (bound3 <= bound2): continue
                for bound4 in boundaries:
                    if (bound4 <= bound3): continue
                    for bound5 in boundaries:
                        if (bound5 <= bound4): continue
                        bsets.append([boundaries[0],bound1,bound2,bound3,bound4,bound5,boundaries[len(boundaries)-1]])

elif (nCat == 7):
    for bound1 in boundaries:
        for bound2 in boundaries:
            if (bound2 <= bound1): continue
            for bound3 in boundaries:
                if (bound3 <= bound2): continue
                for bound4 in boundaries:
                    if (bound4 <= bound3): continue
                    for bound5 in boundaries:
                        if (bound5 <= bound4): continue
                        for bound6 in boundaries:
                            if (bound6 <= bound5): continue        
                            bsets.append([boundaries[0],bound1,bound2,bound3,bound4,bound5,bound6,boundaries[len(boundaries)-1]])

elif (nCat == 8):
    for bound1 in boundaries:
        for bound2 in boundaries:
            if (bound2 <= bound1): continue
            for bound3 in boundaries:
                if (bound3 <= bound2): continue
                for bound4 in boundaries:
                    if (bound4 <= bound3): continue
                    for bound5 in boundaries:
                        if (bound5 <= bound4): continue
                        for bound6 in boundaries:
                            if (bound6 <= bound5): continue   
                            for bound7 in boundaries:
                                if (bound7 <= bound6): continue     
                                bsets.append([boundaries[0],bound1,bound2,bound3,bound4,bound5,bound6,bound7,boundaries[len(boundaries)-1]])
bestSig = 0.
bestWP = ()

#print bsets
for bset in bsets:
    combSens = 0.
    for i in range(len(bset)-1):
        binLow = hs.GetXaxis().FindBin(bset[i])
        binHigh = hs.GetXaxis().FindBin(bset[i+1])
        S = hs.Integral(binLow, binHigh)
        B = hbkg.Integral(binLow, binHigh)
        S_arr[i] = S
        B_arr[i] = B
        bound_arr[i] = bset[i]
        #print binLow, binHigh, S, B
        if (B > 0):
            combSens += pow(S,2)/B
            #combSens += pow(S,2)/pow(sqrt(B)+0.1*B,2)
    combSens = sqrt(combSens)
    Sig[0] = float(combSens)
    bound_arr[len(bset)-1] = bset[len(bset)-1]
    #print ncat[0],S_arr,B_arr,bound_arr,Sig[0]
    otree.Fill()
    if (nCat == 2):
        print hs.GetXaxis().FindBin(bset[1]), combSens
        hSig.Fill(bset[1], combSens)
    elif (nCat == 3):
        hSig2.Fill(bset[1], bset[2], combSens)
    bset.append(combSens)
    if (combSens > bestSig):
        bestSig = combSens
        bestWP = bset 

#for bset in bsets:
#    print bset

#for i in range(1,nbinsNew+1):
#   sensHigh = 0.
#   S_high = hs.Integral(i, nbinsNew)
#   B_high = hbkg.Integral(i, nbinsNew) 
#   if (B_high > 0):
#       sensHigh = S_high/sqrt(B_high)
#
#   sensLow = 0.
#   S_low = 0.
#   B_low = 0.
#   if (i != 1):
#      S_low = hs.Integral(1, i-1)
#      B_low = hbkg.Integral(1, i-1)
#   if (B_low > 0):
#       sensLow = S_low/sqrt(B_low)
#
#   hSig.SetBinContent(i, sqrt(pow(sensLow, 2) + pow(sensHigh, 2)) )
#   print S_low,B_low,S_high,B_high
#   print "%i: %f: %f: %f" % (i, sensLow, sensHigh, sqrt(pow(sensLow, 2) + pow(sensHigh, 2)))

print "best category boundaries are..."
print bestWP
print "with combined sensitivity: %f" % bestSig

if (nCat == 2 or nCat == 3):
    c1 = ROOT.TCanvas("c1","c1")
    if (nCat == 2):
        hSig.Draw("hist")
    if (nCat == 3):
        hSig2.Draw("colZ")
    raw_input()
ofile.cd()
otree.Write()
ofile.Close()

