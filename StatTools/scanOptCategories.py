import ROOT
from math import sqrt,pow
import numpy
import sys
#import subprocess
import os

#ifile = ROOT.TFile("TMVA_13TeV_Dec4_H125Sig_0b1b2bWjetsTTbarBkg_Mjj.root")
#ifile = ROOT.TFile("/uscms_data/d3/sbc01/HbbAnalysis13TeV/PrincetonAnalysisTools/PlottingTools/Nm1Cuts/Dec4_McForOpt/output_allsamples_wBDTS.root")
if (len(sys.argv) != 4 and len(sys.argv) != 5):
    print "Usage: python scanOptCategories.py [ifilename] [nCat] [ofilename] [doWPs=False]"
    sys.exit(1)

ROOT.gROOT.SetBatch(True)

ifilename = sys.argv[1]
ifile = ROOT.TFile(ifilename)
tree = ifile.Get("tree")

useCombine = True # calculate significance using full data card with Combine
bdtname = "CMS_vhbb_BDT_Wln_13TeV"
#bdtname = "BDT_wMass_Dec14_3000_5"
presel = "H_mass>90 && H_mass<150"
#presel = "H_mass>0"

nCat = int(sys.argv[2])

ofile = ROOT.TFile(sys.argv[3],"RECREATE")
otree = ROOT.TTree("tree","tree")
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

# split different sig. eff. ranges by specified bin sizes of sig. eff.
# stepSize, sigEffMax
#binSplitting = [(1./50,0.0,0.1),(1./50,0.1,0.75),(1./100,0.75,0.95),(1./200,0.95,1.0)]
binSplitting = [(1./30,0.0,0.5),(1./60,0.5,1.0)]
#binSplitting = [(1./30,0.0,1.0)]
doWPs = False
if (len(sys.argv) > 4):
    if (sys.argv[4] ==  "True"): doWPs = True
    print sys.argv[4]
    print doWPs

nbins = 1000
if (doWPs):
    nbins = 14

xlow = -1.
xhigh = 1.
if (doWPs):
    xlow = 1
    xhigh = 15

hs = ROOT.TH1F("hs","hs",nbins,xlow,xhigh)
hbkg = ROOT.TH1F("hbkg","hbkg",nbins,xlow,xhigh)
hSig = ROOT.TH1F("hSig","hSig",nbins, xlow, xhigh)
hSig2 = ROOT.TH2F("hSig","hSig",nbins,xlow,xhigh,nbins,xlow,xhigh)

hbkgTot = ROOT.TH1F("hbkgTot","hbkgTot",100,-1.,1.)
hsTot = ROOT.TH1F("hsTot","hsTot",100, -1., 1.)
tree.Draw("%s>>hsTot" % bdtname,"(sampleIndex<0&&(%s))*weight" % presel)
tree.Draw("%s>>hbkgTot" % bdtname,"(sampleIndex>0&&(%s))*weight" % presel)

sDen = hsTot.Integral()
bkgDen = hbkgTot.Integral()

print "Total Denominator Yields: ",sDen,bkgDen
if doWPs:
    hs.Reset()
    hs.SetBinContent(1,sDen - 10.23)
    hs.SetBinContent(2,1.32)
    hs.SetBinContent(3,2.432)
    hs.SetBinContent(4,0.627)
    hs.SetBinContent(5,1.097)
    hs.SetBinContent(6,1.286)
    hs.SetBinContent(7,1.212)
    hs.SetBinContent(8,0.966)
    hs.SetBinContent(9,0.397)
    hs.SetBinContent(10,0.231)
    hs.SetBinContent(11,0.198)
    hs.SetBinContent(12,0.114)
    hs.SetBinContent(13,0.136)
    hs.SetBinContent(14,0.214)
    hs.SetBinContent(15,0.0)

    hbkg.Reset()
    hbkg.SetBinContent(1,bkgDen - 2494.9)
    hbkg.SetBinContent(2,523.6)
    hbkg.SetBinContent(3,1090.86)
    hbkg.SetBinContent(4,172.6)
    hbkg.SetBinContent(5,233.87)
    hbkg.SetBinContent(6,227.57)
    hbkg.SetBinContent(7,140.8)
    hbkg.SetBinContent(8,71.72)
    hbkg.SetBinContent(9,18.37)
    hbkg.SetBinContent(10,6.488)
    hbkg.SetBinContent(11,4.035)
    hbkg.SetBinContent(12,1.23)
    hbkg.SetBinContent(13,1.904)
    hbkg.SetBinContent(14,1.853)
    hbkg.SetBinContent(15,0.0)

else: 
    tree.Draw("%s>>hs" % bdtname,"(sampleIndex<0&&(%s))*weight" % presel)
    tree.Draw("%s>>hbkg" % bdtname,"(sampleIndex>0&&(%s))*weight" % presel)

boundaries= []
subBoundaries = {}
for stepSize, minSigEff, maxSigEff in binSplitting:
    subBoundaries[maxSigEff] = []

S_tot = hs.Integral()
boundaries.append(xlow)
# determine category boundaries based on incremental signal efficiency steps
for i in range(2,nbins+1):
    if (doWPs):
        boundaries.append(hs.GetBinLowEdge(i))
    else:
        lowBin = 1
        for stepSize, minSigEff, maxSigEff in binSplitting:
            S_frac = hs.Integral(lowBin,i-1)/S_tot;
            #print i, S_frac
            if (S_frac >= minSigEff and S_frac < maxSigEff and S_frac >= (minSigEff + ((len(subBoundaries[maxSigEff])+1)*stepSize) ) ):
                subBoundaries[maxSigEff].append(hs.GetBinLowEdge(i))
            
for stepSize, minSigEff, maxSigEff in binSplitting:
    boundaries.extend(subBoundaries[maxSigEff])
boundaries.append(xhigh)
print "boundary bin edges are..."
print boundaries
#print subBoundaries

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
    binLow = 1
    bound_label = ""
    for bound in bset:
        bound_label += "%f_" % bound
    bound_label = bound_label.replace('-','m')
    if (useCombine):
        if not os.path.exists("%i" % (nCat)):
            os.mkdir("%i" % (nCat))
        if not os.path.exists("%i/%s" % (nCat, bound_label)):
            os.mkdir("%i/%s" % (nCat, bound_label))
        print "cd %i/%s" % (nCat, bound_label)
        os.chdir("%i/%s" % (nCat, bound_label))
        os.system("pwd")
    for i in range(len(bset)-1):
        binLow = hs.GetXaxis().FindBin(bset[i])
        binHigh = hs.GetXaxis().FindBin(bset[i+1]) - 1
        #if (binHigh == 0): continue
        #print hs.GetBinLowEdge(binLow), hs.GetBinLowEdge(binHigh)
        #if (binLow == binHigh):
        #    S = 0.
        #    B = 0.
        #else:
        S = hs.Integral(binLow, binHigh)
        B = hbkg.Integral(binLow, binHigh)
        S_arr[i] = S
        B_arr[i] = B
        bound_arr[i] = bset[i]
        #print binLow, binHigh, S, B
        if (B > 0 or B<=0):
            combSens += pow(S,2)/B
            #combSens += pow(S,2)/pow(sqrt(B)+0.1*B,2)
            if (useCombine):
                print "python ../../../splitSamples.py ../../%s WmnCat%i 'isWmunu&&(%s)&&(%s>=%f)&&(%s<%f)' ../../../systematics_Wmn.txt" % (ifilename,i,presel,bdtname,bset[i],bdtname,bset[i+1])
		os.system("python ../../../splitSamples.py ../../%s WmnCat%i 'isWmunu&&(%s)&&(%s>=%f)&&(%s<%f)' ../../../systematics_Wmn.txt" % (ifilename,i,presel,bdtname,bset[i],bdtname,bset[i+1]))
		print "python ../../../splitSamples.py ../../%s WenCat%i ''isWenu&&(%s)&&(%s>=%f)&&(%s<%f)' ../../../systematics_Wen.txt" % (ifilename,i,presel,bdtname,bset[i],bdtname,bset[i+1])
                os.system("python ../../../splitSamples.py ../../%s WenCat%i 'isWenu&&(%s)&&(%s>=%f)&&(%s<%f)' ../../../systematics_Wen.txt" % (ifilename,i,presel,bdtname,bset[i],bdtname,bset[i+1]))
        #binLow = binHigh + 1
    combSens = sqrt(combSens)
    if (useCombine):
        cat_labels_Wmn = ""
        cat_labels_Wen = ""
        for i in range(nCat):
            cat_labels_Wmn += "WmnCat%i," % i
            cat_labels_Wen += "WenCat%i," % i
        bound_label = ""
        print "python ../../../printYields.py dc_Wmn.txt %s ../../../systematics_Wmn.txt" % cat_labels_Wmn
        os.system("python ../../../printYields.py dc_Wmn.txt %s ../../../systematics_Wmn.txt" % cat_labels_Wmn)
        print "python ../../../printYields.py dc_Wen.txt %s ../../../systematics_Wen.txt" % cat_labels_Wen
        os.system("python ../../../printYields.py dc_Wen.txt %s ../../../systematics_Wen.txt" % cat_labels_Wen)
        print "combineCards.py Wmn=dc_Wmn.txt Wen=dc_Wen.txt > dc.txt"  
        os.system("combineCards.py Wmn=dc_Wmn.txt Wen=dc_Wen.txt > dc.txt") 
        print "combine -M ProfileLikelihood --significance dc.txt -t -1 --expectSignal=1 > combine_output.txt"
        os.system("combine -M ProfileLikelihood --significance dc.txt -t -1 --expectSignal=1 > combine_output.txt")
        combOut = open("combine_output.txt","r")
        for line in combOut:
            if (line.find("Significance") != -1):
                combSens = float(line.split()[1])
                print combSens
        print "cd -"
        os.chdir("../..")
       
    Sig[0] = float(combSens)
    bound_arr[len(bset)-1] = bset[len(bset)-1]
    #print ncat[0],S_arr,B_arr,bound_arr,Sig[0]
    otree.Fill()
    if (nCat == 2):
        #print hs.GetXaxis().FindBin(bset[1]), combSens
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
        hSig.SetTitle("")
        if (doWPs):
            hSig.GetXaxis().SetTitle("WP Split Point")
        else: hSig.GetXaxis().SetTitle("BDT Score Split Point")
        hSig.GetYaxis().SetTitle("Combined Significance")
        hSig.Draw("hist")
    if (nCat == 3):
        hSig2.SetTitle("")
        if (doWPs):
            hSig2.GetXaxis().SetTitle("WP Split Point 1")
            hSig2.GetXaxis().SetTitle("WP Split Point 2")
        else: 
            hSig2.GetXaxis().SetTitle("BDT Score Split Point 1")
            hSig2.GetXaxis().SetTitle("BDT Score Split Point 2")
        hSig2.Draw("colZ")
    raw_input()
ofile.cd()
otree.Write()
ofile.Close()

