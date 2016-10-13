import ROOT
from math import sqrt,pow
import numpy
import sys
#import subprocess
import os

#ifile = ROOT.TFile("TMVA_13TeV_Dec4_H125Sig_0b1b2bWjetsTTbarBkg_Mjj.root")
#ifile = ROOT.TFile("/uscms_data/d3/sbc01/HbbAnalysis13TeV/PrincetonAnalysisTools/PlottingTools/Nm1Cuts/Dec4_McForOpt/output_allsamples_wBDTS.root")
if (len(sys.argv) != 4 and len(sys.argv) != 5):
    print "Usage: python scanOptCategories.py [ifilename] [nCat] [ofilename] [procStep] [doWPs=False]"
    sys.exit(1)

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

ifilename = sys.argv[1]
ifile = ROOT.TFile(ifilename)
tree = ifile.Get("tree")

useCombine = True # calculate significance using full data card with Combine
#bdtname = "V_pt"
#bdtname = "CMS_vhbb_BDT_Wln_13TeV"
#bdtname = "BDT_May30_noMbb_v2"
bdtname = "V_pt"
#presel = ""
presel = "Pass_nominal==1&&V_pt>100&&(Vtype==2||Vtype==3)"
#presel = "H_mass>0"
xlow = 100
xhigh = 1000
#xlow = -1.
#xhigh = 1.

#bdtname_f = "H_mass"
#bdtname_f = "BDT_May30_350_5"
#bdtname_f = "CMS_vhbb_BDT_Wln_13TeV"
bdtname_f = "BDT_V24_Sep20_fullTTStats_500_4"
xlow_f = -1
xhigh_f = 1
#xlow_f = 50
#xhigh_f = 250

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

procStep = int(sys.argv[4])

# split different sig. eff. ranges by specified bin sizes of sig. eff.
# stepSize, sigEffMax
#binSplitting = [(1./50,0.0,0.1),(1./50,0.1,0.75),(1./100,0.75,0.95),(1./200,0.95,1.0)]
#binSplitting = [(1.,0.0,0.6),(1./20,0.6,0.8),(1./40,0.8,1.0)]
#binSplitting = [(1./4,0.0,0.5),(1./16,0.5,1.0)]
binSplitting = [(1./15,0.0,1.0)]
#binSplitting = [(1./50,0.0,1.0)]
doWPs = False
if (len(sys.argv) > 5):
    if (sys.argv[5] ==  "True"): doWPs = True
    print sys.argv[5]
    print doWPs

nbins = 1000
if (doWPs):
    nbins = 13
    bdtname = "WPInt"

if (doWPs):
    xlow = 1
    xhigh = 14

hs = ROOT.TH1F("hs","hs",nbins,xlow,xhigh)
hbkg = ROOT.TH1F("hbkg","hbkg",nbins,xlow,xhigh)
hSig = ROOT.TH1F("hSig","hSig",nbins, xlow, xhigh)
hSig2 = ROOT.TH2F("hSig","hSig",nbins,xlow,xhigh,nbins,xlow,xhigh)

#hbkgTot = ROOT.TH1F("hbkgTot","hbkgTot",100,xlow,xhigh)
#hsTot = ROOT.TH1F("hsTot","hsTot",100, xlow, xhigh)
#tree.Draw("%s>>hsTot" % bdtname,"(sampleIndex==-12501&&(%s))*weight*weight_PU*bTagWeight*CS_SF*weight_ptQCD*weight_ptEWK*selLeptons_SF_IsoTight[lepInd]*(selLeptons_SF_IdCutTight[lepInd]*isWmunu + selLeptons_SF_IdMVATight[lepInd]*isWenu)*weight_PU*bTagWeight*CS_SF*weight_ptQCD*weight_ptEWK*selLeptons_SF_IsoTight[lepInd]*(selLeptons_SF_IdCutTight[lepInd]*isWmunu + selLeptons_SF_IdMVATight[lepInd]*isWenu)" % presel)
#tree.Draw("%s>>hbkgTot" % bdtname,"(sampleIndex>0&&sampleIndex!=12&&(%s))*weight*weight_PU*bTagWeight*CS_SF*weight_ptQCD*weight_ptEWK*selLeptons_SF_IsoTight[lepInd]*(selLeptons_SF_IdCutTight[lepInd]*isWmunu + selLeptons_SF_IdMVATight[lepInd]*isWenu)*weight_PU*bTagWeight*CS_SF*weight_ptQCD*weight_ptEWK*selLeptons_SF_IsoTight[lepInd]*(selLeptons_SF_IdCutTight[lepInd]*isWmunu + selLeptons_SF_IdMVATight[lepInd]*isWenu)" % presel)
#
#sDen = hsTot.Integral()
#bkgDen = hbkgTot.Integral()
#
#print "Total Denominator Yields: ",sDen,bkgDen
#print "Overal S/sqrt(B): %f" % (sDen/sqrt(bkgDen))

tree.Draw("%s>>hs" % bdtname,"(sampleIndex==-12501&&(%s))*weight" % presel)
tree.Draw("%s>>hbkg" % bdtname,"(sampleIndex>0&&sampleIndex!=50&&sampleIndex!=51&&sampleIndex!=52&&(%s))*weight" % presel)
#tree.Draw("%s>>hs" % bdtname,"(sampleIndex==-12501&&(%s))*weight*weight_PU*bTagWeight*CS_SF*weight_ptQCD*weight_ptEWK*selLeptons_SF_IsoTight[lepInd]*(selLeptons_SF_IdCutTight[lepInd]*isWmunu + selLeptons_SF_IdMVATight[lepInd]*isWenu)" % presel)
#tree.Draw("%s>>hbkg" % bdtname,"(sampleIndex>0&&sampleIndex!=50&&sampleIndex!=51&&sampleIndex!=52&&(%s))*weight*weight_PU*bTagWeight*CS_SF*weight_ptQCD*weight_ptEWK*selLeptons_SF_IsoTight[lepInd]*(selLeptons_SF_IdCutTight[lepInd]*isWmunu + selLeptons_SF_IdMVATight[lepInd]*isWenu)" % presel)

boundaries= []
boundaryEff = {}
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
                boundaryEff[hs.GetBinLowEdge(i)] = S_frac
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
    overlappingCats = False # throw out divisions where a category has 0 
    for i in range(len(bset)-1):
        if (bset[i] == bset[i+1]):
             overlappingCats = True
    if (overlappingCats): continue
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
        #print "cd %i/%s" % (nCat, bound_label)
        os.chdir("%i/%s" % (nCat, bound_label))
        os.system("pwd")
    for i in range(len(bset)-1):
        if (doWPs and i==0): continue # for the moment let's have the lower edge be one of the cut-based WP's instead of presel
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
        if (B > 0):
            combSens += pow(S,2)/B
            #combSens += pow(S,2)/pow(sqrt(B)+0.1*B,2)
            if (useCombine):
                #cardnbin = 40
                #if (i == 2): cardnbin = 20
                #if (i == 3): cardnbin = 15
                cardnbin = 15
                cwd = os.getcwd()
                script_text = '''export ORIG_DIR=$PWD
cd %s
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
''' % cwd
                #script_text += "python ../../../splitSamples.py -i ../../%s -c WmnCat%i -p 'isWmunu&&(%s)&&(%s>=%f)&&(%s<%f)&&evt%%2==0' -s ../../../systematics_Wmn.txt -w '2.0,8.26' -t ../../%s -tt ../../%s -v %s -xl %f -xh %f -o $ORIG_DIR/hists_WmnCat%i.root -b $ORIG_DIR/binStats_WmnCat%i.txt -n %i" % (ifilename,i,presel,bdtname,bset[i],bdtname,bset[i+1],ifilename.replace("output_mc.root","output_data.root"),ifilename.replace("output_mc.root","output_ttpowheg.root"),bdtname_f,xlow_f,xhigh_f,i,i,cardnbin)
                script_text += "python ../../../splitSamples.py -i ../../%s -c WmnCat%i -p 'run<=276811&&isWmunu&&(%s)&&(%s>=%f)&&(%s<%f)&&(evt%%2==0||sampleIndex==0)&&Vtype==2&&H_mass>90&&H_mass<150&&selLeptons_relIso_0[lepInd]<0.06' -s ../../../systematics_Wmn_CS.txt -w '2.0' -t ../../%s -tt ../../%s -v %s -xl %f -xh %f -o $ORIG_DIR/hists_WmnCat%i.root -b $ORIG_DIR/binStats_WmnCat%i.txt -n %i -d 1" % (ifilename,i,presel,bdtname,bset[i],bdtname,bset[i+1],ifilename.replace("output_mc.root","output_data.root"),ifilename.replace("output_mc.root","output_ttpowheg.root"),bdtname_f,xlow_f,xhigh_f,i,i,cardnbin)
#                script_text += "python ../../../splitSamples.py -i ../../%s -c WmnCat%i -p 'isWmunu&&(%s)&&(%s>=%f)&&(%s<%f)&&(H_pt>=%f)&&(H_pt<%f)&&evt%%2==0' -s ../../../systematics_Wmn_CS.txt -w '2.0' -t ../../%s -tt ../../%s -v %s -xl %f -xh %f -o $ORIG_DIR/hists_WmnCat%i.root -b $ORIG_DIR/binStats_WmnCat%i.txt -n %i" % (ifilename,i,presel,bdtname,bset[i],bdtname,bset[i+1],bset[i],bset[i+1],ifilename.replace("output_mc.root","output_data.root"),ifilename.replace("output_mc.root","output_ttpowheg.root"),bdtname_f,xlow_f,xhigh_f,i,i,cardnbin)
                #script_text += "python ../../../splitSamples.py -i ../../%s -c WmnCat%i -p 'isWmunu&&(%s)&&(%s>=%f)&&(%s<%f)&&(H_pt>=%f)&&(H_pt<%f)&&evt%%2==0' -s ../../../systematics_Wmn.txt -w 'selLeptons_SF_IdCutTight[lepInd],selLeptons_SF_IsoTight[lepInd],weight_PU,weight_ptEWK,weight_ptQCD,CS_SF,2.0,8.26' -t ../../%s -tt ../../%s -v %s -xl %f -xh %f -o $ORIG_DIR/hists_WmnCat%i.root -b $ORIG_DIR/binStats_WmnCat%i.txt -n %i" % (ifilename,i,presel,bdtname,bset[i],bdtname,bset[i+1],bset[i],bset[i+1],ifilename.replace("output_mc.root","output_data.root"),ifilename.replace("output_mc.root","output_ttpowheg.root"),bdtname_f,xlow_f,xhigh_f,i,i,cardnbin)
                #script_text += "python ../../../splitSamples.py -i ../../%s -c WmnCat%i -p 'isWmunu&&(%s)&&(%s>=%f)&&(%s<%f)&&evt%%2==0' -s ../../../systematics_Wmn.txt -w 'selLeptons_SF_IdCutTight[lepInd],selLeptons_SF_IsoTight[lepInd],weight_PU,weight_ptEWK,weight_ptQCD,CS_SF,2.0,8.26' -t ../../%s -tt ../../%s -v %s -xl %f -xh %f -o $ORIG_DIR/hists_WmnCat%i.root -b $ORIG_DIR/binStats_WmnCat%i.txt -n %i" % (ifilename,i,presel,bdtname,bset[i],bdtname,bset[i+1],ifilename.replace("output_mc.root","output_data.root"),ifilename.replace("output_mc.root","output_ttpowheg.root"),bdtname_f,xlow_f,xhigh_f,i,i,cardnbin)
                script_text += "\n"
                #os.system("python ../../../splitSamples.py -i ../../%s -c WmnCat%i -p 'isWmunu&&(%s)&&(%s>=%f)&&(%s<%f)&&evt%%2==0' -s ../../../systematics_Wmn.txt -w 'selLeptons_SF_IdCutTight[lepInd],selLeptons_SF_IsoTight[lepInd],weight_PU,weight_ptEWK,weight_ptQCD,CS_SF,2.0,8.26' -t ../../%s -tt ../../%s -v %s -xl %f -xh %f-o $ORIG_DIR/hists_WmnCat%i.root" % (ifilename,i,presel,bdtname,bset[i],bdtname,bset[i+1],ifilename.replace("output_mc.root","output_data.root"),ifilename.replace("output_mc.root","output_ttpowheg.root"),bdtname_f,xlow_f,xhigh_f,i))
		#script_text +=  "python ../../../splitSamples.py -i ../../%s -c WenCat%i -p 'isWenu&&(%s)&&(%s>=%f)&&(%s<%f)&&evt%%2==0' -s ../../../systematics_Wen.txt -w '2.0,8.26' -t ../../%s -tt ../../%s -v %s -xl %f -xh %f -o $ORIG_DIR/hists_WenCat%i.root -b $ORIG_DIR/binStats_WenCat%i.txt -n %i" % (ifilename,i,presel,bdtname,bset[i],bdtname,bset[i+1],ifilename.replace("output_mc.root","output_data.root"),ifilename.replace("output_mc.root","output_ttpowheg.root"),bdtname_f,xlow_f,xhigh_f,i,i,cardnbin)
		script_text +=  "python ../../../splitSamples.py -i ../../%s -c WenCat%i -p 'run<=276811&&isWenu&&(%s)&&(%s>=%f)&&(%s<%f)&&(evt%%2==0||sampleIndex==0)&&Vtype==3&&H_mass>90&&H_mass<150' -s ../../../systematics_Wen_CS.txt -w '2.0' -t ../../%s -tt ../../%s -v %s -xl %f -xh %f -o $ORIG_DIR/hists_WenCat%i.root -b $ORIG_DIR/binStats_WenCat%i.txt -n %i -d 1" % (ifilename,i,presel,bdtname,bset[i],bdtname,bset[i+1],ifilename.replace("output_mc.root","output_data.root"),ifilename.replace("output_mc.root","output_ttpowheg.root"),bdtname_f,xlow_f,xhigh_f,i,i,cardnbin)
#		script_text +=  "python ../../../splitSamples.py -i ../../%s -c WenCat%i -p 'isWenu&&(%s)&&(%s>=%f)&&(%s<%f)&&(H_pt>=%f)&&(H_pt<%f)&&evt%%2==0' -s ../../../systematics_Wen_CS.txt -w '2.0' -t ../../%s -tt ../../%s -v %s -xl %f -xh %f -o $ORIG_DIR/hists_WenCat%i.root -b $ORIG_DIR/binStats_WenCat%i.txt -n %i" % (ifilename,i,presel,bdtname,bset[i],bdtname,bset[i+1],bset[i],bset[i+1],ifilename.replace("output_mc.root","output_data.root"),ifilename.replace("output_mc.root","output_ttpowheg.root"),bdtname_f,xlow_f,xhigh_f,i,i,cardnbin)
		#script_text +=  "python ../../../splitSamples.py -i ../../%s -c WenCat%i -p 'isWenu&&(%s)&&(%s>=%f)&&(%s<%f)&&(H_pt>=%f)&&(H_pt<%f)&&evt%%2==0' -s ../../../systematics_Wen.txt -w 'selLeptons_SF_IdMVATight[lepInd],weight_PU,weight_ptEWK,weight_ptQCD,CS_SF,2.0,8.26' -t ../../%s -tt ../../%s -v %s -xl %f -xh %f -o $ORIG_DIR/hists_WenCat%i.root -b $ORIG_DIR/binStats_WenCat%i.txt -n %i" % (ifilename,i,presel,bdtname,bset[i],bdtname,bset[i+1],bset[i],bset[i+1],ifilename.replace("output_mc.root","output_data.root"),ifilename.replace("output_mc.root","output_ttpowheg.root"),bdtname_f,xlow_f,xhigh_f,i,i,cardnbin)
		#script_text +=  "python ../../../splitSamples.py -i ../../%s -c WenCat%i -p 'isWenu&&(%s)&&(%s>=%f)&&(%s<%f)&&evt%%2==0' -s ../../../systematics_Wen.txt -w 'selLeptons_SF_IdMVATight[lepInd],weight_PU,weight_ptEWK,weight_ptQCD,CS_SF,2.0,8.26' -t ../../%s -tt ../../%s -v %s -xl %f -xh %f -o $ORIG_DIR/hists_WenCat%i.root -b $ORIG_DIR/binStats_WenCat%i.txt -n %i" % (ifilename,i,presel,bdtname,bset[i],bdtname,bset[i+1],ifilename.replace("output_mc.root","output_data.root"),ifilename.replace("output_mc.root","output_ttpowheg.root"),bdtname_f,xlow_f,xhigh_f,i,i,cardnbin)
                script_text += "\n"
                condor_script = open("condor_runscript_%i.sh" %i,"w")
                condor_script.write(script_text)
                condor_script.close()
                submit_text = '''universe = vanilla
Executable = condor_runscript_%i.sh
Should_Transfer_Files = YES
Output = stdout_%i
Error  = stderr_%i
Log    = log_%i
Notification = never
WhenToTransferOutput=On_Exit
Queue  1
                ''' % (i,i,i,i) 
                submitfile = open("submit_%i.txt" % i,"w")
                submitfile.write(submit_text)
                submitfile.close()
                if (procStep == 0):
                    os.system("condor_submit submit_%i.txt" % i)
    combSens = sqrt(combSens)
    if (useCombine):
        cat_labels_Wmn = ""
        cat_labels_Wen = ""
        for i in range(nCat):
            cat_labels_Wmn += "WmnCat%i," % i
            cat_labels_Wen += "WenCat%i," % i
        cat_labels_Wmn = cat_labels_Wmn[:-1]
        cat_labels_Wen = cat_labels_Wen[:-1]
        if (doWPs):
            cat_labels_Wmn = cat_labels_Wmn[cat_labels_Wmn.find("WmnCat0")+8:]
            cat_labels_Wen = cat_labels_Wen[cat_labels_Wen.find("WenCat0")+8:]
        bound_label = ""
        script_text_cards = '''export ORIG_DIR=$PWD
cd /uscms_data/d3/sbc01/HbbAnalysis13TeV/CMSSW_7_1_5/src
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
cd %s
''' % os.getcwd()
        script_text_cards += "\n"
        combineCards_text = ""
        for i in range(len(cat_labels_Wmn.split(','))):
            script_text_cards += "python ../../../printYields.py -o $ORIG_DIR/dc_WmnCat%i.txt -c WmnCat%i -s ../../../systematics_Wmn_CS.txt -r TT_Wln,Wj0b_Wln,Wj1b_Wln,Wj2b_Wln\n" % (i,i)
            script_text_cards += "python ../../../printYields.py -o $ORIG_DIR/dc_WenCat%i.txt -c WenCat%i -s ../../../systematics_Wen_CS.txt -r TT_Wln,Wj0b_Wln,Wj1b_Wln,Wj2b_Wln\n" % (i,i)
            combineCards_text += " WmnCat%i=dc_WmnCat%i.txt" % (i,i)
            combineCards_text += " WenCat%i=dc_WenCat%i.txt" % (i,i)
        script_text_cards += "cp hists_*.root $ORIG_DIR\n"
        script_text_cards += "cd $ORIG_DIR\n"
        script_text_cards += "cp /uscms_data/d3/sbc01/HbbAnalysis13TeV/CMSSW_7_6_3_patch2/src/PrincetonAnalysisTools/StatTools/Jan21/V24_CR_12p9fb/vhbb*.txt $ORIG_DIR\n"
        script_text_cards += "cp /uscms_data/d3/sbc01/HbbAnalysis13TeV/CMSSW_7_6_3_patch2/src/PrincetonAnalysisTools/StatTools/Jan21/V24_CR_12p9fb/hists*.root $ORIG_DIR\n"
        #script_text_cards += "combineCards.py Wmn=dc_Wmn.txt Wen=dc_Wen.txt > dc.txt\n"
        #script_text_cards += "combineCards.py ttWmn=vhbb_ttWmn_13TeV.txt ttWen=vhbb_ttWen_13TeV.txt whfWmn=vhbb_whfWmn_13TeV.txt whfWen=vhbb_whfWen_13TeV.txt wlfWmn=vhbb_wlfWmn_13TeV.txt wlfWen=vhbb_wlfWen_13TeV.txt WmnCat0=dc_WmnCat0.txt WmnCat1=dc_WmnCat1.txt WenCat0=dc_WenCat0.txt WenCat1=dc_WenCat1.txt > dc.txt\n"
        script_text_cards += "combineCards.py ttWmn=vhbb_ttWmn_13TeV.txt ttWen=vhbb_ttWen_13TeV.txt whfWmn=vhbb_whfWmn_13TeV.txt whfWen=vhbb_whfWen_13TeV.txt wlfWmn=vhbb_wlfWmn_13TeV.txt wlfWen=vhbb_wlfWen_13TeV.txt %s > dc.txt\n" % combineCards_text
        script_text_cards += "combine -M ProfileLikelihood --significance dc.txt -t -1 --expectSignal=1 > combine_output.txt\n"
        condor_runscript_cards = open("condor_runscript_cards.sh","w")
        condor_runscript_cards.write(script_text_cards)
        condor_runscript_cards.close()
        submit_text_cards = '''universe = vanilla
Executable = condor_runscript_cards.sh
Should_Transfer_Files = YES
Output = stdout_cards
Error  = stderr_cards
Log    = log_cards
Notification = never
WhenToTransferOutput=On_Exit
Queue  1
        '''
        submitfile_cards = open("submit_cards.txt","w")
        submitfile_cards.write(submit_text_cards)
        submitfile_cards.close()
        if (procStep == 1):
            os.system("condor_submit submit_cards.txt")
        if (procStep == 2):
            combOut = open("combine_output.txt","r")
            for line in combOut:
                 if (line.find("Significance") != -1):
                     combSens = float(line.split()[1])
                     print combSens
        os.chdir("../..")
       
    Sig[0] = float(combSens)
    bound_arr[len(bset)-1] = bset[len(bset)-1]
    print ncat[0],S_arr,B_arr,bound_arr,Sig[0]
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

print "boundary bin edges were"
print boundaries
print "best category boundaries are..."
print bestWP
if not doWPs:
    print "The corresponding signal efficiencies of each boundary edge is:"
    sigEffs = []
    for bound in bestWP:
        if (bound == xlow or bound == xhigh or bound == bestWP[len(bestWP)-1]): continue
        else: sigEffs.append(boundaryEff[bound])
    print sigEffs
print "with combined sensitivity: %f" % bestSig

if (nCat == 2 or nCat == 3):
    c1 = ROOT.TCanvas("c1","c1")
    if (nCat == 2):
        hSig.SetTitle("")
        if (doWPs):
            hSig.GetXaxis().SetTitle("WP Split Point")
        else: hSig.GetXaxis().SetTitle("V pT Split Point")
        hSig.GetYaxis().SetTitle("Combined Significance")
        hSig.Fill(100,0.734816) 
        hSig.GetXaxis().SetRangeUser(100,300)
        hSig.GetYaxis().SetRangeUser(0.5,0.8)
        hSig.Draw("hist")
        #hSig.Draw("pc")
    if (nCat == 3):
        hSig2.SetTitle("")
        if (doWPs):
            hSig2.GetXaxis().SetTitle("WP Split Point 1")
            hSig2.GetXaxis().SetTitle("WP Split Point 2")
        else: 
            hSig2.GetXaxis().SetTitle("V pT Split Point 1")
            hSig2.GetXaxis().SetTitle("V pT Split Point 2")
        hSig2.Draw("text")
    #raw_input()
    c1.SaveAs(sys.argv[3].replace(".root",".pdf"))
ofile.cd()
otree.Write()
ofile.Close()

