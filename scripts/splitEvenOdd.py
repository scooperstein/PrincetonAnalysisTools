import ROOT
import sys

ifilename = sys.argv[1]
ofiledir = sys.argv[2]
itree = ROOT.TChain("tree")
for ifilename in sys.argv[1].split(','):
    itree.Add(ifilename)

doVV = False
if (len(sys.argv) > 3):
    if (sys.argv[3] == "1"):
        doVV = True
if not doVV:
    print "copying signal even tree..."
    ofile_even_sig = ROOT.TFile("%s/ofile_even_sig.root" % ofiledir,"RECREATE")
    itree_even_sig = itree.CopyTree("Pass_nominal==1 && evt%2==0 && (sampleIndex==-12501 || sampleIndex==-12500)")
    ofile_even_sig.cd()
    itree_even_sig.Write()

    print "copying even bkg tree..."
    ofile_even_bkg = ROOT.TFile("%s/ofile_even_bkg.root" % ofiledir,"RECREATE")
    itree_even_bkg = itree.CopyTree("Pass_nominal==1 && evt%2==0 && sampleIndex>0")
    ofile_even_bkg.cd()
    itree_even_bkg.Write()

    print "copying odd sig tree..."
    ofile_odd_sig = ROOT.TFile("%s/ofile_odd_sig.root" % ofiledir,"RECREATE")
    itree_odd_sig = itree.CopyTree("Pass_nominal==1 && evt%2==1 && (sampleIndex==-12501 || sampleIndex==-12500)")
    ofile_odd_sig.cd()
    itree_odd_sig.Write()

    print "copying odd bkg tree..."
    ofile_odd_bkg = ROOT.TFile("%s/ofile_odd_bkg.root" % ofiledir,"RECREATE")
    itree_odd_bkg = itree.CopyTree("Pass_nominal==1 && evt%2==1 && sampleIndex>0")
    ofile_odd_bkg.cd()
    itree_odd_bkg.Write()

    print "done!"
    #ifile.Close()
    ofile_even_sig.Close()
    ofile_odd_sig.Close()
    ofile_even_bkg.Close()
    ofile_odd_bkg.Close()
else:
    # doVV
    print "copying vv even tree..."
    ofile_even_sig = ROOT.TFile("%s/ofile_even_vv.root" % ofiledir,"RECREATE")
    itree_even_sig = itree.CopyTree("Pass_nominal==1 && evt%2==0 && (sampleIndex>=3400&&sampleIndex<=3702)")
    ofile_even_sig.cd()
    itree_even_sig.Write()

    print "copying even bkg tree..."
    ofile_even_bkg = ROOT.TFile("%s/ofile_even_bkgvv.root" % ofiledir,"RECREATE")
    itree_even_bkg = itree.CopyTree("Pass_nominal==1 && evt%2==0 && (sampleIndex<3400||sampleIndex>3702)")
    ofile_even_bkg.cd()
    itree_even_bkg.Write()

    print "copying vv odd tree..."
    ofile_odd_sig = ROOT.TFile("%s/ofile_odd_vv.root" % ofiledir,"RECREATE")
    itree_odd_sig = itree.CopyTree("Pass_nominal==1 && evt%2==1 && (sampleIndex>=3400&&sampleIndex<=3702)")
    ofile_odd_sig.cd()
    itree_odd_sig.Write()

    print "copying odd bkg tree..."
    ofile_odd_bkg = ROOT.TFile("%s/ofile_odd_bkgvv.root" % ofiledir,"RECREATE")
    itree_odd_bkg = itree.CopyTree("Pass_nominal==1 && evt%2==1 && (sampleIndex<3400||sampleIndex>3702)")
    ofile_odd_bkg.cd()
    itree_odd_bkg.Write()

    print "done!"
    #ifile.Close()
    ofile_even_sig.Close()
    ofile_odd_sig.Close()
    ofile_even_bkg.Close()
    ofile_odd_bkg.Close()
