import ROOT

line = "rate                                                        " # yields in datacard format

ifile = ROOT.TFile("hists.root")

for cat in ["WenLowPt","WenHighPt","WmnLowPt","WmnMidPt","WmnHighPt"]:
    print "Calculating Yields for Category: %s" % cat
    nBkgTot = 0.0;
    for sample in ["ZH","WH","s_Top","Zj1b","TT","Zj0b","Wj0b","Wj1b","Wj2b","Zj2b"]:
        nyield = ifile.Get("%s/%s" % (cat,sample)).Integral()
        print "%s = %.4f" % (sample, nyield )
        if (sample != "WH" and sample != "ZH"): nBkgTot += nyield
        if (nyield < 10):
            line += "%.4f      " % nyield
        elif (nyield < 100):
            line += "%.4f     " % nyield
        elif (nyield < 1000):
            line += "%.4f    " % nyield
        else:
            line += "%.4f   " % nyield
    print "Total Background Yield: %.4f" % nBkgTot

print line
ifile.Close()
