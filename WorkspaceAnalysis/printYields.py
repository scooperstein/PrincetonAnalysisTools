import ROOT

rates = "rate                                                        " # yields in datacard format
dc_string = "" # write datacard

ifile = ROOT.TFile("hists.root")

cats = ["WenLowPt","WenHighPt","WmnLowPt","WmnMidPt","WmnHighPt"]
cat_labels = ["ch1_Wenu","ch1_Wenu3","ch2_Wmunu","ch2_Wmunu2","ch2_Wmunu3"]
samples = ["ZH","WH","s_Top","Zj1b","TT","Zj0b","Wj0b","Wj1b","Wj2b","Zj2b"]

for cat in cats:
    print "Calculating Yields for Category: %s" % cat
    nBkgTot = 0.0;
    for sample in samples:
        nyield = ifile.Get("%s/%s" % (cat,sample)).Integral()
        print "%s = %.4f" % (sample, nyield )
        if (sample != "WH" and sample != "ZH"): nBkgTot += nyield
        # for readable aligment
        if (nyield < 10):
            rates += "%.4f      " % nyield
        elif (nyield < 100):
            rates += "%.4f     " % nyield
        elif (nyield < 1000):
            rates += "%.4f    " % nyield
        else:
            rates += "%.4f   " % nyield
    print "Total Background Yield: %.4f" % nBkgTot

print rates

dc_string += "imax %i number of bins\n" % len(cats)
dc_string += "jmax %i number of processes minus 1\n" % (len(samples) - 1)
dc_string += "kmax 0 number of nuisance parameters\n"
dc_string += "----------------------------------------------------------------------------------------------------------------------------------\n"
for i in range(len(cats)):
    dc_string += "shapes *           %s    workspace_%s.root %s:$PROCESS %s:$PROCESS_$SYSTEMATIC\n" % (cat_labels[i],cats[i],cats[i],cats[i])
dc_string += "----------------------------------------------------------------------------------------------------------------------------------\n"
dc_string += "bin          "
for label in cat_labels:
    dc_string += "%s    " % label
dc_string += "\n"
dc_string += "observation  "
for i in range(len(cats)):
    dc_string += "0.0000      "
dc_string += "\n"
dc_string += "bin                                                         "
for label in cat_labels:
    for sample in samples:
        dc_string += label
        for i in range(12 - len(label)):
            dc_string += " "
dc_string += "\n"
dc_string += "process                                                     "
for label in cat_labels:
    for sample in samples:
        dc_string += sample
        for i in range(12 - len(sample)):
            dc_string += " "
dc_string += "\n"
dc_string += "process                                                     "
for label in cat_labels:
    for i in range(len(samples)):
        dc_string += str(i)
        if (i < 10): dc_string += "           "
        else: dc_string += "          "
dc_string += "\n"
dc_string += rates

ofile = open("dc.txt","write")
ofile.write(dc_string)
ofile.close()
print dc_string
ifile.Close()
