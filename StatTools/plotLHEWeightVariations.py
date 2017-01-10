import ROOT
import sys
from numpy import mean, sqrt, square, array, zeros

#ROOT.gROOT.SetBatch(True)
#ROOT.gStyle.SetOptStat(0)

channel = sys.argv[2]
ifile = ROOT.TFile(sys.argv[1])
#ifile = ROOT.TFile("hists_%s.root" % channel)
#samples = ["s_Top","VVHF","VVLF","Zj0b","Zj1b","Zj2b","QCD"]
samples = ["TT","Wj0b","Wj1b","Wj2b","Zj0b","Zj1b","Zj2b","VVLF","VVHF","ZH","WH"]
#samples = ["TT"]
#samples = ["WH"]
nBins = -1

stack_pdfs = {}
for sample in samples:
    print sample
    nPdfs = 101
    stack_pdfs[sample] = []
    for i in range(nPdfs):
        print i
        #stack = ROOT.THStack()
        madeStack = False
        stack = ifile.Get("BDT_%s_WH" % channel )
        if not madeStack:
                stack = ifile.Get("BDT_%s_%s_CMS_vhbb_LHE_weights_pdf_%iUp" % (channel,sample,i) )
                print "BDT_%s_%s_CMS_vhbb_LHE_weights_pdf_%iUp" % (channel,sample,i)
                madeStack = True
                nBins = stack.GetNbinsX()
                print "First time, ",stack.Integral()
        else:
                hist = ifile.Get("BDT_%s_%s_CMS_vhbb_LHE_weights_pdf_%iUp" % (channel,sample,i) )
                print hist.GetName()
                print stack.Integral(),hist.Integral()
                stack.Add(hist)
                print stack.Integral()
        if (stack.Integral() == 0.):
            print "Be careful! The integral is 0.0 for sample %s and pdf # %i" % (sample,i)
        stack_pdfs[sample].append(stack)

print "time to make plots"
canv = ROOT.TCanvas("canv","canv")
ofile = ROOT.TFile("LHEPDF_rms_over_mean.root","RECREATE")
ifile.cd()
for sample in samples:
    graph = ROOT.TGraphErrors(sample,sample)
    graph2 = ROOT.TGraphErrors(sample+"_2",sample+"_2")
    for i in range(1,nBins+1):
        binMean = 0
        binRMS = 0
        vals = []
        xval = 0.
        for stack in stack_pdfs[sample]:
            vals.append(stack.GetBinContent(i))
            xval = stack.GetBinCenter(i)
        print vals
        vals_arr = array(vals)
        binMean = mean(vals_arr)
        for index in range(len(vals_arr)):
            vals_arr[index] = vals_arr[index] - binMean
        print vals_arr
        binRMS = sqrt(mean(square(vals_arr)))
        graph.SetPoint(i-1,xval, binMean)
        graph.SetPointError(i-1,0.5*stack.GetBinWidth(i), binRMS)
        if (binMean!=0):
            graph2.SetPoint(i-1,xval,binRMS/binMean)
        else:
            graph2.SetPoint(i-1,xval,0.0)
        print xval,binMean,binRMS

    hist_nominal = ifile.Get("BDT_%s_%s" % (channel,sample))
    #for i in range(len(samples)):
    #    sample = samples[i]
    #    if (i == 0):
    #        hist_nominal = ifile.Get("BDT_%s_%s" % (channel,sample))
    #    else:
    #        hist_tmp = ifile.Get("BDT_%s_%s" % (channel,sample))
    #        hist_nominal.Add(hist_tmp)

    hist_nominal.SetTitle("")
    hist_nominal.GetXaxis().SetTitle("BDT Score")

    graph.SetTitle("")
    graph.GetXaxis().SetTitle("BDT Score")
    graph2.SetTitle("")
    graph2.GetXaxis().SetTitle("BDT Score")

    hist_nominal.SetLineColor(ROOT.kBlue)
    hist_nominal.SetLineWidth(2)
    graph.SetMarkerColor(ROOT.kRed)
    graph.SetMarkerStyle(22)
    graph.SetFillStyle(3005)
    graph.SetFillColor(ROOT.kBlack)
    graph2.SetMarkerColor(ROOT.kRed)
    graph2.SetMarkerStyle(22)
    graph2.SetFillStyle(3005)
    graph2.SetFillColor(ROOT.kBlack)

    #hist_nominal.SetMinimum(0.01)
    #hist_nominal.SetMinimum(5)

    #hist_nominal.Draw("hist")
    #hist_nominal.Draw("ep same")
    #graph.Draw("ep same")
    #graph.Draw("ep2 same")

    graph2.Draw("ape")
    #graph.Draw("ape")
    print stack.Integral(),hist_nominal.Integral()
    ofile.cd()
    hist = hist_nominal.Clone("PDF_rms_over_mean_%s_%s" % (channel,sample))
    hist.Reset()
    for i in range(nBins):
        x = zeros(1, dtype=float)
        y = zeros(1, dtype=float)
        graph2.GetPoint(i,x,y)
        hist.Fill(x,y)
    #hist.Write("PDF_rms_over_mean_%s_%s" % (channel,sample))
    hist.Write()    

    #leg = ROOT.TLegend(0.15,0.65,0.4,0.9)
    leg = ROOT.TLegend(0.65,0.65,0.9,0.9)
    #leg.AddEntry(hist_nominal,"NOMINAL (WHbb Signal)")
    #leg.AddEntry(hist_nominal,"NOMINAL (Sum Bkg.)")
    #leg.AddEntry(graph, "Mean of PDF variations","p")
    #leg.AddEntry(graph, "RMS of PDF variations","f")
    leg.AddEntry(graph2, "PDF Variation RMS / Mean (Sum Bkg.)","p")
    #leg.AddEntry(graph, "PDF Variation RMS / Mean (WHbb signal)","p")
    leg.Draw("same")

    canv.SetLogy(False)
    canv.SaveAs("PDF_variations_%s_%s_rms_over_mean.pdf" % (channel,sample))
    canv.SaveAs("PDF_variations_%s_%s_rms_over_mean.png" % (channel,sample))
    canv.SetLogy(True)
    canv.SaveAs("PDF_variations_%s_%s_rms_over_mean_log.pdf" % (channel,sample))
    canv.SaveAs("PDF_variations_%s_%s_rms_over_mean_log.png" % (channel,sample))

    hist_nominal.Draw("hist")
    hist_nominal.Draw("ep same")
    #graph.Draw("ep same")
    graph.Draw("ep2 same")
    leg2 = ROOT.TLegend(0.65,0.65,0.9,0.9)
    leg2.AddEntry(hist_nominal,"NOMINAL %s" % sample)
    leg2.AddEntry(graph, "Mean of PDF variations","p")
    leg2.AddEntry(graph, "RMS of PDF variations","f")
    leg2.Draw("same")
    canv.SetLogy(False)
    canv.SaveAs("PDF_variations_%s_%s.pdf" % (channel,sample))
    canv.SaveAs("PDF_variations_%s_%s.png" % (channel,sample))
    canv.SetLogy(True)
    canv.SaveAs("PDF_variations_%s_%s_log.pdf" % (channel,sample))
    canv.SaveAs("PDF_variations_%s_%s_log.png" % (channel,sample))

    #raw_input()
    #canv.Close()
    #graph.Clear()


