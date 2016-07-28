from ROOT import *
import sys

gStyle.SetOptStat(0)

ifile = TFile(sys.argv[1])
#channels = ["ttWmn","ttWen","whfWmn","whfWen","wlfWmn","wlfWen","WmnHighPt","WenHighPt"]
#channels = ["ttWmn","ttWen","whfWmn","whfWen","wlfWmn","wlfWen"]
samples = ["Wj2b","Wj1b","Wj0b","TT","s_Top"]
#channels = ["whf","wlf","tt"]
#channels = ["whfHighPt","wlfHighPt","ttHighPt"]
#channels = ["whfLowPt","wlfLowPt","ttLowPt"]
#channels = ["whfLowPt","whfHighPt","wlfLowPt","wlfHighPt","ttLowPt","ttHighPt"]
channels = ["HighPt","whf","wlf","tt"]
#channels = ["HighPt"]

bkgYield_pre = {}
bkgYield_post = {}
for chan in channels:
    bkgYield_pre[chan] = 0.
    bkgYield_post[chan] = 0.

for sample in samples:
    prefit_integral = 0.
    postfit_integral = 0.
    for chan in channels:
        prefit_integral_em = 0.
        postfit_integral_em = 0.
        for cat in ["e","m"]:
            channel = chan + "W" + cat + "n"
            #print sample,channel
            if (chan.find("HighPt")!=-1):
                #channel = "W" + cat + "n" + chan 
                index = chan.find("HighPt")
                channel = chan[:index] + "W" + cat + "n" + chan[index:] 
            if (chan.find("LowPt")!=-1):
                #channel = "W" + cat + "n" + chan 
                index = chan.find("LowPt")
                channel = chan[:index] + "W" + cat + "n" + chan[index:] 
            #if (channel == "ttWmnHighPt"): 
            #    channel = "ttWmnHigh"
            #print sample,channel
            h_pre = ifile.Get("shapes_prefit/%s/%s" % (channel,sample))
            h_post = ifile.Get("shapes_fit_s/%s/%s" % (channel,sample))

            prefit_integral += h_pre.Integral()
            postfit_integral += h_post.Integral()
            #prefit_integral_em += h_pre.Integral()
            #postfit_integral_em += h_post.Integral()
            prefit_integral_em += h_pre.Integral(14,15)
            postfit_integral_em += h_post.Integral(14,15)
            #print "%s %s: %f" % (sample, channel, (h_post.Integral() / h_pre.Integral()) )
            
            #bkgYield_pre[chan] += h_pre.Integral()
            #bkgYield_post[chan] += h_post.Integral() 
            bkgYield_pre[chan] += h_pre.Integral(14,15)
            bkgYield_post[chan] += h_post.Integral(14,15) 
        #prefit_integral += h_pre.Integral()
        #postfit_integral += h_post.Integral()
        print "%s %s: %f" % (sample, chan, (postfit_integral_em) )
        #print "%s %s: %f" % (sample, chan, (postfit_integral_em / prefit_integral_em) )
    print "%s: %f" % (sample, (postfit_integral / prefit_integral) )
print "Total background yield post/pre-fit..."
total_pre = 0.
total_post = 0.
for chan in channels: 
    total_pre += bkgYield_pre[chan]
    total_post += bkgYield_post[chan]
    #print "%s: %f" % (chan, (bkgYield_post[chan]/bkgYield_pre[chan]))
    print "%s: %f" % (chan, (bkgYield_post[chan]))
    #print "%s: %f" % (chan, (bkgYield_pre[chan]))
print "Overall: %f" % (total_post/total_pre)


