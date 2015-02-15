#!/usr/bin/env python
from CMGTools.TTHAnalysis.plotter.mcAnalysis import *
import re, sys, os, os.path
systs = {}

from optparse import OptionParser
parser = OptionParser(usage="%prog [options] mc.txt cuts.txt var bins systs.txt ")
addMCAnalysisOptions(parser)
parser.add_option("-o",   "--out",    dest="outname", type="string", default="all", help="output name") 
parser.add_option("--od", "--outdir", dest="outdir", type="string", default=None, help="output name") 
parser.add_option("-v", "--verbose",  dest="verbose",  default=0,  type="int",    help="Verbosity level (0 = quiet, 1 = verbose, 2+ = more)")
parser.add_option("-O", "--poly-order", dest="polyorder", default=5, type="int", help="Polynomial order")
parser.add_option("--rebin", dest="rebin", default=1, type="int", help="Rebin factor")

(options, args) = parser.parse_args()
options.weight = True
options.final  = True

mca  = MCAnalysis(args[0],options)
cuts = CutsFile(args[1],options)

binname = options.outname
outdir  = options.outdir+"/" if options.outdir else ""

print "Running on dataset"
report = mca.getPlotsRaw("x", args[2], args[3], cuts.allCuts())

print "Creating workspace"
workspace = ROOT.RooWorkspace("w","w");
workspace.imp = getattr(workspace,'import')
workspace.factory("MH[18.5]")
if args[3][0] != "[":
    nbins,xmin,xmax = args[3].split(",")
    workspace.factory("m4l[%s,%s]" %(xmin,xmax))
    workspace.var("m4l").setBins(int(nbins))
else:
    raise RuntimeError, "Variable binning not supported"
workspace.factory("Gaussian::sig_Shape(m4l,MH,sigma[0.15])")
workspace.factory("Bernstein::bkg_Shape(m4l,{%s})" % (",".join("bkgCoeff%d[0,10]" % (i+1) for i in xrange(options.polyorder-1))))
#workspace.factory("bkg_Shape_norm[%g,0,%g]"%(report['data'].Integral(),2*report['data'].Integral()))

obs = ROOT.RooDataHist("data_obs","data_obs",ROOT.RooArgList(workspace.var("m4l")),report['data'])
workspace.imp(obs)

workspace.Print("")

workspace.writeToFile(binname+".input.root")
workspace.pdf("bkg_Shape").fitTo(obs)
print "Wrote to ",binname+".input.root"

card = open(binname+".card.txt","w")
card.write("""# test card
imax 1
jmax 1
kmax *
-------------------------
shapes *        * {card}.input.root w:$PROCESS_Shape
shapes data_obs * {card}.input.root w:data_obs
-------------------------
bin {card}
observation {events}
-------------------------
bin      {card}   {card}
process  sig      bkg
process  0        1
rate     100      {events}
-------------------------
bn lnU   -        2 
""".format(card=binname,events=report['data'].Integral()))
card.close()
print "Wrote to ",binname+".card.txt"
os.system("cat %s.card.txt | sed 's/^/\\t/'" % binname)

ROOT.gROOT.ProcessLine(".x tdrstyle.cc")

drebin = report['data'].Clone("rebinned")
drebin.Rebin(options.rebin)
obsRebin = ROOT.RooDataHist("data_obs","data_obs",ROOT.RooArgList(workspace.var("m4l")),drebin)
print "bins: ",report['data'].GetNbinsX()/options.rebin
c1 = ROOT.TCanvas("c1","c1")
plot = workspace.var("m4l").frame(ROOT.RooFit.Bins(report['data'].GetNbinsX()/options.rebin))
obsRebin.plotOn(plot)
workspace.pdf("bkg_Shape").plotOn(plot)
plot.Draw()
c1.Print(binname+".fit_b.png")


workspace.factory("SUM::model_sb(f_s[0.01,-0.1,0.1]*sig_Shape,bkg_Shape)");
workspace.var("MH").setConstant(False)
workspace.var("MH").setMin(17)
workspace.var("MH").setMax(20)
workspace.pdf("model_sb").fitTo(obs)
drebin2 = report['data'].Clone("rebinned2")
drebin2.Rebin(max(1,options.rebin/2))
obsRebin = ROOT.RooDataHist("data_obs","data_obs",ROOT.RooArgList(workspace.var("m4l")),drebin2)
plot = workspace.var("m4l").frame(ROOT.RooFit.Bins(drebin2.GetNbinsX()))
obsRebin.plotOn(plot)
workspace.pdf("model_sb").plotOn(plot, ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Components("sig_Shape"), ROOT.RooFit.LineColor(ROOT.kRed))
workspace.pdf("model_sb").plotOn(plot, ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Components("bkg_Shape"), ROOT.RooFit.LineStyle(2))
workspace.pdf("model_sb").plotOn(plot, ROOT.RooFit.LineColor(ROOT.kBlue))
plot.Draw()
c1.Print(binname+".fit_sb.png")
c1.Print(binname+".fit_sb.root")

