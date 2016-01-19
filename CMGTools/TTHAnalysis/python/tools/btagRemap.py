#!/usr/bin/env python
from CMGTools.TTHAnalysis.treeReAnalyzer import ROOT, EventLoop, Module, Collection
import os.path, types
from array import array
from math import log, exp 

class BTagRemapper:
    def __init__(self,effs,csv,measForB="mujets",measForC="mujets",measForU="comb",shiftForBC="central",shiftForU="central"):
        self._flavs = ("b","c","l")
        self._wps   = ("L","M","T")
        self._loadMC(effs)
        self._makeData(csv,measForB,measForC,measForU,shiftForBC,shiftForU)
        hmc = ROOT.TH1F("hmc","hmc",4,array('f',[0,0.605,0.89,0.97,1.0]))
        hdt = ROOT.TH1F("hdt","hdt",4,array('f',[0,0.605,0.89,0.97,1.0]))
        remaps = {}
        weights = {}
        weightSplines = {}
        for flav in self._flavs:
            h2d = self._hflavs[flav]
            h2ds = [ (i+1,wp,self._heffs[(wp,flav)],self._heffsData[(wp,flav)]) for (i,wp) in enumerate(self._wps) ]
            for ix in xrange(1,h2d.GetNbinsX()+1):
                for iy in xrange(1,h2d.GetNbinsY()+1):
                    totmc, totdt = 0, 0
                    for iwp,wp,heff,heffData in reversed(h2ds):
                        effmc = min(max(0, heff.GetBinContent(ix,iy)),     1)
                        effdt = min(max(0, heffData.GetBinContent(ix,iy)), 1)
                        hmc.SetBinContent(iwp+1, max(0,effmc - totmc))
                        hdt.SetBinContent(iwp+1, max(0,effdt - totdt))
                        totmc += max(0,effmc - totmc)
                        totdt += max(0,effdt - totdt)
                        #if (ix == 1 and iy in (3,4,5)):
                        #    print "bin %2d, %2d, flav %s, wp %s (%d): heff mc %.3f, data %.3f, deff mc %.3f, data %.3f, sf %.3f, dsf %.3f" % (ix,iy,flav,wp,iwp,effmc,effdt,hmc.GetBinContent(iwp+1),hdt.GetBinContent(iwp+1),effdt/effmc if effmc else 1, hdt.GetBinContent(iwp+1)/hmc.GetBinContent(iwp+1) if hmc.GetBinContent(iwp+1) else 1)
                    hmc.SetBinContent(1, max(0,1-sum(hmc.GetBinContent(i) for i in (2,3,4))))
                    hdt.SetBinContent(1, max(0,1-sum(hdt.GetBinContent(i) for i in (2,3,4))))
                    hweight = hdt.Clone(""); hweight.SetDirectory(None)
                    for i in 1,2,3,4:
                        hweight.SetBinContent(i, hdt.GetBinContent(i)/hmc.GetBinContent(i) if  hmc.GetBinContent(i) else 1)
                        weights[(flav,ix,iy)] = hweight
                        weightSplines[(flav,ix,iy)] = ROOT.TSpline3(hweight)
                    for i in 1,2,3,4:
                        if hmc.GetBinContent(i) == 0: print "Warning: empty MC   wp bin %d for %2d, %2d, flav %s " % (i, ix, iy, flav)
                        if hdt.GetBinContent(i) == 0: print "Warning: empty Data wp bin %d for %2d, %2d, flav %s " % (i, ix, iy, flav)
                    #if (ix == 1 and iy in (3,4,5)):
                    #    print "bin %2d, %2d, flav %s, wp %s (%d): heff mc %.3f, data %.3f, deff mc %.3f, data %.3f, sf %.3f, dsf %.3f" % (ix,iy,flav,'F',0,1,1,hmc.GetBinContent(1),hdt.GetBinContent(1), 1, hdt.GetBinContent(1)/hmc.GetBinContent(1) if  hmc.GetBinContent(i) else 1)
                    remaps[(flav,ix,iy)] = ROOT.DistributionRemapper(hmc,hdt,True,False,True,False)#(ix == 1 and iy in (3,4,5)))
                    if False and (ix == 1 and iy in (3,4,5)):
                        func = remaps[(flav,ix,iy)]
                        hist   = weights[(flav,ix,iy)]
                        spline = weightSplines[(flav,ix,iy)]
                        for x in [0, 0.2, 0.4, 0.6, 0.7, 0.8, 0.85, 0.9, 0.94, 0.96, 0.98, 1.00]:
                            print "\t\t%.3f -> %.4f  (h %.4f, s %.4f)" % (x,func.Eval(x),hist.GetBinContent(min(max(1, hist.FindBin(x)), 4)),spline.Eval(x))
        self._remaps = remaps
        self._weights = weights
        self._weightSplines = weightSplines
    def _loadMC(self,effs):
        self._heffs  = {}
        self._hflavs = {}
        tfile = ROOT.TFile.Open(os.path.expandvars(effs))
        for flav in self._flavs:
            for wp in self._wps:
                hist = tfile.Get("CSV%s_vs_ptAEta_TT%s" % (wp,flav))
                if not hist: raise RuntimeError, "Missing CSV%s_vs_ptAEta_TT%s" % (wp,flav)
                if hist.ClassName() == "TProfile2D":
                    hist = hist.ProjectionXY()
                hist = hist.Clone(); hist.SetDirectory(None)
                self._heffs[(wp,flav)] = hist
                self._hflavs[flav] = hist
        tfile.Close()
    def _makeData(self,csv,measForB,measForC,measForU,shiftForBC,shiftForU):
        csvfile = open(os.path.expandvars(csv),'r')
        header = None
        self._heffsData = {}
        for key,val in self._heffs.iteritems():
            hdata = val.Clone(val.GetName()+"Data") 
            hdata.SetDirectory(None)
            self._heffsData[key] = hdata
        done = []
        for op in self._wps:
            for flav in self._flavs:
                h2d = self._heffsData[(op,flav)]
                done += [ ((op,flav,ix,iy),0) for ix in xrange(1,h2d.GetNbinsX()+1) for iy in xrange(1,h2d.GetNbinsY()+1) ]
        done = dict(done)
        for line in csvfile:
            fields = [ x.strip() for x in line.split(",") ]
            if not header:
                header = fields
                continue
            (op,meas,syst,flav,etal,etah,ptl,pth,dl,dh,formula) = fields
            (op,flav) = map(int, (op,flav))
            (etal,etah,ptl,pth) = map(float, (etal,etah,ptl,pth))
            if flav == 2:
                if meas != measForU or syst != shiftForU: continue
            elif flav == 1:
                if meas != measForC or syst != shiftForBC: continue
            elif flav == 0:
                if meas != measForB or syst != shiftForBC: continue
            if op >= len(self._wps): continue
            if flav >= len(self._flavs): continue
            sfpt = eval("lambda x : "+formula.replace('"',''))
            h2d = self._heffsData[(self._wps[op],self._flavs[flav])]
            for ix in xrange(1,h2d.GetNbinsX()+1):
                eta = h2d.GetXaxis().GetBinCenter(ix)
                if not(etal <= eta and eta < etah): continue
                for iy in xrange(1,h2d.GetNbinsY()+1):
                    pt  = h2d.GetYaxis().GetBinCenter(iy)
                    if not(ptl <= pt and pt < pth): continue
                    dataeff = sfpt(pt)*h2d.GetBinContent(ix,iy)
                    #if (ix == 1 and iy == 3):
                    #    print "for %s, %s, pt %7.0f, eta %5.2f, sf = %.3f (%s), mc = %.3f, data = %.3f" % (self._wps[op],self._flavs[flav],pt,eta,sfpt(pt), formula, h2d.GetBinContent(ix,iy), dataeff)
                    h2d.SetBinContent(ix,iy,dataeff)
                    done[(self._wps[op],self._flavs[flav],ix,iy)] += 1
        for (op,flav,ix,iy),n in done.iteritems():
            if n == 1: continue
            h2d = self._heffsData[(op,flav)]
            eta = h2d.GetXaxis().GetBinCenter(ix)
            pt  = h2d.GetYaxis().GetBinCenter(iy)
            if h2d.GetBinContent(ix,iy) < 1e-4:
                print "Warning, no EFF for %s, %s, %s, %s, %.6f" % (op,flav,pt,eta, h2d.GetBinContent(ix,iy))
            #if n == 0: print "Warning, no SF for %s, %s, %s, %s" % (op,flav,pt,eta)
            if n >  1: 
                print "ERROR, duplicate SFs for %s, %s, %s, %s" % (op,flav,pt,eta)
                raise RuntimeError
    def remap(self,flav,pt,eta,disc):
        h2d = self._hflavs[flav]
        ieta = min(max(1, h2d.GetXaxis().FindBin(abs(eta))), h2d.GetNbinsX())
        ipt  = min(max(1, h2d.GetYaxis().FindBin(pt)      ), h2d.GetNbinsY())
        return self._remaps[(flav,ieta,ipt)].Eval(min(max(0,disc),1.))
    def sf(self,flav,pt,eta,wp):
        h2d = self._heffsData[(wp,flav)]
        h2m = self._heffs[(wp,flav)]
        ieta = min(max(1, h2d.GetXaxis().FindBin(abs(eta))), h2d.GetNbinsX())
        ipt  = min(max(1, h2d.GetYaxis().FindBin(pt)      ), h2d.GetNbinsY())
        return h2d.GetBinContent(ieta,ipt)/h2m.GetBinContent(ieta,ipt)
    def reweight(self,flav,pt,eta,disc,spline=False):
        h2d = self._hflavs[flav]
        ieta = min(max(1, h2d.GetXaxis().FindBin(abs(eta))), h2d.GetNbinsX())
        ipt  = min(max(1, h2d.GetYaxis().FindBin(pt)      ), h2d.GetNbinsY())
        if spline:
            return self._weightSplines[(flav,ieta,ipt)].Eval(min(max(0.3,disc),1.))
        else:
            hist = self._weights[(flav,ieta,ipt)]
            return hist.GetBinContent(min(max(1,hist.FindBin(disc)),4))
            
BTagRemap74X = lambda : BTagRemapper("$CMSSW_BASE/src/CMGTools/TTHAnalysis/data/btag/bTagEff_TTMC_Spring15.root",
                                     "$CMSSW_BASE/src/CMGTools/TTHAnalysis/data/btag/CSVv2_74X.csv")
class BTagRemapFriend:
    def __init__(self,remap,jets="Jet",inlabel="btagCSV",outlabel="btagCSVRemap",maxjets=20,mcOnly=True):
        self.jets = jets
        self.label = outlabel
        self.blabel = inlabel
        self.maxJets = maxjets
        self.mcOnly = mcOnly
        self._remapper = (remap() if type(remap) == types.FunctionType else remap)
    def listBranches(self):
        return [ ("n"+self.jets,"I"), (self.jets+"_"+self.label,"F",self.maxJets,"n"+self.jets) ]
    def __call__(self,event):
        jets = Collection(event,self.jets)
        ret = { 'n%s'%self.jets : len(jets) }
        ret[self.jets+"_"+self.label] = [ self.remap(j) for j in jets ] 
        return ret
    def remap(self,jet):
        btag = getattr(jet,self.blabel)
        if self.mcOnly and jet.mcPt <= 0: return btag
        if   abs(jet.mcFlavour) == 5: flav = "b"
        elif abs(jet.mcFlavour) == 4: flav = "c"
        else                        : flav = "l"
        return self._remapper.reweght(flav,jet.pt,jet.eta,btag)

class BTagMyReweightFriend:
    def __init__(self,remap,jets="Jet",inlabel="btagCSV",outlabel="btagCSVMyW",maxjets=20,mcOnly=True,spline=True):
        self.jets = jets
        self.label = outlabel
        self.blabel = inlabel
        self.maxJets = maxjets
        self.mcOnly = mcOnly
        self._remapper = (remap() if type(remap) == types.FunctionType else remap)
        self._spline = spline
    def listBranches(self):
        return [ ("n"+self.jets,"I"), (self.jets+"_"+self.label,"F",self.maxJets,"n"+self.jets) ]
    def __call__(self,event):
        jets = Collection(event,self.jets)
        ret = { 'n%s'%self.jets : len(jets) }
        ret[self.jets+"_"+self.label] = [ self.weight(j) for j in jets ] 
        return ret
    def weight(self,jet):
        btag = getattr(jet,self.blabel)
        if self.mcOnly and jet.mcPt <= 0: return btag
        if   abs(jet.mcFlavour) == 5: flav = "b"
        elif abs(jet.mcFlavour) == 4: flav = "c"
        else                        : flav = "l"
        return self._remapper.reweight(flav,jet.pt,jet.eta,btag,spline=self._spline)

class BTagMyLepReweightFriend(BTagMyReweightFriend):
    def __init__(self,remap,jets="LepGood",inlabel="jetBTagCSV",outlabel="jetBTagCSVMyW",maxjets=20,spline=True):
        BTagMyReweightFriend.__init__(self,remap,jets=jets,inlabel=inlabel,outlabel=outlabel,maxjets=maxjets,mcOnly=True,spline=spline)
    def weight(self,lep):
        flav = abs(lep.mcMatchAny)
        if   flav == 5: flav = 'b'
        elif flav == 4: flav = 'c'
        else:           flav = 'l'
        jetpt = lep.pt/lep.jetPtRatiov2
        return self._remapper.reweight(flav, jetpt, lep.eta, getattr(lep,self.blabel), spline=self._spline)
 
if __name__ == '__main__':
    from sys import argv
    file = ROOT.TFile(argv[1])
    tree = file.Get("tree")
    tree.vectorTree = True
    class Tester(Module):
        def __init__(self, name):
            Module.__init__(self,name,None)
            self.sf = BTagRemapFriend(BTagRemap74X)
            self._wps = [ 0.605, 0.89, 0.97 ]
            self._pcountA = 0
            self._pcount0 = [0 for w in self._wps]
            self._pcount1 = [0 for w in self._wps]
            self._weightedSums = dict((s,[0.,0.,0.,0.]) for s in (True,False))
            self._sfAvg   = [[],[],[]]
        def analyze(self,ev):
            if self._pcountA < 10:
                print "\nrun %6d lumi %4d event %d: jets %d" % (ev.run, ev.lumi, ev.evt, ev.nJet25)
            ret = self.sf(ev)["Jet_btagCSVRemap"]
            jets = Collection(ev,"Jet")
            for i,j in enumerate(jets):
                d = min(max(0,j.btagCSV),1)
                sf  = [ self.sf._remapper.sf('b',j.pt,j.eta,wp) for wp in self.sf._remapper._wps ] 
                w_h = self.sf._remapper.reweight('b',j.pt,j.eta,d,spline=False)
                w_s = self.sf._remapper.reweight('b',j.pt,j.eta,d,spline=True)
                if self._pcountA < 10:
                    print "\t%8.2f %+5.2f %+3d %.3f -> %.3f (%.3f, %.3f, %.3f) (%.3f, %.3f)" % (
                                    j.pt, j.eta, j.mcFlavour, d, ret[i],
                                    sf[0], sf[1], sf[2],
                                    w_h, w_s)
                if j.pt > 30 and abs(j.eta) < 2.4 and abs(j.mcFlavour) == 5:
                    self._pcountA += 1.
                    for iw,w in enumerate(self._wps):
                        if j.btagCSV > w: self._pcount0[iw] += 1.
                        if ret[i]    > w: self._pcount1[iw] += 1.
                        self._sfAvg[iw].append(sf[iw])
                    for s in True, False:
                        self._weightedSums[s][-1] += w_s if s else w_h
                        for iw,w in enumerate(self._wps):
                            if d > w: self._weightedSums[s][iw] += w_s if s else w_h
            if self._pcountA < 10:
                print ""
        def done(self):
            for iw,w in enumerate(self._wps):
                eff_pre  = self._pcount0[iw]/self._pcountA
                eff_post = self._pcount1[iw]/self._pcountA 
                sfavg    = sum(self._sfAvg[iw])/len(self._sfAvg[iw])
                eff_w    = self._weightedSums[False][iw]/self._weightedSums[False][-1]
                eff_s    = self._weightedSums[True ][iw]/self._weightedSums[True ][-1]
                print " for WP %.3f, eff(pre) = %.3f, eff(post) = %.3f, SF = %.3f  (real avg %.3f), sf(w) = %.3f, sf(s) = %.3f " % (w, 
                            eff_pre, eff_post, eff_post/eff_pre,
                            sfavg, 
                            eff_w/eff_pre, eff_s/eff_pre)
    T = Tester("tester")
    el = EventLoop([ T ])
    el.loop([tree], maxEvents = 2500)  
    T.done()
