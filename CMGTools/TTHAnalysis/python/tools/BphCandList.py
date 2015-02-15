#!/usr/bin/env python
from CMGTools.TTHAnalysis.treeReAnalyzer import *
from array import array
from glob import glob
import os.path

import os, ROOT


class CandList:
    def __init__(self,evlist,tolerance=0.25):
        self.tolerance = tolerance
        self.candidates = {}
        for line in open(evlist):
            fields = line.split()
            if len(fields) != 4: continue
            (run,lumi,event) = int(fields[0]),int(fields[1]),int(fields[2])
            mass = float(fields[3])
            key = (run,lumi,event)
            if key not in self.candidates:
                self.candidates[key] = [mass]
            else:
                self.candidates[key].append(mass)
        print "Loaded from %s: %d candidates in %d events" % (evlist,sum(len(v) for v in self.candidates.itervalues()),len(self.candidates))
    def listBranches(self):
        return [ ("nCand","I"), ("Cand_isKai","I",80,"nCand") ]
    def __call__(self,event):
        cand = Collection(event,"Cand","nCand")
        rle = (event.run,event.lumi,event.evt)
        isKai = [0 for c in cand]
        ret = { 'nCand' : len(cand) }
        if rle in self.candidates:
            masses = self.candidates[rle]
            for i,c in enumerate(cand):
                if min(abs(m-c.mass) for m in masses) < self.tolerance:
                    isKai[i] = 1
        ret['Cand_isKai'] = isKai
        return ret

class CandList4Mu:
    def __init__(self,evlist,mtolerance=0.25,xyztolerance=0.25,verbose=False,vname="isKai"):
        self.mtolerance = mtolerance
        self.xyztolerance = xyztolerance
        self.verb = verbose
        self.vname = vname
        self.candidates = {}
        self.evlist = evlist
        for line in open(evlist):
            fields = line.split()
            if len(fields) != 3+4*3+1: continue
            (run,lumi,event) = int(fields[0]),int(fields[1]),int(fields[2])
            mass = float(fields[-1])
            pxyz = []
            for i in xrange(4):
                px,py,pz = float(fields[3+i*3]),float(fields[4+i*3]),float(fields[5+i*3])
                pxyz.append((px,py,pz)) 
            key = (run,lumi,event)
            if key not in self.candidates:
                self.candidates[key] = [(mass,pxyz)]
            else:
                self.candidates[key].append((mass,pxyz))
        self.ncands = sum(len(v) for v in self.candidates.itervalues())
        self.nevts  = len(self.candidates)
        self.fcands = 0
        self.pcands = 0
        self.fevts  = 0
        print "Loaded from %s: %d candidates in %d events" % (evlist,sum(len(v) for v in self.candidates.itervalues()),len(self.candidates))
    def listBranches(self):
        return [ ("nCand","I"), ("Cand_"+self.vname,"I",2000,"nCand") ]
    def __call__(self,event):
        cand = Collection(event,"Cand","nCand")
        lep = Collection(event,"Lep","nLep")
        rle = (event.run,event.lumi,event.evt)
        isKai = [0 for c in cand]
        ret = { 'nCand' : len(cand) }
        if rle in self.candidates:
            self.fevts += 1
            nkai = len(self.candidates[rle])
            self.pcands += nkai
            #print "event %s with %d cands (kai %d)" % (rle,event.nCand,nkai)
            revkai = [0 for ikai in xrange(nkai)]
            for i,c in enumerate(cand):
                pxyz = []
                for imu in xrange(4):
                    p4 = lep[int(getattr(c,"l%dh"%(imu+1)))].p4()
                    pxyz.append((p4.X(),p4.Y(),p4.Z()))
                if self.verb: print "Our candidate mass %6.2f, muons: %s" % (c.mass, ", ".join("{ %+6.2f  %+6.2f  %+6.2f }" % p3 for p3 in pxyz))
                for ikai,(kmass,kpxyz) in enumerate(self.candidates[rle]): 
                    if abs(kmass-c.mass) < self.mtolerance*max(c.mass/10,1):
                        if self.verb: print "     Kai cand mass %6.2f, muons: %s" % (kmass, ", ".join("{ %+6.2f  %+6.2f  %+6.2f }" % p3 for p3 in kpxyz))
                        for i0 in xrange(4):
                            if self.verb: print "        dist for cand[0] vs kai[%d]: %f " % (i0, max(abs(a-b) for (a,b) in zip(pxyz[0],kpxyz[i0])))
                            if max(abs(a-b) for (a,b) in zip(pxyz[0],kpxyz[i0])) > self.xyztolerance: continue
                            for i1 in xrange(4):
                                if i1 == i0: continue 
                                if self.verb: print "           dist for cand[1] vs kai[%d]: %f " % (i1, max(abs(a-b) for (a,b) in zip(pxyz[1],kpxyz[i1])))
                                if max(abs(a-b) for (a,b) in zip(pxyz[1],kpxyz[i1])) > self.xyztolerance: continue
                                for i2 in xrange(4):
                                    if i2 == i0 or i2 == i1: continue 
                                    i3 = (0+1+2+3)-(i0+i1+i2)
                                    if self.verb: print "               dist for cand[2] vs kai[%d]: %f " % (i2, max(abs(a-b) for (a,b) in zip(pxyz[2],kpxyz[i2])))
                                    if self.verb: print "               dist for cand[3] vs kai[%d]: %f " % (i3, max(abs(a-b) for (a,b) in zip(pxyz[3],kpxyz[i3])))
                                    if max(abs(a-b) for (a,b) in zip(pxyz[2],kpxyz[i2])) > self.xyztolerance: continue
                                    if max(abs(a-b) for (a,b) in zip(pxyz[3],kpxyz[i3])) > self.xyztolerance: continue
                                    if self.verb: print "               ---> matched, dmax = %.6f dm = %.6f " % ( max(
                                                                max(abs(a-b) for (a,b) in zip(pxyz[0],kpxyz[i0])),
                                                                max(abs(a-b) for (a,b) in zip(pxyz[1],kpxyz[i1])),
                                                                max(abs(a-b) for (a,b) in zip(pxyz[2],kpxyz[i2])),
                                                                max(abs(a-b) for (a,b) in zip(pxyz[3],kpxyz[i3]))),
                                                                abs(kmass-c.mass)/max(c.mass/10,1)
                                                          )
                                    isKai[i] = 1
                                    revkai[ikai] += 1
                                    if revkai[ikai] == 1: self.fcands += 1
                                    break
                                if isKai[i]: break
                            if isKai[i]: break
                        if isKai[i]: break
                    if isKai[i]: break
            nmatch = sum(isKai)
            #print "In event with %d cands, kai has %d, we match %d (%s, %s)" % ( event.nCand, nkai, nmatch, isKai, revkai)
            if not self.verb and ((nmatch != nkai) or any([nuse for nuse in revkai if nuse>1])):
                print "\nAttention, in event with %d cands, kai has %d, we match %d (%s, %s)" % ( event.nCand, nkai, nmatch, isKai, revkai)
                self.verb = True
                self(event)
                print "\n\n"
                self.verb = False
        ret['Cand_'+self.vname] = isKai
        return ret
    def __del__(self):
        print "From %s: matched %d/%d/%d candidates in %d/%d events" % (self.evlist,self.fcands,self.pcands,self.ncands,self.fevts,self.nevts)

if __name__ == '__main__':
    from sys import argv
    file = ROOT.TFile(argv[1])
    tree = file.Get("treeProducerBphFourMuon")
    tree.vectorTree = True
    class Tester(Module):
        def __init__(self, name):
            Module.__init__(self,name,None)
            #self.sf = CandList("/afs/cern.ch/user/g/gpetrucc/scratch0/cmgprod/CMSSW_5_3_22/src/CMGTools/TTHAnalysis/python/plotter/2012MuOniaforCristina.runlist")
            self.sf = CandList4Mu("/afs/cern.ch/user/g/gpetrucc/scratch0/cmgprod/CMSSW_5_3_22/src/CMGTools/TTHAnalysis/python/plotter/MuOnia2012fourmuonsignal4Cristina.txt",
                                  mtolerance=0.2,xyztolerance=0.0005)
            self.sf2 = CandList4Mu("/afs/cern.ch/user/g/gpetrucc/scratch0/cmgprod/CMSSW_5_3_22/src/CMGTools/TTHAnalysis/python/plotter/MuOnia2012fourmuonsignal4Cristina.txt",
                                  xyztolerance=0.2,mtolerance=1.0,verbose=True)
        def analyze(self,ev):
            ret = self.sf(ev)
            if (ev.run, ev.lumi, ev.evt) in self.sf.candidates and sum(ret['Cand_isKai']) != len(self.sf.candidates[(ev.run,ev.lumi,ev.evt)]) and ev.nCand < 100:
                print "\nrun %6d lumi %4d event %d: cands %d" % (ev.run, ev.lumi, ev.evt, ev.nCand)
                ret = self.sf2(ev)
                cand = Collection(ev,"Cand","nCand")
                for i,c in enumerate(cand):
                    print "\tCand mass %6.3f, kai: %d" % (c.mass,ret['Cand_isKai'][i])
                print "\tKai's candidates: %s" % ( self.sf.candidates[ (ev.run, ev.lumi, ev.evt) ] )
    el = EventLoop([ Tester("tester") ])
    el.loop([tree], maxEvents = -1)

        
