import operator 
import itertools
import copy
from math import *

from ROOT import TLorentzVector, TVectorD, TMath

from CMGTools.RootTools.fwlite.Analyzer import Analyzer
from CMGTools.RootTools.fwlite.Event import Event
from CMGTools.RootTools.statistics.Counter import Counter, Counters
from CMGTools.RootTools.fwlite.AutoHandle import AutoHandle
from CMGTools.RootTools.physicsobjects.Lepton import Lepton
from CMGTools.RootTools.physicsobjects.Photon import Photon
from CMGTools.RootTools.physicsobjects.Electron import Electron
from CMGTools.RootTools.physicsobjects.Muon import Muon
from CMGTools.RootTools.physicsobjects.Jet import Jet

from CMGTools.RootTools.utils.DeltaR import * 
from CMGTools.TTHAnalysis.leptonMVA import LeptonMVA
from CMGTools.TTHAnalysis.signedSip import twoTrackChi2, fourTrackChi2
import os

class bphDilepton:
    def __init__(self, idx1, l1, idx2, l2):
        if l1.charge() < 0:
            (idx1,l1,idx2,l2) = (idx2,l2,idx1,l1)
        self.idx1 = idx1
        self.l1 = l1
        self.idx2 = idx2
        self.l2 = l2
        self.p4  = (l1.p4()+l2.p4())
        self.mll = (l1.p4()+l2.p4()).M()
        self.vtx2l = twoTrackChi2(l1,l2)
        self.vtx2l_prob = TMath.Prob(self.vtx2l[0],int(self.vtx2l[1]))
        self.leps = [l1,l2]

class bphFourlepton:
    def __init__(self, dil11, dil12, dil21, dil22):
        self.dil11  = dil11
        self.dil12  = dil12
        self.dil21  = dil21
        self.dil22  = dil22
        self.leps = dil11.leps + dil12.leps 
        self.p4   = (dil11.p4+dil12.p4)
        self.m4l  = (dil11.p4+dil12.p4).M()
        self.vtx4l      = fourTrackChi2(self.leps[0], self.leps[1], self.leps[2], self.leps[3])
        self.vtx4l_prob = TMath.Prob(self.vtx4l[0],int(self.vtx4l[1]))
 
class bphFourMuonAnalyzer( Analyzer ):
    def __init__(self, cfg_ana, cfg_comp, looperName ):
        super(bphFourMuonAnalyzer,self).__init__(cfg_ana,cfg_comp,looperName)

    def declareHandles(self):
        super(bphFourMuonAnalyzer, self).declareHandles()

    def beginLoop(self):
        super(bphFourMuonAnalyzer,self).beginLoop()
        self.counters.addCounter('events')
        count = self.counters.counter('events')
        count.register('all events')

    def makeDileps(self, event):
        event.dileptons = {}
        nlep = len(event.selectedLeptons)
        for i,l1 in enumerate(event.selectedLeptons):
            for j in range(i+1,nlep):
                l2 = event.selectedLeptons[j]    
                if l1.charge() + l2.charge() == 0:
                    dilepton = bphDilepton(i,l1,j,l2) 
                    event.dileptons[ (dilepton.idx1,dilepton.idx2) ] = dilepton 

    def makeFourleps(self, event):
        event.fourleptons = []
        for (i1,i2),dil11 in event.dileptons.iteritems():
            for (i3,i4),dil12 in event.dileptons.iteritems():
                if i1 == i3 or i2 == i4:  continue
                if dil11.mll < dil12.mll: continue 
                dil21 = event.dileptons[(i1,i4)]
                dil22 = event.dileptons[(i3,i2)]
                if dil21.mll < dil22.mll:
                    (dil21,dil22) = (dil22,dil21)
                if dil21.mll > dil11.mll: continue
                event.fourleptons.append( bphFourlepton(dil11,dil12,dil21,dil22) )


    def process(self, iEvent, event):
        self.readCollections( iEvent )
        self.makeDileps( event )
        print "found %d dilepton candidates (possibly partially overlapping)" % len(event.dileptons)
        self.makeFourleps( event )
        print "found %d four-lepton candidates (possibly partially overlapping)" % len(event.fourleptons)
        return len(event.fourleptons) > 0
