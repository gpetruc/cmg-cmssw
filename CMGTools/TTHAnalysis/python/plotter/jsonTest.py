#!/usr/bin/env python
from math import *
import re, string
from CMGTools.TTHAnalysis.treeReAnalyzer import *
from CMGTools.TTHAnalysis.plotter.tree2yield import CutsFile, scalarToVector

import json,array

def json2map(file):
    jsonmap = {}
    J = json.load(open(file, 'r'))
    nrun, nlumi = 0,0
    for r,lumiranges in J.iteritems():
        lastlumi = max(end for (start,end) in lumiranges)
        lumiscore = array.array('b',[0]*(lastlumi+1))
        nrun += 1
        for (start,end) in lumiranges:
            for lumi in xrange(start,end+1): 
                lumiscore[lumi]=1
                nlumi += 1
        jsonmap[int(r)] = lumiscore
    print "Loaded %d runs, %d lumis from %s " % (nrun,nlumi,file)
    return jsonmap

def root2map(file):
    tfile = ROOT.TFile.Open(file)
    tree = tfile.Get("RLTInfo")
    print "Reading %d entries out of %s" % (tree.GetEntries(),file)
    jsonind = {}
    for e in xrange(tree.GetEntries()):
        tree.GetEntry(e)
        run,lumi = tree.run, tree.lumi
        if run not in jsonind:
            jsonind[run] = [lumi]
        else:
            jsonind[run].append(lumi)
    jsonmap = {}
    for r,lumis in jsonind.iteritems():
        lastlumi = max(lumis)
        lumiscore = array.array('b',[0]*(lastlumi+1))
        for lumi in lumis:
            lumiscore[lumi]=1
        jsonmap[r] = lumiscore
    return jsonmap

where = {
'Golden': json2map('Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt'),
'MuonPhys': json2map('Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON_MuonPhys.txt'),
'MuOniaAB': root2map('/afs/cern.ch/work/b/botta/TREES_53X_060215_bphData/MuOniaAB/JSONAnalyzer/RLTInfo.root'),
'MuOniaC': root2map('/afs/cern.ch/work/b/botta/TREES_53X_060215_bphData/MuOniaC/JSONAnalyzer/RLTInfo.root'),
'MuOniaD': root2map('/afs/cern.ch/work/b/botta/TREES_53X_060215_bphData/MuOniaD/JSONAnalyzer/RLTInfo.root'),
'DoubleMuAB': root2map('/afs/cern.ch/work/b/botta/TREES_53X_060215_bphData/DoubleMuAB/JSONAnalyzer/RLTInfo.root'),
'DoubleMuC': root2map('/afs/cern.ch/work/b/botta/TREES_53X_060215_bphData/DoubleMuC/JSONAnalyzer/RLTInfo.root'),
'DoubleMuD': root2map('/afs/cern.ch/work/b/botta/TREES_53X_060215_bphData/DoubleMuD/JSONAnalyzer/RLTInfo.root'),
}

while True:
    line = sys.stdin.readline()
    if not line or "exit" in line: break
    fields = line.split()
    if len(fields) < 2: continue
    run = int(fields[0])
    lumi = int(fields[1])
    print "Looking for run %d, lumi %d: " % (run,lumi),
    for key,lumimap in where.iteritems():
        if run in lumimap and lumimap[run][lumi]: print key,
    print ""
    print ""
      
