#!/usr/bin/env python
from CMGTools.TTHAnalysis.treeReAnalyzer import *
from glob import glob
import os.path

MODULES = []

#from CMGTools.TTHAnalysis.tools.eventVars_2lss import EventVars2LSS 
#MODULES.append( ('2lss', EventVars2LSS()) )
#from CMGTools.TTHAnalysis.tools.finalMVA_2lss import FinalMVA_2LSS
#MODULES.append( ('2lss_mva', FinalMVA_2LSS()) )
#from CMGTools.TTHAnalysis.tools.finalMVA_3l import FinalMVA_3L
#MODULES.append( ('3l_mva', FinalMVA_3L()) )
#from CMGTools.TTHAnalysis.tools.bbvars import bbVars
#MODULES.append( ('bbvars', bbVars()) )
#from CMGTools.TTHAnalysis.tools.finalMVA_2lss_2 import FinalMVA_2LSS_2
#MODULES.append( ('finalMVA_2lss_2', FinalMVA_2LSS_2()) )
#from CMGTools.TTHAnalysis.tools.ttbarEventReco_2lss import TTEventReco_MC
#MODULES.append( ('ttreco_mc', TTEventReco_MC()) )
#from CMGTools.TTHAnalysis.tools.ttbarEventReco_2lss import TTEventReco
#MODULES.append( ('ttreco', TTEventReco(sortersToUse={"BestGuess":"", "BestBySum4NoTJJb":"_bySum4"})) )
#MODULES.append( ('ttreco', TTEventReco(sortersToUse={"BestGuess":"","ByGuessLL2B":"_byLL"})) )
#from CMGTools.TTHAnalysis.tools.LepMVAFriend import LepMVAFriend
#MODULES.append( ('LepMVAFriend', LepMVAFriend("/afs/cern.ch/user/b/botta/CMGToolsGit/newRecipe/CMSSW_5_3_19/src/CMGTools/TTHAnalysis/macros/leptons")) )
from CMGTools.TTHAnalysis.tools.BphCandList import CandList4Mu
MODULES.append( ('Kai', CandList4Mu("/afs/cern.ch/user/g/gpetrucc/scratch0/cmgprod/CMSSW_5_3_22/src/CMGTools/TTHAnalysis/python/plotter/MuOnia2012fourmuonsignal4Cristina.txt",
                                  mtolerance=0.2,xyztolerance=0.0005)) )
MODULES.append( ('Kai', CandList4Mu("/afs/cern.ch/user/g/gpetrucc/scratch0/cmgprod/CMSSW_5_3_22/src/CMGTools/TTHAnalysis/python/plotter/2012MuOniaforGiovanniTMOST.runlist",
                                  mtolerance=0.2,xyztolerance=0.0005,vname="isKai12")) )

class VariableProducer(Module):
    def __init__(self,name,booker,modules):
        Module.__init__(self,name,booker)
        self._modules = modules
    def beginJob(self):
        self.t = PyTree(self.book("TTree","t","t"))
        self.branches = {}
        for name,mod in self._modules:
            for B in mod.listBranches():
                # don't add the same branch twice
                if B in self.branches: 
                    print "Will not add branch %s twice" % (B,)
                    continue
                self.branches[B] = True
                if type(B) == tuple:
                    if len(B) == 2:
                        self.t.branch(B[0],B[1])
                    elif len(B) == 4:
                        self.t.branch(B[0],B[1],n=B[2],lenVar=B[3])
                else:
                    self.t.branch(B ,"F")
    def analyze(self,event):
        for name,mod in self._modules:
            keyvals = mod(event)
            for B,V in keyvals.iteritems():
                setattr(self.t, B, V)
                setattr(event,  B, V)
        self.t.fill()

import os, itertools

from optparse import OptionParser
parser = OptionParser(usage="%prog [options] <TREE_DIR> <OUT>")
parser.add_option("-d", "--dataset", dest="datasets",  type="string", default=[], action="append", help="Process only this dataset (or dataset if specified multiple times)");
parser.add_option("-c", "--chunk",   dest="chunks",    type="int",    default=[], action="append", help="Process only these chunks (works only if a single dataset is selected with -d)");
parser.add_option("-N", "--events",  dest="chunkSize", type="int",    default=500000, help="Default chunk size when splitting trees");
parser.add_option("-j", "--jobs",    dest="jobs",      type="int",    default=1, help="Use N threads");
parser.add_option("-p", "--pretend", dest="pretend",   action="store_true", default=False, help="Don't run anything");
parser.add_option("-T", "--tree-dir",   dest="treeDir",     type="string", default="sf", help="Directory of the friend tree in the file (default: 'sf')");
parser.add_option("-q", "--queue",   dest="queue",     type="string", default=None, help="Run jobs on lxbatch instead of locally");
parser.add_option("-t", "--tree",    dest="tree",      default='ttHLepTreeProducerTTH', help="Pattern for tree name");
parser.add_option("-V", "--vector",  dest="vectorTree", action="store_true", default=True, help="Input tree is a vector");
parser.add_option("-F", "--add-friend",    dest="friendTrees",  action="append", default=[], nargs=2, help="Add a friend tree (treename, filename). Can use {name}, {cname} patterns in the treename") 
parser.add_option("--FMC", "--add-friend-mc",    dest="friendTreesMC",  action="append", default=[], nargs=2, help="Add a friend tree (treename, filename) to MC only. Can use {name}, {cname} patterns in the treename") 
parser.add_option("--FD", "--add-friend-data",    dest="friendTreesData",  action="append", default=[], nargs=2, help="Add a friend tree (treename, filename) to data trees only. Can use {name}, {cname} patterns in the treename") 
(options, args) = parser.parse_args()

if len(args) != 2 or not os.path.isdir(args[0]) or not os.path.isdir(args[1]): 
    print "Usage: program <TREE_DIR> <OUT>"
    exit()
if len(options.chunks) != 0 and len(options.datasets) != 1:
    print "must specify a single dataset with -d if using -c to select chunks"
    exit()

jobs = []
for D in glob(args[0]+"/*"):
    fname = D+"/"+options.tree+"/"+options.tree+"_tree.root"
    if os.path.exists(fname):
        short = os.path.basename(D)
        if options.datasets != []:
            if short not in options.datasets: continue
        data = ("DoubleMu" in short or "MuEG" in short or "DoubleElectron" in short or "SingleMu" in short)
        f = ROOT.TFile.Open(fname);
        t = f.Get(options.tree)
        entries = t.GetEntries()
        f.Close()
        chunk = options.chunkSize
        if entries < chunk:
            print "  ",os.path.basename(D),("  DATA" if data else "  MC")," single chunk"
            jobs.append((short,fname,"%s/evVarFriend_%s.root" % (args[1],short),data,xrange(entries),-1))
        else:
            nchunk = int(ceil(entries/float(chunk)))
            print "  ",os.path.basename(D),("  DATA" if data else "  MC")," %d chunks" % nchunk
            for i in xrange(nchunk):
                if options.chunks != []:
                    if i not in options.chunks: continue
                r = xrange(int(i*chunk),min(int((i+1)*chunk),entries))
                jobs.append((short,fname,"%s/evVarFriend_%s.chunk%d.root" % (args[1],short,i),data,r,i))
print "\n"
print "I have %d taks to process" % len(jobs)

if options.queue:
    import os, sys
    basecmd = "bsub -q {queue} {dir}/lxbatch_runner.sh {dir} {cmssw} python {self} -N {chunkSize} -T '{tdir}' -t {tree} {data} {output}".format(
                queue = options.queue, dir = os.getcwd(), cmssw = os.environ['CMSSW_BASE'], 
                self=sys.argv[0], chunkSize=options.chunkSize, tdir=options.treeDir, tree=options.tree, data=args[0], output=args[1]
            )
    if options.vectorTree: basecmd += " --vector "
    friendPost =  "".join(["  -F  %s %s " % (fn,ft) for fn,ft in options.friendTrees])
    friendPost += "".join([" --FM %s %s " % (fn,ft) for fn,ft in options.friendTreesMC])
    friendPost += "".join([" --FD %s %s " % (fn,ft) for fn,ft in options.friendTreesData])
    for (name,fin,fout,data,range,chunk) in jobs:
        if chunk != -1:
            print "{base} -d {data} -c {chunk} {post}".format(base=basecmd, data=name, chunk=chunk, post=friendPost)
        else:
            print "{base} -d {data} {post}".format(base=basecmd, data=name, chunk=chunk, post=friendPost)
        
    exit()

maintimer = ROOT.TStopwatch()
def _runIt(myargs):
    (name,fin,fout,data,range,chunk) = myargs
    timer = ROOT.TStopwatch()
    fb = ROOT.TFile(fin)
    tb = fb.Get(options.tree)
    if options.vectorTree:
        tb.vectorTree = True
    else:
        tb.vectorTree = False
    friends = options.friendTrees[:]
    friends += (options.friendTreesData if data else options.friendTreesMC)
    friends_ = [] # to make sure pyroot does not delete them
    for tf_tree,tf_file in friends:
        tf = tb.AddFriend(tf_tree, tf_file.format(name=name, cname=name)),
        friends_.append(tf) # to make sure pyroot does not delete them
    nev = tb.GetEntries()
    if options.pretend:
        print "==== pretending to run %s (%d entries, %s) ====" % (name, nev, fout)
        return (name,(nev,0))
    print "==== %s starting (%d entries) ====" % (name, nev)
    booker = Booker(fout)
    el = EventLoop([ VariableProducer(options.treeDir,booker,MODULES), ])
    el.loop([tb], eventRange=range)
    booker.done()
    fb.Close()
    time = timer.RealTime()
    print "=== %s done (%d entries, %.0f s, %.0f e/s) ====" % ( name, nev, time,(nev/time) )
    return (name,(nev,time))

if options.jobs > 0:
    from multiprocessing import Pool
    pool = Pool(options.jobs)
    ret  = dict(pool.map(_runIt, jobs)) if options.jobs > 0 else dict([_runIt(j) for j in jobs])
else:
    ret = dict(map(_runIt, jobs))
fulltime = maintimer.RealTime()
totev   = sum([ev   for (ev,time) in ret.itervalues()])
tottime = sum([time for (ev,time) in ret.itervalues()])
print "Done %d tasks in %.1f min (%d entries, %.1f min)" % (len(jobs),fulltime/60.,totev,tottime/60.)
