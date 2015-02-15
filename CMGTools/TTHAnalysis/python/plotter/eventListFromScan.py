import ROOT
ROOT.gROOT.SetBatch(True)

from CMGTools.TTHAnalysis.treeReAnalyzer import *
from CMGTools.TTHAnalysis.plotter.tree2yield import CutsFile, scalarToVector

from optparse import OptionParser

parser = OptionParser(usage="usage: %prog [options] rootfile [what] \nrun with --help to get list of options")
parser.add_option("-c", "--cut-file",  dest="cutfile", default=None, type="string", help="Cut file to apply")
parser.add_option("-n", "--maxEvents",  dest="maxEvents", default=1000000000, type="int", help="Max events")
parser.add_option("-f", "--format",   dest="fmt",  default="run:lumi:evt", type="string",  help="Print this format string")
parser.add_option("-t", "--tree",          dest="tree", default='ttHLepTreeProducerTTH', help="Pattern for tree name");
parser.add_option("-F", "--add-friend",    dest="friendTrees",  action="append", default=[], nargs=2, help="Add a friend tree (treename, filename)") 

### CUT-file options
parser.add_option("-S", "--start-at-cut",   dest="startCut",   type="string", help="Run selection starting at the cut matched by this regexp, included.") 
parser.add_option("-U", "--up-to-cut",      dest="upToCut",   type="string", help="Run selection only up to the cut matched by this regexp, included.") 
parser.add_option("-X", "--exclude-cut", dest="cutsToExclude", action="append", default=[], help="Cuts to exclude (regexp matching cut name), can specify multiple times.") 
parser.add_option("-I", "--invert-cut",  dest="cutsToInvert",  action="append", default=[], help="Cuts to invert (regexp matching cut name), can specify multiple times.") 
parser.add_option("-R", "--replace-cut", dest="cutsToReplace", action="append", default=[], nargs=3, help="Cuts to invert (regexp of old cut name, new name, new cut); can specify multiple times.") 
parser.add_option("-A", "--add-cut",     dest="cutsToAdd",     action="append", default=[], nargs=3, help="Cuts to insert (regexp of cut name after which this cut should go, new name, new cut); can specify multiple times.") 
parser.add_option("--s2v", "--scalar2vector",     dest="doS2V",    action="store_true", default=False, help="Do scalar to vector conversion") 
 
(options, args) = parser.parse_args()

cut = ""
if options.cutfile:
    cut = CutsFile(options.cutfile,options).allCuts()
if options.doS2V:
    cut = scalarToVector(cut)

file = ROOT.TFile.Open(args[0])
treename = options.tree
tree = file.Get(treename)
for tf_tree,tf_file in options.friendTrees:
    tf = tree.AddFriend(tf_tree, tf_file),
tree.SetScanField(0)
tree.Scan(options.fmt,cut,"precision=10 colsize=15")

