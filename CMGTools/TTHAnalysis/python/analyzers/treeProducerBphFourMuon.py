from CMGTools.TTHAnalysis.analyzers.ttHLepTreeProducerNew import *


bphFourLepType = NTupleObjectType("bphFourLep", variables = [
    NTupleVariable("mass", lambda x : x.p4.M()),
    NTupleVariable("vtx4",   lambda x : x.vtx4l_prob),
    NTupleVariable("m11",  lambda x : x.dil11.mll),
    NTupleVariable("m12",  lambda x : x.dil12.mll),
    NTupleVariable("m21",  lambda x : x.dil21.mll),
    NTupleVariable("m22",  lambda x : x.dil22.mll),
    NTupleVariable("vtx11",  lambda x : x.dil11.vtx2l_prob),
    NTupleVariable("vtx12",  lambda x : x.dil12.vtx2l_prob),
    NTupleVariable("vtx21",  lambda x : x.dil21.vtx2l_prob),
    NTupleVariable("vtx22",  lambda x : x.dil22.vtx2l_prob),
    NTupleVariable("y11",  lambda x : x.dil11.y2l),
    NTupleVariable("y12",  lambda x : x.dil12.y2l),
    NTupleVariable("y21",  lambda x : x.dil21.y2l),
    NTupleVariable("y22",  lambda x : x.dil22.y2l),
    NTupleVariable("nT11",  lambda x : x.dil11.nT),
    NTupleVariable("nT12",  lambda x : x.dil12.nT),
    NTupleVariable("nT21",  lambda x : x.dil21.nT),
    NTupleVariable("nT22",  lambda x : x.dil22.nT),
    NTupleVariable("l1h",  lambda x : x.dil11.idx1),
    NTupleVariable("l2h",  lambda x : x.dil11.idx2),
    NTupleVariable("l3h",  lambda x : x.dil12.idx1),
    NTupleVariable("l4h",  lambda x : x.dil12.idx2),
    NTupleVariable("l1l",  lambda x : x.dil21.idx1),
    NTupleVariable("l2l",  lambda x : x.dil21.idx2),
    NTupleVariable("l3l",  lambda x : x.dil22.idx1),
    NTupleVariable("l4l",  lambda x : x.dil22.idx2),
])

bphDiLepType = NTupleObjectType("bphDiLep", variables = [
    NTupleVariable("mass", lambda x : x.mll),
    NTupleVariable("vtx2",   lambda x : x.vtx2l_prob),
    NTupleVariable("y",  lambda x : x.y2l),
    NTupleVariable("l1",  lambda x : x.idx1),
    NTupleVariable("l2",  lambda x : x.idx2),
])

class treeProducerBphFourMuon( ttHLepTreeProducerNew ):

    #-----------------------------------
    # CORE TREE PRODUCER FOR THE SUSY ANALYSES
    # defines the core variables that will be present in the trees of all final states
    #-----------------------------------
    def __init__(self, cfg_ana, cfg_comp, looperName):
        super(treeProducerBphFourMuon,self).__init__(cfg_ana, cfg_comp, looperName)

        ## Declare what we want to fill
        self.globalVariables = [
            NTupleVariable("nVert",  lambda ev: len(ev.goodVertices), int, help="Number of good vertices"),
            ##--------------------------------------------------
            NTupleVariable("nLepTight", lambda ev: sum([l.tightId() for l in ev.selectedLeptons]), int, help="Number of muons passing tight id"),
        ]

        self.globalObjects = {
        }

        self.collections = {
            "selectedLeptons" : NTupleCollection("Lep",  leptonTypeFull, 100, help="Leptons after the preselection"),
            "otherLeptons"    : NTupleCollection("LepOther", leptonTypeFull, 100, help="Leptons after the preselection"),
            "fourleptons"     : NTupleCollection("Cand", bphFourLepType, 2000, help="Four-lepton candidates"),
            "dileptonList"    : NTupleCollection("DiLep",    bphDiLepType,   1000, help="Di-lepton candidates"),
        }

        ## Book the variables, but only if we're called explicitly and not through a base class
        if cfg_ana.name == "treeProducerBphFourMuon":
            self.initDone = True
            self.declareVariables()
