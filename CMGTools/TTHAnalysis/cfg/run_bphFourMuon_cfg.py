##########################################################
##       CONFIGURATION FOR SUSY MULTILEPTON TREES       ##
## skim condition: >= 2 loose leptons, no pt cuts or id ##
##########################################################

import CMGTools.RootTools.fwlite.Config as cfg
from CMGTools.RootTools.fwlite.Config import printComps
from CMGTools.RootTools.RootTools import *

#Load all analyzers
from CMGTools.TTHAnalysis.analyzers.susyCore_modules_cff import * 

# Redefine what I need

# --- LEPTON SKIMMING ---
ttHLepSkim.minLeptons = 1
ttHLepSkim.maxLeptons = 999
#ttHLepSkim.idCut  = ""
#ttHLepSkim.ptCuts = []

# inclusive very loose muon selection
ttHLepAna.inclusive_muon_id  = ""
ttHLepAna.inclusive_muon_pt  = 2
ttHLepAna.inclusive_muon_eta = 2.4
ttHLepAna.inclusive_muon_dxy = 3.0
ttHLepAna.inclusive_muon_dz  = 30.0
# loose muon selection
ttHLepAna.loose_muon_id     = "" #"POG_ID_SoftNew"
ttHLepAna.loose_muon_pt     = 2
ttHLepAna.loose_muon_eta    = 2.4
ttHLepAna.loose_muon_dxy    = 3.0
ttHLepAna.loose_muon_dz     = 30.0
ttHLepAna.loose_muon_relIso = 9e9
# inclusive very loose electron selection
ttHLepAna.inclusive_electron_id  = ""
ttHLepAna.inclusive_electron_pt  = 9e9
ttHLepAna.loose_electron_id     = "POG_MVA_ID_NonTrig"
ttHLepAna.loose_electron_pt     = 9e9


# Event Analyzer for susy multi-lepton (at the moment, it's the TTH one)
bphFourMuAna = cfg.Analyzer(
    'bphFourMuonAnalyzer',
    )

trigger_3mu_jpsi    = "HLT_Dimuon0_Jpsi_Muon_v*"
trigger_3mu_upsilon = "HLT_Dimuon0_Upsilon_Muon_v*"
trigger_3mu_five    = "HLT_TripleMu5_v*"
triggers_3mu_onia = [ trigger_3mu_jpsi, trigger_3mu_upsilon ]
triggers_3mu = [ trigger_3mu_jpsi, trigger_3mu_upsilon, trigger_3mu_five ]
# Tree Producer
treeProducer = cfg.Analyzer(
    'treeProducerBphFourMuon',
    vectorTree = True,
    saveTLorentzVectors = False,  # can set to True to get also the TLorentzVectors, but trees will be bigger
    PDFWeights = PDFWeights,
    triggerBits = {
            'Jpsi'    : [ trigger_3mu_jpsi ],
            'Upsilon' : [ trigger_3mu_upsilon ],
            'TriMu5'  : [ trigger_3mu_five ],
            'PDMuOnia' : [ trigger_3mu_jpsi, trigger_3mu_upsilon ],
        }
    )


#-------- SAMPLES AND TRIGGERS -----------
from CMGTools.TTHAnalysis.samples.samples_8TeV_v517 import * 

for mc in mcSamples+mcSamplesAll:
    mc.triggers = triggers_3mu

for data in dataSamplesMuOnia:
    data.triggers = triggers_3mu_onia
for data in dataSamplesMu:
    data.triggers = trigger_3mu_five
    data.vetoTriggers= triggers_3mu_onia

#selectedComponents = mcSamplesAll + dataSamplesAll
selectedComponents = dataSamplesMu + dataSamplesMuOnia

#-------- SEQUENCE

sequence = cfg.Sequence(susyCoreSequence+[
    bphFourMuAna,
    treeProducer,
    ])


#-------- HOW TO RUN
test = 2
if test==1:
    # test a single component, using a single thread.
    comp = MuOniaAB
    selectedComponents = [comp]
    comp.splitFactor = 60
    eventSelector.toSelect = [
        (190705,50,47481560),
        (190707,169,181561506),
        (190906,231,225619657),
        (191226,89,50899112),
        (191837,33,30725745),
        (193336,242,185660276),
        (193336,501,353871183),
        (193621,142,112638382),
        (193621,1469,1131495656),
        (193621,521,472972796),
    ]
elif test==2:
    # test a single component, using a single thread.
    selectedComponents = [MuOniaAB,MuOniaC,MuOniaD]
    for comp in selectedComponents:
        comp.splitFactor = 60
    eventSelector.toSelect = set()
    kai = open("/afs/cern.ch/user/g/gpetrucc/scratch0/cmgprod/CMSSW_5_3_22/src/CMGTools/TTHAnalysis/python/plotter/2012MuOniaforCristina.runlist");
    for line in kai:
        fields = line.split()
        if len(fields) < 3: continue
        eventSelector.toSelect.add( (int(fields[0]),int(fields[1]),int(fields[2])) )
    print "trying to select ",len(eventSelector.toSelect),"files"




config = cfg.Config( components = selectedComponents,
                     sequence = [eventSelector]+sequence )

printComps(config.components, True)
