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
ttHLepSkim.minLeptons = 4
ttHLepSkim.maxLeptons = 999
#ttHLepSkim.idCut  = ""
#ttHLepSkim.ptCuts = []

# inclusive very loose muon selection
ttHLepAna.inclusive_muon_id  = "POG_ID_SoftNew"
ttHLepAna.inclusive_muon_pt  = 2
ttHLepAna.inclusive_muon_eta = 2.4
ttHLepAna.inclusive_muon_dxy = 3.0
ttHLepAna.inclusive_muon_dz  = 30.0
# loose muon selection
ttHLepAna.loose_muon_id     = "POG_ID_SoftNew"
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
test = 1
if test==1:
    # test a single component, using a single thread.
    comp = MuOniaD
    comp.files = comp.files[:1]
    selectedComponents = [comp]
    comp.splitFactor = 1


config = cfg.Config( components = selectedComponents,
                     sequence = sequence )

printComps(config.components, True)
