import FWCore.ParameterSet.Config as cms

tauMuSVFit = cms.EDProducer(
    "TauMuWithSVFitProducer",
    diTauSrc = cms.InputTag("cmgTauMuTauPtSel"),
    SVFitVersion =  cms.int32(2), # 1 for 2011 version, 2 for new 2012 (slow) version
    fitAlgo = cms.string('VEGAS'),
    verbose = cms.untracked.bool(False)
    )
