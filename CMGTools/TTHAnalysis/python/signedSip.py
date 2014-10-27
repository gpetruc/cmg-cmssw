import ROOT

#ROOT.gSystem.Load("libCMGToolsTTHAnalysis")
SignedImpactParameterComputer = ROOT.SignedImpactParameter()

def signedSip3D(lepton, vertex=None):
    if vertex is None:
        vertex = lepton.associatedVertex
    dir   = lepton.jet.momentum()         if hasattr(lepton,'jet')      else lepton.momentum()
    track = lepton.sourcePtr().gsfTrack() if abs(lepton.pdgId()) == 11  else lepton.sourcePtr().track()
    meas = SignedImpactParameterComputer.signedIP3D(track.get(), vertex, dir)
    return meas.significance()


def signedIp3D(lepton, vertex=None):
    if vertex is None:
        vertex = lepton.associatedVertex
    dir   = lepton.jet.momentum()         if hasattr(lepton,'jet')      else lepton.momentum()
    track = lepton.sourcePtr().gsfTrack() if abs(lepton.pdgId()) == 11  else lepton.sourcePtr().track()
    meas = SignedImpactParameterComputer.signedIP3D(track.get(), vertex, dir)
    return meas.value()

def twoTrackChi2(lepton1,lepton2):
    track1 = lepton1.sourcePtr().gsfTrack() if abs(lepton1.pdgId()) == 11  else lepton1.sourcePtr().track()
    track2 = lepton2.sourcePtr().gsfTrack() if abs(lepton2.pdgId()) == 11  else lepton2.sourcePtr().track()
    pair = SignedImpactParameterComputer.twoTrackChi2(track1.get(),track2.get())
    return (pair.first,pair.second)

def fourTrackChi2(lepton1,lepton2,lepton3,lepton4):
    track1 = lepton1.sourcePtr().gsfTrack() if abs(lepton1.pdgId()) == 11  else lepton1.sourcePtr().track()
    track2 = lepton2.sourcePtr().gsfTrack() if abs(lepton2.pdgId()) == 11  else lepton2.sourcePtr().track()
    track3 = lepton3.sourcePtr().gsfTrack() if abs(lepton3.pdgId()) == 11  else lepton3.sourcePtr().track()
    track4 = lepton4.sourcePtr().gsfTrack() if abs(lepton4.pdgId()) == 11  else lepton4.sourcePtr().track()
    pair = SignedImpactParameterComputer.fourTrackChi2(track1.get(),track2.get(),track3.get(),track4.get())
    return (pair.first,pair.second)
