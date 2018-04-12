#include "PhysicsTools/Heppy/interface/PATTauDiscriminationByMVAIsolationRun2FWlite.h"
#include "RecoTauTag/RecoTau/interface/PFRecoTauClusterVariables.h"
#include <TFile.h>
#include <cassert>

heppy::PATTauDiscriminationByMVAIsolationRun2FWlite::PATTauDiscriminationByMVAIsolationRun2FWlite(const std::string & fileName, const std::string & mvaName, const std::string &mvaKind)
{
    TFile *f = TFile::Open(fileName.c_str());
    assert(f);
    mvaReader_ = (GBRForest*) f->Get(mvaName.c_str());
    assert(mvaReader_);

    if      ( mvaKind == "oldDMwoLT" ) mvaOpt_ = kOldDMwoLT;
    else if ( mvaKind == "oldDMwLT"  ) mvaOpt_ = kOldDMwLT;
    else if ( mvaKind == "newDMwoLT" ) mvaOpt_ = kNewDMwoLT;
    else if ( mvaKind == "newDMwLT"  ) mvaOpt_ = kNewDMwLT;
    else if ( mvaKind == "DBoldDMwLT"  ) mvaOpt_ = kDBoldDMwLT;
    else if ( mvaKind == "DBnewDMwLT"  ) mvaOpt_ = kDBnewDMwLT;
    else if ( mvaKind == "PWoldDMwLT"  ) mvaOpt_ = kPWoldDMwLT;
    else if ( mvaKind == "PWnewDMwLT"  ) mvaOpt_ = kPWnewDMwLT;
    else if ( mvaKind == "DBoldDMwLTwGJ" ) mvaOpt_ = kDBoldDMwLTwGJ;
    else if ( mvaKind == "DBnewDMwLTwGJ" ) mvaOpt_ = kDBnewDMwLTwGJ;
    else assert(0);

    if      ( mvaOpt_ == kOldDMwoLT || mvaOpt_ == kNewDMwoLT ) mvaInput_ = new float[6];
    else if ( mvaOpt_ == kOldDMwLT  || mvaOpt_ == kNewDMwLT  ) mvaInput_ = new float[12];
    else if ( mvaOpt_ == kDBoldDMwLT || mvaOpt_ == kDBnewDMwLT ||
            mvaOpt_ == kPWoldDMwLT || mvaOpt_ == kPWnewDMwLT ||
            mvaOpt_ == kDBoldDMwLTwGJ || mvaOpt_ == kDBnewDMwLTwGJ) mvaInput_ = new float[23];
    else assert(0);   

    // these are trivial for the moment
    chargedIsoPtSums_ = "chargedIsoPtSum";
    neutralIsoPtSums_ = "neutralIsoPtSum";
    puCorrPtSums_ = "puCorrPtSum";
    footprintCorrection_ = "footprintCorrection";
    photonPtSumOutsideSignalCone_ = "photonPtSumOutsideSignalCone";
}

heppy::PATTauDiscriminationByMVAIsolationRun2FWlite::~PATTauDiscriminationByMVAIsolationRun2FWlite()
{
    delete mvaReader_;
}

float heppy::PATTauDiscriminationByMVAIsolationRun2FWlite::operator()(const pat::Tau *tau) const 
{
    TauIdMVAAuxiliaries clusterVariables_;

    // CV: computation of MVA value requires presence of leading charged hadron
    if ( tau->leadChargedHadrCand().isNull() ) return 0.;

    int tauDecayMode = tau->decayMode();

    if ( ((mvaOpt_ == kOldDMwoLT || mvaOpt_ == kOldDMwLT || mvaOpt_ == kDBoldDMwLT || mvaOpt_ == kPWoldDMwLT || mvaOpt_ == kDBoldDMwLTwGJ)
                && (tauDecayMode == 0 || tauDecayMode == 1 || tauDecayMode == 2 || tauDecayMode == 10))
            ||
            ((mvaOpt_ == kNewDMwoLT || mvaOpt_ == kNewDMwLT || mvaOpt_ == kDBnewDMwLT || mvaOpt_ == kPWnewDMwLT || mvaOpt_ == kDBnewDMwLTwGJ)
             && (tauDecayMode == 0 || tauDecayMode == 1 || tauDecayMode == 2 || tauDecayMode == 5 || tauDecayMode == 6 || tauDecayMode == 10 || tauDecayMode == 11))
       ) {

        float chargedIsoPtSum = tau->tauID(chargedIsoPtSums_);
        float neutralIsoPtSum = tau->tauID(neutralIsoPtSums_);
        float puCorrPtSum     = tau->tauID(puCorrPtSums_);
        float photonPtSumOutsideSignalCone = tau->tauID(photonPtSumOutsideSignalCone_);
        float footprintCorrection = tau->tauID(footprintCorrection_);

        float decayDistX = tau->flightLength().x();
        float decayDistY = tau->flightLength().y();
        float decayDistZ = tau->flightLength().z();
        float decayDistMag = std::sqrt(decayDistX*decayDistX + decayDistY*decayDistY + decayDistZ*decayDistZ);

        // --- The following 5 variables differ slightly between AOD & MiniAOD
        //     because they are recomputed using packedCandidates saved in the tau
        float nPhoton = (float)clusterVariables_.tau_n_photons_total(*tau);
        float ptWeightedDetaStrip = clusterVariables_.tau_pt_weighted_deta_strip(*tau, tauDecayMode);
        float ptWeightedDphiStrip = clusterVariables_.tau_pt_weighted_dphi_strip(*tau, tauDecayMode);
        float ptWeightedDrSignal = clusterVariables_.tau_pt_weighted_dr_signal(*tau, tauDecayMode);
        float ptWeightedDrIsolation = clusterVariables_.tau_pt_weighted_dr_iso(*tau, tauDecayMode);
        // ---
        float leadingTrackChi2 = tau->leadingTrackNormChi2();
        float eRatio = clusterVariables_.tau_Eratio(*tau);

        // Difference between measured and maximally allowed Gottfried-Jackson angle
        float gjAngleDiff = -999;
        if ( tauDecayMode == 10 ) {
            double mTau = 1.77682;
            double mAOne = tau->p4().M();
            double pAOneMag = tau->p();
            double argumentThetaGJmax = (std::pow(mTau,2) - std::pow(mAOne,2) ) / ( 2 * mTau * pAOneMag );
            double argumentThetaGJmeasured = ( tau->p4().px() * decayDistX + tau->p4().py() * decayDistY + tau->p4().pz() * decayDistZ ) / ( pAOneMag * decayDistMag );
            if ( std::abs(argumentThetaGJmax) <= 1. && std::abs(argumentThetaGJmeasured) <= 1. ) {
                double thetaGJmax = std::asin( argumentThetaGJmax );
                double thetaGJmeasured = std::acos( argumentThetaGJmeasured );
                gjAngleDiff = thetaGJmeasured - thetaGJmax;
            }
        }

        if ( mvaOpt_ == kOldDMwoLT || mvaOpt_ == kNewDMwoLT ) {
            mvaInput_[0]  = std::log(std::max(1.f, (float)tau->pt()));
            mvaInput_[1]  = std::abs((float)tau->eta());
            mvaInput_[2]  = std::log(std::max(1.e-2f, chargedIsoPtSum));
            mvaInput_[3]  = std::log(std::max(1.e-2f, neutralIsoPtSum - 0.125f*puCorrPtSum));
            mvaInput_[4]  = std::log(std::max(1.e-2f, puCorrPtSum));
            mvaInput_[5]  = tauDecayMode;
        } else if ( mvaOpt_ == kOldDMwLT || mvaOpt_ == kNewDMwLT  ) {
            mvaInput_[0]  = std::log(std::max(1.f, (float)tau->pt()));
            mvaInput_[1]  = std::abs((float)tau->eta());
            mvaInput_[2]  = std::log(std::max(1.e-2f, chargedIsoPtSum));
            mvaInput_[3]  = std::log(std::max(1.e-2f, neutralIsoPtSum - 0.125f*puCorrPtSum));
            mvaInput_[4]  = std::log(std::max(1.e-2f, puCorrPtSum));
            mvaInput_[5]  = tauDecayMode;
            mvaInput_[6]  = std::copysign(+1.f, tau->dxy());
            mvaInput_[7]  = std::sqrt(std::min(1.f, std::abs(tau->dxy())));
            mvaInput_[8]  = std::min(10.f, std::abs(tau->dxy_Sig()));
            mvaInput_[9]  = ( tau->hasSecondaryVertex() ) ? 1. : 0.;
            mvaInput_[10] = std::sqrt(decayDistMag);
            mvaInput_[11] = std::min(10.f, tau->flightLengthSig());
        } else if ( mvaOpt_ == kDBoldDMwLT || mvaOpt_ == kDBnewDMwLT ) {
            mvaInput_[0]  = std::log(std::max(1.f, (float)tau->pt()));
            mvaInput_[1]  = std::abs((float)tau->eta());
            mvaInput_[2]  = std::log(std::max(1.e-2f, chargedIsoPtSum));
            mvaInput_[3]  = std::log(std::max(1.e-2f, neutralIsoPtSum));
            mvaInput_[4]  = std::log(std::max(1.e-2f, puCorrPtSum));
            mvaInput_[5]  = std::log(std::max(1.e-2f, photonPtSumOutsideSignalCone));
            mvaInput_[6]  = tauDecayMode;
            mvaInput_[7]  = std::min(30.f, nPhoton);
            mvaInput_[8]  = std::min(0.5f, ptWeightedDetaStrip);
            mvaInput_[9]  = std::min(0.5f, ptWeightedDphiStrip);
            mvaInput_[10] = std::min(0.5f, ptWeightedDrSignal);
            mvaInput_[11] = std::min(0.5f, ptWeightedDrIsolation);
            mvaInput_[12] = std::min(100.f, leadingTrackChi2);
            mvaInput_[13] = std::min(1.f, eRatio);
            mvaInput_[14]  = std::copysign(+1.f, tau->dxy());
            mvaInput_[15]  = std::sqrt(std::min(1.f, std::abs(tau->dxy())));
            mvaInput_[16]  = std::min(10.f, std::abs(tau->dxy_Sig()));
            mvaInput_[17]  = std::copysign(+1.f, tau->ip3d());
            mvaInput_[18]  = std::sqrt(std::min(1.f, std::abs(tau->ip3d())));
            mvaInput_[19]  = std::min(10.f, std::abs(tau->ip3d_Sig()));
            mvaInput_[20]  = ( tau->hasSecondaryVertex() ) ? 1. : 0.;
            mvaInput_[21] = std::sqrt(decayDistMag);
            mvaInput_[22] = std::min(10.f, tau->flightLengthSig());
        } else if ( mvaOpt_ == kPWoldDMwLT || mvaOpt_ == kPWnewDMwLT ) {
            mvaInput_[0]  = std::log(std::max(1.f, (float)tau->pt()));
            mvaInput_[1]  = std::abs((float)tau->eta());
            mvaInput_[2]  = std::log(std::max(1.e-2f, chargedIsoPtSum));
            mvaInput_[3]  = std::log(std::max(1.e-2f, neutralIsoPtSum));
            mvaInput_[4]  = std::log(std::max(1.e-2f, footprintCorrection));
            mvaInput_[5]  = std::log(std::max(1.e-2f, photonPtSumOutsideSignalCone));
            mvaInput_[6]  = tauDecayMode;
            mvaInput_[7]  = std::min(30.f, nPhoton);
            mvaInput_[8]  = std::min(0.5f, ptWeightedDetaStrip);
            mvaInput_[9]  = std::min(0.5f, ptWeightedDphiStrip);
            mvaInput_[10] = std::min(0.5f, ptWeightedDrSignal);
            mvaInput_[11] = std::min(0.5f, ptWeightedDrIsolation);
            mvaInput_[12] = std::min(100.f, leadingTrackChi2);
            mvaInput_[13] = std::min(1.f, eRatio);
            mvaInput_[14]  = std::copysign(+1.f, tau->dxy());
            mvaInput_[15]  = std::sqrt(std::min(1.f, std::abs(tau->dxy())));
            mvaInput_[16]  = std::min(10.f, std::abs(tau->dxy_Sig()));
            mvaInput_[17]  = std::copysign(+1.f, tau->ip3d());
            mvaInput_[18]  = std::sqrt(std::min(1.f, std::abs(tau->ip3d())));
            mvaInput_[19]  = std::min(10.f, std::abs(tau->ip3d_Sig()));
            mvaInput_[20]  = ( tau->hasSecondaryVertex() ) ? 1. : 0.;
            mvaInput_[21] = std::sqrt(decayDistMag);
            mvaInput_[22] = std::min(10.f, tau->flightLengthSig());
        } else if ( mvaOpt_ == kDBoldDMwLTwGJ || mvaOpt_ == kDBnewDMwLTwGJ ) {
            mvaInput_[0]  = std::log(std::max(1.f, (float)tau->pt()));
            mvaInput_[1]  = std::abs((float)tau->eta());
            mvaInput_[2]  = std::log(std::max(1.e-2f, chargedIsoPtSum));
            mvaInput_[3]  = std::log(std::max(1.e-2f, neutralIsoPtSum));
            mvaInput_[4]  = std::log(std::max(1.e-2f, puCorrPtSum));
            mvaInput_[5]  = std::log(std::max(1.e-2f, photonPtSumOutsideSignalCone));
            mvaInput_[6]  = tauDecayMode;
            mvaInput_[7]  = std::min(30.f, nPhoton);
            mvaInput_[8]  = std::min(0.5f, ptWeightedDetaStrip);
            mvaInput_[9]  = std::min(0.5f, ptWeightedDphiStrip);
            mvaInput_[10] = std::min(0.5f, ptWeightedDrSignal);
            mvaInput_[11] = std::min(0.5f, ptWeightedDrIsolation);
            mvaInput_[12] = std::min(1.f, eRatio);
            mvaInput_[13]  = std::copysign(+1.f, tau->dxy());
            mvaInput_[14]  = std::sqrt(std::min(1.f, std::abs(tau->dxy())));
            mvaInput_[15]  = std::min(10.f, std::abs(tau->dxy_Sig()));
            mvaInput_[16]  = std::copysign(+1.f, tau->ip3d());
            mvaInput_[17]  = std::sqrt(std::min(1.f, std::abs(tau->ip3d())));
            mvaInput_[18]  = std::min(10.f, std::abs(tau->ip3d_Sig()));
            mvaInput_[19]  = ( tau->hasSecondaryVertex() ) ? 1. : 0.;
            mvaInput_[20] = std::sqrt(decayDistMag);
            mvaInput_[21] = std::min(10.f, tau->flightLengthSig());
            mvaInput_[22] = std::max(-1.f, gjAngleDiff);
        }

        return mvaReader_->GetClassifier(mvaInput_);
    } else {
        return -1.;
    }
}
