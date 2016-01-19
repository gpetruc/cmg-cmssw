################################
#  use mcEfficiencies.py to make plots of efficiencies and fake rates
################################

#T="/afs/cern.ch/work/b/botta/TREES_74X_100815_MiniIso_25nsMC"
#T="/tmp/gpetrucc/TREES_74X_100815_MiniIso_25nsMC"
#T="/afs/cern.ch/user/g/gpetrucc/w/TREES_74X_25ns_ForLepMVATrain"
#T="/afs/cern.ch/work/p/peruzzi/ra5trees/TREES_74X_191115_MiniIso_v5"
T="/data/p/peruzzi/TREES_74X_070116_MiniIso_forLepMVA"
BCORE=" --s2v --tree treeProducerSusyMultilepton ttH-multilepton/lepton-mca.txt object-studies/lepton-perlep.txt -P ${T} "
#BCORE="${BCORE}  --Fs {P}/0_leptonMVA_v1_NW --Fs {P}/0_leptonMVA_v2pt_NW "
#BCORE="${BCORE}  --Fs {P}/37_gioLepMVATests_v1 "
BCORE="${BCORE}  --Fs {P}/0_lepMVA_forMoriond16_v1 "
PBASE="plots/74X/lepMVA/v3.0"
BG=" -j 8 "
DEN=" LepGood_sip3d < 8  && LepGood_mediumMuonId > 0 "
if [[ "$1" == "-b" ]]; then BG=" & "; shift; fi

what=$1;
case $what in
rocs)
    mkdir -p $PBASE/$what 2> /dev/null && cp /afs/cern.ch/user/g/gpetrucc/php/index.php $PBASE/$what/
    ROCS0="rocCurves.py -p TT_red,TTH_true ${BCORE} ttH-multilepton/lepton-sels.txt ttH-multilepton/lepton-xvars.txt --splitSig 1 ";
    ROCS="${ROCS0}  --logx --grid --xrange 0.007 1.0 --yrange 0.80 1 --max-entries 5000000 "
    echo "( python $ROCS  -o $PBASE/$what/mu_pt_25_inf.root  -R pt20 pt '${DEN} && LepGood_pt > 25'  --sP 'mvaTTH.*,miniRelIso,multiIso' ${BG} )"
    ROCS="${ROCS0}  --logx --grid --xrange 0.007 1.0 --yrange 0.50 1 --max-entries 5000000 "
    echo "( python $ROCS  -o $PBASE/$what/mu_pt_10_25.root   -R pt20 pt '${DEN} && LepGood_pt < 25 && LepGood_pt > 10'  --sP 'mvaTTH.*,miniRelIso,multiIso' ${BG} )"
    ROCS="${ROCS0}  --logx --grid --xrange 0.0007 1.0 --yrange 0.00 1 --max-entries 5000000 "
    echo "( python $ROCS  -o $PBASE/$what/mu_pt_5_10.root    -R pt20 pt '${DEN} && LepGood_pt < 10'  --sP 'mvaTTH.*,miniRelIso,hzz' ${BG} )"

    ROCS="${ROCS0}  --logx --grid --xrange 0.007 1.0 --yrange 0.80 1 --max-entries 5000000 "
    echo "( python $ROCS  -o $PBASE/$what/el_pt_25_inf.root  -R pt20 pt '${DEN} && LepGood_pt > 25' -I mu  --sP 'mvaTTH.*,miniRelIso,multiIso' ${BG} )"
    ROCS="${ROCS0}  --logx --grid --xrange 0.007 1.0 --yrange 0.50 1 --max-entries 5000000 "
    echo "( python $ROCS  -o $PBASE/$what/el_pt_10_25.root   -R pt20 pt '${DEN} && LepGood_pt < 25 && LepGood_pt > 10' -I mu  --sP 'mvaTTH.*,miniRelIso,multiIso' ${BG} )"
    #ROCS="${ROCS0}  --logx --grid --xrange 0.0007 1.0 --yrange 0.00 1 --max-entries 500000 "
    #echo "( python $ROCS  -o $PBASE/$what/el_pt_5_10.root    -R pt20 pt 'LepGood_pt < 10'  -I mu --sP 'mvaTTH.*,miniRelIso,hzz' ${BG} )"
    ;;
effs)
    EFFS="mcEfficiencies.py -p TTH_true,TTH_true_tau,TTH_true_top ${BCORE} ttH-multilepton/lepton-sels.txt ttH-multilepton/lepton-xvars.txt --groupBy process,cut --legend=BR";
    echo "( python $EFFS --yrange 0 1.  -o $PBASE/$what/mu_eta_00_12_comp.root  -R pt20 eta '${DEN} && abs(LepGood_eta)<1.2' --sP pt_fine --sP mvaTTHMoriond16_[TM],multiIso ${BG} )"
    echo "( python $EFFS --yrange 0 1.  -o $PBASE/$what/mu_eta_12_24_comp.root  -R pt20 eta '${DEN} && abs(LepGood_eta)>1.2' --sP pt_fine --sP mvaTTHMoriond16_[TM],multiIso ${BG} )"
    echo "( python $EFFS --yrange 0 1.  -o $PBASE/$what/el_EB_comp.root -I mu -R pt20 eta '${DEN} && abs(LepGood_eta)<1.479' --sP pt_fine --sP mvaTTHMoriond16_[TM],multiIso ${BG} )"
    echo "( python $EFFS --yrange 0 1.  -o $PBASE/$what/el_EE_comp.root -I mu -R pt20 eta '${DEN} && abs(LepGood_eta)>1.479' --sP pt_fine --sP mvaTTHMoriond16_[TM],multiIso ${BG} )"
    echo "( python $EFFS --yrange 0.7 1.  -o $PBASE/$what/mu_pt25_comp.root  -R pt20 eta '${DEN} && LepGood_pt > 25' --sP eta_fine --sP mvaTTHMoriond16_[TM],multiIso ${BG} )"
    echo "( python $EFFS --yrange 0.7 1.  -o $PBASE/$what/el_pt25_comp.root -I mu -R pt20 eta '${DEN} && LepGood_pt > 25' --sP eta_fine --sP mvaTTHMoriond16_[TM],multiIso ${BG} )"
    echo "( python $EFFS --yrange 0.0 1.  -o $PBASE/$what/mu_pt10_25_comp.root  -R pt20 eta '${DEN} && LepGood_pt < 25 && LepGood_pt > 10' --sP eta_fine --sP mvaTTHMoriond16_[TM],multiIso ${BG} )"
    echo "( python $EFFS --yrange 0.0 1.  -o $PBASE/$what/el_pt10_25_comp.root -I mu -R pt20 eta '${DEN} && LepGood_pt < 25 && LepGood_pt > 10' --sP eta_fine --sP mvaTTHMoriond16_[TM],multiIso ${BG} )"
    echo "( python $EFFS --yrange 0 1.  -o $PBASE/$what/mu_eta_00_12_hist_T.root  -R pt20 eta '${DEN} && abs(LepGood_eta)<1.2' --sP pt_fine --sP mvaTTH.*_T ${BG} )"
    echo "( python $EFFS --yrange 0 1.  -o $PBASE/$what/mu_eta_12_24_hist_T.root  -R pt20 eta '${DEN} && abs(LepGood_eta)>1.2' --sP pt_fine --sP mvaTTH.*_T ${BG} )"
    echo "( python $EFFS --yrange 0 1.  -o $PBASE/$what/el_EB_hist_T.root -I mu -R pt20 eta '${DEN} && abs(LepGood_eta)<1.479' --sP pt_fine --sP mvaTTH.*_T ${BG} )"
    echo "( python $EFFS --yrange 0 1.  -o $PBASE/$what/el_EE_hist_T.root -I mu -R pt20 eta '${DEN} && abs(LepGood_eta)>1.479' --sP pt_fine --sP mvaTTH.*_T ${BG} )"
    echo "( python $EFFS --yrange 0.7 1.  -o $PBASE/$what/mu_pt25_hist_T.root  -R pt20 eta '${DEN} && LepGood_pt > 25' --sP eta_fine --sP mvaTT.*_T ${BG} )"
    echo "( python $EFFS --yrange 0.7 1.  -o $PBASE/$what/el_pt25_hist_T.root -I mu -R pt20 eta '${DEN} && LepGood_pt > 25' --sP eta_fine --sP mvaTT.*_T ${BG} )"
    echo "( python $EFFS --yrange 0.0 1.  -o $PBASE/$what/mu_pt10_25_hist_T.root  -R pt20 eta '${DEN} && LepGood_pt < 25 && LepGood_pt > 10' --sP eta_fine --sP mvaTT.*_T ${BG} )"
    echo "( python $EFFS --yrange 0.0 1.  -o $PBASE/$what/el_pt10_25_hist_T.root -I mu -R pt20 eta '${DEN} && LepGood_pt < 25 && LepGood_pt > 10' --sP eta_fine --sP mvaTT.*_T ${BG} )"
    ;;
fakes)
    EFFS="mcEfficiencies.py -p TT_red,TT_bjets ${BCORE} ttH-multilepton/lepton-sels.txt ttH-multilepton/lepton-xvars.txt --groupBy process,cut --legend=TL";
    echo "( python $EFFS --yrange 0 0.12  -o $PBASE/$what/mu_eta_00_12_comp.root  -R pt20 eta '${DEN} && abs(LepGood_eta)<1.2' --sP pt_fcoarse --sP mvaTTHMoriond16_[TM],multiIso ${BG} )"
    echo "( python $EFFS --yrange 0 0.12  -o $PBASE/$what/mu_eta_12_24_comp.root  -R pt20 eta '${DEN} && abs(LepGood_eta)>1.2' --sP pt_fcoarse --sP mvaTTHMoriond16_[TM],multiIso ${BG} )"
    echo "( python $EFFS --yrange 0 0.12  -o $PBASE/$what/mu_pt25_comp.root  -R pt20 eta '${DEN} && LepGood_pt > 25' --sP eta_fcoarse --sP mvaTTHMoriond16_[TM],multiIso ${BG} )"
    EFFS="mcEfficiencies.py -p TT_fake ${BCORE} ttH-multilepton/lepton-sels.txt ttH-multilepton/lepton-xvars.txt --groupBy process,cut --legend=TL";
    echo "( python $EFFS --yrange 0 0.35  -o $PBASE/$what/mu_eta_00_12_comp.root  -R pt20 eta '${DEN} && abs(LepGood_eta)<1.2' --sP pt_fvcoarse --sP mvaTTHMoriond16_[TM],multiIso ${BG} )"
    echo "( python $EFFS --yrange 0 0.35  -o $PBASE/$what/mu_eta_12_24_comp.root  -R pt20 eta '${DEN} && abs(LepGood_eta)>1.2' --sP pt_fvcoarse --sP mvaTTHMoriond16_[TM],multiIso ${BG} )"
    echo "( python $EFFS --yrange 0 0.35  -o $PBASE/$what/mu_pt25_comp.root  -R pt20 eta '${DEN} && LepGood_pt > 25' --sP eta_fvcoarse --sP mvaTTHMoriond16_[TM],multiIso ${BG} )"
    EFFS="mcEfficiencies.py -p TT_red,TT_bjets,TT_fake ${BCORE} ttH-multilepton/lepton-sels.txt ttH-multilepton/lepton-xvars.txt --groupBy process,cut --legend=TL";
    echo "( python $EFFS --yrange 0 0.1  -o $PBASE/$what/el_EB_comp.root -I mu -R pt20 eta '${DEN} && abs(LepGood_eta)<1.479' --sP pt_fcoarse --sP mvaTTHMoriond16_[TM],multiIso ${BG} )"
    echo "( python $EFFS --yrange 0 0.1  -o $PBASE/$what/el_EE_comp.root -I mu -R pt20 eta '${DEN} && abs(LepGood_eta)>1.479' --sP pt_fcoarse --sP mvaTTHMoriond16_[TM],multiIso ${BG} )"
    echo "( python $EFFS --yrange 0 0.1  -o $PBASE/$what/el_pt25_comp.root -I mu -R pt20 eta '${DEN} && LepGood_pt > 25' --sP eta_fcoarse --sP mvaTTHMoriond16_[TM],multiIso ${BG} )"

    EFFS="mcEfficiencies.py -p TT_red,TT_bjets ${BCORE} ttH-multilepton/lepton-sels.txt ttH-multilepton/lepton-xvars.txt --groupBy process,cut --legend=TL";
    echo "( python $EFFS --yrange 0 0.12  -o $PBASE/$what/mu_eta_00_12_hist_T.root  -R pt20 eta '${DEN} && abs(LepGood_eta)<1.2' --sP pt_fcoarse --sP mvaTTH.*_T ${BG} )"
    echo "( python $EFFS --yrange 0 0.12  -o $PBASE/$what/mu_eta_12_24_hist_T.root  -R pt20 eta '${DEN} && abs(LepGood_eta)>1.2' --sP pt_fcoarse --sP mvaTTH.*_T ${BG} )"
    echo "( python $EFFS --yrange 0 0.12  -o $PBASE/$what/mu_pt25_hist_T.root  -R pt20 eta '${DEN} && LepGood_pt > 25' --sP eta_fcoarse --sP mvaTTH.*_T ${BG} )"
    EFFS="mcEfficiencies.py -p TT_fake ${BCORE} ttH-multilepton/lepton-sels.txt ttH-multilepton/lepton-xvars.txt --groupBy process,cut --legend=TL";
    echo "( python $EFFS --yrange 0 0.35  -o $PBASE/$what/mu_eta_00_12_hist_T.root  -R pt20 eta '${DEN} && abs(LepGood_eta)<1.2' --sP pt_fvcoarse --sP mvaTTH.*_T ${BG} )"
    echo "( python $EFFS --yrange 0 0.35  -o $PBASE/$what/mu_eta_12_24_hist_T.root  -R pt20 eta '${DEN} && abs(LepGood_eta)>1.2' --sP pt_fvcoarse --sP mvaTTH.*_T ${BG} )"
    echo "( python $EFFS --yrange 0 0.35  -o $PBASE/$what/mu_pt25_hist_T.root  -R pt20 eta '${DEN} && LepGood_pt > 25' --sP eta_fvcoarse --sP mvaTTH.*_T ${BG} )"
    EFFS="mcEfficiencies.py -p TT_red,TT_bjets,TT_fake ${BCORE} ttH-multilepton/lepton-sels.txt ttH-multilepton/lepton-xvars.txt --groupBy process,cut --legend=TL";
    echo "( python $EFFS --yrange 0 0.1  -o $PBASE/$what/el_EB_hist_T.root -I mu -R pt20 eta '${DEN} && abs(LepGood_eta)<1.479' --sP pt_fcoarse --sP mvaTTH.*_T ${BG} )"
    echo "( python $EFFS --yrange 0 0.1  -o $PBASE/$what/el_EE_hist_T.root -I mu -R pt20 eta '${DEN} && abs(LepGood_eta)>1.479' --sP pt_fcoarse --sP mvaTTH.*_T ${BG} )"
    echo "( python $EFFS --yrange 0 0.1  -o $PBASE/$what/el_pt25_hist_T.root -I mu -R pt20 eta '${DEN} && LepGood_pt > 25' --sP eta_fcoarse --sP mvaTTH.*_T ${BG} )"

    ;;
esac

