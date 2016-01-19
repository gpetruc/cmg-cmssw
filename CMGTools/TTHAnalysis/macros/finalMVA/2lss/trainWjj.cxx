void trainWjj(TString name) {
    TTree *dSig = (TTree*) _file0->Get("sf/t");
    TTree *dBg1 = dSig;
    TFile *fOut = new TFile(name+".root","RECREATE");
    TMVA::Factory *factory = new TMVA::Factory(name, fOut, "!V:!Color:Transformations=I");
    
    TString allvars = ""; 
    factory->AddVariable("Wjj_mass", 'D'); 
    factory->AddVariable("Wjj_pt", 'D');  
    //factory->AddVariable("Wjj_j1pt", 'D');
    factory->AddVariable("Wjj_j2pt", 'D');
    //factory->AddVariable("Wjj_dEta", 'D');
    factory->AddVariable("Wjj_maxbtagCSV := max(0,max(Wjj_j1btagCSV,Wjj_j2btagCSV))", 'D');
    //factory->AddVariable("Wjj_minbtagCSV := max(0,min(Wjj_j1btagCSV,Wjj_j2btagCSV))", 'D');
    //factory->AddVariable("Wjj_maxqgl := max(Wjj_j1qgl,Wjj_j2qgl)", 'D');
    factory->AddVariable("Wjj_minqgl := min(Wjj_j1qgl,Wjj_j2qgl)", 'D');

    TCut presel = "Wjj_mass > 30 && Wjj_mass < 130 && max(Wjj_j1btagCSV,Wjj_j2btagCSV) < 0.99 && min(Wjj_j1btagCSV,Wjj_j2btagCSV) < 0.89";

    double wSig = 1.0, wBkg = 1.0;
    factory->AddSignalTree(dSig, wSig);
    factory->AddBackgroundTree(dBg1, wBkg);

    factory->PrepareTrainingAndTestTree( presel+" Wjj_mciW != -1", presel+" Wjj_mciW == -1", "" );

    factory->BookMethod( TMVA::Types::kLD, "LD", "!H:!V:VarTransform=None" );


    factory->BookMethod( TMVA::Types::kLikelihood, "Likelihood", "!H:!V:!TransformOutput" );

    TString BDTGopt = "!H:!V:NTrees=200:BoostType=Grad:Shrinkage=0.10:!UseBaggedGrad:nCuts=2000:nEventsMin=100:NNodesMax=9:UseNvars=9:MaxDepth=5";
    factory->BookMethod( TMVA::Types::kBDT, "BDTG", BDTGopt);
    //TString BDTGopt2 = "!H:!V:NTrees=200:BoostType=Grad:Shrinkage=0.20:!UseBaggedGrad:nCuts=2000:nEventsMin=100:NNodesMax=9:UseNvars=9:MaxDepth=5:CreateMVAPdfs";
    //factory->BookMethod( TMVA::Types::kBDT, "BDTG2", BDTGopt2);

    factory->TrainAllMethods();
    factory->TestAllMethods();
    factory->EvaluateAllMethods();

    fOut->Close();
}
