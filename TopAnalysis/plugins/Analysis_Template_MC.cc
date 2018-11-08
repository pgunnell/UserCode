#include <iostream>
#include <sstream>
#include <istream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <functional>
#include <vector>
#include <cassert>
#include "TMath.h"
#include "TRandom.h"
#include "TFile.h"
#include "TTree.h"

#include "Analysis_Template_MC.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"

//#include "KKousour/TopAnalysis/test/ResolutionFits/RooUnfold-1.1.1/src/RooUnfold.h"
//#include "KKousour/TopAnalysis/test/ResolutionFits/RooUnfold-1.1.1/src/RooUnfoldResponse.h"

using namespace std;

// reader.load(...)     // for FLAV_C
// reader.load(...)     // for FLAV_UDSG


//---------------------------- Constructor Of The Class TriggerTurnOn -------------------------- //
Analysis_Template_MC::Analysis_Template_MC(edm::ParameterSet const& cfg)
{
  mFileName       = cfg.getParameter<std::vector<std::string>>     ("filename");
  mTreeName       = cfg.getParameter<std::string>               ("treename");
  mDirName        = cfg.getParameter<std::string>               ("dirname");
  mMCtype         = cfg.getParameter<std::string>               ("MCtype");
  
  mIsMCarlo       = cfg.getUntrackedParameter<bool>             ("isMCarlo");

  musePDFWeights       = cfg.getUntrackedParameter<bool>             ("usePDFWeights");
  museScaleWeights       = cfg.getUntrackedParameter<bool>             ("useScaleWeights");

  mnumbPDFWeights       = cfg.getUntrackedParameter<int>             ("numbPDFWeight");
  mnumbScaleWeights       = cfg.getUntrackedParameter<int>             ("numbScaleWeight");
  
  mCrossSection   = cfg.getUntrackedParameter<double>             ("CrossSection");
  mIntLumi        = cfg.getUntrackedParameter<double>             ("IntLumi");

  mShift          = cfg.getUntrackedParameter<double>             ("ShiftJEC");
  
  mTriggers       = cfg.getUntrackedParameter<std::vector<std::string>>     ("Triggers");


}

//------------------------------ Declaration Of The Function beginjob() ------------------------//
void Analysis_Template_MC::beginJob()
 {

     //TBranch *branch = mTree->GetBranch("events");

     NEvents=0;
     Neventsdetselection=0;
     Neventsgenselection=0;
     Neventspartonselection=0;
     dijetsafter=0;
     dijetsbefore=0;

     double Massbinning[401];

     for(int i = 0; i < 401; i++) {
       Massbinning[i]=20+2*i;
     }

     int Massbins=400;

     double SubJetMassbinning[22] = {0.001, 30, 50, 70, 75,80,85,90, 100,110,130, 150, 160, 170, 180, 190, 200, 300, 400, 500, 800, 1200};
     
     int SubJetMassbins=21;

     double Ptbinning[81] = {0, 1, 5, 6, 8, 10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,507, 548, 592, 638, 686, 737, 790, 846, 905, 967,1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000,2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832,4037, 4252, 4477, 4713, 4961, 5220, 5492, 5777, 6076, 6389, 6717, 7000};
     
     int Ptbins=80;

     double PtbinningLarge[6] = {395, 548, 737, 967, 1248, 4252};
     int PtbinsLarge=5;

     double PtbinningV2[11] = {395, 471.5, 548, 642.5, 737, 852, 967, 1107.5, 1248, 1588, 4252};
     int PtbinsV2=10;

     double DeltaPhibinning[8] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 2.75, TMath::Pi()};
     int DeltaPhibins=7;
     
     double DeltaPhibinningLarge[15] = {0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.625, 2.75, 2.875, TMath::Pi()};
     int DeltaPhibinsLarge=14;


     //double DeltaPhibinning[6] = {0.0,1.0, 1.5,2.0, 2.5, TMath::Pi()};
     //int DeltaPhibins=5;

     //------------------ Histogram Booking --------------------------- //
     num_of_Vtx     = fs->make<TH1F>("num_of_Vtx","num_of_Vtx",100,0.,100.); num_of_Vtx->Sumw2();
     PileUpInteractions     = fs->make<TH1F>("PileUpInteractions","PileUpInteractions",100,0.,100.); PileUpInteractions->Sumw2();
     num_of_Jets     = fs->make<TH1F>("num_of_Jets","num_of_Jets",100,0.,100.); num_of_Jets->Sumw2();
     num_of_SubJets     = fs->make<TH1F>("num_of_SubJets","num_of_SubJets",100,0.,100.); num_of_SubJets->Sumw2();
     num_of_SubBJets     = fs->make<TH1F>("num_of_SubBJets","num_of_SubBJets",100,0.,100.); num_of_SubBJets->Sumw2();
     num_of_BJets     = fs->make<TH1F>("num_of_BJets","num_of_BJets",100,0.,100.); num_of_BJets->Sumw2();
     num_of_GenJets     = fs->make<TH1F>("num_of_GenJets","num_of_GenJets",100,0.,100.); num_of_GenJets->Sumw2();
     num_of_Leptons     = fs->make<TH1F>("num_of_Leptons","num_of_Leptons",100,0.,100.); num_of_Leptons->Sumw2();
     num_of_GenLeptons     = fs->make<TH1F>("num_of_GenLeptons","num_of_GenLeptons",100,0.,100.); num_of_GenLeptons->Sumw2();

     num_of_VtxSelection     = fs->make<TH1F>("num_of_VtxSelection","num_of_VtxSelection",100,0.,100.); num_of_VtxSelection->Sumw2();

     ptDETJet  = fs->make<TH1F>("ptDETJet","ptDETJet",Ptbins,Ptbinning); ptDETJet->Sumw2();
     yDETJet = fs->make<TH1F>("yDETJet","yDETJet",60,-3.,3.); yDETJet->Sumw2();
     phiDETJet = fs->make<TH1F>("phiDETJet","phiDETJet",60, -TMath::Pi(),TMath::Pi()); phiDETJet->Sumw2();

     ChfDETJet = fs->make<TH1F>("ChfDETJet","ChfDETJet",100,0.,1.); ChfDETJet->Sumw2();
     ElfDETJet = fs->make<TH1F>("ElfDETJet","ElfDETJet",100,0.,1.); ElfDETJet->Sumw2(); 
     PhfDETJet = fs->make<TH1F>("PhfDETJet","PhfDETJet",100,0.,1.); PhfDETJet->Sumw2();
     NhfDETJet = fs->make<TH1F>("NhfDETJet","NhfDETJet",100,0.,1.); NhfDETJet->Sumw2();
     MufDETJet = fs->make<TH1F>("MufDETJet","MufDETJet",100,0.,1.); MufDETJet->Sumw2();

     Tau1DETJet = fs->make<TH1F>("Tau1DETJet","Tau1DETJet",100,0.,1.); Tau1DETJet->Sumw2();
     Tau2DETJet = fs->make<TH1F>("Tau2DETJet","Tau2DETJet",100,0.,1.); Tau2DETJet->Sumw2();
     Tau3DETJet = fs->make<TH1F>("Tau3DETJet","Tau3DETJet",100,0.,1.); Tau3DETJet->Sumw2();
     Tau31DETJet = fs->make<TH1F>("Tau31DETJet","Tau31DETJet",100,0.,1.); Tau31DETJet->Sumw2();
     Tau32DETJet = fs->make<TH1F>("Tau32DETJet","Tau32DETJet",100,0.,1.); Tau32DETJet->Sumw2();

     Tau31GENJet = fs->make<TH1F>("Tau31GENJet","Tau31GENJet",100,0.,1.); Tau31GENJet->Sumw2();
     Tau32GENJet = fs->make<TH1F>("Tau32GENJet","Tau32GENJet",100,0.,1.); Tau32GENJet->Sumw2();

     EnergyDETJet= fs->make<TH1F>("EnergyDETJet","EnergyDETJet",100,0.,600.); EnergyDETJet->Sumw2();
     MassDETJet= fs->make<TH1F>("MassDETJet","MassDETJet",200,0.,2000.); MassDETJet->Sumw2();
     MassSoftDropDETJet= fs->make<TH1F>("MassSoftDropDETJet","MassSoftDropDETJet",200,0.,2000.); MassSoftDropDETJet->Sumw2();
     BtagDETJet= fs->make<TH1F>("BtagDETJet","BtagDETJet",100,0.,1.); BtagDETJet->Sumw2();
     
     PtDETSubJet0  = fs->make<TH1F>("PtDETSubJet0","PtDETSubJet0",Ptbins,Ptbinning); PtDETSubJet0->Sumw2();
     PtDETSubJet1  = fs->make<TH1F>("PtDETSubJet1","PtDETSubJet1",Ptbins,Ptbinning); PtDETSubJet1->Sumw2();
     EtaDETSubJet0 = fs->make<TH1F>("EtaDETSubJet0","EtaDETSubJet0",60,-3.,3.); EtaDETSubJet0->Sumw2();
     EtaDETSubJet1 = fs->make<TH1F>("EtaDETSubJet1","EtaDETSubJet1",60,-3.,3.); EtaDETSubJet1->Sumw2();
     MassDETSubJet0 = fs->make<TH1F>("MassDETSubJet0","MassDETSubJet0",200,0.,2000.); MassDETSubJet0->Sumw2();
     MassDETSubJet1 = fs->make<TH1F>("MassDETSubJet1","MassDETSubJet1",200,0.,2000.); MassDETSubJet1->Sumw2();
     BtagDETSubJet0 = fs->make<TH1F>("BtagDETSubJet0","BtagDETSubJet0",60,0.,1.); BtagDETSubJet0->Sumw2();
     BtagDETSubJet1 = fs->make<TH1F>("BtagDETSubJet1","BtagDETSubJet1",60,0.,1.); BtagDETSubJet1->Sumw2();
     FlavourDETSubJet0 = fs->make<TH1F>("FlavourDETSubJet0","FlavourDETSubJet0",23,-0.5,22.5); FlavourDETSubJet0->Sumw2();
     FlavourDETSubJet1 = fs->make<TH1F>("FlavourDETSubJet1","FlavourDETSubJet1",23,-0.5,22.5); FlavourDETSubJet1->Sumw2();

     PhiDETSubJet0 = fs->make<TH1F>("PhiDETSubJet0","PhiDETSubJet0",60,-TMath::Pi(),TMath::Pi()); PhiDETSubJet0->Sumw2();
     PhiDETSubJet1 = fs->make<TH1F>("PhiDETSubJet1","PhiDETSubJet1",60,-TMath::Pi(),TMath::Pi()); PhiDETSubJet1->Sumw2();

     TaggedAndBJetDETSubJet0NUM = fs->make<TH1F>("TaggedAndBJetDETSubJet0NUM","TaggedAndBJetDETSubJet0NUM",Ptbins,Ptbinning); TaggedAndBJetDETSubJet0NUM->Sumw2();
     TaggedAndBJetDETSubJet1NUM = fs->make<TH1F>("TaggedAndBJetDETSubJet1NUM","TaggedAndBJetDETSubJet1NUM",Ptbins,Ptbinning); TaggedAndBJetDETSubJet1NUM->Sumw2();

     BJetAndTaggedDETSubJet0NUM = fs->make<TH1F>("BJetAndTaggedDETSubJet0NUM","BJetAndTaggedDETSubJet0NUM",Ptbins,Ptbinning); BJetAndTaggedDETSubJet0NUM->Sumw2();
     BJetAndTaggedDETSubJet1NUM = fs->make<TH1F>("BJetAndTaggedDETSubJet1NUM","BJetAndTaggedDETSubJet1NUM",Ptbins,Ptbinning); BJetAndTaggedDETSubJet1NUM->Sumw2();

     TaggedAndBJetDETSubJet0DEN = fs->make<TH1F>("TaggedAndBJetDETSubJet0DEN","TaggedAndBJetDETSubJet0DEN",Ptbins,Ptbinning); TaggedAndBJetDETSubJet0DEN->Sumw2();
     TaggedAndBJetDETSubJet1DEN = fs->make<TH1F>("TaggedAndBJetDETSubJet1DEN","TaggedAndBJetDETSubJet1DEN",Ptbins,Ptbinning); TaggedAndBJetDETSubJet1DEN->Sumw2();

     BJetAndTaggedDETSubJet0DEN = fs->make<TH1F>("BJetAndTaggedDETSubJet0DEN","BJetAndTaggedDETSubJet0DEN",Ptbins,Ptbinning); BJetAndTaggedDETSubJet0DEN->Sumw2();
     BJetAndTaggedDETSubJet1DEN = fs->make<TH1F>("BJetAndTaggedDETSubJet1DEN","BJetAndTaggedDETSubJet1DEN",Ptbins,Ptbinning); BJetAndTaggedDETSubJet1DEN->Sumw2();

     ptGENJet  = fs->make<TH1F>("ptGENJet","ptGENJet",Ptbins,Ptbinning); ptGENJet->Sumw2();
     ptGENJetLeading  = fs->make<TH1F>("ptGENJetLeading","ptGENJetLeading",Ptbins,Ptbinning); ptGENJetLeading->Sumw2();
     yGENJet = fs->make<TH1F>("yGENJet","yGENJet",60,-3.,3.); yGENJet->Sumw2();
     phiGENJet = fs->make<TH1F>("phiGENJet","phiGENJet",60, -TMath::Pi(),TMath::Pi()); phiGENJet->Sumw2();
     EnergyGENJet= fs->make<TH1F>("EnergyGENJet","EnergyGENJet",100,0.,600.); EnergyGENJet->Sumw2();
     MassGENJet= fs->make<TH1F>("MassGENJet","MassGENJet",200,0.,2000.); MassGENJet->Sumw2();
     MassSoftDropGENJet= fs->make<TH1F>("MassSoftDropGENJet","MassSoftDropGENJet",200,0.,2000.); MassSoftDropGENJet->Sumw2();

     JetBTaggedSubjets0Mass_FirstSubjet= fs->make<TH1F>("JetBTaggedSubjets0Mass_FirstSubjet","JetBTaggedSubjets0Mass_FirstSubjet",200,0.,2000.); JetBTaggedSubjets0Mass_FirstSubjet->Sumw2();
     JetBTaggedSubjets1Mass_FirstSubjet = fs->make<TH1F>("JetBTaggedSubjets1Mass_FirstSubjet","JetBTaggedSubjets1Mass_FirstSubjet",200,0.,2000.); JetBTaggedSubjets1Mass_FirstSubjet->Sumw2();
     JetBTaggedSubjets0Mass_SecondSubjet= fs->make<TH1F>("JetBTaggedSubjets0Mass_SecondSubjet","JetBTaggedSubjets0Mass_SecondSubjet",200,0.,2000.); JetBTaggedSubjets0Mass_SecondSubjet->Sumw2();
     JetBTaggedSubjets1Mass_SecondSubjet = fs->make<TH1F>("JetBTaggedSubjets1Mass_SecondSubjet","JetBTaggedSubjets1Mass_SecondSubjet",200,0.,2000.); JetBTaggedSubjets1Mass_SecondSubjet->Sumw2();

     BJetGenFlag = fs->make<TH1F>("BJetGenFlag","BJetGenFlag",2,-0.5,1.5); BJetGenFlag->Sumw2();
     WJetGenFlag = fs->make<TH1F>("WJetGenFlag","WJetGenFlag",2,-0.5,1.5); WJetGenFlag->Sumw2();
     metGEN = fs->make<TH1F>("metGEN","metGEN",100,0.,300); metGEN->Sumw2();
     metSigGEN = fs->make<TH1F>("metSigGEN","metSigGEN",100,0.,1.0); metSigGEN->Sumw2();
     mvaGEN = fs->make<TH1F>("mvaGEN","mvaGEN",100,0.,1.0); mvaGEN->Sumw2();

     metDET = fs->make<TH1F>("metDET","metDET",100,0.,300); metDET->Sumw2();
     metSigDET = fs->make<TH1F>("metSigDET","metSigDET",100,0.,1.0); metSigDET->Sumw2();
     mvaDET = fs->make<TH1F>("mvaDET","mvaDET",100,0.,1.0); mvaDET->Sumw2();

     ptTopParton  = fs->make<TH1F>("ptTopParton","ptTopParton",Ptbins,Ptbinning); ptTopParton->Sumw2();
     yTopParton = fs->make<TH1F>("yTopParton","yTopParton",60,-3.,3.); yTopParton->Sumw2();
     phiTopParton = fs->make<TH1F>("phiTopParton","phiTopParton",60, -TMath::Pi(),TMath::Pi()); phiTopParton->Sumw2();
     EnergyTopParton= fs->make<TH1F>("EnergyTopParton","EnergyTopParton",100,0.,600.); EnergyTopParton->Sumw2();

     ptTopParton_Leading= fs->make<TH1F>("ptTopParton_Leading","ptTopParton_Leading",Ptbins,Ptbinning); ptTopParton_Leading->Sumw2();
     ptTopParton_SubLeading= fs->make<TH1F>("ptTopParton_SubLeading","ptTopParton_SubLeading",Ptbins,Ptbinning); ptTopParton_SubLeading->Sumw2();
     ptTopParton_Leading_AfterMatching= fs->make<TH1F>("ptTopParton_Leading_AfterMatching","ptTopParton_Leading_AfterMatching",Ptbins,Ptbinning); ptTopParton_Leading_AfterMatching->Sumw2();
     ptTopParton_SubLeading_AfterMatching= fs->make<TH1F>("ptTopParton_SubLeading_AfterMatching","ptTopParton_SubLeading_AfterMatching",Ptbins,Ptbinning); ptTopParton_SubLeading_AfterMatching->Sumw2();
     ptTopParton_Leading_AfterMatchingTagging= fs->make<TH1F>("ptTopParton_Leading_AfterMatchingTagging","ptTopParton_Leading_AfterMatchingTagging",Ptbins,Ptbinning); ptTopParton_Leading_AfterMatchingTagging->Sumw2();
     ptTopParton_SubLeading_AfterMatchingTagging= fs->make<TH1F>("ptTopParton_SubLeading_AfterMatchingTagging","ptTopParton_SubLeading_AfterMatchingTagging",Ptbins,Ptbinning); ptTopParton_SubLeading_AfterMatchingTagging->Sumw2();

     ptTopParton_Leading_SubSubLeading = fs->make<TH1F>("ptTopParton_Leading_SubSubLeading","ptTopParton_Leading_SubSubLeading",Ptbins,Ptbinning); ptTopParton_Leading_SubSubLeading->Sumw2();
     ptTopParton_SubLeading_SubSubLeading = fs->make<TH1F>("ptTopParton_SubLeading_SubSubLeading","ptTopParton_SubLeading_SubSubLeading",Ptbins,Ptbinning); ptTopParton_SubLeading_SubSubLeading->Sumw2();

     ptTopPartonLeading  = fs->make<TH1F>("ptTopPartonLeading","ptTopPartonLeading",Ptbins,Ptbinning); ptTopPartonLeading->Sumw2();
     ptTopPartonSubLeading  = fs->make<TH1F>("ptTopPartonSubLeading","ptTopPartonSubLeading",Ptbins,Ptbinning); ptTopPartonSubLeading->Sumw2();
     ptTopPartonMatchedLeading  = fs->make<TH1F>("ptTopPartonMatchedLeading","ptTopPartonMatchedLeading",Ptbins,Ptbinning); ptTopPartonMatchedLeading->Sumw2();
     ptTopPartonMatchedSubLeading  = fs->make<TH1F>("ptTopPartonMatchedSubLeading","ptTopPartonMatchedSubLeading",Ptbins,Ptbinning); ptTopPartonMatchedSubLeading->Sumw2();
     ptTopPartonMatchedWBGenJetLeading  = fs->make<TH1F>("ptTopPartonMatchedWBGenJetLeading","ptTopPartonMatchedWBGenJetLeading",Ptbins,Ptbinning); ptTopPartonMatchedWBGenJetLeading->Sumw2();
     ptTopPartonMatchedWBGenJetSubLeading  = fs->make<TH1F>("ptTopPartonMatchedWBGenJetSubLeading","ptTopPartonMatchedWBGenJetSubLeading",Ptbins,Ptbinning); ptTopPartonMatchedWBGenJetSubLeading->Sumw2();

     ptTopPartonMatchedWSubjetMassBGenJetLeading  = fs->make<TH1F>("ptTopPartonMatchedWSubjetMassBGenJetLeading","ptTopPartonMatchedWSubjetMassBGenJetLeading",Ptbins,Ptbinning); ptTopPartonMatchedWSubjetMassBGenJetLeading->Sumw2();
     ptTopPartonMatchedWSubjetMassBGenJetSubLeading  = fs->make<TH1F>("ptTopPartonMatchedWSubjetMassBGenJetSubLeading","ptTopPartonMatchedWSubjetMassBGenJetSubLeading",Ptbins,Ptbinning); ptTopPartonMatchedWSubjetMassBGenJetSubLeading->Sumw2();

     DeltaR_GenJetPartonLead = fs->make<TH1F>("DeltaR_GenJetPartonLead","DeltaR_GenJetPartonLead",60, 0,2*TMath::Pi()); DeltaR_GenJetPartonLead->Sumw2();
     DeltaR_GenJetPartonSubLead = fs->make<TH1F>("DeltaR_GenJetPartonSubLead","DeltaR_GenJetPartonSubLead",60, 0,2*TMath::Pi()); DeltaR_GenJetPartonSubLead->Sumw2();
     DeltaR_GenJetW = fs->make<TH1F>("DeltaR_GenJetW","DeltaR_GenJetW",60, 0,2*TMath::Pi()); DeltaR_GenJetW->Sumw2();
     DeltaR_DetJetParton = fs->make<TH1F>("DeltaR_DetJetParton","DeltaR_DetJetParton",60, 0,2*TMath::Pi()); DeltaR_DetJetParton->Sumw2();     
     DeltaR_GenJetRecoJet= fs->make<TH1F>("DeltaR_GenJetRecoJet","DeltaR_GenJetRecoJet",60, 0,2*TMath::Pi()); DeltaR_GenJetRecoJet->Sumw2();     

     ResolutionPtPartonGenJet = fs->make<TH1F>("ResolutionPtPartonGenJet","ResolutionPtPartonGenJet",60,-3.,3.); ResolutionPtPartonGenJet->Sumw2();
     ResolutionPtPartonDetJet = fs->make<TH1F>("ResolutionPtPartonDetJet","ResolutionPtPartonDetJet",60,-3.,3.); ResolutionPtPartonDetJet->Sumw2();
     ResolutionPtRecoGenJet = fs->make<TH1F>("ResolutionPtRecoGenJet","ResolutionPtRecoGenJet",60,-3.,3.); ResolutionPtRecoGenJet->Sumw2();

     ResolutionMassRecoGenJet = fs->make<TH1F>("ResolutionMassRecoGenJet","ResolutionMassRecoGenJet",60,-3.,3.); ResolutionMassRecoGenJet->Sumw2();
     ResolutionSDMassRecoGenJet = fs->make<TH1F>("ResolutionSDMassRecoGenJet","ResolutionSDMassRecoGenJet",60,-3.,3.); ResolutionSDMassRecoGenJet->Sumw2();
     ResolutionSubJetMassRecoGenJet = fs->make<TH1F>("ResolutionSubJetMassRecoGenJet","ResolutionSubJetMassRecoGenJet",60,-3.,3.); ResolutionSubJetMassRecoGenJet->Sumw2();

     ResponsePtRecoGenJet = fs->make<TH2F>("ResponsePtRecoGenJet","ResponsePtRecoGenJet",PtbinsLarge,PtbinningLarge,PtbinsLarge,PtbinningLarge); ResponsePtRecoGenJet->Sumw2();
     ResponseRapidityRecoGenJet = fs->make<TH2F>("ResponseRapidityRecoGenJet","ResponseRapidityRecoGenJet",60,-3.,3.,60,-3.,3.); ResponseRapidityRecoGenJet->Sumw2();
     ResponsePtPartonGenJet = fs->make<TH2F>("ResponsePtPartonGenJet","ResponsePtPartonGenJet",PtbinsLarge,PtbinningLarge,PtbinsLarge,PtbinningLarge); ResponsePtPartonGenJet->Sumw2();

     ResolutionPhiPartonGenJet = fs->make<TH1F>("ResolutionPhiPartonGenJet","ResolutionPhiPartonGenJet",60,-3.,3.); ResolutionPhiPartonGenJet->Sumw2();
     ResolutionPhiPartonDetJet = fs->make<TH1F>("ResolutionPhiPartonDetJet","ResolutionPhiPartonDetJet",60,-3.,3.); ResolutionPhiPartonDetJet->Sumw2();
     ResolutionPhiRecoGenJet = fs->make<TH1F>("ResolutionPhiRecoGenJet","ResolutionPhiRecoGenJet",60,-3.,3.); ResolutionPhiRecoGenJet->Sumw2();
     ResolutionRapidityRecoGenJet = fs->make<TH1F>("ResolutionRapidityRecoGenJet","ResolutionRapidityRecoGenJet",60,-3.,3.); ResolutionRapidityRecoGenJet->Sumw2();
     ResolutionDeltaPhiRecoGenJet = fs->make<TH1F>("ResolutionDeltaPhiRecoGenJet","ResolutionDeltaPhiRecoGenJet",2000,-3.,3.); ResolutionDeltaPhiRecoGenJet->Sumw2();

     ResponsePhiPartonGenJet = fs->make<TH2F>("ResponsePhiPartonGenJet","ResponsePhiPartonGenJet",60, 0,2*TMath::Pi(),60, 0,2*TMath::Pi()); ResponsePhiPartonGenJet->Sumw2();
     ResponsePhiRecoGenJet = fs->make<TH2F>("ResponsePhiRecoGenJet","ResponsePhiRecoGenJet",60, 0,2*TMath::Pi(),60, 0,2*TMath::Pi()); ResponsePhiRecoGenJet->Sumw2();
     ResponseDeltaPhiRecoGenJet = fs->make<TH2F>("ResponseDeltaPhiRecoGenJet","ResponseDeltaPhiRecoGenJet",DeltaPhibins,DeltaPhibinning,DeltaPhibins,DeltaPhibinning); ResponseDeltaPhiRecoGenJet->Sumw2();

     ResponseForTUnfoldDeltaPhiRecoGenJet = fs->make<TH2F>("ResponseForTUnfoldDeltaPhiRecoGenJet","ResponseForTUnfoldDeltaPhiRecoGenJet",DeltaPhibinsLarge,DeltaPhibinningLarge,DeltaPhibins,DeltaPhibinning); ResponseForTUnfoldDeltaPhiRecoGenJet->Sumw2();

     PtGENSubJet0  = fs->make<TH1F>("PtGENSubJet0","PtGENSubJet0",Ptbins,Ptbinning); PtGENSubJet0->Sumw2();
     PtGENSubJet1  = fs->make<TH1F>("PtGENSubJet1","PtGENSubJet1",Ptbins,Ptbinning); PtGENSubJet1->Sumw2();
     EtaGENSubJet0 = fs->make<TH1F>("EtaGENSubJet0","EtaGENSubJet0",60,-3.,3.); EtaGENSubJet0->Sumw2();
     EtaGENSubJet1 = fs->make<TH1F>("EtaGENSubJet1","EtaGENSubJet1",60,-3.,3.); EtaGENSubJet1->Sumw2();
     MassGENSubJet0 = fs->make<TH1F>("MassGENSubJet0","MassGENSubJet0",200,0.,2000.); MassGENSubJet0->Sumw2();
     MassGENSubJet1 = fs->make<TH1F>("MassGENSubJet1","MassGENSubJet1",200,0.,2000.); MassGENSubJet1->Sumw2();
     PhiGENSubJet0 = fs->make<TH1F>("PhiGENSubJet0","PhiGENSubJet0",60,-TMath::Pi(),TMath::Pi()); PhiGENSubJet0->Sumw2();
     PhiGENSubJet1 = fs->make<TH1F>("PhiGENSubJet1","PhiGENSubJet1",60,-TMath::Pi(),TMath::Pi()); PhiGENSubJet1->Sumw2();

     DeltaRGENSubJets= fs->make<TH1F>("DeltaRGENSubJets","DeltaRGENSubJets",100,0.,6.28); DeltaRGENSubJets->Sumw2();

     MVADesy = fs->make<TH1F>("MVADesy","MVADesy",200,-2.0,2.0); MVADesy->Sumw2();

     PtMatchedGenJet  = fs->make<TH1F>("PtMatchedGenJet","PtMatchedGenJet",Ptbins,Ptbinning); PtMatchedGenJet->Sumw2();
     PtMatchedGenJetAndBJet  = fs->make<TH1F>("PtMatchedGenJetAndBJet","PtMatchedGenJetAndBJet",Ptbins,Ptbinning); PtMatchedGenJetAndBJet->Sumw2();
     PtMatchedGenJetAndWJet  = fs->make<TH1F>("PtMatchedGenJetAndWJet","PtMatchedGenJetAndWJet",Ptbins,Ptbinning); PtMatchedGenJetAndWJet->Sumw2();
     PtMatchedGenJetAndBWJet  = fs->make<TH1F>("PtMatchedGenJetAndBWJet","PtMatchedGenJetAndBWJet",Ptbins,Ptbinning); PtMatchedGenJetAndBWJet->Sumw2();
     PtMatchedGenJetAndWBosonJet  = fs->make<TH1F>("PtMatchedGenJetAndWBosonJet","PtMatchedGenJetAndWBosonJet",Ptbins,Ptbinning); PtMatchedGenJetAndWBosonJet->Sumw2();
     PtMatchedGenJetAndBWBosonJet  = fs->make<TH1F>("PtMatchedGenJetAndBWBosonJet","PtMatchedGenJetAndBWBosonJet",Ptbins,Ptbinning); PtMatchedGenJetAndBWBosonJet->Sumw2();
     PtMatchedGenJetAndWSubjetMassJet  = fs->make<TH1F>("PtMatchedGenJetAndWSubjetMassJet","PtMatchedGenJetAndWSubjetMassJet",Ptbins,Ptbinning); PtMatchedGenJetAndWSubjetMassJet->Sumw2();
     PtMatchedGenJetAndBWSubjetMassJet  = fs->make<TH1F>("PtMatchedGenJetAndBWSubjetMassJet","PtMatchedGenJetAndBWSubjetMassJet",Ptbins,Ptbinning); PtMatchedGenJetAndBWSubjetMassJet->Sumw2();
     PtMatchedGenJetAndWSubjetMassHighJet  = fs->make<TH1F>("PtMatchedGenJetAndWSubjetMassHighJet","PtMatchedGenJetAndWSubjetMassHighJet",Ptbins,Ptbinning); PtMatchedGenJetAndWSubjetMassHighJet->Sumw2();
     PtMatchedGenJetAndBWSubjetMassHighJet  = fs->make<TH1F>("PtMatchedGenJetAndBWSubjetMassHighJet","PtMatchedGenJetAndBWSubjetMassHighJet",Ptbins,Ptbinning); PtMatchedGenJetAndBWSubjetMassHighJet->Sumw2();

     PtMatchedRecoJet  = fs->make<TH1F>("PtMatchedRecoJet","PtMatchedRecoJet",Ptbins,Ptbinning); PtMatchedRecoJet->Sumw2();
     PtMatchedRecoJetSubBtags  = fs->make<TH1F>("PtMatchedRecoJetSubBtags","PtMatchedRecoJetSubBtags",Ptbins,Ptbinning); PtMatchedRecoJetSubBtags->Sumw2();

     PtMatchedGenJetAndBWJetLeading = fs->make<TH1F>("PtMatchedGenJetAndBWJetLeading","PtMatchedGenJetAndBWJetLeading",Ptbins,Ptbinning); PtMatchedGenJetAndBWJetLeading->Sumw2();
     PtMatchedGenJetAndBWJetLeadingReco = fs->make<TH1F>("PtMatchedGenJetAndBWJetLeadingReco","PtMatchedGenJetAndBWJetLeadingReco",Ptbins,Ptbinning); PtMatchedGenJetAndBWJetLeadingReco->Sumw2();
     PtMatchedGenJetAndBWJetSubLeading = fs->make<TH1F>("PtMatchedGenJetAndBWJetSubLeading","PtMatchedGenJetAndBWJetSubLeading",Ptbins,Ptbinning); PtMatchedGenJetAndBWJetSubLeading->Sumw2();
     PtMatchedGenJetAndBWJetSubLeadingReco = fs->make<TH1F>("PtMatchedGenJetAndBWJetSubLeadingReco","PtMatchedGenJetAndBWJetSubLeadingReco",Ptbins,Ptbinning); PtMatchedGenJetAndBWJetSubLeadingReco->Sumw2();

     discr_ = new BoostedDiscriminatorMVADesy();

     //After prepreselection
     PtDETLeadingPrepreselectionSubJet0  = fs->make<TH1F>("PtDETLeadingPrepreselectionSubJet0","PtDETLeadingPrepreselectionSubJet0",Ptbins,Ptbinning); PtDETLeadingPrepreselectionSubJet0->Sumw2();
     PtDETLeadingPrepreselectionSubJet1  = fs->make<TH1F>("PtDETLeadingPrepreselectionSubJet1","PtDETLeadingPrepreselectionSubJet1",Ptbins,Ptbinning); PtDETLeadingPrepreselectionSubJet1->Sumw2();
     EtaDETLeadingPrepreselectionSubJet0 = fs->make<TH1F>("EtaDETLeadingPrepreselectionSubJet0","EtaDETLeadingPrepreselectionSubJet0",60,-3.,3.); EtaDETLeadingPrepreselectionSubJet0->Sumw2();
     EtaDETLeadingPrepreselectionSubJet1 = fs->make<TH1F>("EtaDETLeadingPrepreselectionSubJet1","EtaDETLeadingPrepreselectionSubJet1",60,-3.,3.); EtaDETLeadingPrepreselectionSubJet1->Sumw2();
     MassDETLeadingPrepreselectionSubJet0 = fs->make<TH1F>("MassDETLeadingPrepreselectionSubJet0","MassDETLeadingPrepreselectionSubJet0",SubJetMassbins,SubJetMassbinning); MassDETLeadingPrepreselectionSubJet0->Sumw2();
     MassDETLeadingPrepreselectionSubJet1 = fs->make<TH1F>("MassDETLeadingPrepreselectionSubJet1","MassDETLeadingPrepreselectionSubJet1",SubJetMassbins,SubJetMassbinning); MassDETLeadingPrepreselectionSubJet1->Sumw2();
     PhiDETLeadingPrepreselectionSubJet0 = fs->make<TH1F>("PhiDETLeadingPrepreselectionSubJet0","PhiDETLeadingPrepreselectionSubJet0",60,-TMath::Pi(),TMath::Pi()); PhiDETLeadingPrepreselectionSubJet0->Sumw2();
     PhiDETLeadingPrepreselectionSubJet1 = fs->make<TH1F>("PhiDETLeadingPrepreselectionSubJet1","PhiDETLeadingPrepreselectionSubJet1",60,-TMath::Pi(),TMath::Pi()); PhiDETLeadingPrepreselectionSubJet1->Sumw2();
     DeltaRDETLeadingPrepreselectionSubJets = fs->make<TH1F>("DeltaRDETLeadingPrepreselectionSubJets","DeltaRDETLeadingPrepreselectionSubJets",60,0,2*TMath::Pi()); DeltaRDETLeadingPrepreselectionSubJets->Sumw2();

     Tau1DETLeadingJetPrepreselection = fs->make<TH1F>("Tau1DETLeadingJetPrepreselection","Tau1DETLeadingJetPrepreselection",100,0.,1.); Tau1DETLeadingJetPrepreselection->Sumw2();
     Tau2DETLeadingJetPrepreselection = fs->make<TH1F>("Tau2DETLeadingJetPrepreselection","Tau2DETLeadingJetPrepreselection",100,0.,1.);Tau2DETLeadingJetPrepreselection->Sumw2();
     Tau3DETLeadingJetPrepreselection = fs->make<TH1F>("Tau3DETLeadingJetPrepreselection","Tau3DETLeadingJetPrepreselection",100,0.,1.);Tau3DETLeadingJetPrepreselection->Sumw2();
     Tau31DETLeadingJetPrepreselection = fs->make<TH1F>("Tau31DETLeadingJetPrepreselection","Tau31DETLeadingJetPrepreselection",100,0.,1.);Tau31DETLeadingJetPrepreselection->Sumw2();
     Tau32DETLeadingJetPrepreselection = fs->make<TH1F>("Tau32DETLeadingJetPrepreselection","Tau32DETLeadingJetPrepreselection",100,0.,1.);Tau32DETLeadingJetPrepreselection->Sumw2();

     PtDETSubLeadingPrepreselectionSubJet0  = fs->make<TH1F>("PtDETSubLeadingPrepreselectionSubJet0","PtDETSubLeadingPrepreselectionSubJet0",Ptbins,Ptbinning); PtDETSubLeadingPrepreselectionSubJet0->Sumw2();
     PtDETSubLeadingPrepreselectionSubJet1  = fs->make<TH1F>("PtDETSubLeadingPrepreselectionSubJet1","PtDETSubLeadingPrepreselectionSubJet1",Ptbins,Ptbinning); PtDETSubLeadingPrepreselectionSubJet1->Sumw2();
     EtaDETSubLeadingPrepreselectionSubJet0 = fs->make<TH1F>("EtaDETSubLeadingPrepreselectionSubJet0","EtaDETSubLeadingPrepreselectionSubJet0",60,-3.,3.); EtaDETSubLeadingPrepreselectionSubJet0->Sumw2();
     EtaDETSubLeadingPrepreselectionSubJet1 = fs->make<TH1F>("EtaDETSubLeadingPrepreselectionSubJet1","EtaDETSubLeadingPrepreselectionSubJet1",60,-3.,3.); EtaDETSubLeadingPrepreselectionSubJet1->Sumw2();
     MassDETSubLeadingPrepreselectionSubJet0 = fs->make<TH1F>("MassDETSubLeadingPrepreselectionSubJet0","MassDETSubLeadingPrepreselectionSubJet0",SubJetMassbins,SubJetMassbinning); MassDETSubLeadingPrepreselectionSubJet0->Sumw2();
     MassDETSubLeadingPrepreselectionSubJet1 = fs->make<TH1F>("MassDETSubLeadingPrepreselectionSubJet1","MassDETSubLeadingPrepreselectionSubJet1",SubJetMassbins,SubJetMassbinning); MassDETSubLeadingPrepreselectionSubJet1->Sumw2();
     PhiDETSubLeadingPrepreselectionSubJet0 = fs->make<TH1F>("PhiDETSubLeadingPrepreselectionSubJet0","PhiDETSubLeadingPrepreselectionSubJet0",60,-TMath::Pi(),TMath::Pi()); PhiDETSubLeadingPrepreselectionSubJet0->Sumw2();
     PhiDETSubLeadingPrepreselectionSubJet1 = fs->make<TH1F>("PhiDETSubLeadingPrepreselectionSubJet1","PhiDETSubLeadingPrepreselectionSubJet1",60,-TMath::Pi(),TMath::Pi()); PhiDETSubLeadingPrepreselectionSubJet1->Sumw2();
     DeltaRDETSubLeadingPrepreselectionSubJets = fs->make<TH1F>("DeltaRDETSubLeadingPrepreselectionSubJets","DeltaRDETSubLeadingPrepreselectionSubJets",60,0,2*TMath::Pi()); DeltaRDETSubLeadingPrepreselectionSubJets->Sumw2();

     Tau1DETSubLeadingJetPrepreselection = fs->make<TH1F>("Tau1DETSubLeadingJetPrepreselection","Tau1DETSubLeadingJetPrepreselection",100,0.,1.);Tau1DETSubLeadingJetPrepreselection->Sumw2();
     Tau2DETSubLeadingJetPrepreselection = fs->make<TH1F>("Tau2DETSubLeadingJetPrepreselection","Tau2DETSubLeadingJetPrepreselection",100,0.,1.);Tau2DETSubLeadingJetPrepreselection->Sumw2();
     Tau3DETSubLeadingJetPrepreselection = fs->make<TH1F>("Tau3DETSubLeadingJetPrepreselection","Tau3DETSubLeadingJetPrepreselection",100,0.,1.);Tau3DETSubLeadingJetPrepreselection->Sumw2();
     Tau31DETSubLeadingJetPrepreselection = fs->make<TH1F>("Tau31DETSubLeadingJetPrepreselection","Tau31DETSubLeadingJetPrepreselection",100,0.,1.);Tau31DETSubLeadingJetPrepreselection->Sumw2();
     Tau32DETSubLeadingJetPrepreselection = fs->make<TH1F>("Tau32DETSubLeadingJetPrepreselection","Tau32DETSubLeadingJetPrepreselection",100,0.,1.);Tau32DETSubLeadingJetPrepreselection->Sumw2();

     PtDETLeadingPrepreselection  = fs->make<TH1F>("PtDETLeadingPrepreselection","PtDETLeadingPrepreselection",Ptbins,Ptbinning); PtDETLeadingPrepreselection->Sumw2();
     PtDETSubLeadingPrepreselection  = fs->make<TH1F>("PtDETSubLeadingPrepreselection","PtDETSubLeadingPrepreselection",Ptbins,Ptbinning); PtDETSubLeadingPrepreselection->Sumw2();
     EtaDETLeadingPrepreselection = fs->make<TH1F>("EtaDETLeadingPrepreselection","EtaDETLeadingPrepreselection",60,-3.,3.); EtaDETLeadingPrepreselection->Sumw2();
     EtaDETSubLeadingPrepreselection = fs->make<TH1F>("EtaDETSubLeadingPrepreselection","EtaDETSubLeadingPrepreselection",60,-3.,3.); EtaDETSubLeadingPrepreselection->Sumw2();
     MassDETLeadingPrepreselection = fs->make<TH1F>("MassDETLeadingPrepreselection","MassDETLeadingPrepreselection",Massbins,Massbinning); MassDETLeadingPrepreselection->Sumw2();
     MassDETSubLeadingPrepreselection = fs->make<TH1F>("MassDETSubLeadingPrepreselection","MassDETSubLeadingPrepreselection",Massbins,Massbinning); MassDETSubLeadingPrepreselection->Sumw2();
     PhiDETLeadingPrepreselection = fs->make<TH1F>("PhiDETLeadingPrepreselection","PhiDETLeadingPrepreselection",60,-TMath::Pi(),TMath::Pi()); PhiDETLeadingPrepreselection->Sumw2();
     PhiDETSubLeadingPrepreselection = fs->make<TH1F>("PhiDETSubLeadingPrepreselection","PhiDETSubLeadingPrepreselection",60,-TMath::Pi(),TMath::Pi()); PhiDETSubLeadingPrepreselection->Sumw2();
     SDMassDETLeadingPrepreselection = fs->make<TH1F>("SDMassDETLeadingPrepreselection","SDMassDETLeadingPrepreselection",Massbins,Massbinning); SDMassDETLeadingPrepreselection->Sumw2();
     SDMassDETSubLeadingPrepreselection = fs->make<TH1F>("SDMassDETSubLeadingPrepreselection","SDMassDETSubLeadingPrepreselection",Massbins,Massbinning); SDMassDETSubLeadingPrepreselection->Sumw2();

     //After preselection
     PtDETLeadingPreselectionSubJet0  = fs->make<TH1F>("PtDETLeadingPreselectionSubJet0","PtDETLeadingPreselectionSubJet0",Ptbins,Ptbinning); PtDETLeadingPreselectionSubJet0->Sumw2();
     PtDETLeadingPreselectionSubJet1  = fs->make<TH1F>("PtDETLeadingPreselectionSubJet1","PtDETLeadingPreselectionSubJet1",Ptbins,Ptbinning); PtDETLeadingPreselectionSubJet1->Sumw2();
     EtaDETLeadingPreselectionSubJet0 = fs->make<TH1F>("EtaDETLeadingPreselectionSubJet0","EtaDETLeadingPreselectionSubJet0",60,-3.,3.); EtaDETLeadingPreselectionSubJet0->Sumw2();
     EtaDETLeadingPreselectionSubJet1 = fs->make<TH1F>("EtaDETLeadingPreselectionSubJet1","EtaDETLeadingPreselectionSubJet1",60,-3.,3.); EtaDETLeadingPreselectionSubJet1->Sumw2();
     MassDETLeadingPreselectionSubJet0 = fs->make<TH1F>("MassDETLeadingPreselectionSubJet0","MassDETLeadingPreselectionSubJet0",SubJetMassbins,SubJetMassbinning); MassDETLeadingPreselectionSubJet0->Sumw2();
     MassDETLeadingPreselectionSubJet1 = fs->make<TH1F>("MassDETLeadingPreselectionSubJet1","MassDETLeadingPreselectionSubJet1",SubJetMassbins,SubJetMassbinning); MassDETLeadingPreselectionSubJet1->Sumw2();
     PhiDETLeadingPreselectionSubJet0 = fs->make<TH1F>("PhiDETLeadingPreselectionSubJet0","PhiDETLeadingPreselectionSubJet0",60,-TMath::Pi(),TMath::Pi()); PhiDETLeadingPreselectionSubJet0->Sumw2();
     PhiDETLeadingPreselectionSubJet1 = fs->make<TH1F>("PhiDETLeadingPreselectionSubJet1","PhiDETLeadingPreselectionSubJet1",60,-TMath::Pi(),TMath::Pi()); PhiDETLeadingPreselectionSubJet1->Sumw2();
     DeltaRDETLeadingPreselectionSubJets = fs->make<TH1F>("DeltaRDETLeadingPreselectionSubJets","DeltaRDETLeadingPreselectionSubJets",60,0,2*TMath::Pi()); DeltaRDETLeadingPreselectionSubJets->Sumw2();

     Tau1DETLeadingJetPreselection = fs->make<TH1F>("Tau1DETLeadingJetPreselection","Tau1DETLeadingJetPreselection",100,0.,1.); Tau1DETLeadingJetPreselection->Sumw2();
     Tau2DETLeadingJetPreselection = fs->make<TH1F>("Tau2DETLeadingJetPreselection","Tau2DETLeadingJetPreselection",100,0.,1.);Tau2DETLeadingJetPreselection->Sumw2();
     Tau3DETLeadingJetPreselection = fs->make<TH1F>("Tau3DETLeadingJetPreselection","Tau3DETLeadingJetPreselection",100,0.,1.);Tau3DETLeadingJetPreselection->Sumw2();
     Tau31DETLeadingJetPreselection = fs->make<TH1F>("Tau31DETLeadingJetPreselection","Tau31DETLeadingJetPreselection",100,0.,1.);Tau31DETLeadingJetPreselection->Sumw2();
     Tau32DETLeadingJetPreselection = fs->make<TH1F>("Tau32DETLeadingJetPreselection","Tau32DETLeadingJetPreselection",100,0.,1.);Tau32DETLeadingJetPreselection->Sumw2();

     PtDETSubLeadingPreselectionSubJet0  = fs->make<TH1F>("PtDETSubLeadingPreselectionSubJet0","PtDETSubLeadingPreselectionSubJet0",Ptbins,Ptbinning); PtDETSubLeadingPreselectionSubJet0->Sumw2();
     PtDETSubLeadingPreselectionSubJet1  = fs->make<TH1F>("PtDETSubLeadingPreselectionSubJet1","PtDETSubLeadingPreselectionSubJet1",Ptbins,Ptbinning); PtDETSubLeadingPreselectionSubJet1->Sumw2();
     EtaDETSubLeadingPreselectionSubJet0 = fs->make<TH1F>("EtaDETSubLeadingPreselectionSubJet0","EtaDETSubLeadingPreselectionSubJet0",60,-3.,3.); EtaDETSubLeadingPreselectionSubJet0->Sumw2();
     EtaDETSubLeadingPreselectionSubJet1 = fs->make<TH1F>("EtaDETSubLeadingPreselectionSubJet1","EtaDETSubLeadingPreselectionSubJet1",60,-3.,3.); EtaDETSubLeadingPreselectionSubJet1->Sumw2();
     MassDETSubLeadingPreselectionSubJet0 = fs->make<TH1F>("MassDETSubLeadingPreselectionSubJet0","MassDETSubLeadingPreselectionSubJet0",SubJetMassbins,SubJetMassbinning); MassDETSubLeadingPreselectionSubJet0->Sumw2();
     MassDETSubLeadingPreselectionSubJet1 = fs->make<TH1F>("MassDETSubLeadingPreselectionSubJet1","MassDETSubLeadingPreselectionSubJet1",SubJetMassbins,SubJetMassbinning); MassDETSubLeadingPreselectionSubJet1->Sumw2();
     PhiDETSubLeadingPreselectionSubJet0 = fs->make<TH1F>("PhiDETSubLeadingPreselectionSubJet0","PhiDETSubLeadingPreselectionSubJet0",60,-TMath::Pi(),TMath::Pi()); PhiDETSubLeadingPreselectionSubJet0->Sumw2();
     PhiDETSubLeadingPreselectionSubJet1 = fs->make<TH1F>("PhiDETSubLeadingPreselectionSubJet1","PhiDETSubLeadingPreselectionSubJet1",60,-TMath::Pi(),TMath::Pi()); PhiDETSubLeadingPreselectionSubJet1->Sumw2();
     DeltaRDETSubLeadingPreselectionSubJets = fs->make<TH1F>("DeltaRDETSubLeadingPreselectionSubJets","DeltaRDETSubLeadingPreselectionSubJets",60,0,2*TMath::Pi()); DeltaRDETSubLeadingPreselectionSubJets->Sumw2();

     Tau1DETSubLeadingJetPreselection = fs->make<TH1F>("Tau1DETSubLeadingJetPreselection","Tau1DETSubLeadingJetPreselection",100,0.,1.);Tau1DETSubLeadingJetPreselection->Sumw2();
     Tau2DETSubLeadingJetPreselection = fs->make<TH1F>("Tau2DETSubLeadingJetPreselection","Tau2DETSubLeadingJetPreselection",100,0.,1.);Tau2DETSubLeadingJetPreselection->Sumw2();
     Tau3DETSubLeadingJetPreselection = fs->make<TH1F>("Tau3DETSubLeadingJetPreselection","Tau3DETSubLeadingJetPreselection",100,0.,1.);Tau3DETSubLeadingJetPreselection->Sumw2();
     Tau31DETSubLeadingJetPreselection = fs->make<TH1F>("Tau31DETSubLeadingJetPreselection","Tau31DETSubLeadingJetPreselection",100,0.,1.);Tau31DETSubLeadingJetPreselection->Sumw2();
     Tau32DETSubLeadingJetPreselection = fs->make<TH1F>("Tau32DETSubLeadingJetPreselection","Tau32DETSubLeadingJetPreselection",100,0.,1.);Tau32DETSubLeadingJetPreselection->Sumw2();

     PtDETLeadingPreselection  = fs->make<TH1F>("PtDETLeadingPreselection","PtDETLeadingPreselection",Ptbins,Ptbinning); PtDETLeadingPreselection->Sumw2();
     PtDETSubLeadingPreselection  = fs->make<TH1F>("PtDETSubLeadingPreselection","PtDETSubLeadingPreselection",Ptbins,Ptbinning); PtDETSubLeadingPreselection->Sumw2();
     EtaDETLeadingPreselection = fs->make<TH1F>("EtaDETLeadingPreselection","EtaDETLeadingPreselection",60,-3.,3.); EtaDETLeadingPreselection->Sumw2();
     EtaDETSubLeadingPreselection = fs->make<TH1F>("EtaDETSubLeadingPreselection","EtaDETSubLeadingPreselection",60,-3.,3.); EtaDETSubLeadingPreselection->Sumw2();
     MassDETLeadingPreselection = fs->make<TH1F>("MassDETLeadingPreselection","MassDETLeadingPreselection",Massbins,Massbinning); MassDETLeadingPreselection->Sumw2();
     MassDETSubLeadingPreselection = fs->make<TH1F>("MassDETSubLeadingPreselection","MassDETSubLeadingPreselection",Massbins,Massbinning); MassDETSubLeadingPreselection->Sumw2();
     PhiDETLeadingPreselection = fs->make<TH1F>("PhiDETLeadingPreselection","PhiDETLeadingPreselection",60,-TMath::Pi(),TMath::Pi()); PhiDETLeadingPreselection->Sumw2();
     PhiDETSubLeadingPreselection = fs->make<TH1F>("PhiDETSubLeadingPreselection","PhiDETSubLeadingPreselection",60,-TMath::Pi(),TMath::Pi()); PhiDETSubLeadingPreselection->Sumw2();
     SDMassDETLeadingPreselection = fs->make<TH1F>("SDMassDETLeadingPreselection","SDMassDETLeadingPreselection",Massbins,Massbinning); SDMassDETLeadingPreselection->Sumw2();
     SDMassDETSubLeadingPreselection = fs->make<TH1F>("SDMassDETSubLeadingPreselection","SDMassDETSubLeadingPreselection",Massbins,Massbinning); SDMassDETSubLeadingPreselection->Sumw2();

     //After selection
     PtDETLeadingSelectionSubJet0  = fs->make<TH1F>("PtDETLeadingSelectionSubJet0","PtDETLeadingSelectionSubJet0",Ptbins,Ptbinning); PtDETLeadingSelectionSubJet0->Sumw2();
     PtDETLeadingSelectionSubJet1  = fs->make<TH1F>("PtDETLeadingSelectionSubJet1","PtDETLeadingSelectionSubJet1",Ptbins,Ptbinning); PtDETLeadingSelectionSubJet1->Sumw2();
     EtaDETLeadingSelectionSubJet0 = fs->make<TH1F>("EtaDETLeadingSelectionSubJet0","EtaDETLeadingSelectionSubJet0",60,-3.,3.); EtaDETLeadingSelectionSubJet0->Sumw2();
     EtaDETLeadingSelectionSubJet1 = fs->make<TH1F>("EtaDETLeadingSelectionSubJet1","EtaDETLeadingSelectionSubJet1",60,-3.,3.); EtaDETLeadingSelectionSubJet1->Sumw2();
     MassDETLeadingSelectionSubJet0 = fs->make<TH1F>("MassDETLeadingSelectionSubJet0","MassDETLeadingSelectionSubJet0",SubJetMassbins,SubJetMassbinning); MassDETLeadingSelectionSubJet0->Sumw2();
     MassDETLeadingSelectionSubJet1 = fs->make<TH1F>("MassDETLeadingSelectionSubJet1","MassDETLeadingSelectionSubJet1",SubJetMassbins,SubJetMassbinning); MassDETLeadingSelectionSubJet1->Sumw2();
     PhiDETLeadingSelectionSubJet0 = fs->make<TH1F>("PhiDETLeadingSelectionSubJet0","PhiDETLeadingSelectionSubJet0",60,-TMath::Pi(),TMath::Pi()); PhiDETLeadingSelectionSubJet0->Sumw2();
     PhiDETLeadingSelectionSubJet1 = fs->make<TH1F>("PhiDETLeadingSelectionSubJet1","PhiDETLeadingSelectionSubJet1",60,-TMath::Pi(),TMath::Pi()); PhiDETLeadingSelectionSubJet1->Sumw2();
     DeltaRDETLeadingSelectionSubJets = fs->make<TH1F>("DeltaRDETLeadingSelectionSubJets","DeltaRDETLeadingSelectionSubJets",60,0,2*TMath::Pi()); DeltaRDETLeadingSelectionSubJets->Sumw2();

     Tau1DETLeadingJetSelection = fs->make<TH1F>("Tau1DETLeadingJetSelection","Tau1DETLeadingJetSelection",100,0.,1.);Tau1DETLeadingJetSelection->Sumw2();
     Tau2DETLeadingJetSelection = fs->make<TH1F>("Tau2DETLeadingJetSelection","Tau2DETLeadingJetSelection",100,0.,1.);Tau2DETLeadingJetSelection->Sumw2();
     Tau3DETLeadingJetSelection = fs->make<TH1F>("Tau3DETLeadingJetSelection","Tau3DETLeadingJetSelection",100,0.,1.);Tau3DETLeadingJetSelection->Sumw2();
     Tau31DETLeadingJetSelection = fs->make<TH1F>("Tau31DETLeadingJetSelection","Tau31DETLeadingJetSelection",100,0.,1.);Tau31DETLeadingJetSelection->Sumw2();
     Tau32DETLeadingJetSelection = fs->make<TH1F>("Tau32DETLeadingJetSelection","Tau32DETLeadingJetSelection",100,0.,1.);Tau32DETLeadingJetSelection->Sumw2();

     PtDETSubLeadingSelectionSubJet0  = fs->make<TH1F>("PtDETSubLeadingSelectionSubJet0","PtDETSubLeadingSelectionSubJet0",Ptbins,Ptbinning); PtDETSubLeadingSelectionSubJet0->Sumw2();
     PtDETSubLeadingSelectionSubJet1  = fs->make<TH1F>("PtDETSubLeadingSelectionSubJet1","PtDETSubLeadingSelectionSubJet1",Ptbins,Ptbinning); PtDETSubLeadingSelectionSubJet1->Sumw2();
     EtaDETSubLeadingSelectionSubJet0 = fs->make<TH1F>("EtaDETSubLeadingSelectionSubJet0","EtaDETSubLeadingSelectionSubJet0",60,-3.,3.); EtaDETSubLeadingSelectionSubJet0->Sumw2();
     EtaDETSubLeadingSelectionSubJet1 = fs->make<TH1F>("EtaDETSubLeadingSelectionSubJet1","EtaDETSubLeadingSelectionSubJet1",60,-3.,3.); EtaDETSubLeadingSelectionSubJet1->Sumw2();
     MassDETSubLeadingSelectionSubJet0 = fs->make<TH1F>("MassDETSubLeadingSelectionSubJet0","MassDETSubLeadingSelectionSubJet0",SubJetMassbins,SubJetMassbinning); MassDETSubLeadingSelectionSubJet0->Sumw2();
     MassDETSubLeadingSelectionSubJet1 = fs->make<TH1F>("MassDETSubLeadingSelectionSubJet1","MassDETSubLeadingSelectionSubJet1",SubJetMassbins,SubJetMassbinning); MassDETSubLeadingSelectionSubJet1->Sumw2();
     PhiDETSubLeadingSelectionSubJet0 = fs->make<TH1F>("PhiDETSubLeadingSelectionSubJet0","PhiDETSubLeadingSelectionSubJet0",60,-TMath::Pi(),TMath::Pi()); PhiDETSubLeadingSelectionSubJet0->Sumw2();
     PhiDETSubLeadingSelectionSubJet1 = fs->make<TH1F>("PhiDETSubLeadingSelectionSubJet1","PhiDETSubLeadingSelectionSubJet1",60,-TMath::Pi(),TMath::Pi()); PhiDETSubLeadingSelectionSubJet1->Sumw2();
     DeltaRDETSubLeadingSelectionSubJets = fs->make<TH1F>("DeltaRDETSubLeadingSelectionSubJets","DeltaRDETSubLeadingSelectionSubJets",60,0,2*TMath::Pi()); DeltaRDETSubLeadingSelectionSubJets->Sumw2();

     Tau1DETSubLeadingJetSelection = fs->make<TH1F>("Tau1DETSubLeadingJetSelection","Tau1DETSubLeadingJetSelection",100,0.,1.);Tau1DETSubLeadingJetSelection->Sumw2();
     Tau2DETSubLeadingJetSelection = fs->make<TH1F>("Tau2DETSubLeadingJetSelection","Tau2DETSubLeadingJetSelection",100,0.,1.);Tau2DETSubLeadingJetSelection->Sumw2();
     Tau3DETSubLeadingJetSelection = fs->make<TH1F>("Tau3DETSubLeadingJetSelection","Tau3DETSubLeadingJetSelection",100,0.,1.);Tau3DETSubLeadingJetSelection->Sumw2();
     Tau31DETSubLeadingJetSelection = fs->make<TH1F>("Tau31DETSubLeadingJetSelection","Tau31DETSubLeadingJetSelection",100,0.,1.);Tau31DETSubLeadingJetSelection->Sumw2();
     Tau32DETSubLeadingJetSelection = fs->make<TH1F>("Tau32DETSubLeadingJetSelection","Tau32DETSubLeadingJetSelection",100,0.,1.);Tau32DETSubLeadingJetSelection->Sumw2();

     FlavourDETSelectionLeadingJet = fs->make<TH1F>("FlavourDETSelectionLeadingJet","FlavourDETSelectionLeadingJet",23,-0.5,22.5); FlavourDETSelectionLeadingJet->Sumw2();
     FlavourDETSelectionSubLeadingJet = fs->make<TH1F>("FlavourDETSelectionSubLeadingJet","FlavourDETSelectionSubLeadingJet",23,-0.5,22.5); FlavourDETSelectionSubLeadingJet->Sumw2();
     FlavourDETPrepreselectionLeadingJet = fs->make<TH1F>("FlavourDETPrepreselectionLeadingJet","FlavourDETPrepreselectionLeadingJet",23,-0.5,22.5); FlavourDETPrepreselectionLeadingJet->Sumw2();
     FlavourDETPrepreselectionSubLeadingJet = fs->make<TH1F>("FlavourDETPrepreselectionSubLeadingJet","FlavourDETPrepreselectionSubLeadingJet",23,-0.5,22.5); FlavourDETPrepreselectionSubLeadingJet->Sumw2();

     PtDETLeadingPrepreselectionBTrue  = fs->make<TH1F>("PtDETLeadingPrepreselectionBTrue","PtDETLeadingPrepreselectionBTrue",Ptbins,Ptbinning); PtDETLeadingPrepreselectionBTrue->Sumw2();
     PtDETSubLeadingPrepreselectionBTrue  = fs->make<TH1F>("PtDETSubLeadingPrepreselectionBTrue","PtDETSubLeadingPrepreselectionBTrue",Ptbins,Ptbinning); PtDETSubLeadingPrepreselectionBTrue->Sumw2();

     PtLeadingGenJet = fs->make<TH1F>("PtLeadingGenJet","PtLeadingGenJet",Ptbins,Ptbinning); PtLeadingGenJet->Sumw2();
     PtSubLeadingGenJet = fs->make<TH1F>("PtSubLeadingGenJet","PtSubLeadingGenJet",Ptbins,Ptbinning); PtSubLeadingGenJet->Sumw2();

     PtLeadingGenJetMatched = fs->make<TH1F>("PtLeadingGenJetMatched","PtLeadingGenJetMatched",Ptbins,Ptbinning); PtLeadingGenJetMatched->Sumw2();
     PtSubLeadingGenJetMatched = fs->make<TH1F>("PtSubLeadingGenJetMatched","PtSubLeadingGenJetMatched",Ptbins,Ptbinning); PtSubLeadingGenJetMatched->Sumw2();

     PtDETLeadingSelection  = fs->make<TH1F>("PtDETLeadingSelection","PtDETLeadingSelection",Ptbins,Ptbinning); PtDETLeadingSelection->Sumw2();
     PtDETSubLeadingSelection  = fs->make<TH1F>("PtDETSubLeadingSelection","PtDETSubLeadingSelection",Ptbins,Ptbinning); PtDETSubLeadingSelection->Sumw2();
     EtaDETLeadingSelection = fs->make<TH1F>("EtaDETLeadingSelection","EtaDETLeadingSelection",60,-3.,3.); EtaDETLeadingSelection->Sumw2();
     EtaDETSubLeadingSelection = fs->make<TH1F>("EtaDETSubLeadingSelection","EtaDETSubLeadingSelection",60,-3.,3.); EtaDETSubLeadingSelection->Sumw2();
     MassDETLeadingSelection = fs->make<TH1F>("MassDETLeadingSelection","MassDETLeadingSelection",Massbins,Massbinning); MassDETLeadingSelection->Sumw2();
     MassDETSubLeadingSelection = fs->make<TH1F>("MassDETSubLeadingSelection","MassDETSubLeadingSelection",Massbins,Massbinning); MassDETSubLeadingSelection->Sumw2();
     PhiDETLeadingSelection = fs->make<TH1F>("PhiDETLeadingSelection","PhiDETLeadingSelection",60,-TMath::Pi(),TMath::Pi()); PhiDETLeadingSelection->Sumw2();
     PhiDETSubLeadingSelection = fs->make<TH1F>("PhiDETSubLeadingSelection","PhiDETSubLeadingSelection",60,-TMath::Pi(),TMath::Pi()); PhiDETSubLeadingSelection->Sumw2();
     SDMassDETLeadingSelection = fs->make<TH1F>("SDMassDETLeadingSelection","SDMassDETLeadingSelection",Massbins,Massbinning); SDMassDETLeadingSelection->Sumw2();
     SDMassDETSubLeadingSelection = fs->make<TH1F>("SDMassDETSubLeadingSelection","SDMassDETSubLeadingSelection",Massbins,Massbinning); SDMassDETSubLeadingSelection->Sumw2();

     //Here I define the SD mass histograms that I will fit in the control region for controlling the background (and then subtracting it from the measured observables)
     TString histo_label_SDdeltaphi;
     TString histo_title_SDdeltaphi; 

     TString histo_label_SDdeltaphiSub;
     TString histo_title_SDdeltaphiSub; 

     TString histo_label_SDPt;
     TString histo_title_SDPt; 

     TString histo_label_SDPtSub;
     TString histo_title_SDPtSub; 
     
     for(int i = 0; i < 6; i++) {

       histo_label_SDdeltaphi = "SDMassDETLeadingSelection_DeltaPhibin_" + TString::Format("%d",i+1);
       histo_title_SDdeltaphi = "SDMassDETLeadingSelection_DeltaPhibin_" + TString::Format("%d",i+1);
       
       SDMassDETLeadingSelection_DeltaPhibins[i] = fs->make<TH1F>(histo_label_SDdeltaphi,histo_title_SDdeltaphi,Massbins,Massbinning); 
       SDMassDETLeadingSelection_DeltaPhibins[i]->Sumw2();

       histo_label_SDdeltaphiSub = "SDMassDETSubLeadingSelection_DeltaPhibin_" + TString::Format("%d",i+1);
       histo_title_SDdeltaphiSub = "SDMassDETSubLeadingSelection_DeltaPhibin_" + TString::Format("%d",i+1);
       
       SDMassDETSubLeadingSelection_DeltaPhibins[i] = fs->make<TH1F>(histo_label_SDdeltaphiSub,histo_title_SDdeltaphiSub,Massbins,Massbinning); 
       SDMassDETSubLeadingSelection_DeltaPhibins[i]->Sumw2();

       histo_label_SDPt = "SDMassDETLeadingSelection_Ptbin_" + TString::Format("%d",i+1);
       histo_title_SDPt = "SDMassDETLeadingSelection_Ptbin_" + TString::Format("%d",i+1);
       
       SDMassDETLeadingSelection_Ptbins[i] = fs->make<TH1F>(histo_label_SDPt,histo_title_SDPt,Massbins,Massbinning); 
       SDMassDETLeadingSelection_Ptbins[i]->Sumw2();

       histo_label_SDPtSub = "SDMassDETSubLeadingSelection_Ptbin_" + TString::Format("%d",i+1);
       histo_title_SDPtSub = "SDMassDETSubLeadingSelection_Ptbin_" + TString::Format("%d",i+1);
       
       SDMassDETSubLeadingSelection_Ptbins[i] = fs->make<TH1F>(histo_label_SDPtSub,histo_title_SDPtSub,Massbins,Massbinning); 
       SDMassDETSubLeadingSelection_Ptbins[i]->Sumw2();

     }
          

     //RooUnfoldResponse resp_DeltaPhi(14,0,3.14159265,14,0,3.14159265, "resp_DeltaPhi", "resp_DeltaPhi");

     //To be measured
     TString histo_label_pt1;
     TString histo_title_pt1;

     TString histo_label_pt2;
     TString histo_title_pt2;

     TString histo_label_genpt1;
     TString histo_title_genpt1;

     TString histo_label_genpt2;
     TString histo_title_genpt2;

     TString histo_label_partonpt1;
     TString histo_title_partonpt1;

     TString histo_label_partonpt2;
     TString histo_title_partonpt2;

     TString ResponseLeadingJetPt;
     TString ResponseSubLeadingJetPt;

     for(int i = 0; i < 5; i++) {
       histo_label_pt1 = "PtDETLeadingSelection_" + TString::Format("%d",i+1) + "bin";
       histo_title_pt1 = "PtDETLeadingSelection_" + TString::Format("%d",i+1) + "bin";
       
       PtDETLeadingSelection_bins[i] = fs->make<TH1F>(histo_label_pt1,histo_title_pt1,PtbinsLarge,PtbinningLarge); 
       PtDETLeadingSelection_bins[i]->Sumw2();

       histo_label_pt2 = "PtDETSubLeadingSelection_" + TString::Format("%d",i+1) + "bin";
       histo_title_pt2 = "PtDETSubLeadingSelection_" + TString::Format("%d",i+1) + "bin";
       
       PtDETSubLeadingSelection_bins[i] = fs->make<TH1F>(histo_label_pt2,histo_title_pt2,PtbinsLarge,PtbinningLarge); 
       PtDETSubLeadingSelection_bins[i]->Sumw2();

       //ForToyFits
       histo_label_pt1 = "PtDETLeadingSelection_" + TString::Format("%d",i+1) + "bin_ForToyFits";
       histo_title_pt1 = "PtDETLeadingSelection_" + TString::Format("%d",i+1) + "bin_ForToyFits";
       
       PtDETLeadingSelectionForToyFits_bins[i] = fs->make<TH1F>(histo_label_pt1,histo_title_pt1,Ptbins,Ptbinning); 
       PtDETLeadingSelectionForToyFits_bins[i]->Sumw2();

       histo_label_pt2 = "PtDETSubLeadingSelection_" + TString::Format("%d",i+1) + "bin_ForToyFits";
       histo_title_pt2 = "PtDETSubLeadingSelection_" + TString::Format("%d",i+1) + "bin_ForToyFits";
       
       PtDETSubLeadingSelectionForToyFits_bins[i] = fs->make<TH1F>(histo_label_pt2,histo_title_pt2,Ptbins,Ptbinning); 
       PtDETSubLeadingSelectionForToyFits_bins[i]->Sumw2();

       //ForDoubleBins
       histo_label_pt1 = "PtDETLeadingSelection_" + TString::Format("%d",i+1) + "bin_ForDoubleBins";
       histo_title_pt1 = "PtDETLeadingSelection_" + TString::Format("%d",i+1) + "bin_ForDoubleBins";
       
       PtDETLeadingSelectionForDoubleBins_bins[i] = fs->make<TH1F>(histo_label_pt1,histo_title_pt1,PtbinsV2,PtbinningV2); 
       PtDETLeadingSelectionForDoubleBins_bins[i]->Sumw2();

       histo_label_pt2 = "PtDETSubLeadingSelection_" + TString::Format("%d",i+1) + "bin_ForDoubleBins";
       histo_title_pt2 = "PtDETSubLeadingSelection_" + TString::Format("%d",i+1) + "bin_ForDoubleBins";
       
       PtDETSubLeadingSelectionForDoubleBins_bins[i] = fs->make<TH1F>(histo_label_pt2,histo_title_pt2,PtbinsV2,PtbinningV2); 
       PtDETSubLeadingSelectionForDoubleBins_bins[i]->Sumw2();


       // GEN level

       histo_label_genpt1 = "PtGENLeadingSelection_" + TString::Format("%d",i+1) + "bin";
       histo_title_genpt1 = "PtGENLeadingSelection_" + TString::Format("%d",i+1) + "bin";
       
       PtGENLeadingSelection_bins[i] = fs->make<TH1F>(histo_label_genpt1,histo_title_genpt1,PtbinsLarge,PtbinningLarge); 
       PtGENLeadingSelection_bins[i]->Sumw2();

       histo_label_genpt2 = "PtGENSubLeadingSelection_" + TString::Format("%d",i+1) + "bin";
       histo_title_genpt2 = "PtGENSubLeadingSelection_" + TString::Format("%d",i+1) + "bin";
       
       PtGENSubLeadingSelection_bins[i] = fs->make<TH1F>(histo_label_genpt2,histo_title_genpt2,PtbinsLarge,PtbinningLarge); 
       PtGENSubLeadingSelection_bins[i]->Sumw2();

       //ForToyFits
       histo_label_pt1 = "PtGENLeadingSelection_" + TString::Format("%d",i+1) + "bin_ForToyFits";
       histo_title_pt1 = "PtGENLeadingSelection_" + TString::Format("%d",i+1) + "bin_ForToyFits";
       
       PtGENLeadingSelectionForToyFits_bins[i] = fs->make<TH1F>(histo_label_pt1,histo_title_pt1,Ptbins,Ptbinning); 
       PtGENLeadingSelectionForToyFits_bins[i]->Sumw2();

       histo_label_pt2 = "PtGENSubLeadingSelection_" + TString::Format("%d",i+1) + "bin_ForToyFits";
       histo_title_pt2 = "PtGENSubLeadingSelection_" + TString::Format("%d",i+1) + "bin_ForToyFits";
       
       PtGENSubLeadingSelectionForToyFits_bins[i] = fs->make<TH1F>(histo_label_pt2,histo_title_pt2,Ptbins,Ptbinning); 
       PtGENSubLeadingSelectionForToyFits_bins[i]->Sumw2();

       //ForDoubleBins
       histo_label_pt1 = "PtGENLeadingSelection_" + TString::Format("%d",i+1) + "bin_ForDoubleBins";
       histo_title_pt1 = "PtGENLeadingSelection_" + TString::Format("%d",i+1) + "bin_ForDoubleBins";
       
       PtGENLeadingSelectionForDoubleBins_bins[i] = fs->make<TH1F>(histo_label_pt1,histo_title_pt1,PtbinsV2,PtbinningV2); 
       PtGENLeadingSelectionForDoubleBins_bins[i]->Sumw2();

       histo_label_pt2 = "PtGENSubLeadingSelection_" + TString::Format("%d",i+1) + "bin_ForDoubleBins";
       histo_title_pt2 = "PtGENSubLeadingSelection_" + TString::Format("%d",i+1) + "bin_ForDoubleBins";
       
       PtGENSubLeadingSelectionForDoubleBins_bins[i] = fs->make<TH1F>(histo_label_pt2,histo_title_pt2,PtbinsV2,PtbinningV2); 
       PtGENSubLeadingSelectionForDoubleBins_bins[i]->Sumw2();

       // PARTON level

       histo_label_partonpt1 = "PtPARTONLeadingSelection_" + TString::Format("%d",i+1) + "bin";
       histo_title_partonpt1 = "PtPARTONLeadingSelection_" + TString::Format("%d",i+1) + "bin";
       
       PtPARTONLeadingSelection_bins[i] = fs->make<TH1F>(histo_label_partonpt1,histo_title_partonpt1,Ptbins,Ptbinning); 
       PtPARTONLeadingSelection_bins[i]->Sumw2();

       histo_label_partonpt2 = "PtPARTONSubLeadingSelection_" + TString::Format("%d",i+1) + "bin";
       histo_title_partonpt2 = "PtPARTONSubLeadingSelection_" + TString::Format("%d",i+1) + "bin";
       
       PtPARTONSubLeadingSelection_bins[i] = fs->make<TH1F>(histo_label_partonpt2,histo_title_partonpt2,Ptbins,Ptbinning); 
       PtPARTONSubLeadingSelection_bins[i]->Sumw2();


       //Response Matrices
       ResponseLeadingJetPt = "ResponseRecoGenLeadingJetPt_" + TString::Format("%d",i+1) + "bin";
       ResponseSubLeadingJetPt = "ResponseRecoGenSubLeadingJetPt_" + TString::Format("%d",i+1) + "bin";

       ResponseLeadingJetPt_bins[i] = fs->make<TH2F>(ResponseLeadingJetPt,ResponseLeadingJetPt,PtbinsLarge,PtbinningLarge,PtbinsLarge,PtbinningLarge); 
       ResponseSubLeadingJetPt_bins[i] = fs->make<TH2F>(ResponseSubLeadingJetPt,ResponseSubLeadingJetPt,PtbinsLarge,PtbinningLarge,PtbinsLarge,PtbinningLarge); 

       //TUnfold
       ResponseLeadingJetPt = "ResponseRecoGenLeadingJetPtForTUnfold_" + TString::Format("%d",i+1) + "bin";
       ResponseSubLeadingJetPt = "ResponseRecoGenSubLeadingJetPtForTUnfold_" + TString::Format("%d",i+1) + "bin";
       
       ResponseLeadingJetPtForTUnfold_bins[i] = fs->make<TH2F>(ResponseLeadingJetPt,ResponseLeadingJetPt,PtbinsV2,PtbinningV2,PtbinsLarge,PtbinningLarge); 
       ResponseSubLeadingJetPtForTUnfold_bins[i] = fs->make<TH2F>(ResponseSubLeadingJetPt,ResponseSubLeadingJetPt,PtbinsV2,PtbinningV2,PtbinsLarge,PtbinningLarge); 

     }

     DeltaPhiDETSelection = fs->make<TH1F>("DeltaPhiDETSelection","DeltaPhiDETSelection",DeltaPhibins,DeltaPhibinning); DeltaPhiDETSelection->Sumw2();
     DeltaPhiGENSelection = fs->make<TH1F>("DeltaPhiGENSelection","DeltaPhiGENSelection",DeltaPhibins,DeltaPhibinning); DeltaPhiGENSelection->Sumw2();
     DeltaPhiDETSelectionFake = fs->make<TH1F>("DeltaPhiDETSelectionFake","DeltaPhiDETSelectionFake",DeltaPhibinsLarge,DeltaPhibinningLarge); DeltaPhiDETSelectionFake->Sumw2();
     DeltaPhiDETSelectionMiss = fs->make<TH1F>("DeltaPhiDETSelectionMiss","DeltaPhiDETSelectionMiss",DeltaPhibins,DeltaPhibinning); DeltaPhiDETSelectionMiss->Sumw2();
     
     DeltaPhiDETSelectionDoubleBins = fs->make<TH1F>("DeltaPhiDETSelectionDoubleBins","DeltaPhiDETSelectionDoubleBins",DeltaPhibinsLarge,DeltaPhibinningLarge); DeltaPhiDETSelectionDoubleBins->Sumw2();
     DeltaPhiGENSelectionDoubleBins = fs->make<TH1F>("DeltaPhiGENSelectionDoubleBins","DeltaPhiGENSelectionDoubleBins",DeltaPhibinsLarge,DeltaPhibinningLarge); DeltaPhiGENSelectionDoubleBins->Sumw2();
     
     DeltaPhiDETSelectionForToyFits = fs->make<TH1F>("DeltaPhiDETSelectionForToyFits","DeltaPhiDETSelectionForToyFits",30,0,TMath::Pi()); DeltaPhiDETSelectionForToyFits->Sumw2();
     DeltaPhiGENSelectionForToyFits = fs->make<TH1F>("DeltaPhiGENSelectionForToyFits","DeltaPhiGENSelectionForToyFits",30,0,TMath::Pi()); DeltaPhiGENSelectionForToyFits->Sumw2();

     MassGENSubJet0LeadingAfterSelection = fs->make<TH1F>("MassGENSubJet0LeadingAfterSelection","MassGENSubJet0LeadingAfterSelection",200,0.,2000.); MassGENSubJet0LeadingAfterSelection->Sumw2();
     MassGENSubJet0SubLeadingAfterSelection = fs->make<TH1F>("MassGENSubJet0SubLeadingAfterSelection","MassGENSubJet0SubLeadingAfterSelection",200,0.,2000.); MassGENSubJet0SubLeadingAfterSelection->Sumw2();

     MassGENJetLeadingAfterSelection = fs->make<TH1F>("MassGENJetLeadingAfterSelection","MassGENJetLeadingAfterSelection",200,0.,2000.); MassGENJetLeadingAfterSelection->Sumw2();
     MassGENJetSubLeadingAfterSelection = fs->make<TH1F>("MassGENJetSubLeadingAfterSelection","MassGENJetSubLeadingAfterSelection",200,0.,2000.); MassGENJetSubLeadingAfterSelection->Sumw2();

     MassGENJetLead_1Selection = fs->make<TH1F>("MassGENJetLead_1Selection","MassGENJetLead_1Selection",200,0.,2000.); MassGENJetLead_1Selection->Sumw2();
     MassGENJetLead_2Selection = fs->make<TH1F>("MassGENJetLead_2Selection","MassGENJetLead_2Selection",200,0.,2000.); MassGENJetLead_2Selection->Sumw2();

     MassGENJetSubLead_1Selection = fs->make<TH1F>("MassGENJetSubLead_1Selection","MassGENJetSubLead_1Selection",200,0.,2000.); MassGENJetSubLead_1Selection->Sumw2();
     MassGENJetSubLead_2Selection = fs->make<TH1F>("MassGENJetSubLead_2Selection","MassGENJetSubLead_2Selection",200,0.,2000.); MassGENJetSubLead_2Selection->Sumw2();

     MassGENJetLead_TopParton = fs->make<TH1F>("MassGENJetLead_TopParton","MassGENJetLead_TopParton",200,0.,2000.); MassGENJetLead_TopParton->Sumw2();
     MassGENJetLead_SubJet1TopParton = fs->make<TH1F>("MassGENJetLead_SubJet1TopParton","MassGENJetLead_SubJet1TopParton",200,0.,2000.); MassGENJetLead_SubJet1TopParton->Sumw2();

     MassGENJetSubLead_TopParton = fs->make<TH1F>("MassGENJetSubLead_TopParton","MassGENJetSubLead_TopParton",200,0.,2000.); MassGENJetSubLead_TopParton->Sumw2();
     MassGENJetSubLead_SubJet1TopParton = fs->make<TH1F>("MassGENJetSubLead_SubJet1TopParton","MassGENJetSubLead_SubJet1TopParton",200,0.,2000.); MassGENJetSubLead_SubJet1TopParton->Sumw2();

     PtGENJetLead_1Selection = fs->make<TH1F>("PtGENJetLead_1Selection","PtGENJetLead_1Selection",Ptbins,Ptbinning); PtGENJetLead_1Selection->Sumw2();
     PtGENJetLead_2Selection = fs->make<TH1F>("PtGENJetLead_2Selection","PtGENJetLead_2Selection",Ptbins,Ptbinning); PtGENJetLead_2Selection->Sumw2();

     PtGENJetSubLead_1Selection = fs->make<TH1F>("PtGENJetSubLead_1Selection","PtGENJetSubLead_1Selection",Ptbins,Ptbinning); PtGENJetSubLead_1Selection->Sumw2();
     PtGENJetSubLead_2Selection = fs->make<TH1F>("PtGENJetSubLead_2Selection","PtGENJetSubLead_2Selection",Ptbins,Ptbinning); PtGENJetSubLead_2Selection->Sumw2();

     SubJet2MassGENJetLeadingAfterSelection = fs->make<TH1F>("SubJet2MassGENJetLeadingAfterSelection","SubJet2MassGENJetLeadingAfterSelection",200,0.,2000.); SubJet2MassGENJetLeadingAfterSelection->Sumw2();
     SubJet2MassGENJetSubLeadingAfterSelection = fs->make<TH1F>("SubJet2MassGENJetSubLeadingAfterSelection","SubJet2MassGENJetSubLeadingAfterSelection",200,0.,2000.); SubJet2MassGENJetSubLeadingAfterSelection->Sumw2();

     SubJet1MassGENJetLeadingAfterSelection = fs->make<TH1F>("SubJet1MassGENJetLeadingAfterSelection","SubJet1MassGENJetLeadingAfterSelection",200,0.,2000.); SubJet1MassGENJetLeadingAfterSelection->Sumw2();
     SubJet1MassGENJetSubLeadingAfterSelection = fs->make<TH1F>("SubJet1MassGENJetSubLeadingAfterSelection","SubJet1MassGENJetSubLeadingAfterSelection",200,0.,2000.); SubJet1MassGENJetSubLeadingAfterSelection->Sumw2();

     
     DeltaPhiPARTONUnmatched = fs->make<TH1F>("DeltaPhiPARTONUnmatched","DeltaPhiPARTONUnmatched",30,0.0,TMath::Pi()); DeltaPhiPARTONUnmatched->Sumw2();

     DeltaPhiPARTONSelection = fs->make<TH1F>("DeltaPhiPARTONSelection","DeltaPhiPARTONSelection",30,0.0,TMath::Pi()); DeltaPhiPARTONSelection->Sumw2();
     DeltaPhiMatchedSelection = fs->make<TH1F>("DeltaPhiMatchedSelection","DeltaPhiMatchedSelection",DeltaPhibinsLarge,DeltaPhibinningLarge); DeltaPhiMatchedSelection->Sumw2();

     TString histo_label;
     TString histo_title;

     for(int i = 0; i < 6; i++) {
       histo_label = "MassDETLeadingSelection_" + TString::Format("%d",i+1) + "bin";
       histo_title = "MassDETLeadingSelection_" + TString::Format("%d",i+1) + "bin";

       MassDETLeadingSelection_bins[i] = fs->make<TH1F>(histo_label,histo_title,Massbins,Massbinning); 
       MassDETLeadingSelection_bins[i]->Sumw2();
     }

     //file_input.open("../test/weights/RunG-scaled.pred");
     //file_input.open("../test/weights/RunG-SignificanceModel-scaled.pred");
     file_input.open("../test/weights/ModelPredictionsRunGData-DNN.txt");


     if(mIsMCarlo){
       file_PUNominal.open("../test/PileUp2016/ReweightFactors-RunG.txt");
       //file_PUNominal.open("../test/PileUp2016/ReweightFactors-RunG-Herwig.txt");
       
       double factor=-10;
       
       while (file_PUNominal >> factor) {
	 ReweightPU.push_back(factor);
       }

       file_PUNominal.close();	
     }

 } // end of function beginJob()


 //------------------------ endjob() function declaration ---------------------- //
 void Analysis_Template_MC::endJob()
 {

   SDMassDETLeadingSelection->Scale(1.00,"width");
   SDMassDETSubLeadingSelection->Scale(1.00,"width");
   MassDETSubLeadingSelection->Scale(1.00,"width");
   MassDETLeadingSelection->Scale(1.00,"width");
   MassDETSubLeadingSelectionSubJet0->Scale(1.00,"width");
   MassDETSubLeadingSelectionSubJet1->Scale(1.00,"width");
   MassDETSubLeadingSelectionSubJet0->Scale(1.00,"width");
   MassDETSubLeadingSelectionSubJet1->Scale(1.00,"width");

   SDMassDETLeadingPreselection->Scale(1.00,"width");
   SDMassDETSubLeadingPreselection->Scale(1.00,"width");
   MassDETSubLeadingPreselection->Scale(1.00,"width");
   MassDETLeadingPreselection->Scale(1.00,"width");
   MassDETSubLeadingPreselectionSubJet0->Scale(1.00,"width");
   MassDETSubLeadingPreselectionSubJet1->Scale(1.00,"width");
   MassDETSubLeadingPreselectionSubJet0->Scale(1.00,"width");
   MassDETSubLeadingPreselectionSubJet1->Scale(1.00,"width");

   SDMassDETLeadingPrepreselection->Scale(1.00,"width");
   SDMassDETSubLeadingPrepreselection->Scale(1.00,"width");
   MassDETSubLeadingPrepreselection->Scale(1.00,"width");
   MassDETLeadingPrepreselection->Scale(1.00,"width");
   MassDETSubLeadingPrepreselectionSubJet0->Scale(1.00,"width");
   MassDETSubLeadingPrepreselectionSubJet1->Scale(1.00,"width");
   MassDETSubLeadingPrepreselectionSubJet0->Scale(1.00,"width");
   MassDETSubLeadingPrepreselectionSubJet1->Scale(1.00,"width");
   
   for(int i = 0; i < 6; i++) {
     MassDETLeadingSelection_bins[i]->Scale(1.00,"width");
   }
   for(int i = 0; i < 5; i++) {
     PtDETLeadingSelectionForToyFits_bins[i]->Scale(1.00,"width");
     PtDETSubLeadingSelectionForToyFits_bins[i]->Scale(1.00,"width");
   } 

   cout<<"Closing file and saving"<<endl;

   cout<<NEvents<<endl;
   
   mInf->Close();

   file_input.close();

 } // closing endJob()





 //--------------------------- analyze() fuction declaration ------------------ //
void Analysis_Template_MC::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup)
 {

   cout<<" Size = "<<mFileName.size()<<endl;

   for(unsigned ifile=0;ifile<mFileName.size();ifile++){

     mInf = TFile::Open(mFileName[ifile].c_str());
     mDir = (TDirectoryFile*)mInf->Get(mDirName.c_str());
     mTree = (TTree*)mDir->Get(mTreeName.c_str());
     
     unsigned NEntries = mTree->GetEntries();
     cout<<"Reading TREE: "<<NEntries<<" events"<<endl;

     TH1D* NumberEvents=(TH1D*) mInf->Get("boosted/TriggerPass");    
     NEvents=NEvents+NumberEvents->GetBinContent(1);
     NEventsQCD=NumberEvents->GetBinContent(1);
     
     int decade = 0 ;
     
     float hweight=1.;  ///Initial value to one
     double hweightAfterFitToSelection=0.655;  //hard-coded value for top cuetp8m2
     //float hweightAfterFitToSelection=0.855;  //hard-coded value for top herwig

     int lumi_=-100; int nvtx_=-100; int nBJets_=-100; int nJets_=-100; int nGenJets_=-100; int nLeptons_=-100; std::vector<int>* jetNSub_=0; std::vector<int>* jetNBSub_=0; 
     std::vector<int>* jetNGenSub_=0;
     float met_=-100; 
     std::vector<float>* jetPt_=0; std::vector<float>* jetEta_=0; std::vector<float>* jetBtag_=0; std::vector<float>* jetPhi_=0; std::vector<float>* jetUnc_=0; std::vector<float>* jethadronflavor_=0; 
     std::vector<float>* jetEnergy_=0; 
     std::vector<float>* jetMass_=0; std::vector<float>* jetMassSoftDrop_=0;
     
     std::vector<float>* jetchf_=0; std::vector<float>* jetnhf_=0; std::vector<float>* jetphf_=0; std::vector<float>* jetmuf_=0; std::vector<float>* jetelf_=0; std::vector<float>* jettau1_=0; std::vector<float>* jettau2_=0; std::vector<float>* jettau3_=0; 
     
     std::vector<float>* GenJetPt_=0; std::vector<float>* GenJetEta_=0; std::vector<float>* GenJetPhi_=0; std::vector<float>* GenJetenergy_=0; std::vector<float>* GenJetmass_=0; std::vector<float>* GenJetMassSoftDrop_=0;
     
     std::vector<float>* GenJetSoftTau31_=0; std::vector<float>* GenJetSoftTau32_=0;

     std::vector<double>* pre_=0;   std::vector<bool>* bit_=0; std::vector<bool>* metbit_=0;
     
     float metGen_=-100; float metSigGen_=-100; 
     std::vector<bool>* isBJetGen_=0;
     std::vector<bool>* isWJetGen_=0;

     bool matchWGenJet[40];

     int decay_=-100;
     int npu_=-100;
     
     std::vector<float>* partonPt_=0; std::vector<float>* partonEta_=0; std::vector<float>* partonPhi_=0; std::vector<float>* partonEnergy_=0; std::vector<float>* partonId_=0;
     std::vector<float>* partonWPt_=0; std::vector<float>* partonWEta_=0; std::vector<float>* partonWPhi_=0; std::vector<float>* partonWEnergy_=0; std::vector<float>* partonWId_=0;
     
     std::vector<float>* jetbtagSub0_=0;
     std::vector<float>* jetbtagSub1_=0;
     std::vector<float>* jetmassSub0_=0;
     std::vector<float>* jetmassSub1_=0;
     std::vector<float>* jetptSub0_=0;
     std::vector<float>* jetptSub1_=0;
     std::vector<float>* jetetaSub0_=0;
     std::vector<float>* jetetaSub1_=0;
     std::vector<float>* jetphiSub0_=0;
     std::vector<float>* jetphiSub1_=0;
     std::vector<int>* jetflavorSub0_=0;
     std::vector<int>* jetflavorSub1_=0;
     
     std::vector<float>* GenSubJet1Pt_=0;
     std::vector<float>* GenSubJet2Pt_=0;
     std::vector<float>* GenSubJet1Eta_=0;
     std::vector<float>* GenSubJet2Eta_=0;
     std::vector<float>* GenSubJet1Phi_=0;
     std::vector<float>* GenSubJet2Phi_=0;
     std::vector<float>* GenSubJet1Mass_=0;
     std::vector<float>* GenSubJet2Mass_=0;
     std::vector<float>* GenSubJetsDeltaR_=0;
     std::vector<float>* GenSubJetsMu_=0;

     std::vector<float>* pdfWeights_=0;
     std::vector<float>* scaleWeights_=0;
     
     mTree->SetBranchAddress("lumi",&lumi_);
     mTree->SetBranchAddress("nvtx",&nvtx_);
     mTree->SetBranchAddress("met",&met_);
     mTree->SetBranchAddress("nBJets",&nBJets_);
     mTree->SetBranchAddress("nJets",&nJets_);
     mTree->SetBranchAddress("nLeptons",&nLeptons_);
     mTree->SetBranchAddress("jetNSub",&jetNSub_);
     mTree->SetBranchAddress("jetNBSub",&jetNBSub_);
     mTree->SetBranchAddress("jetPt",&jetPt_);
     mTree->SetBranchAddress("jetBtag",&jetBtag_);
     mTree->SetBranchAddress("jetEta",&jetEta_);
     mTree->SetBranchAddress("jetPhi",&jetPhi_);
     mTree->SetBranchAddress("jetMass",&jetMass_);
     mTree->SetBranchAddress("jetEnergy",&jetEnergy_);
     mTree->SetBranchAddress("jetMassSoftDrop",&jetMassSoftDrop_);
     
     mTree->SetBranchAddress("jetChf",&jetchf_);
     mTree->SetBranchAddress("jetNhf",&jetnhf_);
     mTree->SetBranchAddress("jetPhf",&jetphf_);
     mTree->SetBranchAddress("jetMuf",&jetmuf_);
     mTree->SetBranchAddress("jetElf",&jetelf_);
     mTree->SetBranchAddress("jetTau1",&jettau1_);
     mTree->SetBranchAddress("jetTau2",&jettau2_);
     mTree->SetBranchAddress("jetTau3",&jettau3_);
     
     mTree->SetBranchAddress("jetBtagSub0",&jetbtagSub0_);
     mTree->SetBranchAddress("jetBtagSub1",&jetbtagSub1_);
     mTree->SetBranchAddress("jetMassSub0",&jetmassSub0_);
     mTree->SetBranchAddress("jetMassSub1",&jetmassSub1_);
     mTree->SetBranchAddress("jetPtSub0",&jetptSub0_);
     mTree->SetBranchAddress("jetPtSub1",&jetptSub1_);
     mTree->SetBranchAddress("jetPhiSub0",&jetphiSub0_);
     mTree->SetBranchAddress("jetPhiSub1",&jetphiSub1_);
     mTree->SetBranchAddress("jetEtaSub0",&jetetaSub0_);
     mTree->SetBranchAddress("jetEtaSub1",&jetetaSub1_);
          
     if(mIsMCarlo){
       mTree->SetBranchAddress("jetFlavorSub0",&jetflavorSub0_);
       mTree->SetBranchAddress("jetFlavorSub1",&jetflavorSub1_);

       mTree->SetBranchAddress("jetFlavorHadron",&jethadronflavor_);
       mTree->SetBranchAddress("jetNSubGen",&jetNGenSub_);
       mTree->SetBranchAddress("nGenJets",&nGenJets_);
       mTree->SetBranchAddress("GenJetpt",&GenJetPt_);
       mTree->SetBranchAddress("GenJeteta",&GenJetEta_);
       mTree->SetBranchAddress("GenJetphi",&GenJetPhi_);
       mTree->SetBranchAddress("GenJetmass",&GenJetmass_);
       mTree->SetBranchAddress("GenJetenergy",&GenJetenergy_);
       mTree->SetBranchAddress("GenSoftDropMass",&GenJetMassSoftDrop_);
       mTree->SetBranchAddress("GenSoftDropTau32",&GenJetSoftTau32_);
       mTree->SetBranchAddress("GenSoftDropTau31",&GenJetSoftTau31_);
       mTree->SetBranchAddress("metGen",&metGen_);
       mTree->SetBranchAddress("decay",&decay_);
       mTree->SetBranchAddress("metGenSig",&metSigGen_);
       mTree->SetBranchAddress("isBJetGen",&isBJetGen_);
       
       mTree->SetBranchAddress("npu",&npu_);

       mTree->SetBranchAddress("partonPt",&partonPt_);
       mTree->SetBranchAddress("partonEta",&partonEta_);
       mTree->SetBranchAddress("partonPhi",&partonPhi_);
       mTree->SetBranchAddress("partonE",&partonEnergy_);
       mTree->SetBranchAddress("partonId",&partonId_);

       if(mMCtype=="TOP" || mMCtype=="TOPHerwig" || mMCtype=="TOP-Old"){
	 mTree->SetBranchAddress("WBosonPt",&partonWPt_);
	 mTree->SetBranchAddress("WBosonEta",&partonWEta_);
	 mTree->SetBranchAddress("WBosonPhi",&partonWPhi_);
	 mTree->SetBranchAddress("WBosonE",&partonWEnergy_);
	 mTree->SetBranchAddress("WBosonId",&partonWId_);
	 mTree->SetBranchAddress("isWJetGen",&isWJetGen_);
       }

       if(mMCtype=="TOP"){
	 mTree->SetBranchAddress("jetUnc",&jetUnc_);

	 mTree->SetBranchAddress("pdfWeights",&pdfWeights_);
	 mTree->SetBranchAddress("scaleWeights",&scaleWeights_);
       }

       mTree->SetBranchAddress("GenSubJet1Pt",&GenSubJet1Pt_);
       mTree->SetBranchAddress("GenSubJet2Pt",&GenSubJet2Pt_);
       mTree->SetBranchAddress("GenSubJet1Eta",&GenSubJet1Eta_);
       mTree->SetBranchAddress("GenSubJet2Eta",&GenSubJet2Eta_);
       mTree->SetBranchAddress("GenSubJet1Phi",&GenSubJet1Phi_);
       mTree->SetBranchAddress("GenSubJet2Phi",&GenSubJet2Phi_);
       mTree->SetBranchAddress("GenSubJet1Mass",&GenSubJet1Mass_);
       mTree->SetBranchAddress("GenSubJet2Mass",&GenSubJet2Mass_);
       mTree->SetBranchAddress("GenSubJetsDeltaR",&GenSubJetsDeltaR_);
       mTree->SetBranchAddress("GenSubJetsMu",&GenSubJetsMu_);
	 
     }
     
     if(!mIsMCarlo){
       mTree->SetBranchAddress("triggerPre",&pre_);
       mTree->SetBranchAddress("triggerBit",&bit_);
       mTree->SetBranchAddress("metBit",&metbit_);
     }
     
     //Trigger infos
     TH1F *hTrigNames = (TH1F*)mInf->Get("boosted/TriggerNames");

     TString HLTJet[10];
     int ihltj[10];
     
     for (unsigned i=0; i<mTriggers.size(); i++){
       HLTJet[i] = mTriggers[i];
     }
     
     for (int i=0; i<hTrigNames->GetNbinsX(); i++){
       ihltj[i] = -1 ;
     }

     cout<<"Finding trigger mapping: "<<endl;
     //----------- loop over the X-axis labels -----------------
     if(!mIsMCarlo){
       for(int ibin=0;ibin<hTrigNames->GetNbinsX();ibin++) {
	 TString ss(hTrigNames->GetXaxis()->GetBinLabel(ibin+1));
	 for (unsigned ii=0; ii<mTriggers.size(); ii++){
	   if (ss == HLTJet[ii]) {
	     ihltj[ii] = ibin;
	     continue;
	   }
	 } 
       }
       
       for (unsigned ij=0; ij<mTriggers.size(); ij++){
	 if (ihltj[ij] == -1) {
	   cout<<"The requested trigger ("<<HLTJet[ij]<<") is not found "<<endl;
	 }
	 else {
	   cout<<HLTJet[ij]<<" --> "<<ihltj[ij]<<endl;
	 }
       }
     }
     
     for(unsigned  l=0; l<NEntries; l++) {
       
       double mvaD=-100;

       double SmearingFactor[2]; SmearingFactor[0] = 1.; SmearingFactor[1]=1.;
       
       //----------- progress report -------------
       double progress = 10.0*l/(1.0*NEntries);
       int k = TMath::FloorNint(progress);
       if (k > decade)
	 cout<<10*k<<" %"<<endl;
       decade = k;
       
       //----------- read the event --------------

       mTree->GetEntry(l);

       if(nGenJets_>1){
	 dijetsbefore++;
       }
       
       bool triggerfired(false);
       
       if(mIsMCarlo){
	 triggerfired=true;
	 hweight=mCrossSection/NEvents;
	 	 
	 if(mMCtype=="QCD"){	 
	   NEventsQCD=NumberEvents->GetBinContent(1);
	   if(ifile==0) hweight=31630./NEventsQCD;
	   if(ifile==1) hweight=6802./NEventsQCD;
	   if(ifile==2) hweight=1206./NEventsQCD;
	   if(ifile==3) hweight=120.4/NEventsQCD;
	   if(ifile==4) hweight=25.25/NEventsQCD;
	 }

	 if(mMCtype=="TOP" || mMCtype=="TOP-Old"){
	   hweight=1.;//TOP events
	   hweightAfterFitToSelection=0.655;//TOP events
	   if(musePDFWeights){
	     hweight=hweight*pdfWeights_->at(mnumbPDFWeights);
	   }

	   if(museScaleWeights){
	     hweight=hweight*scaleWeights_->at(mnumbScaleWeights);
	   }

	   //for(unsigned int i=0;i<scaleWeights_->size(); i++){
	   //cout<<scaleWeights_->at(i)<<" scale i"<<endl;
	   //cout<<pdfWeights_->size()<<" pdf size"<<endl;
	   //}
	 }

	 if(mMCtype=="TOPHerwig"){
	   hweight=1.;//TOP events
	   hweightAfterFitToSelection=0.855;//TOP events
	 }	 

	 PileUpInteractions->Fill(npu_,hweight);
	 hweight=hweight*ReweightPU.at(npu_);
	 hweightAfterFitToSelection=hweightAfterFitToSelection*ReweightPU.at(npu_);	 
       }
       
       if(!mIsMCarlo){
      
	 double prescale=1;
	 bool metfilters=true;

	 for(unsigned imet=0;imet<metbit_->size();imet++){
	  metfilters=metfilters || metbit_->at(imet);
	 }	 

	 if(mMCtype!="RUN") {
	   if(metfilters && bit_->at(3) && jetPt_->at(0)>=400 && jetPt_->at(0)<550) {triggerfired=true; prescale=pre_->at(3);}//3 is the trigger of interest, 4 is the trigger 140 used for tests
	   //if(metfilters && bit_->at(5) && jetPt_->at(0)>=400 && jetPt_->at(0)<550) {triggerfired=true; prescale=pre_->at(5);}//for control region in data-driven (320 single jet trigger)
	   if(metfilters && bit_->at(7) && jetPt_->at(0)>=550) {triggerfired=true; prescale=pre_->at(7);}
	 }

	 if(mMCtype=="RUN") {
	   //if(metfilters && bit_->at(0) && jetPt_->at(0)>=400 && jetPt_->at(0)<550) {triggerfired=true; prescale=pre_->at(0);}
	   if(metfilters && bit_->at(2) && jetPt_->at(0)>=400 && jetPt_->at(0)<550) {triggerfired=true; prescale=pre_->at(2);}
	   if(metfilters && bit_->at(1) && jetPt_->at(0)>=550) {triggerfired=true; prescale=pre_->at(1);}
	 }

	 hweight=prescale/mIntLumi;
       }      
       
       if(triggerfired==false) continue;
       
       int NJetsPreselection=0;
       bool prepreselection=false;
       bool selection=false;
     
       for(int j=0; j< nJets_; j++){

	 double ptaftersmear=jetPt_->at(j);

	 if(mIsMCarlo){
	   for(int i=0; i< nGenJets_; i++){
	     ptaftersmear=PtSmeared(fabs(jetEta_->at(j)),jetPhi_->at(j),jetPt_->at(j),fabs(GenJetEta_->at(i)),GenJetPhi_->at(i),GenJetPt_->at(i));
	   }
	 }
	 
	 if(mShift==0.){
	   if(ptaftersmear>400 && fabs(jetEta_->at(j))<2.4){//mass drop+b-tag
	     NJetsPreselection++;
	   }
	 }

	 if(mShift!=0.){
	   if((ptaftersmear*(1+mShift*jetUnc_->at(j)))>400 && fabs(jetEta_->at(j))<2.4 && jetMassSoftDrop_->at(j)>50){//mass drop+b-tag
	     NJetsPreselection++;
	   }
	 }
       }
       
       if(NJetsPreselection>1) prepreselection=true;

       if(prepreselection){

	 //determining the smearing due to JER for leading and subleading
	 if(mIsMCarlo){
	   double smearLeading=jetPt_->at(0); double smearSubLeading=jetPt_->at(1);
	   for(int j=0; j< nGenJets_; j++){
	     smearLeading=PtSmeared(fabs(jetEta_->at(0)),jetPhi_->at(0),jetPt_->at(0),fabs(GenJetEta_->at(j)),GenJetPhi_->at(j),GenJetPt_->at(j));
	     smearSubLeading=PtSmeared(fabs(jetEta_->at(1)),jetPhi_->at(1),jetPt_->at(1),fabs(GenJetEta_->at(j)),GenJetPhi_->at(j),GenJetPt_->at(j));
	   }
	   
	   SmearingFactor[0]=smearLeading/jetPt_->at(0);
	   SmearingFactor[1]=smearSubLeading/jetPt_->at(1);
	 }

	 double tau31_1jet=-1;
         if(jettau1_->at(0)!=0){ tau31_1jet= jettau3_->at(0)/jettau1_->at(0);}
         double tau32_1jet=-1;
         if(jettau2_->at(0)!=0){ tau32_1jet= jettau3_->at(0)/jettau2_->at(0);}
	 
	 double tau31_2jet=-1;
	 if(jettau1_->at(1)!=0){ tau31_2jet= jettau3_->at(1)/jettau1_->at(1);}
         double tau32_2jet=-1;
         if(jettau2_->at(1)!=0){ tau32_2jet= jettau3_->at(1)/jettau2_->at(1);}
	 
	 double DeltaRSubJet_1jet=DeltaR(jetetaSub0_->at(0),jetphiSub0_->at(0),jetetaSub1_->at(0),jetphiSub1_->at(0));
	 double DeltaRSubJet_2jet=DeltaR(jetetaSub0_->at(1),jetphiSub0_->at(1),jetetaSub1_->at(1),jetphiSub1_->at(1));
      
	 //cout<<"QCD:-1 TOP:1 "<<jetMass_->at(0)<<" "<<jetmassSub0_->at(0)<<" "<<jetmassSub1_->at(0)<<" "<<tau31_1jet<<" "<<tau32_1jet<<" "<<jetptSub0_->at(0)<<" "<<jetptSub1_->at(0)<<" "<<DeltaRSubJet_1jet<<" "<<jetMass_->at(1)<<" "<<jetmassSub0_->at(1)<<" "<<jetmassSub1_->at(1)<<" "<<tau31_2jet<<" "<<tau32_2jet<<" "<<jetptSub0_->at(1)<<" "<<jetptSub1_->at(1)<<" "<<DeltaRSubJet_2jet<<endl;
	 
	 //cout<<"QCD:-1 TOP:1 1:"<<tau31_1jet<<" 2:"<<tau32_1jet<<" 3:"<<tau31_2jet<<" 4:"<<tau32_2jet<<" 5:"<<jetMass_->at(0)<<" 6:"<<jetmassSub0_->at(0)<<" 7:"<<jetMass_->at(1)<<" 8:"<<jetmassSub0_->at(1)<<endl;
	 //cout<<"QCD:-1 TOP:1 1:"<<jettau1_->at(0)<<" "<<jettau2_->at(0)<<" "<<jettau3_->at(0)<<" "<<jettau1_->at(1)<<" "<<jettau2_->at(1)<<" "<<jettau3_->at(1)<<endl;;
	 
	 //double tau31Leading = tau31_1jet;
	 //double tau32Leading = tau32_1jet;
	 //double MassLeading = jetMass_->at(0);
	 //double SDmassLeading = jetMassSoftDrop_->at(0);
	 //double MassSubJet0Leading = jetmassSub0_->at(0);
	 //double MassSubJet1Leading = jetmassSub0_->at(1);
	 //double jetSubJetpt0Leading = jetptSub0_->at(0);
	 //double jetSubJetpt1Leading = jetptSub0_->at(1);
	 
	 //double tau31SubLeading = tau31_2jet;
	 //double tau32SubLeading = tau32_2jet;
	 //double MassSubLeading = jetMass_->at(1);
	 //double SDmassSubLeading = jetMassSoftDrop_->at(0);
	 //double MassSubJet0SubLeading = jetmassSub1_->at(0);
	 //double MassSubJet1SubLeading = jetmassSub1_->at(1);
	 //double jetSubJetpt0SubLeading = jetptSub1_->at(0);
	 //double jetSubJetpt1SubLeading = jetptSub1_->at(1);

	 //mvaD = discr_->eval(tau31Leading,tau32Leading,MassLeading,MassSubJet0Leading,tau31SubLeading,tau32SubLeading,MassSubLeading,MassSubJet0SubLeading); //8variables
	 mvaD = discr_->eval(jettau1_->at(0),jettau2_->at(0),jettau3_->at(0),jettau1_->at(1),jettau2_->at(1),jettau3_->at(1));
	 //mvaD = discr_->eval(MassLeading,MassSubJet0Leading,MassSubLeading,MassSubJet0SubLeading); // Daniela
	 

	 MVADesy->Fill(mvaD,hweight);

	 //Add histograms
	 PtDETLeadingPrepreselectionSubJet0->Fill(jetptSub0_->at(0),hweight);
	 PtDETLeadingPrepreselectionSubJet1->Fill(jetptSub1_->at(0),hweight);
	 EtaDETLeadingPrepreselectionSubJet0->Fill(jetetaSub0_->at(0),hweight);
	 EtaDETLeadingPrepreselectionSubJet1->Fill(jetetaSub1_->at(0),hweight);
	 MassDETLeadingPrepreselectionSubJet0->Fill(jetmassSub0_->at(0),hweight);
	 MassDETLeadingPrepreselectionSubJet1->Fill(jetmassSub1_->at(0),hweight);
	 PhiDETLeadingPrepreselectionSubJet0->Fill(jetphiSub0_->at(0),hweight);
	 PhiDETLeadingPrepreselectionSubJet1->Fill(jetphiSub1_->at(0),hweight);
	 DeltaRDETLeadingPrepreselectionSubJets->Fill(DeltaRSubJet_1jet,hweight);
	 
	 Tau1DETLeadingJetPrepreselection->Fill(jettau1_->at(0),hweight);
	 Tau2DETLeadingJetPrepreselection->Fill(jettau2_->at(0),hweight);
	 Tau3DETLeadingJetPrepreselection->Fill(jettau3_->at(0),hweight);
	 Tau31DETLeadingJetPrepreselection->Fill(tau31_1jet,hweight);
	 Tau32DETLeadingJetPrepreselection->Fill(tau32_1jet,hweight);

	 PtDETSubLeadingPrepreselectionSubJet0->Fill(jetptSub0_->at(1),hweight);
	 PtDETSubLeadingPrepreselectionSubJet1->Fill(jetptSub1_->at(1),hweight);
	 EtaDETSubLeadingPrepreselectionSubJet0->Fill(jetetaSub0_->at(1),hweight);
	 EtaDETSubLeadingPrepreselectionSubJet1->Fill(jetetaSub1_->at(1),hweight);
	 MassDETSubLeadingPrepreselectionSubJet0->Fill(jetmassSub0_->at(1),hweight);
	 MassDETSubLeadingPrepreselectionSubJet1->Fill(jetmassSub1_->at(1),hweight);
	 PhiDETSubLeadingPrepreselectionSubJet0->Fill(jetphiSub0_->at(1),hweight);
	 PhiDETSubLeadingPrepreselectionSubJet1->Fill(jetphiSub1_->at(1),hweight);
	 DeltaRDETSubLeadingPrepreselectionSubJets->Fill(DeltaRSubJet_2jet,hweight);
	 
	 Tau1DETSubLeadingJetPrepreselection->Fill(jettau1_->at(1),hweight);
	 Tau2DETSubLeadingJetPrepreselection->Fill(jettau2_->at(1),hweight);
	 Tau3DETSubLeadingJetPrepreselection->Fill(jettau3_->at(1),hweight);
	 Tau31DETSubLeadingJetPrepreselection->Fill(tau31_2jet,hweight);
	 Tau32DETSubLeadingJetPrepreselection->Fill(tau32_2jet,hweight);
	 
	 PtDETLeadingPrepreselection->Fill(SmearingFactor[0]*jetPt_->at(0),hweight);
	 PtDETSubLeadingPrepreselection->Fill(SmearingFactor[1]*jetPt_->at(1),hweight);
	 EtaDETLeadingPrepreselection->Fill(jetEta_->at(0),hweight);
	 EtaDETSubLeadingPrepreselection->Fill(jetEta_->at(1),hweight);
	 MassDETLeadingPrepreselection->Fill(jetMass_->at(0),hweight);
	 MassDETSubLeadingPrepreselection->Fill(jetMass_->at(1),hweight);
	 PhiDETLeadingPrepreselection->Fill(jetPhi_->at(0),hweight);
	 PhiDETSubLeadingPrepreselection->Fill(jetPhi_->at(1),hweight);
	 SDMassDETLeadingPrepreselection->Fill(jetMassSoftDrop_->at(0),hweight);
	 SDMassDETSubLeadingPrepreselection->Fill(jetMassSoftDrop_->at(1),hweight);

	 num_of_Vtx->Fill(nvtx_,hweight);
	 num_of_Jets->Fill(nJets_,hweight);
	 
	 metDET->Fill(met_,hweight);
	 
	 if(mIsMCarlo){ num_of_GenJets->Fill(nGenJets_,hweight); 
	   metGEN->Fill(metGen_,hweight);
	   metSigGEN->Fill(metSigGen_,hweight);
	   FlavourDETPrepreselectionLeadingJet->Fill(jethadronflavor_->at(0),hweight);
	   FlavourDETPrepreselectionSubLeadingJet->Fill(jethadronflavor_->at(1),hweight);
	 
	   if(jethadronflavor_->at(0)==5) {
	     PtDETLeadingPrepreselectionBTrue->Fill(SmearingFactor[0]*jetPt_->at(0),hweight);
	     if(jethadronflavor_->at(1)==5) {
	       PtDETSubLeadingPrepreselectionBTrue->Fill(SmearingFactor[1]*jetPt_->at(1),hweight);
	     }
	   }
	 }
	 
	 num_of_BJets->Fill(nBJets_,hweight);
	 num_of_Leptons->Fill(nLeptons_,hweight);
	 
	 for(int j=0; j< nJets_; j++){
	   
	   
	   if(jetPt_->at(j)>400 && fabs(jetEta_->at(j))<2.4){
	     
	     num_of_SubJets->Fill(jetNSub_->at(j),hweight);
	     num_of_SubBJets->Fill(jetNBSub_->at(j),hweight);
	     
	     ptDETJet->Fill(jetPt_->at(j),hweight);
	     yDETJet->Fill(jetEta_->at(j),hweight);
	     phiDETJet->Fill(jetPhi_->at(j),hweight);
	     EnergyDETJet->Fill(jetEnergy_->at(j),hweight);
	     MassDETJet->Fill(jetMass_->at(j),hweight);
	     MassSoftDropDETJet->Fill(jetMassSoftDrop_->at(j),hweight);
	     BtagDETJet->Fill(jetBtag_->at(j),hweight);
	     Tau1DETJet->Fill(jettau1_->at(j),hweight);
	     Tau2DETJet->Fill(jettau2_->at(j),hweight);
	     Tau3DETJet->Fill(jettau3_->at(j),hweight);
	     double tau31=-1;
	     if(jettau1_->at(j)!=0){ tau31= jettau3_->at(j)/jettau1_->at(j);}
	     double tau32=-1;
	     if(jettau2_->at(j)!=0){ tau32= jettau3_->at(j)/jettau2_->at(j);}
	     Tau31DETJet->Fill(tau31,hweight);
	     Tau32DETJet->Fill(tau32,hweight);
	     ChfDETJet->Fill(jetchf_->at(j),hweight);
	     NhfDETJet->Fill(jetnhf_->at(j),hweight);
	     PhfDETJet->Fill(jetphf_->at(j),hweight);
	     MufDETJet->Fill(jetmuf_->at(j),hweight);
	     ElfDETJet->Fill(jetelf_->at(j),hweight);
     
	     //subjet information
	     BtagDETSubJet0->Fill(jetbtagSub0_->at(j),hweight);
	     BtagDETSubJet1->Fill(jetbtagSub1_->at(j),hweight);
	     MassDETSubJet0->Fill(jetmassSub0_->at(j),hweight);
	     MassDETSubJet1->Fill(jetmassSub1_->at(j),hweight);
	     PtDETSubJet0->Fill(jetptSub0_->at(j),hweight);
	     PtDETSubJet1->Fill(jetptSub1_->at(j),hweight);
	     EtaDETSubJet0->Fill(jetetaSub0_->at(j),hweight);
	     EtaDETSubJet1->Fill(jetetaSub1_->at(j),hweight);
	     EtaDETSubJet0->Fill(jetphiSub0_->at(j),hweight);
	     EtaDETSubJet1->Fill(jetphiSub1_->at(j),hweight);
	     if(mIsMCarlo) FlavourDETSubJet0->Fill(jetflavorSub0_->at(j),hweight);
	     if(mIsMCarlo) FlavourDETSubJet1->Fill(jetflavorSub1_->at(j),hweight);
	     
	     if(mIsMCarlo) {
	       if(jetbtagSub0_->at(j)>0.8484){
		 TaggedAndBJetDETSubJet0DEN->Fill(jetptSub0_->at(j),hweight);
		 if(fabs(jetflavorSub0_->at(j))==5){
		   TaggedAndBJetDETSubJet0NUM->Fill(jetptSub0_->at(j),hweight);
		 }
	       }

	       if(jetbtagSub1_->at(j)>0.8484){
		 TaggedAndBJetDETSubJet1DEN->Fill(jetptSub1_->at(j),hweight);
		 if(fabs(jetflavorSub1_->at(j))==5){
		   TaggedAndBJetDETSubJet1NUM->Fill(jetptSub1_->at(j),hweight);
		 }
	       }
	     
	       if(fabs(jetflavorSub0_->at(j))==5){
		 BJetAndTaggedDETSubJet0DEN->Fill(jetptSub0_->at(j),hweight);
		 
		 if(jetbtagSub0_->at(j)>0.8484){
		   BJetAndTaggedDETSubJet0NUM->Fill(jetptSub0_->at(j),hweight);
		 }
	       }
	       
	       if(fabs(jetflavorSub1_->at(j))==5){
		 BJetAndTaggedDETSubJet1DEN->Fill(jetptSub1_->at(j),hweight);
		 if(jetbtagSub1_->at(j)>0.8484){	   
		   BJetAndTaggedDETSubJet1NUM->Fill(jetptSub1_->at(j),hweight);
		 }
	       }
	     }

	     if(false){
	       printf("Number of PFJets=%i\n",nJets_);
	       printf("j=%2i  pt=%8.3f  y=%6.3f  phi=%6.3f\n",j,jetPt_->at(j),jetEta_->at(j),jetPhi_->at(j));
	     }
	     
	   }
	 }
       }
       
       bool genpreselection=false;
       bool partonpreselection=false;

       if(mIsMCarlo){
	 //Gen jet part
	 
	 for(int j=0; j< nGenJets_; j++){

	   ptGENJet->Fill(GenJetPt_->at(j),hweight);
	   yGENJet->Fill(GenJetEta_->at(j),hweight);
	   phiGENJet->Fill(GenJetPhi_->at(j),hweight);
	   EnergyGENJet->Fill(GenJetenergy_->at(j),hweight);
	   MassGENJet->Fill(GenJetmass_->at(j),hweight);
	   MassSoftDropGENJet->Fill(GenJetMassSoftDrop_->at(j),hweight);
	   BJetGenFlag->Fill(isBJetGen_->at(j),hweight);
	   if(mMCtype=="TOP" || mMCtype=="SingleTOP" || mMCtype=="TOPHerwig" || mMCtype=="TOP-Old") WJetGenFlag->Fill(isWJetGen_->at(j),hweight);
	   
	   //Matching with W-jet

	   double deltaR_jWMin=10;
	   if(mMCtype=="TOP" || mMCtype=="TOPHerwig" || mMCtype=="TOP-Old"){
	     for(unsigned i=0; i<partonWPt_->size();i++){
	       
	       double deltaR_jW=DeltaR(partonWEta_->at(i),partonWPhi_->at(i),GenJetEta_->at(j),GenJetPhi_->at(j));
	       if(deltaR_jW<deltaR_jWMin){
		 deltaR_jWMin=deltaR_jW;
	       }
	     }
	   }

	   DeltaR_GenJetW->Fill(deltaR_jWMin,hweight);
	   
	   if(deltaR_jWMin<0.4) matchWGenJet[j]=true;
	   if(deltaR_jWMin>=0.4) matchWGenJet[j]=false;
		       	 	    
	   Tau32GENJet->Fill(GenJetSoftTau32_->at(j),hweight);
	   Tau31GENJet->Fill(GenJetSoftTau31_->at(j),hweight);

	   if(jetNGenSub_->at(j)>1){
	     PtGENSubJet0->Fill(GenSubJet1Pt_->at(j),hweight);
	     PtGENSubJet1->Fill(GenSubJet2Pt_->at(j),hweight);
	     EtaGENSubJet0->Fill(GenSubJet1Eta_->at(j),hweight);
	     EtaGENSubJet1->Fill(GenSubJet2Eta_->at(j),hweight);
	     PhiGENSubJet0->Fill(GenSubJet1Phi_->at(j),hweight);
	     PhiGENSubJet1->Fill(GenSubJet2Phi_->at(j),hweight);
	     MassGENSubJet0->Fill(GenSubJet1Mass_->at(j),hweight);
	     MassGENSubJet1->Fill(GenSubJet2Mass_->at(j),hweight);
	   }

	   DeltaRGENSubJets->Fill(GenSubJetsDeltaR_->at(j),hweight);
	   
	 }

	 //Preselection at GEN level
	 bool toptagjets=false;
	 if(nGenJets_>0){
	   if(GenJetPt_->at(0)>950 && fabs(GenJetEta_->at(0))<2.4){// && GenJetMassSoftDrop_->at(0)>50){//mass drop+b-tag
	     ptGENJetLeading->Fill(GenJetPt_->at(0),hweight);
	   }
	 }	     
       
	 if(nGenJets_>1){
	   dijetsafter++;
	   if(GenJetPt_->at(0)>400 && fabs(GenJetEta_->at(0))<2.4 && GenJetMassSoftDrop_->at(0)>50){//mass drop+b-tag
	     if(GenJetPt_->at(1)>400 && fabs(GenJetEta_->at(1))<2.4 && GenJetMassSoftDrop_->at(1)>50){//mass drop+b-tag
	       MassGENJetLead_1Selection->Fill(GenJetMassSoftDrop_->at(0),hweight);
	       MassGENJetSubLead_1Selection->Fill(GenJetMassSoftDrop_->at(1),hweight);

	       PtGENJetLead_1Selection->Fill(GenJetPt_->at(0),hweight);
	       PtGENJetSubLead_1Selection->Fill(GenJetPt_->at(1),hweight);
	       MassGENJetSubLead_1Selection->Fill(GenJetMassSoftDrop_->at(1),hweight);
	       //if(GenSubJet1Mass_->at(0) > 70 && isBJetGen_->at(0)){ //hadron level def
	       //if(matchWGenJet[0] && isBJetGen_->at(0)){ //hadron level def
	       if(isBJetGen_->at(0)){ //hadron level def
		 PtGENJetLead_2Selection->Fill(GenJetPt_->at(0),hweight);
		 //if(GenSubJet1Mass_->at(1) > 70 && isBJetGen_->at(1)){ //hadron level def
		 //if(matchWGenJet[1] && isBJetGen_->at(1)){ //hadron level def
		 if(isBJetGen_->at(1)){ //hadron level def
		   PtGENJetSubLead_2Selection->Fill(GenJetPt_->at(1),hweight);
		   MassGENJetLead_2Selection->Fill(GenJetMassSoftDrop_->at(0),hweight);
		   MassGENJetSubLead_2Selection->Fill(GenJetMassSoftDrop_->at(1),hweight);
		   toptagjets=true;
		 }
	       }
	     }
	   }
	 }
       
	 if(toptagjets) {
	   genpreselection=true;	     
	 }

	 //preselection at PARTON level
	 bool toppartons=false;
	 if(partonPt_->size()>1 && decay_==0){
	   if(partonPt_->size()>2) cout<<"Parton size greater than 2"<<endl;
	   //if(partonPt_->at(0)>400 && fabs(partonEta_->at(0))<2.4){
	   //if(partonPt_->at(1)>400 && fabs(partonEta_->at(1))<2.4){
	   if(partonPt_->at(0)>400){
	     if(partonPt_->at(1)>400){
	       toppartons=true;
	     }
	   } 
	 }
	 
	 if(toppartons) partonpreselection=true;	     


	 if(nGenJets_>1 && (mMCtype=="TOP" || mMCtype=="TOPHerwig" || mMCtype=="TOP-Old")){

	   bool genleadingmatched=false;
	   
	   if(matchWGenJet[0] && isBJetGen_->at(0)) {
	     PtMatchedGenJetAndBWJetLeading->Fill(GenJetPt_->at(0),hweight);
	   
	     for(int i=0; i< nJets_; i++){
	       if(jetPt_->at(i)>400){
		 double deltaR_jt=DeltaR(GenJetEta_->at(0),GenJetPhi_->at(0),jetEta_->at(i),jetPhi_->at(i));
		 if(deltaR_jt<0.4){
		   genleadingmatched=true;
		   PtMatchedGenJetAndBWJetLeadingReco->Fill(GenJetPt_->at(0),hweight);
		 }
	       }
	     }
	   }
	 
	   if(matchWGenJet[1] && isBJetGen_->at(1) && genleadingmatched) {
	     PtMatchedGenJetAndBWJetSubLeading->Fill(GenJetPt_->at(1),hweight);
	   
	     for(int i=0; i< nJets_; i++){
	       if(jetPt_->at(i)>400){
		 double deltaR_jt=DeltaR(GenJetEta_->at(1),GenJetPhi_->at(1),jetEta_->at(i),jetPhi_->at(i));
		 if(deltaR_jt<0.4){
		   PtMatchedGenJetAndBWJetSubLeadingReco->Fill(GenJetPt_->at(1),hweight);
		 }
	       }
	     }
	   }
	 }
       
	 //Match with partons
	 if((mMCtype=="TOP" || mMCtype=="TOPHerwig" || mMCtype=="TOP-Old") && decay_==0){

	   //if(nGenJets_>1 && partonPt_->at(0)>400 && partonPt_->at(1)>400){ //add here more requirements!! 
	   if(nGenJets_>1){ //add here more requirements!!
	     //PtLeadingGenJet->Fill(GenJetPt_->at(0),hweight);			     
	     //PtSubLeadingGenJet->Fill(GenJetPt_->at(1),hweight);	
	     //PtLeadingGenJet->Fill(partonPt_->at(0),hweight);			     
	     //PtSubLeadingGenJet->Fill(partonPt_->at(1),hweight);	

	     if(GenJetPt_->at(0)>400 && fabs(GenJetEta_->at(0))<2.4 && GenJetMassSoftDrop_->at(0)>50 && isBJetGen_->at(0)){//massdrop+b-tag

	     //if(GenJetPt_->at(0)>400 && fabs(GenJetEta_->at(0))<2.4){//mass drop+b-tag
	       if(GenJetPt_->at(1)>400 && fabs(GenJetEta_->at(1))<2.4 && GenJetMassSoftDrop_->at(1)>50 && isBJetGen_->at(1)){//mass drop+b-tag
		 //if(GenJetPt_->at(1)>400 && fabs(GenJetEta_->at(1))<2.4){//mass drop+b-tag
		 
		 PtLeadingGenJet->Fill(GenJetPt_->at(0),hweight);			     
		 PtSubLeadingGenJet->Fill(GenJetPt_->at(1),hweight);	

		 double deltaR_jt_1_min=10;
		 double deltaR_jt_2_min=10;
		 
		 for(unsigned j=0; j<partonPt_->size();j++){
		   double deltaR_jt_1=DeltaR(partonEta_->at(j),partonPhi_->at(j),GenJetEta_->at(0),GenJetPhi_->at(0));
		   if(deltaR_jt_1<deltaR_jt_1_min) deltaR_jt_1_min=deltaR_jt_1;
		   
		   double deltaR_jt_2=DeltaR(partonEta_->at(j),partonPhi_->at(j),GenJetEta_->at(1),GenJetPhi_->at(1));
		   if(deltaR_jt_2<deltaR_jt_2_min) deltaR_jt_2_min=deltaR_jt_2;
		 }
		 
		 if(deltaR_jt_1_min<0.8) PtLeadingGenJetMatched->Fill(GenJetPt_->at(0),hweight);			     
		 if(deltaR_jt_2_min<0.8) PtSubLeadingGenJetMatched->Fill(GenJetPt_->at(1),hweight);			     


		 /*double deltaR_jt1_1jet=DeltaR(partonEta_->at(0),partonPhi_->at(0),GenJetEta_->at(0),GenJetPhi_->at(0));
		 double deltaR_jt1_2jet=DeltaR(partonEta_->at(0),partonPhi_->at(0),GenJetEta_->at(1),GenJetPhi_->at(1));
		 
		 if(deltaR_jt1_1jet<deltaR_jt1_2jet) {deltaR_jt_1_min=deltaR_jt1_1jet;}
		 else {deltaR_jt_1_min=deltaR_jt1_2jet;}

		 if(deltaR_jt_1_min<0.8) PtLeadingGenJetMatched->Fill(partonPt_->at(0),hweight);			     
		 //consider the second parton and the first two jets
		 double deltaR_jt2_1jet=DeltaR(partonEta_->at(1),partonPhi_->at(1),GenJetEta_->at(0),GenJetPhi_->at(0));
		 double deltaR_jt2_2jet=DeltaR(partonEta_->at(1),partonPhi_->at(1),GenJetEta_->at(1),GenJetPhi_->at(1));
	       
		 if(deltaR_jt2_1jet<deltaR_jt2_2jet) {deltaR_jt_2_min=deltaR_jt2_1jet;}
		 else {deltaR_jt_2_min=deltaR_jt2_2jet;}

		 if(deltaR_jt_2_min<0.8) PtSubLeadingGenJetMatched->Fill(partonPt_->at(1),hweight);			   */  
		 
		 if(deltaR_jt_1_min>0.8) {MassGENJetLead_TopParton->Fill(GenJetMassSoftDrop_->at(0),hweight); MassGENJetLead_SubJet1TopParton->Fill(GenSubJet1Mass_->at(0),hweight);

		   if(partonPt_->at(0)>400 &&partonPt_->at(1)>400){

		     double deltaPhiTopJets=DeltaPhi(partonPhi_->at(0),partonPhi_->at(1));
		     DeltaPhiPARTONUnmatched->Fill(deltaPhiTopJets,hweight);
		   }

		   if(nGenJets_>2){
		     
		     //if(GenJetPt_->at(2)>400 && fabs(GenJetEta_->at(2))<2.4 && GenJetMassSoftDrop_->at(2)>50 && isBJetGen_->at(2)){
		     double deltaR_jt2_1jet=DeltaR(partonEta_->at(0),partonPhi_->at(0),GenJetEta_->at(2),GenJetPhi_->at(2));
		     double deltaR_jt2_2jet=DeltaR(partonEta_->at(1),partonPhi_->at(1),GenJetEta_->at(2),GenJetPhi_->at(2));
		     
		     if(deltaR_jt2_1jet<0.8){
		       ptTopParton_Leading_SubSubLeading->Fill(partonPt_->at(0),hweight);
		     }
		     
		     if(deltaR_jt2_2jet<0.8){
		       ptTopParton_SubLeading_SubSubLeading->Fill(partonPt_->at(1),hweight);
		     }
		   }
		 }

		 if(deltaR_jt_2_min>0.8)  {MassGENJetSubLead_TopParton->Fill(GenJetMassSoftDrop_->at(1),hweight); MassGENJetSubLead_SubJet1TopParton->Fill(GenSubJet1Mass_->at(1),hweight); 		  
		   if(partonPt_->at(0)>400 &&partonPt_->at(1)>400){
		     double deltaPhiTopJets=DeltaPhi(partonPhi_->at(0),partonPhi_->at(1));
		     DeltaPhiPARTONUnmatched->Fill(deltaPhiTopJets,hweight);
		   }
		 }

		 DeltaR_GenJetPartonLead->Fill(deltaR_jt_1_min,hweight);
		 DeltaR_GenJetPartonSubLead->Fill(deltaR_jt_2_min,hweight);
	       }
	     }
	   }
	 	 
	   if(partonPt_->size()>1){

	     //ptTopParton_Leading->Fill(partonPt_->at(0),hweight);
	     //ptTopParton_SubLeading->Fill(partonPt_->at(1),hweight);
	     	     
	     if(nGenJets_>1){ 

	       if(GenJetPt_->at(0)>400 && fabs(GenJetEta_->at(0))<2.4 && GenJetPt_->at(1)>400 && fabs(GenJetEta_->at(1))<2.4){
		 
		 ptTopParton_Leading->Fill(GenJetPt_->at(0),hweight);
		 ptTopParton_SubLeading->Fill(GenJetPt_->at(1),hweight);
		 
		 //consider the first parton and the first two jets
		 double deltaR_jt1_1jet=DeltaR(partonEta_->at(0),partonPhi_->at(0),GenJetEta_->at(0),GenJetPhi_->at(0));
		 double deltaR_jt1_2jet=DeltaR(partonEta_->at(0),partonPhi_->at(0),GenJetEta_->at(1),GenJetPhi_->at(1));
		 
		 //consider the second parton and the first two jets
		 double deltaR_jt2_1jet=DeltaR(partonEta_->at(1),partonPhi_->at(1),GenJetEta_->at(0),GenJetPhi_->at(0));
		 double deltaR_jt2_2jet=DeltaR(partonEta_->at(1),partonPhi_->at(1),GenJetEta_->at(1),GenJetPhi_->at(1));
		 
		 if(deltaR_jt1_1jet<0.8 || deltaR_jt2_1jet<0.8){
		   ptTopParton_Leading_AfterMatching->Fill(GenJetPt_->at(0),hweight);
		 }

		 if(deltaR_jt1_2jet<0.8 || deltaR_jt2_2jet<0.8){
		   ptTopParton_SubLeading_AfterMatching->Fill(GenJetPt_->at(1),hweight);
		 }
  	 
		 if((deltaR_jt1_1jet<0.8 || deltaR_jt2_1jet<0.8) && (deltaR_jt1_2jet<0.8 || deltaR_jt2_2jet<0.8)){ //the tops are matched to jets
		   //ptTopParton_Leading_AfterMatching->Fill(partonPt_->at(0),hweight);
		   //ptTopParton_SubLeading_AfterMatching->Fill(partonPt_->at(1),hweight);
		   //ptTopParton_Leading_AfterMatching->Fill(GenJetPt_->at(0),hweight);
		   //ptTopParton_SubLeading_AfterMatching->Fill(GenJetPt_->at(1),hweight);
		 		 
		   if(GenJetMassSoftDrop_->at(0)>50 && isBJetGen_->at(0) && GenJetMassSoftDrop_->at(1)>50 && isBJetGen_->at(1)){
		     //if(isBJetGen_->at(0) && isBJetGen_->at(1)){
		     //ptTopParton_Leading_AfterMatchingTagging->Fill(partonPt_->at(0),hweight);
		     //ptTopParton_SubLeading_AfterMatchingTagging->Fill(partonPt_->at(1),hweight);
		     ptTopParton_Leading_AfterMatchingTagging->Fill(GenJetPt_->at(0),hweight);
		     ptTopParton_SubLeading_AfterMatchingTagging->Fill(GenJetPt_->at(1),hweight);
		   }
		 		    
		 }
	       }
	     }

	   }

	   for(unsigned j=0; j<partonPt_->size();j++){
	     ptTopParton->Fill(partonPt_->at(j),hweight);
	     yTopParton->Fill(partonEta_->at(j),hweight);
	     phiTopParton->Fill(partonPhi_->at(j),hweight);
	     EnergyTopParton->Fill(partonEnergy_->at(j),hweight);
	     
	     for(int i=0; i< nGenJets_; i++){
	       if(partonPt_->at(j)>100 && GenJetPt_->at(i)>400){
		 double deltaR_jt=DeltaR(partonEta_->at(j),partonPhi_->at(j),GenJetEta_->at(i),GenJetPhi_->at(i));

		 if(deltaR_jt<0.4){//only positive top selected
		   double resolutionPt=(GenJetPt_->at(i)-partonPt_->at(j))/partonPt_->at(j);
		   double resolutionPhi=GenJetPhi_->at(i)-partonPhi_->at(j);

		   ResolutionPtPartonGenJet->Fill(resolutionPt,hweight);
		   ResolutionPhiPartonGenJet->Fill(resolutionPhi,hweight);
		   
		   ResponsePtPartonGenJet->Fill(GenJetPt_->at(i),partonPt_->at(j),hweight);
		   ResponsePhiPartonGenJet->Fill(GenJetPhi_->at(i),partonPhi_->at(j),hweight);
		   
		   // check the fractions of the top which are b-tagged and W-tagged
		   PtMatchedGenJet->Fill(GenJetPt_->at(i),hweight);
		   if(isBJetGen_->at(i)) PtMatchedGenJetAndBJet->Fill(GenJetPt_->at(i),hweight);
		   if(isWJetGen_->at(i) && mMCtype=="TOP") PtMatchedGenJetAndWJet->Fill(GenJetPt_->at(i),hweight);
		   if(isWJetGen_->at(i) && isBJetGen_->at(i) && mMCtype=="TOP") PtMatchedGenJetAndBWJet->Fill(GenJetPt_->at(i),hweight);
		   
		   if(matchWGenJet[i]) PtMatchedGenJetAndWBosonJet->Fill(GenJetPt_->at(i),hweight);
		   if(matchWGenJet[i] && isBJetGen_->at(i)) { PtMatchedGenJetAndBWBosonJet->Fill(GenJetPt_->at(i),hweight);}		 

		   if(GenSubJet1Mass_->at(i) > 30) PtMatchedGenJetAndWSubjetMassJet->Fill(GenJetPt_->at(i),hweight);
		   if(GenSubJet1Mass_->at(i) > 30 && isBJetGen_->at(i)) { PtMatchedGenJetAndBWSubjetMassJet->Fill(GenJetPt_->at(i),hweight);}		 
		   		   
		   if(GenSubJet1Mass_->at(i) > 70) PtMatchedGenJetAndWSubjetMassHighJet->Fill(GenJetPt_->at(i),hweight);
		   if(GenSubJet1Mass_->at(i) > 70) { PtMatchedGenJetAndBWSubjetMassHighJet->Fill(GenJetPt_->at(i),hweight);}		 
		 }
	       }
	     }
	     
	     for(int i=0; i< nJets_; i++){
	       if(partonPt_->at(j)>100 && jetPt_->at(i)>400){
		 double deltaR_jt=DeltaR(partonEta_->at(j),partonPhi_->at(j),jetEta_->at(i),jetPhi_->at(i));
		 DeltaR_DetJetParton->Fill(deltaR_jt,hweight);
		 if(deltaR_jt<0.4){
		   double resolutionPt=(jetPt_->at(i)-partonPt_->at(j))/partonPt_->at(j);
		   double resolutionPhi=jetPhi_->at(i)-partonPhi_->at(j);
		   
		   ResolutionPtPartonDetJet->Fill(resolutionPt,hweight);
		   ResolutionPhiPartonDetJet->Fill(resolutionPhi,hweight);
		   
		   // check the fractions of the top which are b-tagged for the subjets
		   PtMatchedRecoJet->Fill(jetPt_->at(i),hweight);
		   if(jetbtagSub0_->at(i)>0.8484 || jetbtagSub1_->at(i)>0.8484) PtMatchedRecoJetSubBtags->Fill(jetPt_->at(i),hweight);
		   
		 }
	       }
	     }
	   }
	 
	   //Study of the partonic reconstruction
	   if(partonPt_->size()>1 && decay_==0){
	     ptTopPartonLeading->Fill(partonPt_->at(0),hweight);
	     
	     bool leadingmatched=false;
	     for(int i=0; i< nGenJets_; i++){
	       if(GenJetPt_->at(i)>400){
		 double deltaR_jt_lead=DeltaR(partonEta_->at(0),partonPhi_->at(0),GenJetEta_->at(i),GenJetPhi_->at(i));
		 if(deltaR_jt_lead<0.4){
		   leadingmatched=true;
		   ptTopPartonMatchedLeading->Fill(partonPt_->at(0),hweight);
		   if(matchWGenJet[i] && isBJetGen_->at(i)) { ptTopPartonMatchedWBGenJetLeading->Fill(partonPt_->at(0),hweight);}		 
		   if(GenJetMassSoftDrop_->at(i) && isBJetGen_->at(i)) { ptTopPartonMatchedWSubjetMassBGenJetLeading->Fill(partonPt_->at(0),hweight);}		 
		 } 	   
	       }
	     }
	     
	     if(leadingmatched){
	       ptTopPartonSubLeading->Fill(partonPt_->at(1),hweight);
	       for(int i=0; i< nGenJets_; i++){
		 if(GenJetPt_->at(i)>400){
		   double deltaR_jt_sublead=DeltaR(partonEta_->at(1),partonPhi_->at(1),GenJetEta_->at(i),GenJetPhi_->at(i));
		   if(deltaR_jt_sublead<0.4){
		     ptTopPartonMatchedSubLeading->Fill(partonPt_->at(1),hweight);
		     if(matchWGenJet[i] && isBJetGen_->at(i)) { ptTopPartonMatchedWBGenJetSubLeading->Fill(partonPt_->at(1),hweight);}		 
		     if(GenJetMassSoftDrop_->at(i) && isBJetGen_->at(i)) { ptTopPartonMatchedWSubjetMassBGenJetSubLeading->Fill(partonPt_->at(1),hweight);}		 
		   } 	   
		   
		 }
	       }
	     }
	   }
	 
	 //Study of the gen-jet reconstruction
	    
	 }
       }

       //Preselection at RECO level
     
       //SVM attempt
       //int SVMoutput=0;
       //file_input>>SVMoutput;
       
       //DNN attempt
       //double DNNoutput=0;
       //file_input>>DNNoutput;

       bool preselection=false;
       if(prepreselection){
	 //if(SVMoutput==1){//SVM selection
	 //if(DNNoutput>0.9){//SVM selection
	   //if(mvaD>-0.35){    //LD cut 4 variables!!
	 //if(mvaD>0.65){      //MLP cut 8 variables!!
	   //if(mvaD>0.45){      //BDT cut Daniela variables!!
	 if(mvaD>0.45){      //MLP 6 variables!!
	   preselection=true;
	 }
       }
	

       // for QCD control region
       /*if(jetPt_->at(0)>500 && fabs(jetEta_->at(0))<2.4 && jetMassSoftDrop_->at(0)>50){//mass drop+b-tag
	 if(jetbtagSub0_->at(0)<0.8484 && jetbtagSub1_->at(0)<0.8484){ //b-tagged subjets
	   if(jetPt_->at(1)>500 && fabs(jetEta_->at(1))<2.4 && jetMassSoftDrop_->at(1)>50){//mass drop+b-tag
	     if(jetbtagSub0_->at(1)<0.8484 && jetbtagSub1_->at(1)<0.8484){ //b-tagged subjets
	       btagjets=true;
	     }
	   }
	 }
	 }*/


       if(preselection){

	 double tau31_1jet=-1;
         if(jettau1_->at(0)!=0){ tau31_1jet= jettau3_->at(0)/jettau1_->at(0);}
         double tau32_1jet=-1;
         if(jettau2_->at(0)!=0){ tau32_1jet= jettau3_->at(0)/jettau2_->at(0);}
	 
	 double tau31_2jet=-1;
	 if(jettau1_->at(1)!=0){ tau31_2jet= jettau3_->at(1)/jettau1_->at(1);}
         double tau32_2jet=-1;
         if(jettau2_->at(1)!=0){ tau32_2jet= jettau3_->at(1)/jettau2_->at(1);}
	 
	 double DeltaRSubJet_1jet=DeltaR(jetetaSub0_->at(0),jetphiSub0_->at(0),jetetaSub1_->at(0),jetphiSub1_->at(0));
	 double DeltaRSubJet_2jet=DeltaR(jetetaSub0_->at(1),jetphiSub0_->at(1),jetetaSub1_->at(1),jetphiSub1_->at(1));
      
	 //cout<<"QCD:-1 TOP:1 "<<jetMass_->at(0)<<" "<<jetmassSub0_->at(0)<<" "<<jetmassSub1_->at(0)<<" "<<tau31_1jet<<" "<<tau32_1jet<<" "<<jetptSub0_->at(0)<<" "<<jetptSub1_->at(0)<<" "<<DeltaRSubJet_1jet<<" "<<jetMass_->at(1)<<" "<<jetmassSub0_->at(1)<<" "<<jetmassSub1_->at(1)<<" "<<tau31_2jet<<" "<<tau32_2jet<<" "<<jetptSub0_->at(1)<<" "<<jetptSub1_->at(1)<<" "<<DeltaRSubJet_2jet<<endl;
	 
	 //double tau31Leading = tau31_1jet;
	 //double tau32Leading = tau32_1jet;
	 //double MassLeading = jetMass_->at(0);
	 //double SDmassLeading = jetMassSoftDrop_->at(0);
	 //double MassSubJet0Leading = jetmassSub0_->at(0);
	 //double MassSubJet1Leading = jetmassSub0_->at(1);
	 //double jetSubJetpt0Leading = jetptSub0_->at(0);
	 //double jetSubJetpt1Leading = jetptSub0_->at(1);
	 
	 //double tau31SubLeading = tau31_2jet;
	 //double tau32SubLeading = tau32_2jet;
	 //double MassSubLeading = jetMass_->at(1);
	 //double SDmassSubLeading = jetMassSoftDrop_->at(0);
	 //double MassSubJet0SubLeading = jetmassSub1_->at(0);
	 //double MassSubJet1SubLeading = jetmassSub1_->at(1);
	 //double jetSubJetpt0SubLeading = jetptSub1_->at(0);
	 //double jetSubJetpt1SubLeading = jetptSub1_->at(1);
	 
	 //double mvaD = discr_->eval(tau31Leading,tau32Leading,MassLeading,MassSubJet0Leading,tau31SubLeading,tau32SubLeading,MassSubLeading,MassSubJet0SubLeading);
	 //double mvaD = discr_->eval(tau31Leading,tau32Leading,MassLeading,MassSubJet0Leading);

	 double scale_factor=1.;

	 bool applyscale_factors=false;

	 if(mIsMCarlo && applyscale_factors){
	   double jet_scalefactor_1subjet=1.;
	   if(jetbtagSub0_->at(0)>0.8484){//b- and c-quark scale factors are the same for subjets!
	     if(fabs(jetflavorSub0_->at(0))==5) jet_scalefactor_1subjet = readerHF.eval_auto_bounds("central",BTagEntry::FLAV_B, jetetaSub0_->at(0), jetptSub0_->at(0)); 
	     if(fabs(jetflavorSub0_->at(0))==4) jet_scalefactor_1subjet = readerHF.eval_auto_bounds("central",BTagEntry::FLAV_B, jetetaSub0_->at(0), jetptSub0_->at(0)); 
	     if(fabs(jetflavorSub0_->at(0))!=4 && fabs(jetflavorSub0_->at(0))!=5) jet_scalefactor_1subjet = reader.eval_auto_bounds("central",BTagEntry::FLAV_UDSG, jetetaSub0_->at(0), jetptSub0_->at(0)); 
	     
	   }
	   
	   if(jetbtagSub0_->at(1)>0.8484){
	     if(fabs(jetflavorSub0_->at(1))==5) jet_scalefactor_1subjet = readerHF.eval_auto_bounds("central",BTagEntry::FLAV_B, jetetaSub0_->at(1), jetptSub0_->at(1)); 
	     if(fabs(jetflavorSub0_->at(1))==4) jet_scalefactor_1subjet = readerHF.eval_auto_bounds("central",BTagEntry::FLAV_B, jetetaSub0_->at(1), jetptSub0_->at(1)); 
	     if(fabs(jetflavorSub0_->at(1))!=4 && fabs(jetflavorSub0_->at(1))!=5) jet_scalefactor_1subjet = reader.eval_auto_bounds("central",BTagEntry::FLAV_UDSG, jetetaSub0_->at(1), jetptSub0_->at(1)); 
	   }


	   double jet_scalefactor_2subjet=1.;
	   if(jetbtagSub1_->at(0)>0.8484){
	     if(fabs(jetflavorSub1_->at(0))==5) jet_scalefactor_2subjet = readerHF.eval_auto_bounds("central",BTagEntry::FLAV_B, jetetaSub1_->at(0), jetptSub1_->at(0)); 
	     if(fabs(jetflavorSub1_->at(0))==4) jet_scalefactor_2subjet = readerHF.eval_auto_bounds("central",BTagEntry::FLAV_B, jetetaSub1_->at(0), jetptSub1_->at(0)); 
	     if(fabs(jetflavorSub1_->at(0))!=4 && fabs(jetflavorSub1_->at(0))!=5) jet_scalefactor_2subjet = reader.eval_auto_bounds("central",BTagEntry::FLAV_UDSG, jetetaSub1_->at(0), jetptSub1_->at(0)); 
	   }
	   
	   if(jetbtagSub1_->at(1)>0.8484){
	     if(fabs(jetflavorSub1_->at(1))==5) jet_scalefactor_2subjet = readerHF.eval_auto_bounds("central",BTagEntry::FLAV_B, jetetaSub1_->at(1), jetptSub1_->at(1)); 
	     if(fabs(jetflavorSub1_->at(1))==4) jet_scalefactor_2subjet = readerHF.eval_auto_bounds("central",BTagEntry::FLAV_B, jetetaSub1_->at(1), jetptSub1_->at(1)); 
	     if(fabs(jetflavorSub1_->at(1))!=4 && fabs(jetflavorSub1_->at(1))!=5) jet_scalefactor_2subjet = reader.eval_auto_bounds("central",BTagEntry::FLAV_UDSG, jetetaSub1_->at(1), jetptSub1_->at(1)); 
	   }

	   if(jet_scalefactor_1subjet!=0 && jet_scalefactor_2subjet!=0){
	     scale_factor=jet_scalefactor_1subjet*jet_scalefactor_2subjet;
	   }
	   //if(scale_factor<0.01){scale_factor=1.;}
	 }

	 hweight=hweight*scale_factor;
	 hweightAfterFitToSelection=hweightAfterFitToSelection*scale_factor;
	 if(hweightAfterFitToSelection<0.) hweightAfterFitToSelection=1.; 
       
	 //Add histograms
	 PtDETLeadingPreselectionSubJet0->Fill(jetptSub0_->at(0),hweight);
	 PtDETLeadingPreselectionSubJet1->Fill(jetptSub1_->at(0),hweight);
	 EtaDETLeadingPreselectionSubJet0->Fill(jetetaSub0_->at(0),hweight);
	 EtaDETLeadingPreselectionSubJet1->Fill(jetetaSub1_->at(0),hweight);
	 MassDETLeadingPreselectionSubJet0->Fill(jetmassSub0_->at(0),hweight);
	 MassDETLeadingPreselectionSubJet1->Fill(jetmassSub1_->at(0),hweight);
	 PhiDETLeadingPreselectionSubJet0->Fill(jetphiSub0_->at(0),hweight);
	 PhiDETLeadingPreselectionSubJet1->Fill(jetphiSub1_->at(0),hweight);
	 DeltaRDETLeadingPreselectionSubJets->Fill(DeltaRSubJet_1jet,hweight);
	   
	 Tau1DETLeadingJetPreselection->Fill(jettau1_->at(0),hweight);
	 Tau2DETLeadingJetPreselection->Fill(jettau2_->at(0),hweight);
	 Tau3DETLeadingJetPreselection->Fill(jettau3_->at(0),hweight);
	 Tau31DETLeadingJetPreselection->Fill(tau31_1jet,hweight);
	 Tau32DETLeadingJetPreselection->Fill(tau32_1jet,hweight);
	 
	 PtDETSubLeadingPreselectionSubJet0->Fill(jetptSub0_->at(1),hweight);
	 PtDETSubLeadingPreselectionSubJet1->Fill(jetptSub1_->at(1),hweight);
	 EtaDETSubLeadingPreselectionSubJet0->Fill(jetetaSub0_->at(1),hweight);
	 EtaDETSubLeadingPreselectionSubJet1->Fill(jetetaSub1_->at(1),hweight);
	 MassDETSubLeadingPreselectionSubJet0->Fill(jetmassSub0_->at(1),hweight);
	 MassDETSubLeadingPreselectionSubJet1->Fill(jetmassSub1_->at(1),hweight);
	 PhiDETSubLeadingPreselectionSubJet0->Fill(jetphiSub0_->at(1),hweight);
	 PhiDETSubLeadingPreselectionSubJet1->Fill(jetphiSub1_->at(1),hweight);
	 DeltaRDETSubLeadingPreselectionSubJets->Fill(DeltaRSubJet_2jet,hweight);
	 
	 Tau1DETSubLeadingJetPreselection->Fill(jettau1_->at(1),hweight);
	 Tau2DETSubLeadingJetPreselection->Fill(jettau2_->at(1),hweight);
	 Tau3DETSubLeadingJetPreselection->Fill(jettau3_->at(1),hweight);
	 Tau31DETSubLeadingJetPreselection->Fill(tau31_2jet,hweight);
	 Tau32DETSubLeadingJetPreselection->Fill(tau32_2jet,hweight);
	 
	 PtDETLeadingPreselection->Fill(SmearingFactor[0]*jetPt_->at(0),hweight);
	 PtDETSubLeadingPreselection->Fill(SmearingFactor[1]*jetPt_->at(1),hweight);
	 EtaDETLeadingPreselection->Fill(jetEta_->at(0),hweight);
	 EtaDETSubLeadingPreselection->Fill(jetEta_->at(1),hweight);
	 MassDETLeadingPreselection->Fill(jetMass_->at(0),hweight);
	 MassDETSubLeadingPreselection->Fill(jetMass_->at(1),hweight);
	 PhiDETLeadingPreselection->Fill(jetPhi_->at(0),hweight);
	 PhiDETSubLeadingPreselection->Fill(jetPhi_->at(1),hweight);
	 SDMassDETLeadingPreselection->Fill(jetMassSoftDrop_->at(0),hweight);
	 SDMassDETSubLeadingPreselection->Fill(jetMassSoftDrop_->at(1),hweight);
	 
	 //Actual selection after MVA analysis
	 
	 //if(mvaD>0.01){ //TMVA selection
	 
	 //if((jetbtagSub0_->at(0)>0.8484 || jetbtagSub1_->at(0)>0.8484) && (jetbtagSub0_->at(1)>0.8484 || jetbtagSub1_->at(1)>0.8484)){ //b-tagged subjets
	 if(jetbtagSub0_->at(0)< 0.5426 && jetbtagSub1_->at(0)<  0.5426  && jetbtagSub0_->at(1)< 0.5426  && jetbtagSub1_->at(1)< 0.5426 ){ // for QCD control region  
	 //if(jetbtagSub0_->at(0)< 0.8484 && jetbtagSub1_->at(0)<  0.8484  && jetbtagSub0_->at(1)< 0.8484  && jetbtagSub1_->at(1)< 0.8484 ){ // for QCD control region  
	   selection=true;

	   num_of_VtxSelection->Fill(nvtx_,hweight);

	   if(mIsMCarlo){
	     FlavourDETSelectionLeadingJet->Fill(jethadronflavor_->at(0),hweight);
	     FlavourDETSelectionSubLeadingJet->Fill(jethadronflavor_->at(1),hweight);
	   

	     if(jetflavorSub1_->at(0)==5){
	       JetBTaggedSubjets0Mass_SecondSubjet->Fill(jetmassSub0_->at(0),hweight);
	       JetBTaggedSubjets1Mass_SecondSubjet->Fill(jetmassSub1_->at(0),hweight);
	     }

	     if(jetflavorSub0_->at(0)==5){
	       JetBTaggedSubjets0Mass_FirstSubjet->Fill(jetmassSub0_->at(0),hweight);
	       JetBTaggedSubjets1Mass_FirstSubjet->Fill(jetmassSub1_->at(0),hweight);
	     }
	   }

	   PtDETLeadingSelectionSubJet0->Fill(jetptSub0_->at(0),hweight);
	   PtDETLeadingSelectionSubJet1->Fill(jetptSub1_->at(0),hweight);
	   EtaDETLeadingSelectionSubJet0->Fill(jetetaSub0_->at(0),hweight);
	   EtaDETLeadingSelectionSubJet1->Fill(jetetaSub1_->at(0),hweight);
	   MassDETLeadingSelectionSubJet0->Fill(jetmassSub0_->at(0),hweight);
	   MassDETLeadingSelectionSubJet1->Fill(jetmassSub1_->at(0),hweight);
	   PhiDETLeadingSelectionSubJet0->Fill(jetphiSub0_->at(0),hweight);
	   PhiDETLeadingSelectionSubJet1->Fill(jetphiSub1_->at(0),hweight);
	   DeltaRDETLeadingSelectionSubJets->Fill(DeltaRSubJet_1jet,hweight);
	   
	   Tau1DETLeadingJetSelection->Fill(jettau1_->at(0),hweight);
	   Tau2DETLeadingJetSelection->Fill(jettau2_->at(0),hweight);
	   Tau3DETLeadingJetSelection->Fill(jettau3_->at(0),hweight);
	   Tau31DETLeadingJetSelection->Fill(tau31_1jet,hweight);
	   Tau32DETLeadingJetSelection->Fill(tau32_1jet,hweight);
	   
	   PtDETSubLeadingSelectionSubJet0->Fill(jetptSub0_->at(1),hweight);
	   PtDETSubLeadingSelectionSubJet1->Fill(jetptSub1_->at(1),hweight);
	   EtaDETSubLeadingSelectionSubJet0->Fill(jetetaSub0_->at(1),hweight);
	   EtaDETSubLeadingSelectionSubJet1->Fill(jetetaSub1_->at(1),hweight);
	   MassDETSubLeadingSelectionSubJet0->Fill(jetmassSub0_->at(1),hweight);
	   MassDETSubLeadingSelectionSubJet1->Fill(jetmassSub1_->at(1),hweight);
	   PhiDETSubLeadingSelectionSubJet0->Fill(jetphiSub0_->at(1),hweight);
	   PhiDETSubLeadingSelectionSubJet1->Fill(jetphiSub1_->at(1),hweight);
	   DeltaRDETSubLeadingSelectionSubJets->Fill(DeltaRSubJet_2jet,hweight);
	   
	   Tau1DETSubLeadingJetSelection->Fill(jettau1_->at(1),hweight);
	   Tau2DETSubLeadingJetSelection->Fill(jettau2_->at(1),hweight);
	   Tau3DETSubLeadingJetSelection->Fill(jettau3_->at(1),hweight);
	   Tau31DETSubLeadingJetSelection->Fill(tau31_2jet,hweight);
	   Tau32DETSubLeadingJetSelection->Fill(tau32_2jet,hweight);
	   
	   if(mShift!=0.) PtDETLeadingSelection->Fill(SmearingFactor[0]*jetPt_->at(0)*(1+mShift*jetUnc_->at(0)),hweight);
	   if(mShift!=0.) PtDETSubLeadingSelection->Fill(SmearingFactor[1]*jetPt_->at(1)*(1+mShift*jetUnc_->at(1)),hweight);
	   if(mShift==0.) PtDETLeadingSelection->Fill(SmearingFactor[0]*jetPt_->at(0),hweight);
	   if(mShift==0.) PtDETSubLeadingSelection->Fill(SmearingFactor[1]*jetPt_->at(1),hweight);
	   EtaDETLeadingSelection->Fill(jetEta_->at(0),hweight);
	   EtaDETSubLeadingSelection->Fill(jetEta_->at(1),hweight);
	   MassDETLeadingSelection->Fill(jetMass_->at(0),hweight);
	   MassDETSubLeadingSelection->Fill(jetMass_->at(1),hweight);
	   PhiDETLeadingSelection->Fill(jetPhi_->at(0),hweight);
	   PhiDETSubLeadingSelection->Fill(jetPhi_->at(1),hweight);
	   SDMassDETLeadingSelection->Fill(jetMassSoftDrop_->at(0),hweight);
	   SDMassDETSubLeadingSelection->Fill(jetMassSoftDrop_->at(1),hweight);

	   double deltaPhiTopJets=DeltaPhi(jetPhi_->at(0),jetPhi_->at(1));
	   DeltaPhiDETSelection->Fill(deltaPhiTopJets,hweight);
	   DeltaPhiDETSelectionDoubleBins->Fill(deltaPhiTopJets,hweight);
	   DeltaPhiDETSelectionForToyFits->Fill(deltaPhiTopJets,hweight);

	   Neventsdetselection++;

	   double ptedges[7]={400,500,600,800,1000,1500,2000};

	   for(int iy = 0; iy < 6; iy++) {
	     if(SmearingFactor[0]*jetPt_->at(0) > ptedges[iy] && SmearingFactor[0]*jetPt_->at(0) <= ptedges[iy+1]) {
	       MassDETLeadingSelection_bins[iy]->Fill(jetMass_->at(0),hweight);
	     }
	   }

	   double yedges[6]={0,0.5,1.0,1.5,2.0,2.5};

	   for(int iy = 0; iy < 5; iy++) {
	     if(fabs(jetEta_->at(0)) > yedges[iy] && fabs(jetEta_->at(0)) <= yedges[iy+1]) {
	       if(mShift!=0.) PtDETLeadingSelection_bins[iy]->Fill(SmearingFactor[0]*jetPt_->at(0)*(1+mShift*jetUnc_->at(0)),hweight);
	       if(mShift==0.) PtDETLeadingSelection_bins[iy]->Fill(SmearingFactor[0]*jetPt_->at(0),hweight);
	       PtDETLeadingSelectionForDoubleBins_bins[iy]->Fill(SmearingFactor[0]*jetPt_->at(0),hweight);
	       PtDETLeadingSelectionForToyFits_bins[iy]->Fill(SmearingFactor[0]*jetPt_->at(0),hweight);
	     }
	     
	     if(fabs(jetEta_->at(1)) > yedges[iy] && fabs(jetEta_->at(1)) <= yedges[iy+1]) {
	       if(mShift!=0.) PtDETSubLeadingSelection_bins[iy]->Fill(SmearingFactor[1]*jetPt_->at(1)*(1+mShift*jetUnc_->at(1)),hweight);
	       if(mShift==0.) PtDETSubLeadingSelection_bins[iy]->Fill(SmearingFactor[1]*jetPt_->at(1),hweight);
	       PtDETSubLeadingSelectionForToyFits_bins[iy]->Fill(SmearingFactor[1]*jetPt_->at(1),hweight);
	       PtDETSubLeadingSelectionForDoubleBins_bins[iy]->Fill(SmearingFactor[1]*jetPt_->at(1),hweight);
	     }
	   }
	   
	   double DeltaPhiSlices[]={0.0,1.0,1.75,2.25,2.625,2.875,3.14};
	   double PtSlices[]={395, 471.5, 548, 642.5, 737, 852., 4252.};
     
	   //filling the SD masses for the different DeltaPhi bins
	   for(int iy = 0; iy < 6; iy++) {
	     if(deltaPhiTopJets<DeltaPhiSlices[iy+1] && deltaPhiTopJets>DeltaPhiSlices[iy]){
	       SDMassDETLeadingSelection_DeltaPhibins[iy]->Fill(jetMassSoftDrop_->at(0),hweight);
	       SDMassDETSubLeadingSelection_DeltaPhibins[iy]->Fill(jetMassSoftDrop_->at(1),hweight);
	     }
	   }
	   
	   for(int iy = 0; iy < 6; iy++) {
	     if((SmearingFactor[0]*jetPt_->at(0))<PtSlices[iy+1] && (SmearingFactor[0]*jetPt_->at(0))>PtSlices[iy]){
	       SDMassDETLeadingSelection_Ptbins[iy]->Fill(jetMassSoftDrop_->at(0),hweight);
	     }
	     if((SmearingFactor[1]*jetPt_->at(1))<PtSlices[iy+1] && (SmearingFactor[1]*jetPt_->at(1))>PtSlices[iy]){
	       SDMassDETSubLeadingSelection_Ptbins[iy]->Fill(jetMassSoftDrop_->at(1),hweight);
	     }
	   }
	 }
       }
	
       if(mIsMCarlo){
	 //Match with RecoJets
	 
	 if(preselection && genpreselection){
	     
	   //Here you don't require the matching in eta and phi of the gen-reco jets
	   
	   double resolutionPt=(SmearingFactor[0]*jetPt_->at(0)-GenJetPt_->at(0))/GenJetPt_->at(0);
	   double resolutionPhi=jetPhi_->at(0)-GenJetPhi_->at(0);	       
	   ResolutionPtRecoGenJet->Fill(resolutionPt,hweight);
	   ResolutionPhiRecoGenJet->Fill(resolutionPhi,hweight);
	   ResponsePtRecoGenJet->Fill(GenJetPt_->at(0),SmearingFactor[0]*jetPt_->at(0),hweight);
	   ResponsePhiRecoGenJet->Fill(GenJetPhi_->at(0),jetPhi_->at(0),hweight);
	   double resolutionMass=(jetMass_->at(0)-GenJetmass_->at(0))/GenJetmass_->at(0);
	   ResolutionMassRecoGenJet->Fill(resolutionMass,hweight);
	   double resolutionSDMass=(jetMassSoftDrop_->at(0)-GenJetMassSoftDrop_->at(0))/GenJetMassSoftDrop_->at(0);
	   ResolutionSDMassRecoGenJet->Fill(resolutionSDMass,hweight);
	   
	   resolutionPt=(SmearingFactor[1]*jetPt_->at(1)-GenJetPt_->at(1))/GenJetPt_->at(1);
	   resolutionPhi=jetPhi_->at(1)-GenJetPhi_->at(1);	       
	   ResolutionPtRecoGenJet->Fill(resolutionPt,hweight);
	   ResolutionPhiRecoGenJet->Fill(resolutionPhi,hweight);
	   ResponsePtRecoGenJet->Fill(GenJetPt_->at(1),SmearingFactor[1]*jetPt_->at(1),hweight);
	   ResponsePhiRecoGenJet->Fill(GenJetPhi_->at(1),jetPhi_->at(1),hweight);
	   resolutionMass=(jetMass_->at(1)-GenJetmass_->at(1))/GenJetmass_->at(1);
	   ResolutionMassRecoGenJet->Fill(resolutionMass,hweight);
	   resolutionSDMass=(jetMassSoftDrop_->at(1)-GenJetMassSoftDrop_->at(1))/GenJetMassSoftDrop_->at(1);
	   ResolutionSDMassRecoGenJet->Fill(resolutionSDMass,hweight);
	   
	   double resolutionRapidity=jetEta_->at(0)-GenJetEta_->at(0);	       
	   ResolutionRapidityRecoGenJet->Fill(resolutionRapidity,hweight);
	   ResponseRapidityRecoGenJet->Fill(GenJetEta_->at(0),jetEta_->at(0),hweight);
	   
	   resolutionRapidity=jetEta_->at(1)-GenJetEta_->at(1);	       
	   ResolutionRapidityRecoGenJet->Fill(resolutionRapidity,hweight);
	   ResponseRapidityRecoGenJet->Fill(GenJetEta_->at(1),jetEta_->at(1),hweight);
	   
	 }          
	 
	 
	 if(genpreselection && !selection) {
	   double deltaPhiGen=DeltaPhi(GenJetPhi_->at(1),GenJetPhi_->at(0));
	   double deltaPhiReco=-0.5;
	   //ResponseForTUnfoldDeltaPhiRecoGenJet->Fill(deltaPhiReco,deltaPhiGen,hweight);
	   ResponseForTUnfoldDeltaPhiRecoGenJet->Fill(deltaPhiReco,deltaPhiGen,hweightAfterFitToSelection);

	   DeltaPhiDETSelectionMiss->Fill(deltaPhiGen,hweight);
	 
	   double yedges[6]={0,0.5,1.0,1.5,2.0,2.5};
	   
	   double RecoPt=300;
	   for(int iy = 0; iy < 5; iy++) {
	     if(fabs(GenJetEta_->at(0)) > yedges[iy] && fabs(GenJetEta_->at(0)) <= yedges[iy+1]) {
	       //ResponseLeadingJetPtForTUnfold_bins[iy]->Fill(RecoPt,GenJetPt_->at(0),hweight);
	       ResponseLeadingJetPtForTUnfold_bins[iy]->Fill(RecoPt,GenJetPt_->at(0),hweightAfterFitToSelection);
	     }
	     
	     if(fabs(GenJetEta_->at(1)) > yedges[iy] && fabs(GenJetEta_->at(1)) <= yedges[iy+1]) {
	       //ResponseSubLeadingJetPtForTUnfold_bins[iy]->Fill(RecoPt,GenJetPt_->at(1),hweight);
	       ResponseSubLeadingJetPtForTUnfold_bins[iy]->Fill(RecoPt,GenJetPt_->at(1),hweightAfterFitToSelection);
	     }
	     
	   }

	 }
	 
	 if(!genpreselection && selection) {
	   double deltaPhiGen=-0.5;
	   double deltaPhiReco=DeltaPhi(jetPhi_->at(1),jetPhi_->at(0));
	   ResponseForTUnfoldDeltaPhiRecoGenJet->Fill(deltaPhiReco,deltaPhiGen,hweightAfterFitToSelection);
	   
	   DeltaPhiDETSelectionFake->Fill(deltaPhiReco,hweightAfterFitToSelection);
	 
	   double yedges[6]={0,0.5,1.0,1.5,2.0,2.5};
	   
	   double GenPt=300;

	   for(int iy = 0; iy < 5; iy++) {
	     if(fabs(jetEta_->at(0)) > yedges[iy] && fabs(jetEta_->at(0)) <= yedges[iy+1]) {
	       ResponseLeadingJetPtForTUnfold_bins[iy]->Fill(SmearingFactor[0]*jetPt_->at(0),GenPt,hweightAfterFitToSelection);
	     }
	     
	     if(fabs(jetEta_->at(1)) > yedges[iy] && fabs(jetEta_->at(1)) <= yedges[iy+1]) {
	       ResponseSubLeadingJetPtForTUnfold_bins[iy]->Fill(SmearingFactor[1]*jetPt_->at(1),GenPt,hweightAfterFitToSelection);
	     }
	     
	   }

	 }
	 
	 if(genpreselection && selection) { 
	   double deltaPhiTopJets=DeltaPhi(GenJetPhi_->at(0),GenJetPhi_->at(1));
	   DeltaPhiMatchedSelection->Fill(deltaPhiTopJets,hweightAfterFitToSelection);
	   
	   double deltaPhiGen=DeltaPhi(GenJetPhi_->at(1),GenJetPhi_->at(0));
	   double deltaPhiReco=DeltaPhi(jetPhi_->at(1),jetPhi_->at(0));
	   
	   double resolutionDeltaPhi=deltaPhiGen-deltaPhiReco;	       
	   ResolutionDeltaPhiRecoGenJet->Fill(resolutionDeltaPhi,hweightAfterFitToSelection);
	   ResponseDeltaPhiRecoGenJet->Fill(deltaPhiReco,deltaPhiGen,hweightAfterFitToSelection);
	   ResponseForTUnfoldDeltaPhiRecoGenJet->Fill(deltaPhiReco,deltaPhiGen,hweightAfterFitToSelection);

	   double hweightAfterFitToSelectionForGenMisses=1-hweightAfterFitToSelection;
	   if(hweightAfterFitToSelectionForGenMisses<0) hweightAfterFitToSelectionForGenMisses=0.;
	   //Compensating the SF applied to RECO for GEN level
	   ResponseForTUnfoldDeltaPhiRecoGenJet->Fill(-0.5,deltaPhiGen,hweightAfterFitToSelectionForGenMisses);
	   
	   double yedges[6]={0,0.5,1.0,1.5,2.0,2.5};
	   
	   for(int iy = 0; iy < 5; iy++) {
	     if(fabs(jetEta_->at(0)) > yedges[iy] && fabs(jetEta_->at(0)) <= yedges[iy+1]) {
	       //PtDETLeadingSelection_bins[iy]->Fill(SmearingFactor[0]*jetPt_->at(0)*(1+mShift*jetUnc_->at(0)),hweight);
	       ResponseLeadingJetPt_bins[iy]->Fill(SmearingFactor[0]*jetPt_->at(0),GenJetPt_->at(0),hweightAfterFitToSelection);
	       ResponseLeadingJetPtForTUnfold_bins[iy]->Fill(SmearingFactor[0]*jetPt_->at(0),GenJetPt_->at(0),hweightAfterFitToSelection);

	       //Compensating the SF applied to RECO for GEN level
	       ResponseLeadingJetPtForTUnfold_bins[iy]->Fill(300,GenJetPt_->at(0),hweightAfterFitToSelectionForGenMisses);
	       
	     }
	     
	     if(fabs(jetEta_->at(1)) > yedges[iy] && fabs(jetEta_->at(1)) <= yedges[iy+1]) {
	       //PtDETSubLeadingSelection_bins[iy]->Fill(SmearingFactor[0]*jetPt_->at(1)*(1+mShift*jetUnc_->at(1)),hweight);
	       ResponseSubLeadingJetPt_bins[iy]->Fill(SmearingFactor[1]*jetPt_->at(1),GenJetPt_->at(1),hweightAfterFitToSelection);
	       ResponseSubLeadingJetPtForTUnfold_bins[iy]->Fill(SmearingFactor[1]*jetPt_->at(1),GenJetPt_->at(1),hweightAfterFitToSelection);

	       //Compensating the SF applied to RECO for GEN level
	       ResponseSubLeadingJetPtForTUnfold_bins[iy]->Fill(300,GenJetPt_->at(1),hweightAfterFitToSelectionForGenMisses);
	     }
	     
	   }
	 }
       
	 if(genpreselection && decay_==0){ //hadronic decays
	   //if(genpreselection && decay_==2){ //leptonic decays
	   
	   Neventsgenselection++;
	   
	   double yedges[6]={0,0.5,1.0,1.5,2.0,2.5};
	   
	   for(int iy = 0; iy < 5; iy++) {
	     if(fabs(GenJetEta_->at(0)) > yedges[iy] && fabs(GenJetEta_->at(0)) <= yedges[iy+1]) {
	       PtGENLeadingSelection_bins[iy]->Fill(GenJetPt_->at(0),hweight);
	       PtGENLeadingSelectionForToyFits_bins[iy]->Fill(GenJetPt_->at(0),hweight);
	       PtGENLeadingSelectionForDoubleBins_bins[iy]->Fill(GenJetPt_->at(0),hweight);
	     }
	     
	     if(fabs(GenJetEta_->at(1)) > yedges[iy] && fabs(GenJetEta_->at(1)) <= yedges[iy+1]) {
	       PtGENSubLeadingSelection_bins[iy]->Fill(GenJetPt_->at(1),hweight);
	       PtGENSubLeadingSelectionForToyFits_bins[iy]->Fill(GenJetPt_->at(1),hweight);
	       PtGENSubLeadingSelectionForDoubleBins_bins[iy]->Fill(GenJetPt_->at(1),hweight);
	     }
	   }
	     
	   MassGENSubJet0LeadingAfterSelection->Fill(GenSubJet1Mass_->at(0),hweight);
	   MassGENSubJet0SubLeadingAfterSelection->Fill(GenSubJet1Mass_->at(1),hweight);
	   
	   MassGENJetLeadingAfterSelection->Fill(GenJetMassSoftDrop_->at(0),hweight);
	   MassGENJetSubLeadingAfterSelection->Fill(GenJetMassSoftDrop_->at(1),hweight);
	   
	   SubJet1MassGENJetLeadingAfterSelection->Fill(GenSubJet1Mass_->at(0),hweight);
	   SubJet2MassGENJetLeadingAfterSelection->Fill(GenSubJet2Mass_->at(0),hweight);
	   SubJet1MassGENJetSubLeadingAfterSelection->Fill(GenSubJet1Mass_->at(1),hweight);
	   SubJet2MassGENJetSubLeadingAfterSelection->Fill(GenSubJet2Mass_->at(1),hweight);

	   double deltaPhiTopJets=DeltaPhi(GenJetPhi_->at(0),GenJetPhi_->at(1));
	   DeltaPhiGENSelection->Fill(deltaPhiTopJets,hweight);
	   DeltaPhiGENSelectionDoubleBins->Fill(deltaPhiTopJets,hweight);
	   DeltaPhiGENSelectionForToyFits->Fill(deltaPhiTopJets,hweight);
	 }
       
     
	 if(partonpreselection){
	   
	   Neventspartonselection++;
	   
	   double yedges[6]={0,0.5,1.0,1.5,2.0,2.5};
	   
	   for(int iy = 0; iy < 5; iy++) {
	     if(fabs(partonEta_->at(0)) > yedges[iy] && fabs(partonEta_->at(0)) <= yedges[iy+1]) {
	       PtPARTONLeadingSelection_bins[iy]->Fill(partonPt_->at(0),hweight);
	     }
	     
	     if(fabs(partonEta_->at(1)) > yedges[iy] && fabs(partonEta_->at(1)) <= yedges[iy+1]) {
	       PtPARTONSubLeadingSelection_bins[iy]->Fill(partonPt_->at(1),hweight);
	     }
	   }
	   
	   if(partonPt_->at(0)> 400 && partonPt_->at(1)>400){
	   
	     double deltaPhiTopJets=DeltaPhi(partonPhi_->at(0),partonPhi_->at(1));
	     DeltaPhiPARTONSelection->Fill(deltaPhiTopJets,hweight);
	   }
	 }
       }
     } // end of event loop  
   }//end of file loop

   cout<<" Number of selected events " << NEvents<<endl;
   cout<<" Number of selected det events " <<Neventsdetselection<<endl;
   cout<<" Number of selected gen events " <<Neventsgenselection<<endl;
   cout<<" Number of selected parton events " <<Neventspartonselection<<endl;
   cout<<" Number of selected dijets (after and before) " <<dijetsafter<<" "<<dijetsbefore<<endl;
 
 } // closing analyze() function

Analysis_Template_MC::~Analysis_Template_MC()
{
}


DEFINE_FWK_MODULE(Analysis_Template_MC);

