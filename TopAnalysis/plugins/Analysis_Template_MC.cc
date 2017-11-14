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



using namespace std;

//---------------------------- Constructor Of The Class TriggerTurnOn -------------------------- //
Analysis_Template_MC::Analysis_Template_MC(edm::ParameterSet const& cfg)
{
     mFileName       = cfg.getParameter<std::string>               ("filename");
     mTreeName       = cfg.getParameter<std::string>               ("treename");
     mDirName        = cfg.getParameter<std::string>               ("dirname");
     
     mIsMCarlo       = cfg.getUntrackedParameter<bool>             ("isMCarlo");

     mCrossSection   = cfg.getUntrackedParameter<double>             ("CrossSection");
     mIntLumi        = cfg.getUntrackedParameter<double>             ("IntLumi");

     mTriggers       = cfg.getUntrackedParameter<std::vector<std::string>>     ("Triggers");
}

//------------------------------ Declaration Of The Function beginjob() ------------------------//
void Analysis_Template_MC::beginJob()
 {

     mInf = TFile::Open(mFileName.c_str());
     mDir = (TDirectoryFile*)mInf->Get(mDirName.c_str());
     mTree = (TTree*)mDir->Get(mTreeName.c_str());
     
     //TBranch *branch = mTree->GetBranch("events");
     numerator1=0;
     numerator2=0;
     numerator3=0;
     denominator=0;

     numeratorParton1=0;
     numeratorParton2=0;
     numeratorParton3=0;
     denominatorParton=0;
     
     //------------------ Histogram Booking --------------------------- //
     num_of_Vtx     = fs->make<TH1F>("num_of_Vtx","num_of_Vtx",100,0.,100.);
     num_of_Jets     = fs->make<TH1F>("num_of_Jets","num_of_Jets",100,0.,100.);
     num_of_SubJets     = fs->make<TH1F>("num_of_SubJets","num_of_SubJets",100,0.,100.);
     num_of_SubBJets     = fs->make<TH1F>("num_of_SubBJets","num_of_SubBJets",100,0.,100.);
     num_of_BJets     = fs->make<TH1F>("num_of_BJets","num_of_BJets",100,0.,100.);
     num_of_GenJets     = fs->make<TH1F>("num_of_GenJets","num_of_GenJets",100,0.,100.);
     num_of_Leptons     = fs->make<TH1F>("num_of_Leptons","num_of_Leptons",100,0.,100.);
     num_of_GenLeptons     = fs->make<TH1F>("num_of_GenLeptons","num_of_GenLeptons",100,0.,100.);

     ptDETJet  = fs->make<TH1F>("ptDETJet","ptDETJet",200,0.,2000.); ptDETJet->Sumw2();
     yDETJet = fs->make<TH1F>("yDETJet","yDETJet",60,-3.,3.); yDETJet->Sumw2();
     phiDETJet = fs->make<TH1F>("phiDETJet","phiDETJet",60, -TMath::Pi(),TMath::Pi()); phiDETJet->Sumw2();

     ChfDETJet = fs->make<TH1F>("ChfDETJet","ChfDETJet",100,0.,1.);
     ElfDETJet = fs->make<TH1F>("ElfDETJet","ElfDETJet",100,0.,1.);
     PhfDETJet = fs->make<TH1F>("PhfDETJet","PhfDETJet",100,0.,1.);
     NhfDETJet = fs->make<TH1F>("NhfDETJet","NhfDETJet",100,0.,1.);
     MufDETJet = fs->make<TH1F>("MufDETJet","MufDETJet",100,0.,1.);

     Tau1DETJet = fs->make<TH1F>("Tau1DETJet","Tau1DETJet",100,0.,1.);
     Tau2DETJet = fs->make<TH1F>("Tau2DETJet","Tau2DETJet",100,0.,1.);
     Tau3DETJet = fs->make<TH1F>("Tau3DETJet","Tau3DETJet",100,0.,1.);
     Tau21DETJet = fs->make<TH1F>("Tau21DETJet","Tau21DETJet",100,0.,1.);
     Tau32DETJet = fs->make<TH1F>("Tau32DETJet","Tau32DETJet",100,0.,1.);

     EnergyDETJet= fs->make<TH1F>("EnergyDETJet","EnergyDETJet",100,0.,600.);
     MassDETJet= fs->make<TH1F>("MassDETJet","MassDETJet",100,0.,600.);
     MassSoftDropDETJet= fs->make<TH1F>("MassSoftDropDETJet","MassSoftDropDETJet",100,0.,600.);
     BtagDETJet= fs->make<TH1F>("BtagDETJet","BtagDETJet",100,0.,1.);
     
     PtDETSubJet0  = fs->make<TH1F>("PtDETSubJet0","PtDETSubJet0",200,0.,2000.); PtDETSubJet0->Sumw2();
     PtDETSubJet1  = fs->make<TH1F>("PtDETSubJet1","PtDETSubJet1",200,0.,2000.); PtDETSubJet1->Sumw2();
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

     ptGENJet  = fs->make<TH1F>("ptGENJet","ptGENJet",200,0.,2000.); ptGENJet->Sumw2();
     yGENJet = fs->make<TH1F>("yGENJet","yGENJet",60,-3.,3.); yGENJet->Sumw2();
     phiGENJet = fs->make<TH1F>("phiGENJet","phiGENJet",60, -TMath::Pi(),TMath::Pi()); phiGENJet->Sumw2();
     Tau1GENJet = fs->make<TH1F>("Tau1GENJet","Tau1GENJet",100,0.,1.);
     Tau2GENJet = fs->make<TH1F>("Tau2GENJet","Tau2GENJet",100,0.,1.);
     Tau3GENJet = fs->make<TH1F>("Tau3GENJet","Tau3GENJet",100,0.,1.);
     EnergyGENJet= fs->make<TH1F>("EnergyGENJet","EnergyGENJet",100,0.,600.);
     MassGENJet= fs->make<TH1F>("MassGENJet","MassGENJet",100,0.,600.);
     MassSoftDropGENJet= fs->make<TH1F>("MassSoftDropGENJet","MassSoftDropGENJet",100,0.,600.);

     BJetGenFlag = fs->make<TH1F>("BJetGenFlag","BJetGenFlag",2,-0.5,1.5);
     metGEN = fs->make<TH1F>("metGEN","metGEN",100,0.,300);
     metSigGEN = fs->make<TH1F>("metSigGEN","metSigGEN",100,0.,1.0);
     mvaGEN = fs->make<TH1F>("mvaGEN","mvaGEN",100,0.,1.0);

     metDET = fs->make<TH1F>("metDET","metDET",100,0.,300);
     metSigDET = fs->make<TH1F>("metSigDET","metSigDET",100,0.,1.0);
     mvaDET = fs->make<TH1F>("mvaDET","mvaDET",100,0.,1.0);

     ptTopParton  = fs->make<TH1F>("ptTopParton","ptTopParton",200,0.,2000.); ptTopParton->Sumw2();
     yTopParton = fs->make<TH1F>("yTopParton","yTopParton",60,-3.,3.); yTopParton->Sumw2();
     phiTopParton = fs->make<TH1F>("phiTopParton","phiTopParton",60, -TMath::Pi(),TMath::Pi()); phiTopParton->Sumw2();
     EnergyTopParton= fs->make<TH1F>("EnergyTopParton","EnergyTopParton",100,0.,600.);

     DeltaR_GenJetParton = fs->make<TH1F>("DeltaR_GenJetParton","DeltaR_GenJetParton",60, 0,2*TMath::Pi()); DeltaR_GenJetParton->Sumw2();
     DeltaR_DetJetParton = fs->make<TH1F>("DeltaR_DetJetParton","DeltaR_DetJetParton",60, 0,2*TMath::Pi()); DeltaR_DetJetParton->Sumw2();     
     DeltaR_GenJetRecoJet= fs->make<TH1F>("DeltaR_GenJetRecoJet","DeltaR_GenJetRecoJet",60, 0,2*TMath::Pi()); DeltaR_GenJetRecoJet->Sumw2();     

     ResolutionPtPartonGenJet = fs->make<TH1F>("ResolutionPtPartonGenJet","ResolutionPtPartonGenJet",60,-3.,3.); ResolutionPtPartonGenJet->Sumw2();
     ResolutionPtPartonDetJet = fs->make<TH1F>("ResolutionPtPartonDetJet","ResolutionPtPartonDetJet",60,-3.,3.); ResolutionPtPartonDetJet->Sumw2();
     ResolutionPtRecoGenJet = fs->make<TH1F>("ResolutionPtRecoGenJet","ResolutionPtRecoGenJet",60,-3.,3.); ResolutionPtRecoGenJet->Sumw2();

     ResolutionMassRecoGenJet = fs->make<TH1F>("ResolutionMassRecoGenJet","ResolutionMassRecoGenJet",60,-3.,3.); ResolutionMassRecoGenJet->Sumw2();
     ResolutionSDMassRecoGenJet = fs->make<TH1F>("ResolutionSDMassRecoGenJet","ResolutionSDMassRecoGenJet",60,-3.,3.); ResolutionSDMassRecoGenJet->Sumw2();
     ResolutionSubJetMassRecoGenJet = fs->make<TH1F>("ResolutionSubJetMassRecoGenJet","ResolutionSubJetMassRecoGenJet",60,-3.,3.); ResolutionSubJetMassRecoGenJet->Sumw2();

     ResponsePtRecoGenJet = fs->make<TH2F>("ResponsePtRecoGenJet","ResponsePtRecoGenJet",40,0.,2000.,40,0.,2000.); ResponsePtRecoGenJet->Sumw2();
     ResponsePtPartonGenJet = fs->make<TH2F>("ResponsePtPartonGenJet","ResponsePtPartonGenJet",40,0.,2000.,40,0.,2000.); ResponsePtPartonGenJet->Sumw2();

     ResolutionPhiPartonGenJet = fs->make<TH1F>("ResolutionPhiPartonGenJet","ResolutionPhiPartonGenJet",60,-3.,3.); ResolutionPhiPartonGenJet->Sumw2();
     ResolutionPhiPartonDetJet = fs->make<TH1F>("ResolutionPhiPartonDetJet","ResolutionPhiPartonDetJet",60,-3.,3.); ResolutionPhiPartonDetJet->Sumw2();
     ResolutionPhiRecoGenJet = fs->make<TH1F>("ResolutionPhiRecoGenJet","ResolutionPhiRecoGenJet",60,-3.,3.); ResolutionPhiRecoGenJet->Sumw2();

     ResponsePhiPartonGenJet = fs->make<TH2F>("ResponsePhiPartonGenJet","ResponsePhiPartonGenJet",60, 0,2*TMath::Pi(),60, 0,2*TMath::Pi()); ResponsePhiPartonGenJet->Sumw2();
     ResponsePhiRecoGenJet = fs->make<TH2F>("ResponsePhiRecoGenJet","ResponsePhiRecoGenJet",60, 0,2*TMath::Pi(),60, 0,2*TMath::Pi()); ResponsePhiRecoGenJet->Sumw2();

     PtGENSubJet0  = fs->make<TH1F>("PtGENSubJet0","PtGENSubJet0",200,0.,2000.); PtGENSubJet0->Sumw2();
     PtGENSubJet1  = fs->make<TH1F>("PtGENSubJet1","PtGENSubJet1",200,0.,2000.); PtGENSubJet1->Sumw2();
     EtaGENSubJet0 = fs->make<TH1F>("EtaGENSubJet0","EtaGENSubJet0",60,-3.,3.); EtaGENSubJet0->Sumw2();
     EtaGENSubJet1 = fs->make<TH1F>("EtaGENSubJet1","EtaGENSubJet1",60,-3.,3.); EtaGENSubJet1->Sumw2();
     MassGENSubJet0 = fs->make<TH1F>("MassGENSubJet0","MassGENSubJet0",200,0.,2000.); MassGENSubJet0->Sumw2();
     MassGENSubJet1 = fs->make<TH1F>("MassGENSubJet1","MassGENSubJet1",200,0.,2000.); MassGENSubJet1->Sumw2();
     PhiGENSubJet0 = fs->make<TH1F>("PhiGENSubJet0","PhiGENSubJet0",60,-TMath::Pi(),TMath::Pi()); PhiGENSubJet0->Sumw2();
     PhiGENSubJet1 = fs->make<TH1F>("PhiGENSubJet1","PhiGENSubJet1",60,-TMath::Pi(),TMath::Pi()); PhiGENSubJet1->Sumw2();

     DeltaRGENSubJets= fs->make<TH1F>("DeltaRGENSubJets","DeltaRGENSubJets",100,0.,6.28);
     MuGENSubJets = fs->make<TH1F>("MuGENSubJets","MuGENSubJets",100,0.,1.0);

 } // end of function beginJob()





 //------------------------ endjob() function declaration ---------------------- //
 void Analysis_Template_MC::endJob()
 {
   mInf->Close();
   
   cout<<numerator1<<" "<<numerator2<<" "<<numerator3<<" "<<denominator<<endl;
   cout<<numeratorParton1<<" "<<numeratorParton2<<" "<<numeratorParton3<<" "<<denominatorParton<<endl;

   //9295 7567 4182 10827
   //110433 66105 22777 121986


 } // closing endJob()





 //--------------------------- analyze() fuction declaration ------------------ //
void Analysis_Template_MC::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup)
 {

   unsigned NEntries = mTree->GetEntries();
   cout<<"Reading TREE: "<<NEntries<<" events"<<endl;
   
   int decade = 0 ;
   
   float hweight=1.;  ///Initial value to one

   int lumi_=-100; int nvtx_=-100; int nBJets_=-100; int nJets_=-100; int nGenJets_=-100; int nLeptons_=-100; std::vector<int>* jetNSub_=0; std::vector<int>* jetNBSub_=0; int nGenLeptons_=-100; std::vector<int>* jetNGenSub_=0;
   float met_=-100; float mva_=-100; float ht_=-100; float mJJ_=-100; float yJJ_=-100; float ptJJ_=-100; float dRJJ_=-100; float metSig_=-100;
   std::vector<float>* jetPt_=0; std::vector<float>* jetEta_=0; std::vector<float>* jetBtag_=0; std::vector<float>* jetPhi_=0; std::vector<float>* jetEnergy_=0; std::vector<float>* jetMass_=0; std::vector<float>* jetMassSoftDrop_=0;

   std::vector<float>* jetchf_=0; std::vector<float>* jetnhf_=0; std::vector<float>* jetphf_=0; std::vector<float>* jetmuf_=0; std::vector<float>* jetelf_=0; std::vector<float>* jettau1_=0; std::vector<float>* jettau2_=0; std::vector<float>* jettau3_=0; 

   std::vector<float>* GenJetPt_=0; std::vector<float>* GenJetEta_=0; std::vector<float>* GenJetPhi_=0; std::vector<float>* GenJetenergy_=0; std::vector<float>* GenJetmass_=0; std::vector<float>* GenJetMassSoftDrop_=0; std::vector<float>* GenJettau1_=0; std::vector<float>* GenJettau2_=0; std::vector<float>* GenJettau3_=0; 

   std::vector<double>* pre_=0;   std::vector<bool>* bit_=0;

   float metGen_=-100; float mvaGen_=-100; float metSigGen_=-100; 
   std::vector<bool>* isBJetGen_=0;

   std::vector<float>* partonPt_=0; std::vector<float>* partonEta_=0; std::vector<float>* partonPhi_=0; std::vector<float>* partonEnergy_=0;

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

   mTree->SetBranchAddress("lumi",&lumi_);
   mTree->SetBranchAddress("nvtx",&nvtx_);
   mTree->SetBranchAddress("met",&met_);
   mTree->SetBranchAddress("metSig",&metSig_);
   mTree->SetBranchAddress("mva",&mva_);
   mTree->SetBranchAddress("ht",&ht_);
   mTree->SetBranchAddress("nBJets",&nBJets_);
   mTree->SetBranchAddress("nJets",&nJets_);
   mTree->SetBranchAddress("nLeptons",&nLeptons_);
   mTree->SetBranchAddress("mJJ",&mJJ_);
   mTree->SetBranchAddress("yJJ",&yJJ_);
   mTree->SetBranchAddress("ptJJ",&ptJJ_);
   mTree->SetBranchAddress("dRJJ",&dRJJ_);
   mTree->SetBranchAddress("jetNSub",&jetNSub_);
   mTree->SetBranchAddress("jetNSubGen",&jetNGenSub_);
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
   mTree->SetBranchAddress("jetFlavorSub0",&jetflavorSub0_);
   mTree->SetBranchAddress("jetFlavorSub1",&jetflavorSub1_);

   if(mIsMCarlo){
     mTree->SetBranchAddress("nGenJets",&nGenJets_);
     mTree->SetBranchAddress("GenJetpt",&GenJetPt_);
     mTree->SetBranchAddress("GenJeteta",&GenJetEta_);
     mTree->SetBranchAddress("GenJetphi",&GenJetPhi_);
     mTree->SetBranchAddress("GenJetmass",&GenJetmass_);
     mTree->SetBranchAddress("GenJetenergy",&GenJetenergy_);
     mTree->SetBranchAddress("GenSoftDropMass",&GenJetMassSoftDrop_);
     mTree->SetBranchAddress("GenJettau1",&GenJettau1_);
     mTree->SetBranchAddress("GenJettau2",&GenJettau2_);
     mTree->SetBranchAddress("GenJettau3",&GenJettau3_);
     mTree->SetBranchAddress("metGen",&metGen_);
     mTree->SetBranchAddress("metGenSig",&metSigGen_);
     mTree->SetBranchAddress("mvaGen",&mvaGen_);
     mTree->SetBranchAddress("nGenLeptons",&nGenLeptons_);
     mTree->SetBranchAddress("isBJetGen",&isBJetGen_);

     mTree->SetBranchAddress("partonPt",&partonPt_);
     mTree->SetBranchAddress("partonEta",&partonEta_);
     mTree->SetBranchAddress("partonPhi",&partonPhi_);
     mTree->SetBranchAddress("partonE",&partonEnergy_);

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
   
   mTree->SetBranchAddress("triggerPre",&pre_);
   mTree->SetBranchAddress("triggerBit",&bit_);

   TH1D* NumberEvents=(TH1D*) mInf->Get("boosted/TriggerPass");

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
   
   for(unsigned  l=0; l<NEntries; l++) {
     
    //----------- progress report -------------
    double progress = 10.0*l/(1.0*NEntries);
    int k = TMath::FloorNint(progress);
    if (k > decade)
      cout<<10*k<<" %"<<endl;
    decade = k;
   
    //----------- read the event --------------
    mTree->GetEntry(l);
    
    double NEvents=NumberEvents->GetBinContent(1);
    bool triggerfired=false;

    if(mIsMCarlo){
      triggerfired=true;
      hweight=mCrossSection/NEvents;
    }

    if(!mIsMCarlo){
      if(bit_->at(5)) triggerfired=true;
      double prescale=pre_->at(5);
      hweight=prescale/mIntLumi;
    }

    if(triggerfired==false) continue;

    num_of_Vtx->Fill(nvtx_,hweight);
    num_of_Jets->Fill(nJets_,hweight);

    metDET->Fill(met_,hweight);
    metSigDET->Fill(metSig_,hweight);
    mvaDET->Fill(mva_,hweight);

    if(mIsMCarlo){ num_of_GenJets->Fill(nGenJets_,hweight); 
      metGEN->Fill(metGen_,hweight);
      metSigGEN->Fill(metSigGen_,hweight);
      mvaGEN->Fill(mvaGen_,hweight);
      num_of_GenLeptons->Fill(nGenLeptons_,hweight); 
    }
    num_of_BJets->Fill(nBJets_,hweight);
    num_of_Leptons->Fill(nLeptons_,hweight);
    
    for(int j=0; j< nJets_; j++){

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
      double tau21=-1;
      if(jettau1_->at(j)!=0){ tau21= jettau2_->at(j)/jettau1_->at(j);}
      double tau32=-1;
      if(jettau2_->at(j)!=0){ tau32= jettau3_->at(j)/jettau2_->at(j);}
      Tau21DETJet->Fill(tau21,hweight);
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
      FlavourDETSubJet0->Fill(jetflavorSub0_->at(j),hweight);
      FlavourDETSubJet1->Fill(jetflavorSub1_->at(j),hweight);
      
      if(false){
	printf("Number of PFJets=%i\n",nJets_);
	printf("j=%2i  pt=%8.3f  y=%6.3f  phi=%6.3f\n",j,jetPt_->at(j),jetEta_->at(j),jetPhi_->at(j));
      }

    }

    if(mIsMCarlo){
      bool PreliminarySubJetMassSwitch=true;
      //Gen jet part
      for(int j=0; j< nGenJets_; j++){
	ptGENJet->Fill(GenJetPt_->at(j),hweight);
	yGENJet->Fill(GenJetEta_->at(j),hweight);
	phiGENJet->Fill(GenJetPhi_->at(j),hweight);
	EnergyGENJet->Fill(GenJetenergy_->at(j),hweight);
	MassGENJet->Fill(GenJetmass_->at(j),hweight);
	MassSoftDropGENJet->Fill(GenJetMassSoftDrop_->at(j),hweight);
	Tau1GENJet->Fill(GenJettau1_->at(j),hweight);
	Tau2GENJet->Fill(GenJettau2_->at(j),hweight);
	Tau3GENJet->Fill(GenJettau3_->at(j),hweight);
	BJetGenFlag->Fill(isBJetGen_->at(j),hweight);

	if(jetNGenSub_->at(j)==2){
	  PtGENSubJet0->Fill(GenSubJet1Pt_->at(j),hweight);
	  PtGENSubJet1->Fill(GenSubJet2Pt_->at(j),hweight);
	  EtaGENSubJet0->Fill(GenSubJet1Eta_->at(j),hweight);
	  EtaGENSubJet1->Fill(GenSubJet2Eta_->at(j),hweight);
	  PhiGENSubJet0->Fill(GenSubJet1Phi_->at(j),hweight);
	  PhiGENSubJet1->Fill(GenSubJet2Phi_->at(j),hweight);
	  MassGENSubJet0->Fill(GenSubJet1Mass_->at(j),hweight);
	  MassGENSubJet1->Fill(GenSubJet2Mass_->at(j),hweight);
	}
	else{ PreliminarySubJetMassSwitch=false;}
	  
	DeltaRGENSubJets->Fill(GenSubJetsDeltaR_->at(j),hweight);
	MuGENSubJets->Fill(GenSubJetsMu_->at(j),hweight);

      }

      //Match with partons
      for(unsigned j=0; j<partonPt_->size();j++){
	ptTopParton->Fill(partonPt_->at(j),hweight);
	yTopParton->Fill(partonEta_->at(j),hweight);
	phiTopParton->Fill(partonPhi_->at(j),hweight);
	EnergyTopParton->Fill(partonEnergy_->at(j),hweight);
	
	for(int i=0; i< nGenJets_; i++){
	  if(partonPt_->at(j)>500 && GenJetPt_->at(i)>350){
	    double deltaR_jt=DeltaR(partonEta_->at(j),partonPhi_->at(j),GenJetEta_->at(i),GenJetPhi_->at(i));
	    DeltaR_GenJetParton->Fill(deltaR_jt,hweight);
	    if(deltaR_jt<0.8){
	      double resolutionPt=(GenJetPt_->at(i)-partonPt_->at(j))/partonPt_->at(j);
	      double resolutionPhi=GenJetPhi_->at(i)-partonPhi_->at(j);
	      
	      ResolutionPtPartonGenJet->Fill(resolutionPt,hweight);
	      ResolutionPhiPartonGenJet->Fill(resolutionPhi,hweight);

	      ResponsePtPartonGenJet->Fill(GenJetPt_->at(i),partonPt_->at(j),hweight);
	      ResponsePhiPartonGenJet->Fill(GenJetPhi_->at(i),partonPhi_->at(j),hweight);
	      
	      denominatorParton++;

	      if(PreliminarySubJetMassSwitch) 
		{
		  if(isBJetGen_->at(i)) numeratorParton1++;
		  if(isBJetGen_->at(i) && GenJetmass_->at(i)>150 && GenJetmass_->at(i)<200) numeratorParton2++;
		  if(isBJetGen_->at(i) && GenJetmass_->at(i)>150 && GenJetmass_->at(i)<200 && GenSubJet1Mass_->at(i) < 90 && GenSubJet1Mass_->at(i)>70) numeratorParton3++;
		}
	      
	    }
	  }
	}
	
	for(int i=0; i< nJets_; i++){
	  if(partonPt_->at(j)>500 && jetPt_->at(i)>350){
	    double deltaR_jt=DeltaR(partonEta_->at(j),partonPhi_->at(j),jetEta_->at(i),jetPhi_->at(i));
	    DeltaR_DetJetParton->Fill(deltaR_jt,hweight);
	    if(deltaR_jt<0.8){
	      double resolutionPt=(jetPt_->at(i)-partonPt_->at(j))/partonPt_->at(j);
	      double resolutionPhi=jetPhi_->at(i)-partonPhi_->at(j);
	      
	      ResolutionPtPartonDetJet->Fill(resolutionPt,hweight);
	      ResolutionPhiPartonDetJet->Fill(resolutionPhi,hweight);
	    }
	  }
	}
      }
      

      //Preselection at RECO level
      int btagjets=0;

      for(int j=0; j< nJets_; j++){
	if(jetPt_->at(j)>400 && fabs(jetEta_->at(j))<2.4 && jetMassSoftDrop_->at(j)>50){//mass drop+b-tag
	  if(jetbtagSub0_->at(j)>0.814 || jetbtagSub1_->at(j)>0.814){ //b-tagged subjets
	    btagjets++;
	  }
	}
      }

      bool preselection=false;
      if(btagjets>1) preselection=true;

      if(preselection){
	double tau21_1jet=-1;
	if(jettau1_->at(0)!=0){ tau21_1jet= jettau2_->at(0)/jettau1_->at(0);}
	double tau32_1jet=-1;
	if(jettau2_->at(0)!=0){ tau32_1jet= jettau3_->at(0)/jettau2_->at(0);}
	
	double tau21_2jet=-1;
	if(jettau1_->at(1)!=0){ tau21_2jet= jettau2_->at(1)/jettau1_->at(1);}
	double tau32_2jet=-1;
	if(jettau2_->at(1)!=0){ tau32_2jet= jettau3_->at(1)/jettau2_->at(1);}
	

	double DeltaRSubJet_1jet=DeltaR(jetetaSub0_->at(0),jetphiSub0_->at(0),jetetaSub1_->at(0),jetphiSub1_->at(0));
	double DeltaRSubJet_2jet=DeltaR(jetetaSub0_->at(1),jetphiSub0_->at(1),jetetaSub1_->at(1),jetphiSub1_->at(1));
	
	//cout<<"QCD:-1 TOP:1 "<<jetPt_->at(0)<<" "<<jetMassSoftDrop_->at(0)<<" "<<jetMass_->at(0)<<" "<<jetmassSub0_->at(0)<<" "<<jetmassSub1_->at(0)<<" "<<jettau1_->at(0)<<" "<<jettau2_->at(0)<<" "<<jettau3_->at(0)<<tau21_1jet<<" "<<tau32_1jet<<" "<<jetptSub0_->at(0)<<" "<<jetptSub1_->at(0)<<" "<<DeltaRSubJet_1jet<<" "<<jetPt_->at(1)<<" "<<jetMassSoftDrop_->at(1)<<" "<<jetMass_->at(1)<<" "<<jetmassSub0_->at(1)<<" "<<jetmassSub1_->at(1)<<" "<<jettau1_->at(1)<<" "<<jettau2_->at(1)<<" "<<jettau3_->at(1)<<tau21_2jet<<" "<<tau32_2jet<<" "<<jetptSub0_->at(1)<<" "<<jetptSub1_->at(1)<<" "<<DeltaRSubJet_2jet<<endl
;

	double tau21Leading = tau21_1jet;
	double tau32Leading = tau32_1jet;
	double SDmassLeading = jetMassSoftDrop_->at(0);
	double MassSubJet0Leading = jetmassSub0_->at(0);
	double MassSubJet1Leading = jetmassSub0_->at(1);
	double jetSubJetpt0Leading = jetptSub0_->at(0);
	double jetSubJetpt1Leading = jetptSub0_->at(1);
	
	double tau21SubLeading = tau21_2jet;
	double tau32SubLeading = tau32_2jet;
	double SDmassSubLeading = jetMassSoftDrop_->at(0);
	double MassSubJet0SubLeading = jetmassSub1_->at(0);
	double MassSubJet1SubLeading = jetmassSub1_->at(1);
	double jetSubJetpt0SubLeading = jetptSub1_->at(0);
	double jetSubJetpt1SubLeading = jetptSub1_->at(1);

      }
	

      //Match with RecoJets
      for(int j=0; j< nJets_; j++){
	for(int i=0; i< nGenJets_; i++){
	  if(preselection){
	    double deltaR_jt=DeltaR(jetEta_->at(j),jetPhi_->at(j),GenJetEta_->at(i),GenJetPhi_->at(i));
	    DeltaR_GenJetRecoJet->Fill(deltaR_jt,hweight);
	    if(deltaR_jt<0.8){
	      double resolutionPt=(jetPt_->at(j)-GenJetPt_->at(i))/GenJetPt_->at(i);
	      double resolutionPhi=jetPhi_->at(j)-GenJetPhi_->at(i);
	      
	      ResolutionPtRecoGenJet->Fill(resolutionPt,hweight);
	      ResolutionPhiRecoGenJet->Fill(resolutionPhi,hweight);

	      ResponsePtRecoGenJet->Fill(GenJetPt_->at(i),jetPt_->at(j),hweight);
	      ResponsePhiRecoGenJet->Fill(GenJetPhi_->at(i),jetPhi_->at(j),hweight);

	      double resolutionMass=(jetMass_->at(j)-GenJetmass_->at(i))/GenJetmass_->at(i);
	      ResolutionMassRecoGenJet->Fill(resolutionMass,hweight);
	      
	      double resolutionSDMass=(jetMassSoftDrop_->at(j)-GenJetMassSoftDrop_->at(i))/GenJetMassSoftDrop_->at(i);
	      ResolutionSDMassRecoGenJet->Fill(resolutionSDMass,hweight);
	      
	      if(PreliminarySubJetMassSwitch){
		double resolutionSubGenJetMass=(jetmassSub0_->at(j)-GenSubJet1Mass_->at(i))/GenSubJet1Mass_->at(i);
		ResolutionSubJetMassRecoGenJet->Fill(resolutionSubGenJetMass,hweight);
	      }
	    }
	    
	    if(PreliminarySubJetMassSwitch){
	      if(jetBtag_->at(j)>0.814 && jetMass_->at(j)>150 && jetMass_->at(j)<200 && jetmassSub0_->at(j) < 90 && jetmassSub0_->at(j)>70){
		denominator++;
		if(isBJetGen_->at(i)) numerator1++;
		if(isBJetGen_->at(i) && GenJetmass_->at(i)>150 && GenJetmass_->at(i)<200) numerator2++;
		if(isBJetGen_->at(i) && GenJetmass_->at(i)>150 && GenJetmass_->at(i)<200 && GenSubJet1Mass_->at(i) < 90 && GenSubJet1Mass_->at(i)>70) numerator3++;
	      }
	    }
	  }
	}
      }
    }
    



   } // end of event loop


} // closing analyze() function

Analysis_Template_MC::~Analysis_Template_MC()
{
}


DEFINE_FWK_MODULE(Analysis_Template_MC);

