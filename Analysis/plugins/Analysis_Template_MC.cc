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
     
     mCrossSection   = cfg.getUntrackedParameter<double>             ("CrossSection");
     mIntLumi        = cfg.getUntrackedParameter<double>             ("IntLumi",1);

}

//------------------------------ Declaration Of The Function beginjob() ------------------------//
void Analysis_Template_MC::beginJob()
 {

     mInf = TFile::Open(mFileName.c_str());
     mDir = (TDirectoryFile*)mInf->Get(mDirName.c_str());
     mTree = (TTree*)mInf->Get(mTreeName.c_str());
     NEvents=0;
     
     //------------------ Histogram Booking --------------------------- //
     num_of_Vtx     = fs->make<TH1F>("num_of_Vtx","num_of_Vtx",100,0.,100.);
     num_of_Jets     = fs->make<TH1F>("num_of_Jets","num_of_Jets",100,0.,100.);
     num_of_GenJets     = fs->make<TH1F>("num_of_GenJets","num_of_GenJets",100,0.,100.);

     ptDETJet  = fs->make<TH1F>("ptDETJet","ptDETJet",200,0.,2000.); ptDETJet->Sumw2();
     yDETJet = fs->make<TH1F>("yDETJet","yDETJet",60,-3.,3.); yDETJet->Sumw2();
     phiDETJet = fs->make<TH1F>("phiDETJet","phiDETJet",60, -TMath::Pi(),TMath::Pi()); phiDETJet->Sumw2();

     DeltaPhiDETJet = fs->make<TH1F>("DeltaPhiDETJet","DeltaPhiDETJet",30, 0 ,TMath::Pi()); DeltaPhiDETJet->Sumw2();

     //Tau1DETJet = fs->make<TH1F>("Tau1DETJet","Tau1DETJet",100,0.,1.);
     //Tau2DETJet = fs->make<TH1F>("Tau2DETJet","Tau2DETJet",100,0.,1.);
     //Tau3DETJet = fs->make<TH1F>("Tau3DETJet","Tau3DETJet",100,0.,1.);
     //Tau21DETJet = fs->make<TH1F>("Tau21DETJet","Tau21DETJet",100,0.,1.);
     //Tau32DETJet = fs->make<TH1F>("Tau32DETJet","Tau32DETJet",100,0.,1.);

     MassDETJet= fs->make<TH1F>("MassDETJet","MassDETJet",100,0.,600.);
     MassSoftDropDETJet= fs->make<TH1F>("MassSoftDropDETJet","MassSoftDropDETJet",100,0.,600.);

 } // end of function beginJob()


 //------------------------ endjob() function declaration ---------------------- //
 void Analysis_Template_MC::endJob()
 {
   mInf->Close();
   
 } // closing endJob()

 //--------------------------- analyze() fuction declaration ------------------ //
void Analysis_Template_MC::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup)
 {

   unsigned NEntries = mTree->GetEntries();
   cout<<"Reading TREE: "<<NEntries<<" events"<<endl;
   
   int decade = 0 ;
   
   float hweight= mCrossSection;  //Histograms already scaled by the cross section

   std::vector<float>* jetPt_=0; std::vector<float>* jetEta_=0; std::vector<float>* jetPhi_=0; std::vector<float>* jetMass_=0; std::vector<float>* jetFlavour_=0;   std::vector<float>* jetBTag_=0; 
   std::vector<float>* AK8jetPt_=0; std::vector<float>* AK8jetEta_=0; std::vector<float>* AK8jetPhi_=0; std::vector<float>* AK8jetMass_=0; std::vector<float>* AK8jetFlavour_=0;   std::vector<float>* AK8jetBTag_=0; 

   std::vector<float>* AK8jetTau1_=0;    std::vector<float>* AK8jetTau2_=0;    std::vector<float>* AK8jetTau3_=0; 
   std::vector<float>* AK8genjetTau1_=0;    std::vector<float>* AK8genjetTau2_=0;    std::vector<float>* AK8genjetTau3_=0; 

   //subjets
   std::vector<float>* AK8jetSD1Pt_=0; std::vector<float>* AK8jetSD1Eta_=0; std::vector<float>* AK8jetSD1Phi_=0; std::vector<float>* AK8jetSD1Mass_=0; 
   std::vector<float>* AK8jetSD2Pt_=0; std::vector<float>* AK8jetSD2Eta_=0; std::vector<float>* AK8jetSD2Phi_=0; std::vector<float>* AK8jetSD2Mass_=0; 

   //gen jets
   std::vector<float>* genjetPt_=0; std::vector<float>* genjetEta_=0; std::vector<float>* genjetPhi_=0; std::vector<float>* genjetMass_=0; std::vector<float>* genjetFlavour_=0; 
   std::vector<float>* AK8genjetPt_=0; std::vector<float>* AK8genjetEta_=0; std::vector<float>* AK8genjetPhi_=0; std::vector<float>* AK8genjetMass_=0; std::vector<float>* AK8genjetFlavour_=0; 


   //jets AK4 CHS
   mTree->SetBranchAddress("jetPt",&jetPt_);
   mTree->SetBranchAddress("jetEta",&jetEta_);
   mTree->SetBranchAddress("jetPhi",&jetPhi_);
   mTree->SetBranchAddress("jetMass",&jetMass_);
   mTree->SetBranchAddress("jetFlavour",&jetFlavour_);
   mTree->SetBranchAddress("jetBTag",&jetBTag_);

   //genjets AK4
   mTree->SetBranchAddress("genjetPt",&genjetPt_);
   mTree->SetBranchAddress("genjetEta",&genjetEta_);
   mTree->SetBranchAddress("genjetPhi",&genjetPhi_);
   mTree->SetBranchAddress("genjetMass",&genjetMass_);
   mTree->SetBranchAddress("genjetFlavour",&genjetFlavour_);

   //genjets AK8
   mTree->SetBranchAddress("AK8genjetPt",&AK8genjetPt_);
   mTree->SetBranchAddress("AK8genjetEta",&AK8genjetEta_);
   mTree->SetBranchAddress("AK8genjetPhi",&AK8genjetPhi_);
   mTree->SetBranchAddress("AK8genjetMass",&AK8genjetMass_);
   mTree->SetBranchAddress("AK8genjetFlavour",&AK8genjetFlavour_);

   //jets AK8 CHS
   mTree->SetBranchAddress("AK8jetPt",&AK8jetPt_);
   mTree->SetBranchAddress("AK8jetEta",&AK8jetEta_);
   mTree->SetBranchAddress("AK8jetPhi",&AK8jetPhi_);
   mTree->SetBranchAddress("AK8jetMass",&AK8jetMass_);
   mTree->SetBranchAddress("AK8jetFlavour",&AK8jetFlavour_);
   mTree->SetBranchAddress("AK8jetBTag",&AK8jetBTag_);

   //subjettiness AK8
   mTree->SetBranchAddress("AK8jetTau1",&AK8jetTau1_);
   mTree->SetBranchAddress("AK8jetTau2",&AK8jetTau2_);
   mTree->SetBranchAddress("AK8jetTau3",&AK8jetTau3_);

   mTree->SetBranchAddress("AK8genjetTau1",&AK8genjetTau1_);
   mTree->SetBranchAddress("AK8genjetTau2",&AK8genjetTau2_);
   mTree->SetBranchAddress("AK8genjetTau3",&AK8genjetTau3_);

   //subjets AK8
   mTree->SetBranchAddress("AK8jetSD1Pt",&AK8jetSD1Pt_); 
   mTree->SetBranchAddress("AK8jetSD1Eta",&AK8jetSD1Eta_);
   mTree->SetBranchAddress("AK8jetSD1Phi",&AK8jetSD1Phi_);
   mTree->SetBranchAddress("AK8jetSD1Mass",&AK8jetSD1Mass_);

   mTree->SetBranchAddress("AK8jetSD2Pt",&AK8jetSD2Pt_); 
   mTree->SetBranchAddress("AK8jetSD2Eta",&AK8jetSD2Eta_);
   mTree->SetBranchAddress("AK8jetSD2Phi",&AK8jetSD2Phi_);
   mTree->SetBranchAddress("AK8jetSD2Mass",&AK8jetSD2Mass_);

   //gen information
   float nVtx_; float scale_; float processID_;

   mTree->SetBranchAddress("nVtx",&nVtx_);
   mTree->SetBranchAddress("scale",&scale_);
   mTree->SetBranchAddress("processID",&processID_);   

   float metEt_;    float metEta_;    float metPhi_;
   mTree->SetBranchAddress("metEt",&metEt_);
   //mTree->SetBranchAddress("metEta",&metEta_); //not there yet
   mTree->SetBranchAddress("metPhi",&metPhi_);
   
   for(unsigned  l=0; l<NEntries; l++) {
     
    //----------- progress report -------------
    double progress = 10.0*l/(1.0*NEntries);
    int k = TMath::FloorNint(progress);
    if (k > decade)
      cout<<10*k<<" %"<<endl;
    decade = k;
   
    //----------- read the event --------------
    mTree->GetEntry(l);
    
    for(unsigned j=0; j< jetPt_->size(); j++){

      ptDETJet->Fill(jetPt_->at(j),hweight);
      yDETJet->Fill(jetEta_->at(j),hweight);
      phiDETJet->Fill(jetPhi_->at(j),hweight);
      MassDETJet->Fill(jetMass_->at(j),hweight);
      //MassSoftDropDETJet->Fill(jetMassSoftDrop_->at(j),hweight);
      //to add this information, you have to add the branches in the analyze method
      //Tau1DETJet->Fill(jettau1_->at(j),hweight);
      //Tau2DETJet->Fill(jettau2_->at(j),hweight);
      //Tau3DETJet->Fill(jettau3_->at(j),hweight);
      //double tau21=-1;
      //if(jettau1_->at(j)!=0){ tau21= jettau2_->at(j)/jettau1_->at(j);}
      //double tau32=-1;
      //if(jettau2_->at(j)!=0){ tau32= jettau3_->at(j)/jettau2_->at(j);}
      //Tau21DETJet->Fill(tau21,hweight);
      //Tau32DETJet->Fill(tau32,hweight);
    }

    //Gen jet information
    //for(int j=0; j< nGenJets_; j++){ //example histograms
	//ptGENJet->Fill(GenJetPt_->at(j),hweight);
	//yGENJet->Fill(GenJetEta_->at(j),hweight);
	//phiGENJet->Fill(GenJetPhi_->at(j),hweight);
	//EnergyGENJet->Fill(GenJetenergy_->at(j),hweight);
	//MassGENJet->Fill(GenJetmass_->at(j),hweight);
	//MassSoftDropGENJet->Fill(GenJetMassSoftDrop_->at(j),hweight);
	//Tau1GENJet->Fill(GenJettau1_->at(j),hweight);
	//Tau2GENJet->Fill(GenJettau2_->at(j),hweight);
	//Tau3GENJet->Fill(GenJettau3_->at(j),hweight);
    //}
    
   } // end of event loop

   if(jetPt_->size()>1){
     double deltaPhiTwoJets = DeltaPhi(jetPhi_->at(0),jetPhi_->at(1));
     DeltaPhiDETJet->Fill(deltaPhiTwoJets,hweight);
   }
   
 } // closing analyze() function

Analysis_Template_MC::~Analysis_Template_MC()
{
}


DEFINE_FWK_MODULE(Analysis_Template_MC);
