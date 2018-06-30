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
     
     mIsMCarlo       = cfg.getUntrackedParameter<bool>             ("isMCarlo",true);

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
   
   float hweight=1.;  //Initial value set to one

   std::vector<float>* jetPt_=0; std::vector<float>* jetEta_=0; std::vector<float>* jetPhi_=0; std::vector<float>* jetMass_=0; std::vector<float>* jetFlavour_=0;
   //std::vector<float>* jetMassSoftDrop_=0; 

   mTree->SetBranchAddress("jetPt",&jetPt_);
   mTree->SetBranchAddress("jetEta",&jetEta_);
   mTree->SetBranchAddress("jetPhi",&jetPhi_);
   mTree->SetBranchAddress("jetMass",&jetMass_);
   mTree->SetBranchAddress("jetFlavour",&jetFlavour_);
   mTree->SetBranchAddress("jetBTag",&jetBTag_);
   //mTree->SetBranchAddress("jetTau1",&jetTau1_);
   //mTree->SetBranchAddress("jetTau2",&jetTau2_);
   //mTree->SetBranchAddress("jetTau3",&jetTau3_);
   
   //To be added when they are inside the Ntuples
   //mTree->SetBranchAddress("jetMassSoftDrop",&jetMassSoftDrop_); 
   
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
   }
   
 } // closing analyze() function

Analysis_Template_MC::~Analysis_Template_MC()
{
}


DEFINE_FWK_MODULE(Analysis_Template_MC);
