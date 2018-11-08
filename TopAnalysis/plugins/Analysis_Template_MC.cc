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
          
     //------------------ Histogram Booking --------------------------- //

     num_of_GenJets     = fs->make<TH1F>("num_of_GenJets","num_of_GenJets",100,0.,100.);
     ptGENJet  = fs->make<TH1F>("ptGENJet","ptGENJet",200,0.,2000.); ptGENJet->Sumw2();
     yGENJet = fs->make<TH1F>("yGENJet","yGENJet",60,-3.,3.); yGENJet->Sumw2();
     phiGENJet = fs->make<TH1F>("phiGENJet","phiGENJet",60, -TMath::Pi(),TMath::Pi()); phiGENJet->Sumw2();

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

   int nGenJets_=-100; 
   std::vector<float>* GenJetPt_=0; std::vector<float>* GenJetEta_=0; std::vector<float>* GenJetPhi_=0; 
	
   mTree->SetBranchAddress("nGenJets",&nGenJets_);
   mTree->SetBranchAddress("GenJetpt",&GenJetPt_);
   mTree->SetBranchAddress("GenJeteta",&GenJetEta_);
   mTree->SetBranchAddress("GenJetphi",&GenJetPhi_);
      
   for(unsigned  l=0; l<NEntries; l++) {
     
     //----------- progress report -------------
     double progress = 10.0*l/(1.0*NEntries);
     int k = TMath::FloorNint(progress);
     if (k > decade)
       cout<<10*k<<" %"<<endl;
     decade = k;
     
     //----------- read the event --------------
     mTree->GetEntry(l);
     
     if(mIsMCarlo){
       hweight=mCrossSection/NEntries;
     }
     
     num_of_GenJets->Fill(nGenJets_,hweight); 
     
     //Gen jet information
     for(int j=0; j< nGenJets_; j++){ //example histograms
       ptGENJet->Fill(GenJetPt_->at(j),hweight);
       yGENJet->Fill(GenJetEta_->at(j),hweight);
       phiGENJet->Fill(GenJetPhi_->at(j),hweight);
     }
    
   } // end of event loop

   
 } // closing analyze() function

Analysis_Template_MC::~Analysis_Template_MC()
{
}


DEFINE_FWK_MODULE(Analysis_Template_MC);
