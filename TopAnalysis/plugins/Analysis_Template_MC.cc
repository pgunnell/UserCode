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
     num_of_Vtx     = fs->make<TH1F>("num_of_Vtx","num_of_Vtx",100,0.,100.);
     num_of_Jets     = fs->make<TH1F>("num_of_Jets","num_of_Jets",100,0.,100.);
     num_of_GenJets     = fs->make<TH1F>("num_of_GenJets","num_of_GenJets",100,0.,100.);

     ptDETJet  = fs->make<TH1F>("ptDETJet","ptDETJet",200,0.,2000.); ptDETJet->Sumw2();
     yDETJet = fs->make<TH1F>("yDETJet","yDETJet",60,-3.,3.); yDETJet->Sumw2();
     phiDETJet = fs->make<TH1F>("phiDETJet","phiDETJet",60, -TMath::Pi(),TMath::Pi()); phiDETJet->Sumw2();

     ChfDETJet = fs->make<TH1F>("ChfDETJet","ChfDETJet",100,0.,1.);
     ElfDETJet = fs->make<TH1F>("ElfDETJet","ElfDETJet",100,0.,1.);
     PhfDETJet = fs->make<TH1F>("PhfDETJet","PhfDETJet",100,0.,1.);
     NhfDETJet = fs->make<TH1F>("NhfDETJet","NhfDETJet",100,0.,1.);
     MufDETJet = fs->make<TH1F>("MufDETJet","MufDETJet",100,0.,1.);

     //Tau1DETJet = fs->make<TH1F>("Tau1DETJet","Tau1DETJet",100,0.,1.);
     //Tau2DETJet = fs->make<TH1F>("Tau2DETJet","Tau2DETJet",100,0.,1.);
     //Tau3DETJet = fs->make<TH1F>("Tau3DETJet","Tau3DETJet",100,0.,1.);
     //Tau21DETJet = fs->make<TH1F>("Tau21DETJet","Tau21DETJet",100,0.,1.);
     //Tau32DETJet = fs->make<TH1F>("Tau32DETJet","Tau32DETJet",100,0.,1.);

     EnergyDETJet= fs->make<TH1F>("EnergyDETJet","EnergyDETJet",100,0.,600.);
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

   int nvtx_=-100; int nJets_=-100; int nGenJets_=-100; 
   std::vector<float>* jetPt_=0; std::vector<float>* jetEta_=0; std::vector<float>* jetPhi_=0; std::vector<float>* jetEnergy_=0; std::vector<float>* jetMass_=0; std::vector<float>* jetMassSoftDrop_=0;

   std::vector<float>* jetchf_=0; std::vector<float>* jetnhf_=0; std::vector<float>* jetphf_=0; std::vector<float>* jetmuf_=0; std::vector<float>* jetelf_=0; std::vector<float>* jettau1_=0; std::vector<float>* jettau2_=0; std::vector<float>* jettau3_=0; 

   mTree->SetBranchAddress("nvtx",&nvtx_);
   mTree->SetBranchAddress("nJets",&nJets_);
   mTree->SetBranchAddress("jetPt",&jetPt_);
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
   
   if(mIsMCarlo){
     mTree->SetBranchAddress("nGenJets",&nGenJets_);
     //mTree->SetBranchAddress("GenJetpt",&GenJetPt_);
     //mTree->SetBranchAddress("GenJeteta",&GenJetEta_);
     //mTree->SetBranchAddress("GenJetphi",&GenJetPhi_);
     //mTree->SetBranchAddress("GenJetmass",&GenJetmass_);
     //mTree->SetBranchAddress("GenJetenergy",&GenJetenergy_);
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
    
    bool triggerfired=false;

    if(mIsMCarlo){
      triggerfired=true;
      hweight=mCrossSection/NEvents;
    }

    if(!mIsMCarlo){
      if(bit_->at(5)) triggerfired=true; //choose your trigger here from the list inside the Ntuples
      double prescale=pre_->at(5);
      hweight=prescale/mIntLumi;
    }

    if(triggerfired==false) continue;

    num_of_Vtx->Fill(nvtx_,hweight);
    num_of_Jets->Fill(nJets_,hweight);

    if(mIsMCarlo){ num_of_GenJets->Fill(nGenJets_,hweight);  }
    
    for(int j=0; j< nJets_; j++){

      ptDETJet->Fill(jetPt_->at(j),hweight);
      yDETJet->Fill(jetEta_->at(j),hweight);
      phiDETJet->Fill(jetPhi_->at(j),hweight);
      EnergyDETJet->Fill(jetEnergy_->at(j),hweight);
      MassDETJet->Fill(jetMass_->at(j),hweight);
      MassSoftDropDETJet->Fill(jetMassSoftDrop_->at(j),hweight);
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
      ChfDETJet->Fill(jetchf_->at(j),hweight);
      NhfDETJet->Fill(jetnhf_->at(j),hweight);
      PhfDETJet->Fill(jetphf_->at(j),hweight);
      MufDETJet->Fill(jetmuf_->at(j),hweight);
      ElfDETJet->Fill(jetelf_->at(j),hweight);

      if(false){
	printf("Number of PFJets=%i\n",nJets_);
	printf("j=%2i  pt=%8.3f  y=%6.3f  phi=%6.3f\n",j,jetPt_->at(j),jetEta_->at(j),jetPhi_->at(j));
      }

    }

    if(mIsMCarlo){
      //Gen jet information
      for(int j=0; j< nGenJets_; j++){ //example histograms
	//ptGENJet->Fill(GenJetPt_->at(j),hweight);
	//yGENJet->Fill(GenJetEta_->at(j),hweight);
	//phiGENJet->Fill(GenJetPhi_->at(j),hweight);
	//EnergyGENJet->Fill(GenJetenergy_->at(j),hweight);
	//MassGENJet->Fill(GenJetmass_->at(j),hweight);
	//MassSoftDropGENJet->Fill(GenJetMassSoftDrop_->at(j),hweight);
	//Tau1GENJet->Fill(GenJettau1_->at(j),hweight);
	//Tau2GENJet->Fill(GenJettau2_->at(j),hweight);
	//Tau3GENJet->Fill(GenJettau3_->at(j),hweight);
      }
    }
    
   } // end of event loop

   
 } // closing analyze() function

Analysis_Template_MC::~Analysis_Template_MC()
{
}


DEFINE_FWK_MODULE(Analysis_Template_MC);
