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
#include "TLorentzVector.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"

#include "Analysis_Template_MCFastSim.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"



using namespace std;

//---------------------------- Constructor Of The Class TriggerTurnOn -------------------------- //
Analysis_Template_MCFastSim::Analysis_Template_MCFastSim(edm::ParameterSet const& cfg)
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
void Analysis_Template_MCFastSim::beginJob()
{
  mInf = TFile::Open(mFileName.c_str());
  mDir = (TDirectoryFile*) mInf->Get(mDirName.c_str());
  mTree = (TTree*) mInf->Get(mTreeName.c_str());

  //------------------ Histogram Booking --------------------------- //

  //DET Jets
  nJets  = fs->make<TH1F>("nJets","nJets",15,0.,15.);
  jetPt  = fs->make<TH1F>("jetPt","jetPt",200,0.,2000.);                   jetPt->Sumw2();
  jetEta = fs->make<TH1F>("jetEta","jetEta",60,-3.,3.);                    jetEta->Sumw2();
  jetPhi = fs->make<TH1F>("jetPhi","jetPhi",60, -TMath::Pi(),TMath::Pi()); jetPhi->Sumw2();

  //GEN JETS
  nGenJets        = fs->make<TH1F>("nGenJets","nGenJets",15,0.,15.);
  GenJetPt        = fs->make<TH1F>("GenJetPt","GenJetPt",200,0.,2000.);                   GenJetPt->Sumw2();
  GenJetPt_total  = fs->make<TH1F>("GenJetPt_total" ,"GenJetPt_total" ,85,0.,1700.);    GenJetPt_total->Sumw2();
  GenJetPt_passed = fs->make<TH1F>("GenJetPt_passed","GenJetPt_passed",85,0.,1700.);    GenJetPt_passed->Sumw2();
  GenJetEta       = fs->make<TH1F>("GenJetEta","GenJetEta",60,-3.,3.);                    GenJetEta->Sumw2();
  GenJetPhi       = fs->make<TH1F>("GenJetPhi","GenJetPhi",60, -TMath::Pi(),TMath::Pi()); GenJetPhi->Sumw2();

  //DET TOPS
  nTopParticles   = fs->make<TH1F>("nTopParticles","nTopParticles",10,0.,10.);
  TopParticlePt   = fs->make<TH1F>("TopParticlePt","TopParticlePt",200,0.,2000.);                  TopParticlePt->Sumw2();
  TopParticleEta  = fs->make<TH1F>("TopParticleEta","TopParticleEta",60,-3.,3);                    TopParticleEta->Sumw2();
  TopParticlePhi  = fs->make<TH1F>("TopParticlePhi","TopParticlePhi",60,-TMath::Pi(),TMath::Pi()); TopParticlePhi->Sumw2();
  TopParticleMass = fs->make<TH1F>("TopParticleMass","TopParticleMass",200,0.,2000.);              TopParticleMass->Sumw2();

  //GEN LQ
  nLeptoQuarkParticles   = fs->make<TH1F>("nLeptoQuarkParticles","nLeptoQuarkParticles",30,0.,30);
  LeptoQuarkParticlePt   = fs->make<TH1F>("LeptoQuarkParticlePt","LeptoQuarkParticlePt",200,0.,2000.);                  LeptoQuarkParticlePt->Sumw2();
  LeptoQuarkParticlePhi  = fs->make<TH1F>("LeptoQuarkParticlePhi","LeptoQuarkParticlePhi",60,-TMath::Pi(),TMath::Pi()); LeptoQuarkParticlePhi->Sumw2();
  LeptoQuarkParticleEta  = fs->make<TH1F>("LeptoQuarkParticleEta","LeptoQuarkParticleEta",60,-3.,3);                    LeptoQuarkParticleEta->Sumw2();
  LeptoQuarkParticleMass = fs->make<TH1F>("LeptoQuarkParticleMass","LeptoQuarkParticleMass",200,0.,2000.);              LeptoQuarkParticleMass->Sumw2();

  //DET Muons
  nMuons  = fs->make<TH1F>("nMuons","nMuons",30,0.,30);
  muonPt  = fs->make<TH1F>("muonPt","muonPt",200,0.,2000.);                  muonPt->Sumw2();
  muonPhi = fs->make<TH1F>("muonPhi","muonPhi",60,-TMath::Pi(),TMath::Pi()); muonPhi->Sumw2();
  muonEta = fs->make<TH1F>("muonEta","muonEta",60,-3.,3);                    muonEta->Sumw2();

  //GEN MUONS
  nMuonParticles  = fs->make<TH1F>("nMuonParticles","nMuonParticles",30,0.,30);
  MuonParticlePt  = fs->make<TH1F>("MuonParticlePt","MuonParticlePt",200,0.,2000.);                  MuonParticlePt->Sumw2();
  MuonParticlePhi = fs->make<TH1F>("MuonParticlePhi","MuonParticlePhi",60,-TMath::Pi(),TMath::Pi()); MuonParticlePhi->Sumw2();
  MuonParticleEta = fs->make<TH1F>("MuonParticleEta","MuonParticleEta",60,-3.,3);                    MuonParticleEta->Sumw2();

  MLQreco_detMu = fs->make<TH1F>("MLQreco_detMu","MLQreco_detMu",200,0.,2000.); MLQreco_detMu->Sumw2();
  MLQreco_genMu = fs->make<TH1F>("MLQreco_genMu","MLQreco_genMu",200,0.,2000.); MLQreco_genMu->Sumw2();
  ST            = fs->make<TH1F>("ST","ST",200,1000.,5000.);                    ST->Sumw2();


  WbosonMass = fs->make<TH1F>("WbosonMass","WbosonMass",200,0.,200.); WbosonMass->Sumw2();
  MLQfullreco_detMu = fs->make<TH1F>("MLQfullreco_detMu","MLQfullreco_detMu",200,0.,4000.); MLQfullreco_detMu->Sumw2();
  
  LeptTopfullreco = fs->make<TH1F>("LeptTopfullreco","LeptTopfullreco",100,0.,1000.); LeptTopfullreco->Sumw2();
  HadrTopfullreco = fs->make<TH1F>("HadrTopfullreco","HadrTopfullreco",100,0.,1000.); HadrTopfullreco->Sumw2();
  chi2        = fs->make<TH1F>("chi2","chi2",150,0.,1500.); chi2->Sumw2(); 
  LQHadrMassreco = fs->make<TH1F>("LQHadrMassreco","LQHadrMassreco",100,0.,4000.); LQHadrMassreco->Sumw2();
  LQLeptMassreco = fs->make<TH1F>("LQLeptMassreco","LQLeptMassreco",100,0.,4000.); LQLeptMassreco->Sumw2();
}


//------------------------ endjob() function declaration ---------------------- //
void Analysis_Template_MCFastSim::endJob()
{
  mInf->Close();
}

//--------------------------- analyze() fuction declaration ------------------ //
void Analysis_Template_MCFastSim::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup)
{
  unsigned NEntries = mTree->GetEntries();
  cout << "Reading TREE: " << NEntries << " events" <<endl;

  int decade = 0 ;

  float hweight = 1.;  //Initial value set to one

  //GEN JETS
  std::vector<float>* GenJetPt_ = 0;
  std::vector<float>* GenJetEta_ = 0;
  std::vector<float>* GenJetPhi_ = 0;
  mTree->SetBranchAddress("GenJetPt",&GenJetPt_);
  mTree->SetBranchAddress("GenJetEta",&GenJetEta_);
  mTree->SetBranchAddress("GenJetPhi",&GenJetPhi_);
  //DET JETS
  std::vector<float>* jetPt_   = 0;
  std::vector<float>* jetEta_  = 0;
  std::vector<float>* jetPhi_  = 0;
  std::vector<float>* jetMass_ = 0;
  std::vector<int>* jetHadronFlavour_ = 0;
  mTree->SetBranchAddress("jetPt",  &jetPt_);
  mTree->SetBranchAddress("jetEta", &jetEta_);
  mTree->SetBranchAddress("jetPhi", &jetPhi_);
  mTree->SetBranchAddress("jetMass",&jetMass_);
  mTree->SetBranchAddress("jetFlavorHadron",&jetHadronFlavour_);

  //DET TOPS
  std::vector<float>* TopParticlePt_   = 0;
  std::vector<float>* TopParticleEta_  = 0;
  std::vector<float>* TopParticlePhi_  = 0;
  std::vector<float>* TopParticleId_   = 0;
  std::vector<float>* TopParticleMass_ = 0;
  mTree->SetBranchAddress("TopParticlePt",  &TopParticlePt_);
  mTree->SetBranchAddress("TopParticleEta", &TopParticleEta_);
  mTree->SetBranchAddress("TopParticlePhi", &TopParticlePhi_);
  mTree->SetBranchAddress("TopParticleId",  &TopParticleId_);
  mTree->SetBranchAddress("TopParticleMass",&TopParticleMass_);

  //GEN LQ
  std::vector<int>*   LeptoQuarkParticleId_   = 0;
  std::vector<float>* LeptoQuarkParticlePt_   = 0;
  std::vector<float>* LeptoQuarkParticlePhi_  = 0;
  std::vector<float>* LeptoQuarkParticleEta_  = 0;
  std::vector<float>* LeptoQuarkParticleMass_ = 0;
  mTree->SetBranchAddress("LeptoQuarkParticleId",  &LeptoQuarkParticleId_);
  mTree->SetBranchAddress("LeptoQuarkParticlePt",  &LeptoQuarkParticlePt_);
  mTree->SetBranchAddress("LeptoQuarkParticlePhi", &LeptoQuarkParticlePhi_);
  mTree->SetBranchAddress("LeptoQuarkParticleEta", &LeptoQuarkParticleEta_);
  mTree->SetBranchAddress("LeptoQuarkParticleMass",&LeptoQuarkParticleMass_);

  //GEN MUONS
  std::vector<int>*   MuonParticleId_     = 0;
  std::vector<int>*   MuonParticleStatus_ = 0;
  std::vector<float>* MuonParticlePt_     = 0;
  std::vector<float>* MuonParticlePhi_    = 0;
  std::vector<float>* MuonParticleEta_    = 0;
  std::vector<float>* MuonParticleMass_   = 0;
  mTree->SetBranchAddress("MuonParticleId",    &MuonParticleId_);
  mTree->SetBranchAddress("MuonParticleStatus",&MuonParticleStatus_);
  mTree->SetBranchAddress("MuonParticlePt",    &MuonParticlePt_);
  mTree->SetBranchAddress("MuonParticlePhi",   &MuonParticlePhi_);
  mTree->SetBranchAddress("MuonParticleEta",   &MuonParticleEta_);
  mTree->SetBranchAddress("MuonParticleMass",  &MuonParticleMass_);

  //DET MUONS
  std::vector<float>* muonPt_     = 0;
  std::vector<float>* muonEta_    = 0;
  std::vector<float>* muonPhi_    = 0;
  std::vector<int>*   muonCharge_ = 0;
  mTree->SetBranchAddress("muonPt",    &muonPt_);
  mTree->SetBranchAddress("muonEta",   &muonEta_);
  mTree->SetBranchAddress("muonPhi",   &muonPhi_);
  mTree->SetBranchAddress("muonCharge",&muonCharge_);
  const double muonMass = 0.1056583715;

  //DET ELECTRONS
  std::vector<float>* elecPt_  = 0;
  std::vector<float>* elecEta_ = 0;
  std::vector<float>* elecPhi_ = 0;
  std::vector<int>* elecCharge_ = 0;
  mTree->SetBranchAddress("elecPt", &elecPt_);
  mTree->SetBranchAddress("elecEta",&elecEta_);
  mTree->SetBranchAddress("elecPhi",&elecPhi_);
  mTree->SetBranchAddress("elecCharge",&elecCharge_);

  //MET
  float MET;
  float METeta;
  float METphi;
  mTree->SetBranchAddress("MET",   &MET);
  mTree->SetBranchAddress("METphi",   &METphi);
  mTree->SetBranchAddress("METeta",   &METeta);

  //DET SCALAR HT
  float ScalarHT    = 0;
  mTree->SetBranchAddress("ScalarHT",   &ScalarHT);

  //counter for events
  int numberOfEvents = 0;

  for(unsigned  l=0; l<NEntries; l++)
  {
    //----------- progress report -------------
    double progress = 10.0*l/(1.0*NEntries);
    int k = TMath::FloorNint(progress);
    if(k > decade)
    {
      cout << 10*k << " %" << endl;
    }
    decade = k;
    //cout<<"new event"<<endl;

    //READ THE EVENT
    mTree->GetEntry(l);
    if(mIsMCarlo)
    {
      hweight=mCrossSection/NEntries;
    }

    //NUMBER OF PARTICLES
    const int nMuons_               = muonPt_              ->size();
    const int nMuonParticles_       = MuonParticlePt_      ->size();
    const int nElecs_               = elecPt_              ->size();
    const int nJets_                = jetPt_               ->size();
    const int nGenJets_             = GenJetPt_            ->size();
    const int nTopParticles_        = TopParticlePt_       ->size();
    const int nLeptoQuarkParticles_ = LeptoQuarkParticlePt_->size();

    //SETUP FOR EVENT SELECTION
    //DET JETS
    int counterJet = 0;
    std::vector<float> selJetEta;
    std::vector<float> selJetPhi;
    for(int i=0; i<nJets_; i++)
    {
      if(jetPt_->at(i) > 30 && abs(jetEta_->at(i)) < 2.4)  //pT, eta selection for jets
      {
        counterJet++;
        selJetEta.push_back(jetEta_->at(i));
        selJetPhi.push_back(jetPhi_->at(i));
      }
    }

    //GEN JETS
    int counterGenJet = 0;
    std::vector<float> selGenJetEta;
    std::vector<float> selGenJetPhi;
    std::vector<float> selGenJetPt;
    for(int i=0; i<nGenJets_; i++)
    {
      if(GenJetPt_->at(i) > 30 && abs(GenJetEta_->at(i)) < 2.4)  //pT, eta selection for jets
      {
        counterGenJet++;
        selGenJetEta.push_back(GenJetEta_->at(i));
        selGenJetPhi.push_back(GenJetPhi_->at(i));
        selGenJetPt .push_back(GenJetPt_ ->at(i));
      }
    }


    //DET MUONS
    double STmu = 0.0;
    int counterSelMu      = 0;
    int counterSelMuPlus  = 0;
    int counterSelMuMinus = 0;
    TLorentzVector leadingMuMinus;
    TLorentzVector leadingMuPlus;
    
    for(int j=0; j<nMuons_; j++)
    {
      if(muonPt_->at(j) > 30. && fabs(muonEta_->at(j)) < 2.4)   //pT, eta selection for muons
      {
        counterSelMu++;
        STmu += muonPt_->at(j);
        if(muonCharge_->at(j) == -1)
        {
          counterSelMuMinus++;
          if(counterSelMuMinus == 1)
          {
            leadingMuMinus.SetPtEtaPhiM(muonPt_->at(j),muonEta_->at(j),muonPhi_->at(j),muonMass);
          }
        }
        else if(muonCharge_->at(j) == +1)
        {
          counterSelMuPlus++;
          if(counterSelMuPlus == 1)
          {
            leadingMuPlus.SetPtEtaPhiM(muonPt_->at(j),muonEta_->at(j),muonPhi_->at(j),muonMass);
          }
        }
      }
    }


    //GEN MUONS
    double STmuPart = 0.0;
    int counterSelMuPart      = 0;
    int counterSelMuPartPlus  = 0;
    int counterSelMuPartMinus = 0;
    TLorentzVector leadingMuPartMinus;
    TLorentzVector leadingMuPartPlus;
    for(int j=0; j<nMuonParticles_; j++)
    {
      if(MuonParticlePt_->at(j) > 30. && abs(MuonParticleEta_->at(j)) < 2.4  && MuonParticleStatus_->at(j) == 23)   //pT, eta selection for muons
      {
        counterSelMuPart++;
        STmuPart += MuonParticlePt_->at(j);
        if(MuonParticleId_->at(j) == 13)
        {
          counterSelMuPartMinus++;
          if(counterSelMuPartMinus == 1)
          {
            leadingMuPartMinus.SetPtEtaPhiM(MuonParticlePt_->at(j),MuonParticleEta_->at(j),MuonParticlePhi_->at(j),MuonParticleMass_->at(j));
          }
        }
        else if(MuonParticleId_->at(j) == -13)
        {
          counterSelMuPartPlus++;
          if(counterSelMuPartPlus == 1)
          {
            leadingMuPartPlus.SetPtEtaPhiM(MuonParticlePt_->at(j),MuonParticleEta_->at(j),MuonParticlePhi_->at(j),MuonParticleMass_->at(j));
          }
        }
      }
    }


    //ELECTRONS
    double STel = 0.0;
    int counterSelEl = 0;
    for(int j=0; j<nElecs_; j++)
    {
      if(elecPt_->at(j) > 30. && abs(elecEta_->at(j)) < 2.4)          //pT, eta selection for electrons
      {
        counterSelEl++;
        STel += elecPt_->at(j);
      }
    }

    //LEPTOQUARK MASS RECONSTUCTION FROM TOP AND DET MUONS: always pair production -> always exactly 2 tops (t + tbar)
    TLorentzVector top;
    TLorentzVector topBar;
    double recoMLQ_top    = 0.0;
    double recoMLQ_topBar = 0.0;
    //double recoMLQ        = 0.0;
    for(int i=0; i<nTopParticles_; i++)
    {
      if(TopParticleId_->at(i) == 6)         //t + mu^minus
      {
        top.SetPtEtaPhiM(TopParticlePt_->at(i),TopParticleEta_->at(i),TopParticlePhi_->at(i),TopParticleMass_->at(i));
        recoMLQ_top = (top + leadingMuMinus).M();

      }
      else if(TopParticleId_->at(i) == -6)   //t^bar + mu^plus
      {
        topBar.SetPtEtaPhiM(TopParticlePt_->at(i),TopParticleEta_->at(i),TopParticlePhi_->at(i),TopParticleMass_->at(i));
        recoMLQ_topBar = (topBar + leadingMuPlus).M();
      }
    }
    //recoMLQ = (recoMLQ_top + recoMLQ_topBar) / 2.;  //averaged


    //GEN MUONS RECONSTRUCTION
    double recoMLQ_top_GenMu    = 0.0;
    double recoMLQ_topBar_GenMu = 0.0;
    //double recoMLQ_GenMu        = 0.0;
    for(int i=0; i<nTopParticles_; i++)
    {
      if(TopParticleId_->at(i) == 6)         //t + mu^minus
      {
        top.SetPtEtaPhiM(TopParticlePt_->at(i),TopParticleEta_->at(i),TopParticlePhi_->at(i),TopParticleMass_->at(i));
        recoMLQ_top_GenMu = (top + leadingMuPartMinus).M();
      }
      else if(TopParticleId_->at(i) == -6)   //t^bar + mu^plus
      {
        topBar.SetPtEtaPhiM(TopParticlePt_->at(i),TopParticleEta_->at(i),TopParticlePhi_->at(i),TopParticleMass_->at(i));
        recoMLQ_topBar_GenMu = (topBar + leadingMuPartPlus).M();
      }
    }
    //recoMLQ_GenMu = (recoMLQ_top_GenMu + recoMLQ_topBar_GenMu) / 2.;  //averaged


    bool selectionPass=false;

    //FULL SELECTION: DETECTOR LEVEL
    if(counterJet >= 2)    // >= 2 jets
    {
      if(counterSelMuPlus >= 1 && counterSelMuMinus >= 1)     // >=2 selected muons
      {
        if(STmu+STel > 200.)     // STlep > 200 GeV
        {
          if(ScalarHT > 350.)     // ST > 350 GeV
          {
	    
	    selectionPass=true;
            numberOfEvents++;
            //FILL DET MU HISTOGRAMS
            nMuons    ->Fill(counterSelMu,hweight);
            for(int j=0; j<nMuons_; j++)
            {
              if(muonPt_->at(j) > 30. && abs(muonEta_->at(j)) < 2.4)
              {
                muonPt ->Fill(muonPt_ ->at(j),hweight);
                muonPhi->Fill(muonPhi_->at(j),hweight);
                muonEta->Fill(muonEta_->at(j),hweight);
              }
            }
            //FILL JET HISTOGRAMS
            nJets->Fill(nJets_,hweight);
            for(int i=0; i<nJets_; i++)
            {
              if(jetPt_->at(i) > 30 && abs(jetEta_->at(i)) < 2.4)
              {
                jetPt-> Fill(jetPt_ ->at(i),hweight);
                jetEta->Fill(jetEta_->at(i),hweight);
                jetPhi->Fill(jetPhi_->at(i),hweight);
              }
            }
            //FILL TOP HISTOGRAMS
            nTopParticles->Fill(nTopParticles_, hweight);
            for(int i=0; i<nTopParticles_; i++)
            {
              TopParticlePt  ->Fill(TopParticlePt_->at(i)  ,hweight);
              TopParticlePhi ->Fill(TopParticlePhi_->at(i) ,hweight);
              TopParticleEta ->Fill(TopParticleEta_->at(i) ,hweight);
              TopParticleMass->Fill(TopParticleMass_->at(i),hweight);
            }

            //FILL LQ HISTOGRAMS
            nLeptoQuarkParticles->Fill(nLeptoQuarkParticles_,hweight);
            for(int i=0; i<nLeptoQuarkParticles_; i++)
            {
              LeptoQuarkParticlePt->  Fill(LeptoQuarkParticlePt_  ->at(i),hweight);
              LeptoQuarkParticlePhi-> Fill(LeptoQuarkParticlePhi_ ->at(i),hweight);
              LeptoQuarkParticleEta-> Fill(LeptoQuarkParticleEta_ ->at(i),hweight);
              LeptoQuarkParticleMass->Fill(LeptoQuarkParticleMass_->at(i),hweight);
            }

            //LQ RECO MASS & S_T
            MLQreco_detMu->Fill(recoMLQ_top,hweight);
            MLQreco_detMu->Fill(recoMLQ_topBar,hweight);
            ST           ->Fill(ScalarHT,hweight);
          }
        }
      }
    }

    //FULL SELECTION: GENERATOR LEVEL
    if(counterGenJet >= 2)    // >= 2 jets
    {
      if(counterSelMuPartPlus >= 1 && counterSelMuPartMinus >= 1)     // >=2 selected muons
      {
        if(STmuPart+STel > 200.)     // STlep > 200 GeV
        {
          if(ScalarHT > 350.)     // ST > 350 GeV
          {
            //FILL GEN MU HISTOGRAMS
            nMuonParticles->Fill(counterSelMuPart,hweight);
            for(int j=0; j<nMuonParticles_; j++)
            {
              if(MuonParticlePt_->at(j) > 30. && abs(MuonParticleEta_->at(j)) < 2.4)
              {
                MuonParticlePt ->Fill(MuonParticlePt_ ->at(j),hweight);
                MuonParticlePhi->Fill(MuonParticlePhi_->at(j),hweight);
                MuonParticleEta->Fill(MuonParticleEta_->at(j),hweight);
              }
            }

            if(counterSelMuPartMinus >= 1 && counterSelMuPartPlus >= 1)   //&& counterSelMu+counterSelEl >= 3
            {
              MLQreco_genMu->Fill(recoMLQ_top_GenMu   ,hweight);
              MLQreco_genMu->Fill(recoMLQ_topBar_GenMu,hweight);
            }
          }
        }
      }
    }

    //JET MACHTING EFFICIENCY
    nGenJets->Fill(nGenJets_,hweight);
    bool passMatching;
    for(int i=0; i<nGenJets_; i++)
    {
      GenJetPt->Fill(GenJetPt_->at(i),hweight);
      GenJetPt_total ->Fill(GenJetPt_ ->at(i));         //no weight for efficiency calculation
      GenJetEta->Fill(GenJetEta_->at(i),hweight);
      GenJetPhi->Fill(GenJetPhi_->at(i),hweight);
      passMatching = false;
      for(int k=0; k<nJets_; k++)
      {
        if(DeltaRFS(GenJetEta_->at(i),GenJetPhi_->at(i),jetEta_->at(k),jetPhi_->at(k)) < 0.2)
        {
          passMatching = true;
        }
      }
      if(passMatching == true)
      {
        GenJetPt_passed->Fill(GenJetPt_->at(i)); //no weight for efficiency calculation
      }
    }

    vector <TLorentzVector> Jet_v4;
    
    for(int j=0; j<nJets_; j++){
      if(j==7) break;
      if(jetPt_->at(j)>30 && fabs(jetEta_->at(j))<2.4){
	TLorentzVector jet_v4_temp;
	jet_v4_temp.SetPtEtaPhiM(jetPt_->at(j),jetEta_->at(j),jetPhi_->at(j),jetMass_->at(j));
	Jet_v4.push_back(jet_v4_temp);
      }
    }

    vector <TLorentzVector> Jet_v4_bflavour;
    vector <TLorentzVector> Jet_v4_otherflavour;
    
    for(int j=0; j<nJets_; j++){
      if(j==7) break;
      if(jetPt_->at(j)>30 && fabs(jetEta_->at(j))<2.4 && fabs(jetHadronFlavour_->at(j))==5){
	TLorentzVector jet_v4_temp;
	jet_v4_temp.SetPtEtaPhiM(jetPt_->at(j),jetEta_->at(j),jetPhi_->at(j),jetMass_->at(j));
	Jet_v4_bflavour.push_back(jet_v4_temp);
      }
      if(jetPt_->at(j)>30 && fabs(jetEta_->at(j))<2.4 && fabs(jetHadronFlavour_->at(j))!=5){
	TLorentzVector jet_v4_temp;
	jet_v4_temp.SetPtEtaPhiM(jetPt_->at(j),jetEta_->at(j),jetPhi_->at(j),jetMass_->at(j));
	Jet_v4_otherflavour.push_back(jet_v4_temp);
      }
    }

    if(Jet_v4.size()<4) continue;

    if(!selectionPass) continue;

    vector <TLorentzVector> Tops_Lept_Reco_v4;
    vector <TLorentzVector> Tops_Hadr_Reco_v4;

    vector <TLorentzVector> LQRecos_Hadr_v4;
    vector <TLorentzVector> LQRecos_Lept_v4; 

    //calculate the best combination for having the two reco-top (two hadronic)  
    if(counterSelMu==2 && counterSelEl==0 && selectionPass){
      
      //perhaps here one can use the b-hadron flavour

      //cout<<"size "<<Jet_v4_bflavour.size()<<" "<<Jet_v4_otherflavour.size()<<endl;

      for (unsigned int first=0; first<Jet_v4_bflavour.size(); first++){
	for (unsigned int second=0; second<Jet_v4_otherflavour.size(); second++){
	  for (unsigned int third=second+1; third<Jet_v4_otherflavour.size(); third++){
	    Tops_Hadr_Reco_v4.push_back(Jet_v4_bflavour.at(first)+Jet_v4_otherflavour.at(second)+Jet_v4_otherflavour.at(third));//for three jets
	    Tops_Lept_Reco_v4.push_back(Jet_v4_bflavour.at(first)+Jet_v4_otherflavour.at(second)+Jet_v4_otherflavour.at(third));//for three jets
	    
	    for (unsigned int fourth=second+1; fourth<Jet_v4_otherflavour.size(); fourth++){
	      Tops_Hadr_Reco_v4.push_back(Jet_v4_bflavour.at(first)+Jet_v4_otherflavour.at(second)+Jet_v4_otherflavour.at(third)+Jet_v4_otherflavour.at(fourth));//for fourth jets
	      Tops_Lept_Reco_v4.push_back(Jet_v4_bflavour.at(first)+Jet_v4_otherflavour.at(second)+Jet_v4_otherflavour.at(third)+Jet_v4_otherflavour.at(fourth));//for fourth jets
	    }
	  }
	}
      
      }
      
      //we have now all possible tops and we can calculate the reconstructed LQ mass
      //the constrain is that the muon and the electron needs to be oppositely charged
      
      //LQ from leptonic top
      for (unsigned int j=0; j<Tops_Lept_Reco_v4.size(); j++){
	TLorentzVector muon_to_associate_v4_temp_1;
	muon_to_associate_v4_temp_1.SetPtEtaPhiM(muonPt_->at(0),muonEta_->at(0),muonPhi_->at(0),0.105658);
	TLorentzVector LQRecos_v4_temp_1 = Tops_Lept_Reco_v4.at(j)+muon_to_associate_v4_temp_1;
	LQRecos_Lept_v4.push_back(LQRecos_v4_temp_1);

	TLorentzVector muon_to_associate_v4_temp_2;
	muon_to_associate_v4_temp_2.SetPtEtaPhiM(muonPt_->at(1),muonEta_->at(1),muonPhi_->at(1),0.105658);
	TLorentzVector LQRecos_v4_temp_2 = Tops_Lept_Reco_v4.at(j)+muon_to_associate_v4_temp_2;
	LQRecos_Lept_v4.push_back(LQRecos_v4_temp_2);

	TLorentzVector LQRecos_v4_temp_3 = Tops_Hadr_Reco_v4.at(j)+muon_to_associate_v4_temp_1;
	LQRecos_Hadr_v4.push_back(LQRecos_v4_temp_1);

	TLorentzVector LQRecos_v4_temp_4 = Tops_Hadr_Reco_v4.at(j)+muon_to_associate_v4_temp_2;
	LQRecos_Hadr_v4.push_back(LQRecos_v4_temp_2);
      
      }
    
    }

    //END OF HADRONIC PART

    bool electronChannel=false;
    vector <TLorentzVector> wboson_v4;
    vector <TLorentzVector> neutrino_v4;
    TLorentzVector electron_v4;
  
    if(counterSelMu==2 && counterSelEl>=1 && selectionPass){

      //if(!selectionPass || nElecs_<1) continue;
      electronChannel=true;

      //Top mass reconstruction
      //simple case: electron in the final state
      
      //select an electron and you have the four momentum
    
      if(elecPt_->at(0)<30 || fabs(elecEta_->at(0))>2.4) continue;
      electron_v4.SetPtEtaPhiM(elecPt_->at(0),elecEta_->at(0),elecPhi_->at(0),0.000511);
      
      //define the neutrino from the MET
      float metPx = MET * cos(METphi);
      float metPy = MET * sin(METphi); 
      
      const float mass_w = 80.399;
      float mu = mass_w * mass_w / 2 + electron_v4.Px() * metPx + electron_v4.Py() * metPy;//scalar product between lepton and neutrino
      float A = - (electron_v4.E() * electron_v4.E());
      float B = mu * electron_v4.Pz();
      float C = mu * mu - electron_v4.E() * electron_v4.E() * (MET * MET);
      float discriminant = B * B - A * C;
      
      if (0 >= discriminant) {
	// Take only real part of the solution for pz:
	TLorentzVector solution;
	solution.SetPxPyPzE(metPx,metPy,-B / A,0);
	solution.SetE(solution.P());
	//neutrino_v4.push_back(solution);
	solution.SetPxPyPzE(metPx,metPy,-B / A,solution.E());
	neutrino_v4.push_back(solution);
      }
      else {
	discriminant = sqrt(discriminant);
	TLorentzVector solution;
	solution.SetPxPyPzE(metPx,metPy,(-B - discriminant) / A,0);
	solution.SetE(solution.P());
	solution.SetPxPyPzE(metPx,metPy,(-B - discriminant) / A,solution.E());
	neutrino_v4.push_back(solution);
	
	TLorentzVector solution2;
	solution2.SetPxPyPzE(metPx,metPy,(-B + discriminant) / A,0);
	solution2.SetE(solution2.P());
	solution2.SetPxPyPzE(metPx,metPy,(-B + discriminant) / A,solution2.E());
	neutrino_v4.push_back(solution2);

      }

      //neutrino_v4.SetPtEtaPhiM(MET,METeta,METphi,0.);
      
      //define the W system
      for(unsigned int j=0; j<neutrino_v4.size(); j++){
	wboson_v4.push_back(electron_v4+neutrino_v4.at(j));
      }
    
      //cout<<wboson_v4.M()<<endl;

      //now you need to associate the W system to the different jets in the system
      //define the four vectors of all jets
  
      //calculate the best combination for having the two reco-top (one leptonic and one hadronic)  
      vector <TLorentzVector> Tops_Lept_Reco_v4;
      vector <TLorentzVector> Tops_Hadr_Reco_v4;
      
      
      //for(unsigned int i=0; i<wboson_v4.size(); i++){
      //cout<<"wmass "<<i<<" "<<wboson_v4.at(i).M()<<endl;
      //}

      for (unsigned int first=0; first<Jet_v4.size(); first++){
	//Leptonic top
	for(unsigned int i=0; i<wboson_v4.size(); i++){
	  
	  TLorentzVector top_v4_temp=wboson_v4.at(i)+Jet_v4.at(first);//for single jets
	  Tops_Lept_Reco_v4.push_back(top_v4_temp);
	}

	for (unsigned int second=first+1; second<Jet_v4.size(); second++){
	  for(unsigned int i=0; i<wboson_v4.size(); i++){
	    Tops_Lept_Reco_v4.push_back(wboson_v4.at(i)+Jet_v4.at(first)+Jet_v4.at(second));//for dijets
	  }
	  //Hadronic top
	  for (unsigned int third=second+1; third<Jet_v4.size(); third++){
	    Tops_Hadr_Reco_v4.push_back(Jet_v4.at(first)+Jet_v4.at(second)+Jet_v4.at(third));//for three jets
	    
	    for (unsigned int fourth=third+1; fourth<Jet_v4.size(); fourth++){
	      Tops_Hadr_Reco_v4.push_back(Jet_v4.at(first)+Jet_v4.at(second)+Jet_v4.at(third)+Jet_v4.at(fourth));//for fourth jets
	    }
	  }
	}
	
      }

      //we have now all possible tops and we can calculate the reconstructed LQ mass
      //the constrain is that the muon and the electron needs to be oppositely charged
      vector <TLorentzVector> LQRecos_Hadr_v4;
      vector <TLorentzVector> LQRecos_Lept_v4;
      
      for (int i = 0; i<2; i++){
	//LQ from leptonic top
	if(muonCharge_->at(i)!=elecCharge_->at(0)){     
	  for (unsigned int j=0; j<Tops_Lept_Reco_v4.size(); j++){
	    TLorentzVector muon_to_associate_v4_temp;
	    muon_to_associate_v4_temp.SetPtEtaPhiM(muonPt_->at(i),muonEta_->at(i),muonPhi_->at(i),0.105658);
	    TLorentzVector LQRecos_v4_temp = Tops_Lept_Reco_v4.at(j)+muon_to_associate_v4_temp;
	    LQRecos_Lept_v4.push_back(LQRecos_v4_temp);
	  }
	}
	
	//LQ from hadronic top
	else{
	  for (unsigned int j=0; j<Tops_Hadr_Reco_v4.size(); j++){
	    TLorentzVector muon_to_associate_v4_temp;
	    muon_to_associate_v4_temp.SetPtEtaPhiM(muonPt_->at(i),muonEta_->at(i),muonPhi_->at(i),0.105658);
	    TLorentzVector LQRecos_v4_temp = Tops_Hadr_Reco_v4.at(j)+muon_to_associate_v4_temp;
	    LQRecos_Hadr_v4.push_back(LQRecos_v4_temp);
	  }
	}
      }
    }

    //END OF LEPTONIC PART
  
    //CASE OF THREE MUONS (muonic decay of one of the tops)     

    //defintion of the chi2
    const double mass_thad = 174.;
    const double mass_thad_sigma = 16.;
    const double mass_tlep = 173.;
    const double mass_tlep_sigma = 22.;
    const double mass_LQ_diff_rel = -0.013;
    const double mass_LQ_diff_rel_sigma = 0.088;
    const double pt_LQ_diff = 0.65;
    const double pt_LQ_diff_sigma = 46.;
    const double pt_ratio = 1;
    const double pt_ratio_sigma = 0.15;
    const double dphi_LQ = 3.14;
    const double dphi_LQ_sigma = 0.09;
    double chi2_total=100000.;
    const double mass_LQ_mean_rec=1200.;

    double LQRecos_Hadr_rec_mass=-10;
    double LQRecos_Lept_rec_mass=-10;

    double Top_Hadr_rec_mass=-10;
    double Top_Lept_rec_mass=-10;

    double LQHadrMass=-10;
    double LQLeptMass=-10;
  
    for (unsigned int j=0; j<Tops_Lept_Reco_v4.size(); j++){
      double chi2_top_lept = pow((Tops_Lept_Reco_v4.at(j).M() - mass_tlep) / mass_tlep_sigma, 2);

      if(fabs(Tops_Lept_Reco_v4.at(j).M() - mass_tlep) > 1000) continue;

      for (unsigned int i=0; i<Tops_Hadr_Reco_v4.size(); i++){
	double chi2_thad = pow((Tops_Hadr_Reco_v4.at(i).M() - mass_thad) / mass_thad_sigma, 2);
	if(fabs(Tops_Hadr_Reco_v4.at(i).M() - mass_thad) > 1000) continue;

	for (unsigned int t=0; t<LQRecos_Hadr_v4.size(); t++){
	  for (unsigned int g=0; g<LQRecos_Lept_v4.size(); g++){
	    double chi2_MLQdiff_rel = pow(((LQRecos_Hadr_v4.at(t).M() - LQRecos_Lept_v4.at(g).M())/mass_LQ_mean_rec - mass_LQ_diff_rel) / mass_LQ_diff_rel_sigma, 2);

	    double chi2_total_temp = chi2_thad+chi2_MLQdiff_rel+chi2_top_lept;
	    if (chi2_total_temp<chi2_total){
	      chi2_total=chi2_total_temp;
	      LQHadrMass=LQRecos_Hadr_v4.at(t).M();
	      LQLeptMass=LQRecos_Lept_v4.at(g).M();
	      LQRecos_Hadr_rec_mass=LQRecos_Hadr_v4.at(t).M();
	      LQRecos_Lept_rec_mass=LQRecos_Lept_v4.at(g).M();
	      Top_Hadr_rec_mass=Tops_Hadr_Reco_v4.at(i).M();
	      Top_Lept_rec_mass=Tops_Lept_Reco_v4.at(j).M();
	    }
	  }
	}
      }
    }

    //plot top mass and chi2
  
    double AverageMass= (LQRecos_Hadr_rec_mass+ LQRecos_Lept_rec_mass)/2;
    MLQfullreco_detMu->Fill(AverageMass,hweight);
    LQHadrMassreco->Fill(LQHadrMass,hweight);
    LQLeptMassreco->Fill(LQLeptMass,hweight);
    
    if(electronChannel){
      double Wmass=wboson_v4.at(0).M();
      if(wboson_v4.size()>1){
	if(fabs(Wmass-80.399)>fabs(wboson_v4.at(1).M()-80.399)){
	  Wmass=wboson_v4.at(1).M();
	}
      }

      WbosonMass->Fill(Wmass,hweight);
    }

    HadrTopfullreco->Fill(Top_Hadr_rec_mass,hweight);
    LeptTopfullreco->Fill(Top_Lept_rec_mass,hweight);

    chi2->Fill(chi2_total,hweight);

    //clear all vectors
    LQRecos_Hadr_v4.clear();
    LQRecos_Lept_v4.clear();
    Tops_Hadr_Reco_v4.clear();
    Tops_Lept_Reco_v4.clear();
    Jet_v4.clear();
    neutrino_v4.clear();
    wboson_v4.clear();

    //cout << endl;
  } // end of event loop
  cout << "numberOfEventsAfterSelection: " << numberOfEvents << endl;



  //JET MATCHING EFFICIENCY PLOT
  //If the denominator becomes 0 or pass > total, the corresponding bin is skipped -> change plot range or binning!
  TGraphAsymmErrors jetMatchEff = TGraphAsymmErrors(GenJetPt_passed,GenJetPt_total);
  cout << "jetMatchEff successfully computed" << endl;
  jetMatchEff.Write("jetMatchEff");

  



}  // closing analyze() function

Analysis_Template_MCFastSim::~Analysis_Template_MCFastSim()
{
}

DEFINE_FWK_MODULE(Analysis_Template_MCFastSim);
