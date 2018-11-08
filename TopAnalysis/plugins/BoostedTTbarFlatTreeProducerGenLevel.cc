#include <string>
#include <cmath>
#include <functional>
#include <vector>
#include <cassert>
#include "Math/SpecFuncMathMore.h"
#include "TMath.h"
#include "TVectorD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"

#include "UserCode/TopAnalysis/plugins/BoostedTTbarFlatTreeProducerGenLevel.h"

using namespace std;
using namespace reco;
using namespace fastjet;

BoostedTTbarFlatTreeProducerGenLevel::BoostedTTbarFlatTreeProducerGenLevel(edm::ParameterSet const& cfg)
{ 
  genjetsToken          = consumes<GenJetCollection>(cfg.getUntrackedParameter<edm::InputTag>("genjets",edm::InputTag("")));
  genParticlesToken     = consumes<edm::View<reco::GenParticle> >(edm::InputTag("genParticles"));
  GenptMin_             = cfg.getUntrackedParameter<double>("GenptMin");
  GenetaMax_            = cfg.getUntrackedParameter<double>("GenetaMax");
  genEvtInfoToken = consumes<GenEventInfoProduct>(edm::InputTag("generator"));
  jetFlavourInfosToken_ = consumes<reco::JetFlavourInfoMatchingCollection>( cfg.getParameter<edm::InputTag>("jetFlavourInfos"));
  isHigherOrder_ = cfg.getUntrackedParameter<bool>("isHigherOrder",false);
  pupInfoToken          = consumes<edm::View<PileupSummaryInfo> >(edm::InputTag("slimmedAddPileupInfo"));
  lheEvtInfoToken = consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer"));

  //Gen Jet information
  fAKJetDef = new fastjet::JetDefinition(fastjet::antikt_algorithm, 0.8);
  int activeAreaRepeats = 1;
  double ghostArea      = 0.01;
  double ghostEtaMax    = 7.0;
  fActiveArea           = new fastjet::ActiveAreaSpec (ghostEtaMax,activeAreaRepeats,ghostArea);
  fAreaDefinition       = new fastjet::AreaDefinition (fastjet::active_area_explicit_ghosts, *fActiveArea );
  sd = new fastjet::contrib::SoftDrop(0.0,0.1,0.8);//beta_, zCut_, R0 );
}
//////////////////////////////////////////////////////////////////////////////////////////
void BoostedTTbarFlatTreeProducerGenLevel::beginJob() 
{

  //--- book the tree -----------------------
  outTree_ = fs_->make<TTree>("events","events");
  EventsHisto_ = fs_->make<TH1F>("EventHisto","EventHisto",1,0.5,1.5);
  
  outTree_->Branch("nGenJets"             ,&nGenJets_          ,"nGenJets_/I");
   
  //------------------- MC ---------------------------------

  outTree_->Branch("decay"                ,&decay_             ,"decay_/I");
  outTree_->Branch("weight"                ,&weight_             ,"weight_/D");
  
  ///////
  WBosonId_       = new std::vector<int>;
  WBosonSt_       = new std::vector<int>;
  WBosonPt_       = new std::vector<float>;
  WBosonEta_      = new std::vector<float>;
  WBosonPhi_      = new std::vector<float>;
  WBosonE_        = new std::vector<float>;
  outTree_->Branch("WBosonId"       ,"vector<int>"   ,&WBosonId_);
  outTree_->Branch("WBosonSt"       ,"vector<int>"   ,&WBosonSt_);
  outTree_->Branch("WBosonPt"       ,"vector<float>" ,&WBosonPt_);
  outTree_->Branch("WBosonEta"      ,"vector<float>" ,&WBosonEta_);
  outTree_->Branch("WBosonPhi"      ,"vector<float>" ,&WBosonPhi_);
  outTree_->Branch("WBosonE"        ,"vector<float>" ,&WBosonE_);
  
  //gen jets
  GenSoftDropMass_ = new std::vector<float>;
  GenSoftDropTau32_  = new std::vector<float>;
  GenSoftDropTau31_  = new std::vector<float>;
  GenSoftDropTau3_  = new std::vector<float>;
  GenSoftDropTau2_  = new std::vector<float>;
  GenSoftDropTau1_  = new std::vector<float>;
  GenJetpt_   = new std::vector<float>;
  GenJetphi_  = new std::vector<float>;
  GenJeteta_  = new std::vector<float>;
  GenJetenergy_   = new std::vector<float>;
  GenJetmass_   = new std::vector<float>;
  isBJetGen_   = new std::vector<bool>;
  isWJetGen_   = new std::vector<bool>;
  
  outTree_->Branch("GenSoftDropMass"       ,"vector<float>"   ,&GenSoftDropMass_);
  outTree_->Branch("GenSoftDropTau32"       ,"vector<float>"   ,&GenSoftDropTau32_);
  outTree_->Branch("GenSoftDropTau31"       ,"vector<float>"   ,&GenSoftDropTau31_);
  outTree_->Branch("GenSoftDropTau3"       ,"vector<float>"   ,&GenSoftDropTau3_);
  outTree_->Branch("GenSoftDropTau2"       ,"vector<float>"   ,&GenSoftDropTau2_);
  outTree_->Branch("GenSoftDropTau1"       ,"vector<float>"   ,&GenSoftDropTau1_);
  outTree_->Branch("GenJetpt"       ,"vector<float>"   ,&GenJetpt_);
  outTree_->Branch("GenJeteta"       ,"vector<float>"   ,&GenJeteta_);
  outTree_->Branch("GenJetphi"       ,"vector<float>"   ,&GenJetphi_);
  outTree_->Branch("GenJetenergy"       ,"vector<float>"   ,&GenJetenergy_);
  outTree_->Branch("GenJetmass"       ,"vector<float>"   ,&GenJetmass_);   
  outTree_->Branch("isBJetGen"       ,"vector<bool>"   ,&isBJetGen_);   
  outTree_->Branch("isWJetGen"       ,"vector<bool>"   ,&isWJetGen_);   
  
  GenSubJet1Pt_= new std::vector<float>;
  GenSubJet2Pt_= new std::vector<float>;
  GenSubJet1Eta_= new std::vector<float>;
  GenSubJet2Eta_= new std::vector<float>;
  GenSubJet1Phi_= new std::vector<float>;
  GenSubJet2Phi_= new std::vector<float>;
  GenSubJet1Mass_= new std::vector<float>;
  GenSubJet2Mass_= new std::vector<float>;
  
  outTree_->Branch("GenSubJet1Pt"       ,"vector<float>"   ,&GenSubJet1Pt_);
  outTree_->Branch("GenSubJet2Pt"       ,"vector<float>"   ,&GenSubJet2Pt_);
  outTree_->Branch("GenSubJet1Eta"      ,"vector<float>"   ,&GenSubJet1Eta_);
  outTree_->Branch("GenSubJet2Eta"      ,"vector<float>"   ,&GenSubJet2Eta_);
  outTree_->Branch("GenSubJet1Phi"      ,"vector<float>"   ,&GenSubJet1Phi_);
  outTree_->Branch("GenSubJet2Phi"      ,"vector<float>"   ,&GenSubJet2Phi_);
  outTree_->Branch("GenSubJet1Mass"     ,"vector<float>"   ,&GenSubJet1Mass_);
  outTree_->Branch("GenSubJet2Mass"     ,"vector<float>"   ,&GenSubJet2Mass_);

  outTree_->Branch("npu"                  ,&npu_               ,"npu_/I");
  outTree_->Branch("genEvtWeight"         ,&genEvtWeight_      ,"genEvtWeight_/F");
  outTree_->Branch("lheOriginalXWGTUP"    ,&lheOriginalXWGTUP_ ,"lheOriginalXWGTUP_/F");
  if (isHigherOrder_) {
    scaleWeights_  = new std::vector<float>;
    pdfWeights_  = new std::vector<float>;
    outTree_->Branch("scaleWeights"         ,"vector<float>"     ,&scaleWeights_);
    outTree_->Branch("pdfWeights"           ,"vector<float>"     ,&pdfWeights_);
  }

  cout<<"Begin job finished"<<endl;
}
//////////////////////////////////////////////////////////////////////////////////////////
void BoostedTTbarFlatTreeProducerGenLevel::endJob() 
{  
    
  delete WBosonSt_;
  delete WBosonId_;
  delete WBosonPt_;
  delete WBosonEta_;
  delete WBosonPhi_;
  delete WBosonE_; 

  delete GenSoftDropMass_;
  delete GenSoftDropTau32_;
  delete GenSoftDropTau3_;
  delete GenSoftDropTau2_;
  delete GenSoftDropTau1_;
  delete GenSoftDropTau31_;
  delete GenJetpt_;
  delete GenJetphi_;
  delete GenJeteta_;
  delete GenJetenergy_;
  delete GenJetmass_;
  delete isBJetGen_;
  delete isWJetGen_;

  delete GenSubJet1Pt_;
  delete GenSubJet2Pt_;
  delete GenSubJet1Eta_;
  delete GenSubJet2Eta_;
  delete GenSubJet1Phi_;
  delete GenSubJet2Phi_;
  delete GenSubJet1Mass_;
  delete GenSubJet2Mass_;

  if (isHigherOrder_) {
    delete scaleWeights_;
    delete pdfWeights_;
  }

}
//////////////////////////////////////////////////////////////////////////////////////////
void BoostedTTbarFlatTreeProducerGenLevel::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) 
{
  cout<<"Begin run finished"<<endl;
}
//////////////////////////////////////////////////////////////////////////////////////////
void BoostedTTbarFlatTreeProducerGenLevel::endRun(edm::Run const& iRun, edm::EventSetup const& iSetup) 
{
}
//////////////////////////////////////////////////////////////////////////////////////////
void BoostedTTbarFlatTreeProducerGenLevel::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) 
{
  initialize();

  //---------- mc -----------------------
  bool cut_GEN(true);

  if (!iEvent.isRealData()) { 
    iEvent.getByToken(genEvtInfoToken,genEvtInfo);
    iEvent.getByToken(genParticlesToken,genParticles);
    iEvent.getByToken(genjetsToken,genjets);

    weight_ = genEvtInfo->weight();

    LorentzVector p4T(0,0,0,0),p4Tbar(0,0,0,0);
    bool WPlusLep(false),WMinusLep(false);
    
    for(unsigned ip = 0; ip < genParticles->size(); ++ ip) {
      const GenParticle &p = (*genParticles)[ip];
      if (p.pdgId() == 24) {
	for(unsigned k = 0; k < p.numberOfDaughters(); k++) {
	  int daughterID = p.daughter(k)->pdgId();
	  if (daughterID == -11 || daughterID == -13 || daughterID == -15) {
	    WPlusLep = true;
	  }
	}
      }
      if (p.pdgId() == -24) {
	for(unsigned k = 0; k < p.numberOfDaughters(); k++) {
	  int daughterID = p.daughter(k)->pdgId();
	  if (daughterID == 11 || daughterID == 13 || daughterID == 15) {
	    WMinusLep = true;
	  }
	}
      }
      
      if (fabs(p.pdgId()) == 24) {
	WBosonId_ ->push_back(p.pdgId());
	WBosonSt_ ->push_back(p.status()); 
	WBosonPt_ ->push_back(p.pt());
	WBosonEta_->push_back(p.eta());
	WBosonPhi_->push_back(p.phi());
	WBosonE_  ->push_back(p.energy());
      }
    }// end of particle loop
      
    if (WPlusLep && WMinusLep)   decay_ = 2;
    if (WPlusLep && !WMinusLep)  decay_ = 1;
    if (!WPlusLep && WMinusLep)  decay_ = 1;
    if (!WPlusLep && !WMinusLep) decay_ = 0;
      
    iEvent.getByToken(genjetsToken,genjets);

    if (isHigherOrder_) {
      iEvent.getByToken(lheEvtInfoToken,lheEvtInfo);
      iEvent.getByToken(pupInfoToken,pupInfo);

      edm::View<PileupSummaryInfo>::const_iterator PUI;
      for(PUI = pupInfo->begin(); PUI != pupInfo->end(); ++PUI) {
	if (PUI->getBunchCrossing() == 0) {
	  npu_ = PUI->getTrueNumInteractions();
	}
      }
      
      genEvtWeight_ = genEvtInfo->weight();
      lheOriginalXWGTUP_ = lheEvtInfo->originalXWGTUP();

      for(unsigned i=0;i<lheEvtInfo->weights().size();i++) {
        string wtid(lheEvtInfo->weights()[i].id);
        float wgt(lheEvtInfo->weights()[i].wgt);
        if (wtid == "1002" || wtid == "2") scaleWeights_->push_back(wgt/lheOriginalXWGTUP_);
        if (wtid == "1003" || wtid == "3") scaleWeights_->push_back(wgt/lheOriginalXWGTUP_);
        if (wtid == "1004" || wtid == "4") scaleWeights_->push_back(wgt/lheOriginalXWGTUP_);
        if (wtid == "1005" || wtid == "5") scaleWeights_->push_back(wgt/lheOriginalXWGTUP_);
        if (wtid == "1007" || wtid == "7") scaleWeights_->push_back(wgt/lheOriginalXWGTUP_);
        if (wtid == "1009" || wtid == "9") scaleWeights_->push_back(wgt/lheOriginalXWGTUP_); 

        if ((stoi(wtid) > 2000 && stoi(wtid) <= 2102) || (stoi(wtid) > 10 && stoi(wtid) <= 110)) {
          pdfWeights_->push_back(wgt/lheOriginalXWGTUP_);
        }
      }
    }
    
    vector<LorentzVector> vP4Gen;
  
    for(GenJetCollection::const_iterator i_gen = genjets->begin(); i_gen != genjets->end(); i_gen++) {
      if (i_gen->pt() > GenptMin_ && fabs(i_gen->y()) < GenetaMax_ ) {
	
	nGenJets_++;
	
	//W jet flag
	for(unsigned ip = 0; ip < genParticles->size(); ++ ip) {
	  const GenParticle &p = (*genParticles)[ip];
	  if (fabs(p.pdgId()) == 24) {//bug it should be the absolute value!!
            double deltaEta=i_gen->eta()-p.momentum().eta();
            double deltaPhiJP=deltaPhi(i_gen->phi(),p.momentum().phi());
            double deltaRJPMin=sqrt(pow(deltaEta,2)+pow(deltaPhiJP,2));
            if(deltaRJPMin<0.6) {isWJetGen_->push_back(true);}
	    else { isWJetGen_->push_back(false); }
          }
	}
	
	GenJetpt_->push_back(i_gen->pt());
	GenJetphi_->push_back(i_gen->phi());
	GenJeteta_->push_back(i_gen->y());
	GenJetenergy_->push_back(i_gen->energy());
	GenJetmass_->push_back(i_gen->mass());

        std::vector<fastjet::PseudoJet>  lClusterParticles;
        for(unsigned int ic=0; ic<i_gen->numberOfDaughters(); ic++) {
          const reco::Candidate* gencand = i_gen->daughter(ic);
          fastjet::PseudoJet   pPart(gencand->px(),gencand->py(),gencand->pz(),gencand->energy());
          lClusterParticles.push_back(pPart);
        }
	
        fClustering = new fastjet::ClusterSequenceArea(lClusterParticles, *fAKJetDef, *fAreaDefinition);
        std::vector<fastjet::PseudoJet>  lOutJets = fClustering->inclusive_jets(20.0);
	
        if(lOutJets.size() == 0) {
          delete fClustering;
          return;
        }
	
	fastjet::PseudoJet pT1JetSD = (*sd)( lOutJets[0]);
        float iMSoftDrop = pT1JetSD.m();
		
        fastjet::contrib::NormalizedMeasure          normalizedMeasure        (1.0,0.8);
        fastjet::contrib::Njettiness routine(fastjet::contrib::Njettiness::onepass_kt_axes,normalizedMeasure);
        float iTau1 = routine.getTau(1.,lClusterParticles);                                                                                                                   
        float iTau2 = routine.getTau(2.,lClusterParticles);
        float iTau3 = routine.getTau(3.,lClusterParticles);

	float Tau32=-1;
	float Tau31=-1;

	GenSoftDropTau3_->push_back(iTau3);
	GenSoftDropTau2_->push_back(iTau2);
	GenSoftDropTau1_->push_back(iTau1);

	if(iTau2!=0) Tau32=iTau3/iTau2;
	if(iTau1!=0) Tau31=iTau3/iTau1;
	GenSoftDropTau32_->push_back(Tau32);
	GenSoftDropTau31_->push_back(Tau31);
        GenSoftDropMass_->push_back(iMSoftDrop);

	lOutJets.clear();
        delete fClustering;
	lClusterParticles.clear();
	
      }
    }

    edm::Handle<reco::JetFlavourInfoMatchingCollection> theJetFlavourInfos;
    iEvent.getByToken(jetFlavourInfosToken_, theJetFlavourInfos );
    
    for ( reco::JetFlavourInfoMatchingCollection::const_iterator j  = theJetFlavourInfos->begin();j != theJetFlavourInfos->end();++j ) {
      const reco::Jet *aJet = (*j).first.get();
      reco::JetFlavourInfo aInfo = (*j).second;

      if (aJet->pt() > GenptMin_ && fabs(aJet->y()) < GenetaMax_ ) {
        
	int FlavourGenHadron = aInfo.getHadronFlavour();
	if(abs(FlavourGenHadron)==5){ 
	  isBJetGen_->push_back(true);
	}
	else isBJetGen_->push_back(false);
      }
    }

    cut_GEN = (nGenJets_>1);
    
  }//--- end of MC -------

  EventsHisto_->Fill(1.,1.);
  
  if (cut_GEN) {
    outTree_->Fill();
  }
}
//////////////////////////////////////////////////////////////////////////////////////////
void BoostedTTbarFlatTreeProducerGenLevel::initialize()
{

  //cout<<"Ini starts"<<endl;
  nGenJets_       = 0;

  //----- MC -------
  decay_ = 10;
  weight_=1;

  WBosonSt_->clear();
  WBosonId_->clear();
  WBosonPt_->clear();
  WBosonEta_->clear();
  WBosonPhi_->clear();
  WBosonE_->clear(); 
         
  isBJetGen_->clear();
  isWJetGen_->clear();
  GenSoftDropMass_->clear();
  GenSoftDropTau32_->clear();
  GenSoftDropTau31_->clear();
  GenSoftDropTau3_->clear();
  GenSoftDropTau2_->clear();
  GenSoftDropTau1_->clear();
  GenJetpt_->clear();
  GenJetphi_->clear();
  GenJeteta_->clear();
  GenJetenergy_->clear();
  GenJetmass_->clear();
  
  GenSubJet1Pt_->clear();
  GenSubJet2Pt_->clear();
  GenSubJet1Eta_->clear();
  GenSubJet2Eta_->clear();
  GenSubJet1Phi_->clear();
  GenSubJet2Phi_->clear();
  GenSubJet1Mass_->clear();
  GenSubJet2Mass_->clear();

  if (isHigherOrder_) {
    scaleWeights_->clear();
    pdfWeights_->clear();
  }

  npu_               = -1;
  genEvtWeight_      = -999;
  lheOriginalXWGTUP_ = -999;
    
}
//////////////////////////////////////////////////////////////////////////////////////////
BoostedTTbarFlatTreeProducerGenLevel::~BoostedTTbarFlatTreeProducerGenLevel() 
{
}

DEFINE_FWK_MODULE(BoostedTTbarFlatTreeProducerGenLevel);















