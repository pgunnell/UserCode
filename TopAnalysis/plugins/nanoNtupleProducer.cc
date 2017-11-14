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

#include "KKousour/TopAnalysis/plugins/nanoNtupleProducer.h"

using namespace std;
using namespace reco;
using namespace fastjet;

nanoNtupleProducer::nanoNtupleProducer(edm::ParameterSet const& cfg)
{ 
  jetsToken             = consumes<pat::JetCollection>(cfg.getParameter<edm::InputTag>("jets"));
  genjetsToken          = consumes<GenJetCollection>(cfg.getUntrackedParameter<edm::InputTag>("genjets",edm::InputTag("")));
  metToken              = consumes<pat::METCollection>(cfg.getParameter<edm::InputTag>("met"));
  candsToken            = consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("candidates"));
  rhoToken              = consumes<double>(cfg.getParameter<edm::InputTag>("rho"));
  recVtxsToken          = consumes<reco::VertexCollection>(cfg.getParameter<edm::InputTag>("vertices"));
  triggerResultsToken   = consumes<edm::TriggerResults>(cfg.getParameter<edm::InputTag>("triggerResults"));
  triggerPrescalesToken = consumes<pat::PackedTriggerPrescales>(cfg.getParameter<edm::InputTag>("triggerPrescales"));
  pupInfoToken          = consumes<edm::View<PileupSummaryInfo> >(edm::InputTag("slimmedAddPileupInfo"));
  genEvtInfoToken       = consumes<GenEventInfoProduct>(edm::InputTag("generator"));
  genParticlesToken     = consumes<edm::View<reco::GenParticle> >(edm::InputTag("prunedGenParticles"));
  lheEvtInfoToken       = consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer"));
  runInfoToken          = consumes<LHERunInfoProduct>(edm::InputTag("externalLHEProducer"));
  triggerNames_         = cfg.getParameter<std::vector<std::string> >("triggerNames");
  srcBtag_ = cfg.getParameter<std::string>("btagger");
  etaMax_               = cfg.getParameter<double>("etaMax");
  ptMin_                = cfg.getParameter<double>("ptMin");
  ptMinLeading_         = cfg.getParameter<double>("ptMinLeading");
  massMin_              = cfg.getParameter<double>("massMin");
  btagMin_              = cfg.getParameter<double>("btagMin");
  isMC_                 = cfg.getUntrackedParameter<bool>("isMC",false);
  isPrint_              = cfg.getUntrackedParameter<bool>("isPrint",false);
  saveWeights_          = cfg.getUntrackedParameter<bool>("saveWeights",true);
  debug_                = cfg.getUntrackedParameter<bool>("debug",false);
  GenptMin_             = cfg.getUntrackedParameter<double>("GenptMin");
  GenetaMax_            = cfg.getUntrackedParameter<double>("GenetaMax");

  jetFlavourInfosToken_ = consumes<reco::JetFlavourInfoMatchingCollection>( cfg.getParameter<edm::InputTag>("jetFlavourInfos"));
}
//////////////////////////////////////////////////////////////////////////////////////////
void nanoNtupleProducer::beginJob() 
{
  //--- book the trigger histograms ---------
  triggerNamesHisto_ = fs_->make<TH1F>("TriggerNames","TriggerNames",1,0,1);
  triggerNamesHisto_->SetCanExtend(TH1::kAllAxes);
  for(unsigned i=0;i<triggerNames_.size();i++) {
    triggerNamesHisto_->Fill(triggerNames_[i].c_str(),1);
  }
  triggerPassHisto_ = fs_->make<TH1F>("TriggerPass","TriggerPass",1,0,1);
  triggerPassHisto_->SetCanExtend(TH1::kAllAxes);

  cutFlowHisto_ = fs_->make<TH1F>("CutFlow","CutFlow",1,0,1);
  cutFlowHisto_->SetCanExtend(TH1::kAllAxes);
 
  //--- book the tree -----------------------
  outTree_ = fs_->make<TTree>("events","events");
  outTree_->Branch("runNo"                ,&run_               ,"run_/I");
  outTree_->Branch("evtNo"                ,&evt_               ,"evt_/I");
  outTree_->Branch("lumi"                 ,&lumi_              ,"lumi_/I");
  outTree_->Branch("nvtx"                 ,&nVtx_              ,"nVtx_/I");
  outTree_->Branch("nJets"                ,&nJets_             ,"nJets_/I");
  if (isMC_) outTree_->Branch("nGenJets"             ,&nGenJets_          ,"nGenJets_/I");
  outTree_->Branch("nBJets"               ,&nBJets_            ,"nBJets_/I");
  outTree_->Branch("pvRho"                ,&pvRho_             ,"pvRho_/F");
  outTree_->Branch("pvz"                  ,&pvz_               ,"pvz_/F");
  outTree_->Branch("pvchi2"               ,&pvchi2_            ,"pvchi2_/F");
  outTree_->Branch("pvndof"               ,&pvndof_            ,"pvndof_/F");
  outTree_->Branch("rho"                  ,&rho_               ,"rho_/F");
  outTree_->Branch("ht"                   ,&ht_                ,"ht_/F");
  outTree_->Branch("met"                  ,&met_               ,"met_/F");
  if (isMC_) outTree_->Branch("metGen"               ,&metGen_               ,"metGen_/F");
  outTree_->Branch("metSig"               ,&metSig_            ,"metSig_/F");
  if (isMC_) outTree_->Branch("metGenSig"            ,&metGenSig_            ,"metGenSig_/F");
   
  //------------------------------------------------------------------
  isBtag_         = new std::vector<bool>;
  flavor_         = new std::vector<int>;
  flavorHadron_   = new std::vector<int>;
  pt_             = new std::vector<float>;
  unc_            = new std::vector<float>;
  cor_            = new std::vector<float>;
  btag_           = new std::vector<float>;  
  eta_            = new std::vector<float>;
  phi_            = new std::vector<float>;
  mass_           = new std::vector<float>;
  energy_         = new std::vector<float>;
  chf_            = new std::vector<float>;
  nhf_            = new std::vector<float>;
  phf_            = new std::vector<float>;
  muf_            = new std::vector<float>;
  elf_            = new std::vector<float>;

  outTree_->Branch("jetIsBtag"            ,"vector<bool>"      ,&isBtag_); 
  if (isMC_) outTree_->Branch("jetFlavor"            ,"vector<int>"       ,&flavor_);
  if (isMC_) outTree_->Branch("jetFlavorHadron"      ,"vector<int>"       ,&flavorHadron_);
  outTree_->Branch("jetPt"                ,"vector<float>"     ,&pt_);
  outTree_->Branch("jetCorr"              ,"vector<float>"     ,&cor_);
  outTree_->Branch("jetUnc"               ,"vector<float>"     ,&unc_);
  outTree_->Branch("jetBtag"              ,"vector<float>"     ,&btag_);  
  outTree_->Branch("jetEta"               ,"vector<float>"     ,&eta_);
  outTree_->Branch("jetPhi"               ,"vector<float>"     ,&phi_);
  outTree_->Branch("jetMass"              ,"vector<float>"     ,&mass_);
  outTree_->Branch("jetEnergy"            ,"vector<float>"     ,&energy_);
  outTree_->Branch("jetChf"               ,"vector<float>"     ,&chf_);
  outTree_->Branch("jetNhf"               ,"vector<float>"     ,&nhf_);
  outTree_->Branch("jetPhf"               ,"vector<float>"     ,&phf_);
  outTree_->Branch("jetMuf"               ,"vector<float>"     ,&muf_);
  outTree_->Branch("jetElf"               ,"vector<float>"     ,&elf_);

  //------------------------------------------------------------------
  triggerBit_ = new std::vector<bool>;
  triggerPre_ = new std::vector<int>;
  outTree_->Branch("triggerBit"           ,"vector<bool>"      ,&triggerBit_);
  outTree_->Branch("triggerPre"           ,"vector<int>"       ,&triggerPre_);
  //------------------- MC ---------------------------------
  if (isMC_) {
    outTree_->Branch("npu"                  ,&npu_               ,"npu_/I");
    outTree_->Branch("genEvtWeight"         ,&genEvtWeight_      ,"genEvtWeight_/F");
    outTree_->Branch("lheOriginalXWGTUP"    ,&lheOriginalXWGTUP_ ,"lheOriginalXWGTUP_/F");
    if (saveWeights_) {
      scaleWeights_  = new std::vector<float>;
      pdfWeights_  = new std::vector<float>;
      outTree_->Branch("scaleWeights"         ,"vector<float>"     ,&scaleWeights_);
      outTree_->Branch("pdfWeights"           ,"vector<float>"     ,&pdfWeights_);
    }

    GenJetpt_   = new std::vector<float>;
    GenJetphi_  = new std::vector<float>;
    GenJeteta_  = new std::vector<float>;
    GenJetenergy_   = new std::vector<float>;
    GenJetmass_   = new std::vector<float>;
    isBJetGen_   = new std::vector<bool>;

    outTree_->Branch("GenJetpt"       ,"vector<float>"   ,&GenJetpt_);
    outTree_->Branch("GenJeteta"       ,"vector<float>"   ,&GenJeteta_);
    outTree_->Branch("GenJetphi"       ,"vector<float>"   ,&GenJetphi_);
    outTree_->Branch("GenJetenergy"       ,"vector<float>"   ,&GenJetenergy_);
    outTree_->Branch("GenJetmass"       ,"vector<float>"   ,&GenJetmass_);   
    outTree_->Branch("isBJetGen"       ,"vector<bool>"   ,&isBJetGen_);   

  }
  cout<<"Begin job finished"<<endl;
}
//////////////////////////////////////////////////////////////////////////////////////////
void nanoNtupleProducer::endJob() 
{  
  delete isBtag_;
  delete flavor_;
  delete flavorHadron_;
  delete pt_;
  delete cor_;
  delete unc_;
  delete btag_;
  delete eta_;
  delete phi_;
  delete mass_;
  delete energy_;
  delete chf_;
  delete nhf_;
  delete phf_;
  delete muf_;
  delete elf_;
  delete triggerBit_;
  delete triggerPre_;
  if (isMC_) {
    if (saveWeights_) {
      delete scaleWeights_;
      delete pdfWeights_;
    }

    delete GenJetpt_;
    delete GenJetphi_;
    delete GenJeteta_;
    delete GenJetenergy_;
    delete GenJetmass_;
    delete isBJetGen_;

  }
}
//////////////////////////////////////////////////////////////////////////////////////////
void nanoNtupleProducer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) 
{
  if (isMC_ && debug_) {
    iRun.getByToken(runInfoToken,runInfo);
    for(vector<LHERunInfoProduct::Header>::const_iterator it = runInfo->headers_begin();it != runInfo->headers_end(); it++) {
      cout<<it->tag()<<endl;
      vector<string> lines = it->lines();
      for(unsigned int iLine = 0; iLine < lines.size(); iLine++) {
        cout<< lines.at(iLine);
      }
    }
  }
}
//////////////////////////////////////////////////////////////////////////////////////////
void nanoNtupleProducer::endRun(edm::Run const& iRun, edm::EventSetup const& iSetup) 
{
}
//////////////////////////////////////////////////////////////////////////////////////////
bool nanoNtupleProducer::isGoodJet(const pat::Jet &jet)
{
  bool res  = true; // by default is good, unless fails a cut bellow
  float chf = jet.chargedHadronEnergyFraction();
  float nhf = jet.neutralHadronEnergyFraction();
  float phf = jet.photonEnergyFraction();
  float muf = jet.muonEnergyFraction();
  float elf = jet.electronEnergyFraction();
  int chm   = jet.chargedHadronMultiplicity();
  int npr   = jet.neutralMultiplicity()+jet.chargedMultiplicity();
  float eta = fabs(jet.eta());
  float pt  = jet.pt();
  bool idL  = (npr>1 && phf<0.99 && nhf<0.99);
  //bool idM = (idL && ((eta<=2.4 && nhf<0.9 && phf<0.9 && elf<0.99 && muf<0.99 && chf>0 && chm>0) || eta>2.4));
  bool idT = (idL && ((eta<=2.4 && nhf<0.9 && phf<0.9 && elf<0.9 && muf<0.9 && chf>0 && chm>0) || eta>2.4));
  if (!idT) res = false;
  if (pt < ptMin_) res = false;
  if (eta > etaMax_) res = false;
  return res;
}

//////////////////////////////////////////////////////////////////////////////////////////
void nanoNtupleProducer::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) 
{
  initialize();

  if(isPrint_) cout<<"**** EVENT ****"<<endl;

  iEvent.getByToken(jetsToken,jets);
  iEvent.getByToken(metToken,met);
  iEvent.getByToken(candsToken,cands);
  iEvent.getByToken(rhoToken,rho);
  iEvent.getByToken(recVtxsToken,recVtxs);  
  iEvent.getByToken(triggerResultsToken,triggerResults);  
  iEvent.getByToken(triggerPrescalesToken,triggerPrescales); 

  //-------------- Trigger Info -----------------------------------
  triggerPassHisto_->Fill("totalEvents",1);
  const edm::TriggerNames &names = iEvent.triggerNames(*triggerResults);  
  bool passTrigger(false);
  for(unsigned int k=0;k<triggerNames_.size();k++) {
    bool bit(false);
    int pre(1);
    for(unsigned int itrig=0;itrig<triggerResults->size();itrig++) {
      string trigger_name = string(names.triggerName(itrig));
      //--- erase the last character, i.e. the version number----
      trigger_name.pop_back();
      if (trigger_name == triggerNames_[k]) {
        bit = triggerResults->accept(itrig); 
        pre = triggerPrescales->getPrescaleForIndex(itrig);
        if (bit) {
          triggerPassHisto_->Fill(triggerNames_[k].c_str(),1);
        } 
      }
    }
    //--- if at least one monitored trigger has fired passTrigger becomes true
    passTrigger += bit;
    triggerBit_->push_back(bit); 
    triggerPre_->push_back(pre);   
  }   
  vector<const reco::Candidate *> myLeptons;
  //----- at least one good vertex -----------
  bool cut_vtx = (recVtxs->size() > 0);
  if (cut_vtx) {
    pvRho_    = (*recVtxs)[0].position().Rho();
    pvz_    = (*recVtxs)[0].z();
    pvndof_ = (*recVtxs)[0].ndof();
    pvchi2_ = (*recVtxs)[0].chi2();
  }// if vtx
  //----- PF jets ------------------------------
  nJets_  = 0;
  nGenJets_  = 0;
  nBJets_ = 0;
  ht_     = 0.0;

  //edm::ESHandle<JetCorrectorParametersCollection> PFJetCorParCollCHS;
  //iSetup.get<JetCorrectionsRecord>().get("AK8PFchs",PFJetCorParCollCHS);
  //JetCorrectorParameters const& PFJetCorParCHS = (*PFJetCorParCollCHS)["Uncertainty"];
  //mPFUncCHS = new JetCorrectionUncertainty(PFJetCorParCHS);//"Summer16_23Sep2016V4_MC_Uncertainty_AK8PFchs.txt");

  double unc=0.0;
  
  vector<LorentzVector> vP4; 
  for(pat::JetCollection::const_iterator ijet =jets->begin();ijet != jets->end(); ++ijet) {  
    if (isGoodJet(*ijet)) {
      float btag= ijet->bDiscriminator(srcBtag_.c_str());
      bool isBtag = (btag >= btagMin_);
      if(isPrint_) {if (isBtag==true) cout<<"RECO "<<ijet->pt()<<" "<<ijet->eta()<<" "<<ijet->phi()<<endl;}
      flavor_        ->push_back(ijet->partonFlavour());
      flavorHadron_  ->push_back(ijet->hadronFlavour());
      chf_           ->push_back(ijet->chargedHadronEnergyFraction());
      nhf_           ->push_back(ijet->neutralHadronEnergyFraction());
      phf_           ->push_back(ijet->photonEnergyFraction());
      elf_           ->push_back(ijet->electronEnergyFraction());
      muf_           ->push_back(ijet->muonEnergyFraction());
      pt_            ->push_back(ijet->pt());
      
      //mPFUncCHS->setJetEta(ijet->eta());
      //mPFUncCHS->setJetPt(ijet->pt()); // here you must use the CORRECTED jet pt
      //unc = mPFUncCHS->getUncertainty(true);
      cor_           ->push_back(1+unc);
      unc_           ->push_back(unc);
      phi_           ->push_back(ijet->phi());
      eta_           ->push_back(ijet->eta());
      mass_          ->push_back(ijet->mass());
      energy_        ->push_back(ijet->energy());
      btag_          ->push_back(btag);
      isBtag_        ->push_back(isBtag);
      vP4.push_back(ijet->p4());
      ht_ += ijet->pt();
      nJets_++;
      if (isBtag) {
	nBJets_++;
      } 
    }// if good jet
  }// jet loop       
  rho_    = *rho;
  met_    = (*met)[0].et();
  if ((*met)[0].sumEt() > 0) {
    metSig_ = (*met)[0].et()/(*met)[0].sumEt();
  }

  if(isMC_){
    metGen_ = (*met)[0].genMET()->et();
    if ((*met)[0].genMET()->sumEt() > 0) {
      metGenSig_ = (*met)[0].genMET()->et()/(*met)[0].genMET()->sumEt();
    }
  }

  if(isPrint_) cout<<"MET "<<metGen_<<" "<<met_<<endl;

  nVtx_   = recVtxs->size();
  run_    = iEvent.id().run();
  evt_    = iEvent.id().event();
  lumi_   = iEvent.id().luminosityBlock();
  
  //---------- mc -----------------------
  bool cut_GEN(true);
  
  if (!iEvent.isRealData()) { 
    iEvent.getByToken(genEvtInfoToken,genEvtInfo);
    iEvent.getByToken(lheEvtInfoToken,lheEvtInfo);
    iEvent.getByToken(genParticlesToken,genParticles);
    iEvent.getByToken(pupInfoToken,pupInfo);
    iEvent.getByToken(genjetsToken,genjets);

    genEvtWeight_ = genEvtInfo->weight();
    lheOriginalXWGTUP_ = lheEvtInfo->originalXWGTUP();

    if (saveWeights_) {
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
 
        
    edm::View<PileupSummaryInfo>::const_iterator PUI;
    for(PUI = pupInfo->begin(); PUI != pupInfo->end(); ++PUI) {
      if (PUI->getBunchCrossing() == 0) {
        npu_ = PUI->getTrueNumInteractions();
      }
    }

    vector<LorentzVector> vP4Gen;

    for(GenJetCollection::const_iterator i_gen = genjets->begin(); i_gen != genjets->end(); i_gen++) {
      if (i_gen->pt() > GenptMin_ && fabs(i_gen->y()) < GenetaMax_ ) {

	nGenJets_++;
	
	GenJetpt_->push_back(i_gen->pt());
	GenJetphi_->push_back(i_gen->phi());
	GenJeteta_->push_back(i_gen->y());
	GenJetenergy_->push_back(i_gen->energy());
	GenJetmass_->push_back(i_gen->mass());

	vP4Gen.push_back(i_gen->p4());
        
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
	  if(isPrint_) cout<<"GEN "<<aJet->pt()<<" "<<aJet->y()<<" "<<aJet->phi()<<endl;	
	}
	else isBJetGen_->push_back(false);
      }
    }

    cut_GEN = (nGenJets_>1);

  }//--- end of MC -------

  bool cut_RECO = (nJets_ > 0);  
 
  cutFlowHisto_->Fill("All",1);
  if (iEvent.isRealData()) {
    if (cut_RECO) {
      cutFlowHisto_->Fill("Two jets",1);
      if (passTrigger) {
        cutFlowHisto_->Fill("Trigger",1);
        outTree_->Fill();
      }
    }
  } 
  else {
    if (cut_RECO || cut_GEN) {
      cutFlowHisto_->Fill("Two jets (reco || gen)",1);
      outTree_->Fill();
    }
  }

}
//////////////////////////////////////////////////////////////////////////////////////////
void nanoNtupleProducer::initialize()
{
  run_            = -1;
  evt_            = -1;
  lumi_           = -1;
  nVtx_           = -1;
  nJets_          = -1;
  nGenJets_       = -1;
  nBJets_         = -1;
  rho_            = -1;
  met_            = -1;
  metGen_         = -1;
  metSig_         = -1;
  metGenSig_      = -1;
  ht_             = -1;
  pvRho_          = -999;
  pvz_            = -999;
  pvndof_         = -999;
  pvchi2_         = -999;
  flavor_         ->clear();
  flavorHadron_   ->clear();
  pt_             ->clear();
  cor_             ->clear();
  unc_             ->clear();
  eta_            ->clear();
  phi_            ->clear();
  mass_           ->clear();
  energy_         ->clear();
  chf_            ->clear();
  nhf_            ->clear();
  phf_            ->clear();
  elf_            ->clear();
  muf_            ->clear();
  btag_           ->clear();
  isBtag_         ->clear();
  triggerBit_     ->clear();
  triggerPre_     ->clear();
    
  //----- MC -------
  if (isMC_) {
    npu_ = -1;
    genEvtWeight_ = -999;
    lheOriginalXWGTUP_ = -999;
    if (saveWeights_) {
      scaleWeights_->clear();
      pdfWeights_->clear();
    }
    
    isBJetGen_->clear();
    
    GenJetpt_->clear();
    GenJetphi_->clear();
    GenJeteta_->clear();
    GenJetenergy_->clear();
    GenJetmass_->clear();
    
  }
}
//////////////////////////////////////////////////////////////////////////////////////////
nanoNtupleProducer::~nanoNtupleProducer() 
{
}

DEFINE_FWK_MODULE(nanoNtupleProducer);
















