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
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

#include "KKousour/TopAnalysis/plugins/nanoNtupleProducerCandidates.h"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/contrib/SoftDrop.hh"

using namespace std;
using namespace reco;
using namespace fastjet;

nanoNtupleProducerCandidates::nanoNtupleProducerCandidates(edm::ParameterSet const& cfg)
{ 
  jetsToken             = consumes<pat::JetCollection>(cfg.getParameter<edm::InputTag>("jets"));
  genjetsToken          = consumes<GenJetCollection>(cfg.getUntrackedParameter<edm::InputTag>("genjets",edm::InputTag("")));
  metToken              = consumes<pat::METCollection>(cfg.getParameter<edm::InputTag>("met"));
  candsToken            = consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("candidates"));
  rhoToken              = consumes<double>(cfg.getParameter<edm::InputTag>("rho"));
  recVtxsToken          = consumes<reco::VertexCollection>(cfg.getParameter<edm::InputTag>("vertices"));
  triggerResultsToken   = consumes<edm::TriggerResults>(cfg.getParameter<edm::InputTag>("triggerResults"));
  triggerPrescalesToken = consumes<pat::PackedTriggerPrescales>(cfg.getParameter<edm::InputTag>("triggerPrescales"));
  triggerObjects_ = consumes<pat::TriggerObjectStandAloneCollection>(cfg.getParameter<edm::InputTag> ("triggerObjects"));
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
  isHiggs_              = cfg.getUntrackedParameter<bool>("isHiggs",false);
  mHigherOrder_         = cfg.getUntrackedParameter<bool>("HigherOrder",false);
  isPrint_              = cfg.getUntrackedParameter<bool>("isPrint",false);
  saveWeights_          = cfg.getUntrackedParameter<bool>("saveWeights",true);
  debug_                = cfg.getUntrackedParameter<bool>("debug",false);
  GenptMin_             = cfg.getUntrackedParameter<double>("GenptMin");
  GenetaMax_            = cfg.getUntrackedParameter<double>("GenetaMax");
  muonsToken            = consumes<pat::MuonCollection>(cfg.getParameter<edm::InputTag>("muons"));
  electronsToken = consumes<pat::ElectronCollection>(cfg.getParameter<edm::InputTag>("electrons"));

  jetFlavourInfosToken_ = consumes<reco::JetFlavourInfoMatchingCollection>( cfg.getParameter<edm::InputTag>("jetFlavourInfos"));

  //define jet radius and jet definition for fast jets
  jetRadius_              = cfg.getParameter<double>("jetRadius");

  //Gen Jet information
  fAKJetDef = new fastjet::JetDefinition(fastjet::antikt_algorithm, jetRadius_);
  int activeAreaRepeats = 1;
  double ghostArea      = 0.01;
  double ghostEtaMax    = 7.0;
  fActiveArea           = new fastjet::ActiveAreaSpec (ghostEtaMax,activeAreaRepeats,ghostArea);
  fAreaDefinition       = new fastjet::AreaDefinition (fastjet::active_area_explicit_ghosts, *fActiveArea );
  sd = new fastjet::contrib::SoftDrop(0.0,0.1,jetRadius_);//beta_, zCut_, R0 );
}
//////////////////////////////////////////////////////////////////////////////////////////
void nanoNtupleProducerCandidates::beginJob() 
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
  MatchedLeptons_ = new std::vector<bool>;
  MatchedHiggs_   = new std::vector<bool>;
  WPlusLep_       = new std::vector<bool>;
  WMinusLep_       = new std::vector<bool>;
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
  chm_            = new std::vector<float>;
  npr_            = new std::vector<float>;
  nhf_            = new std::vector<float>;
  phf_            = new std::vector<float>;
  muf_            = new std::vector<float>;
  elf_            = new std::vector<float>;

  outTree_->Branch("jetIsBtag"            ,"vector<bool>"      ,&isBtag_); 
  if (isMC_) outTree_->Branch("jetFlavor"            ,"vector<int>"       ,&flavor_);
  if (isMC_) outTree_->Branch("jetFlavorHadron"      ,"vector<int>"       ,&flavorHadron_);
  if (isMC_) outTree_->Branch("MatchedLeptons"      ,"vector<bool>"       ,&MatchedLeptons_);
  if (isMC_) outTree_->Branch("MatchedHiggs"      ,"vector<bool>"       ,&MatchedHiggs_);
  if (isMC_) outTree_->Branch("WMinusLep"      ,"vector<bool>"       ,&WMinusLep_);
  if (isMC_) outTree_->Branch("WPlusLep"      ,"vector<bool>"       ,&WPlusLep_);

  outTree_->Branch("jetPt"                ,"vector<float>"     ,&pt_);
  outTree_->Branch("jetCorr"              ,"vector<float>"     ,&cor_);
  outTree_->Branch("jetUnc"               ,"vector<float>"     ,&unc_);
  outTree_->Branch("jetBtag"              ,"vector<float>"     ,&btag_);  
  outTree_->Branch("jetEta"               ,"vector<float>"     ,&eta_);
  outTree_->Branch("jetPhi"               ,"vector<float>"     ,&phi_);
  outTree_->Branch("jetMass"              ,"vector<float>"     ,&mass_);
  outTree_->Branch("jetEnergy"            ,"vector<float>"     ,&energy_);
  outTree_->Branch("jetChf"               ,"vector<float>"     ,&chf_);
  outTree_->Branch("jetChm"               ,"vector<float>"     ,&chm_);
  outTree_->Branch("jetNpr"               ,"vector<float>"     ,&npr_);
  outTree_->Branch("jetNhf"               ,"vector<float>"     ,&nhf_);
  outTree_->Branch("jetPhf"               ,"vector<float>"     ,&phf_);
  outTree_->Branch("jetMuf"               ,"vector<float>"     ,&muf_);
  outTree_->Branch("jetElf"               ,"vector<float>"     ,&elf_);

  massSoftDrop_   = new std::vector<float>;
  tau1_           = new std::vector<float>;
  tau2_           = new std::vector<float>;
  tau3_           = new std::vector<float>;
  tau4_           = new std::vector<float>;

  outTree_->Branch("jetMassSoftDrop"      ,"vector<float>"     ,&massSoftDrop_);
  outTree_->Branch("jetEnergy"            ,"vector<float>"     ,&energy_);
  outTree_->Branch("jetChf"               ,"vector<float>"     ,&chf_);
  outTree_->Branch("jetNhf"               ,"vector<float>"     ,&nhf_);
  outTree_->Branch("jetPhf"               ,"vector<float>"     ,&phf_);
  outTree_->Branch("jetMuf"               ,"vector<float>"     ,&muf_);
  outTree_->Branch("jetElf"               ,"vector<float>"     ,&elf_);
  outTree_->Branch("jetTau1"              ,"vector<float>"     ,&tau1_);
  outTree_->Branch("jetTau2"              ,"vector<float>"     ,&tau2_);
  outTree_->Branch("jetTau3"              ,"vector<float>"     ,&tau3_);
  outTree_->Branch("jetTau4"              ,"vector<float>"     ,&tau4_);

  //Candidate information
  candPdgId_ = new std::vector<float>;
  candPhi_ = new std::vector<float>;
  candEta_ = new std::vector<float>;
  candPt_ = new std::vector<float>;
  candMass_ = new std::vector<float>;
  candPuppiWeight_ = new std::vector<float>;
  candPvAssociationQuality_ = new std::vector<float>;
  candDXY_ = new std::vector<float>;
  candDZ_ = new std::vector<float>;
  candDZAssociatedPV_ = new std::vector<float>;
  candPx_ = new std::vector<float>;
  candPy_ = new std::vector<float>;
  candPz_ = new std::vector<float>;
  candEnergy_ = new std::vector<float>;
  candSubJetPart_ = new std::vector<float>;
  ncandidates_ = new std::vector<int>;
  nGencandidates_ = new std::vector<int>;

  btagSub0_= new std::vector<float>;
  btagSub1_= new std::vector<float>;
  massSub0_= new std::vector<float>;
  massSub1_= new std::vector<float>;
  ptSub0_= new std::vector<float>;
  ptSub1_= new std::vector<float>;
  etaSub0_= new std::vector<float>;
  etaSub1_= new std::vector<float>;
  phiSub0_= new std::vector<float>;
  phiSub1_= new std::vector<float>;
  flavorSub0_= new std::vector<int>;
  flavorSub1_= new std::vector<int>;
  flavorHadronSub0_= new std::vector<int>;
  flavorHadronSub1_= new std::vector<int>;
  nSubJets_= new std::vector<int>;
  nBSubJets_= new std::vector<int>;
  
  outTree_->Branch("btagSub0"        ,"vector<float>"       ,&btagSub0_);
  outTree_->Branch("btagSub1"        ,"vector<float>"       ,&btagSub1_);
  outTree_->Branch("massSub0"        ,"vector<float>"       ,&massSub0_);
  outTree_->Branch("massSub1"        ,"vector<float>"       ,&massSub1_);
  outTree_->Branch("ptSub0"        ,"vector<float>"       ,&ptSub0_);
  outTree_->Branch("ptSub1"        ,"vector<float>"       ,&ptSub1_);
  outTree_->Branch("etaSub0"        ,"vector<float>"       ,&etaSub0_);
  outTree_->Branch("etaSub1"        ,"vector<float>"       ,&etaSub1_);
  outTree_->Branch("phiSub0"        ,"vector<float>"       ,&phiSub0_);
  outTree_->Branch("phiSub1"        ,"vector<float>"       ,&phiSub1_);
  outTree_->Branch("flavorSub0"        ,"vector<int>"       ,&flavorSub0_);
  outTree_->Branch("flavorSub1"        ,"vector<int>"       ,&flavorSub1_);
  outTree_->Branch("flavorHadronSub0"        ,"vector<int>"       ,&flavorHadronSub0_);
  outTree_->Branch("flavorHadronSub1"        ,"vector<int>"       ,&flavorHadronSub1_);
  outTree_->Branch("nSubJets"        ,"vector<int>"       ,&nSubJets_);
  outTree_->Branch("nBSubJets"        ,"vector<int>"       ,&nBSubJets_);

  outTree_->Branch("candPdgId"        ,"vector<float>"       ,&candPdgId_);
  outTree_->Branch("candPhi"        ,"vector<float>"       ,&candPhi_);
  outTree_->Branch("candEta"        ,"vector<float>"       ,&candEta_);
  outTree_->Branch("candPt"        ,"vector<float>"       ,&candPt_);
  outTree_->Branch("candMass"        ,"vector<float>"       ,&candMass_);
  outTree_->Branch("candPx"        ,"vector<float>"       ,&candPx_);
  outTree_->Branch("candPy"        ,"vector<float>"       ,&candPy_);
  outTree_->Branch("candPz"        ,"vector<float>"       ,&candPz_);
  outTree_->Branch("candEnergy"        ,"vector<float>"       ,&candEnergy_);
  outTree_->Branch("candPuppiWeight"        ,"vector<float>"       ,&candPuppiWeight_);

  outTree_->Branch("candPvAssociationQuality"        ,"vector<float>"       ,&candPvAssociationQuality_);
  outTree_->Branch("candDXY"        ,"vector<float>"       ,&candDXY_);
  outTree_->Branch("candDZ"        ,"vector<float>"       ,&candDZ_);
  outTree_->Branch("candDZAssociatedPV"        ,"vector<float>"       ,&candDZAssociatedPV_);
  outTree_->Branch("candSubJetPart"        ,"vector<float>"       ,&candSubJetPart_);
  outTree_->Branch("ncandidates"        ,"vector<int>"       ,&ncandidates_);
  outTree_->Branch("nGencandidates"        ,"vector<int>"       ,&nGencandidates_);

  //------------------------------------------------------------------
  triggerBit_ = new std::vector<bool>;
  triggerPre_ = new std::vector<int>;
  outTree_->Branch("triggerBit"           ,"vector<bool>"      ,&triggerBit_);
  outTree_->Branch("triggerPre"           ,"vector<int>"       ,&triggerPre_);

  TrigObjpt_            = new std::vector<double>;
  TrigObjeta_            = new std::vector<double>;
  TrigObjphi_            = new std::vector<double>;
  TrigObjcollection_            = new std::vector<string>;

  outTree_->Branch("TrigObjpt"           ,"vector<double>"      ,&TrigObjpt_);
  outTree_->Branch("TrigObjeta"           ,"vector<double>"       ,&TrigObjeta_);
  outTree_->Branch("TrigObjphi"           ,"vector<double>"       ,&TrigObjphi_);
  outTree_->Branch("TrigObjcollection"           ,"vector<string>"       ,&TrigObjcollection_);
  outTree_->Branch("nTriggerObjects"      ,&nTriggerObjects_  ,"nTriggerObjects_/I");

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

    //gen substructure + candidates
    GenJetTau1_ = new std::vector<float>;
    GenJetTau2_ = new std::vector<float>;
    GenJetTau3_ = new std::vector<float>;
    GenJetTau4_ = new std::vector<float>;
    GenSoftDropMass_ = new std::vector<float>;
    
    outTree_->Branch("GenJetTau1"       ,"vector<float>"   ,&GenJetTau1_);
    outTree_->Branch("GenJetTau2"       ,"vector<float>"   ,&GenJetTau2_);
    outTree_->Branch("GenJetTau3"       ,"vector<float>"   ,&GenJetTau3_);
    outTree_->Branch("GenJetTau4"       ,"vector<float>"   ,&GenJetTau4_);
    outTree_->Branch("GenSoftDropMass"       ,"vector<float>"   ,&GenSoftDropMass_);

    GenCandPhi_ = new std::vector<float>;
    GenCandEta_ = new std::vector<float>;
    GenCandPx_ = new std::vector<float>;
    GenCandPy_ = new std::vector<float>;
    GenCandPz_ = new std::vector<float>;
    GenCandPt_ = new std::vector<float>;
    GenCandEnergy_ = new std::vector<float>;
    
    outTree_->Branch("GenCandPhi"        ,"vector<float>"       ,&GenCandPhi_);
    outTree_->Branch("GenCandEta"        ,"vector<float>"       ,&GenCandEta_);
    outTree_->Branch("GenCandPt"        ,"vector<float>"       ,&GenCandPt_);
    outTree_->Branch("GenCandPx"        ,"vector<float>"       ,&GenCandPx_);
    outTree_->Branch("GenCandPy"        ,"vector<float>"       ,&GenCandPy_);
    outTree_->Branch("GenCandPz"        ,"vector<float>"       ,&GenCandPz_);
    outTree_->Branch("GenCandEnergy"        ,"vector<float>"       ,&GenCandEnergy_);
  }

  cout<<"Begin job finished"<<endl;
}
//////////////////////////////////////////////////////////////////////////////////////////
void nanoNtupleProducerCandidates::endJob() 
{  
  //candidate information
  delete candPdgId_;
  delete candPhi_;
  delete candEta_;
  delete candPt_;
  delete candPx_;
  delete candPy_;
  delete candPz_;
  delete candPuppiWeight_;
  delete candMass_;
  delete candPvAssociationQuality_;
  delete candDXY_;
  delete candDZ_;
  delete candDZAssociatedPV_;
  delete candSubJetPart_;
  delete ncandidates_;
  delete nGencandidates_;
  //new ones
  delete massSoftDrop_;
  delete tau1_;
  delete tau2_;
  delete tau3_;
  delete tau4_;

  delete MatchedLeptons_;
  delete MatchedHiggs_;
  delete WPlusLep_;
  delete WMinusLep_;
  //subjets
  delete btagSub0_;
  delete btagSub1_;
  delete massSub0_;
  delete massSub1_;
  delete ptSub0_;
  delete ptSub1_;
  delete etaSub0_;
  delete etaSub1_;
  delete phiSub0_;
  delete phiSub1_;
  delete flavorSub0_;
  delete flavorSub1_;
  delete flavorHadronSub0_;
  delete flavorHadronSub1_;
  delete nSubJets_;
  delete nBSubJets_;

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
  delete chm_;
  delete npr_;
  delete nhf_;
  delete phf_;
  delete muf_;
  delete elf_;
  delete triggerBit_;
  delete triggerPre_;
  delete TrigObjpt_;
  delete TrigObjeta_;
  delete TrigObjphi_;
  delete TrigObjcollection_;
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

    delete GenJetTau1_;
    delete GenJetTau2_;
    delete GenJetTau3_;
    delete GenJetTau4_;
    delete GenSoftDropMass_;

    delete GenCandPhi_;
    delete GenCandEta_;
    delete GenCandPx_;
    delete GenCandPy_;
    delete GenCandPz_;
    delete GenCandPt_;
    delete GenCandEnergy_;


  }
}
//////////////////////////////////////////////////////////////////////////////////////////
void nanoNtupleProducerCandidates::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) 
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
void nanoNtupleProducerCandidates::endRun(edm::Run const& iRun, edm::EventSetup const& iSetup) 
{
}
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
float nanoNtupleProducerCandidates::getPFMiniIsolation(edm::Handle<pat::PackedCandidateCollection> pfcands,const reco::Candidate *cand)
{
  if (cand->pt() < 5.0) return 99999;
  float deadcone_ch(0.0),deadcone_nh(0.0),deadcone_ph(0.0),deadcone_pu(0.0);
  
  if (cand->isElectron()) {
    if (fabs(cand->eta()) > 1.479) {
      deadcone_ch = 0.015;
      deadcone_nh = 0.0;
      deadcone_ph = 0.08;
      deadcone_pu = 0.015;
    }
  }
  if (cand->isMuon()) {
    deadcone_ch = 0.0001 ;
    deadcone_nh = 0.01;
    deadcone_ph = 0.01;
    deadcone_pu = 0.01;
  }

  float r_iso(0.2);

  if (cand->pt() > 50 && cand->pt() < 200) {
    r_iso = 10./cand->pt();
  }
  if (cand->pt() >= 200) {
    r_iso = 0.05;
  }

  float iso_ch(0.0),iso_nh(0.0),iso_ph(0.0),iso_pu(0.0);

  for(const pat::PackedCandidate &pfc : *pfcands) {
    if (abs(pfc.pdgId()) < 7) continue;
    float dr = reco::deltaR(pfc,*cand);
    if (dr > r_iso) continue;
    //--- NEUTRALS ----------
    if (pfc.charge() == 0) {
      if (pfc.pt() > 0.5) {
        //--- PHOTONS ------
        if (abs(pfc.pdgId()) == 22) {
          if (dr < deadcone_ph) continue;
          iso_ph += pfc.pt();
        }
        //--- HADRONS ------
        if (abs(pfc.pdgId()) == 130) {
          if (dr < deadcone_nh) continue;
          iso_nh += pfc.pt();
        }  
      } 
    }
    //--- CHARGED FROM PV -----
    else if (pfc.fromPV() > 1) {
      if (fabs(pfc.pdgId()) == 211) {
        if (dr < deadcone_ch) continue;
        iso_ch += pfc.pt();
      }
    }
    //--- CHARGED FROM PU -----
    else {
      if (pfc.pt() > 0.5) {
        if (dr < deadcone_pu) continue;
        iso_pu += pfc.pt();
      }
    }
  }
  float iso = iso_ch + std::max(0.0,iso_ph + iso_nh - 0.5*iso_pu);
  return iso;
}
//////////////////////////////////////////////////////////////////////////////////////////
bool nanoNtupleProducerCandidates::isGoodElectron(const pat::Electron &el,const reco::Vertex &vtx,edm::Handle<pat::PackedCandidateCollection> pfcands)
{
  bool res = true; // by default is good, unless fails a cut bellow
  if(el.pt() < 10) res = false;
  if(fabs(el.eta()) > 3.5 && res == true) res = false;
  bool isEB = fabs(el.superCluster()->eta()) < 1.479 ? 1 : 0;
  bool isEE = fabs(el.superCluster()->eta()) > 1.479 ? 1 : 0;
  if(res) {
    float trackMomentumAtVtx = (float)sqrt(el.trackMomentumAtVtx().mag2());
    float ecalEnergy = (float)el.ecalEnergy();
    float full5x5_sigmaIetaIeta = (float)el.full5x5_sigmaIetaIeta();
    float dEtaIn = (float)el.deltaEtaSuperClusterTrackAtVtx();
    float dPhiIn = (float)el.deltaPhiSuperClusterTrackAtVtx();
    float HoE = (float)el.hadronicOverEm();
    float ooEmooP = (float)fabs(1/ecalEnergy - 1/trackMomentumAtVtx);
    float d0 = (float)el.gsfTrack()->dxy(vtx.position());
    float dz = (float)el.gsfTrack()->dz(vtx.position());
    int expectedMissingInnerHits = el.gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS);
    bool passConversionVeto = el.passConversionVeto();
    if(isEB) {// tight working point
      if(res && full5x5_sigmaIetaIeta > 0.0128) res = false;
      if(res && fabs(dEtaIn) > 0.00353) res = false;
      if(res && fabs(dPhiIn) > 0.0499) res = false;
      if(res && HoE > 0.026) res = false;
      if(res && ooEmooP > 0.0361) res = false;
      if(res && fabs(d0) > 0.005) res = false;
      if(res && fabs(dz) > 0.1) res = false;
      if(res && expectedMissingInnerHits >= 2 ) res = false;
      if(res && passConversionVeto == false ) res = false;
    }
    if(isEE) {// tight working point
      if(res && full5x5_sigmaIetaIeta > 0.0445) res = false;
      if(res && fabs(dEtaIn) > 0.00984) res = false;
      if(res && fabs(dPhiIn) > 0.157) res = false;
      if(res && HoE > 0.0615) res = false;
      if(res && ooEmooP > 0.00999) res = false;
      if(res && fabs(d0) > 0.1) res = false;
      if(res && fabs(dz) > 0.2) res = false;
      if(res && expectedMissingInnerHits >= 4 ) res = false;
      if(res && passConversionVeto == false ) res = false;
    }
  }
  if(res && getPFMiniIsolation(pfcands,(reco::Candidate*)&el)/el.pt() > 0.1) res = false;
  return res;
}

//////////////////////////////////////////////////////////////////////////////////////////
bool nanoNtupleProducerCandidates::isGoodJet(const pat::Jet &jet)
{
  bool res  = true; // by default is good, unless fails a cut bellow
  //float chf = jet.chargedHadronEnergyFraction();
  //float nhf = jet.neutralHadronEnergyFraction();
  //float phf = jet.photonEnergyFraction();
  //float muf = jet.muonEnergyFraction();
  //float elf = jet.electronEnergyFraction();
  //int chm   = jet.chargedHadronMultiplicity();
  //int npr   = jet.neutralMultiplicity()+jet.chargedMultiplicity();
  float eta = fabs(jet.eta());
  float pt  = jet.pt();
  //bool idL  = (npr>1 && phf<0.99 && nhf<0.99);
  //bool idM = (idL && ((eta<=2.4 && nhf<0.9 && phf<0.9 && elf<0.99 && muf<0.99 && chf>0 && chm>0) || eta>2.4));
  //bool idT = (idL && ((eta<=2.4 && nhf<0.9 && phf<0.9 && elf<0.9 && muf<0.9 && chf>0 && chm>0) || eta>2.4));
  //if (!idT) res = false;
  if (pt < ptMin_) res = false;
  if (eta > etaMax_) res = false;
  return res;
}

//////////////////////////////////////////////////////////////////////////////////////////
void nanoNtupleProducerCandidates::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) 
{
  initialize();

  iEvent.getByToken(jetsToken,jets);
  iEvent.getByToken(metToken,met);
  iEvent.getByToken(candsToken,cands);
  iEvent.getByToken(rhoToken,rho);
  iEvent.getByToken(recVtxsToken,recVtxs);  

  bool passTrigger(false);

  if(!isMC_){
    iEvent.getByToken(triggerResultsToken,triggerResults);  
    iEvent.getByToken(triggerPrescalesToken,triggerPrescales); 
    
    edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
    iEvent.getByToken(triggerObjects_, triggerObjects);
    
    //-------------- Trigger Info -----------------------------------
    triggerPassHisto_->Fill("totalEvents",1);
    const edm::TriggerNames &names = iEvent.triggerNames(*triggerResults);  
        
    int numTriggerObjects=0;
    
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
    
    for (pat::TriggerObjectStandAlone obj : *triggerObjects) { // note: not "const &" since we want to call unpackPathNames
      obj.unpackPathNames(names);
      //string str2="hltAK4PF";
      //string str3="hltPFJets";
      //string str4="hltAK8";
      bool save_object=false;
      
      for (unsigned h = 0; h < obj.filterIds().size(); ++h) {
	if(fabs(obj.filterIds()[h]) == 85) save_object = true; }
      //if(obj.collection().find(str2)!=string::npos || obj.collection().find(str3)!=string::npos){// || obj.collection().find(str4)!=string::npos){
      
      if(save_object){
	numTriggerObjects++;
	
	TrigObjpt_->push_back(obj.pt());
	TrigObjeta_->push_back(obj.eta());
	TrigObjphi_->push_back(obj.phi());
	//TrigObjcollection_->push_back(obj.collection());
	//cout<<obj.collection()<<" ";
	//}
      }
      
      //for (unsigned h = 0; h < obj.filterIds().size(); ++h) std::cout << " Trigger id " << obj.filterIds()[h]<<" "<<endl;
    }
   
    nTriggerObjects_=numTriggerObjects;
  }

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

  double unc=0.0;
  int i=0;

  iEvent.getByToken(muonsToken,muons);
  iEvent.getByToken(electronsToken,electrons);

  vector<const reco::Candidate *> myLeptons;

  for (const pat::Muon &mu : *muons) {
    if (mu.isTightMuon((*recVtxs)[0]) && mu.pt()>10) myLeptons.push_back(&mu);
  }
  
  for (const pat::Electron &el : *electrons) {
    if (isGoodElectron(el,(*recVtxs)[0],cands)) myLeptons.push_back(&el);
  }
  std::sort(myLeptons.begin(),myLeptons.end(),[](const reco::Candidate *a,const reco::Candidate *b){return a->pt() > b->pt();});
  
  if(isMC_) iEvent.getByToken(genParticlesToken,genParticles);

  vector<LorentzVector> vP4; 
  for(pat::JetCollection::const_iterator ijet =jets->begin();ijet != jets->end(); ++ijet) {  
    if (isGoodJet(*ijet)) {

      float DRmax = jetRadius_;
      bool isHiggsMatched=false;
      bool WPlusLep=false;
      bool WMinusLep=false;
      
      if(isMC_){
	
	for(unsigned ip = 0; ip < genParticles->size(); ++ ip) {
	  const GenParticle &p = (*genParticles)[ip];
	  if (p.pdgId() == 25) {
	    if( deltaR(p.p4().eta(),p.p4().phi(),ijet->eta(),ijet->phi()) < DRmax) isHiggsMatched=true;
	  }
	  
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
	}
      }
    
      bool isLeptonMatched = false;

      for(auto & lep: myLeptons) if( deltaR(lep->eta(),lep->phi(),ijet->eta(),ijet->phi()) < DRmax ) isLeptonMatched = true;
      
      if(!isHiggs_) isHiggsMatched=true; //if it's not the signal, we take the leading jet

      if(isHiggsMatched){
      
	MatchedLeptons_->push_back(isLeptonMatched);
          
	float btag= ijet->bDiscriminator(srcBtag_.c_str());
	bool isBtag = (btag >= btagMin_);
	//if (bittest) cout<<"RECO "<<ijet->pt()<<" "<<ijet->eta()<<" "<<ijet->phi()<<endl;
	flavor_        ->push_back(ijet->partonFlavour());
	flavorHadron_  ->push_back(ijet->hadronFlavour());
	if(ijet->isPFJet()){
	  chf_           ->push_back(ijet->chargedHadronEnergyFraction());
	  nhf_           ->push_back(ijet->neutralHadronEnergyFraction());
	  phf_           ->push_back(ijet->photonEnergyFraction());
	  elf_           ->push_back(ijet->electronEnergyFraction());
	  muf_           ->push_back(ijet->muonEnergyFraction());
	  chm_           ->push_back(ijet->chargedHadronMultiplicity());
	  npr_           ->push_back(ijet->neutralMultiplicity()+ijet->chargedMultiplicity());
	}
      
	MatchedHiggs_->push_back(isHiggsMatched);
	WMinusLep_->push_back(WMinusLep);
	WPlusLep_->push_back(WPlusLep);

	pt_            ->push_back(ijet->pt());
	
	/*massSoftDrop_ ->push_back(ijet->userFloat("ak8PFJetsPuppiSoftDropMass"));
	tau1_          ->push_back(ijet->userFloat("NjettinessAK8Puppi:tau1"));
	tau2_          ->push_back(ijet->userFloat("NjettinessAK8Puppi:tau2"));
	tau3_ ->push_back(ijet->userFloat("NjettinessAK8Puppi:tau3"));
	tau4_ ->push_back(ijet->userFloat("NjettinessAK8Puppi:tau4"));*/
	
	//---- subjets --------------------
      	
	/*int nSub((ijet->subjets("SoftDropPuppi")).size());
	int nBSub(0);
	if (nSub > 0) {
	  btagSub0_->push_back((ijet->subjets("SoftDropPuppi"))[0]->bDiscriminator(srcBtag_.c_str()));
	  massSub0_->push_back((ijet->subjets("SoftDropPuppi"))[0]->mass());
	  ptSub0_->push_back((ijet->subjets("SoftDropPuppi"))[0]->pt());
	  etaSub0_->push_back((ijet->subjets("SoftDropPuppi"))[0]->eta());
	  phiSub0_->push_back((ijet->subjets("SoftDropPuppi"))[0]->phi());
	  if (isMC_) flavorSub0_->push_back((ijet->subjets("SoftDropPuppi"))[0]->partonFlavour());
	  if (isMC_) flavorHadronSub0_->push_back((ijet->subjets("SoftDropPuppi"))[0]->hadronFlavour());
	  if ((ijet->subjets("SoftDropPuppi"))[0]->bDiscriminator(srcBtag_.c_str()) >= btagMin_) {
	    nBSub++;
	  }
	  if (nSub > 1) {
	    btagSub1_->push_back((ijet->subjets("SoftDropPuppi"))[1]->bDiscriminator(srcBtag_.c_str()));
	    massSub1_->push_back((ijet->subjets("SoftDropPuppi"))[1]->mass());
	    ptSub1_->push_back((ijet->subjets("SoftDropPuppi"))[1]->pt());
	    etaSub1_->push_back((ijet->subjets("SoftDropPuppi"))[1]->eta());
	    phiSub1_->push_back((ijet->subjets("SoftDropPuppi"))[1]->phi());
	    if (isMC_) flavorSub1_->push_back((ijet->subjets("SoftDropPuppi"))[1]->partonFlavour());
	    if (isMC_) flavorHadronSub1_->push_back((ijet->subjets("SoftDropPuppi"))[1]->hadronFlavour());
	    if ((ijet->subjets("SoftDropPuppi"))[1]->bDiscriminator(srcBtag_.c_str()) >= btagMin_) {
	      nBSub++;
	    }
	  } 
	}

	nSubJets_->push_back(nSub);
	nBSubJets_->push_back(nBSub);*/
    
	int ncandidatesjet=0;

	//Include constituents of the jets (you have that the first two constituents are the subjets which contain further constituents (the particles that defined the jets) and the particles that are excluded by the two subjets are listed afterwards

	for ( unsigned ida = 0; ida < ijet->numberOfDaughters(); ++ida ) {
	  reco::Candidate const * cand = ijet->daughter(ida);
	  
	  if ( cand->numberOfDaughters() == 0 ){
	    
	    for(const pat::PackedCandidate &pfc : *cands) {
	      if(pfc.pt()==cand->pt() && pfc.eta()==cand->eta() && pfc.phi()==cand->phi()){
		candPvAssociationQuality_->push_back(pfc.pvAssociationQuality());
		candDXY_->push_back(pfc.dxy());
		candDZ_->push_back(pfc.dz());
		candPuppiWeight_->push_back(pfc.puppiWeight());
		candDZAssociatedPV_->push_back(pfc.dzAssociatedPV());
		//break;
		//candVertex_->push_back(cand2.vertex());//it is a 3D vector
	      }
	    }
	      
	    candPx_->push_back(cand->px());
	    candPy_->push_back(cand->py());
	    candPz_->push_back(cand->pz());
	    candEnergy_->push_back(cand->energy());
	    
	    candPdgId_->push_back(cand->pdgId());
	    candPhi_->push_back(cand->phi());
	    candEta_->push_back(cand->eta());
	    candPt_->push_back(cand->pt());
	    candMass_->push_back(cand->mass()); // pfc.p4().mass()
	    candSubJetPart_->push_back(0);
	    ncandidatesjet++;
	  }
	
	  else {
	    for ( unsigned jda = 0; jda < cand->numberOfDaughters(); ++jda ) {
	      reco::Candidate const * cand2 = cand->daughter(jda);
	      candPdgId_->push_back(cand2->pdgId());
	      candPhi_->push_back(cand2->phi());
	      candEta_->push_back(cand2->eta());
	      candPt_->push_back(cand2->pt());
	      candMass_->push_back(cand2->mass());

	      candPx_->push_back(cand2->px());
	      candPy_->push_back(cand2->py());
	      candPz_->push_back(cand2->pz());
	      candEnergy_->push_back(cand2->energy());

	      for(const pat::PackedCandidate &pfc : *cands) {
		if(pfc.pt()==cand2->pt() && pfc.eta()==cand2->eta() && pfc.phi()==cand2->phi()){
		  candPvAssociationQuality_->push_back(pfc.pvAssociationQuality());
		  candDXY_->push_back(pfc.dxy());
		  candDZ_->push_back(pfc.dz());
		  candDZAssociatedPV_->push_back(pfc.dzAssociatedPV());
		  //break;
		  //candVertex_->push_back(cand2.vertex());//it is a 3D vector
		}
	      }

	      candSubJetPart_->push_back(1);
	      ncandidatesjet++;
	    }
	  }
	}
      
	ncandidates_->push_back(ncandidatesjet);
      
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
	i++;
	break;//this is to save only one jet per event (the one non-matched to a lepton)
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
  
  if (!iEvent.isRealData() && mHigherOrder_) { 
    iEvent.getByToken(genEvtInfoToken,genEvtInfo);
    iEvent.getByToken(lheEvtInfoToken,lheEvtInfo);
    iEvent.getByToken(pupInfoToken,pupInfo);

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
        
    if (!iEvent.isRealData()) { 
      edm::View<PileupSummaryInfo>::const_iterator PUI;
      for(PUI = pupInfo->begin(); PUI != pupInfo->end(); ++PUI) {
	if (PUI->getBunchCrossing() == 0) {
	  npu_ = PUI->getTrueNumInteractions();
	}
      }
    }
  }

  vector<LorentzVector> vP4Gen;

  if (isMC_) { 
    iEvent.getByToken(genjetsToken,genjets);
	
    for(GenJetCollection::const_iterator i_gen = genjets->begin(); i_gen != genjets->end(); i_gen++) {

      if (i_gen->pt() > GenptMin_ && fabs(i_gen->y()) < GenetaMax_ ) {

      
	bool isHiggsMatchedGen=false;

	if(!isHiggs_) isHiggsMatchedGen=true; //if it's not the signal, we take the leading jet

	float DRmaxGen = jetRadius_;


	for(unsigned ip = 0; ip < genParticles->size(); ++ ip) {
	  const GenParticle &p = (*genParticles)[ip];
	  if (p.pdgId() == 25) {
	    if( deltaR(p.p4().eta(),p.p4().phi(),i_gen->eta(),i_gen->phi()) < DRmaxGen) isHiggsMatchedGen=true;
	  }
	}
	
	if(isHiggsMatchedGen){

	  nGenJets_++;
	  
	  GenJetpt_->push_back(i_gen->pt());
	  GenJetphi_->push_back(i_gen->phi());
	  GenJeteta_->push_back(i_gen->y());
	  GenJetenergy_->push_back(i_gen->energy());
	  GenJetmass_->push_back(i_gen->mass());

	  vP4Gen.push_back(i_gen->p4());

	  //candidates for gen jets + jet substructure
	  int nGenCandidates=0;
	  
	  std::vector<fastjet::PseudoJet>  lClusterParticles;
	  for(unsigned int ic=0; ic<i_gen->numberOfDaughters(); ic++) {
	    const reco::Candidate* gencand = i_gen->daughter(ic);

	    nGenCandidates++;
	    GenCandPt_->push_back(gencand->pt());
	    GenCandPx_->push_back(gencand->px());
	    GenCandPy_->push_back(gencand->py());
	    GenCandPz_->push_back(gencand->pz());
	    GenCandEnergy_->push_back(gencand->energy());
	    GenCandEta_->push_back(gencand->eta());
	    GenCandPhi_->push_back(gencand->phi());

	    fastjet::PseudoJet   pPart(gencand->px(),gencand->py(),gencand->pz(),gencand->energy());
	    lClusterParticles.push_back(pPart);
	  }

	  nGencandidates_->push_back(nGenCandidates);
	  
	  fClustering = new fastjet::ClusterSequenceArea(lClusterParticles, *fAKJetDef, *fAreaDefinition);
	  std::vector<fastjet::PseudoJet>  lOutJets = fClustering->inclusive_jets(20.0);
	  
	  if(lOutJets.size() == 0) {
	    delete fClustering;
	    return;
	  }
	  
	  fastjet::PseudoJet pT1JetSD = (*sd)( lOutJets[0]);
	  float iMSoftDrop = pT1JetSD.m();

	  fastjet::contrib::NormalizedMeasure normalizedMeasure(1.0,jetRadius_);
	  fastjet::contrib::Njettiness routine(fastjet::contrib::Njettiness::onepass_kt_axes,normalizedMeasure);
          float iTau1 = routine.getTau(1.,lClusterParticles);
          float iTau2 = routine.getTau(2.,lClusterParticles);
          float iTau3 = routine.getTau(3.,lClusterParticles);
	  float iTau4 = routine.getTau(4.,lClusterParticles);
          GenJetTau1_->push_back(iTau1);
          GenJetTau2_->push_back(iTau2);
          GenJetTau3_->push_back(iTau3);
	  GenJetTau4_->push_back(iTau4);
	  GenSoftDropMass_->push_back(iMSoftDrop);
	     
	  lOutJets.clear();
	  delete fClustering;
	  lClusterParticles.clear();

	  //now candidates

	}
      }
      break;//this is to save only one jet per event (the one non-matched to a lepton)
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
void nanoNtupleProducerCandidates::initialize()
{

  run_            = -1;
  evt_            = -1;
  lumi_           = -1;
  nVtx_           = -1;
  nJets_          = -1;
  nSubJets_ -> clear();
  nBSubJets_ -> clear();
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
  chm_            ->clear();
  npr_            ->clear();
  nhf_            ->clear();
  phf_            ->clear();
  elf_            ->clear();
  muf_            ->clear();
  btag_           ->clear();
  isBtag_         ->clear();
  triggerBit_     ->clear();
  triggerPre_     ->clear();

  MatchedLeptons_->clear();
  MatchedHiggs_->clear();

  WPlusLep_->clear();
  WMinusLep_->clear();

  nTriggerObjects_=-1;
  TrigObjeta_     ->clear();
  TrigObjphi_     ->clear();
  TrigObjpt_     ->clear();
  TrigObjcollection_ ->clear();

  //candidate information
  candPdgId_->clear();
  candPhi_->clear();
  candEta_->clear();
  candPt_->clear();
  candPx_->clear();
  candPy_->clear();
  candPz_->clear();
  candEnergy_->clear();
  candPuppiWeight_->clear();
  candMass_->clear();
  candPvAssociationQuality_->clear();
  candDXY_->clear();
  candDZ_->clear();
  candDZAssociatedPV_->clear();
  candSubJetPart_->clear();
  ncandidates_->clear();
  nGencandidates_->clear();

  btagSub0_->clear();
  btagSub1_->clear();
  massSub0_->clear();
  massSub1_->clear();
  ptSub0_->clear();
  ptSub1_->clear();
  etaSub0_->clear();
  etaSub1_->clear();
  phiSub0_->clear();
  phiSub1_->clear();
  flavorSub0_->clear();
  flavorSub1_->clear();
  flavorHadronSub0_->clear();
  flavorHadronSub1_->clear();

  //new ones
  massSoftDrop_->clear();
  tau1_->clear();
  tau2_->clear();
  tau3_->clear();
  tau4_->clear();

  //you need to initialize and delete all the objects of the subjets
  
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

    GenJetTau1_->clear();
    GenJetTau2_->clear();
    GenJetTau3_->clear();
    GenJetTau4_->clear();
    GenSoftDropMass_->clear();
    
    GenCandPhi_->clear();
    GenCandEta_->clear();
    GenCandPx_->clear();
    GenCandPy_->clear();
    GenCandPz_->clear();
    GenCandPt_->clear();
    GenCandEnergy_->clear();   
  }
}
//////////////////////////////////////////////////////////////////////////////////////////
nanoNtupleProducerCandidates::~nanoNtupleProducerCandidates() 
{
}

DEFINE_FWK_MODULE(nanoNtupleProducerCandidates);















