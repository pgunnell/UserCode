#ifndef nanoNtupleProducerCandidates_h
#define nanoNtupleProducerCandidates_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h" 
#include "FWCore/Common/interface/TriggerNames.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionData.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionEvaluator.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionParser.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "TTree.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "SimDataFormats/JetMatching/interface/JetFlavourInfo.h"//add
#include "SimDataFormats/JetMatching/interface/JetFlavourInfoMatching.h"//add

#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/contrib/SoftDrop.hh"

using namespace reco;

class nanoNtupleProducerCandidates : public edm::EDAnalyzer 
{
  public:
    typedef reco::Particle::LorentzVector LorentzVector;
    explicit nanoNtupleProducerCandidates(edm::ParameterSet const& cfg);
    virtual void beginJob();
    virtual void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup);
    virtual void endRun(edm::Run const& iRun, edm::EventSetup const& iSetup);
    virtual void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup);
    virtual void endJob();
    virtual ~nanoNtupleProducerCandidates();

  private:  
    virtual bool isGoodJet(const pat::Jet &jet);
    virtual bool isGoodElectron(const pat::Electron &el,const reco::Vertex &vtx,edm::Handle<pat::PackedCandidateCollection> pfcands);
    float getPFMiniIsolation(edm::Handle<pat::PackedCandidateCollection> pfcands,const reco::Candidate *cand);

    void initialize();
    //---- configurable parameters --------  
    edm::EDGetTokenT<pat::JetCollection> jetsToken;
    edm::EDGetTokenT<GenJetCollection> genjetsToken;
    edm::EDGetTokenT<pat::METCollection> metToken;
    edm::EDGetTokenT<pat::PackedCandidateCollection> candsToken;
    edm::EDGetTokenT<double> rhoToken;
    edm::EDGetTokenT<reco::VertexCollection> recVtxsToken;
    edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken;
    edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescalesToken;
    edm::EDGetTokenT<edm::View<PileupSummaryInfo> > pupInfoToken;
    edm::EDGetTokenT<GenEventInfoProduct> genEvtInfoToken;
    edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesToken;
    edm::EDGetTokenT<LHEEventProduct> lheEvtInfoToken;
    edm::EDGetTokenT<LHERunInfoProduct> runInfoToken;
    edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
    edm::EDGetTokenT<pat::MuonCollection> muonsToken;
    edm::EDGetTokenT<pat::ElectronCollection> electronsToken;

    edm::Handle<pat::MuonCollection> muons;
    edm::Handle<pat::ElectronCollection> electrons;
 
    std::string srcBtag_;
    std::vector<std::string> triggerNames_;
    double etaMax_,ptMin_,ptMinLeading_,massMin_,btagMin_,minMuPt_,minElPt_,GenetaMax_,GenptMin_;
    bool   isMC_,isPrint_,saveWeights_,debug_,mHigherOrder_,isHiggs_; 
    //---------------------------------------------------
    edm::Service<TFileService> fs_;
    TTree *outTree_; 
    TH1F *cutFlowHisto_;
    //---- TRIGGER -------------------------  
    TH1F *triggerPassHisto_,*triggerNamesHisto_;
    //---- output TREE variables ------
    //---- global event variables -----
    int   run_,evt_,nVtx_,lumi_,nJets_,nBJets_,nLeptons_,nGenJets_;
    float rho_,met_,metSig_,ht_,pvRho_,pvz_,pvndof_,pvchi2_,metGenSig_,metGen_;

    int nTriggerObjects_;
    std::vector<double> *TrigObjpt_;
    std::vector<double> *TrigObjeta_;
    std::vector<double> *TrigObjphi_;

    std::vector<string> *TrigObjcollection_;

    std::vector<bool> *triggerBit_;
    std::vector<int>  *triggerPre_;
    //---- jet variables --------------
    std::vector<bool> *MatchedLeptons_;
    std::vector<bool> *MatchedHiggs_;
    std::vector<bool> *WPlusLep_;
    std::vector<bool> *WMinusLep_;
    std::vector<bool>  *isBtag_;
    std::vector<int>   *flavor_,*nSubJets_,*nSubGenJets_,*nBSubJets_,*flavorHadron_;
    std::vector<float> *cor_, *unc_,*pt_,*eta_,*phi_,*mass_,*massSoftDrop_,*energy_,*chf_,*nhf_,*phf_,*elf_,*muf_,*btag_,*tau1_,*tau2_,*tau3_,*tau4_,*npr_,*chm_;
    std::vector<float> *btagSub0_,*btagSub1_,*massSub0_,*massSub1_,*ptSub0_,*ptSub1_,*etaSub0_,*etaSub1_,*phiSub0_,*phiSub1_;
    std::vector<int>   *flavorSub0_,*flavorSub1_,*flavorHadronSub0_,*flavorHadronSub1_;
    //---- MC variables ---------------
    int npu_,decay_;
    float genEvtWeight_,lheOriginalXWGTUP_;
    std::vector<float> *scaleWeights_;
    std::vector<float> *pdfWeights_;

    //candidate information
    std::vector<float> *candPdgId_,*candPhi_,*candEta_,*candPt_,  *candMass_,    *candPvAssociationQuality_,    *candDXY_,    *candDZ_,    *candDZAssociatedPV_,    *candSubJetPart_;
    std::vector<float> *candEnergy_, *candPx_,  *candPy_,  *candPz_, *candPuppiWeight_;
    std::vector<int> *ncandidates_;
    std::vector<int> *nGencandidates_;


    edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection> jetFlavourInfosToken_;

    //gen jets
    std::vector<float> *GenJetpt_;
    std::vector<float> *GenJetphi_;
    std::vector<float> *GenJeteta_;
    std::vector<float> *GenJetenergy_;
    std::vector<float> *GenJetmass_;
    std::vector<bool> *isBJetGen_;

    edm::Handle<pat::JetCollection> jets;
    edm::Handle<GenJetCollection> genjets;
    edm::Handle<pat::METCollection> met;
    edm::Handle<pat::PackedCandidateCollection> cands;
    edm::Handle<double> rho;
    edm::Handle<reco::VertexCollection> recVtxs;
    edm::Handle<edm::TriggerResults> triggerResults;
    edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
    edm::Handle<edm::View<PileupSummaryInfo> > pupInfo;
    edm::Handle<edm::View<reco::GenParticle> > genParticles;
    edm::Handle<GenEventInfoProduct> genEvtInfo;
    edm::Handle<LHEEventProduct> lheEvtInfo;
    edm::Handle<LHERunInfoProduct> runInfo;

    JetCorrectionUncertainty *mPFUncCHS;

    std::vector<float> *GenJetTau1_;
    std::vector<float> *GenJetTau2_;
    std::vector<float> *GenJetTau3_;
    std::vector<float> *GenJetTau4_;
    std::vector<float> *GenSoftDropMass_;

    std::vector<float> *GenCandPhi_;
    std::vector<float> *GenCandEta_;
    std::vector<float> *GenCandPx_;
    std::vector<float> *GenCandPy_;
    std::vector<float> *GenCandPz_;
    std::vector<float> *GenCandPt_;
    std::vector<float> *GenCandEnergy_;

    double jetRadius_;
    fastjet::JetAlgorithm fastjetalgorithm_;

    fastjet::Filter* fTrimmer1;
    fastjet::JetDefinition*       fAKJetDef;
    fastjet::ActiveAreaSpec*      fActiveArea;
    fastjet::AreaDefinition*      fAreaDefinition;
    fastjet::ClusterSequenceArea* fClustering;
    fastjet::contrib::SoftDrop* sd;
};





#endif
