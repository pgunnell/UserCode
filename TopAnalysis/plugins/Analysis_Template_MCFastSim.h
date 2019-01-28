#ifndef AnalysisTemplateFastSim_h
#define AnalysisTemplateFastSim_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TProfile.h"
#include <TMath.h>
using namespace edm;
using namespace std;

float DeltaPhiFS(float phi1, float phi2)
{
	float dPhi = TMath::Abs(phi1 - phi2);
	if(dPhi > TMath::Pi())
	{
		dPhi = TMath::TwoPi() - dPhi;
	}
	return dPhi;
}

float DeltaRFS(float eta1, float phi1, float eta2, float phi2)
{
	float deltaEta = eta1 - eta2;
	float deltaPhi = DeltaPhiFS(phi1, phi2);
	return TMath::Sqrt(deltaEta * deltaEta + deltaPhi * deltaPhi);
}

class Analysis_Template_MCFastSim : public edm::EDAnalyzer
{



public:
	//typedef reco::Particle::LorentzVector LorentzVector;
	explicit Analysis_Template_MCFastSim(edm::ParameterSet const& cfg);
	virtual void beginJob();
	virtual void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup);
	virtual void endJob();
	virtual ~Analysis_Template_MCFastSim();

private:
	//---- configurable parameters --------
	std::string mFileName,mTreeName,mDirName;

	double mMinPt, mYMax;
	int    mJetID;        // looseID==1 tightID==2
	int    mprintOk;       // noPrint=0  Print=1
	bool mIsMCarlo;
	double mCrossSection;
	double mIntLumi;
	vector<std::string> mTriggers;
	edm::Service<TFileService> fs;
	TTree *mTree;
	TFile *mInf, *mPuf;
	TDirectoryFile *mDir;

	//--------- Histogram Declaration --------------------//
	//DET JETS
	TH1F *nJets;
	TH1F *jetPt;
	TH1F *jetEta;
	TH1F *jetPhi;

	//GEN JETS
	TH1F *nGenJets;
	TH1F *GenJetPt;
	TH1F *GenJetPt_passed;
	TH1F *GenJetPt_total;
	TH1F *GenJetEta;
	TH1F *GenJetPhi;

	//GEN LQ
	TH1F *nLeptoQuarkParticles;
	TH1F *LeptoQuarkParticlePt;
	TH1F *LeptoQuarkParticlePhi;
	TH1F *LeptoQuarkParticleEta;
	TH1F *LeptoQuarkParticleMass;

	//DET MUONS
	TH1F *nMuons;
	TH1F *muonPt;
	TH1F *muonPhi;
	TH1F *muonEta;

	//GEN MUONS
	TH1F *nMuonParticles;
	TH1F *MuonParticlePt;
	TH1F *MuonParticlePhi;
	TH1F *MuonParticleEta;

	//DET TOPS
	TH1F *nTopParticles;
	TH1F *TopParticlePt;
	TH1F *TopParticleEta;
	TH1F *TopParticlePhi;
	TH1F *TopParticleMass;

	//Measurements of MLQreco and S_T
	TH1F *MLQreco_detMu;
	TH1F *MLQreco_genMu;
	TH1F *ST;

	//Added
	TH1F *MLQfullreco_detMu;
	TH1F *WbosonMass;

	TH1F *chi2;
	TH1F* LeptTopfullreco;
	TH1F* HadrTopfullreco;

	TH1F* LQHadrMassreco;
	TH1F* LQLeptMassreco;


};

#endif
