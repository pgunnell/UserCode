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

vector <TLorentzVector> WbosonReconstruction(TLorentzVector electron_v4, float metPx, float metPy, float MET){

  vector <TLorentzVector> neutrino_v4;

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

  vector <TLorentzVector> wboson_v4; 

  for(unsigned int j=0; j<neutrino_v4.size(); j++){
    wboson_v4.push_back(electron_v4+neutrino_v4.at(j));
  }

  return wboson_v4;
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
	bool mElSelection;
	bool mMuSelection;
	bool mHadSelection;
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

	TH1F* JetRes;
	TH2F* JetRes2D;
	TH1F* JetResponse;

};

#endif
