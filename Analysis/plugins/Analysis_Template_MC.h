#ifndef AnalysisTemplate_h
#define AnalysisTemplate_h

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

float DeltaR(float eta1,float phi1,float eta2,float phi2)
{
	float deltaPhi = TMath::Abs(phi1-phi2);
	float deltaEta = eta1-eta2;
	if(deltaPhi > TMath::Pi())
	deltaPhi = TMath::TwoPi() - deltaPhi;
	return TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
}


float DeltaPhi(float phi1,float phi2)
{
	float deltaPhi = TMath::Abs(phi1-phi2);
	if(deltaPhi > TMath::Pi())
	deltaPhi = TMath::TwoPi() - deltaPhi;
	return deltaPhi;
}

class Analysis_Template_MC : public edm::EDAnalyzer
 {

  //typedef reco::Particle::LorentzVector LorentzVector;

  public:
    explicit Analysis_Template_MC(edm::ParameterSet const& cfg);
    virtual void beginJob();
    virtual void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup);
    virtual void endJob();
    virtual ~Analysis_Template_MC();


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
    // Vertices
    TH1F *num_of_Vtx;
    TH1F *num_of_Jets;
    TH1F *num_of_GenJets;
    
    int NEvents=0;
    ///Measurement Det jets
    TH1F *ptDETJet;
    TH1F *yDETJet;
    TH1F *phiDETJet;

    TH1F *ChfDETJet;
    TH1F *ElfDETJet;
    TH1F *PhfDETJet;
    TH1F *NhfDETJet;
    TH1F *MufDETJet;

    //TH1F *Tau1DETJet;
    //TH1F *Tau2DETJet;
    //TH1F *Tau3DETJet;
    //TH1F *Tau21DETJet;
    //TH1F *Tau32DETJet;
    
    TH1F *EnergyDETJet;
    TH1F *MassDETJet;
    TH1F *MassSoftDropDETJet;

 };

#endif
