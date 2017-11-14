#ifndef My_azimuthal_MC_h
#define My_azimuthal_MC_h

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
    std::string mFileName,mTreeNameGen,mTreeName,mDirName, mGlobalTag, mjettype;
    double mMinPt, mYMax;
    int    mJetID;        // looseID==1 tightID==2
    int    mprintOk;       // noPrint=0  Print=1
    bool mIsMCarlo;
    double mCrossSection;
    double mIntLumi;
    std::vector<std::string> mJECUncSrcNames;
    std::vector<std::string> mTriggers;

    edm::Service<TFileService> fs;
    TTree *mTree;
    TTree *mTreeGen;
    TFile *mInf, *mPuf;
    TDirectoryFile *mDir;

    //--------- Histogram Declaration --------------------//
    // Vertices
    TH1F *num_of_Vtx;
    TH1F *num_of_Jets;
    TH1F *num_of_GenJets;
    TH1F *num_of_SubJets;
    TH1F *num_of_SubBJets;
    TH1F *num_of_BJets;
    TH1F *num_of_Leptons;
    TH1F *num_of_GenLeptons;

    ///Measurement Gen jets
    TH1F *ptDETJet;
    TH1F *yDETJet;
    TH1F *phiDETJet;

    TH1F *ptGENJet;
    TH1F *yGENJet;
    TH1F *phiGENJet;
    TH1F *Tau1GENJet;
    TH1F *Tau2GENJet;
    TH1F *Tau3GENJet;
    TH1F *EnergyGENJet;
    TH1F *MassGENJet;
    TH1F *MassSoftDropGENJet;

    TH1F *ChfDETJet;
    TH1F *ElfDETJet;
    TH1F *PhfDETJet;
    TH1F *NhfDETJet;
    TH1F *MufDETJet;

    TH1F *Tau1DETJet;
    TH1F *Tau2DETJet;
    TH1F *Tau3DETJet;
    TH1F *Tau21DETJet;
    TH1F *Tau32DETJet;
    
    TH1F *EnergyDETJet;
    TH1F *MassDETJet;
    TH1F *MassSoftDropDETJet;
    TH1F *BtagDETJet;

    TH1F *PtDETSubJet0;
    TH1F *PtDETSubJet1;
    TH1F *EtaDETSubJet0;
    TH1F *EtaDETSubJet1;
    TH1F *PhiDETSubJet0;
    TH1F *PhiDETSubJet1;
    TH1F *MassDETSubJet0;
    TH1F *MassDETSubJet1;
    TH1F *BtagDETSubJet0;
    TH1F *BtagDETSubJet1;
    TH1F *FlavourDETSubJet0;
    TH1F *FlavourDETSubJet1;

    //Last include
    TH1F* BJetGenFlag;
    TH1F* mvaGEN;
    TH1F* metGEN;
    TH1F* metSigGEN;
    TH1F* mvaDET;
    TH1F* metDET;
    TH1F* metSigDET;
    
    TH1F* MVADesy;
    
    TH1F *ptTopParton;
    TH1F *yTopParton;
    TH1F *phiTopParton;
    TH1F *EnergyTopParton;

    TH1F *DeltaR_GenJetParton;
    TH1F *DeltaR_DetJetParton;
    TH1F *DeltaR_GenJetRecoJet;

    TH1F *ResolutionPtPartonGenJet;
    TH1F *ResolutionPtPartonDetJet;
    TH1F *ResolutionPtRecoGenJet;

    TH1F *ResolutionPhiPartonGenJet;
    TH1F *ResolutionPhiPartonDetJet;
    TH1F *ResolutionPhiRecoGenJet;

    TH1F *ResolutionMassRecoGenJet;
    TH1F *ResolutionSDMassRecoGenJet;
    TH1F *ResolutionSubJetMassRecoGenJet;

    TH1F *PtGENSubJet0;
    TH1F *PtGENSubJet1;
    TH1F *EtaGENSubJet0;
    TH1F *EtaGENSubJet1;
    TH1F *PhiGENSubJet0;
    TH1F *PhiGENSubJet1;
    TH1F *MassGENSubJet0;
    TH1F *MassGENSubJet1;
    TH1F* DeltaRGENSubJets;
    TH1F* MuGENSubJets;

    TH2F* ResponsePtRecoGenJet;
    TH2F* ResponsePhiRecoGenJet;

    TH2F* ResponsePtPartonGenJet;
    TH2F* ResponsePhiPartonGenJet;

    int numerator1;
    int numerator2;
    int numerator3;
    int denominator;

    int numeratorParton1;
    int numeratorParton2;
    int numeratorParton3;
    int denominatorParton;

 };

#endif

