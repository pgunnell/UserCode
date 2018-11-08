#ifndef My_azimuthal_MC_h
#define My_azimuthal_MC_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"

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
  float deltaPhi = (phi1-phi2);
  float deltaEta = eta1-eta2;
  while (deltaPhi > M_PI) deltaPhi -= 2*M_PI;
  while (deltaPhi <= -M_PI) deltaPhi += 2*M_PI;
  return TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
}

float DeltaPhi(float phi1,float phi2)
{
  float deltaPhi = TMath::Abs(phi1-phi2);
  if(deltaPhi > TMath::Pi())
    deltaPhi = TMath::TwoPi() - deltaPhi;
  return deltaPhi;
}

double CorrectFactorSmear(double eta, int choice){

  double corrfactor=0;
  if(fabs(eta)<0.5) corrfactor=1.109;
  if(fabs(eta)>=0.5 && fabs(eta)<0.8) corrfactor=1.138;
  if(fabs(eta)>=1.1 && fabs(eta)<1.1) corrfactor=1.114;
  if(fabs(eta)>=1.1 && fabs(eta)<1.3) corrfactor=1.123;
  if(fabs(eta)>=1.3 && fabs(eta)<1.7) corrfactor=1.084;
  if(fabs(eta)>=1.7 && fabs(eta)<1.9) corrfactor=1.082;
  if(fabs(eta)>=1.9 && fabs(eta)<2.1) corrfactor=1.140;
  if(fabs(eta)>=2.1 && fabs(eta)<2.3) corrfactor=1.067;
  if(fabs(eta)>=2.3 && fabs(eta)<2.5) corrfactor=1.177;

  if(choice==1){//up choice
    if(fabs(eta)<0.5) corrfactor=1.109+0.008;
    if(fabs(eta)>=0.5 && fabs(eta)<0.8) corrfactor=1.138+0.013;
    if(fabs(eta)>=1.1 && fabs(eta)<1.1) corrfactor=1.114+0.013;
    if(fabs(eta)>=1.1 && fabs(eta)<1.3) corrfactor=1.123+0.024;
    if(fabs(eta)>=1.3 && fabs(eta)<1.7) corrfactor=1.084+0.011;
    if(fabs(eta)>=1.7 && fabs(eta)<1.9) corrfactor=1.082+0.035;
    if(fabs(eta)>=1.9 && fabs(eta)<2.1) corrfactor=1.140+0.047;
    if(fabs(eta)>=2.1 && fabs(eta)<2.3) corrfactor=1.067+0.053;
    if(fabs(eta)>=2.3 && fabs(eta)<2.5) corrfactor=1.177+0.041;
  }

  if(choice==2){//down choice
    if(fabs(eta)<0.5) corrfactor=1.109-0.008;
    if(fabs(eta)>=0.5 && fabs(eta)<0.8) corrfactor=1.138-0.013;
    if(fabs(eta)>=1.1 && fabs(eta)<1.1) corrfactor=1.114-0.013;
    if(fabs(eta)>=1.1 && fabs(eta)<1.3) corrfactor=1.123-0.024;
    if(fabs(eta)>=1.3 && fabs(eta)<1.7) corrfactor=1.084-0.011;
    if(fabs(eta)>=1.7 && fabs(eta)<1.9) corrfactor=1.082-0.035;
    if(fabs(eta)>=1.9 && fabs(eta)<2.1) corrfactor=1.140-0.047;
    if(fabs(eta)>=2.1 && fabs(eta)<2.3) corrfactor=1.067-0.053;
    if(fabs(eta)>=2.3 && fabs(eta)<2.5) corrfactor=1.177-0.041;
  }

  return corrfactor;
}

double PtSmeared(double PFphi,double PFeta,double PFpt,double GenPhi,double GenEta,double GenPt){

  double deltaPhiMatch1=-100;
  double deltaR1=100;
  double deltaRtest=100;
  double smear1=-100;
  double Pi=3.141592653;
  double newpt1=-100;

  deltaPhiMatch1=PFphi-GenPhi;
  if(deltaPhiMatch1<-Pi) deltaPhiMatch1=deltaPhiMatch1+2*Pi;
  if(deltaPhiMatch1>Pi) deltaPhiMatch1=deltaPhiMatch1-2*Pi;
  deltaPhiMatch1=fabs(deltaPhiMatch1);

  deltaRtest=sqrt(pow(deltaPhiMatch1,2)+pow(PFeta-GenEta,2));

  if(deltaRtest<deltaR1){
    newpt1=GenPt;
    deltaR1=deltaRtest;
  }

  if(deltaR1<=0.4){
    smear1=newpt1+CorrectFactorSmear(PFeta,0)*(PFpt-newpt1);
    return smear1;
  }

  if(deltaR1>0.4) smear1=PFpt;

  return smear1;
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

    BTagCalibration calib;
    BTagCalibrationReader reader;

    BTagCalibration calibHF;
    BTagCalibrationReader readerHF;

    //---- configurable parameters --------
    std::string mTreeNameGen,mTreeName,mDirName, mGlobalTag, mjettype,mMCtype;
    std::vector<std::string> mFileName;
    double mMinPt, mYMax;
    int    mJetID;        // looseID==1 tightID==2
    int    mprintOk;       // noPrint=0  Print=1
    bool mIsMCarlo;
    double mCrossSection;
    double mIntLumi;
    double mShift;
    std::vector<std::string> mJECUncSrcNames;
    std::vector<std::string> mTriggers;

    bool musePDFWeights;
    bool  museScaleWeights;
    
    int mnumbPDFWeights;
    int mnumbScaleWeights;


    edm::Service<TFileService> fs;
    TTree *mTree;
    TTree *mTreeGen;
    TFile *mInf, *mPuf;
    TDirectoryFile *mDir;

    //--------- Histogram Declaration --------------------//
    // Vertices
    TH1F *num_of_Vtx;
    TH1F *PileUpInteractions;
    TH1F *num_of_Jets;
    TH1F *num_of_GenJets;
    TH1F *num_of_SubJets;
    TH1F *num_of_SubBJets;
    TH1F *num_of_BJets;
    TH1F *num_of_Leptons;
    TH1F *num_of_GenLeptons;

    TH1F *num_of_VtxSelection;

    ///Measurement Gen jets
    TH1F *ptDETJet;
    TH1F *yDETJet;
    TH1F *phiDETJet;

    TH1F *ptGENJet;
    TH1F *ptGENJetLeading;
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
    TH1F *Tau31DETJet;
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
    TH1F* WJetGenFlag;
    TH1F* mvaGEN;
    TH1F* metGEN;
    TH1F* metSigGEN;
    TH1F* mvaDET;
    TH1F* metDET;
    TH1F* metSigDET;

    TH1F* Tau32GENJet;
    TH1F* Tau31GENJet;
    
    TH1F *ptTopParton;
    TH1F *yTopParton;
    TH1F *phiTopParton;
    TH1F *EnergyTopParton;

    TH1F *ptTopParton_Leading;
    TH1F *ptTopParton_SubLeading;
    TH1F *ptTopParton_Leading_AfterMatching;
    TH1F *ptTopParton_SubLeading_AfterMatching;
    TH1F *ptTopParton_Leading_AfterMatchingTagging;
    TH1F *ptTopParton_SubLeading_AfterMatchingTagging;

    TH1F *ptTopPartonLeading;
    TH1F *ptTopPartonSubLeading;
    TH1F *ptTopPartonMatchedLeading;
    TH1F *ptTopPartonMatchedSubLeading;
    TH1F *ptTopPartonMatchedWBGenJetLeading;
    TH1F *ptTopPartonMatchedWBGenJetSubLeading;
    TH1F *ptTopPartonMatchedWSubjetMassBGenJetLeading;
    TH1F *ptTopPartonMatchedWSubjetMassBGenJetSubLeading;
    
    TH1F *DeltaR_GenJetPartonLead;
    TH1F *DeltaR_GenJetPartonSubLead;
    TH1F *DeltaR_GenJetW;
    TH1F *DeltaR_DetJetParton;
    TH1F *DeltaR_GenJetRecoJet;

    TH1F *ResolutionPtPartonGenJet;
    TH1F *ResolutionPtPartonDetJet;
    TH1F *ResolutionPtRecoGenJet;

    TH1F *ResolutionPhiPartonGenJet;
    TH1F *ResolutionPhiPartonDetJet;
    TH1F *ResolutionPhiRecoGenJet;
    TH1F *ResolutionRapidityRecoGenJet;
    TH1F *ResolutionDeltaPhiRecoGenJet;

    TH1F *MassGENSubJet0SubLeadingAfterSelection;
    TH1F *MassGENSubJet0LeadingAfterSelection;
    TH1F *MassGENJetLeadingAfterSelection;
    TH1F *MassGENJetSubLeadingAfterSelection;

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
    TH2F* ResponseDeltaPhiRecoGenJet;
    TH2F* ResponseForTUnfoldDeltaPhiRecoGenJet;
    TH2F* ResponseRapidityRecoGenJet;

    TH2F* ResponsePtPartonGenJet;
    TH2F* ResponsePhiPartonGenJet;

    int NEvents=0;
    int Neventsdetselection=0;
    int Neventsgenselection=0;
    int Neventspartonselection=0;
    int NEventsQCD=0;
    int dijetsafter=0;
    int dijetsbefore=0;

    TH1F* MVADesy;

    //Prepreselection

    TH1F* PtDETLeadingPrepreselectionSubJet0;
    TH1F* PtDETLeadingPrepreselectionSubJet1;
    TH1F* EtaDETLeadingPrepreselectionSubJet0;
    TH1F* EtaDETLeadingPrepreselectionSubJet1;
    TH1F* MassDETLeadingPrepreselectionSubJet0;
    TH1F* MassDETLeadingPrepreselectionSubJet1;
    TH1F* PhiDETLeadingPrepreselectionSubJet0;
    TH1F* PhiDETLeadingPrepreselectionSubJet1;
    TH1F* DeltaRDETLeadingPrepreselectionSubJets;

    TH1F* Tau1DETLeadingJetPrepreselection;
    TH1F* Tau2DETLeadingJetPrepreselection;
    TH1F* Tau3DETLeadingJetPrepreselection;
    TH1F* Tau31DETLeadingJetPrepreselection;
    TH1F* Tau32DETLeadingJetPrepreselection;

    TH1F* PtDETSubLeadingPrepreselectionSubJet0;
    TH1F* PtDETSubLeadingPrepreselectionSubJet1;
    TH1F* EtaDETSubLeadingPrepreselectionSubJet0;
    TH1F* EtaDETSubLeadingPrepreselectionSubJet1;
    TH1F* MassDETSubLeadingPrepreselectionSubJet0;
    TH1F* MassDETSubLeadingPrepreselectionSubJet1;
    TH1F* PhiDETSubLeadingPrepreselectionSubJet0;
    TH1F* PhiDETSubLeadingPrepreselectionSubJet1;
    TH1F* DeltaRDETSubLeadingPrepreselectionSubJets;

    TH1F* Tau1DETSubLeadingJetPrepreselection;
    TH1F* Tau2DETSubLeadingJetPrepreselection;
    TH1F* Tau3DETSubLeadingJetPrepreselection;
    TH1F* Tau31DETSubLeadingJetPrepreselection;
    TH1F* Tau32DETSubLeadingJetPrepreselection;

    TH1F* PtDETLeadingPrepreselection;
    TH1F* PtDETSubLeadingPrepreselection;
    TH1F* EtaDETLeadingPrepreselection;
    TH1F* EtaDETSubLeadingPrepreselection;
    TH1F* MassDETLeadingPrepreselection;
    TH1F* MassDETSubLeadingPrepreselection;
    TH1F* PhiDETLeadingPrepreselection;
    TH1F* PhiDETSubLeadingPrepreselection;
    TH1F* SDMassDETLeadingPrepreselection;
    TH1F* SDMassDETSubLeadingPrepreselection;

    //preselection

    TH1F* PtDETLeadingPreselectionSubJet0;
    TH1F* PtDETLeadingPreselectionSubJet1;
    TH1F* EtaDETLeadingPreselectionSubJet0;
    TH1F* EtaDETLeadingPreselectionSubJet1;
    TH1F* MassDETLeadingPreselectionSubJet0;
    TH1F* MassDETLeadingPreselectionSubJet1;
    TH1F* PhiDETLeadingPreselectionSubJet0;
    TH1F* PhiDETLeadingPreselectionSubJet1;
    TH1F* DeltaRDETLeadingPreselectionSubJets;

    TH1F* Tau1DETLeadingJetPreselection;
    TH1F* Tau2DETLeadingJetPreselection;
    TH1F* Tau3DETLeadingJetPreselection;
    TH1F* Tau31DETLeadingJetPreselection;
    TH1F* Tau32DETLeadingJetPreselection;

    TH1F* PtDETSubLeadingPreselectionSubJet0;
    TH1F* PtDETSubLeadingPreselectionSubJet1;
    TH1F* EtaDETSubLeadingPreselectionSubJet0;
    TH1F* EtaDETSubLeadingPreselectionSubJet1;
    TH1F* MassDETSubLeadingPreselectionSubJet0;
    TH1F* MassDETSubLeadingPreselectionSubJet1;
    TH1F* PhiDETSubLeadingPreselectionSubJet0;
    TH1F* PhiDETSubLeadingPreselectionSubJet1;
    TH1F* DeltaRDETSubLeadingPreselectionSubJets;

    TH1F* Tau1DETSubLeadingJetPreselection;
    TH1F* Tau2DETSubLeadingJetPreselection;
    TH1F* Tau3DETSubLeadingJetPreselection;
    TH1F* Tau31DETSubLeadingJetPreselection;
    TH1F* Tau32DETSubLeadingJetPreselection;

    TH1F* PtDETLeadingPreselection;
    TH1F* PtDETSubLeadingPreselection;
    TH1F* EtaDETLeadingPreselection;
    TH1F* EtaDETSubLeadingPreselection;
    TH1F* MassDETLeadingPreselection;
    TH1F* MassDETSubLeadingPreselection;
    TH1F* PhiDETLeadingPreselection;
    TH1F* PhiDETSubLeadingPreselection;
    TH1F* SDMassDETLeadingPreselection;
    TH1F* SDMassDETSubLeadingPreselection;

    //Selection
    TH1F* PtDETLeadingSelectionSubJet0;
    TH1F* PtDETLeadingSelectionSubJet1;
    TH1F* EtaDETLeadingSelectionSubJet0;
    TH1F* EtaDETLeadingSelectionSubJet1;
    TH1F* MassDETLeadingSelectionSubJet0;
    TH1F* MassDETLeadingSelectionSubJet1;
    TH1F* PhiDETLeadingSelectionSubJet0;
    TH1F* PhiDETLeadingSelectionSubJet1;
    TH1F* DeltaRDETLeadingSelectionSubJets;

    TH1F* Tau1DETLeadingJetSelection;
    TH1F* Tau2DETLeadingJetSelection;
    TH1F* Tau3DETLeadingJetSelection;
    TH1F* Tau31DETLeadingJetSelection;
    TH1F* Tau32DETLeadingJetSelection;

    TH1F* PtDETSubLeadingSelectionSubJet0;
    TH1F* PtDETSubLeadingSelectionSubJet1;
    TH1F* EtaDETSubLeadingSelectionSubJet0;
    TH1F* EtaDETSubLeadingSelectionSubJet1;
    TH1F* MassDETSubLeadingSelectionSubJet0;
    TH1F* MassDETSubLeadingSelectionSubJet1;
    TH1F* PhiDETSubLeadingSelectionSubJet0;
    TH1F* PhiDETSubLeadingSelectionSubJet1;
    TH1F* DeltaRDETSubLeadingSelectionSubJets;

    TH1F* Tau1DETSubLeadingJetSelection;
    TH1F* Tau2DETSubLeadingJetSelection;
    TH1F* Tau3DETSubLeadingJetSelection;
    TH1F* Tau31DETSubLeadingJetSelection;
    TH1F* Tau32DETSubLeadingJetSelection;

    TH1F* PtDETLeadingSelection;
    TH1F* PtDETSubLeadingSelection;
    TH1F* EtaDETLeadingSelection;
    TH1F* EtaDETSubLeadingSelection;
    TH1F* MassDETLeadingSelection;
    TH1F* MassDETSubLeadingSelection;
    TH1F* PhiDETLeadingSelection;
    TH1F* PhiDETSubLeadingSelection;
    TH1F* SDMassDETLeadingSelection;
    TH1F* SDMassDETSubLeadingSelection;

    TH1F* MassDETLeadingSelection_bins[6];
    TH1F* PtDETLeadingSelection_bins[5];
    TH1F* PtDETSubLeadingSelection_bins[5];

    TH1F* DeltaPhiPARTONUnmatched;
    TH1F* DeltaPhiDETSelection;
    TH1F* DeltaPhiDETSelectionDoubleBins;
    TH1F* DeltaPhiDETSelectionForToyFits;
    TH1F* DeltaPhiPARTONSelection;

    TH1F* DeltaPhiDETSelectionFake;    TH1F* DeltaPhiDETSelectionMiss;

    TH1F* PtGENLeadingSelection_bins[5];
    TH1F* PtGENSubLeadingSelection_bins[5];

    TH1F* PtPARTONLeadingSelection_bins[5];
    TH1F* PtPARTONSubLeadingSelection_bins[5];

    TH1F* PtGENLeadingSelectionForDoubleBins_bins[5];
    TH1F* PtGENLeadingSelectionForToyFits_bins[5];

    TH1F* PtDETLeadingSelectionForDoubleBins_bins[5];
    TH1F* PtDETLeadingSelectionForToyFits_bins[5];

    TH1F* PtGENSubLeadingSelectionForDoubleBins_bins[5];
    TH1F* PtGENSubLeadingSelectionForToyFits_bins[5];

    TH1F* PtDETSubLeadingSelectionForDoubleBins_bins[5];
    TH1F* PtDETSubLeadingSelectionForToyFits_bins[5];

    TH1F* DeltaPhiGENSelectionDoubleBins;
    TH1F* DeltaPhiGENSelectionForToyFits;
    TH1F* DeltaPhiGENSelection;
    TH1F* DeltaPhiMatchedSelection;

    TH1F* PtMatchedGenJet;
    TH1F* PtMatchedGenJetAndBJet;
    TH1F* PtMatchedGenJetAndWJet;
    TH1F* PtMatchedGenJetAndBWJet;
    TH1F* PtMatchedGenJetAndBWBosonJet;
    TH1F* PtMatchedGenJetAndWBosonJet;
    TH1F* PtMatchedGenJetAndBWSubjetMassJet;
    TH1F* PtMatchedGenJetAndWSubjetMassJet;
    TH1F* PtMatchedGenJetAndBWSubjetMassHighJet;
    TH1F* PtMatchedGenJetAndWSubjetMassHighJet;

    TH1F* PtMatchedRecoJet;
    TH1F* PtMatchedRecoJetSubBtags;

    TH1F* PtMatchedGenJetAndBWJetLeading;
    TH1F* PtMatchedGenJetAndBWJetLeadingReco;

    TH1F* PtMatchedGenJetAndBWJetSubLeading;
    TH1F* PtMatchedGenJetAndBWJetSubLeadingReco;

    ifstream file_input;
    ifstream file_PUNominal;

    TH1F* TaggedAndBJetDETSubJet0NUM;
    TH1F* TaggedAndBJetDETSubJet1NUM;

    TH1F* BJetAndTaggedDETSubJet0NUM;
    TH1F* BJetAndTaggedDETSubJet1NUM;

    TH1F* TaggedAndBJetDETSubJet0DEN;
    TH1F* TaggedAndBJetDETSubJet1DEN;

    TH1F* BJetAndTaggedDETSubJet0DEN;
    TH1F* BJetAndTaggedDETSubJet1DEN;

    TH2F* ResponseLeadingJetPt_bins[5];
    TH2F* ResponseSubLeadingJetPt_bins[5];

    TH2F* ResponseLeadingJetPtForTUnfold_bins[5];
    TH2F* ResponseSubLeadingJetPtForTUnfold_bins[5];

    TH1F* SDMassDETLeadingSelection_DeltaPhibins[6];
    TH1F* SDMassDETSubLeadingSelection_DeltaPhibins[6];

    TH1F* SDMassDETLeadingSelection_Ptbins[6];
    TH1F* SDMassDETSubLeadingSelection_Ptbins[6];

    TH1F* SubJet1MassGENJetLeadingAfterSelection;
    TH1F* SubJet2MassGENJetLeadingAfterSelection;

    TH1F* SubJet1MassGENJetSubLeadingAfterSelection;
    TH1F* SubJet2MassGENJetSubLeadingAfterSelection;

    TH1F* MassGENJetLead_1Selection;
    TH1F* MassGENJetLead_2Selection;

    TH1F* MassGENJetSubLead_1Selection;
    TH1F* MassGENJetSubLead_2Selection;

    TH1F* PtGENJetLead_1Selection;
    TH1F* PtGENJetLead_2Selection;

    TH1F* PtGENJetSubLead_1Selection;
    TH1F* PtGENJetSubLead_2Selection;

    TH1F* PtDETLeadingPrepreselectionBTrue;
    TH1F* PtDETSubLeadingPrepreselectionBTrue;

    TH1F* FlavourDETSelectionLeadingJet;
    TH1F* FlavourDETSelectionSubLeadingJet;
    TH1F* FlavourDETPrepreselectionLeadingJet;
    TH1F* FlavourDETPrepreselectionSubLeadingJet;

    TH1F* PtLeadingGenJetMatched;
    TH1F* PtLeadingGenJet;

    TH1F* PtSubLeadingGenJetMatched;
    TH1F* PtSubLeadingGenJet;

    TH1F* MassGENJetLead_TopParton;
    TH1F* MassGENJetLead_SubJet1TopParton;
    
    TH1F* MassGENJetSubLead_TopParton;
    TH1F* MassGENJetSubLead_SubJet1TopParton;

    TH1F* JetBTaggedSubjets1Mass_FirstSubjet;
    TH1F* JetBTaggedSubjets0Mass_FirstSubjet;

    TH1F* JetBTaggedSubjets1Mass_SecondSubjet;
    TH1F* JetBTaggedSubjets0Mass_SecondSubjet;

    TH1F* ptTopParton_Leading_SubSubLeading;
    TH1F* ptTopParton_SubLeading_SubSubLeading;

    vector<double> ReweightPU;
 };

#endif
