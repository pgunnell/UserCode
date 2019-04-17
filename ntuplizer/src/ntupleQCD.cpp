/*
 * ntupleQCD.cpp
 *
 *  Created on: 24 Aug 2016
 *      Author: jkiesele
 */

#include "interface/ntupleQCD.h"



void ntupleQCD::analyze(size_t childid /* this info can be used for printouts */){

	/*
	 * This skeleton analyser runs directly on the Delphes output.
	 * It can be used to create histograms directly or a skim.
	 * If a skim is created, a new input configuration will be written automatically
	 * and stored in the output directory together with the ntuples.
	 * The skim can contain delphes objects again or can be flat. This is up
	 * to the user.
	 * Examples for both are given here.
	 *
	 * The same skeleton can be used to read the skim. Please refer to the comments
	 * marked with "==SKIM=="
	 *
	 * These parts are commented, since the code is supposed to work as an example without
	 * modifications on Delphes output directly.
	 */



	/*
	 * Define the branches that are to be considered for the analysis
	 * This branch handler (notice the "d")
	 * is used to run directly in Delphes output.
	 * For skimmed ntuples, see below
	 */

	/*
	 * Other branches might be the following
	 * (for a full list, please inspect the Delphes sample root file with root)
	 * For the Delphes class description, see $DELPHES_PATH/classes/DelphesClasses.h
	 */
	//

	d_ana::dBranchHandler<Electron> elecs(tree(),"Electron");
	d_ana::dBranchHandler<HepMCEvent>  event(tree(),"Event");
	d_ana::dBranchHandler<GenParticle> genpart(tree(),"Particle");
	d_ana::dBranchHandler<Jet>         jet(tree(),"Jet");
	d_ana::dBranchHandler<Jet>         genjet(tree(),"GenJet");
	d_ana::dBranchHandler<Muon>        muontight(tree(),"MuonTight");
	d_ana::dBranchHandler<Muon>        muonloose(tree(),"MuonLoose");
	d_ana::dBranchHandler<MissingET>   met(tree(),"MissingET");
	d_ana::dBranchHandler<MissingET>   genmet(tree(),"GenMissingET");
	d_ana::dBranchHandler<ScalarHT>   scalarHt(tree(),"ScalarHT");
	d_ana::dBranchHandler<Vertex>   vertex(tree(),"Vertex");

	/* ==SKIM==
	 *
	 * If a skim of the Delphes outout was created in a way indicated
	 * further below, use the tBranchHandler (please notive the "t")
	 * to access vectors of objects...
	 *
	 */
	// d_ana::tBranchHandler<std::vector<Electron> > electrons(tree(),"Electrons");

	/*==SKIM==
	 *
	 * Or an object directly
	 *
	 */
	//d_ana::tBranchHandler<MissingET> met(tree(),"MET");



	/*
	 * Always use this function to add a new histogram (can also be 2D)!
	 * Histograms created this way are automatically added to the output file
	 */
	//TH1* histo=addPlot(new TH1D("histoname1","histotitle1",100,0,100),"p_{T} [GeV]","N_{e}"); //commented PG


	/*
	 * If (optionally) a skim or a flat ntuple is to be created, please use the following function to initialize
	 * the tree.
	 * The output files will be written automatically, and a config file will be created.
	 */
	TTree* myskim = addTree();
	/*
	 * Add a simple branch to the skim
	 */
	std::vector<float> elecPt; 	std::vector<float> elecPhi; 	std::vector<float> elecEta; 
	std::vector<float> muonPt; 	std::vector<float> muonPhi; 	std::vector<float> muonEta; 
	myskim->Branch("elecPt","vector<float>" ,&elecPt);
	myskim->Branch("elecEta","vector<float>" ,&elecEta);
	myskim->Branch("elecPhi","vector<float>" ,&elecPhi);

	myskim->Branch("muonPt","vector<float>" ,&muonPt);
	myskim->Branch("muonEta","vector<float>" ,&muonEta);
	myskim->Branch("muonPhi","vector<float>" ,&muonPhi);

	std::vector<float> elecIsoVar; 	std::vector<float> elecIsoVarRhoCorr; std::vector<float> elecSumPtCharged; std::vector<float> elecSumPtNeutral; std::vector<float> elecSumPtChargedPU; std::vector<float> elecSumPt; std::vector<int> elecCharge;

	myskim->Branch("elecIsoVar","vector<float>" ,&elecIsoVar);
	myskim->Branch("elecIsoVarRhoCorr","vector<float>" ,&elecIsoVarRhoCorr);
	myskim->Branch("elecSumPtCharged","vector<float>" ,&elecSumPtCharged);
	myskim->Branch("elecSumPtNeutral","vector<float>" ,&elecSumPtNeutral);
	myskim->Branch("elecSumPtChargedPU","vector<float>" ,&elecSumPtChargedPU);
	myskim->Branch("elecSumPt","vector<float>" ,&elecSumPt);
	myskim->Branch("elecCharge","vector<int>" ,&elecCharge);


	std::vector<float> muonIsoVar; 	std::vector<float> muonIsoVarRhoCorr; std::vector<float> muonSumPtCharged; std::vector<float> muonSumPtNeutral; std::vector<float> muonSumPtChargedPU; std::vector<float> muonSumPt; std::vector<int> muonCharge;

	myskim->Branch("muonIsoVar","vector<float>" ,&muonIsoVar);
	myskim->Branch("muonIsoVarRhoCorr","vector<float>" ,&muonIsoVarRhoCorr);
	myskim->Branch("muonSumPtCharged","vector<float>" ,&muonSumPtCharged);
	myskim->Branch("muonSumPtNeutral","vector<float>" ,&muonSumPtNeutral);
	myskim->Branch("muonSumPtChargedPU","vector<float>" ,&muonSumPtChargedPU);
	myskim->Branch("muonSumPt","vector<float>" ,&muonSumPt);
	myskim->Branch("muonCharge","vector<int>" ,&muonCharge);

	//jets

        std::vector<float> pt_;         std::vector<float> btag_ ;  std::vector<float> btagAlgo_ ;         std::vector<float> eta_;           std::vector<float> phi_;         std::vector<float> mass_;      std::vector<float> jetflavour_ ; std::vector<float> jetTauTag_ ;

	myskim->Branch("jetPt"                ,"vector<float>"     ,&pt_);
	myskim->Branch("jetBtag"              ,"vector<float>"     ,&btag_);  
	myskim->Branch("jetBtagAlgo"              ,"vector<float>"     ,&btagAlgo_);  
	myskim->Branch("jetEta"               ,"vector<float>"     ,&eta_);
	myskim->Branch("jetPhi"               ,"vector<float>"     ,&phi_);
	myskim->Branch("jetMass"              ,"vector<float>"     ,&mass_);
	myskim->Branch("jetFlavorHadron" ,"vector<float>" ,&jetflavour_);
	myskim->Branch("jetTauTag" ,"vector<float>" ,&jetTauTag_);

        std::vector<float> GenJetpt_;    std::vector<float> GenJeteta_;           std::vector<float> GenJetphi_;      std::vector<float> GenJetMass_;     

	myskim->Branch("GenJetPt"                ,"vector<float>"     ,&GenJetpt_);
	myskim->Branch("GenJetEta"               ,"vector<float>"     ,&GenJeteta_);
	myskim->Branch("GenJetPhi"               ,"vector<float>"     ,&GenJetphi_); 
	myskim->Branch("GenJetMass"               ,"vector<float>"     ,&GenJetMass_); 

	std::vector<float> GenParticlept_;    std::vector<float> GenParticleeta_;           std::vector<float> GenParticlephi_;     std::vector<int> GenParticleid_;     std::vector<float> GenParticleMass_;     

	myskim->Branch("GenParticlePt"                ,"vector<float>"     ,&GenParticlept_);
	myskim->Branch("GenParticleEta"               ,"vector<float>"     ,&GenParticleeta_);
	myskim->Branch("GenParticlePhi"               ,"vector<float>"     ,&GenParticlephi_);
	myskim->Branch("GenParticleId"               ,"vector<int>"     ,&GenParticleid_);
	myskim->Branch("GenParticleMass"               ,"vector<float>"     ,&GenParticleMass_);

	std::vector<float> MuonParticlept_;    std::vector<float> MuonParticleeta_;           std::vector<float> MuonParticlephi_;     std::vector<int> MuonParticleid_;     std::vector<float> MuonParticleMass_;     std::vector<int> MuonParticleStatus_;

	myskim->Branch("MuonParticlePt"                ,"vector<float>"     ,&MuonParticlept_);
	myskim->Branch("MuonParticleEta"               ,"vector<float>"     ,&MuonParticleeta_);
	myskim->Branch("MuonParticlePhi"               ,"vector<float>"     ,&MuonParticlephi_);
	myskim->Branch("MuonParticleId"               ,"vector<int>"     ,&MuonParticleid_);
	myskim->Branch("MuonParticleMass"               ,"vector<float>"     ,&MuonParticleMass_);
	myskim->Branch("MuonParticleStatus"               ,"vector<int>"     ,&MuonParticleStatus_);

	std::vector<float> LeptoQuarkParticlept_;    std::vector<float> LeptoQuarkParticleeta_;           std::vector<float> LeptoQuarkParticlephi_;     std::vector<int> LeptoQuarkParticleid_;     std::vector<float> LeptoQuarkParticleMass_;     

	myskim->Branch("LeptoQuarkParticlePt"                ,"vector<float>"     ,&LeptoQuarkParticlept_);
	myskim->Branch("LeptoQuarkParticleEta"               ,"vector<float>"     ,&LeptoQuarkParticleeta_);
	myskim->Branch("LeptoQuarkParticlePhi"               ,"vector<float>"     ,&LeptoQuarkParticlephi_);
	myskim->Branch("LeptoQuarkParticleId"               ,"vector<int>"     ,&LeptoQuarkParticleid_);
	myskim->Branch("LeptoQuarkParticleMass"               ,"vector<float>"     ,&LeptoQuarkParticleMass_);

	std::vector<float> TopParticlept_;    std::vector<float> TopParticleeta_;           std::vector<float> TopParticlephi_;     std::vector<int> TopParticleid_;     std::vector<float> TopParticleMass_;     

	myskim->Branch("TopParticlePt"                ,"vector<float>"     ,&TopParticlept_);
	myskim->Branch("TopParticleEta"               ,"vector<float>"     ,&TopParticleeta_);
	myskim->Branch("TopParticlePhi"               ,"vector<float>"     ,&TopParticlephi_);
	myskim->Branch("TopParticleId"               ,"vector<int>"     ,&TopParticleid_);
	myskim->Branch("TopParticleMass"               ,"vector<float>"     ,&TopParticleMass_);

	//MET
	float METvalue_;
	float METeta_;
	float METphi_;

	myskim->Branch("MET"                   ,&METvalue_,"METvalue_/F");
	myskim->Branch("METeta"                ,&METeta_,"METeta_/F" );
	myskim->Branch("METphi"                ,&METphi_,"METphi_/F" );

	float GenMETvalue_;
	float GenMETeta_;
	float GenMETphi_;

	myskim->Branch("GenMET"                   ,&GenMETvalue_,"GenMETvalue_/F");
	myskim->Branch("GenMETeta"                ,&GenMETeta_,"GenMETeta_/F" );
	myskim->Branch("GenMETphi"                ,&GenMETphi_,"GenMETphi_/F" );

	//Scalar HT
	float ScalarHTvalue_;
	myskim->Branch("ScalarHT"                ,&ScalarHTvalue_,"ScalarHTvalue_/F");

	int NVertex_; std::vector<float> VertexX_;    std::vector<float> VertexY_;           std::vector<float> VertexZ_;
	myskim->Branch("nVertex"                ,&NVertex_,"nVertex_/I");
	myskim->Branch("VertexX"                ,"vector<float>"  ,&VertexX_);
	myskim->Branch("VertexY"                ,"vector<float>"  ,&VertexY_);
	myskim->Branch("VertexZ"                ,"vector<float>"  ,&VertexZ_);

	/*
	 * Or store a vector of objects (also possible to store only one object)
	 */
	//std::vector<Electron> skimmedelecs; //commented PG
	//myskim->Branch("Electrons",&skimmedelecs); //commented PG

	size_t nevents=tree()->entries();
	if(isTestMode())
		nevents/=100;
	for(size_t eventno=0;eventno<nevents;eventno++){
	
	  /*
	   * The following two lines report the status and set the event link
	   * Do not remove!
	   */
	  reportStatus(eventno,nevents);
	  tree()->setEntry(eventno);

	  /*
	   * Begin the event-by-event analysis
	   */

	  /*
	   * Or to fill the skim
	   */

	  NVertex_=0;

	  for(size_t i=0;i<vertex.size();i++){
	    NVertex_++;
	    VertexX_.push_back(vertex.at(i)->X);
	    VertexY_.push_back(vertex.at(i)->Y);
	    VertexZ_.push_back(vertex.at(i)->Z);
	  }

	  for(size_t i=0;i<elecs.size();i++){
	    //flat info
	    if (fabs(elecs.at(i)->Eta) > 2.8 || (fabs(elecs.at(i)->Eta) > 1.4442 && fabs(elecs.at(i)->Eta) < 1.5660)) continue;
	    if (elecs.at(i)->PT < 30.) continue;
	    if (elecs.at(i)->IsolationVar > 0.14) continue;
	    //added by PG
	    elecPt.push_back(elecs.at(i)->PT);
	    elecEta.push_back(elecs.at(i)->Eta); 
	    elecPhi.push_back(elecs.at(i)->Phi); 
	    //isolation variables
	    elecIsoVar.push_back(elecs.at(i)->IsolationVar); 
	    elecIsoVarRhoCorr.push_back(elecs.at(i)->IsolationVarRhoCorr); 
	    elecSumPtCharged.push_back(elecs.at(i)->SumPtCharged); 
	    elecSumPtNeutral.push_back(elecs.at(i)->SumPtNeutral); 
	    elecSumPtChargedPU.push_back(elecs.at(i)->SumPtChargedPU); 
	    elecSumPt.push_back(elecs.at(i)->SumPt); 
	    elecCharge.push_back(elecs.at(i)->Charge); 
	  }

	  for(size_t i=0;i<muontight.size();i++){
	    //flat info
	    //added by PG
	    if (fabs(muontight.at(i)->Eta) > 2.8) continue;
            if (muontight.at(i)->PT < 30.) continue;
            if (muontight.at(i)->IsolationVar > 0.14) continue;
	    muonPt.push_back(muontight.at(i)->PT);
	    muonEta.push_back(muontight.at(i)->Eta); 
	    muonPhi.push_back(muontight.at(i)->Phi); 
	    //isolation variables
	    muonIsoVar.push_back(muontight.at(i)->IsolationVar); 
	    muonIsoVarRhoCorr.push_back(muontight.at(i)->IsolationVarRhoCorr); 
	    muonSumPtCharged.push_back(muontight.at(i)->SumPtCharged); 
	    muonSumPtNeutral.push_back(muontight.at(i)->SumPtNeutral); 
	    muonSumPtChargedPU.push_back(muontight.at(i)->SumPtChargedPU); 
	    muonSumPt.push_back(muontight.at(i)->SumPt); 
	    muonCharge.push_back(muontight.at(i)->Charge); 
	  }
	
	  for(size_t i=0;i<jet.size();i++){
	    bool overlaps = false;
	    for (size_t j = 0; j < elecs.size(); j++) {
	      if (fabs(jet.at(i)->PT-elecs.at(j)->PT) < 0.01*elecs.at(j)->PT && elecs.at(j)->P4().DeltaR(jet.at(i)->P4()) < 0.01) {
		overlaps = true;
		break;
	      }
	    }
	    if (overlaps) continue;
	    for (size_t j = 0; j < muontight.size(); j++) {
	      if (fabs(jet.at(i)->PT-muontight.at(j)->PT) < 0.01*muontight.at(j)->PT && muontight.at(j)->P4().DeltaR(jet.at(i)->P4()) < 0.01) {
		overlaps = true;
		break;
	      }
	    }
	    if (overlaps) continue;

	    pt_.push_back(jet.at(i)->PT);
	    eta_.push_back(jet.at(i)->Eta); 
	    phi_.push_back(jet.at(i)->Phi); 
	    mass_.push_back(jet.at(i)->Mass); 
	    jetTauTag_.push_back(jet.at(i)->TauTag); 
	    jetflavour_.push_back(jet.at(i)->Flavor);
	    btag_.push_back(jet.at(i)->BTag);
	    btagAlgo_.push_back(jet.at(i)->BTagAlgo);
	  }

	  for(size_t i=0;i<genjet.size();i++){
	    //flat info
	    //added by PG
	    GenJetpt_.push_back(genjet.at(i)->PT);
	    GenJeteta_.push_back(genjet.at(i)->Eta); 
	    GenJetphi_.push_back(genjet.at(i)->Phi); 
	    GenJetMass_.push_back(genjet.at(i)->Mass); 
	  }

	  for(size_t i=0;i<genpart.size();i++){
	    //flat info
	    if(abs(genpart.at(i)->PID)<=25 && genpart.at(i)->Status<=23 && genpart.at(i)->Status>=22){
	      GenParticlept_.push_back(genpart.at(i)->PT);
	      GenParticleeta_.push_back(genpart.at(i)->Eta); 
	      GenParticlephi_.push_back(genpart.at(i)->Phi); 
	      GenParticleid_.push_back(genpart.at(i)->PID); 
	      GenParticleMass_.push_back(genpart.at(i)->Mass); 

	      //std::cout<<"particle "<<genpart.at(i)->Status<<" "<<genpart.at(i)->PID<<std::endl;
	    
	    }
	    
	    if(abs(genpart.at(i)->PID)==13){
	      //std::cout<<"particle "<<genpart.at(i)->Status<<" "<<genpart.at(i)->Mass<<std::endl;
	      if(genpart.at(i)->Status==23 || genpart.at(i)->Status==51){
		MuonParticlept_.push_back(genpart.at(i)->PT);
		MuonParticleeta_.push_back(genpart.at(i)->Eta); 
		MuonParticlephi_.push_back(genpart.at(i)->Phi); 
		MuonParticleid_.push_back(genpart.at(i)->PID); 
		MuonParticleMass_.push_back(genpart.at(i)->Mass); 
		MuonParticleStatus_.push_back(genpart.at(i)->Status); 
	      }
	    }

	    if((abs(genpart.at(i)->PID)>=9000000 && abs(genpart.at(i)->PID)<=9000009 && genpart.at(i)->Status==22) || (abs(genpart.at(i)->PID)==42 && genpart.at(i)->Status==22)){
	      LeptoQuarkParticlept_.push_back(genpart.at(i)->PT);
	      LeptoQuarkParticleeta_.push_back(genpart.at(i)->Eta); 
	      LeptoQuarkParticlephi_.push_back(genpart.at(i)->Phi); 
	      LeptoQuarkParticleid_.push_back(genpart.at(i)->PID); 
	      LeptoQuarkParticleMass_.push_back(genpart.at(i)->Mass); 
	    }

	    if(abs(genpart.at(i)->PID)==6 && genpart.at(i)->Status<=24 && genpart.at(i)->Status>=22){
	      TopParticlept_.push_back(genpart.at(i)->PT);
	      TopParticleeta_.push_back(genpart.at(i)->Eta); 
	      TopParticlephi_.push_back(genpart.at(i)->Phi); 
	      TopParticleid_.push_back(genpart.at(i)->PID); 
	      TopParticleMass_.push_back(genpart.at(i)->Mass); 
	    }

	  }

	  METvalue_=met.at(0)->MET;
	  METphi_=met.at(0)->Phi;
	  METeta_=met.at(0)->Eta;

	  GenMETvalue_=genmet.at(0)->MET;
	  GenMETphi_=genmet.at(0)->Phi;
	  GenMETeta_=genmet.at(0)->Eta;

	  ScalarHTvalue_=scalarHt.at(0)->HT;

	  /*==SKIM==
	   * Access the branches of the skim
	   */
	  //std::vector<Electron> * skimelecs=electrons.content();
	  //for(size_t i=0;i<skimelecs->size();i++){
	  //	histo->Fill(skimelecs->at(i).PT);
	  //}
		
	  myskim->Fill();	  

	  pt_ .clear();
	  eta_.clear();
	  phi_.clear();

	  GenJetpt_ .clear();
	  GenJeteta_.clear();
	  GenJetphi_.clear();

	  GenParticlept_ .clear();
	  GenParticleeta_.clear();
	  GenParticlephi_.clear();
	  GenParticleid_.clear();
	  GenParticleMass_.clear();

	  MuonParticlept_ .clear();
	  MuonParticleeta_.clear();
	  MuonParticlephi_.clear();
	  MuonParticleid_.clear();
	  MuonParticleMass_.clear();
	  MuonParticleStatus_.clear();

	  LeptoQuarkParticlept_ .clear();
	  LeptoQuarkParticleeta_.clear();
	  LeptoQuarkParticlephi_.clear();
	  LeptoQuarkParticleid_.clear();
	  LeptoQuarkParticleMass_.clear();

	  TopParticlept_ .clear();
	  TopParticleeta_.clear();
	  TopParticlephi_.clear();
	  TopParticleid_.clear();
	  TopParticleMass_.clear();

	  mass_.clear();
	  jetflavour_.clear();
	  jetTauTag_.clear();
	  btag_.clear();
	  btagAlgo_.clear();
	  muonPt.clear();
	  muonEta.clear();
	  muonPhi.clear();
	  elecPt.clear();
	  elecEta.clear();
	  elecPhi.clear();

	  elecIsoVar.clear();
	  elecIsoVarRhoCorr.clear(); 
	  elecSumPtCharged.clear(); 
	  elecSumPtNeutral.clear();
	  elecSumPtChargedPU.clear(); 
	  elecSumPt.clear(); 
	  elecCharge.clear();

	  muonIsoVar.clear();
	  muonIsoVarRhoCorr.clear(); 
	  muonSumPtCharged.clear(); 
	  muonSumPtNeutral.clear();
	  muonSumPtChargedPU.clear(); 
	  muonSumPt.clear(); 
	  muonCharge.clear();

	  VertexX_.clear();
	  VertexY_.clear();
	  VertexZ_.clear();

	}

	
	/*
	 * Must be called in the end, takes care of thread-safe writeout and
	 * call-back to the parent process
	 */
	processEndFunction();
}



void ntupleQCD::postProcess(){
	/*
	 * This function can be used to analyse the output histograms, e.g. extract a signal contribution etc.
	 * The function can also be called directly on an output file with the histograms, if
	 * RunOnOutputOnly = true
	 * is set in the analyser's config file
	 *
	 * This function also represents an example of how the output of the analyser can be
	 * read-back in an external program.
	 * Just include the sampleCollection.h header and follow the procedure below
	 *
	 */

	/*
	 * Here, the input file to the extraction of parameters from the histograms is the output file
	 * of the parallelised analysis.
	 * The sampleCollection class can also be used externally for accessing the output consistently
	 */
  //d_ana::sampleCollection samplecoll;
  //	samplecoll.readFromFile(getOutPath());
  
  //	std::vector<TString> alllegends = samplecoll.listAllLegends();
  
  /*
   * Example how to process the output.
   * Usually, one would define the legendname of the histogram to be used here
	 * by hand, e.g. "signal" or "background".
	 * To make this example run in any case, I am using alllegends.at(0), which
	 * could e.g. be the signal legend.
	 *
	 * So in practise, the following would more look like
	 * samplecoll.getHistos("signal");
	 */
  //if(alllegends.size()>0){
  //d_ana::histoCollection histos=samplecoll.getHistos(alllegends.at(0));

		/*
		 * here, the histogram created in the analyze() function is selected and evaluated
		 * The histoCollection maintains ownership (you don't need to delete the histogram)
		 */
		//const TH1* myplot=histos.getHisto("histoname1"); //commented by PG

		//std::cout << "(example output): the integral is " << myplot->Integral() <<std::endl;

		/*
		 * If the histogram is subject to changes, please clone it and take ownership
		 */

		//TH1* myplot2=histos.cloneHisto("histoname1"); //commented by PG

		/*
		 * do something with the histogram
		 */

		//delete myplot2;
  //	}

	/*
	 * do the extraction here.
	 */



}



