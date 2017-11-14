import FWCore.ParameterSet.Config as cms 
process = cms.Process('myprocess')
process.TFileService=cms.Service("TFileService",fileName=cms.string('flatTreeFileData.root'))
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = '92X_dataRun2_HLT_v7'
##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(20000))
process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
        #"file:10B2302E-DA08-E611-884E-001E67DDD348.root")
        #"file:004A0552-3929-E611-BD44-0025905A48F0.root"),
        "file:006233DA-3599-E711-911B-0025905A6118.root"),
)
#############   Format MessageLogger #################
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

from PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi import selectedPatJets
process.goodJets = selectedPatJets.clone(src='slimmedJets',cut='pt>30 & abs(eta)<2.4')

process.load('RecoJets.JetProducers.QGTagger_cfi')
process.QGTagger.srcJets   = cms.InputTag('slimmedJets')
process.QGTagger.jetsLabel = cms.string('QGL_AK4PFchs')

from RecoJets.Configuration.RecoGenJets_cff import ak4GenJets, ak8GenJets

genParticleCollection = 'prunedGenParticles'

from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
from RecoJets.JetProducers.GenJetParameters_cfi import *

process.ak8GenJetsCustom = ak4GenJets.clone(
    src = genParticleCollection,
    rParam = cms.double(0.8),
    jetAlgorithm = cms.string("AntiKt")
)

genParticleCollection = 'prunedGenParticles'
genJetCollection = 'slimmedGenJetsAK8'

from PhysicsTools.JetMCAlgos.HadronAndPartonSelector_cfi import selectedHadronsAndPartons
process.selectedHadronsAndPartons = selectedHadronsAndPartons.clone(
    particles = genParticleCollection
)

from PhysicsTools.JetMCAlgos.AK4PFJetsMCFlavourInfos_cfi import ak4JetFlavourInfos
process.ak8genJetFlavourInfos = ak4JetFlavourInfos.clone(
    jets = genJetCollection,
    rParam = cms.double(0.8),
)

#only needed if information of the associated b hadrons are required
from PhysicsTools.JetMCAlgos.GenHFHadronMatcher_cff import matchGenBHadron
process.matchGenBHadron = matchGenBHadron.clone(
    genParticles = genParticleCollection,
    jetFlavourInfos = cms.InputTag("ak8genJetFlavourInfos"),
    flavour = cms.int32(5),
    onlyJetClusteredHadrons = cms.bool(True),
    noBBbarResonances = cms.bool(False),
)

##-------------------- User analyzer  --------------------------------
process.boosted = cms.EDAnalyzer('miniNtupleProducer',
  jets             = cms.InputTag('slimmedJetsAK8'),
  muons            = cms.InputTag('slimmedMuons'),
  electrons        = cms.InputTag('slimmedElectrons'),
  met              = cms.InputTag('slimmedMETs'),
  vertices         = cms.InputTag('offlineSlimmedPrimaryVertices'),
  rho              = cms.InputTag('fixedGridRhoFastjetAll'),
  candidates       = cms.InputTag('packedPFCandidates'),
  triggerPrescales = cms.InputTag('patTrigger'),
  xmlFile          = cms.string('boosted_mva_Fisher.weights.xml'),
  nJetsMin         = cms.int32(6),
  nBJetsMin        = cms.int32(2),
  ptMin            = cms.double(50), 
  minMuPt          = cms.double(30),                              
  minElPt          = cms.double(30),                              
  ptMinLeading     = cms.double(50),   
  massMin          = cms.double(50),
  htMin            = cms.double(5),
  etaMax           = cms.double(5.0),
  kinfit           = cms.string('kinFitTtFullHadEvent'),
  btagMinThreshold = cms.double(0.814),
  btagMin          = cms.double(0.1),                     
  btagMaxThreshold = cms.double(1.1),
  btagger          = cms.string('pfCombinedInclusiveSecondaryVertexV2BJetTags'),
  qgtagger         = cms.InputTag('QGTagger','qgLikelihood'),
  pu               = cms.untracked.string("addPileupInfo"),
  genparticles     = cms.untracked.InputTag('prunedGenParticles'),  
  triggerNames     = cms.vstring('HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20_v','HLT_PFJet40_v1','HLT_PFJet60_v1', 'HLT_PFJet80_v1', 'HLT_PFJet140_v', 'HLT_PFJet200_v1', 'HLT_PFJet260_v1','HLT_PFJet320_v1', 'HLT_PFJet400_v1', 'HLT_PFJet450_v1','HLT_PFJet500_v1'),
  triggerResults   = cms.InputTag('TriggerResults','','HLT'),
  isMC             = cms.untracked.bool(False),                              
  genjets          = cms.untracked.InputTag('slimmedGenJetsAK8'),
  GenptMin        = cms.untracked.double(50),
  GenetaMax       = cms.untracked.double(2.5),
  jetFlavourInfos = cms.InputTag("ak8genJetFlavourInfos"),                 
  isPrint         = cms.untracked.bool(False),                           
)

process.load('TopQuarkAnalysis.TopKinFitter.TtFullHadKinFitProducer_cfi')

process.kinFitTtFullHadEvent.jets                = 'goodJets'
process.kinFitTtFullHadEvent.bTagAlgo            = 'pfCombinedInclusiveSecondaryVertexV2BJetTags'
process.kinFitTtFullHadEvent.minBTagValueBJet    = 0.814
process.kinFitTtFullHadEvent.maxBTagValueNonBJet = 0.814
process.kinFitTtFullHadEvent.bTags               = 2
process.kinFitTtFullHadEvent.maxNJets            = 8

process.kinFitTtFullHadEventNoW = process.kinFitTtFullHadEvent.clone(constraints = cms.vuint32(5))
process.boostedNoW = process.boosted.clone(kinfit = 'kinFitTtFullHadEventNoW')

process.p = cms.Path(
   process.goodJets * 
   process.kinFitTtFullHadEvent * 
   process.boosted 
)








