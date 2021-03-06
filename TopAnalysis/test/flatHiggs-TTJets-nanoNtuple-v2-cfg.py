import FWCore.ParameterSet.Config as cms 
process = cms.Process('myprocess')
process.TFileService=cms.Service("TFileService",fileName=cms.string('flatTreeFileHiggs-nano.root'))
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = '94X_mc2017_realistic_v13'
from PhysicsTools.PatAlgos.tools.jetTools import *

##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(2000))
process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
        #"file:10B2302E-DA08-E611-884E-001E67DDD348.root")
        #"file:004A0552-3929-E611-BD44-0025905A48F0.root"),
        'root://cms-xrd-global.cern.ch//store/user/pgunnell/HiggsZProduction/CRAB3_HZ_HiggsToWWZToLL_Pythia8_TuneCP5-DRMiniAOD/180212_213039/0000/SMP-RunIIFall17MiniAOD-00006_1.root'
        ),
)

#create a task
task = cms.Task()

#############   Format MessageLogger #################
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)


from RecoJets.Configuration.RecoPFJets_cff import ak4PFJets, ak4PFJetsCHS
from RecoJets.Configuration.RecoPFJets_cff import ak8PFJetsCHS

#run jets with a different cone size
# PUPPI JETS
process.load('CommonTools/PileupAlgos/Puppi_cff')
## e.g. to run on miniAOD
process.puppi.candName = cms.InputTag('packedPFCandidates')
process.puppi.vertexName = cms.InputTag('offlineSlimmedPrimaryVertices')
process.puppi.clonePackedCands   = cms.bool(True)
task.add(process.puppi)

#create AK15PuppiJets
process.ak15PuppiJets = ak8PFJetsCHS.clone(src = cms.InputTag('puppi'), rParam = 1.5)#, R0 = 1.5, zcut = cms.double(0.1), beta = cms.double(0))
#create CA15PuppiJets
process.ca15PuppiJets = ak8PFJetsCHS.clone(src = cms.InputTag('puppi'), rParam = 1.5, jetAlgorithm = "CambridgeAachen")#, R0 = 1.5, zcut = cms.double(0.1), beta = cms.double(0))

#TestAK8 puppi
process.ak8PuppiJets = ak8PFJetsCHS.clone(src = cms.InputTag('puppi'), doAreaFastjet = True, rParam = 0.8)#, R0 = 1.5, zcut = cms.double(0.1), beta = cms.double(0))

task.add(process.ak15PuppiJets)
task.add(process.ak8PuppiJets)
task.add(process.ca15PuppiJets)

#create the GenJet collection (AK15 and CA15)
from RecoJets.Configuration.RecoGenJets_cff import ak4GenJets, ak8GenJets

from RecoJets.JetProducers.GenJetParameters_cfi import *
from RecoJets.JetProducers.AnomalousCellParameters_cfi import *

process.packedGenParticlesForJetsNoNu = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedGenParticles"), cut = cms.string("abs(pdgId) != 12 && abs(pdgId) != 14 && abs(pdgId) != 16"))
task.add(process.packedGenParticlesForJetsNoNu)

setattr(process,
        "ak15GenJets",
        ak4GenJets.clone(
        src = cms.InputTag("packedGenParticlesForJetsNoNu"),
        rParam       = cms.double(1.5))
        )

task.add(process.ak15GenJets)

setattr(process,
        "ak8GenJetsReClust",
        ak4GenJets.clone(
        src = cms.InputTag("packedGenParticlesForJetsNoNu"),
        rParam       = cms.double(0.8))
        )

task.add(process.ak8GenJetsReClust)

setattr(process,
        "ca15GenJets",
        ak4GenJets.clone(
        src = cms.InputTag("packedGenParticlesForJetsNoNu"),
        rParam       = cms.double(1.5),
        jetAlgorithm = cms.string("CambridgeAachen"),
        ))

task.add(process.ca15GenJets)

task.add(getattr(process, "ca15GenJets"))
task.add(getattr(process, "ak15GenJets"))

#setattr(process, 'ak15GenJets' )
#setattr(process, 'ca15GenJets', )

addJetCollection(process,labelName = 'ak15PuppiJets', 
                 jetSource = cms.InputTag('ak15PuppiJets'), 
                 algo = 'AK', rParam=1.5, 
                 #genJetCollection=cms.InputTag('slimmedGenJetsAK8'), 
                 genJetCollection=cms.InputTag('ak15GenJets'), 
                 jetCorrections = ('AK8PFPuppi', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None'),
                 pfCandidates = cms.InputTag('packedPFCandidates'),
                 pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
                 svSource = cms.InputTag('slimmedSecondaryVertices'),
                 muSource =cms.InputTag( 'slimmedMuons'),
                 elSource = cms.InputTag('slimmedElectrons'),
                 genParticles = cms.InputTag('prunedGenParticles'),
                 getJetMCFlavour=False,
)

addJetCollection(process,labelName = 'ak8PuppiJets', 
                 jetSource = cms.InputTag('ak8PuppiJets'), 
                 algo = 'AK', rParam=0.8, 
                 genJetCollection=cms.InputTag('slimmedGenJetsAK8'), 
                 jetCorrections = ('AK8PFPuppi', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None'),
                 pfCandidates = cms.InputTag('packedPFCandidates'),
                 pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
                 svSource = cms.InputTag('slimmedSecondaryVertices'),
                 muSource =cms.InputTag( 'slimmedMuons'),
                 elSource = cms.InputTag('slimmedElectrons'),
                 genParticles = cms.InputTag('prunedGenParticles'),
                 getJetMCFlavour=False,
)

addJetCollection(process,labelName = 'ca15PuppiJets', 
                 jetSource = cms.InputTag('ca15PuppiJets'), 
                 algo = 'CA', rParam=1.5, 
                 genJetCollection=cms.InputTag('ca15GenJets'), 
                 jetCorrections = ('AK8PFPuppi', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None'),
                 pfCandidates = cms.InputTag('packedPFCandidates'),
                 pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
                 svSource = cms.InputTag('slimmedSecondaryVertices'),
                 muSource =cms.InputTag( 'slimmedMuons'),
                 elSource = cms.InputTag('slimmedElectrons'),
                 genParticles = cms.InputTag('prunedGenParticles'),
                 getJetMCFlavour=False,
)

#patJetsAk15PFJetsPuppi
task.add(getattr(process, 'patJetsAk15PuppiJets' ))
task.add(getattr(process, 'patJetsAk8PuppiJets' ))
task.add(getattr(process, 'patJetsCa15PuppiJets' ))

#adding subjettiness
process.load('RecoJets.JetProducers.nJettinessAdder_cfi')
task.add(process.Njettiness)
process.NjettinessAK15 = process.Njettiness.clone() #probably you need to change something here
task.add(process.NjettinessAK15)
process.NjettinessAK15.src = cms.InputTag("ak15PuppiJets")

getattr(process,"patJetsAk15PuppiJets").userData.userFloats.src += ['NjettinessAK15:tau1','NjettinessAK15:tau2','NjettinessAK15:tau3','NjettinessAK15:tau4']

process.NjettinessCA15 = process.Njettiness.clone() #probably you need to change something here
task.add(process.NjettinessCA15)
process.NjettinessCA15.src = cms.InputTag("ca15PuppiJets")

getattr(process,"patJetsCa15PuppiJets").userData.userFloats.src += ['NjettinessCA15:tau1','NjettinessCA15:tau2','NjettinessCA15:tau3','NjettinessCA15:tau4']

#For doing soft drop algorithm to an existing jet
#From UHH2 software, modified
#one has to create both soft drop jets and soft drop subjets

process.ak15PuppiJetsSoftDrop = process.ak15PuppiJets.clone(
    rParam=1.5,
    beta = cms.double(0.0),
    jetCollInstanceName = cms.string("SubJets"),
    useExplicitGhosts = cms.bool(True),
    useSoftDrop = cms.bool(True),
    writeCompound = cms.bool(True),
    zcut = cms.double(0.1),
    R0=cms.double(1.5)
)
task.add(process.ak15PuppiJetsSoftDrop)

process.ca15PuppiJetsSoftDrop = process.ca15PuppiJets.clone(
    rParam=1.5,
    beta = cms.double(0.0),
    jetCollInstanceName = cms.string("SubJets"),
    useExplicitGhosts = cms.bool(True),
    useSoftDrop = cms.bool(True),
    writeCompound = cms.bool(True),
    zcut = cms.double(0.1),
    R0=cms.double(1.5)
)
task.add(process.ca15PuppiJetsSoftDrop)

#process.ak15PuppiJetsSoftDropforsub = process.ak15PuppiJets.clone(
#    rParam=1.5,
#    beta = cms.double(0.0),
#    jetCollInstanceName = cms.string("SubJets"),
#    useExplicitGhosts = cms.bool(True),
#    useSoftDrop = cms.bool(True),
#    writeCompound = cms.bool(True),
#    zcut = cms.double(0.1),
#    R0=cms.double(1.5)
#)
#task.add(process.ak15PuppiJetsSoftDropforsub)


#Now patifyng both the fat jets and the subjet
addJetCollection(process,labelName = 'ak15PuppiJetsSoftDrop',
                 jetSource = cms.InputTag('ak15PuppiJetsSoftDrop'),
                 algo = 'AK', rParam=1.5,
                 genJetCollection=cms.InputTag('ak15GenJets'),
                 jetCorrections = ('AK8PFPuppi', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None'),
                 pfCandidates = cms.InputTag('packedPFCandidates'),
                 pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
                 svSource = cms.InputTag('slimmedSecondaryVertices'),
                 muSource =cms.InputTag( 'slimmedMuons'),
                 elSource = cms.InputTag('slimmedElectrons'),
                 genParticles = cms.InputTag('prunedGenParticles'),
                 getJetMCFlavour=False,
)

task.add(getattr(process, 'patJetsAk15PuppiJetsSoftDrop' ))

addJetCollection(process,labelName = 'ak15PuppiJetsSoftDropforsub',
                 jetSource = cms.InputTag('ak15PuppiJetsSoftDrop',"SubJets"),
                 algo = 'AK', rParam=1.5,
                 genJetCollection=cms.InputTag('ak15GenJets'),
                 jetCorrections = ('AK8PFPuppi', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None'),
                 pfCandidates = cms.InputTag('packedPFCandidates'),
                 pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
                 svSource = cms.InputTag('slimmedSecondaryVertices'),
                 muSource =cms.InputTag( 'slimmedMuons'),
                 elSource = cms.InputTag('slimmedElectrons'),
                 genParticles = cms.InputTag('prunedGenParticles'),
                 getJetMCFlavour=False,
)

task.add(getattr(process, 'patJetsAk15PuppiJetsSoftDropforsub' ))

#CA15 jets

addJetCollection(process,labelName = 'ca15PuppiJetsSoftDrop',
                 jetSource = cms.InputTag('ca15PuppiJetsSoftDrop'),
                 algo = 'CA', rParam=1.5,
                 genJetCollection=cms.InputTag('ca15GenJets'),
                 jetCorrections = ('AK8PFPuppi', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None'),
                 pfCandidates = cms.InputTag('packedPFCandidates'),
                 pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
                 svSource = cms.InputTag('slimmedSecondaryVertices'),
                 muSource =cms.InputTag( 'slimmedMuons'),
                 elSource = cms.InputTag('slimmedElectrons'),
                 genParticles = cms.InputTag('prunedGenParticles'),
                 getJetMCFlavour=False,
)

task.add(getattr(process, 'patJetsCa15PuppiJetsSoftDrop' ))

addJetCollection(process,labelName = 'ca15PuppiJetsSoftDropforsub',
                 jetSource = cms.InputTag('ca15PuppiJetsSoftDrop',"SubJets"),
                 algo = 'CA', rParam=1.5,
                 genJetCollection=cms.InputTag('ca15GenJets'),
                 jetCorrections = ('AK8PFPuppi', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None'),
                 pfCandidates = cms.InputTag('packedPFCandidates'),
                 pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
                 svSource = cms.InputTag('slimmedSecondaryVertices'),
                 muSource =cms.InputTag( 'slimmedMuons'),
                 elSource = cms.InputTag('slimmedElectrons'),
                 genParticles = cms.InputTag('prunedGenParticles'),
                 getJetMCFlavour=False,
)

task.add(getattr(process, 'patJetsCa15PuppiJetsSoftDropforsub' ))

#Soft drop mass inclusion
#from RecoJets.JetProducers.ak8PFJetsPuppi_groomingValueMaps_cfi import ak8PFJetsPuppiSoftDropMass

#process.ak15PuppiSoftDropMass = ak8PFJetsPuppiSoftDropMass.clone(
#    distMax = cms.double(1.5),
#    matched = cms.InputTag("ak15PuppiJets"),
#    src = cms.InputTag("ak15PuppiJetsSoftDrop"),
#    value = cms.string('mass')
#)
#task.add(process.ak15PuppiSoftDropMass)

#getattr(process,"patJetsAk15PuppiJetsSoftDrop").userData.userFloats.src += ['ak15PuppiJetsSoftDropMass']

#slimming of the subjet collection
#process.slimmedJetsAk15PuppiJetsSoftDropforsub = cms.EDProducer("PATJetSlimmer",
#                                                            src = cms.InputTag("patJetsAk15PuppiJetsSoftDropforsub"),
#                                                            packedPFCandidates = cms.InputTag("packedPFCandidates"),
#                                                            dropJetVars = cms.string("1"),
#                                                            dropDaughters = cms.string("0"),
#                                                            rekeyDaughters = cms.string("1"),
#                                                            dropTrackRefs = cms.string("1"),
#                                                            dropSpecific = cms.string("1"),
#                                                            dropTagInfos = cms.string("1"),
#                                                            modifyJets = cms.bool(True),
#                                                            mixedDaughters = cms.bool(False),
#                                                            modifierConfig = cms.PSet( modifications = cms.VPSet() )
#)

#task.add(process.slimmedJetsAk15PuppiJetsSoftDropforsub)

#merging the two collections in one jet
process.patJetsAk15PuppiJetsSoftDropPacked = cms.EDProducer("BoostedJetMerger",
    jetSrc = cms.InputTag("patJetsAk15PuppiJetsSoftDrop"),
    subjetSrc = cms.InputTag("patJetsAk15PuppiJetsSoftDropforsub")
)

task.add(process.patJetsAk15PuppiJetsSoftDropPacked)

#To add subjets to the ungroomed fat jets - then, you can use them
# From UHH2 software, modified
# Add subjets from groomed fat jet to its corresponding ungroomed fatjet

process.packedPatJetsAk15PuppiJets = cms.EDProducer("JetSubstructurePacker",
    algoLabels = cms.vstring('SoftDropPuppi'),
    algoTags = cms.VInputTag(cms.InputTag("patJetsAk15PuppiJetsSoftDropPacked")),
    distMax = cms.double(1.5),
    fixDaughters = cms.bool(False),
    jetSrc = cms.InputTag("patJetsAk15PuppiJets")
)

task.add(process.packedPatJetsAk15PuppiJets)

#CA15 jets
process.patJetsCa15PuppiJetsSoftDropPacked = cms.EDProducer("BoostedJetMerger",
    jetSrc = cms.InputTag("patJetsCa15PuppiJetsSoftDrop"),
    subjetSrc = cms.InputTag("patJetsCa15PuppiJetsSoftDropforsub")
)

task.add(process.patJetsCa15PuppiJetsSoftDropPacked)

#To add subjets to the ungroomed fat jets - then, you can use them
# From UHH2 software, modified
# Add subjets from groomed fat jet to its corresponding ungroomed fatjet

process.packedPatJetsCa15PuppiJets = cms.EDProducer("JetSubstructurePacker",
    algoLabels = cms.vstring('SoftDropPuppi'),
    algoTags = cms.VInputTag(cms.InputTag("patJetsCa15PuppiJetsSoftDropPacked")),
    distMax = cms.double(1.5),
    fixDaughters = cms.bool(False),
    jetSrc = cms.InputTag("patJetsCa15PuppiJets")
)

task.add(process.packedPatJetsCa15PuppiJets)

#following things
process.load('RecoJets.JetProducers.QGTagger_cfi')
process.QGTagger.srcJets   = cms.InputTag('slimmedJets')
process.QGTagger.jetsLabel = cms.string('QGL_AK4PFchs')

from RecoJets.Configuration.RecoGenJets_cff import ak4GenJets, ak8GenJets

genParticleCollection = 'prunedGenParticles'

from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
from RecoJets.JetProducers.GenJetParameters_cfi import *

process.ak4GenJetsCustom = ak4GenJets.clone(
    src = genParticleCollection,
    rParam = cms.double(0.4),
    jetAlgorithm = cms.string("AntiKt")
)

process.ak8GenJetsCustom = ak4GenJets.clone(
    src = genParticleCollection,
    rParam = cms.double(0.8),
    jetAlgorithm = cms.string("AntiKt")
)

genParticleCollection = 'prunedGenParticles'
genJetCollection = 'slimmedGenJets'

from PhysicsTools.JetMCAlgos.HadronAndPartonSelector_cfi import selectedHadronsAndPartons
process.selectedHadronsAndPartons = selectedHadronsAndPartons.clone(
    particles = genParticleCollection
)

from PhysicsTools.JetMCAlgos.AK4PFJetsMCFlavourInfos_cfi import ak4JetFlavourInfos
process.ak4genJetFlavourInfos = ak4JetFlavourInfos.clone(
    jets = genJetCollection,
    rParam = cms.double(0.4),
)

genJetCollection = 'slimmedGenJetsAK8'

process.ak8genJetFlavourInfos = ak4JetFlavourInfos.clone(
    jets = genJetCollection,
    rParam = cms.double(0.8),
)

#only needed if information of the associated b hadrons are required
from PhysicsTools.JetMCAlgos.GenHFHadronMatcher_cff import matchGenBHadron
process.matchGenBHadron = matchGenBHadron.clone(
    genParticles = genParticleCollection,
    jetFlavourInfos = cms.InputTag("ak4genJetFlavourInfos"),
    flavour = cms.int32(5),
    onlyJetClusteredHadrons = cms.bool(True),
    noBBbarResonances = cms.bool(False),
)

##-------------------- User analyzer  --------------------------------
process.boostedAK4 = cms.EDAnalyzer('nanoNtupleProducerCandidates',
  jets             = cms.InputTag('slimmedJets'),
  met              = cms.InputTag('slimmedMETs'),
  vertices         = cms.InputTag('offlineSlimmedPrimaryVertices'),
  rho              = cms.InputTag('fixedGridRhoFastjetAll'),
  candidates       = cms.InputTag('packedPFCandidates'),
  triggerPrescales = cms.InputTag('patTrigger'),
  triggerObjects = cms.InputTag("slimmedPatTrigger"),
  nJetsMin         = cms.int32(1),
  nBJetsMin        = cms.int32(2),
  muons            = cms.InputTag('slimmedMuons'),                                  
  electrons        = cms.InputTag('slimmedElectrons'),                                  
  ptMin            = cms.double(50), 
  minMuPt          = cms.double(30),                              
  minElPt          = cms.double(30),                              
  ptMinLeading     = cms.double(50),   
  massMin          = cms.double(50),
  htMin            = cms.double(5),
  etaMax           = cms.double(5.0),
  btagger          = cms.string('pfCombinedInclusiveSecondaryVertexV2BJetTags'),
  btagMinThreshold = cms.double(0.814),
  btagMin          = cms.double(0.1),                     
  btagMaxThreshold = cms.double(1.1),
  qgtagger         = cms.InputTag('QGTagger','qgLikelihood'),
  pu               = cms.untracked.string("addPileupInfo"),
  genparticles     = cms.untracked.InputTag('prunedGenParticles'),  
  triggerNames     = cms.vstring(''),
  triggerResults   = cms.InputTag('TriggerResults','','HLT'),
  isMC             = cms.untracked.bool(True),                              
  isHiggs          = cms.untracked.bool(True),                              
  genjets          = cms.untracked.InputTag('slimmedGenJets'),
  GenptMin        = cms.untracked.double(50),
  GenetaMax       = cms.untracked.double(2.5),
  jetFlavourInfos = cms.InputTag("ak4genJetFlavourInfos"),                 
  isPrint         = cms.untracked.bool(False),   
  jetlabel        = cms.string('AK4'),
)

#create also GenJets with different cone sizes

process.boostedAK8 = process.boostedAK4.clone(
    jets             = cms.InputTag('slimmedJetsAK8'),
    genjets          = cms.untracked.InputTag('slimmedGenJetsAK8'),
    #genjets          = cms.untracked.InputTag('ak8GenJetsReClust'),
    jetFlavourInfos = cms.InputTag("ak8genJetFlavourInfos"),                     
    jetRadius       = cms.double(0.8),
    jetlabel        = cms.string('AK8Puppi'),                                  
)

process.boostedAK8puppi = process.boostedAK4.clone(
    jets             = cms.InputTag('patJetsAk8PuppiJets'),
    genjets          = cms.untracked.InputTag('slimmedGenJetsAK8'),
    jetFlavourInfos = cms.InputTag("ak8genJetFlavourInfos"),                    
    jetRadius       = cms.double(0.8),                                   
    jetlabel        = cms.string('AK8'),
)

process.boostedAK15 = process.boostedAK4.clone(
    #jets             = cms.InputTag('patJetsAk15PuppiJets'),
    jets             = cms.InputTag('packedPatJetsAk15PuppiJets'),
    #jets             = cms.InputTag('patJetsAk15PuppiJetsSoftDrop'),
    genjets          = cms.untracked.InputTag('ak15GenJets'),
    jetFlavourInfos = cms.InputTag("ak8genJetFlavourInfos"),                     
    jetRadius       = cms.double(1.5),                                  
    jetlabel        = cms.string('AK15'),
)

process.boostedCA15 = process.boostedAK4.clone(
    jets             = cms.InputTag('packedPatJetsCa15PuppiJets'),
    genjets          = cms.untracked.InputTag('ca15GenJets'),
    jetFlavourInfos = cms.InputTag("ak8genJetFlavourInfos"),                     
    jetRadius       = cms.double(1.5),                                  
    jetlabel        = cms.string('CA15'),
)


process.p = cms.Path(
   #process.boostedAK4 *
   process.boostedAK8*  #*process.boostedAK8puppi
   process.boostedAK15*
   process.boostedCA15
)

process.p.associate(task)
process.p.associate(process.patAlgosToolsTask)
