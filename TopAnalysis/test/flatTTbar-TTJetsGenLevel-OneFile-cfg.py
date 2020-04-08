import FWCore.ParameterSet.Config as cms 
import os

process = cms.Process('myprocess')

process.TFileService=cms.Service("TFileService",fileName=cms.string('GenLevelOutput.root'))

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
       'file:PPD-RunIISummer15wmLHEGS-00008.root')
)
#############   Format MessageLogger #################
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

from RecoJets.Configuration.RecoGenJets_cff import ak4GenJets, ak8GenJets

genParticleCollection = 'genParticles'

from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
from RecoJets.JetProducers.GenJetParameters_cfi import *

##-------------------- User analyzer  --------------------------------
process.boosted = cms.EDAnalyzer('BoostedTTbarFlatTreeProducerGenLevel',
  genparticles     = cms.untracked.InputTag('genParticles'),  
  isMC             = cms.untracked.bool(True),                              
  genjets          = cms.untracked.InputTag('ak4GenJets'),
  GenptMin        = cms.untracked.double(200),
  GenetaMax       = cms.untracked.double(2.5),
  isPrint         = cms.untracked.bool(True),                           
  isHigherOrder   = cms.untracked.bool(True),                           
  EFT_weights     = cms.untracked.vstring("rwgt_SM","rwgt_ctG_min3p0_ctW_2p0")
)

process.p = cms.Path(
    process.boosted 
)







