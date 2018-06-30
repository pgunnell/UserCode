# -*- coding: utf-8 -*-                                                                                                                                                                                      
import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")
process.TFileService=cms.Service("TFileService",fileName=cms.string("Output.root"))

##-------------------- Define the source  ----------------------------                                                                                                                                       
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(1)
        )

process.source = cms.Source("EmptySource")


##-------------------- User analyzer  --------------------------------                                                                                                                                       
process.analysis  = cms.EDAnalyzer('Analysis_Template_MC',
                                   filename        = cms.string('file:example0_ntuple.root'),
                                   treename        = cms.string('Delphes'),
                                   dirname         = cms.string('boostedAK8'),
                                   CrossSection    = cms.untracked.double(100),
                                   )

process.p = cms.Path(process.analysis)


