# -*- coding: utf-8 -*-
import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")
process.TFileService=cms.Service("TFileService",fileName=cms.string("Test-1.root"))

##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(1)
        )

process.source = cms.Source("EmptySource")


##-------------------- User analyzer  --------------------------------
process.analysis  = cms.EDAnalyzer('Analysis_Template_MCFastSim',
                                   filename        = cms.string('file:/nfs/dust/cms/user/gunnep/DelphesAnalysis/ntuplizer/outputNew/LG_1300_S130_ntuple.root'),
                                   treename        = cms.string('Delphes'),
                                   dirname         = cms.string(''),
                                   isMCarlo        = cms.untracked.bool(True),
                                   CrossSection    = cms.untracked.double(0.00009817),
                                   IntLumi         = cms.untracked.double(100),
                                   Triggers        = cms.untracked.vstring(''),
                                   ElSelection     = cms.untracked.bool(False),
                                   MuSelection     = cms.untracked.bool(True),
                                   HadSelection    = cms.untracked.bool(False),
                                   )

process.p = cms.Path(process.analysis)
