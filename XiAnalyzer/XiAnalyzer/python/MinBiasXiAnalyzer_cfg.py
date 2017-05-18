import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")
### standard includes
process.load("Configuration.StandardSequences.Digi_cff")
#process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
#process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("SimGeneral.MixingModule.mixNoPU_cfi")

### To use the TransientTrakcBuilder
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

### GlobalTag
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = 'GR_P_V43F::All'
process.GlobalTag.globaltag = 'GR_R_53_V21A::All'

process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")
process.load("SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(5000)
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    '/store/user/davidlw/PAMinBiasUPC/PA2013_FlowCorr_PromptReco_MB_Gplus_Rereco_ReTracking_v18/25c9a89be536a41c8ccb3c75e9fd9358/pPb_HM_1000_1_Bgt.root',
    '/store/user/davidlw/PAMinBiasUPC/PA2013_FlowCorr_PromptReco_MB_Gplus_Rereco_ReTracking_v18/25c9a89be536a41c8ccb3c75e9fd9358/pPb_HM_1001_1_NY3.root',
    '/store/user/davidlw/PAMinBiasUPC/PA2013_FlowCorr_PromptReco_MB_Gplus_Rereco_ReTracking_v18/25c9a89be536a41c8ccb3c75e9fd9358/pPb_HM_1002_1_Qga.root',
    '/store/user/davidlw/PAMinBiasUPC/PA2013_FlowCorr_PromptReco_MB_Gplus_Rereco_ReTracking_v18/25c9a89be536a41c8ccb3c75e9fd9358/pPb_HM_1003_1_oEd.root',
    '/store/user/davidlw/PAMinBiasUPC/PA2013_FlowCorr_PromptReco_MB_Gplus_Rereco_ReTracking_v18/25c9a89be536a41c8ccb3c75e9fd9358/pPb_HM_1004_1_d8Z.root',
    '/store/user/davidlw/PAMinBiasUPC/PA2013_FlowCorr_PromptReco_MB_Gplus_Rereco_ReTracking_v18/25c9a89be536a41c8ccb3c75e9fd9358/pPb_HM_1005_1_ssw.root'
   )
)


process.ana = cms.EDAnalyzer('XiAnalyzer',
                              vertexSrc = cms.string('offlinePrimaryVertices'),
                              trackSrc = cms.InputTag('generalTracks'),
                            
                              generalV0_xi = cms.InputTag('generalV0CandidatesLowPt:Xi'), 

                              genParticleSrc = cms.InputTag('genParticles'),
                              doGenParticle = cms.untracked.bool(False), 
                              doV0_Xi = cms.untracked.bool(True),
                              doGeneralTracks = cms.untracked.bool(True),
                              DoWeight = cms.untracked.bool(False) ,
                              iG = cms.untracked.int32(0)
)

process.TFileService = cms.Service("TFileService",fileName = cms.string("XiAnalysisPbPb_MinBiasNew.root"))

process.p = cms.Path(process.ana)
