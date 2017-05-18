import FWCore.ParameterSet.Config as cms

import HLTrigger.HLTfilters.hltHighLevel_cfi
hltHM = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
hltHM.HLTPaths = [
                'HLT_PAPixelTracks_Multiplicity100_v*',
                'HLT_PAPixelTracks_Multiplicity130_v*',
                'HLT_PAPixelTracks_Multiplicity160_v*'
                #'HLT_PAPixelTracks_Multiplicity190_v',
                #'HLT_PAPixelTracks_Multiplicity220_v'
            ]

xiCorrelation          = cms.EDAnalyzer('XiCorrelation',
        trkSrc         = cms.InputTag('generalTracks'),
        xiCollection   = cms.InputTag('selectV0CandidatesLowXi:Xi'),
        vertexCollName = cms.InputTag('offlinePrimaryVertices'),
        zVtxHigh       = cms.double(15.0),
        zVtxLow        = cms.double(-15.0),
        ptMax_trg      = cms.double(3.0),
        ptMin_trg      = cms.double(1.0),
        ptMax_ass      = cms.double(3.0),
        ptMin_ass      = cms.double(1.0),
        ptMax_xi       = cms.double(3.0),
        ptMin_xi       = cms.double(1.0),
        etaMax_trg     = cms.double(2.4),
        etaMin_trg     = cms.double(-2.4),
        etaMax_ass     = cms.double(2.4),
        etaMin_ass     = cms.double(-2.4),
        multHigh       = cms.int32(220),
        multLow        = cms.int32(185),
        bkgnum         = cms.double(20.0),
        xiMassLow      = cms.double(1.3118),
        xiMassHigh     = cms.double(1.3326),
        dopTcut        = cms.untracked.bool(True)
        )
