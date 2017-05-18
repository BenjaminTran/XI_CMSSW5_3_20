import FWCore.ParameterSet.Config as cms

selectV0CandidatesLow   = cms.EDProducer('XiSelector',
    vertexCollName     = cms.InputTag('offlinePrimaryVertices'),
    v0CollName         = cms.string("generalV0CandidatesLowPt"),
    v0IDName           = cms.string("Xi"),
    etaCutMin          = cms.double(-2.4),
    etaCutMax          = cms.double(2.4),
    zVertexLow         = cms.double(-15.0),
    zVertexHigh        = cms.double(15.0),
    ptCut1             = cms.double(0.0),
    ptCut2             = cms.double(0.0),
    nHitCut1           = cms.int32(3),
    nHitCut2           = cms.int32(3),
    xi3DIpSigValue     = cms.double(2.5),
    xiPi3DIpSigValue   = cms.double(5),
    VTrkPi3DIpSigValue = cms.double(4),
    VTrkP3DIpSigValue  = cms.double(3),
    xiFlightSigValue   = cms.double(3),
    distanceSigValue   = cms.double(12),
#    dxySigCut_xi_la   = cms.double(1.0), #lambda
#    dzSigCut_xi_la    = cms.double(1.0), #lambda
#    dxySigCut_xi_pi   = cms.double(3.0), #pion transvers impact parameter exceed 3
#    dzSigCut_xi_pi    = cms.double(3.0), #pion longitudinal impact parameter exceed 3
#    dxySigCut_pro     = cms.double(3.0), #lambda decay proton; larger than 3
#    dzSigCut_pro      = cms.double(3.0), #lambda decay proton; larger than 3
#    dxySigCut_la_pi   = cms.double(4.0), #lambda decay pion; larger than 4
#    dxySigCut_la_pi   = cms.double(4.0), #lambda decay pion; larger than 4
#    dxySigCut_xi      = cms.double(2.5),
#    dzSigCut_xi       = cms.double(2.5),
    vtxChi2Cut         = cms.double(10000.0),
    cosThetaCut        = cms.double(0.999),
#    decayLSigCut_xi   = cms.double(3.0),
#    decayLSigCut_la   = cms.double(12.0),
    misIDMassCut       = cms.double(0.010), #Ask Wei what this value should be and what it is
    misIDMassCutEE     = cms.double(0.015)
    )
