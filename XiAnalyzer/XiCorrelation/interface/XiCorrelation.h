// -*- C++ -*-
//
// Package:    XiCorrelation.h
// Class:      XiCorrelation.h
//
/**\class XiCorrelation XiCorrelation.h
 * RiceHIG/V0Analysis/interface/XiCorrelation.h

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author: Benjamin Tran
//
//

//#ifndef RICEHIG__XI_CORRELATION_H
//#define RICEHIG__XI_CORRELATION_H
// For interactive
#ifndef XIANALYZER__XI_CORRELATION_H
#define XIANALYZER__XI_CORRELATION_H

// System include files
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <math.h>
#include <map>
#include <sstream>

// user include files
#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TNtuple.h>
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TRandom.h>
#include <TNtuple.h>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/TrackReco/interface/DeDxData.h"

#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>

//////////////////////////////////////////////
// CMSSW user include files
#include "DataFormats/Common/interface/DetSetAlgorithm.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerLayerIdAccessor.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
//#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"

#include "DataFormats/SiPixelDetId/interface/PixelEndcapName.h"

#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMap.h"

#include "L1Trigger/GlobalTrigger/interface/L1GlobalTrigger.h"

#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

// Heavyion
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"

// Particle Flow
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"

// Vertex significance
#include "RecoBTag/SecondaryVertex/interface/SecondaryVertex.h"

// Root include files
#include "TTree.h"
//
// Track Matching and fake rate calculations     
//#include "RiceHIG/V0Analysis/interface/V0Validator.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

using namespace std;

class XiCorrelation : public edm::EDAnalyzer {
    public :
        explicit XiCorrelation(const edm::ParameterSet&);
        ~XiCorrelation();

    private :
        virtual void beginJob();
        virtual void analyze(const edm::Event&, const edm::EventSetup&);
        virtual void endJob();

        //virtual void beginRun(edm::Run const&, edm::EventSetup const&);
        //virtual void endRun(edm::Run const&, edm::EventSetup const&);
        //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
        //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
        //virtual void fillDescriptions(edm::ConfigurationDescriptions& descriptions const&);

        edm::Service<TFileService> fs;

        edm::InputTag xiCollection_;
        //edm::EDGetToken _xiCollection;
        edm::InputTag vertexCollName_;
        //edm::EDGetToken _vertexCollName;
        edm::InputTag trkSrc_;
        //edm::EDGetToken _trkSrc;

        double zVtxHigh_;
        double zVtxLow_;
        double ptMax_trg_;
        double ptMin_trg_;
        double ptMax_ass_;
        double ptMin_ass_;
        double ptMax_xi_;
        double ptMin_xi_;
        double etaMax_trg_;
        double etaMin_trg_;
        double etaMax_ass_;
        double etaMin_ass_;
        double bkgnum_;
        double xiMassHigh_;
        double xiMassLow_;
        int multHigh_;
        int multLow_;
        bool dopTcut_;

        vector<TVector3> *pepVect_trkhad;
        vector<TVector3> *pepVect_trkass;
        vector<TVector3> *pepVect_Xi;
        vector<double> *xiMass;
        vector<double> *zvtxVect;
        vector< vector<double> > *xiMass2;
        vector< vector<TVector3> > *PepVect2_Xi;
        vector< vector<TVector3> > *PepVect2_ass;
        vector< vector<TVector3> > *PepVect2_had;


        TH1D* InvMassXi;
        TH1D* PtSpectra;
        TH1D* nTrk;
        TH1D* EtaPtCutnTrackHist;
        TH1D* nEvtCut;
        TH1D* XiPerEvt;
        TH1D* XiPerEvtPeak;
        TH1D* XiPerEvtSide;
        TH1D* HadPerEvt;
        TH1D* TrkassPerEvt;
        //TH1D* TrktrgPerEvt;

        TH2D* BackgroundXi;
        TH2D* BackgroundXiPeak;
        TH2D* BackgroundXiSide;
        TH2D* BackgroundHad;
        TH2D* SignalXi;
        TH2D* SignalXiPeak;
        TH2D* SignalXiSide;
        TH2D* SignalHad;
        TH2D* Correlation;
        TH2D* CorrelationPeak;
        TH2D* CorrelationSide;
        TH2D* CorrelationHad;

        TH1D* PtBin12;
        TH1D* PtBin13;
        TH1D* PtBin14;
        TH1D* PtBin15;
        TH1D* PtBin16;
        TH1D* PtBin23;
        TH1D* PtBin24;
        TH1D* PtBin25;
        TH1D* PtBin26;
        TH1D* PtBin34;
        TH1D* PtBin35;
        TH1D* PtBin36;
        TH1D* PtBin45;
        TH1D* PtBin46;
        TH1D* PtBin56;
        TH1D* PtBin1;
        TH1D* PtBin2;
        TH1D* PtBin3;
        TH1D* PtBin4;
        TH1D* PtBin5;
        TH1D* PtBin6;
        TH1D* PtBinOver6;

};

#endif
