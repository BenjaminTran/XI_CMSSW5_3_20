// -*- C++ -*-
//
// Package:    XiCorrelation
// Class:      XiCorrelation
// 
/**\class XiCorrelation XiCorrelation.cc XiAnalyzer/XiCorrelation/src/XiCorrelation.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Benjamin Lucky Tran,,,
//         Created:  Sat Nov 12 08:57:22 CET 2016
// $Id$
//
//

//#include "RiceHIG/V0Analysis/interface/XiCorrelation.h"
#include "XiAnalyzer/XiCorrelation/interface/XiCorrelation.h" // for interactive

#define PI 3.1416

XiCorrelation::XiCorrelation(const edm::ParameterSet& iConfig)

{
    //now do what ever initialization is needed
    using std::string;
    
    TH1::SetDefaultSumw2();
    trkSrc_         = iConfig.getParameter<edm::InputTag>( "trkSrc" );
    xiCollection_   = iConfig.getParameter<edm::InputTag>("xiCollection");
    vertexCollName_ = iConfig.getParameter<edm::InputTag>("vertexCollName");
    zVtxHigh_       = iConfig.getParameter<double>("zVtxHigh");
    zVtxLow_        = iConfig.getParameter<double>("zVtxLow");
    ptMax_trg_      = iConfig.getParameter<double>( "ptMax_trg" );
    ptMin_trg_      = iConfig.getParameter<double>( "ptMin_trg" );
    ptMax_ass_      = iConfig.getParameter<double>( "ptMax_ass" );
    ptMin_ass_      = iConfig.getParameter<double>( "ptMin_ass" );
    ptMax_xi_       = iConfig.getParameter<double>( "ptMax_xi" );
    ptMin_xi_       = iConfig.getParameter<double>( "ptMin_xi" );
    etaMax_trg_     = iConfig.getParameter<double>( "etaMax_trg" );
    etaMin_trg_     = iConfig.getParameter<double>( "etaMin_trg" );
    etaMax_ass_     = iConfig.getParameter<double>( "etaMax_ass" );
    etaMin_ass_     = iConfig.getParameter<double>( "etaMin_ass" );
    bkgnum_         = iConfig.getParameter<double>( "bkgnum" );
    multHigh_       = iConfig.getParameter<int>( "multHigh" );
    multLow_        = iConfig.getParameter<int>( "multLow" );
    //dopTcut_        = iConfig.getUntrackedParameter<bool>( "dopTcut" );
    
}


XiCorrelation::~XiCorrelation()
{
 
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//

// ------------ method called for each event  ------------
void
XiCorrelation::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    using namespace std;
     
    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByLabel(vertexCollName_, vertices);
    double bestvz      = -999, bestvx        = -999, bestvy        = -999;
    double bestvzError = -999.9, bestvxError = -999.9, bestvyError = -999.9;
    const reco::Vertex & vtx = (*vertices)[0];
    bestvx = vtx.x(  );
    bestvy = vtx.y(  );
    bestvz = vtx.z(  );
    bestvxError = vtx.xError(  );
    bestvyError = vtx.yError(  );
    bestvzError = vtx.zError(  );

    if( bestvz > zVtxHigh_ || bestvz < zVtxLow_ ) return;

    pepVect_trkhad = new vector<TVector3>;
    pepVect_trkass = new vector<TVector3>;
    pepVect_Xi     = new vector<TVector3>;

    edm::Handle<reco::VertexCompositeCandidateCollection> xiCollection;
    iEvent.getByLabel(xiCollection_, xiCollection);

    edm::Handle<reco::TrackCollection> tracks;
    iEvent.getByLabel( trkSrc_, tracks );

    if(!xiCollection.isValid()) return;

    int nTracks = 0;
    int EtaPtCutnTracks = 0;

    // Track selection
    for( unsigned it = 0; it < tracks->size(  ); it++ )
    {
        const reco::Track & trk = ( *tracks )[it];
        math::XYZPoint bestvtx( bestvx, bestvy, bestvz );

        double dzvtx    = trk.dz( bestvtx );
        double dxyvtx   = trk.dxy( bestvtx );
        double dzerror  = sqrt( trk.dzError(  )*trk.dzError(  ) + bestvzError*bestvzError );
        double dxyerror = sqrt( trk.d0Error(  )*trk.d0Error(  ) + bestvxError*bestvyError );

        if( !trk.quality( reco::TrackBase::highPurity ) ) continue;
        if( fabs( trk.ptError(  )/trk.pt(  ) > 0.10 ) )   continue;
        if( fabs( dzvtx/dzerror ) > 3 )                   continue;
        if( fabs( dxyvtx/dxyerror ) > 3 )                 continue;

        nTracks++;
        if( fabs( trk.eta(  ) ) > 2.4 || trk.pt(  ) < 0.4 ) continue;
        EtaPtCutnTracks++;
    }

    nTrk->Fill( nTracks );
    EtaPtCutnTrackHist->Fill( EtaPtCutnTracks );

    if( EtaPtCutnTracks >= multLow_ && EtaPtCutnTracks < multHigh_ ){
        nEvtCut->Fill( 1 );
        for(reco::VertexCompositeCandidateCollection::const_iterator xiCand =
                xiCollection->begin(); xiCand != xiCollection->end(); xiCand++ ) {
        
            
            // Make InvMass for Xi
            double mass = xiCand->mass();
            InvMassXi->Fill(mass);
            cout<<"Fill"<<endl;

            // Xi pT Spectrum
            double xi_pT = xiCand->pt();

            //bool dopTcut_ = true;
            //if( dopTcut_ )
            //{
            //    if( xi_pT < ptMin_xi_ || xi_pT > ptMax_xi_ ) continue;
            //}

            PtSpectra->Fill(xi_pT);

            // Xi pT bins
            if( xi_pT < 1.0 ) PtBin1->Fill( mass );
            if( xi_pT < 2.0 ) PtBin2->Fill( mass );
            if( xi_pT < 3.0 ) PtBin3->Fill( mass );
            if( xi_pT < 4.0 ) PtBin4->Fill( mass );
            if( xi_pT < 5.0 ) PtBin5->Fill( mass );
            if( xi_pT < 6.0 ) PtBin6->Fill( mass );
            if( xi_pT > 6.0 ) PtBinOver6->Fill( mass );

            if( xi_pT >= 1 && xi_pT < 2 ) PtBin12->Fill( mass );
            if( xi_pT >= 2 && xi_pT < 3 ) PtBin23->Fill( mass );
            if( xi_pT >= 3 && xi_pT < 4 ) PtBin34->Fill( mass );
            if( xi_pT >= 4 && xi_pT < 5 ) PtBin45->Fill( mass );
            if( xi_pT >= 5 && xi_pT < 6 ) PtBin56->Fill( mass );

            if( xi_pT >= 1 && xi_pT < 3 ) PtBin13->Fill( mass );
            if( xi_pT >= 1 && xi_pT < 4 ) PtBin14->Fill( mass );
            if( xi_pT >= 1 && xi_pT < 5 ) PtBin15->Fill( mass );
            if( xi_pT >= 1 && xi_pT < 6 ) PtBin16->Fill( mass );

            if( xi_pT >= 2 && xi_pT < 4 ) PtBin24->Fill( mass );
            if( xi_pT >= 2 && xi_pT < 5 ) PtBin25->Fill( mass );
            if( xi_pT >= 2 && xi_pT < 6 ) PtBin26->Fill( mass );

            if( xi_pT >= 3 && xi_pT < 5 ) PtBin35->Fill( mass );
            if( xi_pT >= 3 && xi_pT < 6 ) PtBin36->Fill( mass );

            if( xi_pT >= 4 && xi_pT < 6 ) PtBin46->Fill( mass );
            
            // Xi eta and phi
            double xi_eta = xiCand->eta(  );
            double xi_phi = xiCand->phi(  );

            // Make vector of Xi Candidate parameters
            //
            TVector3 xiPEPvector;
            xiPEPvector.SetPtEtaPhi( xi_pT,xi_eta,xi_phi );
            pepVect_Xi->push_back( xiPEPvector );

        }

        // Make vectors for pairing of primary charged tracks
        //
        for( unsigned it = 0; it < tracks->size(  ); it++ )
        {
            const reco::Track & trk = ( *tracks )[it];
            math::XYZPoint bestvtx( bestvx, bestvy, bestvz );

            double dzvtx    = trk.dz( bestvtx );
            double dxyvtx   = trk.dxy( bestvtx );
            double dzerror  = sqrt( trk.dzError(  )*trk.dzError(  ) + bestvzError*bestvzError );
            double dxyerror = sqrt( trk.d0Error(  )*trk.d0Error(  ) + bestvxError*bestvyError );

            if( !trk.quality( reco::TrackBase::highPurity ) ) continue;
            if( fabs( trk.ptError(  )/trk.pt(  ) > 0.10 ) )   continue;
            if( fabs( dzvtx/dzerror ) > 3 )                   continue;
            if( fabs( dxyvtx/dxyerror ) > 3 )                 continue;

            double eta = trk.eta(  );
            double phi = trk.phi(  );
            double pt  = trk.pt(  );


            // Make vector of relevant trk parameters of primary charged particles
            //
            TVector3 primPEPvector;
            primPEPvector.SetPtEtaPhi( pt,eta,phi );
            // This vector is for producing the signal histogram pairing of two
            // charged primary tracks
            if( eta <= etaMax_trg_  && 
                eta >= etaMin_trg_ && 
                pt <= ptMax_trg_ &&
                pt >= ptMin_trg_ 
            ) 
            //contains same information as pepVect_trkass below but created just
            //for sake of clarity
            pepVect_trkhad->push_back( primPEPvector ); 

            // This vector is for the signal histogram pairing of charged
            // primary tracks and with Xi candidates
            if( eta <= etaMax_ass_  &&
                eta >= etaMin_ass_ &&
                pt <= ptMax_ass_ && //pT associated window
                pt >= ptMin_ass_ //pT associated window
            )
            pepVect_trkass->push_back( primPEPvector );

        }

        // Make signal histogram between Xi candidates and charged primary
        // tracks
        //
        int pepVect_Xi_size     = ( int )pepVect_Xi->size(  );
        int pepVect_trkass_size = ( int )pepVect_trkass->size(  );
        XiPerEvt->Fill( pepVect_Xi_size );
        TrkassPerEvt->Fill( pepVect_trkass_size );

        for( int xi_trg = 0; xi_trg < pepVect_Xi_size; xi_trg++ )
        {
            TVector3 pepVect_trg = (*pepVect_Xi)[xi_trg];
            double eta_trg       = pepVect_trg.Eta(  );
            double phi_trg       = pepVect_trg.Phi(  );

            for( int assoc = 0; assoc < pepVect_trkass_size; assoc++ )
            {
                TVector3 pepVect_ass = ( *pepVect_trkass )[assoc];
                double eta_ass       = pepVect_ass.Eta(  );
                double phi_ass       = pepVect_ass.Phi(  );

                double dEta = eta_ass - eta_trg;
                double dPhi = phi_ass - phi_trg;

                if( dPhi > PI )                    
                    dPhi=dPhi-2*PI;
                if( dPhi < -PI )                   
                    dPhi=dPhi+2*PI;
                if( dPhi > -PI && dPhi < -PI/2.0 ) 
                    dPhi=dPhi+2*PI;

                // To reduce jet fragmentation contributions
                if( fabs(dEta) < 0.028 && fabs( dPhi ) < 0.02  ) continue;
                SignalXi->Fill( dEta, dPhi, 1.0/pepVect_Xi_size );
            }
        }
        PepVect2_Xi->push_back( *pepVect_Xi );
        PepVect2_ass->push_back( *pepVect_trkass );
        zvtxVect->push_back( bestvz );

        delete pepVect_Xi;
        delete pepVect_trkass;


        // Make signal histogram for pairing of two charged primary tracks
        //
        int pepVect_trktrg_size = ( int )pepVect_trkhad->size(  );
        HadPerEvt->Fill( pepVect_trktrg_size );

        for( int trktrg1 = 0; trktrg1 < pepVect_trktrg_size; trktrg1++ )
        {
            TVector3 pepVect_trg = ( *pepVect_trkhad )[trktrg1];
            double eta_trg = pepVect_trg.Eta(  );
            double phi_trg = pepVect_trg.Phi(  );

            for( int trktrg2 = 0; trktrg2 < pepVect_trkass_size; trktrg2++ )
            {
                //if( trktrg2 == trktrg1 ) continue;
                TVector3 pepVect_ass = ( *pepVect_trkass )[trktrg2];
                double eta_ass = pepVect_ass.Eta(  );
                double phi_ass = pepVect_ass.Phi(  );

                double dEta = eta_ass - eta_trg;
                double dPhi = phi_ass - phi_trg;

                if( dPhi > PI )                    
                    dPhi=dPhi-2*PI;
                if( dPhi < -PI )                   
                    dPhi=dPhi+2*PI;
                if( dPhi > -PI && dPhi < -PI/2.0 ) 
                    dPhi=dPhi+2*PI;

                // To reduce jet fragmentation contributions
                if( fabs(dEta) < 0.028 && fabs( dPhi ) < 0.02  ) continue;
                SignalHad->Fill( dEta, dPhi, 1.0/pepVect_trktrg_size );
            }
        }

        PepVect2_had->push_back( *pepVect_trkhad );

        delete pepVect_trkhad;

    }
}


// ------------ method called once each job just before starting event loop  ------------
void 
XiCorrelation::beginJob()
{
    InvMassXi          = fs->make<TH1D>("InvMassXi", "InvMass #Xi", 150, 1.25, 1.40 );
    PtSpectra          = fs->make<TH1D>("pT_Xi", "pT Spectrum #Xi", 100, 0, 10 );
    nTrk               = fs->make<TH1D>( "nTrk", "nTrk", 250, 0, 250 );
    nEvtCut = fs->make<TH1D>( "nEvtCut", "nEvtCut", 10, 0, 10 );
    EtaPtCutnTrackHist = fs->make<TH1D>( "EtaPtCutnTrackHist", "EtaPtCutnTrack", 250, 0, 250 );
    XiPerEvt           = fs->make<TH1D>( "XiPerEvent", "#Xi per Event", 1500, 0, 1500 );
    HadPerEvt           = fs->make<TH1D>( "HadPerEvent", "Hadrons per Event", 1500, 0, 1500 );
    TrkassPerEvt   = fs->make<TH1D>( "TrkassPerEvent", "Associated trks per Event", 300,0, 300 );
    //TrktrgPerEvt   = fs->make<TH1D>( "TrktrgPerEvent", "Charged Trg per Event", 300, 0, 300 );
    BackgroundXi   = fs->make<TH2D>( "Background", "Bkg #Delta#eta;#Delta#phi", 33, -4.95, 4.95, 31, -( 0.5 - 1.0/32 )*PI, ( 1.5 - 1.0/32 )*PI );
    SignalXi       = fs->make<TH2D>( "Signal", "Sig #Delta#eta;#Delta#phi", 33, -4.95, 4.95, 31, -( 0.5 - 1.0/32 )*PI, ( 1.5 - 1.0/32 )*PI );
    BackgroundHad  = fs->make<TH2D>( "BackgroundHad", "BkgHad #Delta#eta;#Delta#phi", 33, -4.95, 4.95, 31, -( 0.5 - 1.0/32 )*PI, ( 1.5 - 1.0/32 )*PI );
    SignalHad      = fs->make<TH2D>( "SignalHad", "SigHad #Delta#eta;#Delta#phi", 33, -4.95, 4.95, 31, -( 0.5 - 1.0/32 )*PI, ( 1.5 - 1.0/32 )*PI );
    Correlation    = fs->make<TH2D>( "Correlation", "Correlation", 33, -4.95, 4.95, 31, -( 0.5 - 1.0/32 )*PI, ( 1.5 - 1.0/32 )*PI );
    CorrelationHad = fs->make<TH2D>( "CorrelationHad", "CorrelationHad", 33, -4.95, 4.95, 31, -( 0.5 - 1.0/32 )*PI, ( 1.5 - 1.0/32 )*PI );

    /*
    // For pT 1-3
    BackgroundXi13   = fs->make<TH2D>( "Background13", "Bkg #Delta#eta#Delta#phi 1 < pT < 3", 33, -4.95, 4.95, 31, -( 0.5 - 1.0/32 )*PI, ( 1.5 - 1.0/32 )*PI );
    SignalXi13       = fs->make<TH2D>( "Signal13", "Sig #Delta#eta#Delta#phi 1 < pT < 3", 33, -4.95, 4.95, 31, -( 0.5 - 1.0/32 )*PI, ( 1.5 - 1.0/32 )*PI );
    BackgroundHad13  = fs->make<TH2D>( "BackgroundHad13", "BkgHad #Delta#eta#Delta#phi 1 < pT < 3", 33, -4.95, 4.95, 31, -( 0.5 - 1.0/32 )*PI, ( 1.5 - 1.0/32 )*PI );
    SignalHad13      = fs->make<TH2D>( "SignalHad13", "SigHad #Delta#eta#Delta#phi 1 < pT < 3", 33, -4.95, 4.95, 31, -( 0.5 - 1.0/32 )*PI, ( 1.5 - 1.0/32 )*PI );
    Correlation13    = fs->make<TH2D>( "Correlation13", "Correlation 1 < pT < 3", 33, -4.95, 4.95, 31, -( 0.5 - 1.0/32 )*PI, ( 1.5 - 1.0/32 )*PI );
    CorrelationHad13 = fs->make<TH2D>( "CorrelationHad13", "CorrelationHad 1 < pT < 3", 33, -4.95, 4.95, 31, -( 0.5 - 1.0/32 )*PI, ( 1.5 - 1.0/32 )*PI );
    */

    PtBin1 = fs->make<TH1D>( "pT_1", "pT < 1 GeV #Xi", 150, 1.25, 1.40 );
    PtBin2 = fs->make<TH1D>( "pT_2", "pT < 2 GeV #Xi", 150, 1.25, 1.40 );
    PtBin3 = fs->make<TH1D>( "pT_3", "pT < 3 GeV #Xi", 150, 1.25, 1.40 );
    PtBin4 = fs->make<TH1D>( "pT_4", "pT < 4 GeV #Xi", 150, 1.25, 1.40 );
    PtBin5 = fs->make<TH1D>( "pT_5", "pT < 5 GeV #Xi", 150, 1.25, 1.40 );
    PtBin6 = fs->make<TH1D>( "pT_6", "pT < 6 GeV #Xi", 150, 1.25, 1.40 );
    PtBinOver6 = fs->make<TH1D>( "LargerThan6", "pT > 6 GeV #Xi", 150, 1.25, 1.40 );
    
    PtBin12 = fs->make<TH1D>( "pT_12", "pT 1-2 GeV #Xi", 150, 1.25, 1.40 );
    PtBin23 = fs->make<TH1D>( "pT_23", "pT 2-3 GeV #Xi", 150, 1.25, 1.40 );
    PtBin34 = fs->make<TH1D>( "pT_34", "pT 3-4 GeV #Xi", 150, 1.25, 1.40 );
    PtBin45 = fs->make<TH1D>( "pT_45", "pT 4-5 GeV #Xi", 150, 1.25, 1.40 );
    PtBin56 = fs->make<TH1D>( "pT_56", "pT 5-6 GeV #Xi", 150, 1.25, 1.40 );

    PtBin13 = fs->make<TH1D>( "pT_13", "pT 1-3 GeV #Xi", 150, 1.25, 1.40 );
    PtBin14 = fs->make<TH1D>( "pT_14", "pT 1-4 GeV #Xi", 150, 1.25, 1.40 );
    PtBin15 = fs->make<TH1D>( "pT_15", "pT 1-5 GeV #Xi", 150, 1.25, 1.40 );
    PtBin16 = fs->make<TH1D>( "pT_16", "pT 1-6 GeV #Xi", 150, 1.25, 1.40 );

    PtBin24 = fs->make<TH1D>( "pT_24", "pT 2-4 GeV #Xi", 150, 1.25, 1.40 );
    PtBin25 = fs->make<TH1D>( "pT_25", "pT 2-5 GeV #Xi", 150, 1.25, 1.40 );
    PtBin26 = fs->make<TH1D>( "pT_26", "pT 2-6 GeV #Xi", 150, 1.25, 1.40 );

    PtBin35 = fs->make<TH1D>( "pT_35", "pT 3-5 GeV #Xi", 150, 1.25, 1.40 );
    PtBin36 = fs->make<TH1D>( "pT_36", "pT 3-6 GeV #Xi", 150, 1.25, 1.40 );

    PtBin46 = fs->make<TH1D>( "pT_46", "pT 4-6 GeV #Xi", 150, 1.25, 1.40 );

    // For Background calculations which must be done in the endJob function
    //
    PepVect2_Xi        = new vector< vector<TVector3> >;
    //PepVect2_Primtrg = new vector< vector< TVector3 > >;
    PepVect2_ass       = new vector< vector<TVector3> >;
    PepVect2_had       = new vector< vector<TVector3> > ;
    zvtxVect           = new vector<double>;


}

// ------------ method called once each job just after ending the event loop  ------------
void 
XiCorrelation::endJob() 
{
    // Make background histograms
    // Xi paired with hadron
    int PepVect2_Xi_size = ( int )PepVect2_Xi->size(  );
    int PepVect2_ass_size = ( int )PepVect2_ass->size(  );

    for( int bkgnum = 0; bkgnum<bkgnum_; bkgnum++ )
    {
        int ncount = 0;
        for( int nevt_ass=0; nevt_ass<PepVect2_ass_size; nevt_ass++ )
        {
            int nevt_trg = gRandom->Integer( PepVect2_Xi_size );
            if( nevt_trg == nevt_ass )
            {
                nevt_ass--;
                continue;
            }
            if( fabs( ( *zvtxVect )[nevt_trg] - ( *zvtxVect )[nevt_ass] ) > 0.5 )
            {
                nevt_ass--;
                ncount++;
                if( ncount > 5000 )
                {
                    nevt_ass++;
                    ncount=0;
                }
                continue;
            }

            vector<TVector3> pepVectTmp_trg = ( *PepVect2_Xi )[nevt_trg];
            vector<TVector3> pepVectTmp_ass = ( *PepVect2_ass )[nevt_ass];
            int nMult_trg = pepVectTmp_trg.size(  );
            int nMult_ass = pepVectTmp_ass.size(  );

            for( int ntrg=0; ntrg<nMult_trg; ntrg++ )
            {
                TVector3 pvectorTmp_trg = pepVectTmp_trg[ntrg];
                double eta_trg = pvectorTmp_trg.Eta(  );
                double phi_trg = pvectorTmp_trg.Phi(  );

                for( int nass=0; nass<nMult_ass; nass++ )
                {
                    TVector3 pvectorTmp_ass = pepVectTmp_ass[nass];
                    double eta_ass = pvectorTmp_ass.Eta(  );
                    double phi_ass = pvectorTmp_ass.Phi(  );

                    double dEta = eta_ass - eta_trg;
                    double dPhi = phi_ass - phi_trg;

                    if( dPhi > PI )                    dPhi=dPhi-2*PI; // The bins range from -pi/2 to 3pi/2 so shouldnt it be dPhi > 3*PI/2.0?
                    if( dPhi < -PI )                   dPhi=dPhi+2*PI;
                    if( dPhi > -PI && dPhi < -PI/2.0 ) dPhi=dPhi+2*PI;

                    if( fabs(dPhi) < 0.028 && fabs( dEta ) < 0.02 ) continue;

                    BackgroundXi->Fill( dEta, dPhi, 1.0/nMult_trg );

                }
            }
        }
    }

    // hadron paired with hadron
    int PepVect2_had_size = ( int )PepVect2_had->size(  );

    for( int bkgnum = 0; bkgnum<bkgnum_; bkgnum++ )
    {
        int ncount = 0;
        for( int nevt_ass=0; nevt_ass<PepVect2_had_size; nevt_ass++ )
        {
            int nevt_trg = gRandom->Integer( PepVect2_had_size );
            if( nevt_trg == nevt_ass )
            {
                nevt_ass--;
                continue;
            }
            if( fabs( ( *zvtxVect )[nevt_trg] - ( *zvtxVect )[nevt_ass] ) > 0.5 )
            {
                nevt_ass--;
                ncount++;
                if( ncount > 5000 )
                {
                    nevt_ass++;
                    ncount=0;
                }
                continue;
            }

            vector<TVector3> pepVectTmp_trg = ( *PepVect2_had )[nevt_trg];
            vector<TVector3> pepVectTmp_ass = ( *PepVect2_ass )[nevt_ass];
            int nMult_trg = pepVectTmp_trg.size(  );
            int nMult_ass = pepVectTmp_ass.size(  );

            for( int ntrg=0; ntrg<nMult_trg; ntrg++ )
            {
                TVector3 pvectorTmp_trg = pepVectTmp_trg[ntrg];
                double eta_trg = pvectorTmp_trg.Eta(  );
                double phi_trg = pvectorTmp_trg.Phi(  );

                for( int nass=0; nass<nMult_ass; nass++ )
                {
                    TVector3 pvectorTmp_ass = pepVectTmp_ass[nass];
                    double eta_ass = pvectorTmp_ass.Eta(  );
                    double phi_ass = pvectorTmp_ass.Phi(  );

                    double dEta = eta_ass - eta_trg;
                    double dPhi = phi_ass - phi_trg;

                    if( dPhi > PI )                    dPhi=dPhi-2*PI; // The bins range from -pi/2 to 3pi/2 so shouldnt it be dPhi > 3*PI/2.0?
                    if( dPhi < -PI )                   dPhi=dPhi+2*PI;
                    if( dPhi > -PI && dPhi < -PI/2.0 ) dPhi=dPhi+2*PI;

                    if( fabs(dPhi) < 0.028 && fabs( dEta ) < 0.02 ) continue;

                    BackgroundHad->Fill( dEta, dPhi, 1.0/nMult_trg );

                }
            }
        }
    }

    int nEvent = XiPerEvt->Integral( 3, 10000 ); // XiPerEvt bin 1 and 2 contain the most so why are we starting from 3 instead of 0?
    //double bz = BackgroundXi->GetBinContent( 275 );
    Correlation->Add( SignalXi );
    //Correlation->Scale( bz );
    Correlation->Divide( BackgroundXi );
    Correlation->Scale( 1.0/nEvent );

    int nEventbkg = HadPerEvt->Integral( 3, 10000 ); // XiPerEvt bin 1 and 2 contain the most so why are we starting from 3 instead of 0?
    //double bz = BackgroundXi->GetBinContent( 275 );
    CorrelationHad->Add( SignalHad );
    //Correlation->Scale( bz );
    CorrelationHad->Divide( BackgroundHad );
    CorrelationHad->Scale( 1.0/nEventbkg );
}

