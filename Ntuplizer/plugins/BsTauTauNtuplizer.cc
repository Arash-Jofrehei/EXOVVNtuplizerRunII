#include "../interface/BsTauTauNtuplizer.h"


//===================================================================================================================
BsTauTauNtuplizer::BsTauTauNtuplizer( edm::EDGetTokenT<reco::MuonCollection>  muonToken   ,
				      edm::EDGetTokenT<edm::SortedCollection<CaloTower>> CaloTowerCollection,
              edm::EDGetTokenT<reco::Centrality> tok_centSrc,
				      edm::EDGetTokenT<reco::VertexCollection> verticeToken, 
				      edm::EDGetTokenT<std::vector<reco::PFCandidate>> packedpfcandidatesToken,
				      //edm::EDGetTokenT<std::vector<reco::TrackCollection>> TrackCollectionToken,
				      edm::EDGetTokenT<edm::TriggerResults> triggertoken,
				      edm::EDGetTokenT<reco::GenParticleCollection> genptoken,
				      edm::EDGetTokenT<std::vector<reco::GenJet>> genttoken,
				      std::map< std::string, bool >& runFlags,
				      std::map< std::string, double >& runValues,
				      std::map< std::string, std::string >& runStrings,
				      NtupleBranches* nBranches )
    : CandidateNtuplizer ( nBranches )
    , muonToken_	        ( muonToken )
    , CaloTowerCollection_ ( CaloTowerCollection )
    , tok_centSrc_ ( tok_centSrc )
    , verticeToken_          ( verticeToken )
    , packedpfcandidatesToken_(packedpfcandidatesToken) 
    //, TrackCollectionToken_(TrackCollectionToken) 
    , HLTtriggersToken_	( triggertoken )
    , genParticlesToken_( genptoken )
    , genTauToken_( genttoken )
    , runOnMC_   (runFlags["runOnMC"])
    , c_dz (runValues["dzcut"])
    , c_fsig (runValues["fsigcut"])
    , c_vprob (runValues["vprobcut"])
   
{

  std::cout << "-- (dzcut, fsigcut, vprobcut) = " << c_dz << " " << c_fsig << " " << c_vprob << std::endl;
  
}

//===================================================================================================================
BsTauTauNtuplizer::~BsTauTauNtuplizer( void )
{

}


bool BsTauTauNtuplizer::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){

  
     //if (event.id().event()== 73) std::cout << "---------------- event, run, lumi = " << event.id().event() << " " << event.id().run() << " " << event.id().luminosityBlock() << "----------------" << std::endl;
      
    // MC Signal Madgraph 2015 mu+3prong and triggered and isRight
    //if (event.id().run()!=1 || event.id().luminosityBlock()!=996 || event.id().event()!=2890) return false;
    
    bool isMC = runOnMC_;
    bool mu_3prong_event = false;
    bool tau_mu_check = false;
    bool tau_3prong_check = false;

    if(isMC){
      event.getByToken(genParticlesToken_ , genParticles_); 
      event.getByToken(genTauToken_, genTaus_);
      
      for( unsigned p=0; p < genParticles_->size(); ++p){
        
        if(TMath::Abs((*genParticles_)[p].pdgId())!=15) continue;
        //std::cout << "pdgId - status of the gen tau: " << (*genParticles_)[p].pdgId() << " - " << (*genParticles_)[p].status() << std::endl;
        if(TMath::Abs((*genParticles_)[p].status())!=2) continue;
        /*std::cout << "mother tau: " << (*genParticles_)[p].pdgId() << " (pT, eta, phi) = " 
  		  << (*genParticles_)[p].pt() << " " 
  		  << (*genParticles_)[p].eta() << " " 
  		  << (*genParticles_)[p].phi() << std::endl;*/
        nBranches_->gen_tau_pt.push_back((*genParticles_)[p].pt());
        nBranches_->tau_pT->Fill((*genParticles_)[p].pt());
        nBranches_->gen_tau_eta.push_back((*genParticles_)[p].eta());
        nBranches_->gen_tau_phi.push_back((*genParticles_)[p].phi());
        
        //std::cout << "\t Tau found with # of daughters = " << (*genParticles_)[p].numberOfDaughters() << " with mother = " << (*genParticles_)[p].mother(0)->pdgId() << std::endl;
        
        int neutral_pion_counter = 0;
        int charged_pion_counter = 0;
        int a1_counter = 0;
        int rho_counter = 0;
        int muon_counter = 0;
        int electron_counter = 0;
        int charged_W_counter = 0;
        int gamma_counter = 0;
        
        for(int idd = 0; idd < (int)(*genParticles_)[p].numberOfDaughters(); idd++){ // loop on daughters of gen tau
  	
  	/*std::cout << "\t\t -> " << (*genParticles_)[p].daughter(idd)->pdgId() << " (pT, eta, phi) = " 
  		  << (*genParticles_)[p].daughter(idd)->pt() << " " 
  		  << (*genParticles_)[p].daughter(idd)->eta() << " " 
  		  << (*genParticles_)[p].daughter(idd)->phi() << std::endl;*/
  
          if (TMath::Abs((*genParticles_)[p].daughter(idd)->pdgId())==111) neutral_pion_counter += 1;
          if (TMath::Abs((*genParticles_)[p].daughter(idd)->pdgId())==211) charged_pion_counter += 1;
          if (TMath::Abs((*genParticles_)[p].daughter(idd)->pdgId())==13) muon_counter += 1;
          if (TMath::Abs((*genParticles_)[p].daughter(idd)->pdgId())==11) electron_counter += 1;
          if (TMath::Abs((*genParticles_)[p].daughter(idd)->pdgId())==20213) a1_counter += 1;
          if (TMath::Abs((*genParticles_)[p].daughter(idd)->pdgId())==100213) rho_counter += 1;
          if (TMath::Abs((*genParticles_)[p].daughter(idd)->pdgId())==24) charged_W_counter += 1;
          if (TMath::Abs((*genParticles_)[p].daughter(idd)->pdgId())==21) gamma_counter += 1;
          nBranches_->gen_tau_daughter_pdgId.push_back((*genParticles_)[p].daughter(idd)->pdgId());
          nBranches_->gen_tau_daughter_pt.push_back((*genParticles_)[p].daughter(idd)->pt());
          nBranches_->gen_tau_daughter_eta.push_back((*genParticles_)[p].daughter(idd)->eta());
          nBranches_->gen_tau_daughter_phi.push_back((*genParticles_)[p].daughter(idd)->phi());
          if(
  	        TMath::Abs((*genParticles_)[p].daughter(idd)->pdgId())==12 ||
  	        TMath::Abs((*genParticles_)[p].daughter(idd)->pdgId())==14 || 
  	        TMath::Abs((*genParticles_)[p].daughter(idd)->pdgId())==16
  	        ) continue;
  
        } //loop on daughters of gen tau
  
        if (charged_pion_counter == 3) tau_3prong_check = true;
        if (muon_counter) tau_mu_check = true;
        //cout << "number of pions: " << charged_pion_counter << endl;
        
        //tau decay modes
        nBranches_->tau_decays->Fill(0);
        if (muon_counter == 1) nBranches_->tau_decays->Fill(1);
        else if (electron_counter == 1) nBranches_->tau_decays->Fill(2);
        else if (charged_pion_counter == 1){
          if (neutral_pion_counter == 0) nBranches_->tau_decays->Fill(3);
          if (neutral_pion_counter == 1) nBranches_->tau_decays->Fill(4);
          if (neutral_pion_counter > 1) nBranches_->tau_decays->Fill(5);
        }
        else if (charged_pion_counter == 3) nBranches_->tau_decays->Fill(6);
        else if (a1_counter == 1) nBranches_->tau_decays->Fill(7);
        else if (rho_counter == 1) nBranches_->tau_decays->Fill(8);
        else if (charged_W_counter == 1) nBranches_->tau_decays->Fill(9);
        else if (gamma_counter >= 1) nBranches_->tau_decays->Fill(10);
        else nBranches_->tau_decays->Fill(11);
        
        
      } // loop on gen particles in the event
      
      if (tau_mu_check && tau_3prong_check) mu_3prong_event = true; // Arash
      if (tau_mu_check) nBranches_->cutflow_perevt->Fill(5);
      if (tau_3prong_check) nBranches_->cutflow_perevt->Fill(6);
      if (mu_3prong_event) nBranches_->cutflow_perevt->Fill(7);
      
      nBranches_->gen_tau_to_mu.push_back(tau_mu_check);
      nBranches_->gen_tau_to_3prong.push_back(tau_3prong_check);
      nBranches_->gen_tautau_to_mu3prong.push_back(mu_3prong_event);
      //cout << "    ---->    Event with mu + 3prong tau: " << mu_3prong_event << endl;
    }
    
    
    
    /********************************************************************
     *
     * Step1: check if the J/psi trigger is fired.
     * Namely, HLT_DoubleMu4_JpsiTrk_Displaced_v
     * and  HLT_Dimuon0_Jpsi3p5_Muon2_v
     ********************************************************************/
    event.getByToken(HLTtriggersToken_, HLTtriggers_);
    nBranches_->cutflow_perevt->Fill(0);
    
    // bool isTriggered = false;
    const edm::TriggerNames& trigNames = event.triggerNames(*HLTtriggers_);
    std::vector<std::string> finalTriggerName;
    //    std::string finalTriggerFilterObjName="";

    bool trigMatch = false;
  
    for (unsigned int i = 0, n = HLTtriggers_->size(); i < n; ++i) {

        // if(trigNames.triggerName(i).find("HLT_DoubleMu4_JpsiTrk_Displaced_v")!= std::string::npos || trigNames.triggerName(i).find("HLT_Dimuon0_Jpsi3p5_Muon2_v")!= std::string::npos ){
           
        //     nBranches_->HLT_BPH_isFired[trigNames.triggerName(i)] = HLTtriggers_->accept(i);
// 2018 trigger
      if(trigNames.triggerName(i).find("HLT_HIUPC_SingleMuOpen_NotMBHF2AND")!= std::string::npos){
// 2015 trigger
      //if(trigNames.triggerName(i).find("HLT_HIUPC")!= std::string::npos)cout << trigNames.triggerName(i);
//      if(trigNames.triggerName(i).find("HLT_HIUPCSingleMuNotHF2Pixel_SingleTrack_v")!= std::string::npos){
        //cout << "    is found.";
        nBranches_->HLT_isFired[trigNames.triggerName(i)] = HLTtriggers_->accept(i);
        if(HLTtriggers_->accept(i)){
          //cout << " and accepted.";
          trigMatch = true;
//          std::cout << trigNames.triggerName(i) << std::endl;
          //	  isTriggered = true;
          //	  std::cout << "This trigger is fired:" << trigNames.triggerName(i) << std::endl;
          finalTriggerName.push_back(trigNames.triggerName(i));
          //                finalTriggerFilterObjName="hltJpsiTkVertexFilter";
//           std::cout << "finalTriggerName = "  << finalTriggerName << std::endl;
                
        }
        //cout << endl;
      }
      //else if(trigNames.triggerName(i).find("HLT_HIUPC")!= std::string::npos) cout << endl;
    }

    //if(!trigMatch) return false;
    //nBranches_->cutflow_perevt->Fill(1);
    if(trigMatch) nBranches_->cutflow_perevt->Fill(1);
    nBranches_->triggered.push_back(trigMatch);

    /********************************************************************
     *
     * Step2: pre-select muons for building J/psi candidates ... 
     * For muons, no requirement applied
     *
     ********************************************************************/

    event.getByToken(verticeToken_   , vertices_     );
    event.getByToken(muonToken_	, muons_    );
    //    event.getByToken(triggerObjects_  , triggerObjects);

    std::vector<reco::Muon> muoncollection;
    muoncollection.clear();


     //   std::cout << "#muons: " << muons_->size() << std::endl;
    // evt Triggered
    //    for (std::vector<reco::Muon>::const_iterator muon = muons_->begin(); muon != muons_->end(); ++muon) {
    for(size_t imuon = 0; imuon < muons_->size(); ++ imuon){

        const reco::Muon & muon = (*muons_)[imuon];
        if( (muon.pt() < 3.5 && TMath::Abs(muon.eta()) < 1.2) || (muon.pt() < 2.5 && TMath::Abs(muon.eta()) > 1.2))  continue;
        if(TMath::Abs(muon.eta()) > 2.4) continue;
        if(!(muon.track().isNonnull())) continue;
        //    bool isSoft = muon.isSoftMuon(*firstGoodVertex);
        //    bool isGlobal = muon.isGlobalMuon();
        //    bool isTracker = muon.isTrackerMuon();
        //    bool isLoose = muon.isLooseMuon();
        //    bool isTight =  muon.isTightMuon(*firstGoodVertex);
        //    bool isPF = muon.isPFMuon();
        //    if(!(isSoft && isGlobal)) continue;
        //    if(TMath::Abs(muon.muonBestTrack()->dz(firstGoodVertex->position())) > 0.5) continue;

    
        // Trigger matching

	//	bool trigMatch = false;

//        for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
//    
//            obj.unpackPathNames(trigNames);
//            obj.unpackFilterLabels(event, *HLTtriggers_);
//
//            std::vector<std::string> pathNamesAll  = obj.pathNames(false);
//
//            bool isPathExist = false;
//
//            for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {
//	      
//      	      for(int iname=0; iname < (int)finalTriggerName.size(); iname ++){
//            		if(pathNamesAll[h]==finalTriggerName[iname]) isPathExist = true;
//      	      }
//            }
//      
//            if(!isPathExist) continue;
//
////            bool isFilterExist = false;
////    
////            for (unsigned hh = 0; hh < obj.filterLabels().size(); ++hh){
////	
////                if(obj.filterLabels()[hh].find(finalTriggerFilterObjName) != std::string::npos){
////                    isFilterExist = true;
////                }
////            }
////      
////            if(!isFilterExist) continue;
//      
////            Float_t trigger_dR = reco::deltaR(obj.eta(), obj.phi(),
////                                              muon.eta(), muon.phi());
//      
//	    //	    if(trigger_dR < 0.1) trigMatch = true;
//        }

	//	if(!trigMatch) continue;

	//	std::cout << "imuon = " << imuon << " "  << muon.pt() << std::endl;

        muoncollection.push_back(muon);
    }


    //    std::cout << "number of matched muon = " << muoncollection.size() << std::endl;
    if(!( muoncollection.size() >= 1)){
      //if (mu_3prong_event) return true;
      //else return false;
      return false;
    }
    nBranches_->cutflow_perevt->Fill(2);
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);
    
    //    std::cout << "number of matched muon = " << muoncollection.size() << std::endl;

    const reco::TrackRef track_muon = muoncollection[0].muonBestTrack();
    reco::TransientTrack tt_muon = (*builder).build(track_muon);

    KinematicParticleFactoryFromTransientTrack pFactory;
    KinematicParticleVertexFitter kpvFitter;

    edm::ESHandle<MagneticField> fieldHandle;
    iSetup.get<IdealMagneticFieldRecord>().get(fieldHandle);
    fMagneticField = fieldHandle.product();

    AnalyticalImpactPointExtrapolator extrapolator(fMagneticField);

    reco::Vertex closestVertex; 
    closestVertex = *(vertices_->begin());
    Float_t max_dz = 999;

    for( reco::VertexCollection::const_iterator vtx = vertices_->begin(); vtx != vertices_->end(); ++vtx){
      
      Float_t _dz_ = TMath::Abs(vtx->position().z() - muoncollection[0].vz());
      
	//      bool isFake = (vtx->chi2()==0 && vtx->ndof()==0);
	//      if( !isFake && vtx->ndof()>=4. && vtx->position().Rho()<=2.0 && fabs(vtx->position().Z())<=24.0) {
//	firstGoodVertex = vtx;
//	break;
//      }
      if(_dz_ < max_dz){
      	max_dz = _dz_;
      	closestVertex = *vtx;
      }
    }

    /********************************************************************
     *
     * Step6: Tau selection
     *        Just select highest in pT but there might be better selection ... 
     *
     ********************************************************************/
    
    event.getByToken( packedpfcandidatesToken_               , packedpfcandidates_      ); 
    //event.getByToken( TrackCollectionToken_               , TrackCollection_      ); 
    
    std::vector<reco::PFCandidate> pfcollection; 
    //    std::vector<pat::PackedCandidate> pfcollection_pre; 
    std::vector<reco::TransientTrack> mytracks;
    
//    for( size_t ii = 0; ii < packedpfcandidates_->size() && ii < 2; ++ii ){   
//	
//      pat::PackedCandidate pf = (*packedpfcandidates_)[ii];
//	
//      pfcollection_pre.push_back(pf);
//    }
    
    //Handle<reco::TrackCollection> tracks;
    //iEvent.getByLabel("generalTracks",tracks);   
    

    Int_t npf_qr = 0;

    nBranches_->PVz.push_back(closestVertex.position().z());

    bool print =  mu_3prong_event;
    //std::cout << "                            ---------------- event, run, lumi = " << event.id().event() << " " << event.id().run() << " " << event.id().luminosityBlock() << "----------------" << std::endl;
    
    if (print) std::cout << "---------------- event, run, lumi = " << event.id().event() << " " << event.id().run() << " " << event.id().luminosityBlock() << "----------------" << std::endl;
    
    
    //for(reco::TrackCollection::const_iterator itTrack = TrackCollection_->begin(); itTrack != TrackCollection_->end(); ++itTrack) {
    //if (print) cout << TrackCollection_->size() << endl;
    //for( size_t ii = 0; ii < TrackCollection_->size(); ++ii ){
    //      reco::TrackCollection track = (*TrackCollection_)[ii];
    //      if (print) cout << track.size() << endl;
    //      if (print) cout << "track pt: " << track[0].pt() << "  eta: " << track[0].eta() << " phi: " << track[0].phi() << endl;
    //}
    
    
    int nPions = 0;
    for( size_t ii = 0; ii < packedpfcandidates_->size(); ++ii ){   

      if (packedpfcandidates_->size() > 200.) continue;
      
      reco::PFCandidate pf = (*packedpfcandidates_)[ii];
      
      if (print) std::cout << "pf pdgID, pT , eta, phi = " << pf.pdgId() << " -> " << pf.pt() << " -> " << pf.eta() << " -> " << pf.phi() << std::endl;
      
      if(TMath::Abs(pf.pdgId())!=211) {if (print) cout << "Not a charged pion. pdgID: " << pf.pdgId() << endl; continue;}
      else {if (print) cout << "Charged pion. pdgID: " << pf.pdgId() << endl;}
      
      if(pf.pt() < 0.5) {if (print) cout << "Pion pT less than 0.5 GeV" << endl; /*continue;*/}

      nPions += 1;
      
      /*bool hasTrackDetails = false;
      auto tk = pf.bestTrack();
      if(tk) {
        int numberOfPixelLayers = tk->hitPattern().pixelLayersWithMeasurement();
        int numberOfStripLayers = tk->hitPattern().stripLayersWithMeasurement();
        if (print) cout << "Number of pixel & strip layers: " << numberOfPixelLayers << " - " << numberOfStripLayers << endl;
        int numberOfPixelHits = tk->hitPattern().numberOfValidPixelHits() - numberOfPixelLayers;
        int numberOfStripHits = tk->hitPattern().numberOfValidHits() - numberOfPixelHits - numberOfPixelLayers - numberOfStripLayers;
        hasTrackDetails = numberOfPixelLayers || numberOfStripLayers || numberOfPixelHits || numberOfStripHits;
        if (print) cout << "Number of pixel & strip hits: " << numberOfPixelHits << " - " << numberOfStripHits << endl;
      }
      if (print) cout << "hasTrackDetails: " << hasTrackDetails << endl;*/

      /**************** ARASH **********************************/
      //      if(!pf.hasTrackDetails()) continue;
            //if(pf.trackRef().isNonnull() && pf.trackRef().isAvailable()){
              //if (!pf.trackRef()->TrackBase::numberOfValidHits()) {cout << "Track has no hits" << endl; continue;}
              //if (pf.bestTrack()) cout << "valid hits in best track: " << pf.bestTrack()->hitPattern().numberOfValidHits() << endl; else cout << "No best track" << endl;
              //if (!pf.trackRef()->hasTrackDetails()) {cout << "Track has no details." << endl; continue;}
            //}
      /**************** ARASH **********************************/
      
      
      if(TMath::Abs(pf.eta()) > 2.5) {if (print) cout << "eta outside +- 2.5" << endl; continue;}
      
      Float_t precut_dz = pf.vz() - closestVertex.position().z();
      //if(TMath::Abs(precut_dz) > c_dz) {if (print) cout << "dz precut not satisfied." << endl; continue;}
      if (print) cout << "pf.z: " << pf.vz() << " PV.z: " << closestVertex.position().z() << " dz: " << precut_dz << endl;
      if(TMath::Abs(precut_dz) > c_dz) {/*continue;*/}
      //        TLorentzVector temp;
      //        temp.SetPtEtaPhiM(pf.pt(), pf.eta(), pf.phi(), 0.13957018);
      //        std::cout << " - - - - - - - > " << temp.Z() << " " << pf.vz() << " " << std::endl;
      //        std::cout << " - - - - - > " << temp.Z() - closestVertex.position().z() << " " << std::endl;
      //        std::cout << " - - - - - > " << pf.vz()  - closestVertex.position().z() << " " << std::endl;
      
      /**************** ARASH **********************************/
            bool hpflag = false;
            if(pf.trackRef().isNonnull() && pf.trackRef().isAvailable()) hpflag = pf.trackRef()->quality(reco::TrackBase::highPurity);
            if(!hpflag) {cout << " Impure track " << endl; continue;}
      /**************** ARASH **********************************/

      /**************** ARASH **********************************/
      
            //if(pf.bestTrack()->hitPattern().numberOfValidPixelHits() < 0) continue;
            //cout <<  pf.bestTrack()->hitPattern() << " pixel hits" << endl;
            if(pf.bestTrack()->hitPattern().numberOfValidHits() < 3) continue;
            if(pf.bestTrack()->normalizedChi2() > 100) continue;
            if (print) cout << "checking ratio: " << pf.bestTrack() ->chi2() / pf.bestTrack()->ndof() << endl;
      /**************** ARASH **********************************/
      
      if (print) cout << "ECAL: " << pf.ecalEnergy() << " --> HCAL: " << pf.hcalEnergy() << endl;


      //      Float_t _dR1 = reco::deltaR(pf.eta(), pf.phi(), 
      //				  mu1_fit->eta(), mu1_fit->phi());
      //      
      //      Float_t _dR2 = reco::deltaR(pf.eta(), pf.phi(), 
      //				  mu2_fit->eta(), mu2_fit->phi());
      //
      //
      //      if(_dR1 < 0.1 || _dR2 < 0.1){
      //	if(TMath::Abs(pf.pdgId()) == 13) continue;
      //      }
      
      nBranches_->pion_z.push_back(pf.vz());
      nBranches_->reco_pion_pt.push_back(pf.pt());
      nBranches_->reco_pion_eta.push_back(pf.eta());
      nBranches_->reco_pion_phi.push_back(pf.phi());
      if (print) cout << "Pion added to tau." << endl;
      //if (event.id().event()== 112) cout << "pf pt: " << pf.pt() << "  eta: " << pf.eta() << "  phi: " << pf.phi() << "  pdgId: " << pf.pdgId() << endl;
      pfcollection.push_back(pf);
      reco::TransientTrack  tt_track = (*builder).build(pf.bestTrack());
      mytracks.push_back(tt_track);
      
    }


    // cut on tracks 500 MeV
    for( size_t ii = 0; ii < packedpfcandidates_->size(); ++ii ){   

      if (packedpfcandidates_->size() > 200.) continue;

    	reco::PFCandidate pf = (*packedpfcandidates_)[ii];
    	
	/**************** ARASH **********************************/
	//    	if(pf.pt() > 0.5 && pf.trackHighPurity()) { 
	    	if(pf.pt() > 0.5 && pf.trackRef().isNonnull() && pf.trackRef().isAvailable()) { 
	/**************** ARASH **********************************/

    	if(pf.trackRef()->quality(reco::TrackBase::highPurity)) { 

        nBranches_->BsTauTau_trackPFactivity_pt .push_back(pf.pt());
        nBranches_->BsTauTau_trackPFactivity_eta.push_back(pf.eta()); 
        nBranches_->BsTauTau_trackPFactivity_phi.push_back(pf.phi());
      }
        }

    }

    // ******************************
    // Calo and HF 

    edm::Handle<edm::SortedCollection<CaloTower>> CaloTowerHandle;
    event.getByToken(CaloTowerCollection_, CaloTowerHandle);
 
    double eta;
    double phi;
    double energy;
 
    for (edm::SortedCollection<CaloTower>::const_iterator calo = CaloTowerHandle->begin(); calo != CaloTowerHandle->end(); ++calo) {
      
      energy=calo->energy();

      phi=calo->phi();
      eta=calo->eta();

      nBranches_->BsTauTau_calo_eta.push_back(eta);
      nBranches_->BsTauTau_calo_phi.push_back(phi);
      nBranches_->BsTauTau_calo_energy.push_back(energy);
 
      ECaloType subdetHad = GetTowerSubdetHad(eta);
//      if (TMath::Abs(eta)>3. && TMath::Abs(eta)<5.2) {
//        std::cout << " - - - - - - - - - >> " << energy << " " << eta << " " << subdetHad << " " << kHFp << " " << kHFm << std::endl;
//      }
      if(subdetHad==kHFp) nBranches_->BsTauTau_calo_energyHFp.push_back(energy);
      else nBranches_->BsTauTau_calo_energyHFp.push_back(-1);
      if(subdetHad==kHFm) nBranches_->BsTauTau_calo_energyHFm.push_back(energy);
      else nBranches_->BsTauTau_calo_energyHFm.push_back(-1);
    }
    
    // ZDC
    event.getByToken(tok_centSrc_, cent);
    
    float ZDCPlus  = (cent.isValid() ? cent->zdcSumPlus()  : -1.);
    float ZDCMinus = (cent.isValid() ? cent->zdcSumMinus() : -1.);
    
    nBranches_->BsTauTau_calo_zdcSumPlus.push_back(ZDCPlus);
    nBranches_->BsTauTau_calo_zdcSumMinus.push_back(ZDCMinus);
  

    // ******************************
    // ******************************

    
    nBranches_->BsTauTau_nPions.push_back(nPions);
    Int_t numOfch = (size_t)pfcollection.size();
    if (print) cout << "nPions: " << nPions << endl;
    if (print) cout << "nch: " << numOfch << endl;
    
    
        
    // retrieve gen. information

    std::vector<std::vector<TLorentzVector>> gps;
    std::vector<Int_t> ppdgId;
    std::vector<Int_t> vec_gentaudm;
    std::vector<Int_t> vec_ppdgId;
    std::vector<TLorentzVector> vec_gentaup4;
    std::vector<TLorentzVector> vec_gentau3pp4;
    Int_t isgen3 = 0;
    Int_t isgen3matched = 0;
  
    //bool isMC = runOnMC_;

    if(isMC){
      event.getByToken(genParticlesToken_ , genParticles_); 
      event.getByToken(genTauToken_, genTaus_);
      //bool mu_3prong_event = false;
      //bool tau_mu_check = false;
      //bool tau_3prong_check = false;
      for( unsigned p=0; p < genParticles_->size(); ++p){
        
        if(TMath::Abs((*genParticles_)[p].pdgId())!=15) continue;
        //if((*genParticles_)[p].pdgId()!=15) continue;
        if(TMath::Abs((*genParticles_)[p].status())!=2) continue;
        
        //std::cout << "\t Tau found with # of daughters = " << (*genParticles_)[p].numberOfDaughters() << " with mother = " << (*genParticles_)[p].mother(0)->pdgId() << std::endl;
        
        // calculate visible pt ... 
  
        TLorentzVector genvis;
        std::vector<TLorentzVector> gp;
        std::vector<TLorentzVector> gp_charged_pion;
        //int charged_pion_counter = 0;
        Bool_t matched = true;
        Int_t nprong = 0;
        
        for(int idd = 0; idd < (int)(*genParticles_)[p].numberOfDaughters(); idd++){ // loop on daughters of gen tau
  	
  	if (print) std::cout << "\t\t -> " << (*genParticles_)[p].daughter(idd)->pdgId() << " (pT, eta, phi, z) = " 
  		  << (*genParticles_)[p].daughter(idd)->pt() << " " 
  		  << (*genParticles_)[p].daughter(idd)->eta() << " " 
  		  << (*genParticles_)[p].daughter(idd)->phi() << " " 
  		  << (*genParticles_)[p].daughter(idd)->vz() << std::endl;
  
          //if (TMath::Abs((*genParticles_)[p].daughter(idd)->pdgId())==211) charged_pion_counter += 1;
          //if (TMath::Abs((*genParticles_)[p].daughter(idd)->pdgId())==13) tau_mu_check = true;
          if(
  	        TMath::Abs((*genParticles_)[p].daughter(idd)->pdgId())==12 ||
  	        TMath::Abs((*genParticles_)[p].daughter(idd)->pdgId())==14 || 
  	        TMath::Abs((*genParticles_)[p].daughter(idd)->pdgId())==16
  	        ) continue;
  
  
         	TLorentzVector _genvis_;
        	_genvis_.SetPtEtaPhiM((*genParticles_)[p].daughter(idd)->pt(),
  			      (*genParticles_)[p].daughter(idd)->eta(),
  			      (*genParticles_)[p].daughter(idd)->phi(),
  			      (*genParticles_)[p].daughter(idd)->mass()
  			      );
  	
        	genvis += _genvis_;
        
        	if(TMath::Abs((*genParticles_)[p].daughter(idd)->pdgId())==211){
  
        	  nprong += 1;
        	  
        	  // check matching to reco PF objects
        	  Float_t min_dr = 999;
  	  
        	  for(int kkk = 0; kkk < numOfch; kkk ++){
        	    
        	    reco::PFCandidate _pf = pfcollection[kkk];
        	    
        	    if(_pf.pdgId()!=(*genParticles_)[p].daughter(idd)->pdgId()) continue;
        	    
        	    Float_t _dR = reco::deltaR(
        				       _genvis_.Eta(), _genvis_.Phi(),
        				       _pf.eta(), _pf.phi()
        				       );
        	    if(_dR < min_dr && _dR < 0.015 && _pf.pt()/_genvis_.Pt() < 1.15 && _pf.pt()/_genvis_.Pt() > 0.85){
        	      min_dr = _dR;
        	    }
        	  }
  	  
  //	  Float_t min_dr2 = 999;
  //	  for( size_t iii = 0; iii < packedpfcandidates_->size(); ++iii ){   
  //      
  //	    pat::PackedCandidate pf = (*packedpfcandidates_)[iii];
  //	    
  //	    if(pf.pdgId()!=(*genParticles_)[p].daughter(idd)->pdgId()) continue;
  //	    
  //	    Float_t _dR = reco::deltaR(
  //				       _genvis_.Eta(), _genvis_.Phi(),
  //				       pf.eta(), pf.phi()
  //				       );
  //	    if(_dR < min_dr2 && _dR < 0.1){
  //	      min_dr2 = _dR;
  //	    }
  //	  }
  
  	  //	  if(min_dr2!=999) std::cout << "pf matched !!!" << std::endl;
  
  
  	  //////////////////////////////////////
  	  if(min_dr == 999) matched = false;
  	  //	  else std::cout << "matched!" << std::endl;
  	  //	  else{
  	  gp.push_back(_genvis_); // Arash: gp is a vector. Each element (_genvis_) is the Lorentz vector of a non-neutrino daughter of a gen tau.
  	  //	  }
  	}
        } //loop on daughters of gen tau
  
        //if (charged_pion_counter == 3) tau_3prong_check = true;
        if(nprong==3) isgen3 += 1;
  
        //////////////////////////
        // check decay mod. To do this, take matching with tau-genjet. 
        //////////////////////////
        Float_t min_gendr = 999;
        Int_t taugendm = -999;
  
        for(size_t i = 0; i < genTaus_->size(); ++ i){      
  	
        	const reco::GenJet & TauCand = (*genTaus_)[i];
        	
        	reco::Particle::LorentzVector visibleP4 = ((*genTaus_)[i]).p4();
        
        	TLorentzVector visp4;
        	visp4.SetPtEtaPhiM(visibleP4.pt(),
  		    visibleP4.eta(),
   		    visibleP4.phi(),
   		    visibleP4.mass());
  	
  	      Float_t dRgen = genvis.DeltaR(visp4);
  	
        	if(dRgen < min_gendr && dRgen < 0.1){
        	  min_gendr = dRgen;
        	  taugendm = aux.decaymode_id(JetMCTagUtils::genTauDecayMode(TauCand));
        	}
        }
        
        vec_ppdgId.push_back((*genParticles_)[p].mother(0)->pdgId());
        vec_gentaudm.push_back(taugendm);
        vec_gentaup4.push_back(genvis);
  
        if(gp.size()==3){
        	//std::cout << "\t -----> This has been registered with mother = " << (*genParticles_)[p].mother(0)->pdgId() << std::endl;
        	gps.push_back(gp);
        	ppdgId.push_back((*genParticles_)[p].mother(0)->pdgId());
        	vec_gentau3pp4.push_back(genvis);
        	
        	//	if(TMath::Abs((*genParticles_)[p].mother(0)->pdgId())==541){
        	if(matched) isgen3matched += 1;
            //	}
        }//else{
        //          isgen3matched = false;
        //      }
      } // loop on gen particles in the event
      
      //if (tau_mu_check && tau_3prong_check) mu_3prong_event = true; // Arash
      //if (tau_mu_check) nBranches_->cutflow_perevt->Fill(5);
      //if (tau_3prong_check) nBranches_->cutflow_perevt->Fill(6);
      //if (mu_3prong_event) nBranches_->cutflow_perevt->Fill(7);
      
      //nBranches_->gen_tau_to_mu.push_back(tau_mu_check);
      //nBranches_->gen_tau_to_3prong.push_back(tau_3prong_check);
      //nBranches_->gen_tautau_to_mu3prong.push_back(mu_3prong_event);
      //cout << "    ---->    Event with mu + 3prong tau: " << mu_3prong_event << endl;
  
      //    std::cout << "\t # of gen. taus with 3prong = " << gps.size() << std::endl;
  
    }
    

  //////////////////////////////
//   std::cout << "Starts to build tau candidate out of " << numOfch << " pion candidates" << std::endl;
    
    //Efficiency table:
    float pt_bins[14] = {0, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.1};
    float eta_bins[3] = {0, 0.8, 1.5};
    float pion_rEff[3][14] = {
    {1.00, 1.01, 1.02, 1.03, 1.04, 1.07, 1.08, 1.05, 1.03, 1.02, 1.01, 1.00, 0.97, 0.96},
    {1.03, 1.03, 1.01, 1.04, 1.07, 1.06, 1.06, 1.12, 1.09, 1.10, 1.07, 1.02, 1.10, 1.13},
    {1.00, 0.99, 1.00, 0.99, 1.05, 1.05, 1.10, 1.10, 1.09, 1.11, 1.12, 1.10, 1.05, 1.04}};
    
    float pion_rEff_error[3][14] = {
    {0.03, 0.02, 0.02, 0.02, 0.03, 0.03, 0.03, 0.04, 0.03, 0.04, 0.05, 0.05, 0.06, 0.06},
    {0.05, 0.05, 0.04, 0.08, 0.05, 0.04, 0.09, 0.05, 0.07, 0.06, 0.07, 0.08, 0.10, 0.18},
    {0.00, 0.13, 0.12, 0.07, 0.08, 0.08, 0.09, 0.12, 0.11, 0.13, 0.16, 0.13, 0.18, 0.24}};
    
    
    
    std::vector<taucand> cands;

    for(int iii = 0; iii < numOfch; iii ++){
      
      reco::PFCandidate pf1 = pfcollection[iii];

      for(int jjj = iii+1; jjj < numOfch; jjj ++){
	
      	reco::PFCandidate pf2 = pfcollection[jjj];

      	for(int kkk = jjj+1; kkk < numOfch; kkk ++){

      	  reco::PFCandidate pf3 = pfcollection[kkk];
      
      	  Int_t tau_charge = pf1.charge() + pf2.charge() + pf3.charge(); 
      
      	  if(TMath::Abs(tau_charge)!=1) continue;
           
          int etaIndex[3] = {0,0,0};
          int ptIndex[3] = {0,0,0};
          for (int etabin = 0; etabin < 3; etabin++){
            if (TMath::Abs(pf1.eta()) > eta_bins[etabin]) etaIndex[0] = etabin;
            if (TMath::Abs(pf2.eta()) > eta_bins[etabin]) etaIndex[1] = etabin;
            if (TMath::Abs(pf3.eta()) > eta_bins[etabin]) etaIndex[2] = etabin;
            for (int pTbin = 0; pTbin < 14; pTbin++){
              if (pf1.pt() > pt_bins[pTbin]) ptIndex[0] = pTbin;
              if (pf2.pt() > pt_bins[pTbin]) ptIndex[1] = pTbin;
              if (pf3.pt() > pt_bins[pTbin]) ptIndex[2] = pTbin;
            }            
          }
          float rEff = 1;
          for (int p = 0; p < 3; p++) rEff *= pion_rEff[etaIndex[p]][ptIndex[p]];
          float relErr_temp[3] = {0,0,0};
          for (int p = 0; p < 3; p++){
            relErr_temp[p] = pion_rEff_error[etaIndex[p]][ptIndex[p]] / pion_rEff[etaIndex[p]][ptIndex[p]];
            relErr_temp[p] = relErr_temp[p]*relErr_temp[p];
          }
          float rEff_error = 0;
          for (int p = 0; p < 3; p++){
            rEff_error += relErr_temp[p];
          }
          rEff_error = rEff * TMath::Sqrt(rEff_error);
          

      	  std::vector<RefCountedKinematicParticle> tauParticles;
      
      	  tauParticles.push_back(pFactory.particle(mytracks[iii], aux.pion_mass, chi, ndf, aux.pion_sigma));
      	  tauParticles.push_back(pFactory.particle(mytracks[jjj], aux.pion_mass, chi, ndf, aux.pion_sigma));
      	  tauParticles.push_back(pFactory.particle(mytracks[kkk], aux.pion_mass, chi, ndf, aux.pion_sigma));
      
        
      	  //reconstructing a tau decay
      	  RefCountedKinematicTree tauTree = kpvFitter.fit(tauParticles);
      
      	  if(tauTree->isEmpty() || !tauTree->isValid() || !tauTree->isConsistent()) continue;
      
      	  //getting the J/Psi KinematicParticle
      	  tauTree->movePointerToTheTop();
      
      	  RefCountedKinematicParticle tau_part = tauTree->currentParticle();
      	  if(!tau_part->currentState().isValid()) continue;
      	  RefCountedKinematicVertex tau_vertex = tauTree->currentDecayVertex();
      	  if(!tau_vertex->vertexIsValid()) continue; 
      
      	  // 6.1.2020 commented out
      	  if(TMath::Prob(tau_vertex->chiSquared(), tau_vertex->degreesOfFreedom()) <= c_vprob) continue;
      	  
      	  std::vector< RefCountedKinematicParticle > tau_children = tauTree->finalStateParticles();
	  
      	  math::PtEtaPhiMLorentzVector tau1_fit = aux.daughter_p4(tau_children, 0);
      	  math::PtEtaPhiMLorentzVector tau2_fit = aux.daughter_p4(tau_children, 1);
      	  math::PtEtaPhiMLorentzVector tau3_fit = aux.daughter_p4(tau_children, 2);
      
      	  //	  math::PtEtaPhiMLorentzVector tlv_tau = tau1_fit + tau2_fit + tau3_fit;
      
      	  particle_cand Taucand; 
      	  Taucand = aux.calculateIPvariables(extrapolator, tau_part, tau_vertex, closestVertex);
      
      
      	  // 6.1.2020 commented out
//      	  if(Taucand.fls3d < c_fsig) continue;
      
          std::vector<RefCountedKinematicParticle> allParticles;
      
      	  allParticles.push_back(pFactory.particle(mytracks[iii], aux.pion_mass, chi, ndf, aux.pion_sigma));
      	  allParticles.push_back(pFactory.particle(mytracks[jjj], aux.pion_mass, chi, ndf, aux.pion_sigma));
      	  allParticles.push_back(pFactory.particle(mytracks[kkk], aux.pion_mass, chi, ndf, aux.pion_sigma));
      	  allParticles.push_back(pFactory.particle(tt_muon, aux.muon_mass, chi, ndf, aux.muon_sigma));
      
      	  RefCountedKinematicTree bcTree = kpvFitter.fit(allParticles);
      
      	  if(bcTree->isEmpty() || !bcTree->isValid() || !bcTree->isConsistent()) continue;
      	  
      
      	  RefCountedKinematicParticle bc_part = bcTree->currentParticle();
      	  if(!bc_part->currentState().isValid()) continue;
      
      	  RefCountedKinematicVertex bc_vertex = bcTree->currentDecayVertex();
      	  if(!bc_vertex->vertexIsValid()) continue;
       
      	  particle_cand Bcand; 
      	  Bcand = aux.calculateIPvariables(extrapolator, bc_part, bc_vertex, closestVertex);
      	  
      	  std::vector< RefCountedKinematicParticle > bc_children = bcTree->finalStateParticles();
      
      	  math::PtEtaPhiMLorentzVector tt1_fit = aux.daughter_p4(bc_children, 0);
      	  math::PtEtaPhiMLorentzVector tt2_fit = aux.daughter_p4(bc_children, 1);
      	  math::PtEtaPhiMLorentzVector tt3_fit = aux.daughter_p4(bc_children, 2);
      
      	  math::PtEtaPhiMLorentzVector tlv_tau_fit = tt1_fit + tt2_fit + tt3_fit;
      
      	  if(tlv_tau_fit.Pt() < 2.) continue;
      	  if(!(0.2 < tlv_tau_fit.M() && tlv_tau_fit.M() < 1.5)) continue;

	  // calculation of the isolation 

	  Float_t iso = 0;
	  Int_t ntracks = 0;
	  Float_t iso_mindoca = 999; 
  
	  
	  for(int itrk = 0; itrk < numOfch;  itrk++){
    
	    if(itrk==iii || itrk==jjj || itrk==kkk) continue;

	    iso += pfcollection[itrk].pt();

            TrajectoryStateOnSurface tsos_pf = extrapolator.extrapolate(mytracks[itrk].impactPointState(), bc_vertex->position());
     
    
	    VertexDistance3D a3d_pf;  

            std::pair<bool,Measurement1D> cur3DIP_pf = aux.absoluteImpactParameter(tsos_pf, bc_vertex, a3d_pf);

            Float_t pvip_pf = cur3DIP_pf.second.value();

            //    std::cout << itrk << ": Distance of closest apporach to the bc vertex = " << pvip << std::endl;
    
            if(pvip_pf < 0.03) ntracks+=1;

            if(iso_mindoca > pvip_pf) iso_mindoca = pvip_pf;
        }



	Float_t max_dr_3prong = -1;

	Float_t dR_12 = reco::deltaR(tau1_fit.Eta(), tau1_fit.Phi(), tau2_fit.Eta(), tau2_fit.Phi());
	Float_t dR_13 = reco::deltaR(tau1_fit.Eta(), tau1_fit.Phi(), tau3_fit.Eta(), tau3_fit.Phi());
	Float_t dR_23 = reco::deltaR(tau2_fit.Eta(), tau2_fit.Phi(), tau3_fit.Eta(), tau3_fit.Phi());

	if(max_dr_3prong < dR_12) max_dr_3prong = dR_12;
	if(max_dr_3prong < dR_13) max_dr_3prong = dR_13;
	if(max_dr_3prong < dR_23) max_dr_3prong = dR_23;

	Bool_t isRight = false; 
	Bool_t isRight1 = false; 
	Bool_t isRight2 = false; 
	Bool_t isRight3 = false; 
//	Float_t dr1 = 999;
//	Float_t dr2 = 999;
//	Float_t dr3 = 999;
//	Float_t ptres1 = 999;
//	Float_t ptres2 = 999;
//	Float_t ptres3 = 999;
	Int_t pid = -999;
	Float_t matched_gentaupt = -999;
	
	if(isMC){

	  for(unsigned int mmm=0; mmm < gps.size(); mmm++){
	    
	    Bool_t isRight1_ = false;
	    Bool_t isRight2_ = false;
	    Bool_t isRight3_ = false;
	    
	    std::vector<TLorentzVector> tlvs = gps[mmm];
	    
	    for(unsigned int nnn=0; nnn < tlvs.size(); nnn++){

	      if(
		 reco::deltaR(tau1_fit.Eta(), tau1_fit.Phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi()) < 0.015 &&
		 tau1_fit.Pt()/tlvs[nnn].Pt() > 0.85 && 
		 tau1_fit.Pt()/tlvs[nnn].Pt() < 1.15
		 ){

		isRight1_ = true;
    nBranches_->matched_pion_deltaPt.push_back(tau1_fit.Pt()-tlvs[nnn].Pt());
		//		dr1 = reco::deltaR(tau1_fit.Eta(), tau1_fit.Phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi());
		//		ptres1 = tau1_fit.Pt()/tlvs[nnn].Pt();

	      }	      
	      if(
		 reco::deltaR(tau2_fit.Eta(), tau2_fit.Phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi()) < 0.015 &&
		 tau2_fit.Pt()/tlvs[nnn].Pt() > 0.85 && 
		 tau2_fit.Pt()/tlvs[nnn].Pt() < 1.15
		 ){ 
		isRight2_ = true;
    nBranches_->matched_pion_deltaPt.push_back(tau2_fit.Pt()-tlvs[nnn].Pt());
		//		dr2 = reco::deltaR(tau2_fit.Eta(), tau2_fit.Phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi());
		//		ptres2 = tau2_fit.Pt()/tlvs[nnn].Pt(); 
	      }
	      
	      if(
		 reco::deltaR(tau3_fit.Eta(), tau3_fit.Phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi()) < 0.015 &&
		 tau3_fit.Pt()/tlvs[nnn].Pt() > 0.85 && 
		 tau3_fit.Pt()/tlvs[nnn].Pt() < 1.15
		 ){

		isRight3_ = true;
    nBranches_->matched_pion_deltaPt.push_back(tau3_fit.Pt()-tlvs[nnn].Pt());
		//		dr3 = reco::deltaR(tau3_fit.Eta(), tau3_fit.Phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi());
		//		ptres3 = tau3_fit.Pt()/tlvs[nnn].Pt(); 

	      }
	      
	    }
	    
	    Bool_t isRight_ = isRight1_ && isRight2_ && isRight3_;
	    if(isRight1_) isRight1 = true;
	    if(isRight2_) isRight2 = true;
	    if(isRight3_) isRight3 = true;
	    
	    if(isRight_){
	      isRight = true;
	      pid = ppdgId[mmm];
	      matched_gentaupt = vec_gentau3pp4[mmm].Pt();
	    }
	  }	
	} // isMC


	taucand _cand_ = {
	    iii,
	    jjj,
	    kkk,
	    (Float_t) tlv_tau_fit.Pt(),
	    (Float_t) tlv_tau_fit.Eta(),
	    (Float_t) tlv_tau_fit.Phi(),
	    (Float_t) tlv_tau_fit.M(),
	    //	    (Float_t) tlv_tau.Pt(),
	    //	    (Float_t) tlv_tau.Eta(),
	    //	    (Float_t) tlv_tau.Phi(),
	    //	    (Float_t) tlv_tau.M(),
//	    (Float_t) Taucand.lip, 
//	    (Float_t) Taucand.lips, 
//	    (Float_t) Taucand.pvip, 
//	    (Float_t) Taucand.pvips, 
//	    (Float_t) Taucand.fl3d,
//	    (Float_t) Taucand.fls3d, 
//	    (Float_t) Taucand.alpha,
	    Taucand,
	    (Float_t) TMath::Prob(tau_vertex->chiSquared(), tau_vertex->degreesOfFreedom()),
	    (Float_t) tau_vertex->vertexState().position().x(), 
	    (Float_t) tau_vertex->vertexState().position().y(), 
	    (Float_t) tau_vertex->vertexState().position().z(), 
	    (Float_t) max_dr_3prong, 
	    (Int_t) tau_charge,
	    (Bool_t) isRight,
	    (Bool_t) isRight1,
	    (Bool_t) isRight2,
	    (Bool_t) isRight3,
//	    (Float_t) dr1,
//	    (Float_t) dr2,
//	    (Float_t) dr3,
//	    (Float_t) ptres1,
//	    (Float_t) ptres2,
//	    (Float_t) ptres3,
	    (Int_t) pid,
	    (Float_t) matched_gentaupt, 
	    (Float_t) -1., 
	    (Float_t) -1.,
	    (Float_t) -1.,
	    (Float_t) -1.,
	    (Float_t) TMath::Prob(bc_vertex->chiSquared(), bc_vertex->degreesOfFreedom()),
	    (Float_t) bc_vertex->vertexState().position().x(),
	    (Float_t) bc_vertex->vertexState().position().y(),
	    (Float_t) bc_vertex->vertexState().position().z(),
	    (Float_t) bc_part->currentState().globalMomentum().perp(),
	    (Float_t) bc_part->currentState().globalMomentum().eta(),
	    (Float_t) bc_part->currentState().globalMomentum().phi(),
	    (Float_t) bc_part->currentState().mass(),
	    Bcand,
//	    (Float_t) Bcand.lip, 
//	    (Float_t) Bcand.lips, 
//	    (Float_t) Bcand.pvip, 
//	    (Float_t) Bcand.pvips, 
//	    (Float_t) Bcand.fl3d,
//	    (Float_t) Bcand.fls3d, 
//	    (Float_t) Bcand.alpha,
	    (Float_t) iso,
	    (Float_t) ntracks,
	    (Float_t) iso_mindoca,
	  };
     
    nBranches_->BsTauTau_tau_rEff.push_back(rEff); //relative efficiency (data/MC)
    nBranches_->BsTauTau_tau_rEff_error.push_back(rEff_error); //error on relative efficiency
	  
	  cands.push_back(_cand_);


	}
      }
    }

    sort(cands.begin(), cands.end());

    //    std::vector<Int_t> dict_idx;

//    bool isRight_bS = false;
//    bool isRight_aS = false;
//    int isRight_bS_ith = -1;
//    int isRight_aS_ith = -1;
//    int isRight_bS_n = 0;
//    int isRight_aS_n = 0;

//    std::cout << "# of taus = " << cands.size() << std::endl;

    if(cands.size()==0) return false; // commented temporarily by Arash

    Int_t ncomb = 0;

    //    std::cout << "passed!" << std::endl;

    for(int ic=0; ic < (int)cands.size(); ic++){


      //      if(cands[ic].cand_tau_pt < 2.) continue;
      ncomb += 1;
      
      //      Int_t _idx1 = cands[ic].cand_tau_id1;
      //      Int_t _idx2 = cands[ic].cand_tau_id2;
      //      Int_t _idx3 = cands[ic].cand_tau_id3;
    
//      bool flag_overlap = false;
//      for(int idc=0; idc<(int) dict_idx.size(); idc++){
//	
//	if(_idx1 == dict_idx[idc] || 
//	   _idx2 == dict_idx[idc] || 
//	   _idx3 == dict_idx[idc])
//	  
//	  flag_overlap = true;
//      }
//
//      if(cands[ic].cand_tau_isRight==true){
//	isRight_bS = true;
//	isRight_bS_ith = ic;
//	isRight_bS_n += 1;
//      }
//      
//      
//      if(flag_overlap) continue; 
//
//      if(cands[ic].cand_tau_isRight==true){
//	isRight_aS = true;
//	isRight_aS_ith = ncomb;
//	isRight_aS_n += 1;
//      }
      

      /********************************************************************
       *
       * Step9: Filling normal branches
       *
       ********************************************************************/

      // if the candidate has less than 2 GeV, remove it.
      //      if(cands[ic].cand_tau_fullfit_pt < 2.) continue;

      nBranches_->BsTauTau_tau_pt.push_back(cands[ic].cand_tau_pt);
      nBranches_->BsTauTau_tau_eta.push_back(cands[ic].cand_tau_eta);
      nBranches_->BsTauTau_tau_phi.push_back(cands[ic].cand_tau_phi);
      nBranches_->BsTauTau_tau_mass.push_back(cands[ic].cand_tau_mass);


      std::vector<Float_t> rhomass;
      reco::PFCandidate pf1 = pfcollection[cands[ic].cand_tau_id1];
      reco::PFCandidate pf2 = pfcollection[cands[ic].cand_tau_id2];
      reco::PFCandidate pf3 = pfcollection[cands[ic].cand_tau_id3];
      
      TLorentzVector tlv_pion1; 
      TLorentzVector tlv_pion2;
      TLorentzVector tlv_pion3;
      
      tlv_pion1.SetPtEtaPhiM(pf1.pt(), pf1.eta(), pf1.phi(), pf1.mass());
      tlv_pion2.SetPtEtaPhiM(pf2.pt(), pf2.eta(), pf2.phi(), pf2.mass());
      tlv_pion3.SetPtEtaPhiM(pf3.pt(), pf3.eta(), pf3.phi(), pf3.mass());

	
      if(pf1.charge()*pf2.charge() == -1){
      	TLorentzVector tlv_rho = tlv_pion1 + tlv_pion2;
      	rhomass.push_back(tlv_rho.M());
      }
      
      if(pf1.charge()*pf3.charge() == -1){
      	TLorentzVector tlv_rho = tlv_pion1 + tlv_pion3;
      	rhomass.push_back(tlv_rho.M());
      }
      
      if(pf2.charge()*pf3.charge() == -1){
      	TLorentzVector tlv_rho = tlv_pion2 + tlv_pion3;
      	rhomass.push_back(tlv_rho.M());
      }
      
	//	std::cout << "rho masses size = " << rhomass.size() << std::endl;

      nBranches_->BsTauTau_tau_rhomass1.push_back(rhomass.at(0));
      nBranches_->BsTauTau_tau_rhomass2.push_back(rhomass.at(1));

      nBranches_->BsTauTau_tau_q.push_back(cands[ic].cand_tau_charge);
      nBranches_->BsTauTau_tau_vx.push_back(cands[ic].cand_tau_vx);
      nBranches_->BsTauTau_tau_vy.push_back(cands[ic].cand_tau_vy);
      nBranches_->BsTauTau_tau_vz.push_back(cands[ic].cand_tau_vz);

      nBranches_->BsTauTau_tau_max_dr_3prong.push_back(cands[ic].cand_tau_max_dr_3prong);
      nBranches_->BsTauTau_tau_lip.push_back(cands[ic].cand_tau.lip);
      nBranches_->BsTauTau_tau_lips.push_back(cands[ic].cand_tau.lips);
      nBranches_->BsTauTau_tau_pvip.push_back(cands[ic].cand_tau.pvip);
      nBranches_->BsTauTau_tau_pvips.push_back(cands[ic].cand_tau.pvips);
      nBranches_->BsTauTau_tau_fl3d.push_back(cands[ic].cand_tau.fl3d);
      nBranches_->BsTauTau_tau_fls3d.push_back(cands[ic].cand_tau.fls3d);
      nBranches_->BsTauTau_tau_alpha.push_back(cands[ic].cand_tau.alpha);
      nBranches_->BsTauTau_tau_vprob.push_back(cands[ic].cand_tau_vprob);
      nBranches_->BsTauTau_tau_isRight.push_back(cands[ic].cand_tau_isRight);
      nBranches_->BsTauTau_tau_isRight1.push_back(cands[ic].cand_tau_isRight1);
      nBranches_->BsTauTau_tau_isRight2.push_back(cands[ic].cand_tau_isRight2);
      nBranches_->BsTauTau_tau_isRight3.push_back(cands[ic].cand_tau_isRight3);
//      nBranches_->BsTauTau_tau_dr1.push_back(cands[ic].cand_tau_dr1);
//      nBranches_->BsTauTau_tau_dr2.push_back(cands[ic].cand_tau_dr2);
//      nBranches_->BsTauTau_tau_dr3.push_back(cands[ic].cand_tau_dr3);
//      nBranches_->BsTauTau_tau_ptres1.push_back(cands[ic].cand_tau_ptres1);
//      nBranches_->BsTauTau_tau_ptres2.push_back(cands[ic].cand_tau_ptres2);
//      nBranches_->BsTauTau_tau_ptres3.push_back(cands[ic].cand_tau_ptres3);
      nBranches_->BsTauTau_tau_matched_ppdgId.push_back(cands[ic].cand_tau_matched_ppdgId);
      nBranches_->BsTauTau_tau_matched_gentaupt.push_back(cands[ic].cand_tau_matched_gentaupt);
      nBranches_->BsTauTau_tau_pfidx1.push_back(cands[ic].cand_tau_id1);
      nBranches_->BsTauTau_tau_pfidx2.push_back(cands[ic].cand_tau_id2);
      nBranches_->BsTauTau_tau_pfidx3.push_back(cands[ic].cand_tau_id3);

      nBranches_->BsTauTau_B_pt.push_back(cands[ic].cand_b_pt);
      nBranches_->BsTauTau_B_eta.push_back(cands[ic].cand_b_eta);
      nBranches_->BsTauTau_B_phi.push_back(cands[ic].cand_b_phi);
      nBranches_->BsTauTau_B_mass.push_back(cands[ic].cand_b_mass);
      nBranches_->BsTauTau_B_vprob.push_back(cands[ic].cand_b_vprob);
      nBranches_->BsTauTau_B_lip.push_back(cands[ic].cand_b.lip);
      nBranches_->BsTauTau_B_lips.push_back(cands[ic].cand_b.lips);
      nBranches_->BsTauTau_B_pvip.push_back(cands[ic].cand_b.pvip);
      nBranches_->BsTauTau_B_pvips.push_back(cands[ic].cand_b.pvips);
      nBranches_->BsTauTau_B_fls3d.push_back(cands[ic].cand_b.fls3d);
      nBranches_->BsTauTau_B_fl3d.push_back(cands[ic].cand_b.fl3d);
      nBranches_->BsTauTau_B_alpha.push_back(cands[ic].cand_b.alpha);

      nBranches_->BsTauTau_tau_pi1_pt.push_back(tlv_pion1.Pt());
      nBranches_->BsTauTau_tau_pi1_eta.push_back(tlv_pion1.Eta());
      nBranches_->BsTauTau_tau_pi1_phi.push_back(tlv_pion1.Phi());
      nBranches_->BsTauTau_tau_pi1_mass.push_back(tlv_pion1.M());
      nBranches_->BsTauTau_tau_pi1_x.push_back(pf1.vx());
      nBranches_->BsTauTau_tau_pi1_y.push_back(pf1.vy());
      nBranches_->BsTauTau_tau_pi1_z.push_back(pf1.vz());
      
      nBranches_->BsTauTau_tau_pi2_pt.push_back(tlv_pion2.Pt());
      nBranches_->BsTauTau_tau_pi2_eta.push_back(tlv_pion2.Eta());
      nBranches_->BsTauTau_tau_pi2_phi.push_back(tlv_pion2.Phi());
      nBranches_->BsTauTau_tau_pi2_mass.push_back(tlv_pion2.M());
      nBranches_->BsTauTau_tau_pi2_x.push_back(pf2.vx());
      nBranches_->BsTauTau_tau_pi2_y.push_back(pf2.vy());
      nBranches_->BsTauTau_tau_pi2_z.push_back(pf2.vz());
      
      nBranches_->BsTauTau_tau_pi3_pt.push_back(tlv_pion3.Pt());
      nBranches_->BsTauTau_tau_pi3_eta.push_back(tlv_pion3.Eta());
      nBranches_->BsTauTau_tau_pi3_phi.push_back(tlv_pion3.Phi());
      nBranches_->BsTauTau_tau_pi3_mass.push_back(tlv_pion3.M());
	    nBranches_->BsTauTau_tau_pi3_x.push_back(pf3.vx());
      nBranches_->BsTauTau_tau_pi3_y.push_back(pf3.vy());
      nBranches_->BsTauTau_tau_pi3_z.push_back(pf3.vz());
      
      std::vector<RefCountedKinematicParticle> allParticles4doc;
      
      allParticles4doc.push_back(pFactory.particle(tt_muon, aux.muon_mass, chi, ndf, aux.muon_sigma));
      allParticles4doc.push_back(pFactory.particle(mytracks[cands[ic].cand_tau_id1], aux.pion_mass, chi, ndf, aux.pion_sigma));
      allParticles4doc.push_back(pFactory.particle(mytracks[cands[ic].cand_tau_id2], aux.pion_mass, chi, ndf, aux.pion_sigma));
      allParticles4doc.push_back(pFactory.particle(mytracks[cands[ic].cand_tau_id3], aux.pion_mass, chi, ndf, aux.pion_sigma));

      nBranches_->BsTauTau_B_maxdoca.push_back(aux.getMaxDoca(allParticles4doc));
      nBranches_->BsTauTau_B_mindoca.push_back(aux.getMinDoca(allParticles4doc));

      nBranches_->BsTauTau_B_vx.push_back(cands[ic].cand_b_vx);
      nBranches_->BsTauTau_B_vy.push_back(cands[ic].cand_b_vy);
      nBranches_->BsTauTau_B_vz.push_back(cands[ic].cand_b_vz);

      nBranches_->BsTauTau_B_iso.push_back(cands[ic].cand_b_iso);
      nBranches_->BsTauTau_B_iso_ntracks.push_back(cands[ic].cand_b_iso_ntracks);
      nBranches_->BsTauTau_B_iso_mindoca.push_back(cands[ic].cand_b_iso_mindoca);
      if (print) cout << "pf size: " << (size_t)pfcollection.size() << endl;
    }

      nBranches_->BsTauTau_mu1_pt.push_back(muoncollection[0].pt());
      nBranches_->BsTauTau_mu1_eta.push_back(muoncollection[0].eta());
      nBranches_->BsTauTau_mu1_phi.push_back(muoncollection[0].phi());
      nBranches_->BsTauTau_mu1_mass.push_back(muoncollection[0].mass());
      nBranches_->BsTauTau_mu1_q.push_back(muoncollection[0].charge());
      /**************** ARASH **********************************/
      //      nBranches_->BsTauTau_mu1_isLoose.push_back(muoncollection[0].isLooseMuon());
      //      nBranches_->BsTauTau_mu1_isTight.push_back(muoncollection[0].isTightMuon(closestVertex));
      //nBranches_->BsTauTau_mu1_isSoft.push_back(muoncollection[0].isSoftMuon(closestVertex));
      bool softMu = false;
      auto mu = muoncollection[0];
      bool muID = muon::isGoodMuon(mu, muon::TMOneStationTight);
      if (muID){
        bool layers = mu.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 && mu.innerTrack()->hitPattern().pixelLayersWithMeasurement() > 0;
        bool ishighq = mu.innerTrack()->quality(reco::Track::highPurity);
        bool ip = std::abs(mu.innerTrack()->dxy(closestVertex.position())) < 0.3 && std::abs(mu.innerTrack()->dz(closestVertex.position())) < 20.;
        softMu =  layers && ip && ishighq;
      }
      nBranches_->BsTauTau_mu1_isSoft.push_back(softMu);
      /**************** ARASH **********************************/
      nBranches_->BsTauTau_mu1_isPF.push_back(muoncollection[0].isPFMuon());
      nBranches_->BsTauTau_mu1_isGlobal.push_back(muoncollection[0].isGlobalMuon());
      nBranches_->BsTauTau_mu1_isTracker.push_back(muoncollection[0].isTrackerMuon());
      nBranches_->BsTauTau_mu1_vx.push_back(muoncollection[0].vx());
      nBranches_->BsTauTau_mu1_vy.push_back(muoncollection[0].vy());
      nBranches_->BsTauTau_mu1_vz.push_back(muoncollection[0].vz());
      nBranches_->BsTauTau_mu1_iso.push_back(1.);
      nBranches_->BsTauTau_mu1_dbiso.push_back(aux.MuonPFIso(muoncollection[0]));
  
      nBranches_->BsTauTau_PV_vx.push_back(vertices_->begin()->position().x());
      nBranches_->BsTauTau_PV_vy.push_back(vertices_->begin()->position().y());
      nBranches_->BsTauTau_PV_vz.push_back(vertices_->begin()->position().z());

//      if(myVertex.isValid()){
//	nBranches_->BsTauTau_bbPV_refit_vx.push_back(myVertex.position().x());
//	nBranches_->BsTauTau_bbPV_refit_vy.push_back(myVertex.position().y());
//	nBranches_->BsTauTau_bbPV_refit_vz.push_back(myVertex.position().z());
//      }else{
      nBranches_->BsTauTau_bbPV_refit_vx.push_back(-1);
      nBranches_->BsTauTau_bbPV_refit_vy.push_back(-1);
      nBranches_->BsTauTau_bbPV_refit_vz.push_back(-1);
      //      }

      nBranches_->BsTauTau_bbPV_vx.push_back(closestVertex.position().x());
      nBranches_->BsTauTau_bbPV_vy.push_back(closestVertex.position().y());
      nBranches_->BsTauTau_bbPV_vz.push_back(closestVertex.position().z());
      
      


      //////////////////////////////




    /********************************************************************
     *
     * Step10: check gen-matching and fill them
     *
     ********************************************************************/

    TVector3 genvertex(-9.,-9.,-9.);
  
    std::vector<const reco::Candidate*> gen_nr_mu;
  
    //      std::cout << "passed 2!" << std::endl;

    if(isMC){

    	//	std::cout << "passed 3: " << genParticles_->size() << std::endl;
        
    	for( unsigned p=0; p < genParticles_->size(); ++p){
    
    	  //	  std::cout << "gen: " << (*genParticles_)[p].pdgId() << " " << (*genParticles_)[p].status() << std::endl;
    
    	  // Bc daughters loop
    	  if(TMath::Abs((*genParticles_)[p].pdgId())==531 && (*genParticles_)[p].status()==2){
    
    	    // retrieve production vertex
    	    genvertex = aux.getVertex((*genParticles_)[p]);
    
    	    for(int idd = 0; idd < (int)(*genParticles_)[p].numberOfDaughters(); idd++){
    	      Int_t dpid = (*genParticles_)[p].daughter(idd)->pdgId();
    	      //	      std::cout << "\t -> " << (*genParticles_)[p].daughter(idd)->status()<< std::endl;
    	      if(TMath::Abs(dpid)==15) gen_nr_mu.push_back((*genParticles_)[p].daughter(idd));
    	    }
    	  }
    
          
    	}
    }
  
    // -9 if there is no Bc found 
    nBranches_->BsTauTau_genPV_vx.push_back(genvertex.x());
    nBranches_->BsTauTau_genPV_vy.push_back(genvertex.y());
    nBranches_->BsTauTau_genPV_vz.push_back(genvertex.z());
    nBranches_->BsTauTau_ngenmuons.push_back(gen_nr_mu.size());

    nBranches_->BsTauTau_isgen3.push_back(isgen3);
    nBranches_->BsTauTau_isgen3matched.push_back(isgen3matched);
    nBranches_->BsTauTau_nch.push_back(numOfch);
    nBranches_->BsTauTau_nch_qr.push_back(npf_qr);
    nBranches_->BsTauTau_ngentau3.push_back(gps.size());
    nBranches_->BsTauTau_ngentau.push_back(vec_gentaudm.size());

    if(vec_gentaudm.size() >=1){
      nBranches_->BsTauTau_gentaupt.push_back(vec_gentaup4[0].Pt());
      nBranches_->BsTauTau_gentaudm.push_back(vec_gentaudm[0]);
    }else{
      nBranches_->BsTauTau_gentaupt.push_back(-1);
      nBranches_->BsTauTau_gentaudm.push_back(-1);
    }

    nBranches_->IsBsTauTau.push_back(1.);
    nBranches_->BsTauTau_nCandidates.push_back(ncomb);
    
    nBranches_->cutflow_perevt->Fill(3);
    //cout << "number of taus: " << (size_t)nBranches_->gen_tau_pt.size() << endl;
    return true;


}

BsTauTauNtuplizer::ECaloType

BsTauTauNtuplizer::GetTowerSubdetHad(double&eta) const
{
//  if (TMath::Abs(eta)>3. && TMath::Abs(eta)<5.2)
//    std::cout << " - - - - - - - - - - - - - >> " << eta << " " << -maxEtaHF << " " << -minEtaHF << " " << kHFm << " " << (eta > -maxEtaHF) << " " << (eta < -minEtaHF) << std::endl;
  if(eta > -maxEtaHF && eta < -minEtaHF) return kHFm;

//  if (TMath::Abs(eta)>3. && TMath::Abs(eta)<5.2)
//    std::cout << " - - - - - - - - - - - - - - - >> " << eta << " " << minEtaHF << " " << maxEtaHF << " " << kHFp << " " << (eta >  minEtaHF) << " " << (eta <  maxEtaHF) << std::endl;
  if(eta >  minEtaHF && eta <  maxEtaHF) return kHFp;

  if(fabs(eta) > 0 && fabs(eta) < maxEtaHB) return kHB;

  if(fabs(eta) > minEtaHE && fabs(eta) < maxEtaHE) return kHE;

  return nCaloTypes;
}


//void BsTauTauNtuplizer::printout(const RefCountedKinematicVertex& myVertex){
//    std::cout << "Vertex:" << std::endl;
//    if (myVertex->vertexIsValid()) {
//        std::cout << "\t Decay vertex: " << myVertex->position() << myVertex->chiSquared() << " " << myVertex->degreesOfFreedom()
//                  << std::endl;
//    } else
//        std::cout << "\t Decay vertex Not valid\n";
//}
//
//void BsTauTauNtuplizer::printout(const RefCountedKinematicParticle& myParticle){
//    std::cout << "Particle:" << std::endl;
//    //accessing the reconstructed Bs meson parameters:
//    //SK: uncomment if needed  AlgebraicVector7 bs_par = myParticle->currentState().kinematicParameters().vector();
//
//    //and their joint covariance matrix:
//    //SK:uncomment if needed  AlgebraicSymMatrix77 bs_er = myParticle->currentState().kinematicParametersError().matrix();
//    std::cout << "\t Momentum at vertex: " << myParticle->currentState().globalMomentum() << std::endl;
//    std::cout << "\t Parameters at vertex: " << myParticle->currentState().kinematicParameters().vector() << std::endl;
//}
//
//void BsTauTauNtuplizer::printout(const RefCountedKinematicTree& myTree){
//    if (!myTree->isValid()) {
//        std::cout << "Tree is invalid. Fit failed.\n";
//        return;
//    }
//
//    //accessing the tree components, move pointer to top
//    myTree->movePointerToTheTop();
//
//    //We are now at the top of the decay tree getting the B_s reconstructed KinematicPartlcle
//    RefCountedKinematicParticle b_s = myTree->currentParticle();
//    printout(b_s);
//
//    // The B_s decay vertex
//    RefCountedKinematicVertex b_dec_vertex = myTree->currentDecayVertex();
//    printout(b_dec_vertex);
//
//    // Get all the children of Bs:
//    //In this way, the pointer is not moved
//    std::vector<RefCountedKinematicParticle> bs_children = myTree->finalStateParticles();
//
//    for (unsigned int i = 0; i < bs_children.size(); ++i) {
//        printout(bs_children[i]);
//    }
//
//    std::cout << "\t ------------------------------------------" << std::endl;
//
//    //Now navigating down the tree , pointer is moved:
//    bool child = myTree->movePointerToTheFirstChild();
//
//    if (child)
//        while (myTree->movePointerToTheNextChild()) {
//            RefCountedKinematicParticle aChild = myTree->currentParticle();
//            printout(aChild);
//        }
//}
