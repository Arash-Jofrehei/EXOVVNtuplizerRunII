#include "../interface/NtupleBranches.h"

//===================================================================================================================        
NtupleBranches::NtupleBranches( std::map< std::string, bool >& runFlags, TTree* tree )
  : tree_( tree )
{
  branch( runFlags );
}

//===================================================================================================================
NtupleBranches::~NtupleBranches( void )
{
}

//===================================================================================================================      
void NtupleBranches::branch( std::map< std::string, bool >& runFlags ){


  if ( runFlags["runOnMC"] ){
    if ( runFlags["doGenParticles"] ){
      /** genParticles */
      tree_->Branch( "genParticle_N"	     , &genParticle_N	       );
      tree_->Branch( "genParticle_pt"	     , &genParticle_pt	       ); 
      //      tree_->Branch( "genParticle_px"	     , &genParticle_px	       ); 
      //      tree_->Branch( "genParticle_py"	     , &genParticle_py	       ); 
      //      tree_->Branch( "genParticle_pz"	     , &genParticle_pz	       ); 
      //      tree_->Branch( "genParticle_e" 	     , &genParticle_e	       ); 
      tree_->Branch( "genParticle_eta"	     , &genParticle_eta        ); 
      tree_->Branch( "genParticle_phi"	     , &genParticle_phi        ); 
      tree_->Branch( "genParticle_mass"	     , &genParticle_mass       ); 
      tree_->Branch( "genParticle_pdgId"     , &genParticle_pdgId      );
      tree_->Branch( "genParticle_status"    , &genParticle_status     );
      tree_->Branch( "genParticle_isPrompt"  , &genParticle_isPrompt   );
      tree_->Branch( "genParticle_isDirectPromptTauDecayProduct"  , &genParticle_isDirectPromptTauDecayProduct);
      tree_->Branch( "genParticle_isDirectHardProcessTauDecayProductFinalState"  , &genParticle_isDirectHardProcessTauDecayProductFinalState);
      tree_->Branch( "genParticle_fromHardProcessFinalState"  , &genParticle_fromHardProcessFinalState   );
      tree_->Branch( "genParticle_mother"    , &genParticle_mother     );
      tree_->Branch( "genParticle_mother_pt" , &genParticle_mother_pt  );
      tree_->Branch( "genParticle_nMoth"     , &genParticle_nMoth      );
      tree_->Branch( "genParticle_nDau"	     , &genParticle_nDau       ); 
      tree_->Branch( "genParticle_dau"	     , &genParticle_dau        );
	
	
	
    } //doGenParticles
      
    if ( runFlags["doGenEvent"] ){
      /** generator info */
      tree_->Branch( "lheV_pt"	             , &lheV_pt                ); 
      tree_->Branch( "lheHT"	             , &lheHT                  ); 
      tree_->Branch( "lheNj"	             , &lheNj                  );
      tree_->Branch( "lheNb"	             , &lheNb                  );
      tree_->Branch( "lheNl"	             , &lheNl                  );
      tree_->Branch( "lheV_mass"           , &lheV_mass              ); 
      tree_->Branch( "genWeight"	         , &genWeight              );
      tree_->Branch( "genFacWeightUp"	     , &genFacWeightUp         );
      tree_->Branch( "genFacWeightDown"	   , &genFacWeightDown       );
      tree_->Branch( "genRenWeightUp"	     , &genRenWeightUp         );
      tree_->Branch( "genRenWeightDown"	   , &genRenWeightDown       );
      tree_->Branch( "genFacRenWeightUp"	 , &genFacRenWeightUp      );
      tree_->Branch( "genFacRenWeightDown" , &genFacRenWeightDown    );
      tree_->Branch( "qScale"	             , &qScale                 );
      tree_->Branch( "PDF_rms"	           , &PDF_rms                );
      tree_->Branch( "PDF_x"	             , &PDF_x                  );
      tree_->Branch( "PDF_xPDF"	           , &PDF_xPDF               );
      tree_->Branch( "PDF_id"	             , &PDF_id                 );

    } //doGenEvent
  } //runOnMC
  
  
  if (runFlags["doTriggerDecisions"]) {
    /** HLT trigger decisions */
    tree_->Branch("HLT_isFired", &HLT_isFired );
  }

  
  if (runFlags["doTriggerObjects"]) {
    /** HLT trigger objects */
    tree_->Branch("triggerObject_pt"		, &triggerObject_pt		);
    tree_->Branch("triggerObject_eta"		, &triggerObject_eta		);
    tree_->Branch("triggerObject_phi"		, &triggerObject_phi	        );
    tree_->Branch("triggerObject_mass"		, &triggerObject_mass	        );
    tree_->Branch("triggerObject_lastname"	, &triggerObject_lastname	);
    tree_->Branch("triggerObject_filterLabels"	, &triggerObject_filterLabels	);
    tree_->Branch("triggerObject_firedTrigger"	, &triggerObject_firedTrigger	);
    tree_->Branch("triggerObject_filterIDs"	, &triggerObject_filterIDs	);

  } //doTriggerObjects
  
  if (runFlags["doHltFilters"]) {
    /** HLT filter decisions */
    tree_->Branch("passFilter_HBHE"                 ,&passFilter_HBHE_                ,"passFilter_HBHE_/O");
    tree_->Branch("passFilter_HBHELoose"            ,&passFilter_HBHELoose_	          ,"passFilter_HBHELoose_/O");
    tree_->Branch("passFilter_HBHETight"            ,&passFilter_HBHETight_	          ,"passFilter_HBHETight_/O");
    tree_->Branch("passFilter_HBHEIso"              ,&passFilter_HBHEIso_	            ,"passFilter_HBHEIso_/O");
    tree_->Branch("passFilter_CSCHalo"              ,&passFilter_CSCHalo_             ,"passFilter_CSCHalo_/O");
    tree_->Branch("passFilter_CSCTightHalo2015"     ,&passFilter_CSCTightHalo2015_    ,"passFilter_CSCTightHalo2015_/O");
    tree_->Branch("passFilter_HCALlaser"            ,&passFilter_HCALlaser_           ,"passFilter_HCALlaser_/O");
    tree_->Branch("passFilter_ECALDeadCell"         ,&passFilter_ECALDeadCell_        ,"passFilter_ECALDeadCell_/O");
    tree_->Branch("passFilter_GoodVtx"              ,&passFilter_GoodVtx_             ,"passFilter_GoodVtx_/O");
    tree_->Branch("passFilter_TrkFailure"           ,&passFilter_TrkFailure_          ,"passFilter_TrkFailure_/O");
    tree_->Branch("passFilter_EEBadSc"              ,&passFilter_EEBadSc_             ,"passFilter_EEBadSc_/O");
    tree_->Branch("passFilter_ECALlaser"            ,&passFilter_ECALlaser_           ,"passFilter_ECALlaser_/O");
    tree_->Branch("passFilter_TrkPOG"               ,&passFilter_TrkPOG_              ,"passFilter_TrkPOG_/O");
    tree_->Branch("passFilter_TrkPOG_manystrip"     ,&passFilter_TrkPOG_manystrip_    ,"passFilter_TrkPOG_manystrip_/O");
    tree_->Branch("passFilter_TrkPOG_toomanystrip"  ,&passFilter_TrkPOG_toomanystrip_ ,"passFilter_TrkPOG_toomanystrip_/O");
    tree_->Branch("passFilter_TrkPOG_logError"      ,&passFilter_TrkPOG_logError_     ,"passFilter_TrkPOG_logError_/O");
    tree_->Branch("passFilter_METFilters"           ,&passFilter_METFilters_          ,"passFilter_METFilters_/O");
    
    //NEW FOR ICHEP
    tree_->Branch("passFilter_CSCTightHaloTrkMuUnvetoFilter", &passFilter_CSCTightHaloTrkMuUnvetoFilter_   ,"passFilter_CSCTightHaloTrkMuUnvetoFilter_/O");
    tree_->Branch("passFilter_globalTightHalo2016"          , &passFilter_globalTightHalo2016_             ,"passFilter_globalTightHalo2016_/O");
    tree_->Branch("passFilter_globalSuperTightHalo2016"          , &passFilter_globalSuperTightHalo2016_             ,"passFilter_globalSuperTightHalo2016_/O");
    tree_->Branch("passFilter_HcalStripHalo"                , &passFilter_HcalStripHalo_                   ,"passFilter_HcalStripHalo_/O");
    tree_->Branch("passFilter_chargedHadronTrackResolution" , &passFilter_chargedHadronTrackResolution_    ,"passFilter_chargedHadronTrackResolution_/O");
    tree_->Branch("passFilter_muonBadTrack"                 , &passFilter_muonBadTrack_                    ,"passFilter_muonBadTrack_/O");
    tree_->Branch("flag_badMuons"                 , &flag_badMuons_                    ,"flag_badMuons_/O");
    tree_->Branch("flag_duplicateMuons"                 , &flag_duplicateMuons_                    ,"flag_duplicateMuons_/O");
    tree_->Branch("flag_nobadMuons"                 , &flag_nobadMuons_                    ,"flag_nobadMuons_/O");
    tree_->Branch("passFilter_ecalBadCalib_"    ,&passFilter_ecalBadCalib_, "passFilter_ecalBadCalib_/O");
  } //do HltFilters

  if (runFlags["doMissingEt"]) {
    /** MET */
    tree_->Branch( "rho", &rho );
    tree_->Branch("METraw_et"		        , &METraw_et	     );
    tree_->Branch("METraw_phi"		        , &METraw_phi	     ); 
    tree_->Branch("METraw_sumEt"		, &METraw_sumEt	     );   
    tree_->Branch("MET_corrPx"		        , &MET_corrPx	     ); 
    tree_->Branch("MET_corrPy"		        , &MET_corrPy	     );   
    tree_->Branch("MET_et"	                , &MET_et  	     ); 
    tree_->Branch("MET_phi"	                , &MET_phi           );
    tree_->Branch("MET_puppi_et"	        , &MET_puppi_et      ); 
    tree_->Branch("MET_puppi_phi"               , &MET_puppi_phi     );
    tree_->Branch("MET_sumEt"	                , &MET_sumEt 	     ); 
    tree_->Branch("MET_JetEnUp"	                , &MET_JetEnUp 	     ); 
    tree_->Branch("MET_JetEnDown"	                , &MET_JetEnDown 	     ); 
    tree_->Branch("MET_JetResUp"	                , &MET_JetResUp 	     ); 
    tree_->Branch("MET_JetResDown"	                , &MET_JetResDown 	     ); 
    tree_->Branch("MET_UnclusteredEnUp"	                , &MET_UnclusteredEnUp 	     ); 
    tree_->Branch("MET_UnclusteredEnDown"	                , &MET_UnclusteredEnDown 	     ); 
    
  } //doMissingEt

  if ( runFlags["doMVAMET"] ){
    /** MET SVift*/
    tree_->Branch("MET_Nmva"	                , &MET_Nmva 	     ); 
    tree_->Branch("MET_mva_et"	                , &MET_mva_et        ); 
    tree_->Branch("MET_mva_phi"                 , &MET_mva_phi       );
    tree_->Branch( "MET_mva_cov00"                                        , &MET_mva_cov00 );
    tree_->Branch( "MET_mva_cov10"                                        , &MET_mva_cov10 );
    tree_->Branch( "MET_mva_cov11"                                        , &MET_mva_cov11 );
    tree_->Branch( "MET_mva_recoil_pt"                                        , &MET_mva_recoil_pt );
    tree_->Branch( "MET_mva_recoil_eta"                                        , &MET_mva_recoil_eta );
    tree_->Branch( "MET_mva_recoil_phi"                                        , &MET_mva_recoil_phi );
    tree_->Branch( "MET_mva_recoil_pdgId"                                        , &MET_mva_recoil_pdgId );

  }

  
  /*------------- ------EVENT infos-----------------------------*/
  tree_->Branch("EVENT_event"	 , &EVENT_event     );
  tree_->Branch("EVENT_run"	 , &EVENT_run	    );
  tree_->Branch("EVENT_lumiBlock", &EVENT_lumiBlock );
  
  if (runFlags["runOnMC"]) {
    if (runFlags["doPileUp"]) {
      /*--------------------------PU infos--------------------------*/
      tree_->Branch("nPuVtxTrue", &nPuVtxTrue );	
      tree_->Branch("nPuVtx"    , &nPuVtx     );
      tree_->Branch("bX"	, &bX	      );
    } //doPileUp
  } //runOnMC
  
  if (runFlags["doVertices"]) {  
    /*--------------------------PV infos--------------------------*/
    tree_->Branch("PV_N"     , &PV_N      );
    tree_->Branch("PV_filter", &PV_filter );
    tree_->Branch("PV_chi2"  , &PV_chi2   );
    tree_->Branch("PV_ndof"  , &PV_ndof   );
    tree_->Branch("PV_rho"   , &PV_rho    );
    tree_->Branch("PV_z"     , &PV_z      );
    tree_->Branch("BeamSpot_x0" , &BeamSpot_x0     );
    tree_->Branch("BeamSpot_y0" , &BeamSpot_y0     );
    tree_->Branch("BeamSpot_z0" , &BeamSpot_z0     );


  }



  if (runFlags["doBsTauTau"] || runFlags["dotwoPiBsTauTau"]){
    tree_->Branch("IsBsTauTau", &IsBsTauTau );

    tree_->Branch("BsTauTau_calo_eta",       &BsTauTau_calo_eta       );
    tree_->Branch("BsTauTau_calo_phi",       &BsTauTau_calo_phi       );
    tree_->Branch("BsTauTau_calo_energy",    &BsTauTau_calo_energy    );
    tree_->Branch("BsTauTau_calo_energyHFp", &BsTauTau_calo_energyHFp );
    tree_->Branch("BsTauTau_calo_energyHFm", &BsTauTau_calo_energyHFm );
    tree_->Branch("BsTauTau_calo_zdcSumPlus", &BsTauTau_calo_zdcSumPlus );
    tree_->Branch("BsTauTau_calo_zdcSumMinus", &BsTauTau_calo_zdcSumMinus );

    tree_->Branch("BsTauTau_trackPFactivity_pt",  &BsTauTau_trackPFactivity_pt);
    tree_->Branch("BsTauTau_trackPFactivity_eta", &BsTauTau_trackPFactivity_eta);
    tree_->Branch("BsTauTau_trackPFactivity_phi", &BsTauTau_trackPFactivity_phi);

    tree_->Branch("BsTauTau_nCandidates", &BsTauTau_nCandidates );

    tree_->Branch("BsTauTau_mu1_pt", &BsTauTau_mu1_pt );
    tree_->Branch("BsTauTau_mu1_eta", &BsTauTau_mu1_eta );
    tree_->Branch("BsTauTau_mu1_phi", &BsTauTau_mu1_phi );
    tree_->Branch("BsTauTau_mu1_mass", &BsTauTau_mu1_mass );
    tree_->Branch("BsTauTau_mu1_q", &BsTauTau_mu1_q );
    tree_->Branch("BsTauTau_mu1_isLoose"  , &BsTauTau_mu1_isLoose   );
    tree_->Branch("BsTauTau_mu1_isTight"  , &BsTauTau_mu1_isTight   );
    tree_->Branch("BsTauTau_mu1_isPF"     , &BsTauTau_mu1_isPF      );
    tree_->Branch("BsTauTau_mu1_isGlobal" , &BsTauTau_mu1_isGlobal  );
    tree_->Branch("BsTauTau_mu1_isTracker", &BsTauTau_mu1_isTracker );
    tree_->Branch("BsTauTau_mu1_isSoft"   , &BsTauTau_mu1_isSoft    );
    tree_->Branch("BsTauTau_mu1_vx"   , &BsTauTau_mu1_vx    );
    tree_->Branch("BsTauTau_mu1_vy"   , &BsTauTau_mu1_vy    );
    tree_->Branch("BsTauTau_mu1_vz"   , &BsTauTau_mu1_vz    );
    tree_->Branch("BsTauTau_mu1_iso"   , &BsTauTau_mu1_iso    );
    tree_->Branch("BsTauTau_mu1_dbiso"   , &BsTauTau_mu1_dbiso    );


    tree_->Branch("BsTauTau_tau_pt", &BsTauTau_tau_pt );
    tree_->Branch("BsTauTau_tau_eta", &BsTauTau_tau_eta );
    tree_->Branch("BsTauTau_tau_phi", &BsTauTau_tau_phi );
    tree_->Branch("BsTauTau_tau_mass", &BsTauTau_tau_mass );
    tree_->Branch("BsTauTau_tau_rhomass1", &BsTauTau_tau_rhomass1 );
    tree_->Branch("BsTauTau_tau_rhomass2", &BsTauTau_tau_rhomass2 );
    tree_->Branch("BsTauTau_tau_q", &BsTauTau_tau_q );
    tree_->Branch("BsTauTau_tau_vx"   , &BsTauTau_tau_vx    );
    tree_->Branch("BsTauTau_tau_vy"   , &BsTauTau_tau_vy    );
    tree_->Branch("BsTauTau_tau_vz"   , &BsTauTau_tau_vz    );

    tree_->Branch("BsTauTau_tau_max_dr_3prong", &BsTauTau_tau_max_dr_3prong);
    tree_->Branch("BsTauTau_tau_lip", &BsTauTau_tau_lip);
    tree_->Branch("BsTauTau_tau_lips", &BsTauTau_tau_lips);
    tree_->Branch("BsTauTau_tau_pvip", &BsTauTau_tau_pvip);
    tree_->Branch("BsTauTau_tau_pvips", &BsTauTau_tau_pvips);
    tree_->Branch("BsTauTau_tau_fl3d", &BsTauTau_tau_fl3d);
    tree_->Branch("BsTauTau_tau_fls3d", &BsTauTau_tau_fls3d);
    tree_->Branch("BsTauTau_tau_alpha", &BsTauTau_tau_alpha);
    tree_->Branch("BsTauTau_tau_vprob", &BsTauTau_tau_vprob);
    tree_->Branch("BsTauTau_tau_isRight", &BsTauTau_tau_isRight);
    tree_->Branch("BsTauTau_tau_isRight1", &BsTauTau_tau_isRight1);
    tree_->Branch("BsTauTau_tau_isRight2", &BsTauTau_tau_isRight2);
    tree_->Branch("BsTauTau_tau_isRight3", &BsTauTau_tau_isRight3);
//    tree_->Branch("BsTauTau_tau_dr1", &BsTauTau_tau_dr1);
//    tree_->Branch("BsTauTau_tau_dr2", &BsTauTau_tau_dr2);
//    tree_->Branch("BsTauTau_tau_dr3", &BsTauTau_tau_dr3);
//    tree_->Branch("BsTauTau_tau_ptres1", &BsTauTau_tau_ptres1);
//    tree_->Branch("BsTauTau_tau_ptres2", &BsTauTau_tau_ptres2);
//    tree_->Branch("BsTauTau_tau_ptres3", &BsTauTau_tau_ptres3);
    tree_->Branch("BsTauTau_tau_matched_ppdgId", &BsTauTau_tau_matched_ppdgId);
    tree_->Branch("BsTauTau_tau_matched_gentaupt", &BsTauTau_tau_matched_gentaupt);
    //    tree_->Branch("BsTauTau_tau_gentaupt", &BsTauTau_tau_gentaupt);
    tree_->Branch("BsTauTau_tau_pfidx1", &BsTauTau_tau_pfidx1);
    tree_->Branch("BsTauTau_tau_pfidx2", &BsTauTau_tau_pfidx2);
    tree_->Branch("BsTauTau_tau_pfidx3", &BsTauTau_tau_pfidx3);

    tree_->Branch("reco_pion_pt", &reco_pion_pt );
    tree_->Branch("reco_pion_eta", &reco_pion_eta );
    tree_->Branch("reco_pion_phi", &reco_pion_phi );

    tree_->Branch("BsTauTau_tau_pi1_pt", &BsTauTau_tau_pi1_pt );
    tree_->Branch("BsTauTau_tau_pi1_eta", &BsTauTau_tau_pi1_eta );
    tree_->Branch("BsTauTau_tau_pi1_phi", &BsTauTau_tau_pi1_phi );
    tree_->Branch("BsTauTau_tau_pi1_mass", &BsTauTau_tau_pi1_mass );
    tree_->Branch("BsTauTau_tau_pi1_x", &BsTauTau_tau_pi1_x );
    tree_->Branch("BsTauTau_tau_pi1_y", &BsTauTau_tau_pi1_y );
    tree_->Branch("BsTauTau_tau_pi1_z", &BsTauTau_tau_pi1_z );

    tree_->Branch("BsTauTau_tau_pi2_pt", &BsTauTau_tau_pi2_pt );
    tree_->Branch("BsTauTau_tau_pi2_eta", &BsTauTau_tau_pi2_eta );
    tree_->Branch("BsTauTau_tau_pi2_phi", &BsTauTau_tau_pi2_phi );
    tree_->Branch("BsTauTau_tau_pi2_mass", &BsTauTau_tau_pi2_mass );
    tree_->Branch("BsTauTau_tau_pi2_x", &BsTauTau_tau_pi2_x );
    tree_->Branch("BsTauTau_tau_pi2_y", &BsTauTau_tau_pi2_y );
    tree_->Branch("BsTauTau_tau_pi2_z", &BsTauTau_tau_pi2_z );

    tree_->Branch("BsTauTau_tau_pi3_pt", &BsTauTau_tau_pi3_pt );
    tree_->Branch("BsTauTau_tau_pi3_eta", &BsTauTau_tau_pi3_eta );
    tree_->Branch("BsTauTau_tau_pi3_phi", &BsTauTau_tau_pi3_phi );
    tree_->Branch("BsTauTau_tau_pi3_mass", &BsTauTau_tau_pi3_mass );
    tree_->Branch("BsTauTau_tau_pi3_x", &BsTauTau_tau_pi3_x );
    tree_->Branch("BsTauTau_tau_pi3_y", &BsTauTau_tau_pi3_y );
    tree_->Branch("BsTauTau_tau_pi3_z", &BsTauTau_tau_pi3_z );

    tree_->Branch("BsTauTau_PV_vx", &BsTauTau_PV_vx );
    tree_->Branch("BsTauTau_PV_vy", &BsTauTau_PV_vy );
    tree_->Branch("BsTauTau_PV_vz", &BsTauTau_PV_vz );

    tree_->Branch("BsTauTau_bbPV_vx", &BsTauTau_bbPV_vx );
    tree_->Branch("BsTauTau_bbPV_vy", &BsTauTau_bbPV_vy );
    tree_->Branch("BsTauTau_bbPV_vz", &BsTauTau_bbPV_vz );

    tree_->Branch("BsTauTau_bbPV_refit_vx", &BsTauTau_bbPV_vx );
    tree_->Branch("BsTauTau_bbPV_refit_vy", &BsTauTau_bbPV_vy );
    tree_->Branch("BsTauTau_bbPV_refit_vz", &BsTauTau_bbPV_vz );

    tree_->Branch("BsTauTau_genPV_vx", &BsTauTau_genPV_vx );
    tree_->Branch("BsTauTau_genPV_vy", &BsTauTau_genPV_vy );
    tree_->Branch("BsTauTau_genPV_vz", &BsTauTau_genPV_vz );

    tree_->Branch("BsTauTau_B_pt", &BsTauTau_B_pt );
    tree_->Branch("BsTauTau_B_eta", &BsTauTau_B_eta );
    tree_->Branch("BsTauTau_B_phi", &BsTauTau_B_phi );
    tree_->Branch("BsTauTau_B_mass", &BsTauTau_B_mass );
    tree_->Branch("BsTauTau_B_vprob", &BsTauTau_B_vprob );
    tree_->Branch("BsTauTau_B_lip", &BsTauTau_B_lip);
    tree_->Branch("BsTauTau_B_lips", &BsTauTau_B_lips);
    tree_->Branch("BsTauTau_B_pvip", &BsTauTau_B_pvip);
    tree_->Branch("BsTauTau_B_pvips", &BsTauTau_B_pvips);
    tree_->Branch("BsTauTau_B_fl3d", &BsTauTau_B_fl3d);
    tree_->Branch("BsTauTau_B_fls3d", &BsTauTau_B_fls3d);
    tree_->Branch("BsTauTau_B_alpha", &BsTauTau_B_alpha);
    tree_->Branch("BsTauTau_B_maxdoca", &BsTauTau_B_maxdoca);
    tree_->Branch("BsTauTau_B_mindoca", &BsTauTau_B_mindoca);
    tree_->Branch("BsTauTau_B_vx", &BsTauTau_B_vx );
    tree_->Branch("BsTauTau_B_vy", &BsTauTau_B_vy );
    tree_->Branch("BsTauTau_B_vz", &BsTauTau_B_vz );
    tree_->Branch("BsTauTau_B_iso", &BsTauTau_B_iso);
    tree_->Branch("BsTauTau_B_iso_ntracks", &BsTauTau_B_iso_ntracks );
    tree_->Branch("BsTauTau_B_iso_mindoca", &BsTauTau_B_iso_mindoca );

    tree_->Branch("BsTauTau_ngenmuons", &BsTauTau_ngenmuons);
    tree_->Branch("BsTauTau_isgen3", &BsTauTau_isgen3);
    tree_->Branch("BsTauTau_isgen3matched", &BsTauTau_isgen3matched);
    tree_->Branch("BsTauTau_tau_rEff", &BsTauTau_tau_rEff);
    tree_->Branch("BsTauTau_tau_rEff_error", &BsTauTau_tau_rEff_error);
    tree_->Branch("BsTauTau_nPions", &BsTauTau_nPions);
    tree_->Branch("BsTauTau_nch", &BsTauTau_nch);
    tree_->Branch("BsTauTau_nch_qr", &BsTauTau_nch_qr);
    tree_->Branch("BsTauTau_ngentau3", &BsTauTau_ngentau3); 
    tree_->Branch("BsTauTau_ngentau", &BsTauTau_ngentau);
    tree_->Branch("BsTauTau_gentaupt", &BsTauTau_gentaupt);
    tree_->Branch("BsTauTau_gentaudm", &BsTauTau_gentaudm);
    
    tree_->Branch("PVz",&PVz);
    tree_->Branch("pion_z",&pion_z);
    tree_->Branch("matched_pion_deltaPt",&matched_pion_deltaPt);
    
    tree_->Branch("triggered",&triggered);
    tree_->Branch("gen_tau_pt",&gen_tau_pt);
    tree_->Branch("gen_tau_eta",&gen_tau_eta);
    tree_->Branch("gen_tau_phi",&gen_tau_phi);
    tree_->Branch("gen_tau_daughter_pdgId",&gen_tau_daughter_pdgId);
    tree_->Branch("gen_tau_daughter_pt",&gen_tau_daughter_pt);
    tree_->Branch("gen_tau_daughter_eta",&gen_tau_daughter_eta);
    tree_->Branch("gen_tau_daughter_phi",&gen_tau_daughter_phi);
    tree_->Branch("gen_tau_to_mu",&gen_tau_to_mu);
    tree_->Branch("gen_tau_to_3prong",&gen_tau_to_3prong);
    tree_->Branch("gen_tautau_to_mu3prong",&gen_tautau_to_mu3prong);

  }



}

//=================================================================================================================== 
void NtupleBranches::reset( void ){

  /** genParticle */
  genParticle_N = 0;
  genParticle_pt.clear();
  //  genParticle_px.clear();
  //  genParticle_py.clear();
  //  genParticle_pz.clear();
  //  genParticle_e.clear();
  genParticle_eta.clear();
  genParticle_phi.clear();
  genParticle_mass.clear();
  genParticle_pdgId.clear();
  genParticle_isPrompt.clear();
  genParticle_isDirectPromptTauDecayProduct.clear();
  genParticle_fromHardProcessFinalState.clear();
  genParticle_isDirectHardProcessTauDecayProductFinalState.clear();
  genParticle_status.clear();
  genParticle_mother.clear();
  genParticle_mother_pt.clear();
  genParticle_nMoth.clear();
  genParticle_nDau.clear();
  genParticle_dau.clear();
  
  /** generator info */
  genWeight   = 0;
  qScale      = 0;
  genFacWeightUp       = 0;
  genFacWeightDown     = 0;
  genRenWeightUp       = 0;
  genRenWeightDown     = 0;
  genFacRenWeightUp    = 0;
  genFacRenWeightDown  = 0;
  PDF_rms = 0;
  PDF_id.clear();  
  PDF_x.clear();	
  PDF_xPDF.clear();
  lheV_pt = 0;
  lheHT = 0;
  lheNj = 0;
  lheNb = 0;
  lheV_mass = 0;
    
 

  /** HLT trigger decisions */
  HLT_isFired.clear();

  /** HLT trigger objects */
  triggerObject_pt.clear();
  triggerObject_eta.clear();
  triggerObject_phi.clear();
  triggerObject_mass.clear();
  triggerObject_lastname.clear();
  triggerObject_filterIDs.clear();
  triggerObject_filterLabels.clear();
  triggerObject_firedTrigger.clear();

  /** HLT filter decisions */
  passFilter_HBHE_                  = false;
  passFilter_HBHELoose_             = false;
  passFilter_HBHETight_             = false;
  passFilter_HBHEIso_               = false;
  passFilter_CSCHalo_               = false;
  passFilter_CSCTightHalo2015_      = false;
  passFilter_HCALlaser_             = false;
  passFilter_ECALDeadCell_          = false;
  passFilter_GoodVtx_               = false;
  passFilter_TrkFailure_            = false;
  passFilter_EEBadSc_               = false;
  passFilter_ECALlaser_             = false;
  passFilter_TrkPOG_                = false;
  passFilter_TrkPOG_manystrip_      = false;
  passFilter_TrkPOG_toomanystrip_   = false;
  passFilter_TrkPOG_logError_       = false;
  passFilter_METFilters_            = false;
  //NEW FOR ICHEP
  passFilter_CSCTightHaloTrkMuUnvetoFilter_   = false;
  passFilter_globalTightHalo2016_             = false;
  passFilter_globalSuperTightHalo2016_             = false;
  passFilter_HcalStripHalo_                   = false;
  passFilter_chargedHadronTrackResolution_    = false;
  passFilter_muonBadTrack_                    = false;
  flag_badMuons_                    = false;
  flag_duplicateMuons_              = false;
  flag_nobadMuons_                  = false;

  /** energy density */
  rho = 0;
  
 

  /** MET */
  METraw_et.clear();	 
  METraw_phi.clear();
  METraw_sumEt.clear();
  MET_corrPx.clear();
  MET_corrPy.clear();
  MET_et.clear();
  MET_phi.clear();
  MET_puppi_et.clear();
  MET_puppi_phi.clear();

  MET_sumEt.clear();
  MET_T1Uncertainty.clear();
  
  MET_JetEnUp.clear();
  MET_JetEnDown.clear();
  MET_JetResUp.clear();
  MET_JetResDown.clear();
  MET_UnclusteredEnUp.clear();
  MET_UnclusteredEnDown.clear();

  /** MET SVift*/
  MET_significance.clear();
  MET_cov00.clear();
  MET_cov10.clear();
  MET_cov11.clear();
  MET_mva_et.clear();
  MET_mva_phi.clear();
  MET_mva_cov00.clear();
  MET_mva_cov10.clear();
  MET_mva_cov11.clear();
  MET_mva_recoil_pt.clear();
  MET_mva_recoil_eta.clear();
  MET_mva_recoil_phi.clear();
  MET_mva_recoil_pdgId.clear();
  MET_Nmva.clear();

  /*------------------------EVENT infos-------------------------*/    
  EVENT_event = 0;
  EVENT_run = 0;
  EVENT_lumiBlock = 0;

  /*--------------------------PV infos--------------------------*/
  PV_N = 0;
  PV_filter = true;
  PV_chi2.clear();
  PV_ndof.clear();
  PV_rho.clear();
  PV_z.clear();
  BeamSpot_x0.clear();
  BeamSpot_y0.clear();
  BeamSpot_z0.clear();
  /*--------------------------PU infos--------------------------*/  			       
  nPuVtxTrue.clear();
  nPuVtx.clear();
  bX.clear();

  /*-------------------------JPSI infos--------------------------*/ 
  IsBsTauTau.clear();


  BsTauTau_nCandidates.clear();

  BsTauTau_calo_eta.clear();      
  BsTauTau_calo_phi.clear();      
  BsTauTau_calo_energy.clear();   
  BsTauTau_calo_energyHFp.clear();
  BsTauTau_calo_energyHFm.clear();
  BsTauTau_calo_zdcSumPlus.clear();
  BsTauTau_calo_zdcSumMinus.clear();

  BsTauTau_trackPFactivity_pt.clear();
  BsTauTau_trackPFactivity_eta.clear();
  BsTauTau_trackPFactivity_phi.clear();

  BsTauTau_mu1_pt.clear();
  BsTauTau_mu1_eta.clear();
  BsTauTau_mu1_phi.clear();
  BsTauTau_mu1_mass.clear();
  BsTauTau_mu1_q.clear();
  BsTauTau_mu1_isLoose.clear();
  BsTauTau_mu1_isTight.clear();
  BsTauTau_mu1_isPF.clear();
  BsTauTau_mu1_isGlobal.clear();
  BsTauTau_mu1_isTracker.clear();
  BsTauTau_mu1_isSoft.clear();
  BsTauTau_mu1_vx.clear();
  BsTauTau_mu1_vy.clear();
  BsTauTau_mu1_vz.clear();
  BsTauTau_mu1_iso.clear();
  BsTauTau_mu1_dbiso.clear();

  BsTauTau_tau_pt.clear();
  BsTauTau_tau_eta.clear();
  BsTauTau_tau_phi.clear();
  BsTauTau_tau_mass.clear();
  BsTauTau_tau_rhomass1.clear();
  BsTauTau_tau_rhomass2.clear();
  BsTauTau_tau_q.clear();
  BsTauTau_tau_vx.clear();
  BsTauTau_tau_vy.clear();
  BsTauTau_tau_vz.clear();

  BsTauTau_tau_max_dr_3prong.clear();
  BsTauTau_tau_lip.clear();
  BsTauTau_tau_lips.clear();
  BsTauTau_tau_pvip.clear();
  BsTauTau_tau_pvips.clear();
  BsTauTau_tau_fl3d.clear();
  BsTauTau_tau_fls3d.clear();
  BsTauTau_tau_alpha.clear();
  BsTauTau_tau_vprob.clear();
  BsTauTau_tau_isRight.clear();
  BsTauTau_tau_isRight1.clear();
  BsTauTau_tau_isRight2.clear();
  BsTauTau_tau_isRight3.clear();
//  BsTauTau_tau_dr1.clear();
//  BsTauTau_tau_dr2.clear();
//  BsTauTau_tau_dr3.clear();
//  BsTauTau_tau_ptres1.clear();
//  BsTauTau_tau_ptres2.clear();
//  BsTauTau_tau_ptres3.clear();
  BsTauTau_tau_matched_ppdgId.clear();
  BsTauTau_tau_matched_gentaupt.clear();
  BsTauTau_tau_pfidx1.clear();
  BsTauTau_tau_pfidx2.clear();
  BsTauTau_tau_pfidx3.clear();
  
  reco_pion_pt.clear();
  reco_pion_eta.clear();
  reco_pion_phi.clear();

  BsTauTau_tau_pi1_pt.clear();
  BsTauTau_tau_pi1_eta.clear();
  BsTauTau_tau_pi1_phi.clear();
  BsTauTau_tau_pi1_mass.clear();
  BsTauTau_tau_pi1_x.clear();
  BsTauTau_tau_pi1_y.clear();
  BsTauTau_tau_pi1_z.clear();
  BsTauTau_tau_pi2_pt.clear();
  BsTauTau_tau_pi2_eta.clear();
  BsTauTau_tau_pi2_phi.clear();
  BsTauTau_tau_pi2_mass.clear();
  BsTauTau_tau_pi2_x.clear();
  BsTauTau_tau_pi2_y.clear();
  BsTauTau_tau_pi2_z.clear();
  BsTauTau_tau_pi3_pt.clear();
  BsTauTau_tau_pi3_eta.clear();
  BsTauTau_tau_pi3_phi.clear();
  BsTauTau_tau_pi3_mass.clear();
  BsTauTau_tau_pi3_x.clear();
  BsTauTau_tau_pi3_y.clear();
  BsTauTau_tau_pi3_z.clear();

  BsTauTau_B_pt.clear();
  BsTauTau_B_eta.clear();
  BsTauTau_B_phi.clear();
  BsTauTau_B_mass.clear();
  BsTauTau_B_vprob.clear();
  BsTauTau_B_lip.clear();
  BsTauTau_B_lips.clear();
  BsTauTau_B_pvip.clear();
  BsTauTau_B_pvips.clear();
  BsTauTau_B_fl3d.clear();
  BsTauTau_B_fls3d.clear();
  BsTauTau_B_alpha.clear();
  BsTauTau_B_maxdoca.clear();
  BsTauTau_B_mindoca.clear();
  BsTauTau_B_vx.clear();
  BsTauTau_B_vy.clear();
  BsTauTau_B_vz.clear();
  BsTauTau_B_iso.clear();
  BsTauTau_B_iso_ntracks.clear();
  BsTauTau_B_iso_mindoca.clear();

  BsTauTau_PV_vx.clear();
  BsTauTau_PV_vy.clear();
  BsTauTau_PV_vz.clear();

  BsTauTau_bbPV_vx.clear();
  BsTauTau_bbPV_vy.clear();
  BsTauTau_bbPV_vz.clear();

  BsTauTau_bbPV_refit_vx.clear();
  BsTauTau_bbPV_refit_vy.clear();
  BsTauTau_bbPV_refit_vz.clear();

  BsTauTau_genPV_vx.clear();
  BsTauTau_genPV_vy.clear();
  BsTauTau_genPV_vz.clear();

  BsTauTau_ngenmuons.clear();

  BsTauTau_isgen3.clear();
  BsTauTau_isgen3matched.clear();
  BsTauTau_tau_rEff.clear();
  BsTauTau_tau_rEff_error.clear();
  BsTauTau_nPions.clear();
  BsTauTau_nch.clear();
  BsTauTau_nch_qr.clear();
  BsTauTau_ngentau3.clear();
  BsTauTau_ngentau.clear();
  BsTauTau_gentaupt.clear();
  BsTauTau_gentaudm.clear();
  
  PVz.clear();
  pion_z.clear();
  matched_pion_deltaPt.clear();
  
  triggered.clear();
  gen_tau_pt.clear();
  gen_tau_eta.clear();
  gen_tau_phi.clear();
  gen_tau_daughter_pdgId.clear();
  gen_tau_daughter_pt.clear();
  gen_tau_daughter_eta.clear();
  gen_tau_daughter_phi.clear();
  gen_tau_to_mu.clear();
  gen_tau_to_3prong.clear();
  gen_tautau_to_mu3prong.clear();




  ////////////////////////


 
} 

void NtupleBranches::LabelHistograms( std::map< std::string, bool >& runFlags ){
  std::vector cutflow_bins_string = {"Precut", "Trigger","#tau #mu", "#tau hadron", "---","gen #mu","gen 3prong #tau","gen #mu-3prong event"};
  for(size_t i=0; i< cutflow_bins_string.size(); i++){
    cutflow_perevt->GetXaxis()->SetBinLabel(i+1, cutflow_bins_string[i]);
  }
  std::vector tau_decay_bins_string = {"number of #tau", "#mu", "e", "#pi^{#pm}", "#pi^{#pm} + #pi^{0}", "#pi^{#pm} + multiple #pi^{0}s", "3 prong inclusive", "a1", "#rho^{#pm}", "W^{#pm}", "#gamma", "The rest"};
  for(size_t i=0; i< tau_decay_bins_string.size(); i++){
    tau_decays->GetXaxis()->SetBinLabel(i+1, tau_decay_bins_string[i]);
  }


//  if ( runFlags["useHammer"] ){
//    std::vector bins_hammer = {"den", "num", "num_a0_up", "num_a0_down", "num_a1_up", "num_a1_down","num_a2_up", "num_a2_down", "num_b0_up", "num_b0_down", "num_b1_up", "num_b1_down","num_b2_up", "num_b2_down", "num_c1_up", "num_c1_down", "num_c2_up", "num_c2_down", "num_d0_up", "num_d0_down", "num_d1_up", "num_d1_down","num_d2_up", "num_d2_down"};
//    
//    for(size_t i=0; i< bins_hammer.size(); i++){
//      hammer_width->GetXaxis()->SetBinLabel(i+1, bins_hammer[i]);
//    }
//  }


  if ( runFlags["doGenHist"] ){
    std::vector bins_string = { "#mu", "#pi^{0}", "#pi^{#pm}","#rho^{0}","#rho^{#pm}","#eta","#eta^{`}","#omega","#phi","K^{0}","K^{#pm}","K^{*0}","K^{*#pm}","D^{#pm}","D^{0}","#eta_{c}","#eta_{b}","#Upsilon"};
    for(size_t i=0; i< bins_string.size(); i++){
      genParticle_Bdau_X_id->GetXaxis()->SetBinLabel(i+1, bins_string[i]);
    }       
  }
}
