#ifndef NtupleBranches_H
#define NtupleBranches_H

#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/Tau.h>
#include <DataFormats/MuonReco/interface/Muon.h>
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include <DataFormats/PatCandidates/interface/MET.h>
#include <DataFormats/PatCandidates/interface/Jet.h>
#include <DataFormats/JetReco/interface/Jet.h>
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "FWCore/Framework/interface/Event.h"
#include "TMatrixD.h"
#include "TMatrix.h"

#include "TTree.h"
#include "TLorentzVector.h"

/*here we declare the input and output variables*/

class NtupleBranches {

public:
  NtupleBranches( std::map< std::string, bool >& runFlags, TTree* tree = 0);
  ~NtupleBranches( void );
   
  void branch( std::map< std::string, bool >& runFlags );
  void getEventByLabels( edm::EventBase const & event );
  void reset( void );
  void fillTree( void ){ tree_->Fill(); };
  void LabelHistograms (std::map< std::string, bool >& runFlags );
  //=================================================================================================================== 

  /* output histogram */
  /* Cutflow */
  TH1F* cutflow_perevt ;
  TH1F* tau_decays;
  TH1F* tau_pT;
  TH1F* q2_nocut;
//  TH1F* hammer_width;
  
  /* Histogram for genParticles */
   TH1F* genParticle_Bdau_X_id;
   TH1F* genParticle_Bdau_X_pt; 
   TH1F* genParticle_Bdau_X_eta; 
   TH1F* genParticle_Bdau_X_phi; 
   TH1F* genParticle_Bdau_X_mass; 
   TH1F* genParticle_Bdau_mu1_pt;
   TH1F* genParticle_Bdau_mu1_eta;
   TH1F* genParticle_Bdau_mu1_phi; 
   TH1F* genParticle_Bdau_mu2_pt;
   TH1F* genParticle_Bdau_mu2_eta; 
   TH1F* genParticle_Bdau_mu2_phi; 
   TH1F* genParticle_Bdau_Jpsi_pt; 
   TH1F* genParticle_Bdau_Jpsi_eta; 
   TH1F* genParticle_Bdau_Jpsi_phi; 
   TH1F* genParticle_Bdau_Jpsi_mass;
   TH1F* genParticle_Bvis_pt; 
   TH1F* genParticle_Bvis_eta;
   TH1F* genParticle_Bvis_phi; 
   TH1F* genParticle_Bvis_mass;

  /* output tree variables*/
    
  /** genParticles */
  int                             genParticle_N;
  std::vector<float>              genParticle_pt       ;
//  std::vector<float>              genParticle_px       ;
//  std::vector<float>              genParticle_py       ;
//  std::vector<float>              genParticle_pz       ;
//  std::vector<float>              genParticle_e        ;
  std::vector<float>              genParticle_eta      ;
  std::vector<float>              genParticle_phi      ;
  std::vector<float>              genParticle_mass     ;
  std::vector<int  >              genParticle_pdgId    ;
  std::vector<int  >              genParticle_isPrompt ;
  std::vector<int  >              genParticle_isDirectPromptTauDecayProduct;
  std::vector<int  >              genParticle_fromHardProcessFinalState;
  std::vector<int  >              genParticle_isDirectHardProcessTauDecayProductFinalState;
  std::vector<int  >              genParticle_status   ;
  std::vector<int  >              genParticle_nDau     ;
  std::vector<int  >              genParticle_nMoth    ;
  std::vector<std::vector<int> >  genParticle_mother   ; 
  std::vector<std::vector<float>> genParticle_mother_pt;
  std::vector<std::vector<int> >  genParticle_dau      ;

  /** generator info */
  float                           lheV_pt              ;
  float                           lheHT                ;
  int                             lheNj                ;
  int                             lheNb                ;
  int                             lheNl                ;
  float                           lheV_mass            ;
  float                           genWeight            ;
  float                           genFacWeightUp       ;
  float                           genFacWeightDown     ;
  float                           genRenWeightUp       ;
  float                           genRenWeightDown     ;
  float                           genFacRenWeightUp    ;
  float                           genFacRenWeightDown  ;
  float                           qScale               ;
  float                           PDF_rms              ;
  std::vector<int  >              PDF_id               ;
  std::vector<float>              PDF_x                ;
  std::vector<float>              PDF_xPDF             ;

          
 


  
   
  
  /** HLT trigger decisions */
  std::map<std::string,bool> HLT_isFired;
	 
  /** HLT trigger objects */
  std::vector<float>  		    triggerObject_pt	      ;
  std::vector<float>  		    triggerObject_eta	      ;
  std::vector<float>  		    triggerObject_phi	      ;
  std::vector<float>  		    triggerObject_mass	      ;
  std::vector<std::string>  		    triggerObject_lastname    ;
  std::vector< std::vector<float> > triggerObject_filterIDs   ; // as defined in http://cmslxr.fnal.gov/lxr/source/DataFormats/HLTReco/interface/TriggerTypeDefs.h
  //  std::vector< std::vector<std::string> > triggerObject_filterLabels;
  std::map<std::string, std::vector<std::string> > triggerObject_filterLabels;
  std::vector< std::vector<int> >   triggerObject_firedTrigger; // as defined in plugins/TriggersNtuplizer.cc


  /** HLT filter decisions */
  bool passFilter_HBHE_;
  bool passFilter_HBHELoose_;
  bool passFilter_HBHETight_;
  bool passFilter_HBHEIso_;
  bool passFilter_CSCHalo_;
  bool passFilter_CSCTightHalo2015_;
  bool passFilter_HCALlaser_;
  bool passFilter_ECALDeadCell_;
  bool passFilter_GoodVtx_;
  bool passFilter_TrkFailure_;
  bool passFilter_EEBadSc_;
  bool passFilter_ECALlaser_;
  bool passFilter_TrkPOG_;
  bool passFilter_TrkPOG_manystrip_;
  bool passFilter_TrkPOG_toomanystrip_;
  bool passFilter_TrkPOG_logError_;
  bool passFilter_METFilters_;
  bool passFilter_CSCTightHaloTrkMuUnvetoFilter_   ;
  bool passFilter_globalTightHalo2016_             ;
  bool passFilter_globalSuperTightHalo2016_             ;
  bool passFilter_HcalStripHalo_                   ;
  bool passFilter_chargedHadronTrackResolution_    ;
  bool passFilter_muonBadTrack_                    ;
  bool flag_badMuons_;
  bool flag_duplicateMuons_;
  bool flag_nobadMuons_;
  bool passFilter_ecalBadCalib_;

  /** MET */
  /** energy density */
  float                             rho;
 
  std::vector<float>                METraw_et		      ;	 
  std::vector<float>                METraw_phi		      ;
  std::vector<float>  	            METraw_sumEt	      ;
  std::vector<float>  	      	    MET_corrPx  	      ;
  std::vector<float>  	            MET_corrPy  	      ;
  std::vector<float>  	      	    MET_et		      ;
  std::vector<float>  	      	    MET_phi		      ;
  std::vector<float>  	      	    MET_puppi_et		      ;
  std::vector<float>  	      	    MET_puppi_phi	      ;
  std::vector<float>  	      	    MET_mva_et		      ;
  std::vector<float>  	      	    MET_mva_phi	              ;
  std::vector<float>  	      	    MET_sumEt		      ;
  std::vector<float>  	      	    MET_T1Uncertainty	      ;
  
  std::vector<float>  	      	    MET_JetEnUp	      ;
  std::vector<float>  	      	    MET_JetEnDown	      ;
  std::vector<float>  	      	    MET_JetResUp	      ;
  std::vector<float>  	      	    MET_JetResDown	      ;
  std::vector<float>  	      	    MET_UnclusteredEnUp	      ;
  std::vector<float>  	      	    MET_UnclusteredEnDown	      ;
  

  /** MET SVift*/
  std::vector<float>                MET_significance	      ;	 
  std::vector<float>                MET_cov00	      ;	 
  std::vector<float>                MET_cov10	      ;	 
  std::vector<float>                MET_cov11	      ;	 
  std::vector<float>                MET_mva_cov00	      ;	 
  std::vector<float>                MET_mva_cov10	      ;	 
  std::vector<float>                MET_mva_cov11	      ;	 
  std::vector< std::vector<float> > MET_mva_recoil_pt;
  std::vector< std::vector<float> > MET_mva_recoil_eta;
  std::vector< std::vector<float> > MET_mva_recoil_phi;
  std::vector< std::vector<int> >   MET_mva_recoil_pdgId;
  std::vector<int>  	            MET_Nmva;

  /*------------------------EVENT infos-------------------------*/    
  int                               EVENT_event            ;
  int                               EVENT_run              ;
  int                               EVENT_lumiBlock        ;
  /*-------------------------JPSI VARIABLES---------------------*/
  std::vector<int  >                IsBsTauTau   ;

  std::vector <double>              BsTauTau_calo_eta;
  std::vector <double>              BsTauTau_calo_phi;
  std::vector <double>              BsTauTau_calo_energy;
  std::vector <double>              BsTauTau_calo_energyHFp;
  std::vector <double>              BsTauTau_calo_energyHFm;
  std::vector <double>              BsTauTau_calo_zdcSumPlus;
  std::vector <double>              BsTauTau_calo_zdcSumMinus;

  std::vector <double>              BsTauTau_trackPFactivity_pt;
  std::vector <double>              BsTauTau_trackPFactivity_eta;
  std::vector <double>              BsTauTau_trackPFactivity_phi;

  std::vector<int  >                BsTauTau_nCandidates ;

  std::vector<float>                BsTauTau_mu1_pt      ;
  std::vector<float>                BsTauTau_mu1_eta     ;
  std::vector<float>                BsTauTau_mu1_phi     ;
  std::vector<float>                BsTauTau_mu1_mass     ;
  std::vector<float>                BsTauTau_mu1_unfit_pt      ;
  std::vector<float>                BsTauTau_mu1_unfit_eta     ;
  std::vector<float>                BsTauTau_mu1_unfit_phi     ;
  std::vector<float>                BsTauTau_mu1_unfit_mass     ;
  std::vector<int  >                BsTauTau_mu1_q      ;   
  std::vector<int  >                BsTauTau_mu1_isLoose   ;
  std::vector<int  >                BsTauTau_mu1_isTight   ;
  std::vector<int  >                BsTauTau_mu1_isPF      ;
  std::vector<int  >                BsTauTau_mu1_isGlobal  ;
  std::vector<int  >                BsTauTau_mu1_isTracker ;
  std::vector<int  >                BsTauTau_mu1_isSoft    ;
  std::vector<float>                BsTauTau_mu1_vx       ;
  std::vector<float>                BsTauTau_mu1_vy       ;
  std::vector<float>                BsTauTau_mu1_vz       ;
  std::vector<float>                BsTauTau_mu1_iso       ;
  std::vector<float>                BsTauTau_mu1_dbiso       ;

  std::vector<float>                BsTauTau_tau_pt      ;
  std::vector<float>                BsTauTau_tau_eta     ;
  std::vector<float>                BsTauTau_tau_phi     ;
  std::vector<float>                BsTauTau_tau_mass     ;
  std::vector<float>                BsTauTau_tau_rhomass1     ;
  std::vector<float>                BsTauTau_tau_rhomass2     ;
  std::vector<int  >                BsTauTau_tau_q   ;   
  std::vector<float>                BsTauTau_tau_vx       ;
  std::vector<float>                BsTauTau_tau_vy       ;
  std::vector<float>                BsTauTau_tau_vz       ;

  std::vector<float>       BsTauTau_tau_max_dr_3prong;
  std::vector<float>       BsTauTau_tau_lip;
  std::vector<float>       BsTauTau_tau_lips;
  std::vector<float>       BsTauTau_tau_pvip;
  std::vector<float>       BsTauTau_tau_pvips;
  std::vector<float>       BsTauTau_tau_fl3d;
  std::vector<float>       BsTauTau_tau_fls3d;
  std::vector<float>       BsTauTau_tau_alpha;
  std::vector<float>       BsTauTau_tau_vprob;
  std::vector<bool>        BsTauTau_tau_isRight;
  std::vector<bool>        BsTauTau_tau_isRight1;
  std::vector<bool>        BsTauTau_tau_isRight2;
  std::vector<bool>        BsTauTau_tau_isRight3;
//  std::vector<float>        BsTauTau_tau_dr1;
//  std::vector<float>        BsTauTau_tau_dr2;
//  std::vector<float>        BsTauTau_tau_dr3;
//  std::vector<float>        BsTauTau_tau_ptres1;
//  std::vector<float>        BsTauTau_tau_ptres2;
//  std::vector<float>        BsTauTau_tau_ptres3;
  std::vector<int>         BsTauTau_tau_matched_ppdgId;
  std::vector<float>       BsTauTau_tau_matched_gentaupt;
  std::vector<int>       BsTauTau_tau_pfidx1;
  std::vector<int>       BsTauTau_tau_pfidx2;
  std::vector<int>       BsTauTau_tau_pfidx3;


  std::vector<float>       reco_pion_pt;
  std::vector<float>       reco_pion_eta;
  std::vector<float>       reco_pion_phi;

  std::vector<float>       BsTauTau_tau_pi1_pt;
  std::vector<float>       BsTauTau_tau_pi1_eta;
  std::vector<float>       BsTauTau_tau_pi1_phi;
  std::vector<float>       BsTauTau_tau_pi1_mass;
  std::vector<float>       BsTauTau_tau_pi1_x;
  std::vector<float>       BsTauTau_tau_pi1_y;
  std::vector<float>       BsTauTau_tau_pi1_z;
  std::vector<float>       BsTauTau_tau_pi2_pt;
  std::vector<float>       BsTauTau_tau_pi2_eta;
  std::vector<float>       BsTauTau_tau_pi2_phi;
  std::vector<float>       BsTauTau_tau_pi2_mass;
  std::vector<float>       BsTauTau_tau_pi2_x;
  std::vector<float>       BsTauTau_tau_pi2_y;
  std::vector<float>       BsTauTau_tau_pi2_z;
  std::vector<float>       BsTauTau_tau_pi3_pt;
  std::vector<float>       BsTauTau_tau_pi3_eta;
  std::vector<float>       BsTauTau_tau_pi3_phi;
  std::vector<float>       BsTauTau_tau_pi3_mass;
  std::vector<float>       BsTauTau_tau_pi3_x;
  std::vector<float>       BsTauTau_tau_pi3_y;
  std::vector<float>       BsTauTau_tau_pi3_z;

  std::vector<float>                BsTauTau_B_pt      ;
  std::vector<float>                BsTauTau_B_eta     ;
  std::vector<float>                BsTauTau_B_phi     ;
  std::vector<float>                BsTauTau_B_mass    ;
  std::vector<float>                BsTauTau_B_vprob ;
  std::vector<float>                BsTauTau_B_lip;
  std::vector<float>                BsTauTau_B_lips;
  std::vector<float>                BsTauTau_B_pvip;
  std::vector<float>                BsTauTau_B_pvips;
  std::vector<float>                BsTauTau_B_fl3d;
  std::vector<float>                BsTauTau_B_fls3d;
  std::vector<float>                BsTauTau_B_alpha;
  std::vector<float>                BsTauTau_B_maxdoca;
  std::vector<float>                BsTauTau_B_mindoca;
  std::vector<float>                BsTauTau_B_vx      ;
  std::vector<float>                BsTauTau_B_vy      ;
  std::vector<float>                BsTauTau_B_vz      ;
  std::vector<float>                BsTauTau_B_iso;
  std::vector<int  >                BsTauTau_B_iso_ntracks;
  std::vector<float>                BsTauTau_B_iso_mindoca;
  std::vector<float>                BsTauTau_B_unfit_pt      ;
  std::vector<float>                BsTauTau_B_unfit_mass       ;
  std::vector<float>                BsTauTau_B_unfit_vprob    ;
  std::vector<float>                BsTauTau_B_unfit_vx;
  std::vector<float>                BsTauTau_B_unfit_vy;
  std::vector<float>                BsTauTau_B_unfit_vz;

  std::vector<float>                BsTauTau_PV_vx       ;
  std::vector<float>                BsTauTau_PV_vy       ;
  std::vector<float>                BsTauTau_PV_vz       ;

  std::vector<float>                BsTauTau_bbPV_vx       ;
  std::vector<float>                BsTauTau_bbPV_vy       ;
  std::vector<float>                BsTauTau_bbPV_vz       ;

  std::vector<float>                BsTauTau_bbPV_refit_vx       ;
  std::vector<float>                BsTauTau_bbPV_refit_vy       ;
  std::vector<float>                BsTauTau_bbPV_refit_vz       ;

  std::vector<float>                BsTauTau_genPV_vx       ;
  std::vector<float>                BsTauTau_genPV_vy       ;
  std::vector<float>                BsTauTau_genPV_vz       ;

  std::vector<int  >                BsTauTau_ngenmuons      ;
  std::vector<int  >                BsTauTau_isgen3;
  std::vector<int  >                BsTauTau_isgen3matched;
  std::vector<float> BsTauTau_tau_rEff;
  std::vector<float> BsTauTau_tau_rEff_error;
  std::vector<int> BsTauTau_nPions;
  std::vector<int> BsTauTau_nch;
  std::vector<int> BsTauTau_nch_qr;
  std::vector<int> BsTauTau_ngentau3;
  std::vector<int> BsTauTau_ngentau;
  std::vector<float> BsTauTau_gentaupt;
  std::vector<int> BsTauTau_gentaudm;
  
  std::vector<float> PVz;
  std::vector<float> pion_z;
  std::vector<float> matched_pion_deltaPt;
  
  std::vector<bool> triggered;
  std::vector<float> gen_tau_pt;
  std::vector<float> gen_tau_eta;
  std::vector<float> gen_tau_phi;
  std::vector<int> gen_tau_daughter_pdgId;
  std::vector<float> gen_tau_daughter_pt;
  std::vector<float> gen_tau_daughter_eta;
  std::vector<float> gen_tau_daughter_phi;
  std::vector<bool> gen_tau_to_mu;
  std::vector<bool> gen_tau_to_3prong;
  std::vector<bool> gen_tautau_to_mu3prong;





  //////////////////////////




  /** HLT trigger objects */
  /* std::vector<float>  		    triggerObject_pt	      ; */
  /* std::vector<float>  		    triggerObject_eta	      ; */
  /* std::vector<float>  		    triggerObject_phi	      ; */
  /* std::vector<float>  		    triggerObject_mass	      ; */
  /* std::vector<std::string>  		    triggerObject_lastname    ; */
  /* std::vector< std::vector<float> > triggerObject_filterIDs   ; // as defined in http://cmslxr.fnal.gov/lxr/source/DataFormats/HLTReco/interface/TriggerTypeDefs.h */
  /* //  std::vector< std::vector<std::string> > triggerObject_filterLabels; */
  /* std::map<std::string, std::vector<std::string> > triggerObject_filterLabels; */
  /* std::vector< std::vector<int> >   triggerObject_firedTrigger; // as defined in plugins/TriggersNtuplizer.cc */


  /*--------------------------PV infos--------------------------*/
  int                               PV_N		     ;
  bool                              PV_filter		 ;
  std::vector<float>                PV_chi2          ;
  std::vector<float>                PV_ndof          ;
  std::vector<float>                PV_rho           ;
  std::vector<float>                PV_z             ;
  std::vector<float>                BeamSpot_x0;
  std::vector<float>                BeamSpot_y0;
  std::vector<float>                BeamSpot_z0;
  /*--------------------------PU infos--------------------------*/  			       
  std::vector<float  >                nPuVtxTrue             ;// the *true* mean number of the poisson distribution for this event from which the number of interactions each bunch crossing has been sampled // In MC it can be float 
  std::vector<int  >                nPuVtx                 ;// the number of pileup interactions that have been added to the event in the current bunch crossing
  std::vector<int  >                bX                     ;// to which bunch crossing do these interaction belong?  
  
private:
  TTree* tree_;

};

#endif // NtupleBranches_H
