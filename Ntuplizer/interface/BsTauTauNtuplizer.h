#ifndef HELPER_H
#include "../interface/helper.h"
#endif

#ifndef BsTauTauNtuplizer_H
#define BsTauTauNtuplizer_H

//#include "../interface/CandidateNtuplizer.h"
//#include "../interface/MyStruct.h"
////#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"
////#include "SimDataFormats/JetMatching/interface/JetFlavour.h"
//
//#include "DataFormats/VertexReco/interface/VertexFwd.h"
//#include "DataFormats/VertexReco/interface/Vertex.h"
////#include "DataFormats/PatCandidates/interface/Jet.h"
//#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
//#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
//
//#include "JetMETCorrections/Modules/interface/JetResolution.h"
////#include <CondFormats/JetMETObjects/interface/JetResolutionObject.h>
//
//#include "DataFormats/BeamSpot/interface/BeamSpot.h"
//
//// for vertexing
//#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"
//#include "EgammaAnalysis/ElectronTools/interface/ElectronEffectiveArea.h"
//
//#include "FWCore/Framework/interface/ESHandle.h"
//#include "FWCore/MessageLogger/interface/MessageLogger.h"
//#include "FWCore/Common/interface/TriggerNames.h"
//
//#include "DataFormats/Common/interface/Handle.h"  // edm::Handle
//#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"   
//#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
//#include "DataFormats/TrackReco/interface/DeDxData.h" 
//#include "DataFormats/MuonReco/interface/MuonSelectors.h"
//#include "DataFormats/MuonReco/interface/MuonQuality.h"
//
//#include "DataFormats/Math/interface/deltaR.h"
////#include "DataFormats/Common/interface/TriggerResults.h"
////#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
//#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
//#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
////#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
//
//#include "TrackingTools/IPTools/interface/IPTools.h"
//#include "TrackingTools/TransientTrack/interface/TransientTrack.h"             
//#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
//#include "TrackingTools/Records/interface/TransientTrackRecord.h"
//#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"
//#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
//#include "TrackingTools/PatternTools/interface/TransverseImpactPointExtrapolator.h"
//#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
//#include "MagneticField/Engine/interface/MagneticField.h"
//#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
//
//#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
//#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
//#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
//#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
//#include "RecoVertex/KinematicFit/interface/MultiTrackMassKinematicConstraint.h"
//#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
//#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
//#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
//#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
//#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
//#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
//#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
//
//#include "RecoVertex/PrimaryVertexProducer/interface/VertexHigherPtSquared.h"
//#include "RecoBTag/SecondaryVertex/interface/SecondaryVertex.h"
//#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"
//
//
//// kinematic fit package
//#include "DataFormats/TrackReco/interface/Track.h"
//#include "DataFormats/TrackReco/interface/TrackFwd.h"
//#include "FWCore/Framework/interface/MakerMacros.h"
//
//#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
//#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
//#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
//
//
//#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
//
//#include "PhysicsTools/JetMCUtils/interface/JetMCTag.h"
//
//#include "FWCore/ParameterSet/interface/FileInPath.h"
//#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"
//
//#include <tuple>
//#include <sstream>
//#include <iostream>
//#include <cstdlib>
//#include <fstream>
//#include <string>
//#include <algorithm> 
//#include <stdio.h>
//#include <stdlib.h>
//#include <vector>
//#include <map>

//struct particle_cand 
//{
//
//  // Longitudinal impact parameter and its significance
//  Float_t lip;
//  Float_t lips;
//  
//  // Impact parameter for the PV and its significance
//  Float_t pvip;
//  Float_t pvips;
//
//  // Flight length and its significance
//  Float_t fl3d;
//  Float_t fls3d;
//  
//  // opening angle
//  Float_t alpha;
//
//};

//#include "../interface/helper.h"

class BsTauTauNtuplizer : public CandidateNtuplizer {


 public:
  BsTauTauNtuplizer( edm::EDGetTokenT<reco::MuonCollection>    muonToken   , 
		     edm::EDGetTokenT<edm::SortedCollection<CaloTower>> CaloTowerCollection ,
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
		     NtupleBranches* nBranches );
  
  ~BsTauTauNtuplizer( void );

//  std::tuple<Float_t, TransientVertex> vertexProb( const std::vector<reco::TransientTrack>& tracks);
//
//  particle_cand calculateIPvariables(AnalyticalImpactPointExtrapolator extrapolator,
//				     RefCountedKinematicParticle particle,
//				     RefCountedKinematicVertex vertex,
//				     reco::Vertex wrtVertex);
//  
//  std::pair<bool, Measurement1D> absoluteImpactParameter(const TrajectoryStateOnSurface& tsos,
//							 RefCountedKinematicVertex vertex,
//							 VertexDistance& distanceComputer);
//
//
//  math::PtEtaPhiMLorentzVector daughter_p4(std::vector< RefCountedKinematicParticle > fitted_children, size_t i);
//
//  float MuonPFIso(pat::Muon muon);
  bool fillBranches( edm::Event const & event, const edm::EventSetup& iSetup );
//  Float_t getMaxDoca(std::vector<RefCountedKinematicParticle> &kinParticles);
//  Float_t getMinDoca(std::vector<RefCountedKinematicParticle> &kinParticles);
//  TVector3 getVertex(const reco::GenParticle& part);
//
//  void printout(const RefCountedKinematicVertex& myVertex);
//  void printout(const RefCountedKinematicParticle& myParticle);
//  void printout(const RefCountedKinematicTree& myTree);
//
//  Int_t decaymode_id(std::string str);
  
private:
   edm::EDGetTokenT<reco::MuonCollection>    muonToken_   ;
   edm::EDGetTokenT<edm::SortedCollection<CaloTower>> CaloTowerCollection_;
   edm::EDGetTokenT<reco::Centrality> tok_centSrc_;
   edm::EDGetTokenT<reco::VertexCollection> verticeToken_   ;
   edm::EDGetTokenT<std::vector<reco::PFCandidate>>   		packedpfcandidatesToken_;
   //edm::EDGetTokenT<std::vector<reco::TrackCollection>>   		TrackCollectionToken_;
   edm::EDGetTokenT<edm::TriggerResults> 		     HLTtriggersToken_;
   edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
   edm::EDGetTokenT<std::vector<reco::GenJet>> genTauToken_;

   edm::Handle<reco::Centrality> cent;
   edm::Handle<std::vector<reco::Muon>>      		       muons_		       ;
   edm::Handle< reco::VertexCollection >  vertices_;
   edm::Handle< reco::BeamSpot >  beamspot_;
   edm::Handle< std::vector<reco::PFCandidate> > packedpfcandidates_   ;
   //edm::Handle< std::vector<reco::TrackCollection> > TrackCollection_   ;
   edm::Handle< edm::TriggerResults> 			     HLTtriggers_;
   edm::Handle< reco::GenParticleCollection >  genParticles_;
   edm::Handle< std::vector<reco::GenJet> >  genTaus_;

   edm::ESHandle<TransientTrackBuilder> builder;

   const MagneticField                 *fMagneticField;

//   ParticleMass muon_mass = 0.1056583;
//   ParticleMass jpsi_mass = 3.09687;
//   ParticleMass pion_mass = 0.139571;
//   ParticleMass kaon_mass = 0.493677;
//   
//   float muon_sigma = 0.0000001;
//   float jp_m_sigma = 0.00004;
//   float pion_sigma = 0.000016;
//   float kaon_sigma = 0.000016;
   float chi = 0.;
   float ndf = 0.;

   bool runOnMC_;   
   bool useDNN_;

   float c_dz;
   float c_fsig;
   float c_vprob;
   helper aux;

   std::vector<std::string> detNames = { "EB", "EE", "HB", "HE", "HFm", "HFp", "unknown" };
   std::vector<std::pair<double,double>> detLimits = {
     {0    , 1.4442 }, // EB (Exclude transition region between calo barrel and endcap)
     {1.566, 3.0   }, // EE (Exclude transition region between calo barrel and endcap)
     {0    , 1.3   }, // HB
     {1.3  , 3.0   }, // HE
     {3.0  , 5.2   }, // HFm    // sleontsi
     {3.0  , 5.2   }, // HFp*/
//     {-5.2 ,-3.0   }, // HFm    // original
//     {3.0  , 5.2   }, // HFp*/

   };
   enum ECaloType { kEB, kEE, kHB, kHE, kHFp, kHFm, nCaloTypes };
   ECaloType GetTowerSubdetHad(double&) const;
   ECaloType GetTowerSubdetEm(double&)  const;

   const std::map<ECaloType, std::string> caloName = {
     { kEB  , "EB"  },
     { kEE  , "EE"  },
     { kHB  , "HB"  },
     { kHE  , "HE"  },
     { kHFp , "HFp" },
     { kHFm , "HFm" },
   };
 
   const std::map<ECaloType, double> noiseThreshold = {
     { kEB  , 0.7 },
     { kEE  , 7.5 },
     { kHB  , 2.8  },
     { kHE  , 2.4  },
     { kHFp , 7.2  },
     { kHFm , 7.5  },
   };
 
   const double maxEtaEB = detLimits.at(0).second;
   const double minEtaEE = detLimits.at(1).first;
   const double maxEtaEE = detLimits.at(1).second;
 
   const double maxEtaHB = detLimits.at(2).second;
   const double minEtaHE = detLimits.at(3).first;
   const double maxEtaHE = detLimits.at(3).second;
 
   const double minEtaHF = abs(detLimits.at(4).first);
   const double maxEtaHF = abs(detLimits.at(4).second);

};




#endif // BsTauTauNtuplizer_H

