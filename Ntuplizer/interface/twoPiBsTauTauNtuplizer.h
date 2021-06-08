#ifndef HELPER_H
#include "../interface/helper.h"
#endif

#ifndef twoPiBsTauTauNtuplizer_H
#define twoPiBsTauTauNtuplizer_H

class twoPiBsTauTauNtuplizer : public CandidateNtuplizer {


 public:
  twoPiBsTauTauNtuplizer( edm::EDGetTokenT<reco::MuonCollection>    muonToken   , 
		     edm::EDGetTokenT<edm::SortedCollection<CaloTower>> CaloTowerCollection ,
         edm::EDGetTokenT<reco::Centrality> tok_centSrc,
		     edm::EDGetTokenT<reco::VertexCollection> verticeToken, 
		     edm::EDGetTokenT<std::vector<reco::PFCandidate>> packedpfcandidatesToken,
		     edm::EDGetTokenT<edm::TriggerResults> triggertoken,
		     edm::EDGetTokenT<reco::GenParticleCollection> genptoken, 
		     edm::EDGetTokenT<std::vector<reco::GenJet>> genttoken,
		     std::map< std::string, bool >& runFlags,
		     std::map< std::string, double >& runValues,
		     std::map< std::string, std::string >& runStrings,
		     NtupleBranches* nBranches );
  
  ~twoPiBsTauTauNtuplizer( void );

  bool fillBranches( edm::Event const & event, const edm::EventSetup& iSetup );
  
private:
   edm::EDGetTokenT<reco::MuonCollection>    muonToken_   ;
   edm::EDGetTokenT<edm::SortedCollection<CaloTower>> CaloTowerCollection_;
   edm::EDGetTokenT<reco::Centrality> tok_centSrc_;
   edm::EDGetTokenT<reco::VertexCollection> verticeToken_   ;
   edm::EDGetTokenT<std::vector<reco::PFCandidate>>   		packedpfcandidatesToken_;
   edm::EDGetTokenT<edm::TriggerResults> 		     HLTtriggersToken_;
   edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
   edm::EDGetTokenT<std::vector<reco::GenJet>> genTauToken_;

   edm::Handle<std::vector<reco::Muon>>      		       muons_		       ;
   edm::Handle< reco::VertexCollection >  vertices_;
   edm::Handle< reco::BeamSpot >  beamspot_;
   edm::Handle< std::vector<reco::PFCandidate> > packedpfcandidates_   ;
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




#endif // twoPiBsTauTauNtuplizer_H

