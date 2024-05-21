#ifndef _JPsiTrk_h
#define _JPsiTrk_h

// system include files
#include <memory>

// user include files

#include "FWCore/Framework/interface/Frameworkfwd.h"
//#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"

//For kinematic fit:
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackFromFTSFactory.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "MagneticField/Engine/interface/MagneticField.h"            
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h" // muy importante para MiniAOD
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/V0Candidate/interface/V0Candidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"

#include "RecoVertex/VertexPrimitives/interface/BasicSingleVertexState.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"


#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"

//
// class decleration
//

//class JPsiTrk : public edm::EDAnalyzer {
class JPsiTrk : public edm::one::EDAnalyzer<> {  
public:
  explicit JPsiTrk(const edm::ParameterSet&);
  //~JPsiTrk();
  ~JPsiTrk() override = default;
  void fillPsi(const reco::Candidate& genpsi);
  void fillV0(const reco::Candidate& genv0);
  //int const getMuCat(reco::Muon const& muon) const;
  bool IsTheSame(const pat::GenericParticle& tk, const pat::Muon& mu);
  bool   isAncestor(const reco::Candidate*, const reco::Candidate*);
  bool isAncestor(int, const reco::Candidate * );
  double GetLifetime(TLorentzVector, TVector3, TVector3);
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  void printout(const RefCountedKinematicVertex& myVertex) const;
  void printout(const RefCountedKinematicParticle& myParticle) const;
  void printout(const RefCountedKinematicTree& myTree) const;

  // ----------member data ---------------------------
  //const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttrkToken_;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttrkToken_;
  edm::EDGetTokenT<edm::View<pat::Muon>> dimuon_Label;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate>> trakCollection_label;
  edm::EDGetTokenT<reco::GenParticleCollection> genCands_;
  edm::EDGetTokenT<pat::PackedGenParticleCollection> packedGenToken_;
  edm::EDGetTokenT<reco::VertexCollection> primaryVertices_Label;
  edm::EDGetTokenT<edm::TriggerResults> triggerResults_Label;
  edm::EDGetTokenT<reco::BeamSpot> BSLabel_;

  bool OnlyBest_;
  bool isMC_;
  bool OnlyGen_;
  bool isRes_;
  bool mumuMassConstraint_;
  std::vector<double> mumuMasscut_;
  double Trkmass_;
  int  Trkpdgid_;
  int  Bpdgid_;
  int  Respdgid_;
  std::vector<double> BarebMasscut_;
  std::vector<double> bMasscut_;

  TTree*      tree_;
  
  Double_t    mumC2;
  int         mumNHits, mumNPHits; 
  Double_t    mupC2;
  int         mupNHits, mupNPHits;
  Double_t    mumdxy, mupdxy, mumdz, mupdz;
  Double_t    muon_dca;

  int         tri_Dim25, tri_JpsiTrk_Bc, tri_JpsiTk; 
  int         tri_DMu4_3_LM;//HLT_DoubleMu4_3_LowMass
  int         tri_DMu4_LM_Displaced;//HLT_DoubleMu4_LowMass_Displaced

  bool       mu1soft, mu2soft, mu1tight, mu2tight;  
  bool       mu1PF, mu2PF, mu1loose, mu2loose;  
 
  // *************************************
  unsigned int    nB;
  unsigned int    nMu;
    
  Double_t       B_mass, B_px, B_py, B_pz, B_charge;
  Double_t       B_k_px, B_k_py, B_k_pz,  B_k_charge1; 
  Double_t       B_k_px_track, B_k_py_track, B_k_pz_track;

  Double_t       B_J_mass, B_J_massErr, B_J_px, B_J_py, B_J_pz;

  Double_t       B_J_pt1, B_J_px1, B_J_py1, B_J_pz1;
  Double_t       B_J_pt2, B_J_px2, B_J_py2, B_J_pz2;
  int            B_J_charge1, B_J_charge2;

  // Primary Vertex (PV)
  UInt_t          nVtx;
  Double_t        priVtxX, priVtxY, priVtxZ, priVtxXE, priVtxYE, priVtxZE, priVtxCL;
  Double_t        priVtxXYE, priVtxXZE, priVtxYZE;
  
  // ********************************** ************************************************************************
 
  Double_t      B_chi2, B_J_chi2;
  Double_t      B_Prob, B_J_Prob;

  Double_t      B_DecayVtxX,  B_DecayVtxY,  B_DecayVtxZ;
  Double_t      B_DecayVtxXE, B_DecayVtxYE, B_DecayVtxZE;
  Double_t      B_DecayVtxXYE, B_DecayVtxXZE, B_DecayVtxYZE;

  UInt_t  run;
  ULong64_t event;
  UInt_t lumiblock;

  TLorentzVector gen_bc_p4,gen_jpsi_p4,gen_pion3_p4,gen_muon1_p4,gen_muon2_p4;
  TVector3       gen_bc_vtx,gen_jpsi_vtx;
  Double_t       gen_bc_ct;

};
#endif
