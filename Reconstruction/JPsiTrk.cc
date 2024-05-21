// -*- C++ -*-
//
// Package:    JPsiTrk
// Class:      JPsiTrk
// 

//=================================================
// Original author:  Jhovanny Andres Mejia        |
//         created:  October of 2021              |
//         <jhovanny.andres.mejia.guisao@cern.ch> 
// MC non Reonant Channel Edit: Oscar Isaac Pérez |
//         added: Fevruary 2024                   |
//=================================================

// user include files
#include "myAnalyzers/JPsiKsPAT/src/JPsiTrk.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TH2F.h"

//
// constants, enums and typedefs
//

  typedef math::Error<3>::type CovarianceMatrix;

//
// static data member definitions
//

//
// constructors and destructor
//

JPsiTrk::JPsiTrk(const edm::ParameterSet& iConfig)
  :
  //ttrkToken_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
  ttrkToken_(esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder"))),
  dimuon_Label(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("dimuons"))),
  trakCollection_label(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("Trak"))),
  genCands_(consumes<reco::GenParticleCollection>(iConfig.getParameter < edm::InputTag > ("GenParticles"))), 
  packedGenToken_(consumes<pat::PackedGenParticleCollection>(iConfig.getParameter <edm::InputTag> ("packedGenParticles"))), 
  primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
  triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
  BSLabel_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("bslabel"))),

  OnlyBest_(iConfig.getParameter<bool>("OnlyBest")),
  isMC_(iConfig.getParameter<bool>("isMC")),
  OnlyGen_(iConfig.getParameter<bool>("OnlyGen")),
  isRes_(iConfig.getParameter<bool>("isRes")),
  mumuMassConstraint_(iConfig.getParameter<bool>("mumuMassConstraint")),
  mumuMasscut_(iConfig.getParameter<std::vector<double> >("mumuMasscut")),
  Trkmass_(iConfig.getParameter<double>("Trkmass")),
  Trkpdgid_(iConfig.getParameter<int>("Trkpdgid")),
  Bpdgid_(iConfig.getParameter<int>("Bpdgid")),
  Respdgid_(iConfig.getParameter<int>("Respdgid")),
  BarebMasscut_(iConfig.getParameter<std::vector<double> >("BarebMasscut")),
  bMasscut_(iConfig.getParameter<std::vector<double> >("bMasscut")),
  
  
  tree_(0), 

  mumC2(0), mumNHits(0), mumNPHits(0),
  mupC2(0), mupNHits(0), mupNPHits(0),
  mumdxy(0), mupdxy(0), mumdz(0), mupdz(0),
  muon_dca(0),

  tri_Dim25(0), tri_JpsiTrk_Bc(0), tri_JpsiTk(0),
  tri_DMu4_3_LM(0), tri_DMu4_LM_Displaced(0),
  
  mu1soft(0), mu2soft(0), mu1tight(0), mu2tight(0), 
  mu1PF(0), mu2PF(0), mu1loose(0), mu2loose(0),
 
  // *******************************************************
 
  nB(0), nMu(0),
  B_mass(0), B_px(0), B_py(0), B_pz(0), B_charge(0),
  B_k_px(0), B_k_py(0), B_k_pz(0), B_k_charge1(0),
  B_k_px_track(0), B_k_py_track(0), B_k_pz_track(0),
  
  B_J_mass(0), B_J_massErr(0), B_J_px(0), B_J_py(0), B_J_pz(0),
  B_J_pt1(0), B_J_px1(0), B_J_py1(0), B_J_pz1(0),
  B_J_pt2(0), B_J_px2(0), B_J_py2(0), B_J_pz2(0), 
  B_J_charge1(0), B_J_charge2(0),

  // Primary Vertex (PV)
  nVtx(0),
  priVtxX(0), priVtxY(0), priVtxZ(0), priVtxXE(0), priVtxYE(0), priVtxZE(0), priVtxCL(0),
  priVtxXYE(0), priVtxXZE(0), priVtxYZE(0),
  
  // ************************ ****************************************************

  B_chi2(0), B_J_chi2(0),
  B_Prob(0), B_J_Prob(0), 
 
  B_DecayVtxX(0),     B_DecayVtxY(0),     B_DecayVtxZ(0),
  B_DecayVtxXE(0),    B_DecayVtxYE(0),    B_DecayVtxZE(0),
  B_DecayVtxXYE(0),   B_DecayVtxXZE(0),   B_DecayVtxYZE(0),
 
  run(0), event(0),
  lumiblock(0)

{
   //now do what ever initialization is needed
}

//JPsiTrk::~JPsiTrk(){}


// ------------ method called to for each event  ------------
void JPsiTrk::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using std::vector;
  using namespace edm;
  using namespace reco;
  using namespace std;

  //*********************************
  // Get gen level information
  //*********************************

  edm::Handle<reco::GenParticleCollection> pruned;
  //edm::Handle<pat::PackedGenParticle> pruned; 
  iEvent.getByToken(genCands_, pruned);
  
  edm::Handle<pat::PackedGenParticleCollection> packed;
  iEvent.getByToken(packedGenToken_,packed);
  
  gen_bc_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_jpsi_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_pion3_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_muon1_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_muon2_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_bc_vtx.SetXYZ(0.,0.,0.);
  gen_jpsi_vtx.SetXYZ(0.,0.,0.);
  gen_bc_ct = -9999.;


if ( (isMC_ || OnlyGen_) && pruned.isValid() && isRes_) {
    int foundit = 0;
    for (size_t i=0; i<pruned->size(); i++) {
      foundit = 0;
      const reco::Candidate *dau = &(*pruned)[i];
      if ( (abs(dau->pdgId()) == Bpdgid_) ) { //&& (dau->status() == 2) ) {
		foundit++;
		gen_bc_p4.SetPtEtaPhiM(dau->pt(),dau->eta(),dau->phi(),dau->mass());
		gen_bc_vtx.SetXYZ(dau->vx(),dau->vy(),dau->vz());
		for (size_t k=0; k<dau->numberOfDaughters(); k++) {
	  		const reco::Candidate *gdau = dau->daughter(k);
	  		if (gdau->pdgId()== Respdgid_ ) { //&& gdau->status()==2) {  443-JPsi  100443 - Psi(2S)
	    		foundit++;
	    		gen_jpsi_vtx.SetXYZ(gdau->vx(),gdau->vy(),gdau->vz());
	    		gen_bc_ct = GetLifetime(gen_bc_p4,gen_bc_vtx,gen_jpsi_vtx);
	    		int nm=0;
	    		for (size_t l=0; l<gdau->numberOfDaughters(); l++) {
	      			const reco::Candidate *mm = gdau->daughter(l);
	      			if (mm->pdgId()==13) { foundit++;
						if (mm->status()!=1) {
		  					for (size_t m=0; m<mm->numberOfDaughters(); m++) {
		    					const reco::Candidate *mu = mm->daughter(m);
		    					if (mu->pdgId()==13 ) { //&& mu->status()==1) {
		      						nm++;
		      						gen_muon1_p4.SetPtEtaPhiM(mu->pt(),mu->eta(),mu->phi(),mu->mass());
		      						break;
		    					}
		  					}
						} else {
		  					gen_muon1_p4.SetPtEtaPhiM(mm->pt(),mm->eta(),mm->phi(),mm->mass());
		  					nm++;
						}
	      			}
	      			if (mm->pdgId()==-13) { foundit++;
						if (mm->status()!=1) {
		  					for (size_t m=0; m<mm->numberOfDaughters(); m++) {
		    					const reco::Candidate *mu = mm->daughter(m);
		    					if (mu->pdgId()==-13 ) { //&& mu->status()==1) {
		      						nm++;
		      						gen_muon2_p4.SetPtEtaPhiM(mu->pt(),mu->eta(),mu->phi(),mu->mass());
		      						break;
		    					}
		  					}
						} else {
		  					gen_muon2_p4.SetPtEtaPhiM(mm->pt(),mm->eta(),mm->phi(),mm->mass());
		  					nm++;
						}
	      			}
	    		}
	    		if (nm==2) gen_jpsi_p4.SetPtEtaPhiM(gdau->pt(),gdau->eta(),gdau->phi(),gdau->mass());
	   	 		else foundit-=nm;
	  		}
		} // for (size_t k

		for (size_t lk=0; lk<packed->size(); lk++) {
	  		const reco::Candidate * dauInPrunedColl = (*packed)[lk].mother(0);
	  		int stable_id = (*packed)[lk].pdgId();
	  
	  		if (dauInPrunedColl != nullptr && isAncestor(dau,dauInPrunedColl)) {
	    		if( abs(stable_id) == Trkpdgid_ ) {foundit++;
	      		gen_pion3_p4.SetPtEtaPhiM((*packed)[lk].pt(),(*packed)[lk].eta(),(*packed)[lk].phi(),(*packed)[lk].mass());
	    		}
	  		}
		}

    }   // if (abs(dau->pdgId())==521 )
      if (foundit>=5) break;
    } // for i
    if (foundit!=5) {
      gen_bc_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      gen_jpsi_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      gen_bc_vtx.SetXYZ(0.,0.,0.);
      gen_jpsi_vtx.SetXYZ(0.,0.,0.);
	  gen_pion3_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      gen_bc_ct = -9999.;
      //std::cout << "Does not found the given decay " << run << "," << event << " foundit=" << foundit << std::endl; // sanity check
    }
  }



//Non resonant MC

  if ( (isMC_ || OnlyGen_) && pruned.isValid() && !isRes_ ) {
    int foundit = 0;
    for (size_t i=0; i<pruned->size(); i++) {
		foundit = 0;
     	const reco::Candidate *dau = &(*pruned)[i];
      	if ( (abs(dau->pdgId()) == Bpdgid_) ) { 
			foundit++;
			gen_bc_p4.SetPtEtaPhiM(dau->pt(),dau->eta(),dau->phi(),dau->mass());
			gen_bc_vtx.SetXYZ(dau->vx(),dau->vy(),dau->vz());
			int nm=0;
			for (size_t k=0; k<dau->numberOfDaughters(); k++) {
	  			const reco::Candidate *gdau = dau->daughter(k);
				//:::Finding a positive muon and quiting all the resonant particles
	  			if (gdau->pdgId()==13 && !isAncestor(443,gdau) && !isAncestor(100443, gdau) ) { 
	    			foundit++;
	    			gen_jpsi_vtx.SetXYZ(gdau->vx(),gdau->vy(),gdau->vz());
	    			gen_bc_ct = GetLifetime(gen_bc_p4, gen_bc_vtx, gen_jpsi_vtx);
					if (gdau->status()!=1){
						for(size_t m=0; m<gdau->numberOfDaughters();m++){
							const reco::Candidate *mu=dau->daughter(m);
							//---Verify that this daughter is indeed a muon
							if(mu->pdgId()==13){
								nm++;//---Add a muon to the found muons list
								//---Store the muon p4
								gen_muon1_p4.SetPtEtaPhiM(mu->pt(),mu->eta(),mu->phi(),mu->mass());
								break;
							}//---End of if(muon's daughter is a muon)

						}//---End for muon's daughters
					}
					else{
						nm++;//---Add a muon to the found muons list
						//---Store the muon p4
						gen_muon1_p4.SetPtEtaPhiM(gdau->pt(),gdau->eta(),gdau->phi(),gdau->mass());
					}
				}
				
				
				//Now a negative muon
				if (gdau->pdgId()==-13 && !isAncestor(443, gdau) && !isAncestor(100443, gdau)){
					//---Here we have two cases, the dau-muon we found is stable and so we can store its information
					//or it is unstable and decays as mu->mu+photon, in which case we need to check the dau-muon daughters.
					//if dau->status()==1, the muon is stable, otherwise it is not.
					foundit++;//--We found a negative muon
					//std::cout<<"muon2 has been found "<<i<<" foundit="<<foundit<<std::endl;
					if (gdau->status()!=1){
						//---For in muon's daughters
						for(size_t m=0; m<dau->numberOfDaughters();m++){
							const reco::Candidate *mu=dau->daughter(m);
							//---Verify that this daughter is indeed a muon
							if(mu->pdgId()==-13){
								nm++;//---Add a muon to the found muons list
								//---Store the muon p4
								gen_muon2_p4.SetPtEtaPhiM(mu->pt(),mu->eta(),mu->phi(),mu->mass());
								break;
							}//---End of if(muon's daughter is a muon)

						}//---End for muon's daughters

					}//---End if(muon unstable)
					else{
						nm++;//---Add a muon to the found muons list
						//---Store the muon p4
						gen_muon2_p4.SetPtEtaPhiM(gdau->pt(),gdau->eta(),gdau->phi(),gdau->mass());	                                
					}

				}//---End of if(B's daughter is -muon)

				/*if (abs(gdau->pdgId())==Trackpdgid_){
					foundit++;
					gen_pion3_p4.SetPtEtaPhiM(dau->pt(), dau->eta(), dau->phi(), dau->mass());
					gen_pion3_vtx.SetXYZ(dau->vx(), dau->vy(), dau->vz());
				}*/


			}

			for (size_t lk=0; lk<packed->size(); lk++) {
	  			const reco::Candidate * dauInPrunedColl = (*packed)[lk].mother(0);
	  			int stable_id = (*packed)[lk].pdgId();
	  
	  			if (dauInPrunedColl != nullptr && isAncestor(dau,dauInPrunedColl)) {
	    			if( abs(stable_id) == Trkpdgid_ ) {
						foundit++;
	      				gen_pion3_p4.SetPtEtaPhiM((*packed)[lk].pt(),(*packed)[lk].eta(),(*packed)[lk].phi(),(*packed)[lk].mass());
	    			}
	  			}
			}	

			if (nm==2) gen_jpsi_p4 = gen_muon1_p4 + gen_muon2_p4;
	    	else foundit-=nm;	
	    			
	  	}

		if (foundit >= 4){
			break; //1.-B+-, 2.-mu1, 3.-mu2, 4.-K+
		}
	}


	if(foundit!=4){
		gen_bc_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      	gen_jpsi_p4.SetPtEtaPhiM(0.,0.,0.,0.);
		gen_pion3_p4.SetPtEtaPhiM(0.,0.,0.,0.);
		gen_muon1_p4.SetPtEtaPhiM(0.,0.,0.,0.);
		gen_muon2_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      	gen_bc_vtx.SetXYZ(0.,0.,0.);
      	gen_jpsi_vtx.SetXYZ(0.,0.,0.);
      	gen_bc_ct = -9999.;
	}

  }//End of non resonant MC



	if ( OnlyGen_ ) { 
    	tree_->Fill();
    	return;
  	}

  //*********************************
  // Get event content information
  //*********************************  
 
  // Kinematic fit
  edm::ESHandle<TransientTrackBuilder> theB = iSetup.getHandle(ttrkToken_);
  //edm::ESHandle<TransientTrackBuilder> theB; 
  //iSetup.getHandle<TransientTrackRecord>().get(ttrkToken_,theB);//
  //iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);// old one 
  //auto const& theB = iSetup.getHandle(transientTrackRecordToken_);
  //const edm::ESHandle<TransientTrackBuilder> TTbuilder = setup.getHandle(ttrkToken_);

  edm::Handle< View<pat::PackedCandidate> > thePATTrackHandle;
  iEvent.getByToken(trakCollection_label,thePATTrackHandle);

  edm::Handle< View<pat::Muon> > thePATMuonHandle;
  iEvent.getByToken(dimuon_Label,thePATMuonHandle);
 
  //*********************************
  //Now we get the primary vertex 
  //*********************************

  reco::Vertex bestVtx;
  edm::Handle<reco::VertexCollection> primaryVertices_handle;
  iEvent.getByToken(primaryVertices_Label, primaryVertices_handle);

  // get primary vertex
  bestVtx = *(primaryVertices_handle->begin());

  priVtxX = bestVtx.x();
  priVtxY = bestVtx.y();
  priVtxZ = bestVtx.z();
  priVtxXE = bestVtx.covariance(0, 0);
  priVtxYE = bestVtx.covariance(1, 1);
  priVtxZE = bestVtx.covariance(2, 2);
  priVtxXYE = bestVtx.covariance(0, 1);
  priVtxXZE = bestVtx.covariance(0, 2);
  priVtxYZE = bestVtx.covariance(1, 2);

  priVtxCL = ChiSquaredProbability((double)(bestVtx.chi2()),(double)(bestVtx.ndof())); 
  nVtx = primaryVertices_handle->size(); 
 
  lumiblock = iEvent.id().luminosityBlock();
  run = iEvent.id().run();
  event = iEvent.id().event();

  //*****************************************
  //Let's begin by looking for J/psi+K^+

  unsigned int nMu_tmp = thePATMuonHandle->size();
  //nMu = nMu_tmp;

 for(View<pat::Muon>::const_iterator iMuon1 = thePATMuonHandle->begin(); iMuon1 != thePATMuonHandle->end(); ++iMuon1) 
    {     
      for(View<pat::Muon>::const_iterator iMuon2 = iMuon1+1; iMuon2 != thePATMuonHandle->end(); ++iMuon2) 
	{
	  if(iMuon1==iMuon2) continue;
	  
	  //opposite charge 
	  if( (iMuon1->charge())*(iMuon2->charge()) == 1) continue;

	  TrackRef glbTrackP;	  
	  TrackRef glbTrackM;	  
	  
	  if(iMuon1->charge() == 1){ glbTrackP = iMuon1->track();}
	  if(iMuon1->charge() == -1){ glbTrackM = iMuon1->track();}
	  
	  if(iMuon2->charge() == 1) { glbTrackP = iMuon2->track();}
	  if(iMuon2->charge() == -1){ glbTrackM = iMuon2->track();}
	  
	  if( glbTrackP.isNull() || glbTrackM.isNull() ) 
	    {
	      //std::cout << "continue due to no track ref" << endl;
	      continue;
	    }

	  if(iMuon1->track()->pt()<4.0) continue;
	  if(iMuon2->track()->pt()<4.0) continue;

	  if(!(glbTrackM->quality(reco::TrackBase::highPurity))) continue;
	  if(!(glbTrackP->quality(reco::TrackBase::highPurity))) continue;	 
	  
	  reco::TransientTrack muon1TT((*theB).build(glbTrackP));
	  reco::TransientTrack muon2TT((*theB).build(glbTrackM));

	 // *****  Trajectory states to calculate DCA for the 2 muons *********************
	  FreeTrajectoryState mu1State = muon1TT.impactPointTSCP().theState();
	  FreeTrajectoryState mu2State = muon2TT.impactPointTSCP().theState();

	  if( !muon1TT.impactPointTSCP().isValid() || !muon2TT.impactPointTSCP().isValid() ) continue;

	  // Measure distance between tracks at their closest approach
	  ClosestApproachInRPhi cApp;
	  cApp.calculate(mu1State, mu2State);
	  if( !cApp.status() ) continue;
	  float dca = fabs( cApp.distance() );	  
	  //if (dca < 0. || dca > 0.5) continue;
	  //cout<<" closest approach  "<<dca<<endl;

	  // *****  end DCA for the 2 muons *********************

	  //Let's check the vertex and mass

	  //The mass of a muon and the insignificant mass sigma 
	  //to avoid singularities in the covariance matrix.
	  ParticleMass muon_mass = 0.10565837; //pdg mass
	  float muon_sigma = muon_mass*1.e-6;
	  
	  //Creating a KinematicParticleFactory
	  KinematicParticleFactoryFromTransientTrack pFactory;
	  
	  //initial chi2 and ndf before kinematic fits.
	  float chi = 0.;
	  float ndf = 0.;
	  vector<RefCountedKinematicParticle> muonParticles;
	  try {
	    muonParticles.push_back(pFactory.particle(muon1TT,muon_mass,chi,ndf,muon_sigma));
	    muonParticles.push_back(pFactory.particle(muon2TT,muon_mass,chi,ndf,muon_sigma));
	  }
	  catch(...) { 
	    std::cout<<" Exception caught ... continuing 1 "<<std::endl; 
	    continue;
	  }
	  
	  KinematicParticleVertexFitter fitter;   
	  
	  RefCountedKinematicTree psiVertexFitTree;
	  try {
	    psiVertexFitTree = fitter.fit(muonParticles); 
	  }
	  catch (...) { 
	    std::cout<<" Exception caught ... continuing 2 "<<std::endl; 
	    continue;
	  }
	  
	  if (!psiVertexFitTree->isValid()) 
	    {
	      //std::cout << "caught an exception in the psi vertex fit" << std::endl;
	      continue; 
	    }
	  
	  psiVertexFitTree->movePointerToTheTop();
	  
	  RefCountedKinematicParticle psi_vFit_noMC = psiVertexFitTree->currentParticle();
	  RefCountedKinematicVertex psi_vFit_vertex_noMC = psiVertexFitTree->currentDecayVertex();
	  
	  if( psi_vFit_vertex_noMC->chiSquared() < 0 )
	    {
	      //std::cout << "negative chisq from psi fit" << endl;
	      continue;
	    }
	  
	  //some loose cuts go here	  
	  if(psi_vFit_vertex_noMC->chiSquared()>50.) continue;
	  // ******************************************** JPsi mass  ********************************************
	  //if(psi_vFit_noMC->currentState().mass()<2.9 || psi_vFit_noMC->currentState().mass()>3.3) continue;
	  // ******************************************** mumu mass  ********************************************
	  //if(psi_vFit_noMC->currentState().mass()<0.2 || psi_vFit_noMC->currentState().mass()>4.9) continue;
	  
	  if(psi_vFit_noMC->currentState().mass()<mumuMasscut_[0] || psi_vFit_noMC->currentState().mass()>mumuMasscut_[1]) continue;
	  

	  double J_Prob_tmp   = TMath::Prob(psi_vFit_vertex_noMC->chiSquared(),(int)psi_vFit_vertex_noMC->degreesOfFreedom());
	  if(J_Prob_tmp<0.01)
	    {
	      continue;
	    }
	  
	  //Now that we have a J/psi candidate, we look for Trk candidates
	  
	  for(View<pat::PackedCandidate>::const_iterator iTrack1 = thePATTrackHandle->begin(); 
		   iTrack1 != thePATTrackHandle->end(); ++iTrack1 ) 
		   {
		     //quality cuts
		   if(iTrack1->charge()==0) continue;
		   if(fabs(iTrack1->pdgId())!=211) continue;
		   //if(iTrack1->pt()<1.2) continue;
		   if(iTrack1->pt()<0.95) continue;
		   if(!(iTrack1->trackHighPurity())) continue;
		   if(iTrack1->numberOfPixelHits()<1)continue;
		   if(iTrack1->numberOfHits()<5)continue;
		   
		   if ( IsTheSame(*iTrack1,*iMuon1) || IsTheSame(*iTrack1,*iMuon2) ) continue;
		    		 		   
		   reco::TransientTrack kaonTT((*theB).build(iTrack1->pseudoTrack()));

		   //ParticleMass kaon_mass = 0.493677;
		   //float kaon_sigma = kaon_mass*1.e-6;
		   //ParticleMass pion_mass = 0.13957018;
		   ParticleMass pion_mass = Trkmass_;
		   float pion_sigma = pion_mass*1.e-6;

		   float chi = 0.;
		   float ndf = 0.;		 

		   // ***************************
		   // JpsiTrk invariant mass (before kinematic vertex fit)
		   // ***************************
		   TLorentzVector kaon14V, Jpsi4V; 
		   kaon14V.SetXYZM(iTrack1->px(),iTrack1->py(),iTrack1->pz(),pion_mass);		   
		   Jpsi4V.SetXYZM(psi_vFit_noMC->currentState().globalMomentum().x(),psi_vFit_noMC->currentState().globalMomentum().y(),psi_vFit_noMC->currentState().globalMomentum().z(),psi_vFit_noMC->currentState().mass());
		   //if ( (kaon14V + Jpsi4V).M()<4.2 || (kaon14V + Jpsi4V).M()>6.8 ) continue;
		   if ( (kaon14V + Jpsi4V).M()<BarebMasscut_[0] || (kaon14V + Jpsi4V).M()>BarebMasscut_[1] ) continue;
	   
		   //Now we are ready to combine!
		   // JPsi mass constraint is applied in the final Bplus fit,
		   vector<RefCountedKinematicParticle> vFitMCParticles;
		   vFitMCParticles.push_back(pFactory.particle(muon1TT,muon_mass,chi,ndf,muon_sigma));
		   vFitMCParticles.push_back(pFactory.particle(muon2TT,muon_mass,chi,ndf,muon_sigma));
		   vFitMCParticles.push_back(pFactory.particle(kaonTT,pion_mass ,chi,ndf,pion_sigma));

		   //bool mumuMassConstraint_ = true;
		   RefCountedKinematicTree vertexFitTree;

		   if(mumuMassConstraint_){
		     // **************  JPsi mass constraint is applied in the final B fit ******************************************************
		     ParticleMass psi_mass = 3.096916;
		     //float psi_sigma = psi_mass*1.e-6; 
		     MultiTrackKinematicConstraint *  j_psi_c = new  TwoTrackMassKinematicConstraint(psi_mass);
		     KinematicConstrainedVertexFitter kcvFitter;
		     vertexFitTree = kcvFitter.fit(vFitMCParticles, j_psi_c);
		     if (!vertexFitTree->isValid()) {
		       //std::cout << "caught an exception in the B vertex fit with MC" << std::endl;
		       continue;
		     }
		   }
		   else
		     {
		       // **************  mumu channel do not have mass constraint in the final Bs fit **********************************************                                                           
		       KinematicParticleVertexFitter kcvFitter;
		       //RefCountedKinematicTree vertexFitTree = kcvFitter.fit(vFitMCParticles);
		       vertexFitTree = kcvFitter.fit(vFitMCParticles);
		       if (!vertexFitTree->isValid()) {
			 //std::cout << "caught an exception in the B vertex fit with MC" << std::endl;                                        
			 continue;
		       }
		     }

		   
		   vertexFitTree->movePointerToTheTop();
		 
		   RefCountedKinematicParticle bCandMC = vertexFitTree->currentParticle();
		   RefCountedKinematicVertex bDecayVertexMC = vertexFitTree->currentDecayVertex();
		   if (!bDecayVertexMC->vertexIsValid()){
		     // cout << "B MC fit vertex is not valid" << endl;
		     continue;
		   }
		   
		   //if( (bCandMC->currentState().mass() < 5.0) || (bCandMC->currentState().mass() > 6.0) ) { continue; }
		   if(bCandMC->currentState().mass()<bMasscut_[0] || bCandMC->currentState().mass()>bMasscut_[1]) continue;
		   
		   double B_Prob_tmp  = TMath::Prob(bDecayVertexMC->chiSquared(),(int)bDecayVertexMC->degreesOfFreedom());
		   if(B_Prob_tmp<0.01)
		     {
		       continue;
		     }
		   
		    // get children from final B fit

		    vertexFitTree->movePointerToTheFirstChild();
		    RefCountedKinematicParticle mu1CandMC = vertexFitTree->currentParticle();
		    
		    vertexFitTree->movePointerToTheNextChild();
		    RefCountedKinematicParticle mu2CandMC = vertexFitTree->currentParticle();
		   

		   vertexFitTree->movePointerToTheNextChild();
		   RefCountedKinematicParticle kCandMC = vertexFitTree->currentParticle();		  

		   KinematicParameters psiMu1KP = mu1CandMC->currentState().kinematicParameters();
		   KinematicParameters psiMu2KP = mu2CandMC->currentState().kinematicParameters();
		   KinematicParameters psiMupKP;
		   KinematicParameters psiMumKP;
	       
		   if ( mu1CandMC->currentState().particleCharge() > 0 ) psiMupKP = psiMu1KP;
		   if ( mu1CandMC->currentState().particleCharge() < 0 ) psiMumKP = psiMu1KP;
		   if ( mu2CandMC->currentState().particleCharge() > 0 ) psiMupKP = psiMu2KP;
		   if ( mu2CandMC->currentState().particleCharge() < 0 ) psiMumKP = psiMu2KP;	 

 		   GlobalVector Jp1vec(mu1CandMC->currentState().globalMomentum().x(),
				       mu1CandMC->currentState().globalMomentum().y(),
 				       mu1CandMC->currentState().globalMomentum().z());


 		   GlobalVector Jp2vec(mu2CandMC->currentState().globalMomentum().x(),
				       mu2CandMC->currentState().globalMomentum().y(),
 				       mu2CandMC->currentState().globalMomentum().z());
		   
		   KinematicParameters VCandKP = kCandMC->currentState().kinematicParameters();
		   	       
		   // ************ fill candidate variables now
		   
		   // Only save the first time
		   if(nB==0){	    
		     nMu  = nMu_tmp;
		     // cout<< "*Number of Muons : " << nMu_tmp << endl;
		   } // end nB==0

		   B_mass = bCandMC->currentState().mass();
		   B_px = bCandMC->currentState().globalMomentum().x();
		   B_py = bCandMC->currentState().globalMomentum().y();
		   B_pz = bCandMC->currentState().globalMomentum().z();
		   B_charge = bCandMC->currentState().particleCharge();

		   // You can get the momentum components (for muons and kaon) from the final B childrens or of the original Tracks. Here, a example for the kaon:
		   B_k_px = VCandKP.momentum().x() ;
		   B_k_py = VCandKP.momentum().y() ;
		   B_k_pz = VCandKP.momentum().z() ;
		   B_k_px_track = iTrack1->px() ;
		   B_k_py_track = iTrack1->py() ;
		   B_k_pz_track = iTrack1->pz() ;
		   B_k_charge1 = kCandMC->currentState().particleCharge();
		  
		   B_J_mass =  psi_vFit_noMC->currentState().mass() ;
		   B_J_massErr = sqrt(psi_vFit_noMC->currentState().kinematicParametersError().matrix()(6,6));
		   B_J_px =  psi_vFit_noMC->currentState().globalMomentum().x() ;
		   B_J_py =  psi_vFit_noMC->currentState().globalMomentum().y() ;
		   B_J_pz =  psi_vFit_noMC->currentState().globalMomentum().z() ;

		   B_J_pt1 = Jp1vec.perp();
		   B_J_px1 = psiMu1KP.momentum().x();
		   B_J_py1 = psiMu1KP.momentum().y();
		   B_J_pz1 = psiMu1KP.momentum().z();
		   B_J_charge1 = mu1CandMC->currentState().particleCharge();

		   B_J_pt2 = Jp2vec.perp();
		   B_J_px2 = psiMu2KP.momentum().x();
		   B_J_py2 = psiMu2KP.momentum().y();
		   B_J_pz2 = psiMu2KP.momentum().z();
		   B_J_charge2 = mu2CandMC->currentState().particleCharge();
		  
		   B_J_chi2 = psi_vFit_vertex_noMC->chiSquared();
		   B_chi2 = bDecayVertexMC->chiSquared();
             
		   //double B_Prob_tmp       = TMath::Prob(bDecayVertexMC->chiSquared(),(int)bDecayVertexMC->degreesOfFreedom());
		   //double J_Prob_tmp   = TMath::Prob(psi_vFit_vertex_noMC->chiSquared(),(int)psi_vFit_vertex_noMC->degreesOfFreedom());
		   B_Prob     = B_Prob_tmp;
		   B_J_Prob   = J_Prob_tmp;

		   B_DecayVtxX  = (*bDecayVertexMC).position().x();    
		   B_DecayVtxY  = (*bDecayVertexMC).position().y();
		   B_DecayVtxZ  = (*bDecayVertexMC).position().z();
		   B_DecayVtxXE  = bDecayVertexMC->error().cxx();   
		   B_DecayVtxYE  = bDecayVertexMC->error().cyy();   
		   B_DecayVtxZE  = bDecayVertexMC->error().czz();
		   B_DecayVtxXYE  = bDecayVertexMC->error().cyx();
		   B_DecayVtxXZE  = bDecayVertexMC->error().czx();
		   B_DecayVtxYZE  = bDecayVertexMC->error().czy();
		 

 // ********************* muon-trigger-machint ****************

		   const pat::Muon* muon1 = &(*iMuon1);
		   const pat::Muon* muon2 = &(*iMuon2);
		   
		   int tri_Dim25_tmp = 0, tri_JpsiTk_tmp = 0,  tri_JpsiTrk_Bc_tmp = 0, tri_DMu4_3_LM_tmp = 0, tri_DMu4_LM_Displaced_tmp = 0; 
		   
		   if(muon1->triggerObjectMatchByPath("HLT_Dimuon25_Jpsi_v*")!=nullptr && muon2->triggerObjectMatchByPath("HLT_Dimuon25_Jpsi_v*")!=nullptr) tri_Dim25_tmp = 1;
		   if(muon1->triggerObjectMatchByPath("HLT_DoubleMu4_MuMuTrk_Displaced_v*")!=nullptr && muon2->triggerObjectMatchByPath("HLT_DoubleMu4_MuMuTrk_Displaced_v*")!=nullptr) tri_JpsiTk_tmp = 1;
		   if(muon1->triggerObjectMatchByPath("HLT_DoubleMu4_JpsiTrk_Bc_v*")!=nullptr && muon2->triggerObjectMatchByPath("HLT_DoubleMu4_JpsiTrk_Bc_v*")!=nullptr) tri_JpsiTrk_Bc_tmp = 1;
		   if(muon1->triggerObjectMatchByPath("HLT_DoubleMu4_3_LowMass_v*")!=nullptr && muon2->triggerObjectMatchByPath("HLT_DoubleMu4_3_LowMass_v*")!=nullptr) tri_DMu4_3_LM_tmp = 1;
		   if(muon1->triggerObjectMatchByPath("HLT_DoubleMu4_LowMass_Displaced_v*")!=nullptr && muon2->triggerObjectMatchByPath("HLT_DoubleMu4_LowMass_Displaced_v*")!=nullptr) tri_DMu4_LM_Displaced_tmp = 1;

		   
		   tri_Dim25 =  tri_Dim25_tmp ;	       
		   tri_JpsiTk =  tri_JpsiTk_tmp ;
		   tri_JpsiTrk_Bc =  tri_JpsiTrk_Bc_tmp ;
		   tri_DMu4_3_LM =  tri_DMu4_3_LM_tmp ;
		   tri_DMu4_LM_Displaced =  tri_DMu4_LM_Displaced_tmp ;
		  
		   // ************ Different muons Id, and other properties  ****************

		   mu1soft = iMuon1->isSoftMuon(bestVtx) ;
		   mu2soft = iMuon2->isSoftMuon(bestVtx) ;
		   mu1tight = iMuon1->isTightMuon(bestVtx) ;
		   mu2tight = iMuon2->isTightMuon(bestVtx) ;
		   mu1PF = iMuon1->isPFMuon();
		   mu2PF = iMuon2->isPFMuon();
		   mu1loose = muon::isLooseMuon(*iMuon1);
		   mu2loose = muon::isLooseMuon(*iMuon2);

		   mumC2 =  glbTrackP->normalizedChi2() ;
		   mumNHits =  glbTrackP->numberOfValidHits() ;
		   mumNPHits =  glbTrackP->hitPattern().numberOfValidPixelHits() ;	       
		   mupC2 =  glbTrackM->normalizedChi2() ;
		   mupNHits =  glbTrackM->numberOfValidHits() ;
		   mupNPHits =  glbTrackM->hitPattern().numberOfValidPixelHits() ;
                   mumdxy = glbTrackP->dxy(bestVtx.position()) ;// 
		   mupdxy = glbTrackM->dxy(bestVtx.position()) ;// 
		   mumdz = glbTrackP->dz(bestVtx.position()) ;
		   mupdz = glbTrackM->dz(bestVtx.position()) ;
		   muon_dca = dca;

		   //fill the tree
		   tree_->Fill();

		   nB++;	       
		   muonParticles.clear();
		   vFitMCParticles.clear();

	    }
	  }
      	}
 
 /*
  if (nB > 0 ) 
    {

      //std::cout << "filling tree" << endl;
      tree_->Fill();
    }
 */

   nB = 0; nMu = 0;

   B_charge = 0;
   B_mass = 0;    B_px = 0;    B_py = 0;    B_pz = 0; 
   B_k_px = 0; B_k_py = 0; B_k_pz = 0;  B_k_charge1 = 0;
   B_k_px_track = 0; B_k_py_track = 0; B_k_pz_track = 0;

   B_J_mass = 0;  B_J_massErr = 0;  B_J_px = 0;  B_J_py = 0;  B_J_pz = 0;
 
   B_J_pt1 = 0;  B_J_px1 = 0;  B_J_py1 = 0;  B_J_pz1 = 0; B_J_charge1 = 0;
   B_J_pt2 = 0;  B_J_px2 = 0;  B_J_py2 = 0;  B_J_pz2 = 0; B_J_charge2 = 0;

   B_chi2 = 0; B_J_chi2 = 0; 
   B_Prob = 0; B_J_Prob = 0;

   B_DecayVtxX = 0;     B_DecayVtxY = 0;     B_DecayVtxZ = 0;
   B_DecayVtxXE = 0;    B_DecayVtxYE = 0;    B_DecayVtxZE = 0;
   B_DecayVtxXYE = 0;   B_DecayVtxXZE = 0;   B_DecayVtxYZE = 0;

   nVtx = 0;
   priVtxX = 0;     priVtxY = 0;     priVtxZ = 0; 
   priVtxXE = 0;    priVtxYE = 0;    priVtxZE = 0; priVtxCL = 0;
   priVtxXYE = 0;   priVtxXZE = 0;   priVtxYZE = 0;    

   mumC2 = 0;
   mumNHits = 0; mumNPHits = 0;
   mupC2 = 0;
   mupNHits = 0; mupNPHits = 0;
   mumdxy = 0; mupdxy = 0; mumdz = 0; mupdz = 0; muon_dca = 0;

   tri_Dim25 = 0; tri_JpsiTrk_Bc = 0; tri_JpsiTk = 0;
   tri_DMu4_3_LM = 0;  tri_DMu4_LM_Displaced = 0;
 
   mu1soft = 0; mu2soft = 0; mu1tight = 0; mu2tight = 0;
   mu1PF = 0; mu2PF = 0; mu1loose = 0; mu2loose = 0; 
  
}

bool JPsiTrk::IsTheSame(const pat::GenericParticle& tk, const pat::Muon& mu){
  double DeltaEta = fabs(mu.eta()-tk.eta());
  double DeltaP   = fabs(mu.p()-tk.p());
  if (DeltaEta < 0.02 && DeltaP < 0.02) return true;
  return false;
}

bool JPsiTrk::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
    if (ancestor == particle ) return true;
    for (size_t i=0; i< particle->numberOfMothers(); i++) {
        if (isAncestor(ancestor,particle->mother(i))) return true;
    }
    return false;
}

bool JPsiTrk::isAncestor(int a_pdgId, const reco::Candidate * particle) {
    if (a_pdgId == particle->pdgId() ) return true;
    for (size_t i=0; i< particle->numberOfMothers(); i++) {
        if (isAncestor(a_pdgId,particle->mother(i))) return true;
    }
    return false;
}

double JPsiTrk::GetLifetime(TLorentzVector b_p4, TVector3 production_vtx, TVector3 decay_vtx) {
   TVector3 pv_dv = decay_vtx - production_vtx;
   TVector3 b_p3  = b_p4.Vect();
   pv_dv.SetZ(0.);
   b_p3.SetZ(0.);
   Double_t lxy   = pv_dv.Dot(b_p3)/b_p3.Mag();
   return lxy*b_p4.M()/b_p3.Mag();
}

// ------------ method called once each job just before starting event loop  ------------

void 
JPsiTrk::beginJob()
{

  std::cout << "Beginning analyzer job with value of isMC= " << isMC_ << std::endl;

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("ntuple","B+->J/psiK+ ntuple");

  tree_->Branch("nB",&nB,"nB/i");
  tree_->Branch("nMu",&nMu,"nMu/i");

  tree_->Branch("B_charge", &B_charge);
  tree_->Branch("B_mass", &B_mass);
  tree_->Branch("B_px", &B_px);
  tree_->Branch("B_py", &B_py);
  tree_->Branch("B_pz", &B_pz);

  tree_->Branch("B_k_charge1", &B_k_charge1);
  tree_->Branch("B_k_px", &B_k_px);
  tree_->Branch("B_k_py", &B_k_py);
  tree_->Branch("B_k_pz", &B_k_pz);
  tree_->Branch("B_k_px_track", &B_k_px_track);
  tree_->Branch("B_k_py_track", &B_k_py_track);
  tree_->Branch("B_k_pz_track", &B_k_pz_track);

  tree_->Branch("B_J_mass", &B_J_mass);
  tree_->Branch("B_J_massErr", &B_J_massErr);
  tree_->Branch("B_J_px", &B_J_px);
  tree_->Branch("B_J_py", &B_J_py);
  tree_->Branch("B_J_pz", &B_J_pz);

  tree_->Branch("B_J_pt1", &B_J_pt1);
  tree_->Branch("B_J_px1", &B_J_px1);
  tree_->Branch("B_J_py1", &B_J_py1);
  tree_->Branch("B_J_pz1", &B_J_pz1);
  tree_->Branch("B_J_charge1", &B_J_charge1);

  tree_->Branch("B_J_pt2", &B_J_pt2);
  tree_->Branch("B_J_px2", &B_J_px2);
  tree_->Branch("B_J_py2", &B_J_py2);
  tree_->Branch("B_J_pz2", &B_J_pz2);
  tree_->Branch("B_J_charge2", &B_J_charge2);

  tree_->Branch("B_chi2",    &B_chi2);
  tree_->Branch("B_J_chi2",  &B_J_chi2);

  tree_->Branch("B_Prob",    &B_Prob);
  tree_->Branch("B_J_Prob",  &B_J_Prob);
       
  tree_->Branch("B_DecayVtxX",     &B_DecayVtxX);
  tree_->Branch("B_DecayVtxY",     &B_DecayVtxY);
  tree_->Branch("B_DecayVtxZ",     &B_DecayVtxZ);
  tree_->Branch("B_DecayVtxXE",    &B_DecayVtxXE);
  tree_->Branch("B_DecayVtxYE",    &B_DecayVtxYE);
  tree_->Branch("B_DecayVtxZE",    &B_DecayVtxZE);
  tree_->Branch("B_DecayVtxXYE",   &B_DecayVtxXYE);
  tree_->Branch("B_DecayVtxXZE",   &B_DecayVtxXZE);
  tree_->Branch("B_DecayVtxYZE",   &B_DecayVtxYZE);

  tree_->Branch("priVtxX",&priVtxX, "priVtxX/D");
  tree_->Branch("priVtxY",&priVtxY, "priVtxY/D");
  tree_->Branch("priVtxZ",&priVtxZ, "priVtxZ/D");
  tree_->Branch("priVtxXE",&priVtxXE, "priVtxXE/D");
  tree_->Branch("priVtxYE",&priVtxYE, "priVtxYE/D");
  tree_->Branch("priVtxZE",&priVtxZE, "priVtxZE/D");
  tree_->Branch("priVtxXYE",&priVtxXYE, "priVtxXYE/D");
  tree_->Branch("priVtxXZE",&priVtxXZE, "priVtxXZE/D");
  tree_->Branch("priVtxYZE",&priVtxYZE, "priVtxYZE/D");
  tree_->Branch("priVtxCL",&priVtxCL, "priVtxCL/D");

  tree_->Branch("nVtx",       &nVtx);
  tree_->Branch("run",        &run,       "run/I");
  tree_->Branch("event",        &event,     "event/L");
  tree_->Branch("lumiblock",&lumiblock,"lumiblock/I");
    
  // *************************
 
  tree_->Branch("mumC2",&mumC2);  
  tree_->Branch("mumNHits",&mumNHits);
  tree_->Branch("mumNPHits",&mumNPHits);
  tree_->Branch("mupC2",&mupC2);  
  tree_->Branch("mupNHits",&mupNHits);
  tree_->Branch("mupNPHits",&mupNPHits);
  tree_->Branch("mumdxy",&mumdxy);
  tree_->Branch("mupdxy",&mupdxy);
  tree_->Branch("mumdz",&mumdz);
  tree_->Branch("mupdz",&mupdz);
  tree_->Branch("muon_dca",&muon_dca);

  tree_->Branch("tri_Dim25",&tri_Dim25);
  tree_->Branch("tri_JpsiTrk_Bc",&tri_JpsiTrk_Bc);
  tree_->Branch("tri_JpsiTk",&tri_JpsiTk);
  tree_->Branch("tri_DMu4_3_LM",&tri_DMu4_3_LM);
  tree_->Branch("tri_DMu4_LM_Displaced",&tri_DMu4_LM_Displaced);


  tree_->Branch("mu1soft",&mu1soft);
  tree_->Branch("mu2soft",&mu2soft);
  tree_->Branch("mu1tight",&mu1tight);
  tree_->Branch("mu2tight",&mu2tight);
  tree_->Branch("mu1PF",&mu1PF);
  tree_->Branch("mu2PF",&mu2PF);
  tree_->Branch("mu1loose",&mu1loose);
  tree_->Branch("mu2loose",&mu2loose);

  // gen
  if (isMC_) {
     tree_->Branch("gen_bc_p4",     "TLorentzVector",  &gen_bc_p4);
     tree_->Branch("gen_jpsi_p4",   "TLorentzVector",  &gen_jpsi_p4);
     tree_->Branch("gen_pion3_p4",  "TLorentzVector",  &gen_pion3_p4);
     tree_->Branch("gen_muon1_p4",  "TLorentzVector",  &gen_muon1_p4);
     tree_->Branch("gen_muon2_p4",  "TLorentzVector",  &gen_muon2_p4);
     tree_->Branch("gen_bc_vtx",    "TVector3",        &gen_bc_vtx);
     tree_->Branch("gen_jpsi_vtx",  "TVector3",        &gen_jpsi_vtx);
     tree_->Branch("gen_bc_ct",     &gen_bc_ct,        "gen_bc_ct/D");
  }

}




// ------------ method called once each job just after ending the event loop  ------------
void JPsiTrk::endJob() {
  tree_->GetDirectory()->cd();
  tree_->Write();
}

//define this as a plug-in
DEFINE_FWK_MODULE(JPsiTrk);
