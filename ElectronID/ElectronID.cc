#include "ElectronID.hh"

#include "edm4hep/utils/vector_utils.h"
#include "edm4eic/ClusterCollection.h"

#include <iostream>

#include <Math/LorentzVector.h>
using ROOT::Math::PxPyPzEVector;

ElectronID::ElectronID() {

	mEe = 10.;
	mEh = 100.;
	std::cout << "!!! ElectronID: You have not specified beam energies...defaulting to 10x100 GeV !!!" << std::endl;

	mEoP_min = 0.8;
	mEoP_max = 1.2;

	mDeltaH_min = 0.*mEe;
	mDeltaH_min = 5.*mEe;
		
	mIsoR = 0.4;
	mIsoE = 0.9;

	minTrackPoints = 3;

	boost = LorentzRotation(); // Initialize to identity
}

ElectronID::ElectronID(double Ee, double Eh) {

	mEe = Ee;
	mEh = Eh;

	mEoP_min = 0.8;
	mEoP_max = 1.2;

	mDeltaH_min = 0.*mEe;
	mDeltaH_min = 5.*mEe;
		
	mIsoR = 0.4;
	mIsoE = 0.9;

	minTrackPoints = 3;

	boost = LorentzRotation(); // Initialize to identity
}

ElectronID::~ElectronID() {
}

void ElectronID::SetEvent(const podio::Frame* event) {

	// std::cout << "** Setting event in ElectronID... " << std::endl;
	mEvent = event;
	eScatIndex = -1;
	hfs_dpt.clear();
	hfs_dpz.clear();
	hfs_de.clear();
	hfs_theta.clear();
	e_det.clear();
	jet_e_det.clear();
	pi_det.clear();
	else_det.clear();
	// std::cout << "** Event set in ElectronID. " << std::endl;
	return;
}

edm4hep::MCParticle ElectronID::GetMC(edm4eic::ReconstructedParticle e_rec) {

	// std::cout << "Available collections:" << std::endl;
	// for (const auto& name : mEvent->getAvailableCollections()) {
	// 	std::cout << "  " << name << std::endl;
	// }

	// std::cout << "Getting MC particle associated with reconstructed particle..." << std::endl;

	const auto& RecoMC = mEvent->get<edm4eic::MCRecoParticleAssociationCollection>("ReconstructedParticleAssociations");
	for(const auto& assoc : RecoMC) {
		if(assoc.getRec() == e_rec)
			return assoc.getSim();
	}

	return edm4hep::MCParticle();
}

int ElectronID::Check_eID(edm4eic::ReconstructedParticle e_rec) {

	edm4hep::MCParticleCollection meMC = GetMCElectron();
	if ( meMC.size() == 0 )
		return 86; // No MC electron found

	const auto& RecoMC = mEvent->get<edm4eic::MCRecoParticleAssociationCollection>("ReconstructedParticleAssociations");
	for(const auto& assoc : RecoMC) {
		if(assoc.getRec() == e_rec)
		{
		 	if (assoc.getSim() == meMC[0]) 
				return 0;
			else
				return assoc.getSim().getPDG();
		}
	}

	return 86;
}

void ElectronID::CheckClusters() {

	const auto& EcalEndcapNClusters = mEvent->get<edm4eic::ClusterCollection>("EcalEndcapNClusters");
	const auto& EcalBarrelScFiClusters = mEvent->get<edm4eic::ClusterCollection>("EcalBarrelScFiClusters");
	const auto& EcalEndcapPClusters = mEvent->get<edm4eic::ClusterCollection>("EcalEndcapPClusters");

	std::cout << " Number of clusters in EcalEndcapN: " << EcalEndcapNClusters.size() << std::endl;
	std::cout << " Number of clusters in EcalBarrelScFi: " << EcalBarrelScFiClusters.size() << std::endl;
	std::cout << " Number of clusters in EcalEndcapP: " << EcalEndcapPClusters.size() << std::endl;

	return;
}

edm4eic::ReconstructedParticleCollection ElectronID::FindHadronicFinalState(int object_id) {

	// edm4eic::HadronicFinalStateCollection meRecon;
	edm4eic::ReconstructedParticleCollection meRecon;
	meRecon.setSubsetCollection();

	// auto& rcparts = mEvent->get<edm4eic::HadronicFinalStateCollection>("HadronicFinalState");
	const auto& rcparts = mEvent->get<edm4eic::ReconstructedParticleCollection>("ReconstructedParticles");

	for(const auto& mcp : rcparts) {
		if ( mcp.getObjectID().index != object_id )
			meRecon.push_back(mcp);
	}

	return meRecon;
}

edm4eic::ReconstructedParticleCollection ElectronID::FindScatteredElectron() {

	// std::cout << "\nFinding scattered electron candidates..." << std::endl;
	// CheckClusters();

	// Get all the edm4eic objects needed for electron ID
	const auto& rcparts = mEvent->get<edm4eic::ReconstructedParticleCollection>("ReconstructedParticles");
	
	// Create collection for storing scattered electron candidates
	// (subset collection of ReconstructedParticleCollection)
	edm4eic::ReconstructedParticleCollection scatteredElectronCandidates;
	scatteredElectronCandidates.setSubsetCollection();

	// Loop over all ReconstructedParticles for this event
	for (const auto& reconPart : rcparts) {

		// Require at least one track and one cluster
		if(reconPart.getClusters().size() == 0 || reconPart.getTracks().size() == 0) continue;

		// Require negative particle
		if(reconPart.getCharge() >= 0) continue;

		int n_track_points = reconPart.getTracks()[0].measurements_size();

		// Calculate rcpart_ member variables for this event
		CalculateParticleValues(reconPart, rcparts);

		// Calculate E/p and isolation fraction for this event
		// Note that the rcpart_ variables are set in CalculateParticleValues
		double recon_EoP = rcpart_sum_cluster_E / edm4hep::utils::magnitude(reconPart.getMomentum());
		double recon_isoE = rcpart_sum_cluster_E / rcpart_isolation_E;

		int found_id = Check_eID(reconPart);
		if ( found_id == 0 )
			e_det.push_back({n_track_points, recon_EoP, recon_isoE});
		else if ( found_id == 11 )
			jet_e_det.push_back({n_track_points, recon_EoP, recon_isoE});
		else if ( found_id == -211 )
			pi_det.push_back({n_track_points, recon_EoP, recon_isoE});
		else
			else_det.push_back({n_track_points, recon_EoP, recon_isoE});

		if ( n_track_points < minTrackPoints ) continue;

		// Apply scattered electron ID cuts
		if(recon_EoP < mEoP_min || recon_EoP > mEoP_max) continue;
		if(recon_isoE < mIsoE) continue;

		// If particle passes cuts, add to output collection
		scatteredElectronCandidates.push_back(reconPart);

	}	

	// std::cout << "Total scattered electron candidates found: " << scatteredElectronCandidates.size() << std::endl;

	// if ( scatteredElectronCandidates.size() == 0 )
	// {
	// 	std::cout << "** No scattered electron candidates found! " << std::endl;
	
	// 	const edm4hep::MCParticleCollection meMC = GetMCElectron();
	// 	const auto& RecoMC = mEvent->get<edm4eic::MCRecoParticleAssociationCollection>("ReconstructedParticleAssociations");

	// 	for(const auto& assoc : RecoMC) 
	// 	{
	// 		auto Sim = assoc.getSim();
	// 		auto Rec = assoc.getRec();

	// 		// if ( Sim == meMC[0] )

	// 		if ( Rec.getClusters().size() != 0 && Rec.getTracks().size() == 0 ) 
	// 			cout << Sim << endl;
	// 	}
	// }
	
	return scatteredElectronCandidates;

}

edm4hep::MCParticleCollection ElectronID::GetMCHadronicFinalState() {

	edm4hep::MCParticleCollection mhMC;
	mhMC.setSubsetCollection();

	const auto& mcparts = mEvent->get<edm4hep::MCParticleCollection>("MCParticles");

	std::vector<edm4hep::MCParticle> mc_hadronic;
	edm4hep::MCParticleCollection meMC = GetMCElectron();

	bool found_scattered_e = false; 
	for(const auto& mcp : mcparts) {
		if (mcp.getGeneratorStatus() == 1)
		{
			if ( meMC.size() == 0 )
				mhMC.push_back(mcp);
			else if (mcp.getObjectID().index != meMC[0].getObjectID().index ) 
				mhMC.push_back(mcp);	
		}
	}

	return mhMC;
}

edm4hep::MCParticleCollection ElectronID::GetMCElectron() {

	edm4hep::MCParticleCollection meMC;
	meMC.setSubsetCollection();
	
	const auto& mcparts = mEvent->get<edm4hep::MCParticleCollection>("MCParticles");
	if ( eScatIndex != -1 )
		meMC.push_back(mcparts[eScatIndex]);

	////

	for (const auto& mcp : mcparts) 
	{
		if (mcp.getPDG() != 11 || mcp.getGeneratorStatus() != 4) 
			continue;

		std::vector<edm4hep::MCParticle> stack;
		stack.insert(stack.end(), mcp.getDaughters().begin(), mcp.getDaughters().end());

		int shortest_gen = 999;
		int generations = 0;
		edm4hep::MCParticle meMC_candidates;
		bool found_candidate = false;

		while (!stack.empty() ) {
			generations++;
			auto cur = stack.back();
			stack.pop_back();

			if (cur.getPDG() == 11 && cur.getGeneratorStatus() == 1) {
				
				if ( generations < shortest_gen )
				{
					shortest_gen = generations;
					meMC_candidates = cur;
					found_candidate = true;
				}
				break;
			}

			const auto& kids = cur.getDaughters();
			if (!kids.empty()) {
				stack.insert(stack.end(), kids.begin(), kids.end());
			}
		}

		if ( found_candidate )
			meMC.push_back(meMC_candidates);
	}

	return meMC;
}

edm4eic::ReconstructedParticleCollection ElectronID::GetTruthReconElectron() {

	// cout << "New process " << endl;

	const edm4hep::MCParticleCollection meMC = GetMCElectron();
	edm4eic::ReconstructedParticleCollection meRecon;
	meRecon.setSubsetCollection();

	if ( meMC.size() == 0 )
		return meRecon; // No MC electron found

	const auto& RecoMC = mEvent->get<edm4eic::MCRecoParticleAssociationCollection>("ReconstructedParticleAssociations");

	for(const auto& assoc : RecoMC) 
	{
		auto e_candidat = assoc.getSim();

		if(assoc.getSim() == meMC[0]) {
			meRecon.push_back(assoc.getRec());
			// break;
		}
	}

	const auto& rcparts = mEvent->get<edm4eic::ReconstructedParticleCollection>("ReconstructedParticles");
	edm4hep::Vector3f lead_pos(meMC[0].getEndpoint().x, meMC[0].getEndpoint().y, meMC[0].getEndpoint().z);
	CheckSurroundingClusters(lead_pos, rcparts);

	// std::cout << meRecon.size() << " reconstructed particles associated with the MC electron." << std::endl;

	return meRecon;
}
	
void ElectronID::CalculateParticleValues(const edm4eic::ReconstructedParticle& rcp,
		const edm4eic::ReconstructedParticleCollection& rcparts) {

	rcpart_sum_cluster_E = 0.;
	rcpart_lead_cluster_E = 0.;
	rcpart_isolation_E = 0.;

	const edm4eic::Cluster* lead_cluster = nullptr;

	for (const auto& cluster : rcp.getClusters()) {
		rcpart_sum_cluster_E += cluster.getEnergy();
		if(cluster.getEnergy() > rcpart_lead_cluster_E) {
			lead_cluster = &cluster;
			rcpart_lead_cluster_E = cluster.getEnergy();
		}
	}

	if(!lead_cluster) return;

	CheckSurroundingClusters(lead_cluster->getPosition(), rcparts);

	return;
}

void ElectronID::CheckSurroundingClusters(const edm4hep::Vector3f& lead_pos,
		const edm4eic::ReconstructedParticleCollection& rcparts) {

	rcpart_isolation_E = 0.;
	rcpart_n_clusters = 0;

	double lead_eta = edm4hep::utils::eta(lead_pos);
	double lead_phi = edm4hep::utils::angleAzimuthal(lead_pos);

	for (const auto& other_rcp : rcparts) {
		for (const auto& other_cluster : other_rcp.getClusters()) {

			const auto& other_pos = other_cluster.getPosition();
			double other_eta = edm4hep::utils::eta(other_pos);
			double other_phi = edm4hep::utils::angleAzimuthal(other_pos);

			double d_eta = other_eta - lead_eta;
			double d_phi = other_phi - lead_phi;

			// Adjust d_phi to be in the range (-pi, pi)
			if (d_phi > M_PI) d_phi-=2*M_PI;
			if (d_phi < -M_PI) d_phi+=2*M_PI;

			double dR = std::sqrt(std::pow(d_eta, 2) + std::pow(d_phi, 2));

			// Check if the cluster is within the isolation cone
			if (dR < mIsoR) {
				rcpart_isolation_E += other_cluster.getEnergy();
				rcpart_n_clusters ++;
			}
		}
	}

	return;
}

void ElectronID::GetEminusPzSum(double &TrackEminusPzSum, double &CalEminusPzSum) {

	const auto& rcparts = mEvent->get<edm4eic::ReconstructedParticleCollection>("ReconstructedParticles");

	for (const auto& reconPart : rcparts) {

		// Require at least one track and one cluster
		if(reconPart.getClusters().size() == 0 || reconPart.getTracks().size() == 0) continue;

		int n_track_points = reconPart.getTracks()[0].measurements_size();
		if ( n_track_points < minTrackPoints ) continue;

		PxPyPzEVector vC(reconPart.getMomentum().x, reconPart.getMomentum().y, reconPart.getMomentum().z, GetCalorimeterEnergy(reconPart));
		vC = boost(vC);
		CalEminusPzSum += (vC.E() - vC.Pz());

		PxPyPzEVector vT(reconPart.getMomentum().x, reconPart.getMomentum().y, reconPart.getMomentum().z, reconPart.getEnergy());
		vT = boost(vT);
		TrackEminusPzSum += (vT.E() - vT.Pz());
	}

	// std::cout << " recon E - Pz sum: " << reconEminusPzSum << std::endl;

	return;
}

edm4eic::ReconstructedParticle ElectronID::SelectHighestPT(const edm4eic::ReconstructedParticleCollection& ecandidates) {

	edm4eic::ReconstructedParticle erec;
	double max_pT = 0.;
	
	for(const auto& ecand : ecandidates) {
		double e_pT = edm4hep::utils::magnitudeTransverse(ecand.getMomentum());
		if(e_pT > max_pT) {
			erec = ecand;
			max_pT = e_pT;
		}
	}

	return erec;

}

double ElectronID::GetCalorimeterEnergy(const edm4eic::ReconstructedParticle& rcp) {

	double sum_cluster_E = 0.;
	for (const auto& cluster : rcp.getClusters()) {
		sum_cluster_E += cluster.getEnergy();
	}
	return sum_cluster_E;

}


