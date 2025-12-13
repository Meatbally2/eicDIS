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

	mMCElectronCached = false;
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

	mMCElectronCached = false;
}

ElectronID::~ElectronID() {
}


void ElectronID::SetEvent(const podio::Frame* event) {

	mEvent = event;
	eScatIndex = -1;
	hfs_dpt.clear();
	hfs_dpz.clear();
	hfs_de.clear();
	hfs_theta.clear();
	e_det.clear();
	pi_det.clear();
	else_det.clear();
	mMCElectronCached = false;
	mCachedMCElectron.clear();
}

edm4hep::MCParticle ElectronID::GetMC(const edm4eic::ReconstructedParticle& e_rec) {

	auto& RecoMC = mEvent->get<edm4eic::MCRecoParticleAssociationCollection>("ReconstructedParticleAssociations");
	for(const auto& assoc : RecoMC) {
		if(assoc.getRec() == e_rec)
			return assoc.getSim();
	}

	return edm4hep::MCParticle();
}

int ElectronID::Check_eID(const edm4eic::ReconstructedParticle& e_rec) {

	const auto& meMC = GetMCElectron();
	auto& RecoMC = mEvent->get<edm4eic::MCRecoParticleAssociationCollection>("ReconstructedParticleAssociations");
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

	auto& EcalEndcapNClusters = mEvent->get<edm4eic::ClusterCollection>("EcalEndcapNClusters");
	auto& EcalBarrelScFiClusters = mEvent->get<edm4eic::ClusterCollection>("EcalBarrelScFiClusters");
	auto& EcalEndcapPClusters = mEvent->get<edm4eic::ClusterCollection>("EcalEndcapPClusters");

	std::cout << " Number of clusters in EcalEndcapN: " << EcalEndcapNClusters.size() << std::endl;
	std::cout << " Number of clusters in EcalBarrelScFi: " << EcalBarrelScFiClusters.size() << std::endl;
	std::cout << " Number of clusters in EcalEndcapP: " << EcalEndcapPClusters.size() << std::endl;

	return;
}

edm4eic::ReconstructedParticleCollection ElectronID::FindHadronicFinalState(bool use_mc, int object_id, bool is_print) {

	edm4eic::ReconstructedParticleCollection meRecon;
	meRecon->setSubsetCollection();

	edm4hep::MCParticleCollection meMiss;
	meMiss->setSubsetCollection();

	auto& rcparts = mEvent->get<edm4eic::ReconstructedParticleCollection>("ReconstructedParticles");

	// Reserve capacity for vector members to reduce memory reallocations
	size_t estimated_size = rcparts.size() > 0 ? rcparts.size() - 1 : 0;
	hfs_dpt.reserve(estimated_size);
	hfs_dpz.reserve(estimated_size);
	hfs_de.reserve(estimated_size);
	hfs_theta.reserve(estimated_size);

	if ( is_print )
		cout << " new HFS loop " << endl;

	if ( use_mc )
	{
		const auto& meMC = GetMCElectron();
		auto& RecoMC = mEvent->get<edm4eic::MCRecoParticleAssociationCollection>("ReconstructedParticleAssociations");

		if ( is_print )
		{
			PxPyPzEVector v(meMC[0].getMomentum().x, meMC[0].getMomentum().y, meMC[0].getMomentum().z, meMC[0].getEnergy());
			v = boost(v);
			cout << "scat e .. " << Form("PDG %d, ", meMC[0].getPDG()) <<  Form("pt %f, ", v.Pt()) << Form("pz %f, ", v.Z()) << Form("E %f, ", v.E()) << Form("theta %f", v.Theta()) << endl;
		}

		for(const auto& assoc : RecoMC)
		{
			if(assoc.getSim() != meMC[0] && assoc.getSim().getGeneratorStatus() == 1)
			{
				PxPyPzEVector v(assoc.getSim().getMomentum().x, assoc.getSim().getMomentum().y, assoc.getSim().getMomentum().z, assoc.getSim().getEnergy());
				PxPyPzEVector u(assoc.getRec().getMomentum().x, assoc.getRec().getMomentum().y, assoc.getRec().getMomentum().z, assoc.getRec().getEnergy());
				PxPyPzEVector c(assoc.getRec().getMomentum().x, assoc.getRec().getMomentum().y, assoc.getRec().getMomentum().z, GetCalorimeterEnergy(assoc.getRec()));

				hfs_dpt.push_back((u.Pt()-v.Pt())/v.Pt());
				hfs_dpz.push_back((u.Z()-v.Z())/v.Z());
				hfs_de.push_back((u.E()-v.E())/v.E());
				hfs_theta.push_back(v.Theta()*(180./M_PI));

				v = boost(v);
				u = boost(u);
				c = boost(c);

				meRecon.push_back(assoc.getRec());

				if ( is_print )
				{
					cout << "selected .." << endl;
					cout << " MC info .. " << Form("PDG %d, ", assoc.getSim().getPDG()) <<  Form("pt %f, ", v.Pt()) << Form("pz %f, ", v.Z()) << Form("E %f, ", v.E()) << Form("theta %f", v.Theta()) << endl;
					cout << " REC info .. " << Form("PDG %d, ", assoc.getRec().getPDG()) <<  Form("pt %f, ", u.Pt()) << Form("pz %f, ", u.Z()) << Form("E %f, ", u.E()) << Form("CAL E %f, ", c.E()) << Form("theta %f", u.Theta()) << endl;
				}

			}
		}

		if ( is_print )
		{
			auto& mcparts = mEvent->get<edm4hep::MCParticleCollection>("MCParticles");
			for(const auto& mcp : mcparts) 
			{
				if (mcp.getGeneratorStatus() == 1)
				{
					bool selected = (mcp == meMC[0]);
					if ( !selected )
					{
						for(const auto& assoc : RecoMC)
							if ( mcp == assoc.getSim() )
								selected = true;
					}

					if ( !selected )
					{
						PxPyPzEVector v(mcp.getMomentum().x, mcp.getMomentum().y, mcp.getMomentum().z, mcp.getEnergy());
						v = boost(v);
						cout << "missed .. " << Form("PDG %d, ", mcp.getPDG()) <<  Form("pt %f, ", v.Pt()) << Form("pz %f, ", v.Z()) << Form("E %f, ", v.E()) << Form("theta %f", v.Theta()) << endl;
					}
				}
			}
		}
	}
	else
	{ 
		for(const auto& mcp : rcparts) {
			if ( mcp.getObjectID().index != object_id )
				meRecon.push_back(mcp);
		}
	}

	return meRecon;
}

edm4eic::ReconstructedParticleCollection ElectronID::FindScatteredElectron() {

	// Get all the edm4eic objects needed for electron ID
	auto& rcparts = mEvent->get<edm4eic::ReconstructedParticleCollection>("ReconstructedParticles");

	// Create collection for storing scattered electron candidates
	// (subset collection of ReconstructedParticleCollection)
	edm4eic::ReconstructedParticleCollection scatteredElectronCandidates;
	scatteredElectronCandidates->setSubsetCollection();

	// Reserve capacity for vector members to reduce memory reallocations
	e_det.reserve(rcparts.size());
	pi_det.reserve(rcparts.size());
	else_det.reserve(rcparts.size());

	// Loop over all ReconstructedParticles for this event
	for (const auto& reconPart : rcparts) {

		// Require negative particle
		if(reconPart.getCharge() >= 0) continue;

		// Require at least one track and one cluster
		if(reconPart.getClusters().size() == 0 || reconPart.getTracks().size() == 0) continue;

		// Calculate rcpart_ member variables for this event
		CalculateParticleValues(reconPart, rcparts);

		// Calculate E/p and isolation fraction for this event
		// Note that the rcpart_ variables are set in CalculateParticleValues
		double recon_EoP = rcpart_sum_cluster_E / edm4hep::utils::magnitude(reconPart.getMomentum());
		double recon_isoE = rcpart_sum_cluster_E / rcpart_isolation_E;

		if ( Check_eID(reconPart) == 0 )
			e_det.push_back({recon_EoP, recon_isoE});
		else if ( Check_eID(reconPart) == -211 )
			pi_det.push_back({recon_EoP, recon_isoE});
		else
			else_det.push_back({recon_EoP, recon_isoE});


		// Apply scattered electron ID cuts
		if(recon_EoP < mEoP_min || recon_EoP > mEoP_max) continue;
		if(recon_isoE < mIsoE) continue;

		// If particle passes cuts, add to output collection
		scatteredElectronCandidates.push_back(reconPart);

	}	

	return scatteredElectronCandidates;

}

edm4hep::MCParticleCollection ElectronID::GetMCHadronicFinalState() {

	edm4hep::MCParticleCollection mhMC;
	mhMC->setSubsetCollection();

	auto& mcparts = mEvent->get<edm4hep::MCParticleCollection>("MCParticles");

	const auto& meMC = GetMCElectron();

	for(const auto& mcp : mcparts) {
		if (mcp.getGeneratorStatus() == 1 && mcp.getObjectID().index != meMC[0].getObjectID().index ) 
			mhMC.push_back(mcp);	
	}

	return mhMC;
}

const edm4hep::MCParticleCollection& ElectronID::GetMCElectron() const {

	// Return cached result if available
	if (mMCElectronCached) {
		return mCachedMCElectron;
	}

	mCachedMCElectron.clear();
	mCachedMCElectron->setSubsetCollection();

	auto& mcparts = mEvent->get<edm4hep::MCParticleCollection>("MCParticles");
	if ( eScatIndex != -1 )
		mCachedMCElectron.push_back(mcparts[eScatIndex]);

	////
	for(const auto& mcp : mcparts) {
		if(mcp.getPDG() == 11 && mcp.getGeneratorStatus() == 1) {
			mCachedMCElectron.push_back(mcp); // For now, just take first electron
			break;
		}
	}

	mMCElectronCached = true;
	return mCachedMCElectron;
}

edm4eic::ReconstructedParticleCollection ElectronID::GetTruthReconElectron() {

	const auto& meMC = GetMCElectron();
	edm4eic::ReconstructedParticleCollection meRecon;
	meRecon->setSubsetCollection();

	auto& RecoMC = mEvent->get<edm4eic::MCRecoParticleAssociationCollection>("ReconstructedParticleAssociations");

	for(const auto& assoc : RecoMC) {
		if(assoc.getSim() == meMC[0]) {
			meRecon.push_back(assoc.getRec());
			break;
		}
	}

	return meRecon;
}
	
void ElectronID::CalculateParticleValues(const edm4eic::ReconstructedParticle& rcp,
		const edm4eic::ReconstructedParticleCollection& rcparts) {

	rcpart_sum_cluster_E = 0.;
	rcpart_lead_cluster_E = 0.;
	rcpart_isolation_E = 0.;
	rcpart_deltaH = 0.;

	const edm4eic::Cluster* lead_cluster = nullptr;

	// Find lead cluster and sum cluster energy
	const auto& clusters = rcp.getClusters();
	for (const auto& cluster : clusters) {
		const double cluster_E = cluster.getEnergy();
		rcpart_sum_cluster_E += cluster_E;
		if(cluster_E > rcpart_lead_cluster_E) {
			lead_cluster = &cluster;
			rcpart_lead_cluster_E = cluster_E;
		}
	}

	if(!lead_cluster) return;

	// Precompute lead cluster position values
	const auto& lead_pos = lead_cluster->getPosition();
	const double lead_eta = edm4hep::utils::eta(lead_pos);
	const double lead_phi = edm4hep::utils::angleAzimuthal(lead_pos);
	const double mIsoR2 = mIsoR * mIsoR;  // Compare squared distance to avoid sqrt

	// Calculate isolation energy
	for (const auto& other_rcp : rcparts) {
		if (&other_rcp == &rcp) continue;  // Skip the same particle

		const auto& other_clusters = other_rcp.getClusters();
		for (const auto& other_cluster : other_clusters) {

			const auto& other_pos = other_cluster.getPosition();
			const double other_eta = edm4hep::utils::eta(other_pos);
			const double other_phi = edm4hep::utils::angleAzimuthal(other_pos);

			double d_eta = other_eta - lead_eta;
			double d_phi = other_phi - lead_phi;

			// Adjust d_phi to be in the range (-pi, pi)
			if (d_phi > M_PI) d_phi -= 2*M_PI;
			else if (d_phi < -M_PI) d_phi += 2*M_PI;

			// Use squared distance to avoid expensive sqrt
			const double dR2 = d_eta*d_eta + d_phi*d_phi;

			// Check if the cluster is within the isolation cone
			if (dR2 < mIsoR2) {
				rcpart_isolation_E += other_cluster.getEnergy();
			}
		}
	}
}

void ElectronID::GetEminusPzSum(double &TrackEminusPzSum, double &CalEminusPzSum) {

	auto& rcparts = mEvent->get<edm4eic::ReconstructedParticleCollection>("ReconstructedParticles");

	for (const auto& reconPart : rcparts) {
		PxPyPzEVector vC(reconPart.getMomentum().x, reconPart.getMomentum().y, reconPart.getMomentum().z, GetCalorimeterEnergy(reconPart));
		vC = boost(vC);
		CalEminusPzSum += (vC.E() - vC.Pz());

		PxPyPzEVector vT(reconPart.getMomentum().x, reconPart.getMomentum().y, reconPart.getMomentum().z, reconPart.getEnergy());
		vT = boost(vT);
		TrackEminusPzSum += (vT.E() - vT.Pz());
	}
}

edm4eic::ReconstructedParticle ElectronID::SelectHighestPT(const edm4eic::ReconstructedParticleCollection& ecandidates) {

	edm4eic::ReconstructedParticle erec;
	double max_pT = 0.;
	
	for(auto ecand : ecandidates) {
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