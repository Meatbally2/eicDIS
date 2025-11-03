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

	mEoP_min = 0.9;
	mEoP_max = 1.2;

	mDeltaH_min = 0.*mEe;
	mDeltaH_min = 5.*mEe;
		
	mIsoR = 0.4;
	mIsoE = 0.9;


}

ElectronID::ElectronID(double Ee, double Eh) {

	mEe = Ee;
	mEh = Eh;

	mEoP_min = 0.9;
	mEoP_max = 1.2;

	mDeltaH_min = 0.*mEe;
	mDeltaH_min = 5.*mEe;
		
	mIsoR = 0.4;
	mIsoE = 0.9;

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
}

int ElectronID::Check_eID(edm4eic::ReconstructedParticle e_rec) {

	edm4hep::MCParticleCollection meMC = GetMCElectron();
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

edm4eic::ReconstructedParticleCollection ElectronID::FindHadronicFinalState(bool use_mc, int object_id, bool is_print, LorentzRotation boost) {

	// edm4eic::HadronicFinalStateCollection meRecon;
	edm4eic::ReconstructedParticleCollection meRecon;
	meRecon->setSubsetCollection();

	edm4hep::MCParticleCollection meMiss;
	meMiss->setSubsetCollection();

	// auto& rcparts = mEvent->get<edm4eic::HadronicFinalStateCollection>("HadronicFinalState");
	auto& rcparts = mEvent->get<edm4eic::ReconstructedParticleCollection>("ReconstructedParticles");

	// cout << "ef_rc_id after " << object_id << endl;

	if ( is_print )
		cout << " new HFS loop " << endl;

	if ( use_mc )
	{
		edm4hep::MCParticleCollection meMC = GetMCElectron();
		auto& RecoMC = mEvent->get<edm4eic::MCRecoParticleAssociationCollection>("ReconstructedParticleAssociations");

		if ( is_print )
		{
			// cout << "scat e\n" << meMC[0] << endl;
			PxPyPzEVector v(meMC[0].getMomentum().x, meMC[0].getMomentum().y, meMC[0].getMomentum().z, meMC[0].getEnergy());
			v = boost(v);
			cout << "scat e .. " << Form("PDG %d, ", meMC[0].getPDG()) <<  Form("pt %f, ", v.Pt()) << Form("pz %f, ", v.Z()) << Form("E %f, ", v.E()) << Form("theta %f", v.Theta()) << endl;
		}

		// if ( is_print )
		// 	cout << "hfs selected .. \n" << endl;
			
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

				// if ( (u.Pt()-v.Pt())/v.Pt() > 0.3 || (u.Pt()-v.Pt())/v.Pt() < -0.3 )
				// 	continue; 
				
				// if ( (u.Z()-v.Z())/v.Z() > 0.3 || (u.Z()-v.Z())/v.Z() < -0.3 )
				// 	continue;

				// if ( (u.E()-v.E())/v.E() > 0.3 || (u.E()-v.E())/v.E() < -0.3 )
				// 	continue;

				
				v = boost(v);
				u = boost(u);
				c = boost(c);

				// if ( assoc.getSim().getMomentum().z > 0 && assoc.getRec().getMomentum().z < 0 )
				// 	continue;

				// if ( assoc.getSim().getMomentum().z < 0 && assoc.getRec().getMomentum().z > 0 )
				// 	continue;

				meRecon.push_back(assoc.getRec());

				if ( is_print )
				{
					// cout << "selected ..\n" << assoc.getSim() << endl;
					cout << "selected .." << endl;
					cout << " MC info .. " << Form("PDG %d, ", assoc.getSim().getPDG()) <<  Form("pt %f, ", v.Pt()) << Form("pz %f, ", v.Z()) << Form("E %f, ", v.E()) << Form("theta %f", v.Theta()) << endl;
					cout << " REC info .. " << Form("PDG %d, ", assoc.getRec().getPDG()) <<  Form("pt %f, ", u.Pt()) << Form("pz %f, ", u.Z()) << Form("E %f, ", u.E()) << Form("CAL E %f, ", c.E()) << Form("theta %f", u.Theta()) << endl;
				}
					
			}

			// cout << "recon asso " << assoc.getSim().getPDG() << endl;
		}

		if ( is_print )
		{
			// cout << "hfs missed .. \n" << endl;

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
						// cout << "missed ..\n" << mcp << endl;
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

	// std::cout << "Finding scattered electron candidates..." << std::endl;

	// Get all the edm4eic objects needed for electron ID
	auto& rcparts = mEvent->get<edm4eic::ReconstructedParticleCollection>("ReconstructedParticles");
	
	// Create collection for storing scattered electron candidates
	// (subset collection of ReconstructedParticleCollection)
	edm4eic::ReconstructedParticleCollection scatteredElectronCandidates;
	scatteredElectronCandidates->setSubsetCollection();

	// Loop over all ReconstructedParticles for this event
	for (const auto& reconPart : rcparts) {

		// std::cout << "par id: " << reconPart.getPDG() << " cluster size: " << reconPart.getClusters().size() << ", track size: " << reconPart.getTracks().size() << std::endl;

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

	std::vector<edm4hep::MCParticle> mc_hadronic;
	edm4hep::MCParticleCollection meMC = GetMCElectron();
	
	bool found_scattered_e = false; 
	for(const auto& mcp : mcparts) {
		if (mcp.getGeneratorStatus() == 1 && mcp.getObjectID().index != meMC[0].getObjectID().index ) 
			mhMC.push_back(mcp);	
	}

	return mhMC;
}

edm4hep::MCParticleCollection ElectronID::GetMCElectron() {

	edm4hep::MCParticleCollection meMC;
	meMC->setSubsetCollection();
	
	auto& mcparts = mEvent->get<edm4hep::MCParticleCollection>("MCParticles");
	if ( eScatIndex != -1 )
		meMC.push_back(mcparts[eScatIndex]);

	////
	for(const auto& mcp : mcparts) {
		if(mcp.getPDG() == 11 && mcp.getGeneratorStatus() == 1) {
			meMC.push_back(mcp); // For now, just take first electron
			break;
		}
	}
	////

	// std::vector<edm4hep::MCParticle> mc_electrons;
/*	
	// cout << "new" << endl;
	for(const auto& mcp : mcparts) 
	{
		if ( mcp.getGeneratorStatus() == 4 && mcp.getPDG() == 11 )
        {
			for (auto md = mcp.daughters_begin(), mend = mcp.daughters_end(); md != mend; ++md) 
			{
				cout << "new daugter " << endl;
				int mdIndex = md->getObjectID().index;
				auto daughter = mcparts[mdIndex];
				cout << daughter << endl;

				// edm4hep::MCParticle daughter = mcp.getDaughters(0); // copy of mcp
				for (auto it = daughter.daughters_begin(), end = daughter.daughters_end(); it != end; ++it) 
				{
					int Index = it->getObjectID().index;
					auto dau = mcparts[Index];
					cout << "dau " << endl;
					cout << dau << endl;
					if ( dau.getGeneratorStatus() == 1 && dau.getPDG() == 11)
					{
						meMC.push_back(dau);
						eScatIndex = Index;
						
						// cout << dau << endl;
						// for (auto dit = dau.parents_begin(), dend = dau.parents_end(); dit != dend; ++dit) 
						// {
						// 	int dauIndex = dit->getObjectID().index;
						// 	auto dpar = mcparts[dauIndex];
						// 	cout << dpar << "dau parent" << endl << dpar.getParents(0) << endl;
						// }
						break;
					}
				}
			}
			break;
		}
	}
*/
	// cout << "full e list" << endl;
	// for(const auto& mcp : mcparts) 
	// {
	// 	if ( mcp.getGeneratorStatus() == 1 && mcp.getPDG() == 11)
	// 	{
	// 		cout << mcp << endl;
	// 	}
	// }
	

	return meMC;
}

edm4eic::ReconstructedParticleCollection ElectronID::GetTruthReconElectron() {

	edm4hep::MCParticleCollection meMC = GetMCElectron();
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

	for (const auto& cluster : rcp.getClusters()) {
		rcpart_sum_cluster_E += cluster.getEnergy();
		if(cluster.getEnergy() > rcpart_lead_cluster_E) {
			lead_cluster = &cluster;
			rcpart_lead_cluster_E = cluster.getEnergy();
		}
	}

	if(!lead_cluster) return;

	const auto& lead_pos = lead_cluster->getPosition();
	double lead_eta = edm4hep::utils::eta(lead_pos);
	double lead_phi = edm4hep::utils::angleAzimuthal(lead_pos);

	for (const auto& other_rcp : rcparts) {
		if (&other_rcp == &rcp) continue;  // Skip the same particle

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
			}
		}
	}

	return;
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


