#ifndef ELECTRONID_HH
#define ELECTRONID_HH

#include "podio/Frame.h"

#include "edm4eic/ReconstructedParticleCollection.h"
#include "edm4eic/MCRecoParticleAssociationCollection.h"
#include "edm4hep/MCParticleCollection.h"

#include <Math/LorentzRotation.h>
using ROOT::Math::LorentzRotation;

#include <algorithm>

class ElectronID{

public:

	ElectronID();
	ElectronID(double Ee, double Eh);
       	~ElectronID();

	inline void SetBeamEnergy(double Ee, double Eh) {mEe = Ee; mEh = Eh;}
	inline void SetEoPMin(double eopmin) {mEoP_min = eopmin;}	
	inline void SetDeltaHMin(double deltahmin) {mDeltaH_min = deltahmin;}	
	inline void SetIsolation(double isor, double isoe) {mIsoR = isor; mIsoE = isoe;}	
	inline void SetMinTrackPoints(int minPoints) { minTrackPoints = minPoints; }

	void SetEvent(const podio::Frame* event); 
	void SetBoost(LorentzRotation fboost) { boost = fboost; }
	TCanvas* MakePlots();

	int Check_eID(edm4eic::ReconstructedParticle e_rec);
	edm4hep::MCParticle GetMC(edm4eic::ReconstructedParticle e_rec);
	edm4eic::ReconstructedParticleCollection FindHadronicFinalState(int object_id);
	edm4eic::ReconstructedParticleCollection FindScatteredElectron();	
	edm4eic::ReconstructedParticleCollection GetTruthReconElectron();	
	const edm4hep::MCParticleCollection& GetMCElectron();	
	edm4hep::MCParticleCollection GetMCHadronicFinalState();
	edm4eic::ReconstructedParticle SelectHighestPT(const edm4eic::ReconstructedParticleCollection& rcparts);
	double GetCalorimeterEnergy(const edm4eic::ReconstructedParticle& rcp);
	double GetClusterCone(const edm4eic::ReconstructedParticle& rcp, double frac=1);
	PxPyPzEVector GetMomentumVectorFromCluster(const edm4eic::ReconstructedParticle& rcp, double mass);
	double GetClusterTheta(const edm4eic::ReconstructedParticle& rcp);
	void GetEminusPzSum(double &TrackEminusPzSum, double &CalEminusPzSum);
	void CheckClusters();

	double get_mEoP_min() const { return mEoP_min; }
	double get_mEoP_max() const { return mEoP_max; }
	double get_mDeltaH_min() const { return mDeltaH_min; }
	double get_mDeltaH_max() const { return mDeltaH_max; }
	double get_mIsoR() const { return mIsoR; }
	double get_mIsoE() const { return mIsoE; }
	double get_mEoEH_min() const { return mEoEH_min; }
	double get_mPID_veto() const { return mPID_veto; }
	int GetMinTrackPoints() const { return minTrackPoints; }

	// for HFS QA
	vector<double> hfs_dpt;
	vector<double> hfs_dpz;
	vector<double> hfs_de;
	vector<double> hfs_theta;

	struct DetValues {
		int parType; // 0 for dis electron, rest follow pdg code
		int nTrackPoints;
		double recon_EoP;
		double recon_isoE;
		double recon_EoEH;
		int recon_pID;
		double recon_Le; // likelihood of being an electron
		double recon_Lh; // likelihood of being a hadron
		double eta;
	};
	vector<DetValues> det_val;
	vector<DetValues> mod_val;

	struct PIDMap {
		int index;
		int detID;
		std::string name;
	};
	PIDMap pid_map[4] = {
		{0, 131, "Backward pfRICH"},
		{1, 92, "Barrel TOF"},
		{2, 90, "Barrel hpDIRC"},
		{3, 120, "Forward dRICH"}
	};

	struct PDGMap {
		int index;
		int pdg;
		std::string name;
	};
	PDGMap pdg_map[5] = {
		{0, 11, "electron"},
		{1, -211, "pion"},
		{2, -321, "kaon"},
		{3, -2212, "proton"},
		{0, 0, "electron"}
	};

	TH1F* h_pid_score[4][4][4]; // [MC][detector][pid]
	TCanvas* c_pid_score;

	double rcpart_n_clusters;

	bool meMCValid = false;
	edm4hep::MCParticleCollection meMC;

private:

	const podio::Frame* mEvent;

	double mEe;
	double mEh;
	LorentzRotation boost;

	double mEoP_min;
	double mEoP_max;
	double mDeltaH_min;
	double mDeltaH_max;
	double mIsoR;
	double mIsoE;
	double mEcone;
	double mEoEH_min;
	double mPID_veto;
	int minTrackPoints = 3;
	
	void CalculateParticleValues(const edm4eic::ReconstructedParticle& rcp,
		const edm4eic::ReconstructedParticleCollection& rcparts);
	void CheckSurroundingClusters(const edm4hep::Vector3f& lead_pos,
		const edm4eic::ReconstructedParticleCollection& rcparts);

	double rcpart_sum_cluster_E;
	double rcpart_lead_cluster_E;
	double rcpart_isolation_E;
	double rcpart_sum_cluster_H;
	double rcpart_sum_cluster_E_wider;
	double rcpart_sum_cluster_H_wider;

	int eScatIndex;
};

#endif
