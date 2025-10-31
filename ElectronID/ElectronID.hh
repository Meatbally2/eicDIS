#ifndef ELECTRONID_HH
#define ELECTRONID_HH

#include "podio/Frame.h"

#include "edm4eic/HadronicFinalStateCollection.h"
#include "edm4eic/ReconstructedParticleCollection.h"
#include "edm4hep/MCParticleCollection.h"

#include <Math/LorentzRotation.h>
using ROOT::Math::LorentzRotation;


class ElectronID{

public:

	ElectronID();
	ElectronID(double Ee, double Eh);
       	~ElectronID();

	inline void SetBeamEnergy(double Ee, double Eh) {mEe = Ee; mEh = Eh;}
	inline void SetEoPMin(double eopmin) {mEoP_min = eopmin;}	
	inline void SetDeltaHMin(double deltahmin) {mDeltaH_min = deltahmin;}	
	inline void SetIsolation(double isor, double isoe) {mIsoR = isor; mIsoE = isoe;}	

	void SetEvent(const podio::Frame* event); 

	int Check_eID(edm4eic::ReconstructedParticle e_rec);
	edm4eic::ReconstructedParticleCollection FindHadronicFinalState(bool use_mc, int object_id, bool is_print, LorentzRotation boost);
	edm4eic::ReconstructedParticleCollection FindScatteredElectron();	
	edm4eic::ReconstructedParticleCollection GetTruthReconElectron();	
	edm4hep::MCParticleCollection GetMCElectron();	
	edm4hep::MCParticleCollection GetMCHadronicFinalState();
	edm4eic::ReconstructedParticle SelectHighestPT(const edm4eic::ReconstructedParticleCollection& rcparts);
	double GetCalorimeterEnergy(const edm4eic::ReconstructedParticle& rcp);

	// for HFS QA
	vector<double> hfs_dpt;
	vector<double> hfs_dpz;
	vector<double> hfs_de;
	vector<double> hfs_theta;

private:

	const podio::Frame* mEvent;

	double mEe;
	double mEh;

	double mEoP_min;
	double mEoP_max;
	double mDeltaH_min;
	double mDeltaH_max;
	double mIsoR;
	double mIsoE;
	
	void CalculateParticleValues(const edm4eic::ReconstructedParticle& rcp,
		const edm4eic::ReconstructedParticleCollection& rcparts);

	double rcpart_sum_cluster_E;
	double rcpart_lead_cluster_E;
	double rcpart_isolation_E;
	double rcpart_deltaH;

	int eScatIndex;
};

#endif
