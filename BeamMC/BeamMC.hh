#ifndef BEAMMC_HH
#define BEAMMC_HH

// #include "../formats/constants.hh"
#include "podio/Frame.h"

#include "edm4hep/MCParticleCollection.h"

#include <Math/LorentzVector.h>
using ROOT::Math::PxPyPzEVector;

class BeamMC{

public:

    BeamMC(double Ee_, double Eh_);
    ~BeamMC();

    void SetEvent(const podio::Frame* event); 
    void GetMCinfo(PxPyPzEVector &mc_e,  PxPyPzEVector &mc_p, int &n_pbg);

private:

    const podio::Frame* mEvent;
    std::vector<TLorentzVector> spectator_protons;

    double Ee;
    double Eh;
};

#endif