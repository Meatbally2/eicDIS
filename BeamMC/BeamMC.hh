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
    void GetSpecInfo(std::vector<int> &SpecPBG, std::vector<int> &SpecIndex, std::vector<PxPyPzEVector> &SpecVec, std::vector<int> &OtherPBG, std::vector<int> &OtherIndex, std::vector<PxPyPzEVector> &OtherVec);
    void GetMCinfo(PxPyPzEVector &mc_e,  PxPyPzEVector &mc_p, int &n_pbg);

private:

    const podio::Frame* mEvent;

    std::vector<int> spec_pbg;
    std::vector<int> spec_index;
    std::vector<PxPyPzEVector> spec_vec;

    std::vector<int> other_pbg;
    std::vector<int> other_index;
    std::vector<PxPyPzEVector> other_vec;
    
    double Ee;
    double Eh;
};

#endif