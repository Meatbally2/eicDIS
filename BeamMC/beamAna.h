#ifndef BEAMANA_HH
#define BEAMANA_HH

#include "podio/Frame.h"
#include "podio/ROOTReader.h"

#include "../GlobalUtil/AnaManager.cc"
#include "../GlobalUtil/Constants.hh"

#include "BeamMC.cc"

#include <Math/LorentzVector.h>
using ROOT::Math::PxPyPzEVector;

void beamAna(int Ee, int Eh, int analyse_p, int select_region, int sr, int file0);
void CreateOutputTree(TString outFileName);
void ResetVariables();

void DefineHistograms();

TFile* outFile;
TTree* outTree;

int N_PDG;
PxPyPzEVector vectE;
PxPyPzEVector vectN;

#endif