#ifndef TAGANA_HH
#define TAGANA_HH

#include "edm4eic/TrackerHitCollection.h"
#include "edm4hep/MCParticleCollection.h"

#include "podio/Frame.h"
#include "podio/ROOTReader.h"

#include "../GlobalUtil/AnaManager.cc"
#include "../GlobalUtil/DrawManager.cc"
#include "../GlobalUtil/Constants.hh"

#include "FarForward.cc"

#include <Math/LorentzVector.h>
using ROOT::Math::PxPyPzEVector;

using namespace std;

struct spectator_info
{
    int pbg;
    int mc_index;
    PxPyPzEVector vec;
    int det_hits[2][4];
    bool tagged;
};

enum SpectatorStatus {
    SPECTATOR_UNCLASSIFIED = 0,
    SPECTATOR_EN_ANCESTRY = 1,
    SPECTATOR_EP_INITIAL = 2,
    SPECTATOR_INCOMPLETE = 3,
    SPECTATOR_MAPPING_FAILED = 4,
    SPECTATOR_EN_KINEMATIC_CLEAN = 5,
    SPECTATOR_EN_KINEMATIC_AMBIGUOUS = 6
};
std::vector<spectator_info*> spec;

void tagAna(int Ee, int Eh, int analyse_p, int select_region, int sr, int file0);
void CreateOutputTree(TString outFileName);
void setup_rp();
void setup_omd();
void setup_zdc();
void process_ff(const podio::Frame* event);
void plot_ff();
void find_spectators(const podio::Frame* event);
TLorentzVector boost_beagle_spectator( const edm4hep::MCParticle& spectator, const edm4hep::MCParticleCollection& mcHeadon);

TFile* outFile;
TTree* outTree;
DrawManager* draw_manager;

FarForward* rpFinder;
FarForward* omdFinder;
FarForward* zdcFinder;
vector<FarForward*> ffFinder;
bool is_free_neutron;

bool is_tagged;
int n_proton_tracks;
double zdc_energy;
int zdc_pbg;
bool fZDCn;
int Spec1_rpHits[4];
int Spec1_omdHits[4];
int Spec2_rpHits[4];
int Spec2_omdHits[4];

std::vector<int> SpecPBG;
std::vector<PxPyPzEVector> SpecVec;
std::vector<bool> SpecTag;

std::vector<int> OtherPBG;
std::vector<PxPyPzEVector> OtherVec;
std::vector<bool> OtherTag;

TH1F* h_xq2_pt;
TH1F* h_xq2_pt_tag;
TH2F* h_xq2_pt_theta;
TH1F* h_tag_mul[2];
TH1F* h_spectator_status;

int struck_type;
int spectator_status;
int n_final_protons;
double ion_momentum_per_nucleon;

#endif
