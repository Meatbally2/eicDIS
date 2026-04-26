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
std::vector<spectator_info*> spec;

void tagAna(int Ee, int Eh, int analyse_p, int select_region, int sr, int file0);
void CreateOutputTree(TString outFileName);
void setup_rp();
void setup_omd();
void setup_zdc();
void process_ff(const podio::Frame* event);
void plot_ff();
void find_spectators(const podio::Frame* event);

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

TH1F* h_xq2_pt;
TH1F* h_xq2_pt_tag;
TH2F* h_xq2_pt_theta;
TH1F* h_tag_mul[2];

int struck_type;

#endif