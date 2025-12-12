#ifndef TAGANA_HH
#define TAGANA_HH

#include "edm4eic/TrackerHitCollection.h"

#include "podio/Frame.h"
#include "podio/ROOTReader.h"

#include "../GlobalUtil/AnaManager.cc"
#include "../GlobalUtil/Constants.hh"
#include "../GlobalUtil/ePICStyle.c"

#include "FarForward.cc"

using namespace std;

void tagAna(int Ee, int Eh, int select_region, int sr, int all_file, int analyse_p);
void CreateOutputTree(TString outFileName);
void setup_rp();
void setup_omd();
void setup_zdc();
void process_ff(const podio::Frame* event);
void plot_ff();

TFile* outFile;
TTree* outTree;

FarForward* rpFinder;
FarForward* omdFinder;
FarForward* zdcFinder;
vector<FarForward*> ffFinder;
bool is_free_neutron;

bool is_tagged;
int n_proton_tracks;
double zdc_energy;
int zdc_pbg;

#endif