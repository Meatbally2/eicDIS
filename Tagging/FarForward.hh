#ifndef FARFORWARD_HH
#define FARFORWARD_HH

#include "podio/Frame.h"

#include "edm4hep/Constants.h"
#include "edm4eic/TrackerHitCollection.h"
#include "edm4eic/TrackerHitData.h"
#include "edm4eic/MCRecoTrackerHitAssociation.h"
#include "edm4eic/ClusterCollection.h"
#include "edm4eic/ReconstructedParticleCollection.h"

#include "../GlobalUtil/drawHelper.h"
#include "../GlobalUtil/Constants.hh"

class FarForward{

public:

    FarForward(std::string d_name, std::string b_name);
    ~FarForward();

    edm4eic::TrackerHitCollection* det[4];

    void SetEvent(const podio::Frame* event); 
    void GetHits(std::vector<int> mc_index);
    double GetEnergy();

    void set_rotation(double r);
    void set_xShift(int p, double xs);
    void set_zRange(int p, double min, double max);
    int GetMCHits(int s, int d);
    void ClearHits();

    int get_n_tracks();

    TH2F* h_det_xy[4];
    TH1F* h_cluster[2];
    TH1F* h_cluster_sum;

    TH1F* hx_raw;
    TH1F* hz_raw;

    void define_histograms();
    void fill_hit_histograms();
    void fill_energy_histograms();
    std::vector<TCanvas*> draw_histograms();

    std::vector<float> ZDC_E;
    std::vector<int> ZDC_PBG;
    bool is_ZDC_neutron;

private:

    const podio::Frame* mEvent;

    std::string det_name;
    std::string det_title;
    std::string branch_name;

    double rotate;
    double xShift[4];
    double zRange[4][2];
    std::vector<int> spec_hit[4]; // [detector_index][spec_index]
};

#endif