#ifndef BEAMANA_HH
#define BEAMANA_HH

#include "../GlobalUtil/DrawManager.cc"
#include "../GlobalUtil/getBoost.h"
// #include "../GlobalUtil/Constants.hh"
#include "../GlobalUtil/luminosityTable.h"

#include "kinematics.cc"
#include "reconMethod.hh"

#include <Math/LorentzVector.h>
using ROOT::Math::PxPyPzEVector;

#include <Math/LorentzRotation.h>
using ROOT::Math::LorentzRotation;

enum {EL, JB, DA, SIG, ESIG, MC};

void define_histograms();
void plot_energy_and_track();
void format_pad();
void format_hist(TH1* hist);
void draw_eta_axis(TH1* hist);
void ReverseXAxis(TH1 *h);

DrawManager* draw_manager;

struct EPQA
{
   TH1F* h_dE[4];
   TH1F* h_dP[4];
   TH1F* h_dtheta[4];
   TH1F* h_dphi[4];

   TH2F* h_dEvT[4];
   TH2F* h_dPvT[4];

   TH2F* h_dEvT_sum;
   TH2F* h_dPvT_sum;
};
EPQA trackQA;
EPQA clusterQA;

TH1F* h_dhf[4];
TH1F* h_dpt[4];

TH1F* h_hfs_dpz;
TH1F* h_hfs_dpt;
TH1F* h_hfs_de;
TH2F* h_hfs_dpz_t;
TH2F* h_hfs_dpt_t;
TH2F* h_hfs_de_t;

TH2F* h_tde_tmp[2][4];
TH2F* h_tdp_tmp[2][4];

TH2F* h_xq_gen[5];
TH2F* h_xq_eff[5];

void format_eid_histogram(TH1F* &h)
{
    h->GetXaxis()->CenterTitle();
    h->GetXaxis()->SetTitleOffset(1.2);
   
    h->GetYaxis()->SetTitleOffset(1.8);
    h->GetYaxis()->CenterTitle();

    h->GetXaxis()->SetNdivisions(5);
    h->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    h->GetXaxis()->SetLabelSize(20);
    h->GetXaxis()->SetTitleSize(20);
    h->GetXaxis()->SetTitleFont(43);

    h->GetYaxis()->SetNdivisions(5);
    h->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    h->GetYaxis()->SetLabelSize(20);
    h->GetYaxis()->SetTitleSize(20);
    h->GetYaxis()->SetTitleFont(43);

    return;
}

void format_eid_histogram(TH2F* &h)
{
    h->GetXaxis()->CenterTitle();
    h->GetXaxis()->SetTitleOffset(1.2);
   
    h->GetYaxis()->SetTitleOffset(1.8);
    h->GetYaxis()->CenterTitle();

    h->GetXaxis()->SetNdivisions(505);
    h->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    h->GetXaxis()->SetLabelSize(20);
    h->GetXaxis()->SetTitleSize(20);
    h->GetXaxis()->SetTitleFont(43);

    h->GetYaxis()->SetNdivisions(505);
    h->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    h->GetYaxis()->SetLabelSize(20);
    h->GetYaxis()->SetTitleSize(20);
    h->GetYaxis()->SetTitleFont(43);

    return;
}

#endif