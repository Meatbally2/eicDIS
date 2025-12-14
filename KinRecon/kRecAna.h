#ifndef BEAMANA_HH
#define BEAMANA_HH

#include "../GlobalUtil/DrawManager.cc"
#include "../GlobalUtil/getBoost.h"
// #include "../GlobalUtil/Constants.hh"

#include "kinematics.cc"
#include "reconMethod.hh"

#include <Math/LorentzVector.h>
using ROOT::Math::PxPyPzEVector;

#include <Math/LorentzRotation.h>
using ROOT::Math::LorentzRotation;

enum {EL, JB, DA, SIG, ESIG, MC};

void define_histograms();
void plot_energy_and_track(TCanvas* &c1, TCanvas* &c2, TCanvas* &c3, TCanvas* &c4);

TH1F* h_dEc[4];
TH1F* h_dEt[4];
TH1F* h_dPc[4];
TH1F* h_dPt[4];
TH1F* h_dt[4];
TH1F* h_dp[4];
TH2F* h_tde[2];
TH2F* h_tdp[2];
TH2F* h_tdphi;
TH2F* h_pde;
TH2F* h_ede;
TH2F* h_pt;
TH2F* h_et;

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