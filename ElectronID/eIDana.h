#ifndef EIDANA_H
#define EIDANA_H

#include "podio/Frame.h"
#include "podio/ROOTReader.h"

#include "../GlobalUtil/AnaManager.cc"
#include "../GlobalUtil/getBoost.h"
#include "../GlobalUtil/Constants.hh"
#include "../GlobalUtil/ePICStyle.c"

#include "ElectronID.cc"

// #include <Math/LorentzRotation.h>
// using ROOT::Math::LorentzRotation;

void eIDana(int Ee, int Eh, int select_region, int sr, int is_truth_eID, int file0, int analyse_p);
// void eIDana(int Ee, int Eh, std::string ev_type, int is_truth_eID, int analyse_p);
void CreateOutputTree(TString outFileName);
void ResetVariables();

void CalculateElectronKinematics(double fEe, double fEh, TLorentzVector kf, double& xB, double& Q2, double& W2, double& y, double& nu);
TLorentzVector GetHadronBeam(double fEh);

void DefineHistograms();
void DrawVerticalLine(TCanvas* &c, double x_pos, double y_max);
void DrawParComparison(TCanvas* &c, TH1D* &h1, TH1D* &h2, TH1D* &h3, double& draw_max);
void DrawTCComparison(TCanvas* &c, TH1D* &h1, TH1D* &h2, double& draw_max);

TFile* outFile;
TTree* outTree;

int eID_status;
int mc_PDG;

double mc_xB;
double mc_Q2;
double mc_y;
double mc_W2;
double mc_nu;

double rec_xB;
double rec_Q2;
double rec_W2;
double rec_y;
double rec_nu;

TH1D* h_EoP_e;
TH1D* h_EoP_pi;
TH1D* h_EoP_else;
TH1D* h_isoE_e;
TH1D* h_isoE_pi;
TH1D* h_isoE_else;

TH1D* h_TrackEminusPz;
TH1D* h_CalEminusPz;

TH2D* h_n_clusters_n_tracks;

TH1D* h_cand_mul;
TH1D* h_cand_mul_eHighPt;
TH1D* h_cand_mul_oHighPt;

#endif