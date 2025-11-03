#ifndef EIDANA_H
#define EIDANA_H

#include "podio/Frame.h"
#include "podio/ROOTReader.h"

#include "../GlobalUtil/AnaManager.cc"
#include "../GlobalUtil/Constants.hh"
#include "../GlobalUtil/ePICStyle.c"

#include "ElectronID.cc"

void eIDana(int Ee, int Eh, int select_region, int sr, int is_truth_eID, int all_file, int analyse_p);
void CreateOutputTree(TString outFileName);
void ResetVariables();

void CalculateElectronKinematics(double fEe, double fEh, TLorentzVector kf, double& xB, double& Q2, double& W2, double& y, double& nu);
TLorentzVector GetHadronBeam(double fEh);

void DefineHistograms();
void DrawComparison(TCanvas* c, TH1D* &h1, TH1D* &h2, TH1D* &h3, double& draw_max);

TFile* outFile;
TTree* outTree;
// ElectronID* eFinder;

int eID_status;
int mc_PBG;

double mc_xB;
double mc_Q2;
double mc_y;
double mc_W2;
double mc_nu;

enum { NO_RECON, FOUND_E, FOUND_PI, FOUND_OTHERS, NO_FOUND };

TH1D* h_EoP_e;
TH1D* h_EoP_pi;
TH1D* h_EoP_else;
TH1D* h_isoE_e;
TH1D* h_isoE_pi;
TH1D* h_isoE_else;

TH2D* h_n_clusters_n_tracks;

#endif