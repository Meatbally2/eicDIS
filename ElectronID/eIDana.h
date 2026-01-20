#ifndef EIDANA_H
#define EIDANA_H

#include "podio/Frame.h"
#include "podio/ROOTReader.h"

#include "../GlobalUtil/AnaManager.cc"
#include "../GlobalUtil/DrawManager.cc"
#include "../GlobalUtil/getBoost.h"
#include "../GlobalUtil/Constants.hh"

#include "ElectronID.cc"

// #include <Math/LorentzRotation.h>
// using ROOT::Math::LorentzRotation;

void eIDana(int Ee, int Eh, int analyse_p, int select_region, int sr, int file0);
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
double EminusPz;

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

PxPyPzEVector vMC_e;
PxPyPzEVector vTRACK_e;
PxPyPzEVector vCLUSTER_e;
std::vector<PxPyPzEVector> vMC_hfs;
std::vector<PxPyPzEVector> vREC_hfs;

TH1D* h_nTPts_e;
TH1D* h_nTPts_pi;
TH1D* h_nTPts_else;
TH1D* h_EoP_e;
TH1D* h_EoP_pi;
TH1D* h_EoP_else;
TH1D* h_isoE_e;
TH1D* h_isoE_pi;
TH1D* h_isoE_else;

TH1D* h_TrackEminusPz;
TH1D* h_CalEminusPz;

TH1D* h_n_scat_elec;
TH2D* h_n_clusters_n_tracks;

TH1D* h_cand_mul;
TH1D* h_cand_mul_eHighPt;
TH1D* h_cand_mul_oHighPt;

#endif