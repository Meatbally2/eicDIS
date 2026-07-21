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

void eIDana(int Ee, int Eh, int beam_type, int select_region=0, int sr=0, int file0=-1);
// void eIDana(int Ee, int Eh, std::string ev_type, int is_truth_eID, int analyse_p);
void CreateOutputTree(TString outFileName);
void ResetVariables();

void CalculateElectronKinematics(double fEe, double fEh, TLorentzVector kf, double& xB, double& Q2, double& W2, double& y, double& nu);
TLorentzVector GetHadronBeam(double fEh);

void DefineHistograms();
void DrawVerticalLine(TCanvas* &c, double x_pos, double y_max);
void DrawParComparison(TCanvas* &c, TH1D* &h1, TH1D* &h2, TH1D* &h3, TH1D* &h4, double& draw_max);
void DrawTCComparison(TCanvas* &c, TH1D* &h1, TH1D* &h2, double& draw_max);

TFile* outFile;
TTree* outTree;

int eID_status;
int eRecon_status;
int mc_PDG;
double EminusPz;
double EoP;

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
PxPyPzEVector vMC_rec;
std::vector<PxPyPzEVector> vMC_hfs;
std::vector<PxPyPzEVector> vREC_hfs;

TH1D* h_nTPts_e;
TH1D* h_nTPts_jet_e;
TH1D* h_nTPts_pi;
TH1D* h_nTPts_else;

TH1D* h_EoP_e;
TH1D* h_EoP_jet_e;
TH1D* h_EoP_pi;
TH1D* h_EoP_else;

TH1D* h_EoP_gapB;
TH1D* h_EoP_gapF;
TH1D* h_EoP_gapF_mod;
TH1D* h_EoP_gapB_mod;

TH1D* h_isoE_e;
TH1D* h_isoE_jet_e;
TH1D* h_isoE_pi;
TH1D* h_isoE_else;

TH1D* h_EoEH_e;
TH1D* h_EoEH_jet_e;
TH1D* h_EoEH_pi;
TH1D* h_EoEH_else;

TH1D* h_EoEH_gapB;
TH1D* h_EoEH_gapF;
TH1D* h_EoEH_gapF_mod;
TH1D* h_EoEH_gapB_mod;

TH1D* h_PIDe_e;
TH1D* h_PIDe_jet_e;
TH1D* h_PIDe_pi;
TH1D* h_PIDe_else;

TH1D* h_PIDh_e;
TH1D* h_PIDh_jet_e;
TH1D* h_PIDh_pi;
TH1D* h_PIDh_else;

TEfficiency* h_pID_pur; // if identify, how many are correct
TEfficiency* h_pID_fal; // if identify, how many are scatter e
TEfficiency* h_pID_eff; // efficiency in using it as veto for scattered electron

TH1D* h_TrackEminusPz;
TH1D* h_CalEminusPz;

TH1D* h_n_scat_elec;
TH2D* h_n_clusters_n_tracks;

TH1D* h_cand_mul;
TH1D* h_cand_mul_eHighPt;
TH1D* h_cand_mul_oHighPt;

TH1D* h_n_cluster_in_cone;
TH1D* h_n_cluster_in_cone_found;

TH2D* h_pt_theta;

#endif