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
void CreateOutputTree(TString outFileName);
void ResetVariables();

void CalculateElectronKinematics(double fEe, double fEh, TLorentzVector kf, double& xB, double& Q2, double& W2, double& y, double& nu);
TLorentzVector GetHadronBeam(double fEh);

void DefineHistograms();
void SetupPurityEff(TEfficiency* &h, const char* name, const char* xAxisTitle, int nBins, double xMin, double xMax, int color, int marker);
void DrawPurityCanvas(TCanvas* &c, const char* canvasName,
	TEfficiency* h_base,
	TEfficiency* h_veto, TEfficiency* h_e, TEfficiency* h_pi, TEfficiency* h_K, TEfficiency* h_p,
	DrawManager* draw_manager);
void DrawVerticalLine(TCanvas* &c, double x_pos, double y_max);
void DrawParComparison(TCanvas* &c, TH1D* &h1, TH1D* &h2, TH1D* &h3, TH1D* &h4, double& draw_max);
void DrawTCComparison(TCanvas* &c, TH1D* &h1, TH1D* &h2, double& draw_max);

void FillEidPurity(int reco_pid, int mc_pdg, double eta, double momentum);
void FillNegTrackPurity( int reco_pid, int mc_pdg, double eta, double momentum);
void FillPidPurity( int reco_pid, int par_type, double eta, double momentum);
void FillPidSuccess( int reco_pid, int par_type);

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
TEfficiency* h_pID_suc; // In general, how many are correct
TEfficiency* h_pID_eff; // efficiency in using it as veto for scattered electron

TEfficiency* h_e_eID_pur_eta;
TEfficiency* h_pi_eID_pur_eta;
TEfficiency* h_K_eID_pur_eta;
TEfficiency* h_p_eID_pur_eta;
TEfficiency* h_eVeto_eID_pur_eta;

TEfficiency* h_e_pID_pur_eta;
TEfficiency* h_pi_pID_pur_eta;
TEfficiency* h_K_pID_pur_eta;
TEfficiency* h_p_pID_pur_eta;
TEfficiency* h_eVeto_pID_pur_eta;

TEfficiency* h_e_eID_pur_p;
TEfficiency* h_pi_eID_pur_p;
TEfficiency* h_K_eID_pur_p;
TEfficiency* h_p_eID_pur_p;
TEfficiency* h_eVeto_eID_pur_p;

TEfficiency* h_e_trk_pur_eta;
TEfficiency* h_pi_trk_pur_eta;
TEfficiency* h_K_trk_pur_eta;
TEfficiency* h_p_trk_pur_eta;
TEfficiency* h_eVeto_trk_pur_eta;

TEfficiency* h_e_trk_pur_p;
TEfficiency* h_pi_trk_pur_p;
TEfficiency* h_K_trk_pur_p;
TEfficiency* h_p_trk_pur_p;
TEfficiency* h_eVeto_trk_pur_p;

TEfficiency* h_e_pID_pur_p;
TEfficiency* h_pi_pID_pur_p;
TEfficiency* h_K_pID_pur_p;
TEfficiency* h_p_pID_pur_p;
TEfficiency* h_eVeto_pID_pur_p;

TEfficiency* h_cur_pur_eta;
TEfficiency* h_cur_pur_p;

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