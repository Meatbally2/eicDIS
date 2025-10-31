#ifndef EIDANA_H
#define EIDANA_H

#include "podio/Frame.h"
#include "podio/ROOTReader.h"

#include "../GlobalUtil/AnaManager.cc"
#include "../GlobalUtil/Constants.hh"

#include "ElectronID.cc"

void eIDana(int Ee, int Eh, int select_region, int sr, int is_truth_eID, int all_file, int analyse_p);
void CreateOutputTree(TString outFileName);
void ResetVariables();

void CalculateElectronKinematics(double fEe, double fEh, TLorentzVector kf, double& xB, double& Q2, double& W2, double& y, double& nu);
TLorentzVector GetHadronBeam(double fEh);

TFile* outFile;
TTree* outTree;

int eID_status;
int mc_PBG;

double mc_xB;
double mc_Q2;
double mc_y;
double mc_W2;
double mc_nu;

enum { NO_FOUND, FOUND_E, FOUND_PI, FOUND_OTHERS };

#endif