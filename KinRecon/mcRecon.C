// only reconstruct real kinematics

#include <Math/LorentzVector.h>
using ROOT::Math::PxPyPzEVector;

#include <Math/LorentzRotation.h>
using ROOT::Math::LorentzRotation;

#include "../GlobalUtil/DrawManager.cc"
#include "../GlobalUtil/getBoost.h"
#include "reconMethod.hh"

const std::string group[4] = {"minQ2=1", "minQ2=10", "minQ2=100", "minQ2=1000"};

void sumRuns(int beam_type, int Ee, int Eh)
{
    std::string setting = beam_type == 0 ? Form("eHe3_%dx%d", (int)Ee, (int)Eh) : Form("ep_%dx%d", (int)Ee, (int)Eh);
    int n_group = beam_type == 0 ? 3 : 4;

     // ePIC plotting style setup
    std::string type_title[5] = {"e^{3}He", "ep", "#gammap", "ep w. BeamBG", "ep"};
    std::string energy_title = beam_type ? Form("%dx%d GeV", Ee, Eh) : Form("%dx%d GeV/A", Ee, Eh);
    std::string campaign = beam_type ? "25.10.0" : "25.10.2";
    if ( beam_type == 4 )
        campaign = "25.10.4";
    DrawManager* draw_manager = new DrawManager(type_title[beam_type], energy_title, campaign);
    draw_manager->SetEPIC();

    double text_lumi = 10; // fb^-1
    if ( beam_type == 0 && Ee == 10 && Eh == 166 )
        text_lumi = 8.65; // fb^-1
    if ( beam_type == 4 && Ee == 10 && Eh == 130 )
        text_lumi = 5.33; // fb^-1
    if ( beam_type == 4 && Ee == 10 && Eh == 250 )
        text_lumi = 9.18; // fb^-1

    // get generated lumi
    double total_lumi = beam_type == 0 ? text_lumi*3 : text_lumi; // fb^-1
    std::vector<double> gen_lumi = get_lumi(beam_type, Ee, Eh);
    if ( gen_lumi.empty() )
    {
        std::cout << "** Lumi table not found! " << std::endl;
        return;
    }

    for ( int i = 0; i < n_group; i ++ )
    {
        std::string date = beam_type == 0 ? "en_25_10_2" : "ep_25_10_0";
        if ( beam_type == 4 )
            date = "ep_25_10_4";
        TFile* beamFile = new TFile(Form("../data/%s/root_files/%s_%s_beam_combined.root", date.c_str(), setting.c_str(), group[i].c_str()));
        TFile* tagFile  = new TFile(Form("../data/%s/root_files/%s_%s_tag_combined.root", date.c_str(), setting.c_str(), group[i].c_str()));
        TFile* eidFile  = new TFile(Form("../data/%s/root_files/%s_%s_eid_combined.root", date.c_str(), setting.c_str() , group[i].c_str()));

        TTreeReader beam_reader("T_Beam", beamFile);
        TTreeReaderValue<int> N_PDG(beam_reader, "N_PDG");
        TTreeReaderValue<PxPyPzEVector> beam_e(beam_reader, "vectE");
        TTreeReaderValue<PxPyPzEVector> beam_n(beam_reader, "vectN");

        TTreeReader tag_reader("T_Tag", tagFile);
        TTreeReaderValue<bool> fTagged(tag_reader, "fTagged");
        TTreeReaderValue<int> nTracks(tag_reader, "nTracks");
        TTreeReaderValue<bool> fZDCn(tag_reader, "fZDCn");

        TTreeReader eid_reader("T_eID", eidFile);
        TTreeReaderValue<int> eID_status(eid_reader, "eID_status");
        TTreeReaderValue<double> EminusPz(eid_reader, "EminusPz");
        TTreeReaderValue<PxPyPzEVector> vMCe(eid_reader, "vMC_e");
        TTreeReaderValue<PxPyPzEVector> vTRe(eid_reader, "vTRACK_e");
        TTreeReaderValue<PxPyPzEVector> vCLe(eid_reader, "vCLUSTER_e");
        TTreeReaderValue<std::vector<PxPyPzEVector>> vMCf(eid_reader, "vMC_hfs");
        TTreeReaderValue<std::vector<PxPyPzEVector>> vTRf(eid_reader, "vREC_hfs");
    }

    //
    return;
}