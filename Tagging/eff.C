#include "../GlobalUtil/DrawManager.cc"
#include "../GlobalUtil/Constants.hh"
#include "../GlobalUtil/drawHelper.h"
#include <Math/LorentzVector.h>
using ROOT::Math::PxPyPzEVector;

const std::string group[3] = {"1to10", "10to100", "100to1000"};
const double total_lumi = 8.65*3; // fb^-1
int e_set[3] = {1, 10, 100};
int n_group = 3;

double cs[3][2] = {{0.198440424611563, 0.205327493968226}, {4.04371412707044E-02, 4.41976212963417E-02}, {1.36416909784756E-03, 1.69583242740138E-03}};
double ev[3][2] = {{333675, 666325}, {333694, 666306}, {333365, 666640}};

void eff()
{
    std::string type_title = "e^{3}He";
    std::string energy_title = Form("%dx%d GeV", 10, 166);
    DrawManager* draw_manager = new DrawManager(type_title, energy_title, "25.10.2");
    draw_manager->SetEPIC();
    // draw_manager->SetDISrange(0.01, 0.95, 4, 2);
    draw_manager->SetLumi(8.65);

    TH2F* h_xq2_all = BookTH2(Form("h_xq2_all"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    TH2F* h_xq2_eff = BookTH2(Form("h_xq2_eff"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    TH2F* h_xq2_pur = BookTH2(Form("h_xq2_pur"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    
    for ( int i = 0; i < n_group; i ++ )
    {
        TFile* beamFile = new TFile(Form("../data/en_25_10_2/root_files/10x166_%s_beam_combined.root", group[i].c_str()));
        TFile* tagFile  = new TFile(Form("../data/en_25_10_2/root_files/10x166_%s_tag_combined.root", group[i].c_str()));
        TFile* eidFile  = new TFile(Form("../data/en_25_10_2/root_files/10x166_%s_eIDrecon_combined.root", group[i].c_str()));
        
        TTreeReader beam_reader("T_Beam", beamFile);
        TTreeReaderValue<int> N_PDG(beam_reader, "N_PDG");
        TTreeReaderValue<PxPyPzEVector> beam_e(beam_reader, "vectE");
        TTreeReaderValue<PxPyPzEVector> beam_n(beam_reader, "vectN");

        TTreeReader tag_reader("T_Tag", tagFile);
        TTreeReaderValue<bool> fTagged(tag_reader, "fTagged");

        TTreeReader eid_reader("T_eID", eidFile);
        TTreeReaderValue<int>    status(eid_reader, "eID_status");
        TTreeReaderValue<double> mc_xB(eid_reader, "mc_xB");
        TTreeReaderValue<double> mc_Q2(eid_reader, "mc_Q2");
        TTreeReaderValue<double> mc_y(eid_reader, "mc_y");
        TTreeReaderValue<double> mc_W2(eid_reader, "mc_W2");
        TTreeReaderValue<double> mc_nu(eid_reader, "mc_nu");

        Long64_t nEntries = beam_reader.GetEntries();

        double gen_lumi = ev[i][0]/(cs[i][0]*(1e-34/1e-43)) + ev[i][1]/(cs[i][1]*(1e-34/1e-43));

        TH2F* h_tmp_all = BookTH2(Form("h_tmp_all_%i", i), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
        TH2F* h_tmp_eff = BookTH2(Form("h_tmp_eff_%i", i), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
        TH2F* h_tmp_pur = BookTH2(Form("h_tmp_pur_%i", i), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);

        for ( int ev = 0; ev < nEntries; ev ++ )
        {
            beam_reader.Next();
            tag_reader.Next();
            eid_reader.Next();

            if ( *status == NO_MC )
                continue;

            if (*mc_y < 0.01 || *mc_y > 0.95) 
                continue;

            if (*mc_W2 < 4)
                continue;

            if (*mc_Q2 < 2)
                continue;

            if ( *N_PDG == ID_NEUTRON )
            {
                h_tmp_all->Fill(*mc_xB, *mc_Q2);

                if ( *fTagged == 1 )
                    h_tmp_pur->Fill(*mc_xB, *mc_Q2);
            }

            if ( *fTagged == 1 )
                h_tmp_eff->Fill(*mc_xB, *mc_Q2);
        }

        h_tmp_all->Scale(total_lumi/gen_lumi);
        h_tmp_eff->Scale(total_lumi/gen_lumi);
        h_tmp_pur->Scale(total_lumi/gen_lumi);

        h_xq2_all->Add(h_tmp_all);
        h_xq2_eff->Add(h_tmp_eff);
        h_xq2_pur->Add(h_tmp_pur);

        beamFile->Close();
        tagFile->Close();
        eidFile->Close();
    }

    set_2d_scale(h_xq2_all);
    TCanvas* c_xq2_all = draw_2d_standard(h_xq2_all, "c_xq2_all", "all n events", 700, 600, true, true);
    TCanvas* c_xq2_eff = draw_2d_standard(h_xq2_eff, "c_xq2_eff", "eff events", 700, 600, true, true);
    TCanvas* c_xq2_pur = draw_2d_standard(h_xq2_pur, "c_xq2_pur", "purity events", 700, 600, true, true);

    TH2F* h_xq2_eff_copy = (TH2F*)h_xq2_eff->Clone();
    process_eff_hist(h_xq2_eff_copy, h_xq2_all);
    TCanvas* c_xq2_ff_eff = draw_2d_efficiency(h_xq2_eff_copy, "c_ff_eff", "ff eff", 1400, 600, false, true);

    TH2F* h_xq2_pur_copy = (TH2F*)h_xq2_pur->Clone();
    process_eff_hist(h_xq2_pur_copy, h_xq2_eff);
    TCanvas* c_xq2_ff_pur = draw_2d_efficiency(h_xq2_pur_copy, "c_ff_pur", "ff pur", 1400, 600, false, true);

    draw_manager->LableAndCollect(c_xq2_all);
    draw_manager->LableAndCollect(c_xq2_eff);
    draw_manager->LableAndCollect(c_xq2_pur);
    draw_manager->LableAndCollect(c_xq2_ff_eff);
    draw_manager->LableAndCollect(c_xq2_ff_pur);   

    c_xq2_ff_eff->SaveAs("../data/eID/ff_eff_10x166_e3He_25.10.2.png");
    c_xq2_ff_pur->SaveAs("../data/eID/ff_pur_10x166_e3He_25.10.2.png");
}