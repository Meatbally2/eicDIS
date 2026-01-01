#include "../GlobalUtil/DrawManager.cc"
#include "../GlobalUtil/Constants.hh"
#include "../GlobalUtil/drawHelper.h"
#include "../GlobalUtil/luminosityTable.h"
#include <Math/LorentzVector.h>
using ROOT::Math::PxPyPzEVector;

int n_group = 3;
const double total_lumi = 8.65*3; // fb^-1
const std::string group[4] = {"1", "10", "100", "1000"};

double cs[3][2] = {{0.198440424611563, 0.205327493968226}, {4.04371412707044E-02, 4.41976212963417E-02}, {1.36416909784756E-03, 1.69583242740138E-03}};
double ev[3][2] = {{333675, 666325}, {333694, 666306}, {333365, 666640}};

void eff(int Ee, int Eh)
{
     // ePIC plotting style setup
    std::string type_title = "e^{3}He";
    std::string energy_title = Form("%dx%d GeV", Ee, Eh);
    DrawManager* draw_manager = new DrawManager(type_title, energy_title, "25.10.2");
    draw_manager->SetEPIC();

    double text_lumi = 10; // fb^-1
    if ( Ee == 10 && Eh == 166 )
        text_lumi = 8.65; // fb^-1
    draw_manager->SetLumi(text_lumi);

    // get generated lumi
    double total_lumi = text_lumi*3; // fb^-1
    std::vector<double> gen_lumi = get_lumi(0, Ee, Eh);
    if ( gen_lumi.empty() )
    {
        std::cout << "** Lumi table not found! " << std::endl;
        return;
    }

    // histograms
    TH2F* h_xq2_MCall = BookTH2(Form("h_xq2_MCall"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    TH2F* h_xq2_MCen = BookTH2(Form("h_xq2_MCen"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    TH2F* h_xq2_DTeff = BookTH2(Form("h_xq2_DTeff"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    TH2F* h_xq2_DTpur = BookTH2(Form("h_xq2_DTpur"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    TH2F* h_xq2_ZDCeff = BookTH2(Form("h_xq2_ZDCeff"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    TH2F* h_xq2_ZDCpur = BookTH2(Form("h_xq2_ZDCpur"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    
    for ( int i = 0; i < n_group; i ++ )
    {
        TFile* beamFile = new TFile(Form("../data/en_25_10_2/root_files/eHe3_%dx%d_minQ2=%s_beam_combined.root", Ee, Eh, group[i].c_str()));
        TFile* tagFile  = new TFile(Form("../data/en_25_10_2/root_files/eHe3_%dx%d_minQ2=%s_tag_combined.root", Ee, Eh, group[i].c_str()));
        TFile* eidFile  = new TFile(Form("../data/en_25_10_2/root_files/eHe3_%dx%d_minQ2=%s_eid_combined.root", Ee, Eh, group[i].c_str()));
        
        TTreeReader beam_reader("T_Beam", beamFile);
        TTreeReaderValue<int> N_PDG(beam_reader, "N_PDG");
        TTreeReaderValue<PxPyPzEVector> beam_e(beam_reader, "vectE");
        TTreeReaderValue<PxPyPzEVector> beam_n(beam_reader, "vectN");

        TTreeReader tag_reader("T_Tag", tagFile);
        TTreeReaderValue<bool> fTagged(tag_reader, "fTagged");
        TTreeReaderValue<int> nTracks(tag_reader, "nTracks");
        TTreeReaderValue<bool> fZDCn(tag_reader, "fZDCn");

        TTreeReader eid_reader("T_eID", eidFile);
        TTreeReaderValue<int>    status(eid_reader, "eID_status");
        TTreeReaderValue<double> mc_xB(eid_reader, "mc_xB");
        TTreeReaderValue<double> mc_Q2(eid_reader, "mc_Q2");
        TTreeReaderValue<double> mc_y(eid_reader, "mc_y");
        TTreeReaderValue<double> mc_W2(eid_reader, "mc_W2");
        TTreeReaderValue<double> mc_nu(eid_reader, "mc_nu");

        Long64_t nEntries = beam_reader.GetEntries();

        double gen_lumi = ev[i][0]/(cs[i][0]*(1e-34/1e-43)) + ev[i][1]/(cs[i][1]*(1e-34/1e-43));

        TH2F* h_tmp_MCall = BookTH2(Form("h_tmp_MCall_%i", i), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
        TH2F* h_tmp_MCen = BookTH2(Form("h_tmp_MCen_%i", i), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
        TH2F* h_tmp_DTeff = BookTH2(Form("h_tmp_DTeff_%i", i), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
        TH2F* h_tmp_DTpur = BookTH2(Form("h_tmp_DTpur_%i", i), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
        TH2F* h_tmp_ZDCeff = BookTH2(Form("h_tmp_ZDCeff_%i", i), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
        TH2F* h_tmp_ZDCpur = BookTH2(Form("h_tmp_ZDCpur_%i", i), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);

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

            h_tmp_MCall->Fill(*mc_xB, *mc_Q2);

            if ( *N_PDG == ID_NEUTRON )
            {
                h_tmp_MCen->Fill(*mc_xB, *mc_Q2);

                // if ( *fTagged == 1 )
                if (*nTracks >= 2)
                {
                    h_tmp_DTpur->Fill(*mc_xB, *mc_Q2);

                    if ( *fZDCn == 1 )
                        h_tmp_ZDCpur->Fill(*mc_xB, *mc_Q2);
                }
            }

            // if ( *fTagged == 1 )
            if (*nTracks >= 2)
            {
                h_tmp_DTeff->Fill(*mc_xB, *mc_Q2);

                if ( *fZDCn == 1 )
                    h_tmp_ZDCeff->Fill(*mc_xB, *mc_Q2);
            }    
        }

        h_tmp_MCall->Scale(total_lumi/gen_lumi);
        h_tmp_MCen->Scale(total_lumi/gen_lumi);
        h_tmp_DTeff->Scale(total_lumi/gen_lumi);
        h_tmp_DTpur->Scale(total_lumi/gen_lumi);
        h_tmp_ZDCeff->Scale(total_lumi/gen_lumi);
        h_tmp_ZDCpur->Scale(total_lumi/gen_lumi);

        h_xq2_MCall->Add(h_tmp_MCall);
        h_xq2_MCen->Add(h_tmp_MCen);
        h_xq2_DTeff->Add(h_tmp_DTeff);
        h_xq2_DTpur->Add(h_tmp_DTpur);
        h_xq2_ZDCeff->Add(h_tmp_ZDCeff);
        h_xq2_ZDCpur->Add(h_tmp_ZDCpur);

        beamFile->Close();
        tagFile->Close();
        eidFile->Close();
    }

    set_2d_scale(h_xq2_MCall);
    TCanvas* c_xq2_MCall = draw_2d_standard(h_xq2_MCall, "c_xq2_MCall", "all MC events", 700, 600, true, true);
    TCanvas* c_xq2_MCen = draw_2d_standard(h_xq2_MCen, "c_xq2_MCen", "all n events", 700, 600, true, true);
    TCanvas* c_xq2_DTeff = draw_2d_standard(h_xq2_DTeff, "c_xq2_DTeff", "tagged eff events", 700, 600, true, true);
    TCanvas* c_xq2_DTpur = draw_2d_standard(h_xq2_DTpur, "c_xq2_DTpur", "tagged purity events", 700, 600, true, true);
    TCanvas* c_xq2_ZDCeff = draw_2d_standard(h_xq2_ZDCeff, "c_xq2_ZDCeff", "ZDC eff events", 700, 600, true, true);
    TCanvas* c_xq2_ZDCpur = draw_2d_standard(h_xq2_ZDCpur, "c_xq2_ZDCpur", "ZDC purity events", 700, 600, true, true);

    TH2F* h_xq2_MCep = (TH2F*)h_xq2_MCall->Clone();
    h_xq2_MCep->Add(h_xq2_MCen,-1);

    TH2F* h_xq2_MCen_copy = (TH2F*)h_xq2_MCen->Clone();
    process_eff_hist(h_xq2_MCen_copy, h_xq2_MCall);
    TCanvas* c_xq2_fen = draw_2d_efficiency(h_xq2_MCen_copy, "c_fen", "frac en", 1400, 600, false, true);

    TH2F* h_xq2_DTeff_copy = (TH2F*)h_xq2_DTeff->Clone();
    process_eff_hist(h_xq2_DTeff_copy, h_xq2_MCen);
    TCanvas* c_xq2_ff_eff = draw_2d_efficiency(h_xq2_DTeff_copy, "c_ff_eff", "ff eff", 1400, 600, false, true);

    TH2F* h_xq2_DTpur_copy = (TH2F*)h_xq2_DTpur->Clone();
    process_eff_hist(h_xq2_DTpur_copy, h_xq2_DTeff);
    TCanvas* c_xq2_ff_pur = draw_2d_efficiency(h_xq2_DTpur_copy, "c_ff_pur", "ff pur", 1400, 600, false, true);

    TH2F* h_xq2_ZDCeff_copy = (TH2F*)h_xq2_ZDCeff->Clone();
    process_eff_hist(h_xq2_ZDCeff_copy, h_xq2_DTeff);
    TCanvas* c_xq2_ZDC_eff = draw_2d_efficiency(h_xq2_ZDCeff_copy, "c_ZDC_eff", "ZDC eff", 1400, 600, false, true);

    TH2F* h_xq2_ZDCpur_copy = (TH2F*)h_xq2_ZDCpur->Clone();
    process_eff_hist(h_xq2_ZDCpur_copy, h_xq2_ZDCeff);
    TCanvas* c_xq2_ZDC_pur = draw_2d_efficiency(h_xq2_ZDCpur_copy, "c_ZDC_pur", "ZDC pur", 1400, 600, false, true);

    draw_manager->LableAndCollect(c_xq2_MCall);
    draw_manager->LableAndCollect(c_xq2_MCen);
    draw_manager->LableAndCollect(c_xq2_fen);
    draw_manager->LableAndCollect(c_xq2_DTeff);
    draw_manager->LableAndCollect(c_xq2_DTpur);
    draw_manager->LableAndCollect(c_xq2_ZDCeff);
    draw_manager->LableAndCollect(c_xq2_ZDCpur);
    draw_manager->LableAndCollect(c_xq2_ff_eff);
    draw_manager->LableAndCollect(c_xq2_ff_pur);
    draw_manager->LableAndCollect(c_xq2_ZDC_eff);
    draw_manager->LableAndCollect(c_xq2_ZDC_pur);   

    c_xq2_ff_eff->SaveAs(Form("../data/eID/ff_eff_%dx%d_e3He_25.10.2.png", Ee, Eh));
    c_xq2_ff_pur->SaveAs(Form("../data/eID/ff_pur_%dx%d_e3He_25.10.2.png", Ee, Eh));
}