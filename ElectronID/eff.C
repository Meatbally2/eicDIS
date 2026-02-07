#include "../GlobalUtil/luminosityTable.h"
#include "../GlobalUtil/DrawManager.cc"
#include "../GlobalUtil/Constants.hh"
#include "../GlobalUtil/drawHelper.h"

#include "draw_angle.C"

const std::string group[4] = {"minQ2=1", "minQ2=10", "minQ2=100", "minQ2=1000"};

void eff(int beam_type, int Ee, int Eh)
{
    // ePIC plotting style setup
    std::string type_title[5] = {"e^{3}He", "ep", "#gammap", "ep w. BeamBG", "ep"};
    std::string energy_title = beam_type ? Form("%dx%d GeV", Ee, Eh) : Form("%dx%d GeV/A", Ee, Eh);
    std::string campaign = beam_type ? "25.10.0" : "25.10.2";
    if ( beam_type == 4 )
        campaign = "25.10.4";
    DrawManager* draw_manager = new DrawManager(type_title[beam_type], energy_title, campaign);
    draw_manager->SetEPIC();

    draw_manager->SetDISrange(0.01, 0.95, 4, 2);

    double text_lumi = 10; // fb^-1
    if ( beam_type == 0 && Ee == 10 && Eh == 166 )
        text_lumi = 8.65; // fb^-1
    draw_manager->SetLumi(text_lumi);
     
    // get generated lumi
    
    double total_lumi = beam_type == 0 ? 10*3 : 10; // fb^-1
    std::vector<double> gen_lumi = get_lumi(beam_type, Ee, Eh);
    if ( gen_lumi.empty() )
    {
        std::cout << "** Lumi table not found! " << std::endl;
        return;
    }

    int n_group = beam_type == 0 ? 3 : 4;
    std::string setting = beam_type == 0 ? Form("eHe3_%dx%d", (int)Ee, (int)Eh) : Form("ep_%dx%d", (int)Ee, (int)Eh);

    TH2F* h_xq2_all = BookTH2(Form("h_xq2_all"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    TH2F* h_xq2_acp = BookTH2(Form("h_xq2_acp"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    
    TH2F* h_xq2_eff = BookTH2(Form("h_xq2_eff"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    TH2F* h_xq2_pur = BookTH2(Form("h_xq2_pur"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    
    TH2F* h_xq2_eID = BookTH2(Form("h_xq2_eID"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    TH2F* h_xq2_piID = BookTH2(Form("h_xq2_piID"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);

    for ( int i = 0; i < n_group; i ++ )
    {
        std::string date = beam_type == 0 ? "en_25_10_2" : "ep_25_10_0";
        if ( beam_type == 4 )
            date = "ep_25_10_4";
        TFile* file  = new TFile(Form("../data/%s/root_files/%s_%s_eid_combined.root", date.c_str(), setting.c_str(), group[i].c_str()));
        // TFile* file  = new TFile(Form("../data/%s/%s_%s_eid_combined.root", date.c_str(), setting.c_str(), group[i].c_str()));

        TH2F* h_tmp_all = BookTH2(Form("h_tmp_all_%d", i), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
        TH2F* h_tmp_acp = BookTH2(Form("h_tmp_acp_%d", i), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
        TH2F* h_tmp_eff = BookTH2(Form("h_tmp_eff_%d", i), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
        TH2F* h_tmp_pur = BookTH2(Form("h_tmp_pur_%d", i), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
        TH2F* h_tmp_eID = BookTH2(Form("h_tmp_eID_%d", i), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
        TH2F* h_tmp_piID = BookTH2(Form("h_tmp_piID_%d", i), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);

        TTreeReader reader("T_eID", file);

        TTreeReaderValue<int>    status(reader, "eID_status");

        TTreeReaderValue<double> mc_xB(reader, "mc_xB");
        TTreeReaderValue<double> mc_Q2(reader, "mc_Q2");
        TTreeReaderValue<double> mc_y(reader, "mc_y");
        TTreeReaderValue<double> mc_W2(reader, "mc_W2");
        TTreeReaderValue<double> mc_nu(reader, "mc_nu");

        TTreeReaderValue<double> rec_xB(reader, "rec_xB");
        TTreeReaderValue<double> rec_Q2(reader, "rec_Q2");
        TTreeReaderValue<double> rec_y(reader, "rec_y");
        TTreeReaderValue<double> rec_W2(reader, "rec_W2");
        TTreeReaderValue<double> rec_nu(reader, "rec_nu");

        Long64_t nEntries = reader.GetEntries();

        for( size_t ev = 0; ev < nEntries; ev++ ) 
        {
            reader.Next();

            // if(ev%100==0) 
            // cout << "Analysing file " << i << " event " << ev << "/" << nEntries << "\t\r" << std::flush;

            if ( *status < FOUND_MC )
                continue; 

            if (*mc_y < 0.01 || *mc_y > 0.95) 
                continue;

            if (*mc_W2 < 4)
                continue;

            if (*mc_Q2 < 2)
                continue;

            if ( beam_type == EP ) 
            {
                if ( i == 3 && *mc_Q2 < pow(10,i) ) 
                    continue;
            
                if ( i < 3 && (*mc_Q2 > pow(10,i+1) || *mc_Q2 < pow(10,i)) ) 
                    continue;
            }

            if ( *status >= FOUND_MC )
                h_tmp_all->Fill(*mc_xB, *mc_Q2);

            if ( *status >= FOUND_TRUTH )
                h_tmp_acp->Fill(*mc_xB, *mc_Q2);

            if ( *status >= FOUND_E )
                h_tmp_eff->Fill(*mc_xB, *mc_Q2);   

            if ( *status == FOUND_E )
            {
                h_tmp_pur->Fill(*mc_xB, *mc_Q2);
                h_tmp_eID->Fill(*mc_xB, *mc_Q2);
            }

            if ( *status >= FOUND_PI )
                h_tmp_piID->Fill(*mc_xB, *mc_Q2);
        }

        h_tmp_all->Scale(total_lumi/gen_lumi[i]);
        h_tmp_acp->Scale(total_lumi/gen_lumi[i]);
        h_tmp_eff->Scale(total_lumi/gen_lumi[i]);
        h_tmp_pur->Scale(total_lumi/gen_lumi[i]);
        h_tmp_eID->Scale(total_lumi/gen_lumi[i]);
        h_tmp_piID->Scale(total_lumi/gen_lumi[i]);

        h_xq2_all->Add(h_tmp_all);
        h_xq2_acp->Add(h_tmp_acp);
        h_xq2_eff->Add(h_tmp_eff);
        h_xq2_pur->Add(h_tmp_pur);
        h_xq2_eID->Add(h_tmp_eID);
        h_xq2_piID->Add(h_tmp_piID);

        file->Close();
    }

    set_2d_scale(h_xq2_all);
    TCanvas* c_xq2_all = draw_2d_standard(h_xq2_all, "c_xq2_all", "all events", 700, 600, true, true);
    TCanvas* c_xq2_acp = draw_2d_standard(h_xq2_acp, "c_xq2_acp", "acp events", 700, 600, true, true);
    TCanvas* c_xq2_eID = draw_2d_standard(h_xq2_eID, "c_xq2_eID", "eID events", 700, 600, true, true);
    TCanvas* c_xq2_piID = draw_2d_standard(h_xq2_piID, "c_xq2_piID", "piID events", 700, 600, true, true);

    TH2F* h_xq2_acp_copy = (TH2F*)h_xq2_acp->Clone();
    process_eff_hist(h_xq2_acp_copy, h_xq2_all);
    TCanvas* c_xq2_acp_eff = draw_2d_efficiency(h_xq2_acp_copy, "c_xq2_acp_eff", "xq2 acp eff", 1400, 600, false, true);

    TH2F* h_xq2_eff_copy = (TH2F*)h_xq2_eff->Clone();
    process_eff_hist(h_xq2_eff_copy, h_xq2_all);
    TCanvas* c_xq2_eff_eff = draw_2d_efficiency(h_xq2_eff_copy, "c_xq2_eff_eff", "xq2 eff eff", 1400, 600, false, true);
    draw_angle(c_xq2_eff_eff, h_xq2_eff_copy, 160.0, Ee*1.0, Eh*1.0);
    draw_angle(c_xq2_eff_eff, h_xq2_eff_copy, 130.0, Ee*1.0, Eh*1.0);
    draw_angle(c_xq2_eff_eff, h_xq2_eff_copy, 40.0, Ee*1.0, Eh*1.0);

    TH2F* h_xq2_pur_copy = (TH2F*)h_xq2_pur->Clone();
    process_eff_hist(h_xq2_pur_copy, h_xq2_eff);
    TCanvas* c_xq2_pur_eff = draw_2d_efficiency(h_xq2_pur_copy, "c_xq2_pur_eff", "xq2 pur eff", 1400, 600, false, true);

    TH2F* h_xq2_eID_copy = (TH2F*)h_xq2_eID->Clone();
    process_eff_hist(h_xq2_eID_copy, h_xq2_acp);
    TCanvas* c_xq2_eID_eff = draw_2d_efficiency(h_xq2_eID_copy, "c_xq2_eID_eff", "xq2 eID eff", 1400, 600, false, true);

    TH2F* h_xq2_piID_copy = (TH2F*)h_xq2_piID->Clone();
    process_eff_hist(h_xq2_piID_copy, h_xq2_acp);
    TCanvas* c_xq2_piID_eff = draw_2d_efficiency(h_xq2_piID_copy, "c_xq2_piID_eff", "xq2 piID eff", 1400, 600, false, true);

    draw_manager->LableAndCollect(c_xq2_all);
    draw_manager->LableAndCollect(c_xq2_acp);
    draw_manager->LableAndCollect(c_xq2_eID);
    draw_manager->LableAndCollect(c_xq2_piID);

    draw_manager->LableAndCollect(c_xq2_acp_eff);
    draw_manager->LableAndCollect(c_xq2_eff_eff);
    draw_manager->LableAndCollect(c_xq2_pur_eff);
    draw_manager->LableAndCollect(c_xq2_eID_eff);
    draw_manager->LableAndCollect(c_xq2_piID_eff);

    // gStyle->SetImageScaling(3.);
    c_xq2_acp_eff->SaveAs(Form("../data/eID/%s_eID_acceptance.png", setting.c_str()));
    c_xq2_eff_eff->SaveAs(Form("../data/eID/%s_eID_efficiency.png", setting.c_str()));
    c_xq2_pur_eff->SaveAs(Form("../data/eID/%s_eID_purity.png", setting.c_str()));

    // c_xq2_all->SaveAs("../data/eID/10x166_en_raw.png");

    return;
}