
#include "../GlobalUtil/DrawManager.cc"
#include "../GlobalUtil/Constants.hh"
#include "../GlobalUtil/drawHelper.h"
#include "../GlobalUtil/luminosityTable.h"

#include "draw_angle.C"

const std::string group[4] = {"minQ2=1", "minQ2=10", "minQ2=100", "minQ2=1000"};

void eff_pi_bg(int beam_type, int Ee, int Eh)
{
    // ePIC plotting style setup
    std::string type_title[3] = {"e^{3}He", "ep", "#gammap"};
    std::string energy_title = beam_type ? Form("%dx%d GeV", Ee, Eh) : Form("%dx%d GeV/A", Ee, Eh);
    std::string campaign = beam_type ? "25.10.0" : "25.10.2";
    DrawManager* draw_manager = new DrawManager(type_title[beam_type], energy_title, campaign);
    draw_manager->SetEPIC();

    double text_lumi = 10; // fb^-1
    if ( beam_type == 0 && Ee == 10 && Eh == 166 )
        text_lumi = 8.65; // fb^-1
    draw_manager->SetLumi(text_lumi);
    draw_manager->SetDISrange(0.01, 0.95, 4, 2);
     
    // get generated lumi
    
    double total_lumi = beam_type == 0 ? 10*3 : 10; // fb^-1
    std::vector<double> gen_lumi_ep = get_lumi(beam_type, Ee, Eh);
    std::vector<double> gen_lumi_pi = get_lumi(PI_BG, Ee, Eh);
    if ( gen_lumi_ep.empty() || gen_lumi_pi.empty() )
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
    TH1F* h_EminusPz = new TH1F("h_EminusPz", "E-P_{z} all;E-P_{z} (GeV);Counts", 100, 0, 50);

    for ( int i = 0; i < n_group; i ++ )
    {
        std::string date = beam_type == 0 ? "en_25_10_2" : "ep_25_10_0";
        TFile* file = new TFile(Form("../data/%s/root_files/%s_%s_eID_combined.root", date.c_str(), setting.c_str(), group[i].c_str()));
        double gen_lumi = gen_lumi_ep[i];

        TH2F* h_tmp_all = BookTH2(Form("h_tmp_all_%d", i), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
        TH2F* h_tmp_acp = BookTH2(Form("h_tmp_acp_%d", i), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
        TH2F* h_tmp_eff = BookTH2(Form("h_tmp_eff_%d", i), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
        TH2F* h_tmp_pur = BookTH2(Form("h_tmp_pur_%d", i), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
        TH2F* h_tmp_eID = BookTH2(Form("h_tmp_eID_%d", i), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
        TH2F* h_tmp_piID = BookTH2(Form("h_tmp_piID_%d", i), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
        TH1F* h_tmp_EminusPz = new TH1F(Form("h_tmp_EminusPz_%d", i), "E-P_{z} all;E-P_{z} (GeV);Counts", 100, 0, 50);

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

        TTreeReaderValue<double> EminusPz(reader, "EminusPz");

        Long64_t nEntries = reader.GetEntries();

        for( size_t ev = 0; ev < nEntries; ev++ ) 
        {
            reader.Next();

            // if(ev%100==0) 
            // cout << "Analysing file " << i << " event " << ev << "/" << nEntries << "\t\r" << std::flush;

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

            if ( *EminusPz < 25 )
                continue;
               
            if ( *status >= FOUND_MC )
                h_tmp_all->Fill(*mc_xB, *mc_Q2);

            if ( *status >= FOUND_TRUTH )
                h_tmp_acp->Fill(*mc_xB, *mc_Q2);

            if ( *status >= FOUND_E )
            {
                h_tmp_eff->Fill(*mc_xB, *mc_Q2);   
                h_tmp_EminusPz->Fill(*EminusPz);
            }

            if ( *status == FOUND_E )
            {
                h_tmp_pur->Fill(*mc_xB, *mc_Q2);
                h_tmp_eID->Fill(*mc_xB, *mc_Q2);
            }

            if ( *status >= FOUND_PI )
                h_tmp_piID->Fill(*mc_xB, *mc_Q2);
        }

        h_tmp_all->Scale(total_lumi/gen_lumi);
        h_tmp_acp->Scale(total_lumi/gen_lumi);
        h_tmp_eff->Scale(total_lumi/gen_lumi);
        h_tmp_pur->Scale(total_lumi/gen_lumi);
        h_tmp_eID->Scale(total_lumi/gen_lumi);
        h_tmp_piID->Scale(total_lumi/gen_lumi);
        h_tmp_EminusPz->Scale(total_lumi/gen_lumi);

        h_xq2_all->Add(h_tmp_all);
        h_xq2_acp->Add(h_tmp_acp);
        h_xq2_eff->Add(h_tmp_eff);
        h_xq2_pur->Add(h_tmp_pur);
        h_xq2_eID->Add(h_tmp_eID);
        h_xq2_piID->Add(h_tmp_piID);
        h_EminusPz->Add(h_tmp_EminusPz);

        file->Close();
    }

    std::cout << "Finished processing signal files." << std::endl;

    TFile* pi_bg_file = new TFile(Form("../data/ep_25_10_0/root_files/piBG_18x275_eid_combined.root"));
    if ( pi_bg_file == nullptr || pi_bg_file->IsZombie() ) {
        std::cerr << "Error: Could not open pi bg file!" << std::endl;
        return;
    }

    TH2F* h_xq2_all_bg = BookTH2(Form("h_xq2_all_bg"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    TH2F* h_xq2_acp_bg = BookTH2(Form("h_xq2_acp_bg"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    TH2F* h_xq2_eff_bg = BookTH2(Form("h_xq2_eff_bg"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    TH2F* h_xq2_pur_bg = BookTH2(Form("h_xq2_pur_bg"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    TH2F* h_xq2_eID_bg = BookTH2(Form("h_xq2_eID_bg"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    TH2F* h_xq2_piID_bg = BookTH2(Form("h_xq2_piID_bg"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    TH1F* h_EminusPz_bg = new TH1F("h_EminusPz_bg", "E-P_{z} all;E-P_{z} (GeV);Counts", 100, 0, 50);

    TTreeReader reader("T_eID", pi_bg_file);

    TTreeReaderValue<int> status(reader, "eID_status");

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

    TTreeReaderValue<double> EminusPz(reader, "EminusPz");

    Long64_t nEntries = reader.GetEntries();

    std::cout << "Processing pi bg file with " << nEntries << " entries." << std::endl;

    for( size_t ev = 0; ev < nEntries; ev++ ) 
    {
        reader.Next();

        // std::cout << *mc_xB << " " << *mc_Q2 << " " ;
        // std::cout << *rec_xB << " " << *rec_Q2 << " ";
        // std::cout << *status << std::endl;

        if ( *rec_y < 0.01 || *rec_y > 0.95 )
            continue;

        if (*rec_W2 < 4)
            continue;

        if (*rec_Q2 < 2)
            continue;
        
        if ( *EminusPz < 25 )
            continue;

        if ( *status >= FOUND_MC )
            h_xq2_all_bg->Fill(*mc_xB, *mc_Q2);

        if ( *status >= FOUND_TRUTH )
            h_xq2_acp_bg->Fill(*mc_xB, *mc_Q2);

        if ( *status >= FOUND_E )
        {
            h_xq2_eff_bg->Fill(*rec_xB, *rec_Q2);   
        }

        if ( *status == FOUND_E )
        {
            // std::cout << "Filling pi bg event: rec_xB = " << *rec_xB << ", rec_Q2 = " << *rec_Q2 << std::endl;
            h_xq2_pur_bg->Fill(*rec_xB, *rec_Q2);
            h_xq2_eID_bg->Fill(*rec_xB, *rec_Q2);
        }

        if ( *status >= FOUND_PI )
        {
            h_xq2_piID_bg->Fill(*rec_xB, *rec_Q2);
            h_EminusPz_bg->Fill(*EminusPz);
        }
            
    }
    
    double gen_lumi = gen_lumi_pi[0]; // only one setting for pi bg

    h_xq2_all_bg->Scale(total_lumi/gen_lumi);
    h_xq2_acp_bg->Scale(total_lumi/gen_lumi);
    h_xq2_eff_bg->Scale(total_lumi/gen_lumi);
    h_xq2_pur_bg->Scale(total_lumi/gen_lumi);
    h_xq2_eID_bg->Scale(total_lumi/gen_lumi);
    h_xq2_piID_bg->Scale(total_lumi/gen_lumi);
    h_EminusPz_bg->Scale(total_lumi/gen_lumi);

    std::cout << "Entries in h_xq2_all_bg: " << h_xq2_all_bg->GetEntries() << std::endl;
    // std::cout << "Entries in h_EminusPz_bg: " << h_EminusPz_bg->GetEntries() << std::endl;
    // pi_bg_file->Close();

    std::cout << "Finished processing pi bg file." << std::endl;

    // o0o0o0o0o

    set_2d_scale(h_xq2_all);
    TCanvas* c_xq2_all = draw_2d_standard(h_xq2_all, "c_xq2_all", "all events", 700, 600, true, true);
    TCanvas* c_xq2_acp = draw_2d_standard(h_xq2_acp, "c_xq2_acp", "acp events", 700, 600, true, true);
    TCanvas* c_xq2_eID = draw_2d_standard(h_xq2_eID, "c_xq2_eID", "eID events", 700, 600, true, true);
    TCanvas* c_xq2_piID = draw_2d_standard(h_xq2_piID, "c_xq2_piID", "piID events", 700, 600, true, true);

    TCanvas* c_xq2_eff_bg = draw_2d_standard(h_xq2_eff_bg, "c_xq2_eff_bg", "PI bg events", 700, 600, true, true);
    TCanvas* c_xq2_piID_bg = draw_2d_standard(h_xq2_piID_bg, "c_xq2_piID_bg", "piID bg events", 700, 600, true, true);

    TH2F* h_xq2_acp_copy = (TH2F*)h_xq2_acp->Clone();
    process_eff_hist(h_xq2_acp_copy, h_xq2_all);
    TCanvas* c_xq2_acp_eff = draw_2d_efficiency(h_xq2_acp_copy, "c_xq2_acp_eff", "xq2 acp eff", 1400, 600, false, true);

    // TH2F* h_xq2_eff_copy = (TH2F*)h_xq2_eff->Clone();
    TH2F* h_xq2_eff_copy = (TH2F*)h_xq2_acp->Clone();
    h_xq2_eff_copy->Add(h_xq2_eff_bg);
    TCanvas* c_xq2_eff_copy = draw_2d_standard(h_xq2_eff_copy, "c_xq2_eff_copy", "eff copy events", 700, 600, true, true);
    TH2F* h_xq2_eff_bg_copy = (TH2F*)h_xq2_eff_bg->Clone();
    process_eff_hist(h_xq2_eff_bg_copy, h_xq2_eff_copy);
    TCanvas* c_xq2_contaimin = draw_2d_efficiency(h_xq2_eff_bg_copy, "c_xq2_contaimin", "xq2 contamination", 1400, 600, false, true);
    // draw_angle(c_xq2_eff_eff, h_xq2_eff_copy, 160.0, Ee*1.0, Eh*1.0);
    // draw_angle(c_xq2_eff_eff, h_xq2_eff_copy, 130.0, Ee*1.0, Eh*1.0);
    // draw_angle(c_xq2_eff_eff, h_xq2_eff_copy, 40.0, Ee*1.0, Eh*1.0);

    // TH2F* h_xq2_eff_copy2 = (TH2F*)h_xq2_eff->Clone();
    TH2F* h_xq2_eff_copy2 = (TH2F*)h_xq2_acp->Clone();
    h_xq2_eff_copy2->Add(h_xq2_eff_bg);
    TH2F* h_xq2_pur_copy = (TH2F*)h_xq2_pur->Clone();
    process_eff_hist(h_xq2_pur_copy, h_xq2_eff_copy2);
    TCanvas* c_xq2_pur_eff = draw_2d_efficiency(h_xq2_pur_copy, "c_xq2_pur_eff", "xq2 pur eff", 1400, 600, false, true);

    TH2F* h_xq2_eID_copy = (TH2F*)h_xq2_eID->Clone();
    process_eff_hist(h_xq2_eID_copy, h_xq2_eff_copy2);
    TCanvas* c_xq2_eID_eff = draw_2d_efficiency(h_xq2_eID_copy, "c_xq2_eID_eff", "xq2 eID eff", 1400, 600, false, true);

    TH2F* h_xq2_piID_copy = (TH2F*)h_xq2_piID->Clone();
    process_eff_hist(h_xq2_piID_copy, h_xq2_eff_copy2);
    TCanvas* c_xq2_piID_eff = draw_2d_efficiency(h_xq2_piID_copy, "c_xq2_piID_eff", "xq2 piID eff", 1400, 600, false, true);

    TCanvas* c_EminusPz = new TCanvas("c_EminusPz", "c_EminusPz", 1000, 600);
    h_EminusPz->SetLineColor(kBlue);
    h_EminusPz->SetFillColor(kBlue);
    h_EminusPz->SetFillStyle(3003);
    h_EminusPz->SetLineWidth(2);
    h_EminusPz->GetXaxis()->SetRangeUser(0, 50);
    double draw_max = h_EminusPz->GetMaximum() * 1.2;
    h_EminusPz->SetMaximum(draw_max);
    h_EminusPz->Draw("HIST");
    h_EminusPz_bg->SetLineColor(kRed);
    h_EminusPz_bg->SetFillColor(kRed);
    h_EminusPz_bg->SetFillStyle(3003);
    h_EminusPz_bg->SetLineWidth(2);
    h_EminusPz_bg->Draw("HIST SAME");

    TLegend* leg = new TLegend(0.8, 0.7, 0.95, 0.85);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(h_EminusPz, "e p", "l");
    leg->AddEntry(h_EminusPz_bg, "#gamma p", "l");
    leg->Draw("SAME");

    c_EminusPz->Modified();
    c_EminusPz->Update();

    TLine* line = new TLine(2*Ee, 0, 2*Ee, draw_max);
    line->SetLineColor(kBlack);
    line->SetLineStyle(7);
    line->Draw("SAME");

    draw_manager->LableAndCollect(c_xq2_piID_bg);
    c_xq2_piID_bg->SetFrameLineWidth(0);
    std::cout << "N events in c_xq2_piID_bg: " << h_xq2_piID_bg->Integral() << std::endl;
    // c_xq2_piID_bg->SaveAs(Form("../data/eID/%s_eID_xq2_piID_bg_wCut.png", setting.c_str()));
    // c_xq2_piID_bg->SaveAs(Form("../data/eID/%s_eID_xq2_piID_bg.png", setting.c_str()));

    draw_manager->LableAndCollect(c_EminusPz);
    // c_EminusPz->SaveAs(Form("../data/eID/%s_eID_EminusPz_wPi.png", setting.c_str()));

    draw_manager->LableAndCollect(c_xq2_contaimin);
    c_xq2_contaimin->SaveAs(Form("../data/eID/%s_eID_PiContam.png", setting.c_str()));

    return;
}