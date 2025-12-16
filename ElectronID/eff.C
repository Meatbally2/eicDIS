
#include "../GlobalUtil/DrawManager.cc"
#include "../GlobalUtil/Constants.hh"
#include "../GlobalUtil/drawHelper.h"


const std::string group[4] = {"1", "10", "100", "1000"};
int e_set[4] = {1, 10, 100, 1000};
const double total_lumi = 10; // fb^-1
double lumi[4] = {6.73335E-03, 7.21215E-02, 1.48384E+00, 6.29541E+01}; // fb^-1
int n_group = 4;

// const std::string group[3] = {"1to10", "10to100", "100to1000"};
// const double total_lumi = 8.65*3; // fb^-1
// int e_set[3] = {1, 10, 100};
// int n_group = 3;

// double cs[3][2] = {{0.198440424611563, 0.205327493968226}, {4.04371412707044E-02, 4.41976212963417E-02}, {1.36416909784756E-03, 1.69583242740138E-03}};
// double ev[3][2] = {{333675, 666325}, {333694, 666306}, {333365, 666640}};

// smaller samples
// double cs[3][2] = {{0.19855280068406961, 0.20505655394590794}, {4.0447896030497726E-002, 4.4180993297364052E-002}, {1.3640235794259955E-003, 1.6954903002019237E-003}};
// double ev[3][2] = {{33202, 66798}, {66885, 133115}, {66556, 133444}};

void draw_angle(TCanvas* &c, TH2F* &h, double angle)
{
    double theta_deg = 180-angle; // formula from yellow report assume 0 deg is forward scattering for electron
    double theta = theta_deg * M_PI / 180.0;

    double Ee = 10.0; 
    double En = 166.0;

    // double Ee = 18.0; 
    // double En = 275.0;
    
    TAxis* xa = h->GetXaxis();
    int nx = xa->GetNbins();
    double last_Q2 = 0;

    TGraph* g_theta = new TGraph();
    for (int i = 1; i <= nx; ++i) {
        double low = xa->GetBinLowEdge(i);
        double high = xa->GetBinUpEdge(i);
        double center = xa->GetBinCenter(i);

        double xB[3] = {low, center, high};
        for ( int j = 0; j < 3; ++j) 
        {
            double Eprime = 2*Ee*En*xB[j] / (En*(1+cos(theta))*xB[j] + Ee*(1-cos(theta)));
            double Q2 = 2*Ee*Eprime*(1 - cos(theta));
            g_theta->AddPoint(xB[j], Q2);
            // std::cout << "xB: " << xB[j] << ", Q2: " << Q2 << std::endl;

            last_Q2 = Q2;
        }
    }

    g_theta->SetLineColor(kRed);
    g_theta->SetLineWidth(2);
    g_theta->SetLineStyle(7);

    c->cd();
    g_theta->Draw("L SAME");

    TLatex *lt = new TLatex(1.05, last_Q2, Form("#theta_{e} = %.0f#circ", angle));
    // lt->SetNDC(kTRUE);
    lt->SetTextFont(42);
    lt->SetTextSize(0.03);
    lt->SetTextColor(kRed);
    lt->Draw();


    return;
}

void eff()
{
    // std::string type_title = "e^{3}He";
    // std::string energy_title = Form("%dx%d GeV", 10, 166);
    std::string type_title = "ep";
    std::string energy_title = Form("%dx%d GeV", 18, 275);
    DrawManager* draw_manager = new DrawManager(type_title, energy_title, "25.10.2");
    draw_manager->SetEPIC();
    // draw_manager->SetDISrange(0.01, 0.95, 4, 2);
    draw_manager->SetLumi(8.65);

    TH2F* h_xq2_all = BookTH2(Form("h_xq2_all"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    TH2F* h_xq2_acp = BookTH2(Form("h_xq2_acp"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    
    TH2F* h_xq2_eff = BookTH2(Form("h_xq2_eff"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    TH2F* h_xq2_pur = BookTH2(Form("h_xq2_pur"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    
    TH2F* h_xq2_eID = BookTH2(Form("h_xq2_eID"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    TH2F* h_xq2_piID = BookTH2(Form("h_xq2_piID"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);

    for ( int i = 0; i < n_group; i ++ )
    {
        // TFile* file = new TFile(Form("../data/en_25_10_2/root_files/10x166_%s_eIDrecon_combined.root", group[i].c_str()));
        TFile* file = new TFile(Form("../data/ep_25_10_0/root_files/ep_18x275_minQ2=%i_eID_combined.root", e_set[i]));
        double gen_lumi = lumi[i];

        // data/ep_25_05_0/root_files/18x275_minQ2=1_eIDrecon_combined.root

        // double gen_lumi = ev[i][0]/(cs[i][0]*(1e-34/1e-43)) + ev[i][1]/(cs[i][1]*(1e-34/1e-43));
        // double n_lumi = ev[i][0]/(cs[i][0]*(1e-34/1e-43));
        // double p_lumi = ev[i][1]/(cs[i][1]*(1e-34/1e-43));

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


            if (*mc_y < 0.01 || *mc_y > 0.95) 
                continue;

            if (*mc_W2 < 4)
                continue;

            if (*mc_Q2 < 2)
                continue;

            // if ( e_set[i] == 1000 && * mc_Q2 < e_set[i] ) 
            //     continue;
            
            // if ( e_set[i] < 1000 &&  (*mc_Q2 > e_set[i]*10 || *mc_Q2 < e_set[i]) ) 
            //     continue;

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

        h_tmp_all->Scale(total_lumi/gen_lumi);
        h_tmp_acp->Scale(total_lumi/gen_lumi);
        h_tmp_eff->Scale(total_lumi/gen_lumi);
        h_tmp_pur->Scale(total_lumi/gen_lumi);
        h_tmp_eID->Scale(total_lumi/gen_lumi);
        h_tmp_piID->Scale(total_lumi/gen_lumi);

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
    draw_angle(c_xq2_eff_eff, h_xq2_eff_copy, 160.0);
    draw_angle(c_xq2_eff_eff, h_xq2_eff_copy, 130.0);
    draw_angle(c_xq2_eff_eff, h_xq2_eff_copy, 40.0);

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
    // c_xq2_eff_eff->SaveAs("../data/eID/10x166_en_eID_efficiency.png");
    // c_xq2_pur_eff->SaveAs("../data/eID/10x166_en_eID_purity.png");

    // c_xq2_all->SaveAs("../data/eID/10x166_en_raw.png");

    return;
}