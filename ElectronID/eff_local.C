
#include "../GlobalUtil/ePICStyle.c"
#include "../GlobalUtil/Constants.hh"
#include "../GlobalUtil/drawHelper.h"

double total_lumi = 10; // fb^-1
double gen_lumi = 9.26e-4; // fb^-1

void draw_yline(TCanvas* &c, TH2F* h, double y_val)
{
    double x_min = h->GetXaxis()->GetXmin();
    double x_max = h->GetXaxis()->GetXmax();

    double s = 4*18*275; // GeV^2
    double q_min = (s - MASS_PROTON*MASS_PROTON)*x_min*y_val;
    double q_max = (s - MASS_PROTON*MASS_PROTON)*x_max*y_val;

    TGraph* g_theta = new TGraph();
    g_theta->AddPoint(x_min, q_min);
    g_theta->AddPoint(x_max, q_max);

    g_theta->SetLineColor(kRed);
    g_theta->SetLineWidth(2);
    // g_theta->SetLineStyle(7);

    c->cd();
    c->Update();
    g_theta->Draw("L SAME");

    TLatex *lt = new TLatex(1.05, q_max, Form("y = %.2f", y_val));
    // lt->SetNDC(kTRUE);
    lt->SetTextFont(42);
    lt->SetTextSize(0.03);
    lt->SetTextColor(kRed);
    lt->Draw();
}

void eff_local()
{
    SetePICStyle();

    // TH2F* h_xq2_all = BookTH2(Form("h_xq2_all"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -15, 0, n_q_bin,  -10, 0, kLightTemperature);
    TH2F* h_xq2_all = BookTH2(Form("h_xq2_all"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    TH2F* h_xq2_acp = BookTH2(Form("h_xq2_acp"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    
    TH2F* h_xq2_eff = BookTH2(Form("h_xq2_eff"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    TH2F* h_xq2_pur = BookTH2(Form("h_xq2_pur"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    
    TH2F* h_xq2_eID = BookTH2(Form("h_xq2_eID"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    TH2F* h_xq2_piID = BookTH2(Form("h_xq2_piID"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);

    // TFile* file = new TFile(Form("18x275_eIDreconnoBG.root"));
    TFile* file = new TFile(Form("../data/eID/18x275_eIDreconlowQ_BG.root"));

    TTreeReader reader("T_eID", file);

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

    Long64_t nEntries = reader.GetEntries();

    for( size_t ev = 0; ev < nEntries; ev++ ) 
    {
        reader.Next();
        // if(ev%100==0) 
        // cout << "Analysing file " << i << " event " << ev << "/" << nEntries << "\t\r" << std::flush;

        //  if ( y <= 0.01 || y >= 0.95 )
        //     continue;

        // if ( Q2 < 2 )
        //     continue;

        // if ( W2 < 4 )
        //     continue;

        if ( *status >= FOUND_MC )
        {
            h_xq2_all->Fill(*mc_xB, *mc_Q2);
            // std::cout << "Event " << ev << ": xB = " << *mc_xB << ", Q2 = " << *mc_Q2 << std::endl;
        }

        if ( *status >= FOUND_TRUTH )
            h_xq2_acp->Fill(*mc_xB, *mc_Q2);

        if ( *status >= FOUND_E )
            h_xq2_eff->Fill(*rec_xB, *rec_Q2);   

        if ( *status == FOUND_E )
            h_xq2_eID->Fill(*rec_xB, *rec_Q2);

        if ( *status >= FOUND_PI )
            h_xq2_piID->Fill(*rec_xB, *rec_Q2);
    }

        h_xq2_all->Scale(total_lumi/gen_lumi);
        h_xq2_acp->Scale(total_lumi/gen_lumi);
        h_xq2_eff->Scale(total_lumi/gen_lumi);
        // h_xq2_pur->Scale(total_lumi/gen_lumi);
        h_xq2_eID->Scale(total_lumi/gen_lumi);
        h_xq2_piID->Scale(total_lumi/gen_lumi);

    std::cout << total_lumi/gen_lumi << " " << h_xq2_eff->Integral() << std::endl;

    file->Close();

    // set_2d_scale(h_xq2_all);
    TCanvas* c_xq2_all = draw_2d_standard(h_xq2_all, "c_xq2_all", "all events", 700, 600, true, true);
    TCanvas* c_xq2_acp = draw_2d_standard(h_xq2_acp, "c_xq2_acp", "acp events", 700, 600, true, true);
    TCanvas* c_xq2_eff = draw_2d_standard(h_xq2_eff, "c_xq2_eff", "eff events", 700, 600, true, true);
    TCanvas* c_xq2_eID = draw_2d_standard(h_xq2_eID, "c_xq2_eID", "eID events", 700, 600, true, true);
    TCanvas* c_xq2_piID = draw_2d_standard(h_xq2_piID, "c_xq2_piID", "piID events", 700, 600, true, true);

    draw_yline(c_xq2_eff, h_xq2_eff, 0.95);

    // TH2F* h_xq2_acp_copy = (TH2F*)h_xq2_acp->Clone();
    // process_eff_hist(h_xq2_acp_copy, h_xq2_all);
    // TCanvas* c_xq2_acp_eff = draw_2d_efficiency(h_xq2_acp_copy, "c_xq2_acp_eff", "xq2 acp eff", 1400, 600, false, true);

    // TH2F* h_xq2_eID_copy = (TH2F*)h_xq2_eID->Clone();
    // process_eff_hist(h_xq2_eID_copy, h_xq2_acp);
    // TCanvas* c_xq2_eID_eff = draw_2d_efficiency(h_xq2_eID_copy, "c_xq2_eID_eff", "xq2 eID eff", 1400, 600, false, true);

    // TH2F* h_xq2_piID_copy = (TH2F*)h_xq2_piID->Clone();
    // process_eff_hist(h_xq2_piID_copy, h_xq2_acp);
    // TCanvas* c_xq2_piID_eff = draw_2d_efficiency(h_xq2_piID_copy, "c_xq2_piID_eff", "xq2 piID eff", 1400, 600, false, true);

    return;
}