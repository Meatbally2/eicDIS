#include "../GlobalUtil/luminosityTable.h"
#include "../GlobalUtil/DrawManager.cc"
#include "../GlobalUtil/Constants.hh"
#include "../GlobalUtil/drawHelper.h"

#include <Math/LorentzVector.h>
using ROOT::Math::PxPyPzEVector;

#include "draw_angle.C"

const std::string group[4] = {"minQ2=1", "minQ2=10", "minQ2=100", "minQ2=1000"};

void ReverseXAxis(TH1* h);
void draw_eta_axis(TH1* hist);
void draw_leg_pad(const std::string& type, const std::string& energy, const std::string& campaign);

void eff(int beam_type, int Ee, int Eh)
{
    // ePIC plotting style setup
    std::string type_title[5] = {"e^{3}He", "ep", "#gammap", "ep w. BeamBG", "ep"};
    std::string energy_title = beam_type ? Form("%dx%d GeV", Ee, Eh) : Form("%dx%d GeV/A", Ee, Eh);
    std::string campaign = beam_type ? "26.02.0" : "26.02.0";
    if ( beam_type == 4 )
        campaign = "25.10.4";
    if ( beam_type == 3 )
        campaign = "26.02.0";
    DrawManager* draw_manager = new DrawManager(type_title[beam_type], energy_title, campaign);
    draw_manager->SetEPIC("Work in Progress");

    draw_manager->SetDISrange(0.01, 0.95, 4, 2);

    int n_group = beam_type == 0 ? 3 : 4;
    if ( beam_type == 3 )
        n_group = 1;

    double text_lumi = 10; // fb^-1
    if ( beam_type == 0 && Ee == 10 && Eh == 166 )
        text_lumi = 8.65; // fb^-1
    if ( beam_type == 4 && Ee == 10 && Eh == 130 )
        text_lumi = 5.33; // fb^-1
    if ( beam_type == 4 && Ee == 10 && Eh == 250 )
        text_lumi = 9.18; // fb^-1
    draw_manager->SetLumi(text_lumi);
    
    std::string setting = beam_type == 0 ? Form("eHe3_%dx%d", (int)Ee, (int)Eh) : Form("ep_%dx%d", (int)Ee, (int)Eh);
    if ( beam_type == 3 )
        setting = Form("beamBG_%dx%d", (int)Ee, (int)Eh);

    TH2F* h_xq2_all = BookTH2(Form("h_xq2_all"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    TH2F* h_xq2_acp = BookTH2(Form("h_xq2_acp"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    TH2F* h_xq2_trk = BookTH2(Form("h_xq2_trk"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    TH2F* h_xq2_cal = BookTH2(Form("h_xq2_cal"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    TH2F* h_xq2_eff = BookTH2(Form("h_xq2_eff"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    TH2F* h_xq2_pur = BookTH2(Form("h_xq2_pur"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    TH2F* h_xq2_eID = BookTH2(Form("h_xq2_eID"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    TH2F* h_xq2_piID = BookTH2(Form("h_xq2_piID"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);

    TH1F* h_angular_all = new TH1F("h_angular_all", ";#theta (deg);%", 180, 0, 180);
    TH1F* h_angular_acp = new TH1F("h_angular_acp", ";#theta (deg);%", 180, 0, 180);
    TH1F* h_angular_trk = new TH1F("h_angular_trk", ";#theta (deg);%", 180, 0, 180);
    TH1F* h_angular_cal = new TH1F("h_angular_cal", ";#theta (deg);%", 180, 0, 180);
    TH1F* h_angular_eff = new TH1F("h_angular_eff", ";#theta (deg);%", 180, 0, 180);
    TH1F* h_angular_pur = new TH1F("h_angular_pur", ";#theta (deg);%", 180, 0, 180);
    TH1F* h_angular_eID = new TH1F("h_angular_eID", ";#theta (deg);%", 180, 0, 180);
    TH1F* h_angular_piID = new TH1F("h_angular_piID", ";#theta (deg);%", 180, 0, 180);

    for ( int i = 0; i < n_group; i ++ )
    {
        std::string date = beam_type == 0 ? "en_26_02_2" : "ep_26_02_0";
         if ( beam_type == 3 )
            date = "ep_26_02_0";
        if ( beam_type == 4 )
            date = "ep_25_10_4";
        // TFile* file = beam_type == 2 ? new TFile("/Users/winl/eic/eicDIS/data/beamBG_10x100_eid_combined.root") : new TFile(Form("../data/%s/root_files/%s_%s_eid_combined.root", date.c_str(), setting.c_str(), group[i].c_str()));
        TFile* file  = new TFile(Form("../data/%s/root_files/%s_%s_eid_combined.root", date.c_str(), setting.c_str(), group[i].c_str()));

        TH2F* h_tmp_all = BookTH2(Form("h_tmp_all_%d", i), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
        TH2F* h_tmp_acp = BookTH2(Form("h_tmp_acp_%d", i), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
        TH2F* h_tmp_trk = BookTH2(Form("h_tmp_trk_%d", i), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
        TH2F* h_tmp_cal = BookTH2(Form("h_tmp_cal_%d", i), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
        TH2F* h_tmp_eff = BookTH2(Form("h_tmp_eff_%d", i), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
        TH2F* h_tmp_pur = BookTH2(Form("h_tmp_pur_%d", i), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
        TH2F* h_tmp_eID = BookTH2(Form("h_tmp_eID_%d", i), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
        TH2F* h_tmp_piID = BookTH2(Form("h_tmp_piID_%d", i), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);

        TH1F* h_tmp_angular_all = new TH1F(Form("h_tmp_angular_all_%d", i), ";#theta (deg);%", 180, 0, 180);
        TH1F* h_tmp_angular_acp = new TH1F(Form("h_tmp_angular_acp_%d", i), ";#theta (deg);%", 180, 0, 180);
        TH1F* h_tmp_angular_trk = new TH1F(Form("h_tmp_angular_trk_%d", i), ";#theta (deg);%", 180, 0, 180);
        TH1F* h_tmp_angular_cal = new TH1F(Form("h_tmp_angular_cal_%d", i), ";#theta (deg);%", 180, 0, 180);
        TH1F* h_tmp_angular_eff = new TH1F(Form("h_tmp_angular_eff_%d", i), ";#theta (deg);%", 180, 0, 180);
        TH1F* h_tmp_angular_pur = new TH1F(Form("h_tmp_angular_pur_%d", i), ";#theta (deg);%", 180, 0, 180);
        TH1F* h_tmp_angular_eID = new TH1F(Form("h_tmp_angular_eID_%d", i), ";#theta (deg);%", 180, 0, 180);
        TH1F* h_tmp_angular_piID = new TH1F(Form("h_tmp_angular_piID_%d", i), ";#theta (deg);%", 180, 0, 180);

        TTreeReader reader("T_eID", file);

        TTreeReaderValue<int>    status(reader, "eID_status");
        TTreeReaderValue<int>    rec_status(reader, "eRecon_status");

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

        TTreeReaderValue<PxPyPzEVector> vMCe(reader, "vMC_e");
        TTreeReaderValue<PxPyPzEVector> vTRe(reader, "vTRACK_e");
        TTreeReaderValue<PxPyPzEVector> vCLe(reader, "vCLUSTER_e");

         TTreeReaderValue<double> EoP(reader, "EoP");
        TTreeReaderValue<double> EminusPz(reader, "EminusPz");

        Long64_t nEntries = reader.GetEntries();

        // get generated lumi
        double total_lumi = beam_type == 0 ? text_lumi*3 : text_lumi; // fb^-1
        double gen_lumi = 0;
        gen_lumi = get_lumi(beam_type, Ee, Eh, i, nEntries, 0);
        if ( beam_type == 2 )
            gen_lumi = get_lumi(1, Ee, Eh, i, nEntries, 0);
        if ( gen_lumi == 0 )
        {
            std::cout << "** Lumi table not found! " << std::endl;
            return;
        }

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

            double theta = 180-vMCe->Theta()*(180./M_PI); // flip to make the plot go from 180 to 0 deg
            // double theta = vMCe->Theta()*(180./M_PI); 

            if ( *status >= FOUND_MC )
            {
                h_tmp_all->Fill(*mc_xB, *mc_Q2);
                h_tmp_angular_all->Fill(theta);
            }
                
            if ( *status >= FOUND_TRUTH )
            {
                h_tmp_acp->Fill(*mc_xB, *mc_Q2);
                h_tmp_angular_acp->Fill(theta);
            }

            if ( *rec_status == FOUND_TRACK_ONLY || *rec_status == FOUND_BOTH )
            {
                h_tmp_trk->Fill(*mc_xB, *mc_Q2);
                h_tmp_angular_trk->Fill(theta);
            }

            if ( *rec_status == FOUND_CLUSTER_ONLY || *rec_status == FOUND_BOTH )
            {
                h_tmp_cal->Fill(*mc_xB, *mc_Q2);
                h_tmp_angular_cal->Fill(theta);
            }

            // if ( ( (180-theta > 120 && 180-theta < 150) || (180-theta > 35 && 180-theta < 60) ) && *EminusPz < 5 )
            //         continue;

            // if ( (180-theta > 125 && 180-theta < 155) && *EminusPz < 5 )
            //         continue;

            // if ( *mc_Q2 < 10 && (*EoP < 0.8 || *EoP > 1.2) )
            //     continue;
            
            if ( *status >= FOUND_E )
            {
                h_tmp_eff->Fill(*mc_xB, *mc_Q2);   
                h_tmp_angular_eff->Fill(theta);
            }

            if ( *status == FOUND_E )
            {
                h_tmp_pur->Fill(*mc_xB, *mc_Q2);
                h_tmp_eID->Fill(*mc_xB, *mc_Q2);
                h_tmp_angular_pur->Fill(theta);
                h_tmp_angular_eID->Fill(theta);
            }

            if ( *status == FOUND_PI )
            {
                h_tmp_piID->Fill(*mc_xB, *mc_Q2);
                h_tmp_angular_piID->Fill(theta);
            }
        }

        h_tmp_all->Scale(total_lumi/gen_lumi);
        h_tmp_acp->Scale(total_lumi/gen_lumi);
        h_tmp_eff->Scale(total_lumi/gen_lumi);
        h_tmp_pur->Scale(total_lumi/gen_lumi);
        h_tmp_eID->Scale(total_lumi/gen_lumi);
        h_tmp_piID->Scale(total_lumi/gen_lumi);
        h_tmp_trk->Scale(total_lumi/gen_lumi);
        h_tmp_cal->Scale(total_lumi/gen_lumi);
        h_tmp_angular_trk->Scale(total_lumi/gen_lumi);
        h_tmp_angular_cal->Scale(total_lumi/gen_lumi);
        h_tmp_angular_all->Scale(total_lumi/gen_lumi);
        h_tmp_angular_acp->Scale(total_lumi/gen_lumi);
        h_tmp_angular_eff->Scale(total_lumi/gen_lumi);
        h_tmp_angular_pur->Scale(total_lumi/gen_lumi);
        h_tmp_angular_eID->Scale(total_lumi/gen_lumi);
        h_tmp_angular_piID->Scale(total_lumi/gen_lumi);

        h_xq2_all->Add(h_tmp_all);
        h_xq2_acp->Add(h_tmp_acp);
        h_xq2_eff->Add(h_tmp_eff);
        h_xq2_trk->Add(h_tmp_trk);
        h_xq2_cal->Add(h_tmp_cal);
        h_xq2_pur->Add(h_tmp_pur);
        h_xq2_eID->Add(h_tmp_eID);
        h_xq2_piID->Add(h_tmp_piID);
        h_angular_all->Add(h_tmp_angular_all);
        h_angular_acp->Add(h_tmp_angular_acp);
        h_angular_eff->Add(h_tmp_angular_eff);
        h_angular_pur->Add(h_tmp_angular_pur);
        h_angular_eID->Add(h_tmp_angular_eID);
        h_angular_piID->Add(h_tmp_angular_piID);
        h_angular_trk->Add(h_tmp_angular_trk);
        h_angular_cal->Add(h_tmp_angular_cal);

        file->Close();
    }

    set_2d_scale(h_xq2_all);
    TCanvas* c_xq2_all = draw_2d_standard(h_xq2_all, "c_xq2_all", "all events", 700, 600, true, true);
    // TCanvas* c_xq2_trk = draw_2d_standard(h_xq2_trk, "c_xq2_trk", "track events", 700, 600, true, true);
    // TCanvas* c_xq2_cal = draw_2d_standard(h_xq2_cal, "c_xq2_cal", "cal events", 700, 600, true, true);
    TCanvas* c_xq2_acp = draw_2d_standard(h_xq2_acp, "c_xq2_acp", "acp events", 700, 600, true, true);
    TCanvas* c_xq2_eID = draw_2d_standard(h_xq2_eID, "c_xq2_eID", "eID events", 700, 600, true, true);
    TCanvas* c_xq2_piID = draw_2d_standard(h_xq2_piID, "c_xq2_piID", "piID events", 700, 600, true, true);

    TH2F* h_xq2_acp_copy = (TH2F*)h_xq2_acp->Clone();
    process_eff_hist(h_xq2_acp_copy, h_xq2_all);
    TCanvas* c_xq2_acp_eff = draw_2d_efficiency(h_xq2_acp_copy, "c_xq2_acp_eff", "xq2 acp eff", 1400, 600, false, true);

    process_eff_hist(h_xq2_trk, h_xq2_all);
    TCanvas* c_xq2_trk_eff = draw_2d_efficiency(h_xq2_trk, "c_xq2_trk", "xq2 trk eff", 1400, 600, false, true);

    process_eff_hist(h_xq2_cal, h_xq2_all);
    TCanvas* c_xq2_cal_eff = draw_2d_efficiency(h_xq2_cal, "c_xq2_cal", "xq2 cal eff", 1400, 600, false, true);

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
    draw_manager->LableAndCollect(c_xq2_trk_eff);
    draw_manager->LableAndCollect(c_xq2_cal_eff);
    draw_manager->LableAndCollect(c_xq2_pur_eff);
    draw_manager->LableAndCollect(c_xq2_eID_eff);
    draw_manager->LableAndCollect(c_xq2_piID_eff);

    // gStyle->SetImageScaling(3.);
    c_xq2_acp_eff->SaveAs(Form("../data/eID/%s_eID_acceptance.png", setting.c_str()));
    c_xq2_eff_eff->SaveAs(Form("../data/eID/%s_eID_efficiency.png", setting.c_str()));
    c_xq2_pur_eff->SaveAs(Form("../data/eID/%s_eID_purity.png", setting.c_str()));
    c_xq2_eID_eff->SaveAs(Form("../data/eID/%s_eID_overall.png", setting.c_str()));
    c_xq2_trk_eff->SaveAs(Form("../data/eID/%s_eID_trk.png", setting.c_str()));
    c_xq2_cal_eff->SaveAs(Form("../data/eID/%s_eID_cal.png", setting.c_str()));

    // c_xq2_all->SaveAs("../data/eID/10x166_en_raw.png");

    TCanvas* c_angular = new TCanvas("c_angular", "c_angular", 1400, 600);

    TPad *draw_pad = new TPad("draw_pad", "draw_pad", 0, 0, 0.8, 1);
    draw_pad->SetFillColor(0);
    draw_pad->SetFillStyle(0);
    draw_pad->Draw();
    draw_pad->SetLeftMargin(1);
    draw_pad->SetBottomMargin(0.25);
    draw_pad->cd();

    TH1F* h_angular_formater = (TH1F*)h_angular_acp->Clone("h_angular_formater");
    h_angular_formater->Reset();
    h_angular_formater->Draw();
    h_angular_formater->SetMaximum(1.05);
    h_angular_formater->SetMinimum(0);

    h_angular_formater->GetYaxis()->CenterTitle();
    h_angular_formater->GetYaxis()->SetTitleOffset(1.);

    TEfficiency*  h_acceptance_angular = new TEfficiency(*h_angular_acp,*h_angular_all);
    TEfficiency*  h_trk_angular = new TEfficiency(*h_angular_trk,*h_angular_all);
    TEfficiency*  h_cal_angular = new TEfficiency(*h_angular_cal,*h_angular_all);
    TEfficiency*  h_efficiency_angular = new TEfficiency(*h_angular_eff,*h_angular_all);
    TEfficiency*  h_purity_angular = new TEfficiency(*h_angular_pur,*h_angular_eff);
    TEfficiency*  h_eID_angular = new TEfficiency(*h_angular_eID,*h_angular_all);
    TEfficiency*  h_piID_angular = new TEfficiency(*h_angular_piID,*h_angular_acp);

    // h_acceptance_angular->SetLineColor(kBlack);
    // h_acceptance_angular->SetMarkerColor(kBlack);

    h_trk_angular->SetLineColor(kOrange+1);
    h_trk_angular->SetMarkerColor(kOrange+1);

    h_cal_angular->SetLineColor(kBlack);
    h_cal_angular->SetMarkerColor(kBlack);

    h_efficiency_angular->SetLineColor(kBlue+1);
    h_efficiency_angular->SetMarkerColor(kBlue+1);

    h_purity_angular->SetLineColor(kGreen+2);
    h_purity_angular->SetMarkerColor(kGreen+2);

    h_eID_angular->SetLineColor(kRed+1);
    h_eID_angular->SetMarkerColor(kRed+1);

    h_piID_angular->SetLineColor(kGray+1);
    h_piID_angular->SetMarkerColor(kGray+1);

    h_acceptance_angular->SetStatisticOption(TEfficiency::kFNormal);
    h_efficiency_angular->SetStatisticOption(TEfficiency::kFNormal);
    h_trk_angular->SetStatisticOption(TEfficiency::kFNormal);
    h_cal_angular->SetStatisticOption(TEfficiency::kFNormal);
    h_purity_angular->SetStatisticOption(TEfficiency::kFNormal);
    h_eID_angular->SetStatisticOption(TEfficiency::kFNormal);
    h_piID_angular->SetStatisticOption(TEfficiency::kFNormal);

    // h_acceptance_angular->Draw("SAME");
    h_efficiency_angular->Draw("SAME");
    h_trk_angular->Draw("SAME");
    h_cal_angular->Draw("SAME");
    h_purity_angular->Draw("SAME ");
    h_eID_angular->Draw("SAME");
    h_piID_angular->Draw("SAME");

    ReverseXAxis(h_angular_formater);
    draw_eta_axis(h_angular_formater);

    c_angular->cd();
    draw_leg_pad(type_title[beam_type], energy_title, campaign);

    TLegend* leg = new TLegend(0.1, 0.3, 0.9, 0.65);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.08);
    leg->AddEntry(h_trk_angular, "e track eff.", "LPE");
    leg->AddEntry(h_cal_angular, "e cluster eff.", "LPE");
    // leg->AddEntry(h_acceptance_angular, "Track x Cal. eff", "LPE");
    leg->AddEntry(h_efficiency_angular, "eID efficiency", "LPE");
    leg->AddEntry(h_purity_angular, "eID purity", "LPE");
    leg->AddEntry(h_eID_angular, "Overall eID rate", "LPE");
    leg->AddEntry(h_piID_angular, "#pi contamination", "LPE");
    leg->Draw();
    
    c_angular->SaveAs(Form("../data/eID/%s_eID_angular.png", setting.c_str()));

    return;
}

void ReverseXAxis(TH1* h)
{
   // Remove the current axis
   h->GetXaxis()->SetTitleOffset(999);
   h->GetXaxis()->SetLabelOffset(999);
   h->GetXaxis()->SetTickLength(0);
 
   // Redraw the new axis
   gPad->Update();

   TGaxis *newaxis = new TGaxis(gPad->GetUxmax(), gPad->GetUymin(), gPad->GetUxmin(), gPad->GetUymin(), h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax(), 510,"-");
   newaxis->SetLabelOffset(-0.03);
   newaxis->SetTitleFont(42);
   newaxis->SetLabelFont(42);
   newaxis->SetTitleOffset(1.2);
   newaxis->CenterTitle();
   newaxis->SetTitle("#theta_{mc} (deg)");
   newaxis->Draw();
}

void draw_eta_axis(TH1* hist)
{
    if (!hist || !gPad) return;

    gPad->Update();

    const double xmin = hist->GetXaxis()->GetXmin(); // theta min (deg)
    const double xmax = hist->GetXaxis()->GetXmax(); // theta max (deg)
    const double ymin = gPad->GetUymin();
    const double ymax = gPad->GetUymax();
    const double yr   = ymax - ymin;

    // Axis baseline below the original x-axis
    const double y_axis = ymin - 0.20 * yr;
    const double tick_h = 0.02 * yr;

    // Draw custom eta axis line
    TLine* axis_line = new TLine(xmin, y_axis, xmax, y_axis);
    axis_line->SetLineWidth(2);
    axis_line->Draw();

    // Fixed eta ticks
    const int nEta = 9;
    const double etaVals[nEta] = {-4, -2, -1, -0.5, 0, 0.5, 1, 2, 4}; 

    TLatex latex;
    latex.SetTextAlign(23);     // centered x, top y
    latex.SetTextSize(0.042);

    for (int i = 0; i < nEta; ++i)
    {
        const double eta = etaVals[i];

        // theta(eta) = 2 * atan(exp(-eta)) in radians -> degrees
        const double theta_deg = 2.0 * TMath::ATan(TMath::Exp(-eta)) * TMath::RadToDeg();

        if (theta_deg < xmin || theta_deg > xmax) continue;

        TLine* tick = new TLine(theta_deg, y_axis, theta_deg, y_axis + tick_h);
        tick->Draw();

        latex.DrawLatex(180-theta_deg, y_axis - 0.013 * yr, Form("%.1f", eta));
    }

    // Axis title
    TLatex title;
    title.SetTextAlign(23);
    title.SetTextSize(0.042);
    title.DrawLatex((xmin + xmax) * 0.5, y_axis - 0.07 * yr, "#eta");
}

void draw_leg_pad(const std::string& type, const std::string& energy, const std::string& campaign)
{
    TPad *leg_pad = new TPad("leg_pad", "leg_pad", 0.79, 0, 1, 1);
    leg_pad->SetFillColor(0);
    leg_pad->SetFillStyle(0);
    leg_pad->SetLeftMargin(0);
    leg_pad->Draw();
    leg_pad->cd();

    // ===== Add ePIC logo to the figure ======
    TImage *logo = TImage::Open("../GlobalUtil/EPIC-logo_black_transparent.png");
    logo->SetConstRatio(kTRUE);

    // TPad *logo_pad = new TPad("logo_pad", "logo_pad", 0.0, 0.75, 0.4, 0.9);
    TPad *logo_pad = new TPad("logo_pad", "logo_pad", 0.0, 0.75, 0.35, 0.9);
    logo_pad->SetFillStyle(0);
    logo_pad->SetFillColor(0);
    logo_pad->SetFrameFillStyle(0);
    logo_pad->SetBorderMode(0);
    logo_pad->SetBorderSize(0);
    logo_pad->SetFrameBorderMode(0);
    logo_pad->SetLineWidth(0);
    logo_pad->SetLineColor(0);
    logo_pad->SetFrameLineWidth(0);
    logo_pad->SetFrameLineColor(0);
    logo_pad->SetMargin(0, 0, 0, 0);

    logo_pad->Draw();
    logo_pad->cd();
    logo->Draw();

    // ===== Add type ======
    leg_pad->cd();
    // TLatex *Text_plotType = new TLatex(0.45, 0.815, "Work in Progress");
    TLatex *Text_plotType = new TLatex(0.37, 0.825, "Work in Progress");
    Text_plotType->SetNDC();
    Text_plotType->SetTextSize(0.13);
    Text_plotType->SetTextSize(0.08);
    Text_plotType->SetTextFont(62);
    Text_plotType->SetTextAlign(12);
    Text_plotType->Draw();

    // ===== Add energy ======
    TLatex *Text_energy = new TLatex(0.1, 0.7, Form("%s, %s", type.c_str(), energy.c_str()));
    Text_energy->SetNDC();
    Text_energy->SetTextSize(0.1);
    Text_energy->SetTextFont(52);
    Text_energy->SetTextAlign(12);
    Text_energy->Draw();

    // ===== Add campaign ======
    TLatex *Text_campaign = new TLatex(0.45, 0.1, Form("Simulation campaign: %s", campaign.c_str()));
    Text_campaign->SetNDC();
    Text_campaign->SetTextSize(0.055);
    Text_campaign->SetTextFont(52);
    Text_campaign->SetTextAlign(22);
    Text_campaign->Draw();
}