#include "../GlobalUtil/DrawManager.cc"
#include "../GlobalUtil/HistManager.cc"
#include "../GlobalUtil/Constants.hh"
#include "../GlobalUtil/luminosityTable.h"

#include <Math/LorentzVector.h>
using ROOT::Math::PxPyPzEVector;

#include "draw_angle.C"

using HistGroup1D = HistManager::HistGroup1D;
using HistGroup2D = HistManager::HistGroup2D;

const std::string group[4] = {"minQ2=1", "minQ2=10", "minQ2=100", "minQ2=1000"};

void ReverseXAxis(TH1* h);
void draw_eta_axis(TH1* hist);
void draw_leg_pad(const std::string& type, const std::string& energy, const std::string& campaign);

void eff(int beam_type, int Ee, int Eh)
{
    // ePIC plotting style setup
    std::string type_title[5] = {"e^{3}He", "ep", "#gammap", "ep + 2#mus beamBG", "ep"};
    std::string energy_title = beam_type ? Form("%dx%d GeV", Ee, Eh) : Form("%dx%d GeV/A", Ee, Eh);
    std::string campaign = "26.03.1";
    if ( beam_type == 1 || beam_type == 3 )
        campaign = "26.04.1";
    DrawManager* draw_manager = new DrawManager(type_title[beam_type], energy_title, campaign);
    draw_manager->SetEPIC("Work in Progress");

    draw_manager->SetDISrange(0.01, 0.95, 4, 2);

    int n_group = beam_type == 0 ? 3 : 4;
    if ( beam_type == 2 )
        n_group = 1;

    double text_lumi = 1; // fb^-1
    if ( beam_type == 0 && Ee == 10 && Eh == 166 )
        text_lumi = 1.5; // fb^-1
    if ( beam_type == 4 && Ee == 10 && Eh == 130 )
        text_lumi = 1.0; // fb^-1
    if ( beam_type == 4 && Ee == 10 && Eh == 250 )
        text_lumi = 2.5; // fb^-1
    draw_manager->SetLumi(text_lumi);

    // histograms
    HistManager* hist_manager = new HistManager();
    
    std::string setting = beam_type == 0 ? Form("eHe3_%dx%d", (int)Ee, (int)Eh) : Form("ep_%dx%d", (int)Ee, (int)Eh);
    if ( beam_type == 2 )
        setting = Form("piBG_%dx%d", (int)Ee, (int)Eh);

    auto h_xq2_all = hist_manager->BookHistograms(Form("h_xq2_all"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5);
    auto h_xq2_acp = hist_manager->BookHistograms(Form("h_xq2_acp"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5);
    auto h_xq2_trk = hist_manager->BookHistograms(Form("h_xq2_trk"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5);
    auto h_xq2_cal = hist_manager->BookHistograms(Form("h_xq2_cal"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5);
    auto h_xq2_eff = hist_manager->BookHistograms(Form("h_xq2_eff"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5);
    auto h_xq2_pur = hist_manager->BookHistograms(Form("h_xq2_pur"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5);
    auto h_xq2_eID = hist_manager->BookHistograms(Form("h_xq2_eID"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5);
    auto h_xq2_piID = hist_manager->BookHistograms(Form("h_xq2_piID"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5);

    auto h_angular_all = hist_manager->MakeHistograms(Form("h_angular_all"), ";#theta (deg);%", 180, 0, 180);
    auto h_angular_acp = hist_manager->MakeHistograms(Form("h_angular_acp"), ";#theta (deg);%", 180, 0, 180);
    auto h_angular_trk = hist_manager->MakeHistograms(Form("h_angular_trk"), ";#theta (deg);%", 180, 0, 180);
    auto h_angular_cal = hist_manager->MakeHistograms(Form("h_angular_cal"), ";#theta (deg);%", 180, 0, 180);
    auto h_angular_eff = hist_manager->MakeHistograms(Form("h_angular_eff"), ";#theta (deg);%", 180, 0, 180);
    auto h_angular_pur = hist_manager->MakeHistograms(Form("h_angular_pur"), ";#theta (deg);%", 180, 0, 180);
    auto h_angular_eID = hist_manager->MakeHistograms(Form("h_angular_eID"), ";#theta (deg);%", 180, 0, 180);
    auto h_angular_piID = hist_manager->MakeHistograms(Form("h_angular_piID"), ";#theta (deg);%", 180, 0, 180);

    for ( int i = 0; i < n_group; i ++ )
    {
        std::string date = beam_type == 0 ? "en_26_03_1" : "ep_26_04_1";
        // TFile* file = beam_type == 2 ? new TFile("/Users/winl/eic/eicDIS/data/beamBG_10x100_eid_combined.root") : new TFile(Form("../data/%s/root_files/%s_%s_eid_combined.root", date.c_str(), setting.c_str(), group[i].c_str()));
        if ( beam_type == 3 )
            date = "ep_26_03_1";
        std::string file_path = beam_type == 3 ? Form("../data/%s/root_files/%s_%s_eid_pythia6_combined.root", date.c_str(), setting.c_str(), group[i].c_str()) : Form("../data/%s/root_files/%s_%s_eid_combined.root", date.c_str(), setting.c_str(), group[i].c_str());
        if ( beam_type == 2 )
            file_path = Form("../data/%s/root_files/%s_eid_combine.root", date.c_str(), setting.c_str());
        TFile* file  = new TFile(file_path.c_str());

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

        // TFile* beamFile = new TFile(Form("../data/%s/root_files/%s_%s_beam_combined.root", date.c_str(), setting.c_str(), group[i].c_str()));
        // TTreeReader beam_reader("T_Beam", beamFile);
        // TTreeReaderValue<int> N_PDG(beam_reader, "N_PDG");

        Long64_t nEntries = reader.GetEntries();

        int en_num = 0;
        int ep_num = 0;

        for( size_t ev = 0; ev < nEntries; ev++ ) 
        {
            reader.Next();
            // beam_reader.Next();

            int fill_index = -1;
            // if ( *N_PDG == ID_PROTON )
            // {
                ep_num ++;
                fill_index = 1;
            // }
            // else if ( *N_PDG == ID_NEUTRON )
            // {
            //     en_num ++;
            //     fill_index = 0;
            // }

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
                h_xq2_all.ind[i][fill_index]->Fill(*mc_xB, *mc_Q2);
                h_angular_all.ind[i][fill_index]->Fill(theta);
            }
                
            if ( *status >= FOUND_TRUTH )
            {
                h_xq2_acp.ind[i][fill_index]->Fill(*mc_xB, *mc_Q2);
                h_angular_acp.ind[i][fill_index]->Fill(theta);
            }

            if ( *rec_status == FOUND_TRACK_ONLY || *rec_status == FOUND_BOTH )
            {
                h_xq2_trk.ind[i][fill_index]->Fill(*mc_xB, *mc_Q2);
                h_angular_trk.ind[i][fill_index]->Fill(theta);
            }

            if ( *rec_status == FOUND_CLUSTER_ONLY || *rec_status == FOUND_BOTH )
            {
                h_xq2_cal.ind[i][fill_index]->Fill(*mc_xB, *mc_Q2);
                h_angular_cal.ind[i][fill_index]->Fill(theta);
            }

            // if ( ( (180-theta > 120 && 180-theta < 150) || (180-theta > 35 && 180-theta < 60) ) && *EminusPz < 5 )
            //         continue;

            // if ( (180-theta > 125 && 180-theta < 155) && *EminusPz < 5 )
            //         continue;

            // if ( *mc_Q2 < 10 && (*EoP < 0.8 || *EoP > 1.2) )
            //     continue;
            
            if ( *status >= FOUND_E )
            {
                h_xq2_eff.ind[i][fill_index]->Fill(*mc_xB, *mc_Q2);   
                h_angular_eff.ind[i][fill_index]->Fill(theta);
            }

            if ( *status == FOUND_E )
            {
                h_xq2_pur.ind[i][fill_index]->Fill(*mc_xB, *mc_Q2);
                h_xq2_eID.ind[i][fill_index]->Fill(*mc_xB, *mc_Q2);
                h_angular_pur.ind[i][fill_index]->Fill(theta);
                h_angular_eID.ind[i][fill_index]->Fill(theta);
            }

            if ( *status == FOUND_PI || *status == FOUND_OTHERS )
            {
                h_xq2_piID.ind[i][fill_index]->Fill(*mc_xB, *mc_Q2);
                h_angular_piID.ind[i][fill_index]->Fill(theta);
            }
        }

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

        // if ( beam_type == 0)
        // {
        //     double n_lumi = get_lumi(beam_type, Ee, Eh, i, en_num, 0);
        //     double p_lumi = get_lumi(beam_type, Ee, Eh, i, 0, ep_num);
        //     hist_manager->SumHistograms(i, total_lumi / n_lumi, total_lumi / p_lumi);
        // }
        // else
        // {
            double p_lumi = get_lumi(beam_type, Ee, Eh, i, ep_num, 0);
            hist_manager->SumHistograms(i, 1, total_lumi / p_lumi);
        // }
        
        file->Close();
    }

    set_2d_scale(h_xq2_all.sum);
    TCanvas* c_xq2_all = draw_2d_standard(h_xq2_all.sum, "c_xq2_all", "all events", 700, 600, true, true);
    // TCanvas* c_xq2_trk = draw_2d_standard(h_xq2_trk, "c_xq2_trk", "track events", 700, 600, true, true);
    // TCanvas* c_xq2_cal = draw_2d_standard(h_xq2_cal, "c_xq2_cal", "cal events", 700, 600, true, true);
    TCanvas* c_xq2_acp = draw_2d_standard(h_xq2_acp.sum, "c_xq2_acp", "acp events", 700, 600, true, true);
    TCanvas* c_xq2_eID = draw_2d_standard(h_xq2_eID.sum, "c_xq2_eID", "eID events", 700, 600, true, true);
    TCanvas* c_xq2_piID = draw_2d_standard(h_xq2_piID.sum, "c_xq2_piID", "piID events", 700, 600, true, true);

    TH2F* h_xq2_acp_copy = (TH2F*)h_xq2_acp.sum->Clone();
    process_eff_hist(h_xq2_acp_copy, h_xq2_all.sum);
    TCanvas* c_xq2_acp_eff = draw_2d_efficiency(h_xq2_acp_copy, "c_xq2_acp_eff", "xq2 acp eff", 1400, 600, false, true);

    process_eff_hist(h_xq2_trk.sum, h_xq2_all.sum);
    TCanvas* c_xq2_trk_eff = draw_2d_efficiency(h_xq2_trk.sum, "c_xq2_trk", "xq2 trk eff", 1400, 600, false, true);

    process_eff_hist(h_xq2_cal.sum, h_xq2_all.sum);
    TCanvas* c_xq2_cal_eff = draw_2d_efficiency(h_xq2_cal.sum, "c_xq2_cal", "xq2 cal eff", 1400, 600, false, true);

    TH2F* h_xq2_eff_copy = (TH2F*)h_xq2_eff.sum->Clone();
    process_eff_hist(h_xq2_eff_copy, h_xq2_all.sum);
    TCanvas* c_xq2_eff_eff = draw_2d_efficiency(h_xq2_eff_copy, "c_xq2_eff_eff", "xq2 eff eff", 1400, 600, false, true);
    draw_angle(c_xq2_eff_eff, h_xq2_eff_copy, 160.0, Ee*1.0, Eh*1.0);
    draw_angle(c_xq2_eff_eff, h_xq2_eff_copy, 130.0, Ee*1.0, Eh*1.0);
    draw_angle(c_xq2_eff_eff, h_xq2_eff_copy, 40.0, Ee*1.0, Eh*1.0);

    TH2F* h_xq2_pur_copy = (TH2F*)h_xq2_pur.sum->Clone();
    process_eff_hist(h_xq2_pur_copy, h_xq2_eff.sum);
    TCanvas* c_xq2_pur_eff = draw_2d_efficiency(h_xq2_pur_copy, "c_xq2_pur_eff", "xq2 pur eff", 1400, 600, false, true);

    TH2F* h_xq2_eID_copy = (TH2F*)h_xq2_eID.sum->Clone();
    process_eff_hist(h_xq2_eID_copy, h_xq2_acp.sum);
    TCanvas* c_xq2_eID_eff = draw_2d_efficiency(h_xq2_eID_copy, "c_xq2_eID_eff", "xq2 eID eff", 1400, 600, false, true);

    TH2F* h_xq2_piID_copy = (TH2F*)h_xq2_piID.sum->Clone();
    process_eff_hist(h_xq2_piID_copy, h_xq2_acp.sum);
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
    std::string date = beam_type == 0 ? "en_26_03_1" : "ep_26_04_1";
    if ( beam_type == 3 )
    {
        c_xq2_acp_eff->SaveAs(Form("../data/%s/%s_eID_acceptance_pythia6.png", date.c_str(), setting.c_str()));
        c_xq2_eff_eff->SaveAs(Form("../data/%s/%s_eID_efficiency_pythia6.png", date.c_str(), setting.c_str()));
        c_xq2_pur_eff->SaveAs(Form("../data/%s/%s_eID_purity_pythia6.png", date.c_str(), setting.c_str()));
        c_xq2_eID_eff->SaveAs(Form("../data/%s/%s_eID_overall_pythia6.png", date.c_str(), setting.c_str()));
        c_xq2_trk_eff->SaveAs(Form("../data/%s/%s_eID_trk_pythia6.png", date.c_str(), setting.c_str()));
        c_xq2_cal_eff->SaveAs(Form("../data/%s/%s_eID_cal_pythia6.png", date.c_str(), setting.c_str()));
    }
    else
    {
        c_xq2_acp_eff->SaveAs(Form("../data/%s/%s_eID_acceptance.png", date.c_str(), setting.c_str()));
        c_xq2_eff_eff->SaveAs(Form("../data/%s/%s_eID_efficiency.png", date.c_str(), setting.c_str()));
        c_xq2_pur_eff->SaveAs(Form("../data/%s/%s_eID_purity.png", date.c_str(), setting.c_str()));
        c_xq2_eID_eff->SaveAs(Form("../data/%s/%s_eID_overall.png", date.c_str(), setting.c_str()));
        c_xq2_trk_eff->SaveAs(Form("../data/%s/%s_eID_trk.png", date.c_str(), setting.c_str()));
        c_xq2_cal_eff->SaveAs(Form("../data/%s/%s_eID_cal.png", date.c_str(), setting.c_str()));
    }
    
    // c_xq2_all->SaveAs("../data/eID/10x166_en_raw.png");

    TCanvas* c_angular = new TCanvas("c_angular", "c_angular", 1400, 600);

    TPad *draw_pad = new TPad("draw_pad", "draw_pad", 0, 0, 0.8, 1);
    draw_pad->SetFillColor(0);
    draw_pad->SetFillStyle(0);
    draw_pad->Draw();
    draw_pad->SetLeftMargin(1);
    draw_pad->SetBottomMargin(0.25);
    draw_pad->cd();

    TH1F* h_angular_formater = (TH1F*)h_angular_acp.sum->Clone("h_angular_formater");
    h_angular_formater->Reset();
    h_angular_formater->Draw();
    h_angular_formater->SetMaximum(1.05);
    h_angular_formater->SetMinimum(0);

    h_angular_formater->GetYaxis()->CenterTitle();
    h_angular_formater->GetYaxis()->SetTitleOffset(1.);

    TEfficiency*  h_acceptance_angular = new TEfficiency(*h_angular_acp.sum,*h_angular_all.sum);
    TEfficiency*  h_trk_angular = new TEfficiency(*h_angular_trk.sum,*h_angular_all.sum);
    TEfficiency*  h_cal_angular = new TEfficiency(*h_angular_cal.sum,*h_angular_all.sum);
    TEfficiency*  h_efficiency_angular = new TEfficiency(*h_angular_eff.sum,*h_angular_all.sum);
    TEfficiency*  h_purity_angular = new TEfficiency(*h_angular_pur.sum,*h_angular_eff.sum);
    TEfficiency*  h_eID_angular = new TEfficiency(*h_angular_eID.sum,*h_angular_all.sum);
    TEfficiency*  h_piID_angular = new TEfficiency(*h_angular_piID.sum,*h_angular_acp.sum);

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

    TLegend* leg = new TLegend(0.1, 0.25, 0.9, 0.65);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.08);
    leg->AddEntry(h_trk_angular, "e track eff.", "LPE");
    leg->AddEntry(h_cal_angular, "e cluster eff.", "LPE");
    // leg->AddEntry(h_acceptance_angular, "Track x Cal. eff", "LPE");
    leg->AddEntry(h_efficiency_angular, "eID efficiency", "LPE");
    leg->AddEntry(h_purity_angular, "eID purity", "LPE");
    leg->AddEntry(h_eID_angular, "Overall eID rate", "LPE");
    // leg->AddEntry(h_piID_angular, "#pi contamination", "LPE");
    leg->AddEntry(h_piID_angular, "Contamination", "LPE");
    leg->Draw();
    
    c_angular->SaveAs(Form("../data/eID/%s_eID_angular.png", setting.c_str()));

    std::string ext[5] = {"", "", "_piBG", "_beamBG", "_pythia6"};
    TFile* outFile = new TFile(Form("../data/%s/%s_eID_eff%s.root", date.c_str(), setting.c_str(), ext[beam_type].c_str()), "RECREATE");
    h_xq2_eff_copy->Write();
    h_xq2_pur_copy->Write();
    h_xq2_eID_copy->Write();
    outFile->Close();

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