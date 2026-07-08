
#include "../GlobalUtil/DrawManager.cc"
#include "../GlobalUtil/HistManager.cc"
#include "../GlobalUtil/Constants.hh"
#include "../GlobalUtil/luminosityTable.h"

#include <Math/LorentzVector.h>
using ROOT::Math::PxPyPzEVector;

#include "draw_angle.C"

using HistGroup1D = HistManager::HistGroup1D;
using HistGroup2D = HistManager::HistGroup2D;

void ReverseXAxis(TH1* h);
void draw_eta_axis(TH1* hist);
void draw_leg_pad(const std::string& type, const std::string& energy, const std::string& campaign);

void eff_pi_bg(int Ee, int Eh)
{
    // ePIC plotting style setup
    std::string type_title[3] = {"e^{3}He", "ep", "#gammap"};
    std::string energy_title = Form("%dx%d GeV", Ee, Eh);
    std::string campaign = "26.04.1";
    DrawManager* draw_manager = new DrawManager("ep + #gammap", energy_title, campaign);
    draw_manager->SetEPIC("Work in Progress");

    double text_lumi = 1; // fb^-1
    draw_manager->SetLumi(text_lumi);
    draw_manager->SetDISrange(0.01, 0.95, 4, 2);

    // histograms
    HistManager* eid_hist_manager = new HistManager();
    HistManager* bg_hist_manager = new HistManager();
     
    int n_group = 4;
    std::string setting = Form("ep_%dx%d", (int)Ee, (int)Eh);

    auto h_xq2_all = eid_hist_manager->BookHistograms(Form("h_xq2_all"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5);
    auto h_xq2_acp = eid_hist_manager->BookHistograms(Form("h_xq2_acp"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5);
    auto h_xq2_trk = eid_hist_manager->BookHistograms(Form("h_xq2_trk"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5);
    auto h_xq2_cal = eid_hist_manager->BookHistograms(Form("h_xq2_cal"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5);
    auto h_xq2_eff = eid_hist_manager->BookHistograms(Form("h_xq2_eff"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5);
    auto h_xq2_pur = eid_hist_manager->BookHistograms(Form("h_xq2_pur"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5);
    auto h_xq2_eID = eid_hist_manager->BookHistograms(Form("h_xq2_eID"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5);
    auto h_xq2_piID = eid_hist_manager->BookHistograms(Form("h_xq2_piID"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5);
    auto h_xq2_oID = eid_hist_manager->BookHistograms(Form("h_xq2_oID"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5);
    auto h_EminusPz = eid_hist_manager->MakeHistograms(Form("h_EminusPz"), "E-P_{z} all;E-P_{z} (GeV);Counts", 100, 0, 50);

    auto h_angular_all = eid_hist_manager->MakeHistograms(Form("h_angular_all"), ";#theta (deg);%", 180, 0, 180);
    auto h_angular_acp = eid_hist_manager->MakeHistograms(Form("h_angular_acp"), ";#theta (deg);%", 180, 0, 180);
    auto h_angular_trk = eid_hist_manager->MakeHistograms(Form("h_angular_trk"), ";#theta (deg);%", 180, 0, 180);
    auto h_angular_cal = eid_hist_manager->MakeHistograms(Form("h_angular_cal"), ";#theta (deg);%", 180, 0, 180);
    auto h_angular_eff = eid_hist_manager->MakeHistograms(Form("h_angular_eff"), ";#theta (deg);%", 180, 0, 180);
    auto h_angular_pur = eid_hist_manager->MakeHistograms(Form("h_angular_pur"), ";#theta (deg);%", 180, 0, 180);
    auto h_angular_eID = eid_hist_manager->MakeHistograms(Form("h_angular_eID"), ";#theta (deg);%", 180, 0, 180);
    auto h_angular_piID = eid_hist_manager->MakeHistograms(Form("h_angular_piID"), ";#theta (deg);%", 180, 0, 180);
    auto h_angular_oID = eid_hist_manager->MakeHistograms(Form("h_angular_oID"), ";#theta (deg);%", 180, 0, 180);

    for ( int i = 0; i < n_group; i ++ )
    {
        std::string date = "ep_26_03_1";
        TFile* file = new TFile(Form("../data/%s/root_files/%s_minQ2=%d_eID_pythia6_combined.root", date.c_str(), setting.c_str(), static_cast<int>(pow(10,i))));

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

            if ( i == 3 && *mc_Q2 < pow(10,i) ) 
                continue;
        
            if ( i < 3 && (*mc_Q2 > pow(10,i+1) || *mc_Q2 < pow(10,i)) ) 
                continue;

            if ( *EminusPz < 15 )
                continue;

             double theta = 180-vMCe->Theta()*(180./M_PI); // flip to make the plot go from 180 to 0 deg
               
            if ( *status >= FOUND_MC )
            {
                h_xq2_all.ind[i][1]->Fill(*mc_xB, *mc_Q2);
                h_angular_all.ind[i][1]->Fill(theta);
            }
                
            if ( *status >= FOUND_TRUTH )
            {
                h_xq2_acp.ind[i][1]->Fill(*mc_xB, *mc_Q2);
                h_angular_acp.ind[i][1]->Fill(theta);
            }

            if ( *rec_status == FOUND_TRACK_ONLY || *rec_status == FOUND_BOTH )
            {
                h_xq2_trk.ind[i][1]->Fill(*mc_xB, *mc_Q2);
                h_angular_trk.ind[i][1]->Fill(theta);
            }

            if ( *rec_status == FOUND_CLUSTER_ONLY || *rec_status == FOUND_BOTH )
            {
                h_xq2_cal.ind[i][1]->Fill(*mc_xB, *mc_Q2);
                h_angular_cal.ind[i][1]->Fill(theta);
            }

            if ( *status >= FOUND_E )
            {
                h_xq2_eff.ind[i][1]->Fill(*mc_xB, *mc_Q2);   
                h_angular_eff.ind[i][1]->Fill(theta);
                h_EminusPz.ind[i][1]->Fill(*EminusPz);
            }

            if ( *status == FOUND_E )
            {
                h_xq2_pur.ind[i][1]->Fill(*mc_xB, *mc_Q2);
                h_xq2_eID.ind[i][1]->Fill(*mc_xB, *mc_Q2);
                h_angular_pur.ind[i][1]->Fill(theta);
                h_angular_eID.ind[i][1]->Fill(theta);
            }

            if ( *status == FOUND_PI )
            {
                h_xq2_piID.ind[i][1]->Fill(*mc_xB, *mc_Q2);
                h_angular_piID.ind[i][1]->Fill(theta);
            }

            if ( *status == FOUND_OTHERS )
            {
                h_xq2_oID.ind[i][1]->Fill(*mc_xB, *mc_Q2);
                h_angular_oID.ind[i][1]->Fill(theta);
            }
        }

        // get generated lumi
        double total_lumi = text_lumi; // fb^-1
        double gen_lumi = get_lumi(4, Ee, Eh, i, nEntries, 0);
        if ( gen_lumi == 0 )
        {
            std::cout << "** Lumi table not found! " << std::endl;
            return;
        }

        eid_hist_manager->SumHistograms(i, 1, total_lumi / gen_lumi);
        file->Close();
    }

    std::cout << "Finished processing signal files." << std::endl;

    TFile* pi_bg_file = new TFile(Form("../data/ep_26_03_1/root_files/piBG_%dx%d_eid_combined.root", (int)Ee, (int)Eh));
    if ( pi_bg_file == nullptr || pi_bg_file->IsZombie() ) {
        std::cerr << "Error: Could not open pi bg file!" << std::endl;
        return;
    }

    auto h_xq2_all_bg = bg_hist_manager->BookHistograms(Form("h_xq2_all_bg"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5);
    auto h_xq2_acp_bg = bg_hist_manager->BookHistograms(Form("h_xq2_acp_bg"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5);
    auto h_xq2_trk_bg = bg_hist_manager->BookHistograms(Form("h_xq2_trk_bg"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5);
    auto h_xq2_cal_bg = bg_hist_manager->BookHistograms(Form("h_xq2_cal_bg"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5);
    auto h_xq2_eff_bg = bg_hist_manager->BookHistograms(Form("h_xq2_eff_bg"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5);
    auto h_xq2_pur_bg = bg_hist_manager->BookHistograms(Form("h_xq2_pur_bg"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5);
    auto h_xq2_eID_bg = bg_hist_manager->BookHistograms(Form("h_xq2_eID_bg"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5);
    auto h_xq2_piID_bg = bg_hist_manager->BookHistograms(Form("h_xq2_piID_bg"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5);
    auto h_xq2_oID_bg = bg_hist_manager->BookHistograms(Form("h_xq2_oID_bg"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5);
    auto h_EminusPz_bg = bg_hist_manager->MakeHistograms(Form("h_EminusPz_bg"), "E-P_{z} all;E-P_{z} (GeV);Counts", 100, 0, 50);

    auto h_angular_trk_bg = bg_hist_manager->MakeHistograms(Form("h_angular_trk_bg"), "Angular distribution of tracks;#theta (deg);Counts", 180, 0, 180);
    auto h_angular_cal_bg = bg_hist_manager->MakeHistograms(Form("h_angular_cal_bg"), "Angular distribution of calorimeter;#theta (deg);Counts", 180, 0, 180);
    auto h_angular_all_bg = bg_hist_manager->MakeHistograms(Form("h_angular_all_bg"), "Angular distribution of all;#theta (deg);Counts", 180, 0, 180);
    auto h_angular_acp_bg = bg_hist_manager->MakeHistograms(Form("h_angular_acp_bg"), "Angular distribution of accepted;#theta (deg);Counts", 180, 0, 180);
    auto h_angular_eff_bg = bg_hist_manager->MakeHistograms(Form("h_angular_eff_bg"), "Angular distribution of efficiency;#theta (deg);Counts", 180, 0, 180);
    auto h_angular_pur_bg = bg_hist_manager->MakeHistograms(Form("h_angular_pur_bg"), "Angular distribution of purity;#theta (deg);Counts", 180, 0, 180);
    auto h_angular_eID_bg = bg_hist_manager->MakeHistograms(Form("h_angular_eID_bg"), "Angular distribution of eID;#theta (deg);Counts", 180, 0, 180);
    auto h_angular_piID_bg = bg_hist_manager->MakeHistograms(Form("h_angular_piID_bg"), "Angular distribution of piID;#theta (deg);Counts", 180, 0, 180);
    auto h_angular_oID_bg = bg_hist_manager->MakeHistograms(Form("h_angular_oID_bg"), "Angular distribution of oID;#theta (deg);Counts", 180, 0, 180);

    TTreeReader reader("T_eID", pi_bg_file);

    TTreeReaderValue<int> status(reader, "eID_status");
    TTreeReaderValue<int> rec_status(reader, "eRecon_status");

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
        
        if ( *EminusPz < 15 )
            continue;

        double theta = 180-vMCe->Theta()*(180./M_PI); // flip to make the plot go from 180 to 0 deg

        if ( *status >= FOUND_MC )
        {
            h_xq2_all_bg.ind[0][1]->Fill(*mc_xB, *mc_Q2);
            h_angular_all_bg.ind[0][1]->Fill(theta);
        }
            
        if ( *status >= FOUND_TRUTH )
        {
            h_xq2_acp_bg.ind[0][1]->Fill(*mc_xB, *mc_Q2);
            h_angular_acp_bg.ind[0][1]->Fill(theta);
        }

        if ( *rec_status == FOUND_TRACK_ONLY || *rec_status == FOUND_BOTH )
        {
            h_xq2_trk_bg.ind[0][1]->Fill(*mc_xB, *mc_Q2);
            h_angular_trk_bg.ind[0][1]->Fill(theta);
        }

        if ( *rec_status == FOUND_CLUSTER_ONLY || *rec_status == FOUND_BOTH )
            {
                h_xq2_cal_bg.ind[0][1]->Fill(*mc_xB, *mc_Q2);
                h_angular_cal_bg.ind[0][1]->Fill(theta);
            }

        if ( *status >= FOUND_E )
        {
            h_xq2_eff_bg.ind[0][1]->Fill(*rec_xB, *rec_Q2);   
            h_angular_eff_bg.ind[0][1]->Fill(theta);
        }

        if ( *status == FOUND_E )
        {
            // std::cout << "Filling pi bg event: rec_xB = " << *rec_xB << ", rec_Q2 = " << *rec_Q2 << std::endl;
            h_xq2_pur_bg.ind[0][1]->Fill(*rec_xB, *rec_Q2);
            h_xq2_eID_bg.ind[0][1]->Fill(*rec_xB, *rec_Q2);
            h_angular_pur_bg.ind[0][1]->Fill(theta);
            h_angular_eID_bg.ind[0][1]->Fill(theta);
        }

        if ( *status == FOUND_PI )
        {
            h_xq2_piID_bg.ind[0][1]->Fill(*rec_xB, *rec_Q2);
            h_EminusPz_bg.ind[0][1]->Fill(*EminusPz);
            h_angular_piID_bg.ind[0][1]->Fill(theta);
        }
            
        if ( *status == FOUND_OTHERS )
        {
            h_xq2_oID_bg.ind[0][1]->Fill(*rec_xB, *rec_Q2);
            h_angular_oID_bg.ind[0][1]->Fill(theta);
        }
    }
    
    double total_lumi = text_lumi; // fb^-1
    double gen_lumi = get_lumi(2, Ee, Eh, 0, nEntries, 0);
    if ( gen_lumi == 0 )
    {
        std::cout << "** Lumi table not found! " << std::endl;
        return;
    }
    bg_hist_manager->SumHistograms(0, 1, total_lumi / gen_lumi);
    
    std::cout << "Entries in h_xq2_all_bg: " << h_xq2_all_bg.sum->GetEntries() << std::endl;
    // std::cout << "Entries in h_EminusPz_bg: " << h_EminusPz_bg->GetEntries() << std::endl;
    // pi_bg_file->Close();

    std::cout << "Finished processing pi bg file." << std::endl;

    // o0o0o0o0o

    set_2d_scale(h_xq2_all.sum);
    TCanvas* c_xq2_all = draw_2d_standard(h_xq2_all.sum, "c_xq2_all", "all events", 700, 600, true, true);
    TCanvas* c_xq2_acp = draw_2d_standard(h_xq2_acp.sum, "c_xq2_acp", "acp events", 700, 600, true, true);
    TCanvas* c_xq2_eID = draw_2d_standard(h_xq2_eID.sum, "c_xq2_eID", "eID events", 700, 600, true, true);
    TCanvas* c_xq2_piID = draw_2d_standard(h_xq2_piID.sum, "c_xq2_piID", "piID events", 700, 600, true, true);

    TCanvas* c_xq2_eff_bg = draw_2d_standard(h_xq2_eff_bg.sum, "c_xq2_eff_bg", "PI bg events", 700, 600, true, true);
    TCanvas* c_xq2_piID_bg = draw_2d_standard(h_xq2_piID_bg.sum, "c_xq2_piID_bg", "piID bg events", 700, 600, true, true);
    TCanvas* c_xq2_oID_bg = draw_2d_standard(h_xq2_oID_bg.sum, "c_xq2_oID_bg", "oID bg events", 700, 600, true, true);

    TH2F* h_xq2_acp_copy = (TH2F*)h_xq2_acp.sum->Clone();
    process_eff_hist(h_xq2_acp_copy, h_xq2_all.sum);
    TCanvas* c_xq2_acp_eff = draw_2d_efficiency(h_xq2_acp_copy, "c_xq2_acp_eff", "xq2 acp eff", 1400, 600, false, true);

    TH2F* h_xq2_eff_copy = (TH2F*)h_xq2_eff.sum->Clone();
    // TH2F* h_xq2_eff_copy = (TH2F*)h_xq2_acp->Clone();
    h_xq2_eff_copy->Add(h_xq2_eff_bg.sum);
    // TCanvas* c_xq2_eff_copy = draw_2d_standard(h_xq2_eff_copy, "c_xq2_eff_copy", "eff copy events", 700, 600, true, true);
    // TH2F* h_xq2_eff_bg_copy = (TH2F*)h_xq2_eff_bg->Clone();
    // process_eff_hist(h_xq2_eff_bg_copy, h_xq2_eff_copy);
    // TCanvas* c_xq2_contaimin = draw_2d_efficiency(h_xq2_eff_bg_copy, "c_xq2_contaimin", "xq2 contamination", 1400, 600, false, true);
    TH2F* h_xq2_all_copy = (TH2F*)h_xq2_all.sum->Clone();
    h_xq2_all_copy->Add(h_xq2_all_bg.sum);
    process_eff_hist(h_xq2_eff_copy, h_xq2_all_copy);
    TCanvas* c_xq2_eff_eff = draw_2d_efficiency(h_xq2_eff_copy, "c_xq2_eff_eff", "xq2 eff eff", 1400, 600, false, true);
    // draw_angle(c_xq2_eff_eff, h_xq2_eff_copy, 160.0, Ee*1.0, Eh*1.0);
    // draw_angle(c_xq2_eff_eff, h_xq2_eff_copy, 130.0, Ee*1.0, Eh*1.0);
    // draw_angle(c_xq2_eff_eff, h_xq2_eff_copy, 40.0, Ee*1.0, Eh*1.0);

    TH2F* h_xq2_eff_copy2 = (TH2F*)h_xq2_eff.sum->Clone();
    // TH2F* h_xq2_eff_copy2 = (TH2F*)h_xq2_acp->Clone();
    h_xq2_eff_copy2->Add(h_xq2_eff_bg.sum);
    TH2F* h_xq2_pur_copy = (TH2F*)h_xq2_pur.sum->Clone();
    process_eff_hist(h_xq2_pur_copy, h_xq2_eff_copy2);
    TCanvas* c_xq2_pur_eff = draw_2d_efficiency(h_xq2_pur_copy, "c_xq2_pur_eff", "xq2 pur eff", 1400, 600, false, true);

    TH2F* h_xq2_eID_copy = (TH2F*)h_xq2_eID.sum->Clone();
    process_eff_hist(h_xq2_eID_copy, h_xq2_eff_copy2);
    TCanvas* c_xq2_eID_eff = draw_2d_efficiency(h_xq2_eID_copy, "c_xq2_eID_eff", "xq2 eID eff", 1400, 600, false, true);

    TH2F* h_xq2_piID_copy = (TH2F*)h_xq2_piID.sum->Clone();
    h_xq2_piID_copy->Add(h_xq2_piID_bg.sum);
    process_eff_hist(h_xq2_piID_copy, h_xq2_eff_copy2);
    TCanvas* c_xq2_piID_eff = draw_2d_efficiency(h_xq2_piID_copy, "c_xq2_piID_eff", "xq2 piID eff", 1400, 600, false, true);

    TH2F* h_xq2_oID_copy = (TH2F*)h_xq2_oID.sum->Clone();
    h_xq2_oID_copy->Add(h_xq2_oID_bg.sum);
    h_xq2_oID_copy->Add(h_xq2_piID_bg.sum);
    process_eff_hist(h_xq2_oID_copy, h_xq2_eff_copy2);
    TCanvas* c_xq2_oID_eff = draw_2d_efficiency(h_xq2_oID_copy, "c_xq2_oID_eff", "xq2 oID eff", 1400, 600, false, true);

    TCanvas* c_EminusPz = new TCanvas("c_EminusPz", "c_EminusPz", 1000, 600);
    h_EminusPz.sum->SetLineColor(kBlue);
    h_EminusPz.sum->SetFillColor(kBlue);
    h_EminusPz.sum->SetFillStyle(3003);
    h_EminusPz.sum->SetLineWidth(2);
    h_EminusPz.sum->GetXaxis()->SetRangeUser(0, 50);
    double draw_max = h_EminusPz.sum->GetMaximum() * 1.2;
    h_EminusPz.sum->SetMaximum(draw_max);
    h_EminusPz.sum->Draw("HIST");
    h_EminusPz_bg.sum->SetLineColor(kRed);
    h_EminusPz_bg.sum->SetFillColor(kRed);
    h_EminusPz_bg.sum->SetFillStyle(3003);
    h_EminusPz_bg.sum->SetLineWidth(2);
    h_EminusPz_bg.sum->Draw("HIST SAME");

    TLegend* leg = new TLegend(0.8, 0.7, 0.95, 0.85);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(h_EminusPz.sum, "e p", "l");
    leg->AddEntry(h_EminusPz_bg.sum, "#gamma p", "l");
    leg->Draw("SAME");

    c_EminusPz->Modified();
    c_EminusPz->Update();

    TLine* line = new TLine(2*Ee, 0, 2*Ee, draw_max);
    line->SetLineColor(kBlack);
    line->SetLineStyle(7);
    line->Draw("SAME");

    std::string date = "ep_26_03_1";

    draw_manager->LableAndCollect(c_xq2_piID_bg);
    c_xq2_piID_bg->SetFrameLineWidth(0);
    std::cout << "N events in c_xq2_piID_bg: " << h_xq2_piID_bg.sum->Integral() << std::endl;
    // c_xq2_piID_bg->SaveAs(Form("../data/eID/%s_eID_xq2_piID_bg_wCut.png", setting.c_str()));
    c_xq2_piID_bg->SaveAs(Form("../data/%s/%s_eID_xq2_piID_bg.png", date.c_str(), setting.c_str()));

    draw_manager->LableAndCollect(c_EminusPz);
    c_EminusPz->SaveAs(Form("../data/%s/%s_eID_EminusPz_PiBG.png", date.c_str(), setting.c_str()));

    draw_manager->LableAndCollect(c_xq2_eff_eff);
    c_xq2_eff_eff->SaveAs(Form("../data/%s/%s_eID_eff_PiBG.png", date.c_str(), setting.c_str()));

    draw_manager->LableAndCollect(c_xq2_pur_eff);
    c_xq2_pur_eff->SaveAs(Form("../data/%s/%s_eID_pur_PiBG.png", date.c_str(), setting.c_str()));

    // 00o0o0o0o

    h_angular_all.sum->Add(h_angular_all_bg.sum);
    h_angular_acp.sum->Add(h_angular_acp_bg.sum);
    h_angular_eff.sum->Add(h_angular_eff_bg.sum);
    h_angular_pur.sum->Add(h_angular_pur_bg.sum);
    h_angular_eID.sum->Add(h_angular_eID_bg.sum);
    h_angular_piID.sum->Add(h_angular_piID_bg.sum);
    h_angular_oID.sum->Add(h_angular_oID_bg.sum);
    h_angular_trk.sum->Add(h_angular_trk_bg.sum);
    h_angular_cal.sum->Add(h_angular_cal_bg.sum);

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
    TEfficiency*  h_piID_angular = new TEfficiency(*h_angular_piID.sum,*h_angular_all.sum);
    TEfficiency*  h_oID_angular = new TEfficiency(*h_angular_oID.sum,*h_angular_all.sum);

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

    h_oID_angular->SetLineColor(kCyan-6);
    h_oID_angular->SetMarkerColor(kCyan-6);

    h_acceptance_angular->SetStatisticOption(TEfficiency::kFNormal);
    h_efficiency_angular->SetStatisticOption(TEfficiency::kFNormal);
    h_trk_angular->SetStatisticOption(TEfficiency::kFNormal);
    h_cal_angular->SetStatisticOption(TEfficiency::kFNormal);
    h_purity_angular->SetStatisticOption(TEfficiency::kFNormal);
    h_eID_angular->SetStatisticOption(TEfficiency::kFNormal);
    h_piID_angular->SetStatisticOption(TEfficiency::kFNormal);
    h_oID_angular->SetStatisticOption(TEfficiency::kFNormal);

    TLine* unity_line = new TLine(h_angular_formater->GetXaxis()->GetXmin(), 1.0, h_angular_formater->GetXaxis()->GetXmax(), 1.0);
    unity_line->SetLineColor(kGray+2);
    unity_line->SetLineStyle(2);
    unity_line->SetLineWidth(2);
    unity_line->Draw("SAME");

    // h_acceptance_angular->Draw("SAME");
    h_efficiency_angular->Draw("SAME");
    h_trk_angular->Draw("SAME");
    h_cal_angular->Draw("SAME");
    h_purity_angular->Draw("SAME ");
    h_eID_angular->Draw("SAME");
    h_piID_angular->Draw("SAME");
    h_oID_angular->Draw("SAME");

    ReverseXAxis(h_angular_formater);
    draw_eta_axis(h_angular_formater);

    c_angular->cd();
    draw_leg_pad(type_title[1], energy_title, campaign);

    TLegend* leg_angle = new TLegend(0.1, 0.3, 0.9, 0.65);
    leg_angle->SetFillColor(0);
    leg_angle->SetBorderSize(0);
    leg_angle->SetTextSize(0.08);
    leg_angle->AddEntry(h_trk_angular, "e track eff.", "LPE");
    leg_angle->AddEntry(h_cal_angular, "e cluster eff.", "LPE");
    // leg_angle->AddEntry(h_acceptance_angular, "Track x Cal. eff", "LPE");
    leg_angle->AddEntry(h_efficiency_angular, "eID efficiency", "LPE");
    leg_angle->AddEntry(h_purity_angular, "eID purity", "LPE");
    leg_angle->AddEntry(h_eID_angular, "Overall eID rate", "LPE");
    leg_angle->AddEntry(h_piID_angular, "#pi contamination", "LPE");
    leg_angle->AddEntry(h_oID_angular, "Others", "LPE");
    leg_angle->Draw();

    c_angular->Modified();
    c_angular->Update();
    c_angular->SaveAs(Form("../data/%s/%s_eID_angular_eff_PiBG.png", date.c_str(), setting.c_str()));

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