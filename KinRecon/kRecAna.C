#include "kRecAna.h"

double Ee = 10;
double Eh = 166;
const std::string group[3] = {"1to10", "10to100", "100to1000"};
const double total_lumi = 8.65*3; // fb^-1
int e_set[3] = {1, 10, 100};
const int n_group = 3;

double cs[3][2] = {{0.198440424611563, 0.205327493968226}, {4.04371412707044E-02, 4.41976212963417E-02}, {1.36416909784756E-03, 1.69583242740138E-03}};
double ev[3][2] = {{333675, 666325}, {333694, 666306}, {333365, 666640}};

void kRecAna()
{
    std::string type_title = "e^{3}He";
    std::string energy_title = Form("%dx%d GeV", 10, 166);
    DrawManager* draw_manager = new DrawManager(type_title, energy_title, "25.10.2");
    draw_manager->SetEPIC();
    draw_manager->SetLumi(8.65);

    define_histograms();

    vector<Kinematics> algorithm[n_group];

    for ( int i = 0; i < n_group; i ++ )
    {
        // load files and trees
        TFile* beamFile = new TFile(Form("../data/en_25_10_2/root_files/10x166_%s_beam_combined.root", group[i].c_str()));
        TFile* eidFile  = new TFile(Form("../data/en_25_10_2/root_files/10x166_%s_eIDrecon_combined.root", group[i].c_str()));
        
        TTreeReader beam_reader("T_Beam", beamFile);
        TTreeReaderValue<int> N_PDG(beam_reader, "N_PDG");
        TTreeReaderValue<PxPyPzEVector> beam_e(beam_reader, "vectE");
        TTreeReaderValue<PxPyPzEVector> beam_n(beam_reader, "vectN");

        TTreeReader eid_reader("T_eID", eidFile);
        TTreeReaderValue<int>    status(eid_reader, "eID_status");
        TTreeReaderValue<double> mc_xB(eid_reader, "mc_xB");
        TTreeReaderValue<double> mc_Q2(eid_reader, "mc_Q2");
        TTreeReaderValue<double> mc_y(eid_reader, "mc_y");
        TTreeReaderValue<double> mc_W2(eid_reader, "mc_W2");
        TTreeReaderValue<double> mc_nu(eid_reader, "mc_nu");

        TTreeReaderValue<PxPyPzEVector> vMCe(eid_reader, "vMC_e");
        TTreeReaderValue<PxPyPzEVector> vTRe(eid_reader, "vTRACK_e");
        TTreeReaderValue<PxPyPzEVector> vCLe(eid_reader, "vCLUSTER_e");
        TTreeReaderValue<std::vector<PxPyPzEVector>> vMCf(eid_reader, "vMC_hfs");
        TTreeReaderValue<std::vector<PxPyPzEVector>> vTRf(eid_reader, "vREC_hfs");

        // setup for scaling
        double gen_lumi = ev[i][0]/(cs[i][0]*(1e-34/1e-43)) + ev[i][1]/(cs[i][1]*(1e-34/1e-43));
        
        algorithm[i].push_back(Kinematics("E_Method", 2)); // electron method
        algorithm[i].push_back(Kinematics("JB_Method", 4)); // Jacquet-Blondel method
        algorithm[i].push_back(Kinematics("DA_Method", 1)); // Double-Angle method
        algorithm[i].push_back(Kinematics("Sig_Method", kOrange+1)); // Sigma method
        algorithm[i].push_back(Kinematics("ESig_Method", kGreen+3)); // E-Sigma method
        algorithm[i].push_back(Kinematics("MC_Method", kGray)); // MC information

        Long64_t nEntries = beam_reader.GetEntries();
        for ( int ev = 0; ev < nEntries; ev ++ )
        {
            beam_reader.Next();
            eid_reader.Next();
            
            if ( *status < FOUND_E || *N_PDG == 0 )
                continue;

            if (*mc_y < 0.01 || *mc_y > 0.95) 
                continue;

            if (*mc_W2 < 4)
                continue;

            if (*mc_Q2 < 2)
                continue;
            
            LorentzRotation boost = *N_PDG == ID_PROTON ? getBoost( Ee, Eh, MASS_ELECTRON, MASS_PROTON) : getBoost( Ee, Eh, MASS_ELECTRON, MASS_NEUTRON);
            
            PxPyPzEVector boost_beam_e = boost(*beam_e);
            PxPyPzEVector boost_beam_n = boost(*beam_n);
            
            PxPyPzEVector boosted_vMCe = boost(*vMCe);
            PxPyPzEVector boosted_vTRe = boost(*vTRe);
            PxPyPzEVector boosted_vCLe = boost(*vCLe);

            double pxsum = 0;
            double pysum = 0;
            double pzsum = 0;
            double Esum  = 0;
            for ( const auto& p : *vMCf )
            {
                PxPyPzEVector boosted_p = boost(p);
                pxsum += boosted_p.Px();
                pysum += boosted_p.Py();
                pzsum += boosted_p.Pz();
                Esum  += boosted_p.E();
            }
            double mc_pt_had = sqrt(pxsum*pxsum + pysum*pysum);
            double mc_sigma_h = Esum - pzsum;  

            pxsum = 0;
            pysum = 0;
            pzsum = 0;
            Esum  = 0;
            for ( const auto& p : *vTRf )
            {
                PxPyPzEVector boosted_p = boost(p);
                pxsum += boosted_p.Px();
                pysum += boosted_p.Py();
                pzsum += boosted_p.Pz();
                Esum  += boosted_p.E();
            }
            double pt_had = sqrt(pxsum*pxsum + pysum*pysum);
            double sigma_h = Esum - pzsum;

            double Mp = *N_PDG == ID_PROTON ? MASS_PROTON : MASS_NEUTRON;
            algorithm[i][MC].save_variables(calc_elec_method(boosted_vMCe.E(), boosted_vMCe.theta(), mc_pt_had, mc_sigma_h, boost_beam_e.E(), boost_beam_n.E()), Mp);
            algorithm[i][EL].save_variables(calc_elec_method(boosted_vTRe.E(), boosted_vTRe.theta(), pt_had, sigma_h, Ee, Eh), Mp);

            if ( algorithm[i][EL].Q2 < 10 || algorithm[i][EL].Q2 >= 1000 )
                algorithm[i][EL].save_variables(calc_elec_method(boosted_vCLe.E(), boosted_vCLe.theta(), pt_had, sigma_h, Ee, Eh), Mp);
        
            if ( vTRf->size() != 0 ) 
            {
                algorithm[i][JB].save_variables(calc_jb_method(boosted_vTRe.E(), boosted_vTRe.theta(), pt_had, sigma_h, Ee, Eh), Mp);
                algorithm[i][DA].save_variables(calc_da_method(boosted_vTRe.E(), boosted_vTRe.theta(), pt_had, sigma_h, Ee, Eh), Mp);
                algorithm[i][SIG].save_variables(calc_sig_method(boosted_vTRe.E(), boosted_vTRe.theta(), pt_had, sigma_h, Ee, Eh), Mp);
                algorithm[i][ESIG].save_variables(calc_esig_method(boosted_vTRe.E(), boosted_vTRe.theta(), pt_had, sigma_h, Ee, Eh), Mp);
            }

            for ( int j = 0; j < 5; j ++ )
                algorithm[i][j].Process_and_QA(algorithm[i][MC].xB, algorithm[i][MC].y, algorithm[i][MC].Q2);

            int Qrange = 0;
            if ( algorithm[i][MC].Q2 > 1000 )
                Qrange = 3;
            else if ( algorithm[i][MC].Q2 > 100 )
                Qrange = 2;
            else if ( algorithm[i][MC].Q2 > 10 )
                Qrange = 1;

           for ( int j = 0; j < 6; j ++ )
                algorithm[i][j].reset();

            h_dEt[Qrange]->Fill((vTRe->E()-vMCe->E())/vMCe->E());
            h_dPt[Qrange]->Fill((vTRe->P()-vMCe->P())/vMCe->P());
            h_dt[Qrange]->Fill((vTRe->Theta()-vMCe->Theta())/vMCe->Theta());
            h_dp[Qrange]->Fill((vTRe->Phi()-vMCe->Phi())/vMCe->Phi());

            // h_tde[1]->Fill(vMCe->Theta()*(180./M_PI), (vTRe->E()-vMCe->E())/vMCe->E());
            // h_tdp[1]->Fill(vMCe->Theta()*(180./M_PI), (vTRe->P()-vMCe->P())/vMCe->P());
            h_tde_tmp[1][Qrange]->Fill(vMCe->Theta()*(180./M_PI), (vTRe->E()-vMCe->E())/vMCe->E());
            h_tdp_tmp[1][Qrange]->Fill(vMCe->Theta()*(180./M_PI), (vTRe->P()-vMCe->P())/vMCe->P());

            h_tdphi->Fill(vMCe->Theta()*(180./M_PI), (vTRe->Phi()-vMCe->Phi())/vMCe->Phi());
            h_pde->Fill(vMCe->Pt(), (vTRe->E()-vMCe->E())/vMCe->E());
            h_pt->Fill(vMCe->Theta()*(180./M_PI),vMCe->Pt());

            h_dEc[Qrange]->Fill((vCLe->E()-vMCe->E())/vMCe->E());
            h_dPc[Qrange]->Fill((vCLe->P()-vMCe->P())/vMCe->P());
            // h_tde[0]->Fill(vMCe->Theta()*(180./M_PI), (vCLe->E()-vMCe->E())/vMCe->E());
            // h_tdp[0]->Fill(vMCe->Theta()*(180./M_PI), (vCLe->P()-vMCe->P())/vMCe->P());

            h_tde_tmp[0][Qrange]->Fill(vMCe->Theta()*(180./M_PI), (vTRe->E()-vMCe->E())/vMCe->E());
            h_tdp_tmp[0][Qrange]->Fill(vMCe->Theta()*(180./M_PI), (vTRe->P()-vMCe->P())/vMCe->P());

            h_ede->Fill(vMCe->E(), (vCLe->E()-vMCe->E())/vMCe->E());
            h_et->Fill(vMCe->Theta()*(180./M_PI),vMCe->E());
        }

        for ( int j = 0; j < 6; j ++ )
        {
            algorithm[i][j].h_qvq->Scale( total_lumi / gen_lumi );
            algorithm[i][j].h_xvx->Scale( total_lumi / gen_lumi );
            algorithm[i][j].h_yvy->Scale( total_lumi / gen_lumi );

            algorithm[i][j].h_xq->Scale( total_lumi / gen_lumi );
            algorithm[i][j].h_xq_eff->Scale( total_lumi / gen_lumi );
        }

        for ( int i = 0; i < 2; i ++ )
            for ( int j = 0; j < 4; j ++ )
            {
                h_tde_tmp[i][j]->Scale( total_lumi / gen_lumi );
                h_tdp_tmp[i][j]->Scale( total_lumi / gen_lumi );

                h_tde[i]->Add(h_tde_tmp[i][j]);
                h_tdp[i]->Add(h_tdp_tmp[i][j]);
            }
        

        // beamFile->Close();
        // eidFile->Close();
    }

    const int n_methods = algorithm[0].size()-1;

    TH2F* h_qvq[n_methods]; 
    TH2F* h_xvx[n_methods]; 
    TH2F* h_yvy[n_methods];

    TCanvas* c_reco_v_true = new TCanvas("c_reco_v_true", "RECO v TRUE", 1800, 1000);
    c_reco_v_true->Divide(n_methods, 3);

    for ( int i = 0; i < n_methods; i ++ )
    {
        h_qvq[i] = BookTH2(Form("H_%s_QVQ_all",algorithm[0][i].method_name.c_str()), Form("%s;Q_{true}^{2} (GeV/c^{2})^{2};Q_{reco}^{2} (GeV/c^{2})^{2}", algorithm[0][i].method_name.c_str()), 200, 0, 4, 200,  0, 4, kLightTemperature);
        h_xvx[i] = BookTH2(Form("H_%s_XVX_all", algorithm[0][i].method_name.c_str()), Form("%s;x_{true};x_{reco}", algorithm[0][i].method_name.c_str()), 100, -4, 0, 100,  -4, 0, kLightTemperature);
        h_yvy[i] = new TH2F(Form("H_%s_YVY_all", algorithm[0][i].method_name.c_str()), Form("%s;y_{true};y_{reco}", algorithm[0][i].method_name.c_str()), 100, 0, 1, 100,  0, 1);

        h_xq_gen[i] = BookTH2(Form("H_%s_Events_all", algorithm[0][i].method_name.c_str()), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
        h_xq_eff[i] = BookTH2(Form("H_%s_Efficiency_all", algorithm[0][i].method_name.c_str()), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);

        algorithm[0][i].StyleHist(h_qvq[i]);
        algorithm[0][i].StyleHist(h_xvx[i]);
        algorithm[0][i].StyleHist(h_yvy[i]);
        algorithm[0][i].StyleHist(h_xq_gen[i]);
        algorithm[0][i].StyleHist(h_xq_eff[i]);

        for ( int r = 0; r < n_group; r ++ )
        {
            h_qvq[i]->Add(algorithm[r][i].h_qvq);
            h_xvx[i]->Add(algorithm[r][i].h_xvx);
            h_yvy[i]->Add(algorithm[r][i].h_yvy);
            h_xq_gen[i]->Add(algorithm[r][i].h_xq);
            h_xq_eff[i]->Add(algorithm[r][i].h_xq_eff);
        }

        c_reco_v_true->cd(i+1);
        h_qvq[i]->Draw("COLZ");

        c_reco_v_true->cd(i+1+n_methods);
        h_xvx[i]->Draw("COLZ");

        c_reco_v_true->cd(i+1+2*n_methods);
        h_yvy[i]->Draw("COLZ");
    }

    for ( int i = 0; i < n_methods*3; i ++ )
    {
        c_reco_v_true->cd(i+1);
        gPad->SetLogz();
        if ( i < n_methods*2 )
        {
            gPad->SetLogx();
            gPad->SetLogy();
        } 
        gPad->SetTickx(1);
        gPad->SetTicky(1);
        gPad->SetLeftMargin(0.2);
        gPad->SetBottomMargin(0.2);
    }

    c_reco_v_true->cd(1);
    gPad->Update();
	TLine* l_q_diagonal = new TLine(1, 1, 1e4, 1e4);
	l_q_diagonal->SetLineColor(kBlack);
	l_q_diagonal->SetLineStyle(7);

    c_reco_v_true->cd(1+n_methods);
    gPad->Update();
	TLine* l_x_diagonal = new TLine(1e-4, 1e-4, 1, 1);
	l_x_diagonal->SetLineColor(kBlack);
	l_x_diagonal->SetLineStyle(7);

    c_reco_v_true->cd(1+2*n_methods);
    gPad->Update();
	TLine* l_y_diagonal = new TLine(0, 0, 1, 1);
	l_y_diagonal->SetLineColor(kBlack);
	l_y_diagonal->SetLineStyle(7);
    
    for ( int i = 0; i < n_methods; i ++ )
    {
        c_reco_v_true->cd(i+1);
        l_q_diagonal->Draw("SAME");

        c_reco_v_true->cd(i+1+n_methods);
        l_x_diagonal->Draw("SAME");

        c_reco_v_true->cd(i+1+2*n_methods);
        l_y_diagonal->Draw("SAME");
    }

    draw_manager->LableAndCollectSpecial(c_reco_v_true);
    c_reco_v_true->SaveAs(Form("../data/en_25_10_2/reco_v_true_%dx%d.png", 10, 166));

    TCanvas* c_energy_track_qa1;
    TCanvas* c_energy_track_qa2;
    TCanvas* c_energy_track_qa3;
    TCanvas* c_energy_track_qa4;
    plot_energy_and_track(c_energy_track_qa1, c_energy_track_qa2, c_energy_track_qa3, c_energy_track_qa4);
    draw_manager->LableAndCollectSpecial(c_energy_track_qa1);
    // c_energy_track_qa1->SaveAs(Form("../data/en_25_10_2/track_vs_cluster_%dx%d.png", 10, 166));

    draw_manager->LableAndCollect(c_energy_track_qa2);
    // c_energy_track_qa2->SaveAs(Form("../data/en_25_10_2/energy_and_theta_qa_%dx%d.png", 10, 166));

    set_2d_scale(h_xq_gen[0]);
    TCanvas* c_eRecon_eff[n_methods];
    for ( int i = 0; i < n_methods; i ++ )
    {
        process_eff_hist(h_xq_eff[i], h_xq_gen[i]);
        c_eRecon_eff[i] = draw_2d_efficiency(h_xq_eff[i], Form("c_%s_eff", algorithm[0][i].method_name.c_str()), Form("%s Efficiency", algorithm[0][i].method_name.c_str()), 1400, 600, false, true);
    }

    TFile*outFile = new TFile("eHe3_10x166_ReconStats.root", "RECREATE");
    outFile->cd();
    for ( int i = 0; i < n_methods; i ++ )
        c_eRecon_eff[i]->Write();
}

void define_histograms()
{
    std::string QrangeName[4] = {"lowQ", "midQ", "highQ", "ultraQ"};
    std::string Qrange[4] = {"Q^{2} #leq 10 (GeV/c^{2})^{2}", "10 (GeV/c^{2})^{2} < Q^{2} #leq 100 (GeV/c^{2})^{2}", "100 (GeV/c^{2})^{2} < Q^{2} #leq 1000 (GeV/c^{2})^{2}", "Q^{2} > 1000 (GeV/c^{2})^{2}"};
    for ( int i = 0; i < 4; i ++ )
    {
        h_dEc[i] = new TH1F(Form("H_%s_dE_cal",  QrangeName[i].c_str()), Form("%s;(Ee_{reco}-Ee_{true})/Ee_{true};Counts", Qrange[i].c_str()), 100, -1, 1);
        h_dEt[i] = new TH1F(Form("H_%s_dE_track",  QrangeName[i].c_str()), Form("%s;(Ee_{reco}-Ee_{true})/Ee_{true};Counts", Qrange[i].c_str()), 100, -1, 1);
        h_dPc[i] = new TH1F(Form("H_%s_dP_cal",  QrangeName[i].c_str()), Form("%s;(Pe_{reco}-Pe_{true})/Pe_{true};Counts", Qrange[i].c_str()), 100, -1, 1);
        h_dPt[i] = new TH1F(Form("H_%s_dP_track",  QrangeName[i].c_str()), Form("%s;(Pe_{reco}-Pe_{true})/Pe_{true};Counts", Qrange[i].c_str()), 100, -1, 1);
        h_dt[i] = new TH1F(Form("H_%s_dtheta",  QrangeName[i].c_str()), Form("%s;(e_{reco}-e_{true})/e_{true};Counts", Qrange[i].c_str()), 40, -0.01, 0.01);
        h_dp[i] = new TH1F(Form("H_%s_dphi",  QrangeName[i].c_str()), Form("%s;(e_{reco}-e_{true})/e_{true};Counts", Qrange[i].c_str()), 40, -0.01, 0.01);

        h_dhf[i] = new TH1F(Form("H_%s_dhf",  QrangeName[i].c_str()), Form("%s;(#delta_{h, reco}-#delta_{h, true})/#delta_{h, true};Counts", Qrange[i].c_str()), 100, -2, 2);
        h_dpt[i] = new TH1F(Form("H_%s_dpt",  QrangeName[i].c_str()), Form("%s;(p_{t, reco}-p_{t, true})/p_{t, true};Counts", Qrange[i].c_str()), 100, -2, 2);
    }

    std::string etype[2] = {"cal", "track"};
    for ( int i = 0; i < 2; i ++ )
    {
        h_tde[i] = new TH2F(Form("H_dEvT_%s", etype[i].c_str()), Form("dE_{%s} vs #theta_{mc};#theta_{e};dE/E", etype[i].c_str()), 180, 0, 180, 100,  -1, 1);
        h_tdp[i] = new TH2F(Form("H_dPvT_%s", etype[i].c_str()), Form("dP_{%s} vs #theta_{mc};#theta_{e};dP/P", etype[i].c_str()), 180, 0, 180, 100,  -1, 1);
        
        for ( int j = 0; j < 4; j ++ )
        {
            h_tde_tmp[i][j] = new TH2F(Form("H_dEvT_%s_%d", etype[i].c_str(),j), Form("dE_{%s} vs #theta_{mc};#theta_{e};dE/E", etype[i].c_str()), 180, 0, 180, 100,  -1, 1);
            h_tdp_tmp[i][j] = new TH2F(Form("H_dPvT_%s_%d", etype[i].c_str(),j), Form("dP_{%s} vs #theta_{mc};#theta_{e};dP/P", etype[i].c_str()), 180, 0, 180, 100,  -1, 1);
        }
    }
    h_tdphi = new TH2F(Form("H_dPhivT"), Form("d#phi vs #theta_{mc};#theta_{e};d#phi/#phi"), 180, 0, 180, 40, -0.01, 0.01);
    h_pde = new TH2F(Form("H_dEvP_%s", etype[1].c_str()), Form("dE_{%s} vs p_{mc};p_{t,e};dE/E", etype[1].c_str()), 180, 0, 60, 100,  -0.2, 0.2);
    h_ede = new TH2F(Form("H_dEvE_%s", etype[0].c_str()), Form("dE_{%s} vs E_{mc};E_{e};dE/E", etype[0].c_str()), 180, 0, 60, 100,  -0.2, 0.2);
    h_pt = new TH2F(Form("H_PvT_%s", etype[1].c_str()), Form("p_{t,%s} vs #theta_{mc};#theta_{e};p_{t}", etype[1].c_str()), 180, 0, 180, 240, 0, 80);
    h_et = new TH2F(Form("H_EvT_%s", etype[0].c_str()), Form("E_{%s} vs #theta_{mc};#theta_{e};E", etype[0].c_str()), 180, 0, 180, 240, 0, 80);

    h_hfs_dpz = new TH1F("H_HFS_dpz", Form(";(p_{z, reco}-p_{z, true})/p_{z, true};Counts"), 100, -1, 1);
    h_hfs_dpt = new TH1F("H_HFS_dpt", Form(";(p_{t, reco}-p_{t, true})/p_{t, true};Counts"), 100, -1, 1);
    h_hfs_de = new TH1F("H_HFS_de", Form(";(E_{reco}-E_{true})/E_{true};Counts"), 100, -1, 1);
    h_hfs_dpz_t = new TH2F("H_HFS_dpz_vT", Form(";#theta_{mc};(p_{z, reco}-p_{z, true})/p_{z, true}"), 180, 0, 180, 100, -1, 1);
    h_hfs_dpt_t = new TH2F("H_HFS_dpt_vT", Form(";#theta_{mc};(p_{t, reco}-p_{t, true})/p_{t, true}"), 180, 0, 180, 100, -1, 1);
    h_hfs_de_t = new TH2F("H_HFS_de_vT", Form(";#theta_{mc};(E_{reco}-E_{true})/E_{true}"), 180, 0, 180, 100, -1, 1);

}

void plot_energy_and_track(TCanvas* &c1, TCanvas* &c2, TCanvas* &c3, TCanvas* &c4)
{
    c1 = new TCanvas("c_cal_and_track_1", "Energy and Theta QA", 1600, 800);
    c1->Divide(4, 2);

    for ( int i = 0; i < 8; i ++ )
    {
        c1->cd(i+1);
        gPad->SetTickx(1);
        gPad->SetTicky(1);
        gPad->SetLeftMargin(0.2);
        gPad->SetBottomMargin(0.2);
    }
    
    TLegend* leg = new TLegend(0.60, 0.70, 0.88, 0.85);
    leg->SetTextSize(0.04);
	leg->SetBorderSize(0);
	leg->SetFillColor(0);
    leg->SetFillStyle(0);

    leg->AddEntry(h_dEt[0], "Track", "L");
    leg->AddEntry(h_dEc[0], "Calorimeter", "L");

    TLegend* leg_a = new TLegend(0.70, 0.5, 0.88, 0.85);
    leg_a->SetTextSize(0.04);
	leg_a->SetBorderSize(0);
	leg_a->SetFillColor(0);
    leg_a->SetFillStyle(0);

    leg_a->AddEntry(h_dp[0], "#phi_{e}", "L");
    leg_a->AddEntry(h_dt[0], "#theta_{e}", "L");
    

    for ( int i = 0; i < 4; i ++ )
    {
        c1->cd(i+1);
        h_dEt[i]->Draw("HIST SAME");
        h_dEt[i]->SetStats(0);
        h_dEt[i]->SetLineColor(kRed);
        format_eid_histogram(h_dEt[i]);

        h_dEc[i]->Draw("HIST SAME");
        h_dEc[i]->SetStats(0);
        h_dEc[i]->SetLineColor(kBlue);
        format_eid_histogram(h_dEc[i]);

        h_dEt[i]->GetYaxis()->SetRangeUser(0, fmax(h_dEt[i]->GetMaximum(), h_dEc[i]->GetMaximum())*1.1);
        h_dEt[i]->SetStats(0);
    
        leg->Draw();

        c1->cd(i+1+4);
        h_dt[i]->Draw("HIST SAME");
        h_dt[i]->SetStats(0);
        format_eid_histogram(h_dt[i]);

        h_dp[i]->Draw("HIST SAME");
        h_dp[i]->SetStats(0);
        h_dp[i]->SetLineColor(kRed);

        leg_a->Draw();
    }

    c2 = new TCanvas("c_cal_and_track_2", "Energy and Theta QA", 1200, 800);
    c2->Divide(3, 3);

    for ( int i = 0; i < 9; i ++ )
    {
        c2->cd(i+1);
        gPad->SetTickx(1);
        gPad->SetTicky(1);
        gPad->SetLogz();
        gPad->SetLeftMargin(0.2);
        gPad->SetBottomMargin(0.2);
    }

    c2->cd(1);
    h_tde[0]->Draw("COLZ");
    h_tde[0]->SetStats(0);
    format_eid_histogram(h_tde[0]);

    c2->cd(2);
    h_tde[1]->Draw("COLZ");
    h_tde[1]->SetStats(0);
    format_eid_histogram(h_tde[1]);

    c2->cd(3);
    h_tdphi->Draw("COLZ");
    h_tdphi->SetStats(0);
    format_eid_histogram(h_tdphi);

    c2->cd(4);
    h_et->Draw("COLZ");
    h_et->SetStats(0);
    format_eid_histogram(h_et);

    c2->cd(5);
    h_pt->Draw("COLZ");
    h_pt->SetStats(0);
    format_eid_histogram(h_pt);

    c2->cd(6);
    h_tdp[0]->Draw("COLZ");
    h_tdp[0]->SetStats(0);
    format_eid_histogram(h_tdp[0]);

    c2->cd(7);
    h_ede->Draw("COLZ");
    h_ede->SetStats(0);
    format_eid_histogram(h_ede);

    c2->cd(8);
    h_pde->Draw("COLZ");
    h_pde->SetStats(0);
    format_eid_histogram(h_pde);

    c2->cd(9);
    h_tdp[1]->Draw("COLZ");
    h_tdp[1]->SetStats(1);
    format_eid_histogram(h_tdp[1]);

    c3 = new TCanvas("c_cal_and_track_3", "Hadronic final state QA", 1600, 800);
    c3->Divide(4, 2);

    for ( int i = 0; i < 8; i ++ )
    {
        c3->cd(i+1);
        gPad->SetTickx(1);
        gPad->SetTicky(1);
        gPad->SetLeftMargin(0.2);
        gPad->SetBottomMargin(0.2);
    }

    for ( int i = 0; i < 4; i ++ )
    {
        c3->cd(i+1);
        h_dhf[i]->Draw("HIST SAME");
        h_dhf[i]->SetStats(0);
        // h_dhf[i]->SetLineColor(kRed);
        format_eid_histogram(h_dhf[i]);

        c3->cd(i+1+4);
        h_dpt[i]->Draw("HIST SAME");
        h_dpt[i]->SetStats(0);
        format_eid_histogram(h_dpt[i]);
    }

    c4 = new TCanvas("c_cal_and_track_4", "Hadronic final state QA 2", 1600, 800);
    c4->Divide(3, 2);

    for ( int i = 0; i < 6; i ++ )
    {
        c4->cd(i+1);
        gPad->SetTickx(1);
        gPad->SetTicky(1);
        gPad->SetLogz();
        gPad->SetLeftMargin(0.2);
        gPad->SetBottomMargin(0.2);
    }

    c4->cd(1);
    h_hfs_dpt->Draw("HIST");
    h_hfs_dpt->SetStats(0);
    format_eid_histogram(h_hfs_dpt);

    c4->cd(2);
    h_hfs_dpz->Draw("HIST");
    h_hfs_dpz->SetStats(0);
    format_eid_histogram(h_hfs_dpz);

    c4->cd(3);
    h_hfs_de->Draw("HIST");
    h_hfs_de->SetStats(0);
    format_eid_histogram(h_hfs_de);

    c4->cd(4);
    h_hfs_dpt_t->Draw("COLZ");
    h_hfs_dpt_t->SetStats(0);
    format_eid_histogram(h_hfs_dpt_t);

    c4->cd(5);
    h_hfs_dpt_t->Draw("COLZ");
    h_hfs_dpt_t->SetStats(0);
    format_eid_histogram(h_hfs_dpt_t);

    c4->cd(6);
    h_hfs_de_t->Draw("COLZ");  
    h_hfs_de_t->SetStats(0); 
    format_eid_histogram(h_hfs_de_t); 

    return;
}
