#include "kRecAna.h"

int n_group = 0;
const std::string group[4] = {"minQ2=1", "minQ2=10", "minQ2=100", "minQ2=1000"};

void kRecAna(int beam_type, int Ee, int Eh)
{
    // ECCE binning
    // set_ECCE_bin();

    // ePIC plotting style setup
    std::string type_title[5] = {"e^{3}He", "ep", "#gammap", "ep w. BeamBG", "ep"};
    std::string energy_title = beam_type ? Form("%dx%d GeV", Ee, Eh) : Form("%dx%d GeV/A", Ee, Eh);
    std::string campaign = "26.03.1";
    if ( beam_type == 4 )
        campaign = "25.10.4";
    draw_manager = new DrawManager(type_title[beam_type], energy_title, campaign);
    draw_manager->SetEPIC("Performance");

    double text_lumi = 1; // fb^-1
    if ( beam_type == 0 && Ee == 10 && Eh == 166 )
        text_lumi = 1.5; // fb^-1
    if ( beam_type == 4 && Ee == 10 && Eh == 130 )
        text_lumi = 1.0; // fb^-1
    if ( beam_type == 4 && Ee == 10 && Eh == 250 )
        text_lumi = 2.5; // fb^-1
    draw_manager->SetLumi(text_lumi);

    //

    std::string setting = beam_type == 0 ? Form("eHe3_%dx%d", (int)Ee, (int)Eh) : Form("ep_%dx%d", (int)Ee, (int)Eh);
    n_group = beam_type == 0 ? 3 : 4;

    define_histograms();

    vector<Kinematics> algorithm[n_group];

    for ( int i = 0; i < n_group; i ++ )
    {
        // load files and trees
        std::string date = "en_26_03_1";
        if ( beam_type == 4 )
            date = "ep_25_10_4";
        TFile* beamFile = new TFile(Form("../data/%s/root_files/%s_%s_beam_combined.root", date.c_str(), setting.c_str(), group[i].c_str()));
        TFile* eidFile  = new TFile(Form("../data/%s/root_files/%s_%s_eid_combined.root", date.c_str(), setting.c_str(), group[i].c_str()));

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

        // setup kinematics algorithms
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

            double scatE = boosted_vTRe.E();
            if ( algorithm[i][MC].Q2 < 10 && boosted_vCLe.E() > 0 )
                scatE = boosted_vCLe.E();
            
            algorithm[i][EL].save_variables(calc_elec_method(scatE, boosted_vTRe.theta(), pt_had, sigma_h, Ee, Eh), Mp);

            // if ( algorithm[i][EL].Q2 < 10 || algorithm[i][EL].Q2 >= 1000 )
            //     algorithm[i][EL].save_variables(calc_elec_method(boosted_vCLe.E(), boosted_vCLe.theta(), pt_had, sigma_h, Ee, Eh), Mp);
        
            if ( vTRf->size() != 0 ) 
            {
                algorithm[i][JB].save_variables(calc_jb_method(scatE, boosted_vTRe.theta(), pt_had, sigma_h, Ee, Eh), Mp);
                algorithm[i][DA].save_variables(calc_da_method(scatE, boosted_vTRe.theta(), pt_had, sigma_h, Ee, Eh), Mp);
                algorithm[i][SIG].save_variables(calc_sig_method(scatE, boosted_vTRe.theta(), pt_had, sigma_h, Ee, Eh), Mp);
                algorithm[i][ESIG].save_variables(calc_esig_method(scatE, boosted_vTRe.theta(), pt_had, sigma_h, Ee, Eh), Mp);
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


            // track resolution -- lab frame
            
            trackQA.h_dE[Qrange]->Fill((vTRe->E()-vMCe->E())/vMCe->E());
            trackQA.h_dP[Qrange]->Fill((vTRe->P()-vMCe->P())/vMCe->P());
            trackQA.h_dtheta[Qrange]->Fill((vTRe->Theta()-vMCe->Theta())/vMCe->Theta());
            trackQA.h_dphi[Qrange]->Fill((vTRe->Phi()-vMCe->Phi())/vMCe->Phi());

            trackQA.h_dEvT[Qrange]->Fill(180-vMCe->Theta()*(180./M_PI), (vTRe->E()-vMCe->E())/vMCe->E());
            trackQA.h_dPvT[Qrange]->Fill(180-vMCe->Theta()*(180./M_PI), (vTRe->P()-vMCe->P())/vMCe->P());

            // cluster resolution -- lab frame

            double cluster_P = sqrt(vCLe->E()*vCLe->E() - MASS_ELECTRON*MASS_ELECTRON);

            clusterQA.h_dE[Qrange]->Fill((vCLe->E()-vMCe->E())/vMCe->E());
            clusterQA.h_dP[Qrange]->Fill((cluster_P-vMCe->P())/vMCe->P());
            
            clusterQA.h_dEvT[Qrange]->Fill(180-vMCe->Theta()*(180./M_PI), (vCLe->E()-vMCe->E())/vMCe->E());
            clusterQA.h_dPvT[Qrange]->Fill(180-vMCe->Theta()*(180./M_PI), (cluster_P-vMCe->P())/vMCe->P());
        }

        double total_lumi = beam_type == 0 ? text_lumi*3 : text_lumi; // fb^-1
        double gen_lumi = 0;
        gen_lumi = get_lumi(beam_type, Ee, Eh, i, nEntries, 0);
        if ( beam_type == 2 )
            gen_lumi = get_lumi(1, Ee, Eh, i, nEntries, 0);

        for ( int j = 0; j < 6; j ++ )
        {
            algorithm[i][j].h_qvq->Scale( total_lumi / gen_lumi );
            algorithm[i][j].h_xvx->Scale( total_lumi / gen_lumi );
            algorithm[i][j].h_yvy->Scale( total_lumi / gen_lumi );

            algorithm[i][j].h_xq->Scale( total_lumi / gen_lumi );
            algorithm[i][j].h_xq_eff->Scale( total_lumi / gen_lumi );
        } 

        trackQA.h_dEvT[i]->Scale( total_lumi / gen_lumi );
        trackQA.h_dPvT[i]->Scale( total_lumi / gen_lumi );
        clusterQA.h_dEvT[i]->Scale( total_lumi / gen_lumi );
        clusterQA.h_dPvT[i]->Scale( total_lumi / gen_lumi );

        trackQA.h_dEvT_sum->Add(trackQA.h_dEvT[i]);
        trackQA.h_dPvT_sum->Add(trackQA.h_dPvT[i]);
        clusterQA.h_dEvT_sum->Add(clusterQA.h_dEvT[i]);
        clusterQA.h_dPvT_sum->Add(clusterQA.h_dPvT[i]);
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

    draw_manager->LableAndCollect(c_reco_v_true);
    // c_reco_v_true->SaveAs(Form("../data/en_25_10_2/reco_v_true_%dx%d.png", 10, 166));

    plot_energy_and_track();

    set_2d_scale(h_xq_gen[0]);
    TCanvas* c_eRecon_eff[n_methods];
    for ( int i = 0; i < n_methods; i ++ )
    {
        process_eff_hist(h_xq_eff[i], h_xq_gen[i]);
        c_eRecon_eff[i] = draw_2d_efficiency(h_xq_eff[i], Form("c_%s_eff", algorithm[0][i].method_name.c_str()), Form("%s Efficiency", algorithm[0][i].method_name.c_str()), 1400, 600, false, true);
        draw_manager->LableAndCollect(c_eRecon_eff[i]);
    }

    TFile*outFile = new TFile(Form("%s_ReconStats.root", setting.c_str()), "RECREATE");
    // TFile*outFile = new TFile(Form("%s_ReconStats_ECCE.root", setting.c_str()), "RECREATE");
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
        trackQA.h_dE[i]= new TH1F(Form("H_%s_dE_track",  QrangeName[i].c_str()), Form("%s;(E_{reco}-E_{true})/E_{true};Counts", Qrange[i].c_str()), 100, -1, 1);
        trackQA.h_dP[i] = new TH1F(Form("H_%s_dP_track",  QrangeName[i].c_str()), Form("%s;(P_{reco}-P_{true})/P_{true};Counts", Qrange[i].c_str()), 100, -1, 1);
        trackQA.h_dtheta[i] = new TH1F(Form("H_%s_dtheta_track",  QrangeName[i].c_str()), Form("%s;(#delta_{reco}-#delta_{true})/#delta_{true};Counts", Qrange[i].c_str()), 40, -0.01, 0.01);
        trackQA.h_dphi[i] = new TH1F(Form("H_%s_dphi_track",  QrangeName[i].c_str()), Form("%s;(#phi_{reco}-#phi_{true})/#phi_{true};Counts", Qrange[i].c_str()), 40, -0.01, 0.01);

        clusterQA.h_dE[i] = new TH1F(Form("H_%s_dE_cluster",  QrangeName[i].c_str()), Form("%s;(E_{reco}-E_{true})/E_{true};Counts", Qrange[i].c_str()), 100, -1, 1);
        clusterQA.h_dP[i] = new TH1F(Form("H_%s_dP_cluster",  QrangeName[i].c_str()), Form("%s;(P_{reco}-P_{true})/P_{true};Counts", Qrange[i].c_str()), 100, -1, 1);

        trackQA.h_dEvT[i] = new TH2F(Form("H_dEvT_track_%s", QrangeName[i].c_str()), Form("%s;#theta_{mc} (deg);(E_{reco}-E_{true})/E_{true}", Qrange[i].c_str()), 180, 0, 180, 100, -1, 1);
        trackQA.h_dPvT[i] = new TH2F(Form("H_dPvT_track_%s", QrangeName[i].c_str()), Form("%s;#theta_{mc} (deg);(P_{reco}-P_{true})/P_{true}", Qrange[i].c_str()), 180, 0, 180, 100, -1, 1);
        clusterQA.h_dEvT[i] = new TH2F(Form("H_dEvT_cluster_%s", QrangeName[i].c_str()), Form("%s;#theta_{mc} (deg);(E_{reco}-E_{true})/E_{true}", Qrange[i].c_str()), 180, 0, 180, 100, -1, 1);    
        clusterQA.h_dPvT[i] = new TH2F(Form("H_dPvT_cluster_%s", QrangeName[i].c_str()), Form("%s;#theta_{mc} (deg);(P_{reco}-P_{true})/P_{true}", Qrange[i].c_str()), 180, 0, 180, 100, -1, 1);

        h_dhf[i] = new TH1F(Form("H_%s_dhf",  QrangeName[i].c_str()), Form("%s;(#delta_{h, reco}-#delta_{h, true})/#delta_{h, true};Counts", Qrange[i].c_str()), 100, -2, 2);
        h_dpt[i] = new TH1F(Form("H_%s_dpt",  QrangeName[i].c_str()), Form("%s;(p_{t, reco}-p_{t, true})/p_{t, true};Counts", Qrange[i].c_str()), 100, -2, 2);
    }

    trackQA.h_dEvT_sum = new TH2F(Form("H_dEvT_track_sum"), Form(";#theta_{mc} (deg);(E_{track}-E_{true})/E_{true}"), 180, 0, 180, 100, -1, 1);
    trackQA.h_dPvT_sum = new TH2F(Form("H_dPvT_track_sum"), Form(";#theta_{mc} (deg);(P_{track}-P_{true})/P_{true}"), 180, 0, 180, 100, -1, 1);
    clusterQA.h_dEvT_sum = new TH2F(Form("H_dEvT_cluster_sum"), Form(";#theta_{mc} (deg);(E_{cluster}-E_{true})/E_{true}"), 180, 0, 180, 100, -1, 1);    
    clusterQA.h_dPvT_sum = new TH2F(Form("H_dPvT_cluster_sum"), Form(";#theta_{mc} (deg);(P_{cluster}-P_{true})/P_{true}"), 180, 0, 180, 100, -1, 1);

    h_hfs_dpz = new TH1F("H_HFS_dpz", Form(";(p_{z, reco}-p_{z, true})/p_{z, true};Counts"), 100, -1, 1);
    h_hfs_dpt = new TH1F("H_HFS_dpt", Form(";(p_{t, reco}-p_{t, true})/p_{t, true};Counts"), 100, -1, 1);
    h_hfs_de = new TH1F("H_HFS_de", Form(";(E_{reco}-E_{true})/E_{true};Counts"), 100, -1, 1);
    h_hfs_dpz_t = new TH2F("H_HFS_dpz_vT", Form(";#theta_{mc};(p_{z, reco}-p_{z, true})/p_{z, true}"), 180, 0, 180, 100, -1, 1);
    h_hfs_dpt_t = new TH2F("H_HFS_dpt_vT", Form(";#theta_{mc};(p_{t, reco}-p_{t, true})/p_{t, true}"), 180, 0, 180, 100, -1, 1);
    h_hfs_de_t = new TH2F("H_HFS_de_vT", Form(";#theta_{mc};(E_{reco}-E_{true})/E_{true}"), 180, 0, 180, 100, -1, 1);

}

void plot_energy_and_track()
{
    TCanvas* c1 = new TCanvas("c_EP_resolution", "Energy and Momentum QA", 1400, 1000);
    c1->Divide(2, 2);

    for ( int i = 0; i < 4; i ++ )
    {
        c1->cd(i+1);
        format_pad();
    }
        
    c1->cd(1);
    trackQA.h_dEvT_sum->Draw("COLZ");
    format_hist(trackQA.h_dEvT_sum);
    ReverseXAxis(trackQA.h_dEvT_sum);
    draw_eta_axis(trackQA.h_dEvT_sum);

    const double xmin = gPad->GetUxmin();
    const double xmax = gPad->GetUxmax();
    const double ymin = gPad->GetUymin();
    const double ymax = gPad->GetUymax();

    TLine* l_zero = new TLine(xmin, 0, xmax, 0);
	l_zero->SetLineColor(kBlack);
	l_zero->SetLineStyle(7);
    l_zero->SetLineWidth(2);
	l_zero->Draw();

    TLine* l_gap1 = new TLine(180-25, ymin, 180-25, ymax);
	l_gap1->SetLineColor(kBlack);
	l_gap1->SetLineStyle(7);
    l_gap1->SetLineWidth(2);
	// l_gap1->Draw();

    TLine* l_gap2 = new TLine(180-160, ymin, 180-160, ymax);
	l_gap2->SetLineColor(kBlack);
	l_gap2->SetLineStyle(7);
    l_gap2->SetLineWidth(2);
	// l_gap2->Draw();

    c1->cd(2);
    clusterQA.h_dEvT_sum->Draw("COLZ");
    format_hist(clusterQA.h_dEvT_sum);
    ReverseXAxis(clusterQA.h_dEvT_sum);
    draw_eta_axis(clusterQA.h_dEvT_sum);

    l_zero->Draw();
    // l_gap1->Draw();
    // l_gap2->Draw();

    c1->cd(3);
    trackQA.h_dPvT_sum->Draw("COLZ");
    format_hist(trackQA.h_dPvT_sum);
    ReverseXAxis(trackQA.h_dPvT_sum);
    draw_eta_axis(trackQA.h_dPvT_sum);

    l_zero->Draw();
    // l_gap1->Draw();
    // l_gap2->Draw();

    c1->cd(4);
    clusterQA.h_dPvT_sum->Draw("COLZ");
    format_hist(clusterQA.h_dPvT_sum);
    ReverseXAxis(clusterQA.h_dPvT_sum);
    draw_eta_axis(clusterQA.h_dPvT_sum);

    l_zero->Draw();
    // l_gap1->Draw();
    // l_gap2->Draw();

    draw_manager->LableAndCollect(c1,1);
    c1->SaveAs(Form("../data/eID/26.03_energy_and_track_resolution.png"));

    TCanvas* c2 = new TCanvas("c_Q2range_resolution", "Resolution comparison for Q2", 1400, 400);
    c2->Divide(n_group, 1);
    for ( int i = 0; i < n_group; i ++ )
    {
        c2->cd(i+1);
        clusterQA.h_dE[i]->Draw("HIST SAME");
        clusterQA.h_dE[i]->SetLineColor(kBlue);

        trackQA.h_dE[i]->Draw("HIST SAME");
        trackQA.h_dE[i]->SetLineColor(kRed);
    }
    draw_manager->LableAndCollect(c2);

    return;
}

void format_pad()
{
    gPad->SetLogz();
    gPad->SetTickx(1);
    gPad->SetTicky(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.25);

    gPad->Modified();
    gPad->Update();
}

void format_hist(TH1* hist)
{
    hist->GetXaxis()->SetTitleOffset(1.2);
    hist->GetXaxis()->CenterTitle();
    hist->GetXaxis()->SetLabelFont(42);
    hist->GetXaxis()->SetTitleFont(42);

    hist->GetYaxis()->SetTitleOffset(1.4);
    hist->GetYaxis()->CenterTitle();
    hist->GetYaxis()->SetTitleFont(42);
    hist->GetYaxis()->SetLabelFont(42);
    hist->GetYaxis()->SetNdivisions(505);
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
    const int nEta = 5;
    const double etaVals[nEta] = {-3.5, -0.88, 0, 0.88, 3.5}; 

    TLatex latex;
    latex.SetTextAlign(23);     // centered x, top y
    latex.SetTextSize(0.043);

    for (int i = 0; i < nEta; ++i)
    {
        const double eta = etaVals[i];

        // theta(eta) = 2 * atan(exp(-eta)) in radians -> degrees
        const double theta_deg = 2.0 * TMath::ATan(TMath::Exp(-eta)) * TMath::RadToDeg();

        if (theta_deg < xmin || theta_deg > xmax) continue;

        TLine* tick = new TLine(theta_deg, y_axis, theta_deg, y_axis + tick_h);
        tick->Draw();

        latex.DrawLatex(180-theta_deg, y_axis - 0.013 * yr, Form("%.2f", eta));
    }

    // Axis title
    TLatex title;
    title.SetTextAlign(23);
    title.SetTextSize(0.043);
    title.DrawLatex((xmin + xmax) * 0.5, y_axis - 0.066 * yr, "#eta");
}

void ReverseXAxis(TH1 *h)
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