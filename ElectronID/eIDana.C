// Find inclusive scattered electrons

#include "eIDana.h"

// void eIDana(int Ee, int Eh, int select_region, int sr, int is_truth_eID, int file0, int analyse_p)
void eIDana(int Ee, int Eh, std::string ev_type, int is_truth_eID, int analyse_p)
{
    std::cout << "** Analysing inclusive electrons, energy is set to: " << Ee << "x" << Eh << std::endl;

    // Standard setup

    SetePICStyle();

    std::string eID_type;
    if  ( is_truth_eID == 1 )
        eID_type = "truth";
    else if ( is_truth_eID == 2 )
        eID_type = "mc";
    else
        eID_type = "recon"; 

    AnaManager* ana_manager = new AnaManager("eID" + eID_type + ev_type);
    // ana_manager->Initialize(select_region, sr, file0, analyse_p);
    ana_manager->InitializeForLocal(ev_type);

    // .. input setup
    auto reader = podio::ROOTReader();
    // reader.openFiles(ana_manager->GetInputNames());
    reader.openFiles(ana_manager->GetLocalInputNames());

    // .. output setup;
    CreateOutputTree(ana_manager->GetOutputName()); 

    // .. ElectronID setup
    ElectronID* eFinder = new ElectronID(Ee, Eh);
    LorentzRotation boost = analyse_p ? getBoost( Ee, Eh, MASS_ELECTRON, MASS_PROTON) : getBoost( Ee, Eh, MASS_ELECTRON, MASS_NEUTRON);
    eFinder->SetBoost(boost);

    DefineHistograms();

    // Analysis loop

    for( size_t ev = 0; ev < reader.getEntries("events"); ev++ )
    {
        auto raw = reader.readNextEntry("events");
        if(!raw) 
        {
            std::cerr << "readNextEntry returned null at event " << ev << "\n";
            break;
        }
        podio::Frame event(std::move(raw));
        eFinder->SetEvent(&event);

        if(ev%100==0) 
        cout << "Analysing event " << ev << "/" << reader.getEntries("events") << std::endl;

        // Generator information (mcID)
        edm4hep::MCParticleCollection e_mc = eFinder->GetMCElectron();
        if(e_mc.size() > 0) 
        {
            eID_status = FOUND_MC;
            // Calculate kinematic variables using MC electron
            TLorentzVector kprime;
            kprime.SetXYZM(e_mc[0].getMomentum().x, e_mc[0].getMomentum().y, e_mc[0].getMomentum().z, MASS_ELECTRON);
            CalculateElectronKinematics(Ee, Eh, kprime, mc_xB, mc_Q2, mc_W2, mc_y, mc_nu);
        }

        // Use MC to find reconstructed electron (TruthID)
        auto e_truth = eFinder->GetTruthReconElectron();
        if(e_truth.size() > 0) 
        {
            eID_status = FOUND_TRUTH;
            h_n_clusters_n_tracks->Fill( e_truth[0].getTracks().size(), e_truth[0].getClusters().size());
        }
            
        // Find scattered electrons (reconID)
        auto e_candidates = eFinder->FindScatteredElectron();
        edm4eic::ReconstructedParticle e_rec;

        double TrackEminusPzSum = 0;
        double CalEminusPzSum = 0;
        eFinder->GetEminusPzSum(TrackEminusPzSum, CalEminusPzSum);
        h_TrackEminusPz->Fill(TrackEminusPzSum);
        h_CalEminusPz->Fill(CalEminusPzSum);

        // If there are multiple candidates, select one with highest pT
        if(e_candidates.size() > 0) 
        {			
            e_rec = eFinder->SelectHighestPT(e_candidates);
            mc_PDG = eFinder->Check_eID(e_rec);

            if ( mc_PDG == 0 )
                eID_status = FOUND_E;
            else if ( mc_PDG == -211 )
                eID_status = FOUND_PI;
            else
                eID_status = FOUND_OTHERS;

            auto recoMC = eFinder->GetMC(e_rec);
            if ( recoMC.isAvailable() )
            {
                TLorentzVector recokprime;
                recokprime.SetXYZM(recoMC.getMomentum().x, recoMC.getMomentum().y, recoMC.getMomentum().z, MASS_ELECTRON);
                CalculateElectronKinematics(Ee, Eh, recokprime, rec_xB, rec_Q2, rec_W2, rec_y, rec_nu);
            }

            h_cand_mul->Fill(e_candidates.size());
            if ( mc_PDG == 0 )
                h_cand_mul_eHighPt->Fill(e_candidates.size());
            else
                h_cand_mul_oHighPt->Fill(e_candidates.size());
        }

        // Fill histograms
        for ( const auto& det_val : eFinder->e_det )
        {
            h_EoP_e->Fill(det_val.recon_EoP);
            h_isoE_e->Fill(det_val.recon_isoE);
        }
        for ( const auto& det_val : eFinder->pi_det )
        {
            h_EoP_pi->Fill(det_val.recon_EoP);
            h_isoE_pi->Fill(det_val.recon_isoE);
        }
        for ( const auto& det_val : eFinder->else_det )
        {
            h_EoP_else->Fill(det_val.recon_EoP);
            h_isoE_else->Fill(det_val.recon_isoE);
        }

        outTree->Fill();
        ResetVariables();
    }

    // Canvas
    double draw_max = 0.;

    TCanvas* c_EoP = new TCanvas("c_EoP", "c_EoP", 1000, 600);
    c_EoP->SetLogy();

    DrawParComparison(c_EoP, h_EoP_e, h_EoP_pi, h_EoP_else, draw_max);
    DrawVerticalLine(c_EoP, eFinder->get_mEoP_min(), draw_max);
    DrawVerticalLine(c_EoP, eFinder->get_mEoP_max(), draw_max);

    TCanvas* c_isoE = new TCanvas("c_isoE", "c_isoE", 1000, 600);
    c_isoE->SetLogy();

    DrawParComparison(c_isoE, h_isoE_e, h_isoE_pi, h_isoE_else, draw_max);
    DrawVerticalLine(c_isoE, eFinder->get_mIsoE(), draw_max);

    TCanvas* c_EminusPz = new TCanvas("c_EminusPz", "c_EminusPz", 1000, 600);

    DrawTCComparison(c_EminusPz, h_TrackEminusPz, h_CalEminusPz, draw_max);
    DrawVerticalLine(c_EminusPz, 2*Ee, draw_max);

    TCanvas* c_reco_mul = new TCanvas("c_reco_mul", "c_reco_mul", 1000, 600);

    h_cand_mul->Draw("HIST");
    h_cand_mul->SetLineColor(kGray+2);

    h_cand_mul_eHighPt->Draw("HIST SAME");
    h_cand_mul_eHighPt->SetLineColor(kBlue);

    h_cand_mul_oHighPt->Draw("HIST SAME");
    h_cand_mul_oHighPt->SetLineColor(kOrange+7);

    TLegend* leg_mul = new TLegend(0.4, 0.6, 0.6, 0.88);
    leg_mul->SetBorderSize(0);
    leg_mul->SetFillStyle(0);
    leg_mul->AddEntry(h_cand_mul, "All candidates", "L");
    leg_mul->AddEntry(h_cand_mul_eHighPt, "Scat. e has highest p_{T}", "L");
    leg_mul->AddEntry(h_cand_mul_oHighPt, "Others have highest p_{T}", "L");
    leg_mul->Draw();

    TCanvas* c_n_clusters_n_tracks = new TCanvas("c_n_clusters_n_tracks", "c_n_clusters_n_tracks", 1000, 600);
    h_n_clusters_n_tracks->Scale(1.0/h_n_clusters_n_tracks->GetEntries());
    h_n_clusters_n_tracks->Draw("COLZ TEXT");

    // Save

    outFile->cd();
    outTree->Write(outTree->GetName(), 2);

    c_EoP->Write(c_EoP->GetName(), 2);
    c_isoE->Write(c_isoE->GetName(), 2);
    c_n_clusters_n_tracks->Write(c_n_clusters_n_tracks->GetName(), 2);
    c_EminusPz->Write(c_EminusPz->GetName(), 2);
    c_reco_mul->Write(c_reco_mul->GetName(), 2);

    c_EoP->SaveAs(Form("%dx%d_%s_EoP.pdf", 18, 275, ev_type.c_str()));
    c_isoE->SaveAs(Form("%dx%d_%s_isoE.pdf", 18, 275, ev_type.c_str()));
    c_EminusPz->SaveAs(Form("%dx%d_%s_EminusPz.pdf", 18, 275, ev_type.c_str()));
    c_reco_mul->SaveAs(Form("%dx%d_%s_reco_mul.pdf", 18, 275, ev_type.c_str()));
    c_n_clusters_n_tracks->SaveAs(Form("%dx%d_%s_n_clusters_n_tracks.pdf", 18, 275, ev_type.c_str()));

    return;
}

void DefineHistograms() {

    h_EoP_e = new TH1D("h_EoP_e", "EoP e; E/p; Counts", 100, 0., 2.);
    h_EoP_pi = new TH1D("h_EoP_pi", "EoP pi; E/p; Counts", 100, 0., 2.);
    h_EoP_else = new TH1D("h_EoP_else", "EoP; E/p; Counts", 100, 0., 2.);

    h_isoE_e = new TH1D("h_isoE_e", "Isolation Energy; Iso. E; Counts", 100, 0., 2.);
    h_isoE_pi = new TH1D("h_isoE_pi", "Isolation Energy; Iso. E; Counts", 100, 0., 2.);
    h_isoE_else = new TH1D("h_isoE_else", "Isolation Energy; Iso. E; Counts", 100, 0., 2.);

    h_TrackEminusPz = new TH1D("h_TrackEminusPz", "#Sigma(E - Pz); #Sigma(E - Pz); Counts", 200, 0., 50.);
    h_CalEminusPz = new TH1D("h_CalEminusPz", "#Sigma(E - Pz); #Sigma(E - Pz); Counts", 200, 0., 50.);

    h_n_clusters_n_tracks = new TH2D("h_n_clusters_n_tracks", "Number of clusters vs number of tracks; N_{tracks}; N_{clusters}", 5, -0.5, 4.5, 5, -0.5, 4.5);

    h_cand_mul = new TH1D("h_cand_mul", "Scattered electron candidates multiplicity; N_{candidates}; Counts", 10, -0.5, 9.5);
    h_cand_mul_eHighPt = new TH1D("h_cand_mul_eHighPt", "Scattered electron candidates multiplicity (high p_{T,e}); N_{candidates}; Counts", 10, -0.5, 9.5);
    h_cand_mul_oHighPt = new TH1D("h_cand_mul_oHighPt", "Scattered electron candidates multiplicity (high p_{T,others}); N_{candidates}; Counts", 10, -0.5, 9.5);

    return;
}

void DrawVerticalLine(TCanvas* &c, double x_pos, double y_max) {

    c->cd();
    c->Modified();
    c->Update();

    TLine* line = new TLine(x_pos, 0, x_pos, y_max);
    line->SetLineColor(kBlack);
    line->SetLineStyle(7);
    line->Draw("SAME");

    return;
}

void DrawTCComparison(TCanvas* &c, TH1D* &ht, TH1D* &hc, double &draw_max) {

    c->cd();

    hc->SetLineColor(kGray);
    hc->SetFillColor(kGray);
    hc->SetFillStyle(3003);
    hc->Draw("HIST");
    draw_max = 1.2*std::max({hc->GetMaximum(), ht->GetMaximum()});
    hc->SetMaximum(draw_max);

    ht->SetLineColor(kBlue);
    ht->SetFillColor(kBlue);
    ht->SetFillStyle(3003);
    ht->Draw("HIST SAME");

    TLegend* leg = new TLegend(0.2, 0.6, 0.4, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(ht, "Using E_{Track}", "L");
    leg->AddEntry(hc, "Using E_{Cluster}", "L");
    leg->Draw();

    return;
}

void DrawParComparison(TCanvas* &c, TH1D* &h1, TH1D* &h2, TH1D* &h3, double &draw_max) {

    c->cd();

    h3->Draw("HIST");
    h3->SetLineColor(kGray+2);
    // h3->SetFillColor(kGray);
    // h3->SetFillStyle(3003);
    draw_max = 1.2*std::max({h1->GetMaximum(), h2->GetMaximum(), h3->GetMaximum()});
    h3->SetMaximum(draw_max);

    h2->Draw("HIST SAME");
    h2->SetLineColor(kBlue);
    // h2->SetFillColor(kBlue);
    // h2->SetFillStyle(3003);

    h1->Draw("HIST SAME");
    h1->SetLineWidth(2);
    h1->SetLineColor(kRed);
    h1->SetFillColor(kRed);
    h1->SetFillStyle(3003);

    TLegend* leg = new TLegend(0.6, 0.6, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(h1, "Electrons", "L");
    leg->AddEntry(h2, "Pions", "L");
    leg->AddEntry(h3, "Others", "L");
    leg->Draw();

    return;
}

void CreateOutputTree(TString outFileName) {

	outFile = new TFile(outFileName, "RECREATE");
	outTree = new TTree("T_eID", "T_eID");

    outTree->Branch("eID_status", &eID_status);
    outTree->Branch("mc_PDG", &mc_PDG);

	outTree->Branch("mc_xB", &mc_xB);
	outTree->Branch("mc_Q2", &mc_Q2);
	outTree->Branch("mc_W2", &mc_W2);
	outTree->Branch("mc_y",	 &mc_y);
	outTree->Branch("mc_nu", &mc_nu);
    
    return;
}

void ResetVariables() {

	eID_status = NO_MC;
    mc_PDG = 0;

	mc_xB = -999;
	mc_Q2 = -999;
	mc_W2 = -999;
	mc_y = -999;
	mc_nu = -999;

    rec_xB = -999;
	rec_Q2 = -999;
	rec_W2 = -999;
	rec_y = -999;
	rec_nu = -999;

    return;
}

void CalculateElectronKinematics(double fEe, double fEh, TLorentzVector kf, double& xB, double& Q2, double& W2, double& y, double& nu) {

		TLorentzVector ki; ki.SetXYZM(0., 0., -fEe, MASS_ELECTRON);
		TLorentzVector P = GetHadronBeam(fEh);
		TLorentzVector q = ki - kf;
		Q2 = -(q.Dot(q));
		nu = (q.Dot(P))/MASS_PROTON;
		xB = Q2/(2.*MASS_PROTON*nu);
		y  = (q.Dot(P))/(ki.Dot(P));
		W2  = MASS_PROTON*MASS_PROTON + (2.*MASS_PROTON*nu) - Q2;		
}

TLorentzVector GetHadronBeam(double fEh) {
 
	TLorentzVector hadron_beam;
	hadron_beam.SetX(fEh*sin(CROSSING_ANGLE));
	hadron_beam.SetY(0.);
	hadron_beam.SetZ(fEh*cos(CROSSING_ANGLE));
	hadron_beam.SetE(std::hypot(fEh, MASS_PROTON));
	return hadron_beam;

}