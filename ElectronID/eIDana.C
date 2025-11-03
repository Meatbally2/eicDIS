// Find inclusive scattered electrons

#include "eIDana.h"

void eIDana(int Ee, int Eh, int select_region, int sr, int is_truth_eID, int all_file, int analyse_p)
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

    AnaManager* ana_manager = new AnaManager("eID" + eID_type);
    ana_manager->Initialize(select_region, sr, all_file, analyse_p);

    // .. input setup
    auto reader = podio::ROOTReader();
    reader.openFiles(ana_manager->GetInputNames());

    // .. output setup;
    CreateOutputTree(ana_manager->GetOutputName()); 

    // .. ElectronID setup
    ElectronID* eFinder = new ElectronID(Ee, Eh);

    DefineHistograms();

    // Analysis loop

    for( size_t ev = 0; ev < reader.getEntries("events"); ev++ )
    {
        const auto event = podio::Frame(reader.readNextEntry("events"));
        eFinder->SetEvent(&event);

        if(ev%100==0) 
        cout << "Analysing event " << ev << "/" << reader.getEntries("events") << std::endl;

        // Generator information (mcID)
        edm4hep::MCParticleCollection e_mc = eFinder->GetMCElectron();
        if(e_mc.size() == 0) 
            continue;

        // Use MC to find reconstructed electron (TruthID)
        auto e_truth = eFinder->GetTruthReconElectron();
        if(e_truth.size() > 0) 
        {
            eID_status = NO_FOUND;
            h_n_clusters_n_tracks->Fill( e_truth[0].getTracks().size(), e_truth[0].getClusters().size());
        }
            
        // Find scattered electrons (reconID)
        auto e_candidates = eFinder->FindScatteredElectron();
        edm4eic::ReconstructedParticle e_rec;
        
        // If there are multiple candidates, select one with highest pT
        if(e_candidates.size() > 0) 
        {			
            e_rec = eFinder->SelectHighestPT(e_candidates);
            mc_PBG = eFinder->Check_eID(e_rec);

            if ( mc_PBG == 0 )
                eID_status = FOUND_E;
            else if ( mc_PBG == -211 )
                eID_status = FOUND_PI;
            else
                eID_status = FOUND_OTHERS;
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

        // Calculate kinematic variables using MC electron
		TLorentzVector kprime;
		kprime.SetXYZM(e_mc[0].getMomentum().x, e_mc[0].getMomentum().y, e_mc[0].getMomentum().z, MASS_ELECTRON);
		CalculateElectronKinematics(Ee, Eh, kprime, mc_xB, mc_Q2, mc_W2, mc_y, mc_nu);

        outTree->Fill();
        ResetVariables();
    }

    // Canvas
    double draw_max = 0.;

    TCanvas* c_EoP = new TCanvas("c_EoP", "c_EoP", 1000, 600);
    c_EoP->SetLogy();

    DrawComparison(c_EoP, h_EoP_e, h_EoP_pi, h_EoP_else, draw_max);

    c_EoP->cd();
    c_EoP->Update();

    TLine* line_EoP_min = new TLine(eFinder->get_mEoP_min(), 0, eFinder->get_mEoP_min(), draw_max);
    line_EoP_min->SetLineColor(kBlack);
    line_EoP_min->SetLineStyle(7);
    line_EoP_min->Draw("SAME");
    TLine* line_EoP_max = new TLine(eFinder->get_mEoP_max(), 0, eFinder->get_mEoP_max(), draw_max);
    line_EoP_max->SetLineColor(kBlack);
    line_EoP_max->SetLineStyle(7);
    line_EoP_max->Draw("SAME");

    TCanvas* c_isoE = new TCanvas("c_isoE", "c_isoE", 1000, 600);
    c_isoE->SetLogy();

    DrawComparison(c_isoE, h_isoE_e, h_isoE_pi, h_isoE_else, draw_max);
    c_isoE->cd();
    c_isoE->Update();

    TLine* line_isoE_min = new TLine(eFinder->get_mIsoE(), 0, eFinder->get_mIsoE(), draw_max);
    line_isoE_min->SetLineColor(kBlack);
    line_isoE_min->SetLineStyle(7);
    line_isoE_min->Draw("SAME");

    TCanvas* c_n_clusters_n_tracks = new TCanvas("c_n_clusters_n_tracks", "c_n_clusters_n_tracks", 1000, 600);
    h_n_clusters_n_tracks->Scale(1.0/h_n_clusters_n_tracks->GetEntries());
    h_n_clusters_n_tracks->Draw("COLZ TEXT");

    // Save

    outFile->cd();
    outTree->Write(outTree->GetName(), 2);

    c_EoP->Write(c_EoP->GetName(), 2);
    c_isoE->Write(c_isoE->GetName(), 2);
    c_n_clusters_n_tracks->Write(c_n_clusters_n_tracks->GetName(), 2);

    return;
}

void DefineHistograms() {

    h_EoP_e = new TH1D("h_EoP_e", "EoP e; E/p; Counts", 100, 0., 2.);
    h_EoP_pi = new TH1D("h_EoP_pi", "EoP pi; E/p; Counts", 100, 0., 2.);
    h_EoP_else = new TH1D("h_EoP_else", "EoP; E/p; Counts", 100, 0., 2.);

    h_isoE_e = new TH1D("h_isoE_e", "Isolation Energy; Iso. E; Counts", 100, 0., 2.);
    h_isoE_pi = new TH1D("h_isoE_pi", "Isolation Energy; Iso. E; Counts", 100, 0., 2.);
    h_isoE_else = new TH1D("h_isoE_else", "Isolation Energy; Iso. E; Counts", 100, 0., 2.);

    h_n_clusters_n_tracks = new TH2D("h_n_clusters_n_tracks", "Number of clusters vs number of tracks; N_{tracks}; N_{clusters}", 5, -0.5, 4.5, 5, -0.5, 4.5);

    return;
}

void DrawComparison(TCanvas* c, TH1D* &h1, TH1D* &h2, TH1D* &h3, double &draw_max) {

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
    leg->AddEntry(h1, "Electrons", "f");
    leg->AddEntry(h2, "Pions", "f");
    leg->AddEntry(h3, "Others", "f");
    leg->Draw();

    return;
}

void CreateOutputTree(TString outFileName) {

	outFile = new TFile(outFileName, "RECREATE");
	outTree = new TTree("T_eID", "T_eID");

    outTree->Branch("eID_status", &eID_status);
    outTree->Branch("mc_PBG", &mc_PBG);

	outTree->Branch("mc_xB", &mc_xB);
	outTree->Branch("mc_Q2", &mc_Q2);
	outTree->Branch("mc_W2", &mc_W2);
	outTree->Branch("mc_y",	 &mc_y);
	outTree->Branch("mc_nu", &mc_nu);
    
    return;
}

void ResetVariables() {

	eID_status = NO_RECON;
    mc_PBG = 0;

	mc_xB = -999;
	mc_Q2 = -999;
	mc_W2 = -999;
	mc_y = -999;
	mc_nu = -999;

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