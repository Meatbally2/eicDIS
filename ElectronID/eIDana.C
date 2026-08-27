// Find inclusive scattered electrons

#include "../GlobalUtil/preLoadLib.hh" // load lib first otherwise the newest eicshell will not work with the code
#include "eIDana.h"

void eIDana(int Ee, int Eh, int beam_type, int select_region=0, int sr=0, int file0=-1)
{
    std::cout << "** Analysing inclusive electrons, energy is set to: " << Ee << "x" << Eh << std::endl;

    // Standard setup
    AnaManager* ana_manager = new AnaManager("eid");
    ana_manager->SetBeamEnergy(Ee, Eh);
    ana_manager->Initialize(select_region, sr, file0, beam_type);

    std::string type_title[6] = {"e^{3}He", "ep", "#gammap", "ep w. BeamBG", "ep", "ep"};
    std::string energy_title = beam_type ? Form("%dx%d GeV", Ee, Eh) : Form("%dx%d GeV/A", Ee, Eh);
    DrawManager* draw_manager = new DrawManager(type_title[beam_type], energy_title, ana_manager->campaign);
    draw_manager->SetEPIC();

    if (select_region) // no sr for pi_bg runs
    {
        if ( beam_type )
            draw_manager->SetQ2min(pow(10,sr));
        else
            draw_manager->SetQ2range(pow(10,sr), pow(10,sr+1));
    }

    // .. input setup
    podio::ROOTReader* reader = new podio::ROOTReader();
    reader->openFiles(ana_manager->GetInputNames());

    std::cout << "** Input files loaded. Setting up analysis... " << std::endl;

    // .. output setup;
    CreateOutputTree(ana_manager->GetOutputName()); 

    // .. ElectronID setup
    ElectronID* eFinder = new ElectronID(Ee, Eh);
    eFinder->SetMinTrackPoints(4);
    // techinically need to add a check for types of nucleon, but good enough for a quick check of x Q2. precise recon will be done in kin recon.
    LorentzRotation boost = beam_type ? getBoost( Ee, Eh, MASS_ELECTRON, MASS_PROTON) : getBoost( Ee, Eh, MASS_ELECTRON, MASS_NEUTRON); 
    eFinder->SetBoost(boost);

    DefineHistograms();

    // Analysis loop

    std::cout << "** Starting analysis loop... " << std::endl;

    int countMCe = 0;
    int countReconE = 0;

    for( size_t ev = 0; ev < reader->getEntries("events"); ev++ )
    {
        // std::cout << "** Processing event " << ev << " ... " << std::endl;
        // auto raw = reader->readNextEntry("events");
        auto raw = reader->readEntry("events", ev);  // Use readEntry with index instead
            
        if(!raw) 
        {
            std::cerr << "readNextEntry returned null at event " << ev << "\n";
            break;
        }
        // else
        //     std::cout << "** Event " << ev << " read. " << std::endl;
        podio::Frame event(std::move(raw));
        eFinder->SetEvent(&event);

        if(ev%100==0) 
            cout << "Analysing event " << ev << "/" << reader->getEntries("events") << endl;

        // Generator information (mcID)
        // eFinder->GetMCElectron();
        const auto& e_mc = eFinder->GetMCElectron();
        // if (e_mc.size() != 1)
        //     continue;
        
        h_n_scat_elec->Fill(e_mc.size());
        if(e_mc.size() > 0) 
        {
            eID_status = FOUND_MC;
            // Calculate kinematic variables using MC electron
            TLorentzVector kprime;
            kprime.SetXYZM(e_mc[0].getMomentum().x, e_mc[0].getMomentum().y, e_mc[0].getMomentum().z, MASS_ELECTRON);
            CalculateElectronKinematics(Ee, Eh, kprime, mc_xB, mc_Q2, mc_W2, mc_y, mc_nu);
            vMC_e.SetPxPyPzE(e_mc[0].getMomentum().x, e_mc[0].getMomentum().y, e_mc[0].getMomentum().z, e_mc[0].getEnergy());
            // vMC_e = boost(vMC_e);
            // if ( (vMC_e.Theta()*(180./M_PI) > 150) && (vMC_e.Theta()*(180./M_PI) < 155) )
            countMCe += e_mc.size();

            h_pt_theta->Fill(vMC_e.Theta()*(180./M_PI), abs(vMC_e.Pt()));
        }

        // Use MC to find reconstructed electron (TruthID)
        auto e_truth = eFinder->GetTruthReconElectron();
        if(e_truth.size() > 0) 
        {
            eID_status = FOUND_TRUTH;
            h_n_clusters_n_tracks->Fill( e_truth[0].getTracks().size(), e_truth[0].getClusters().size());

            h_n_cluster_in_cone->Fill(eFinder->rcpart_n_clusters);
            if ( e_truth[0].getClusters().size() > 0 )
                h_n_cluster_in_cone_found->Fill(eFinder->rcpart_n_clusters);

            if ( e_truth[0].getTracks().size() > 0 && e_truth[0].getClusters().size() > 0 )
                eRecon_status = FOUND_BOTH;
            else if ( e_truth[0].getTracks().size() > 0 )
                eRecon_status = FOUND_TRACK_ONLY;
            else if ( e_truth[0].getClusters().size() > 0 )
                eRecon_status = FOUND_CLUSTER_ONLY;
        }
            
        // Find scattered electrons (reconID)
        auto e_candidates = eFinder->FindScatteredElectron();
        const auto& summaries = eFinder->GetCandidateSummary();
        edm4eic::ReconstructedParticle e_rec;

        double TrackEminusPzSum = 0;
        double CalEminusPzSum = 0;
        eFinder->GetEminusPzSum(TrackEminusPzSum, CalEminusPzSum);
        EminusPz = TrackEminusPzSum;

        // If there are multiple candidates, select one with highest pT
        if(e_candidates.size() > 0) 
        {			
            h_TrackEminusPz->Fill(TrackEminusPzSum);
            h_CalEminusPz->Fill(CalEminusPzSum);

            e_rec = eFinder->SelectHighestPT(e_candidates);
            EoP = eFinder->GetCalorimeterEnergy(e_rec) / edm4hep::utils::magnitude(e_rec.getMomentum());

            vTRACK_e.SetPxPyPzE(e_rec.getMomentum().x, e_rec.getMomentum().y, e_rec.getMomentum().z, e_rec.getEnergy());
            vCLUSTER_e = eFinder->GetMomentumVectorFromCluster(e_rec, MASS_ELECTRON);
            // vCLUSTER_e.SetPxPyPzE(e_rec.getMomentum().x, e_rec.getMomentum().y, e_rec.getMomentum().z, eFinder->GetCalorimeterEnergy(e_rec));
            // vTRACK_e = boost(vTRACK_e);
            // vCLUSTER_e = boost(vCLUSTER_e);

            mc_PDG = eFinder->Check_eID(e_rec);

            h_cur_pur_eta->Fill(mc_PDG == 0, edm4hep::utils::eta(e_rec.getMomentum()));
            h_cur_pur_p->Fill(mc_PDG == 0, edm4hep::utils::magnitude(e_rec.getMomentum()));

            if ( mc_PDG == 0 )
                eID_status = FOUND_E;
            else if ( mc_PDG == -211 )
                eID_status = FOUND_PI;
            else
                eID_status = FOUND_OTHERS;

            // if ( (vMC_e.Theta()*(180./M_PI) > 150) && (vMC_e.Theta()*(180./M_PI) < 155) )
            if (eID_status == FOUND_E)
                countReconE++;

            auto recoMC = eFinder->GetMC(e_rec);
            if ( recoMC.isAvailable() )
            {
                TLorentzVector recokprime;
                recokprime.SetXYZM(recoMC.getMomentum().x, recoMC.getMomentum().y, recoMC.getMomentum().z, MASS_ELECTRON);
                CalculateElectronKinematics(Ee, Eh, recokprime, rec_xB, rec_Q2, rec_W2, rec_y, rec_nu);
                vMC_rec.SetPxPyPzE(recoMC.getMomentum().x, recoMC.getMomentum().y, recoMC.getMomentum().z, recoMC.getEnergy());
            }

            h_cand_mul->Fill(e_candidates.size());
            if ( mc_PDG == 0 )
                h_cand_mul_eHighPt->Fill(e_candidates.size());
            else
                h_cand_mul_oHighPt->Fill(e_candidates.size());

            // Reconstruct HFS
            auto mc_hfsCollection = eFinder->GetMCHadronicFinalState();
            for (const auto p : mc_hfsCollection) {
                PxPyPzEVector hf(p.getMomentum().x, p.getMomentum().y, p.getMomentum().z, p.getEnergy());
                // hf = boost(hf);
                vMC_hfs.push_back(hf);
            }

            auto hfsCollection = eFinder->FindHadronicFinalState(e_rec.getObjectID().index);
            for (const auto p : hfsCollection) {;
                PxPyPzEVector hf(p.getMomentum().x, p.getMomentum().y, p.getMomentum().z, p.getEnergy());
                // hf = boost(hf);
                vREC_hfs.push_back(hf);
            }

            if ( e_rec.getPDG() != 0 )
            {
                int reco_pid = abs(e_rec.getPDG());
                double reco_eta = edm4hep::utils::eta(recoMC.getMomentum());
                double reco_p = edm4hep::utils::magnitude(recoMC.getMomentum());
                FillEidPurity(reco_pid, mc_PDG, reco_eta, reco_p);
            }
        }

        // Additional purity set: all particles with negative tracks (no cluster/E/p requirement).
        for (const auto& s : summaries)
        {
            if (s.n_tracks <= 0 || s.charge >= 0)
                continue;

            FillNegTrackPurity(abs(s.reco_pdg), s.mc_pdg, s.det.eta, s.det.p);
        }

        // Fill histograms
        for ( const auto& s : summaries )
        {
            const auto& d = s.mod;
            if ( d.parType == 0 )
            {
                if ( s.track_theta > 158 && s.track_theta < 162 )
                {
                    h_EoP_gapF_mod->Fill(d.recon_EoP);
                    h_EoEH_gapF_mod->Fill(d.recon_EoEH);
                }
                    
                if ( s.cluster_theta > 22 && s.cluster_theta < 33 )
                {
                    h_EoP_gapB_mod->Fill(d.recon_EoP);
                    h_EoEH_gapB_mod->Fill(d.recon_EoEH);
                }
            }
        }

        for ( const auto& s : summaries )
        {
            if ( s.n_tracks <= 0 || s.charge >= 0 )
                continue;
                
            const auto& d = s.det;
            if ( d.parType == 0 )
            {
                h_nTPts_e->Fill(d.nTrackPoints);
                h_EoP_e->Fill(d.recon_EoP);
                h_isoE_e->Fill(d.recon_isoE);
                h_EoEH_e->Fill(d.recon_EoEH);
                h_PIDe_e->Fill(d.recon_Le);
                h_PIDh_e->Fill(d.recon_Lh);

                if ( vMC_e.Theta()*(180./M_PI) > 158 && vMC_e.Theta()*(180./M_PI) < 162 ) 
                {
                    h_EoP_gapF->Fill(d.recon_EoP);
                    h_EoEH_gapF->Fill(d.recon_EoEH);
                }
                    
                if ( vMC_e.Theta()*(180./M_PI) > 22 && vMC_e.Theta()*(180./M_PI) < 33 )
                {
                    h_EoP_gapB->Fill(d.recon_EoP);
                    h_EoEH_gapB->Fill(d.recon_EoEH);
                }
            }
            else if ( d.parType == -211 )
            {
                h_nTPts_pi->Fill(d.nTrackPoints);
                h_EoP_pi->Fill(d.recon_EoP);
                h_isoE_pi->Fill(d.recon_isoE);
                h_EoEH_pi->Fill(d.recon_EoEH);
                h_PIDe_pi->Fill(d.recon_Le);
                h_PIDh_pi->Fill(d.recon_Lh);
            }
            else if ( abs(d.parType) == 11 )
            {
                h_nTPts_jet_e->Fill(d.nTrackPoints);
                h_EoP_jet_e->Fill(d.recon_EoP);
                h_isoE_jet_e->Fill(d.recon_isoE);
                h_EoEH_jet_e->Fill(d.recon_EoEH);
                h_PIDe_jet_e->Fill(d.recon_Le);
                h_PIDh_jet_e->Fill(d.recon_Lh);
            }
            else
            {
                h_nTPts_else->Fill(d.nTrackPoints);
                h_EoP_else->Fill(d.recon_EoP);
                h_isoE_else->Fill(d.recon_isoE);
                h_EoEH_else->Fill(d.recon_EoEH);
                h_PIDe_else->Fill(d.recon_Le);
                h_PIDh_else->Fill(d.recon_Lh);
            }

            // How well is PID matching true particle overall if PID is reporting a particle
            FillPidPurity(d.recon_pID, d.parType, d.eta, d.p);

            // how often is PID reporting true particle in general
            FillPidSuccess(d.recon_pID, d.parType);
        }
  
        
        outTree->Fill();
        ResetVariables();
    }

    cout << "** Analysis finished. " << std::endl;
    cout << "DIS eID rate: " << (double)countReconE/countMCe * 100 << "%" << std::endl;

    // Canvas
    double draw_max = 0.;

    // TCanvas* c_nScatElec = new TCanvas("c_nScatElec", "c_nScatElec", 1000, 600);
    // h_n_scat_elec->SetLineColor(kBlue);
    // h_n_scat_elec->Draw("HIST");
    // draw_manager->LableAndCollect(c_nScatElec,2);

    TCanvas* c_nTPts = new TCanvas("c_nTPts", "c_nTPts", 1000, 600);
    DrawParComparison(c_nTPts, h_nTPts_e, h_nTPts_jet_e, h_nTPts_pi, h_nTPts_else, draw_max);
    DrawVerticalLine(c_nTPts, eFinder->GetMinTrackPoints()-0.5, draw_max);
    draw_manager->LableAndCollect(c_nTPts);

    TCanvas* c_EoP = new TCanvas("c_EoP", "c_EoP", 1000, 600);
    c_EoP->SetLogy();

    DrawParComparison(c_EoP, h_EoP_e, h_EoP_jet_e, h_EoP_pi, h_EoP_else, draw_max);
    DrawVerticalLine(c_EoP, eFinder->get_mEoP_min(), draw_max);
    DrawVerticalLine(c_EoP, eFinder->get_mEoP_max(), draw_max);
    draw_manager->LableAndCollect(c_EoP);

    TCanvas* c_gap = new TCanvas("c_gap", "c_gap", 1400, 800);
    c_gap->Divide(2,2);

    c_gap->cd(1);
    gPad->SetLogy();
    h_EoP_gapB->Draw("HIST");
    h_EoP_gapB->SetLineColor(kGray+2);
    h_EoP_gapB_mod->Draw("HIST SAME");
    h_EoP_gapB_mod->SetLineColor(kRed);

    c_gap->cd(2);
    gPad->SetLogy();
    h_EoP_gapF->Draw("HIST");
    h_EoP_gapF->SetLineColor(kGray+2);
    h_EoP_gapF_mod->Draw("HIST SAME");
    h_EoP_gapF_mod->SetLineColor(kRed);

    c_gap->cd(3);
    gPad->SetLogy();
    h_EoEH_gapB->Draw("HIST");
    h_EoEH_gapB->SetLineColor(kGray+2);
    h_EoEH_gapB_mod->Draw("HIST SAME");
    h_EoEH_gapB_mod->SetLineColor(kRed);

    c_gap->cd(4);
    gPad->SetLogy();
    h_EoEH_gapF->Draw("HIST");
    h_EoEH_gapF->SetLineColor(kGray+2);
    h_EoEH_gapF_mod->Draw("HIST SAME");
    h_EoEH_gapF_mod->SetLineColor(kRed);

    draw_manager->LableAndCollect(c_gap);

    TCanvas* c_isoE = new TCanvas("c_isoE", "c_isoE", 1000, 600);
    c_isoE->SetLogy();

    DrawParComparison(c_isoE, h_isoE_e, h_isoE_jet_e, h_isoE_pi, h_isoE_else, draw_max);
    DrawVerticalLine(c_isoE, eFinder->get_mIsoE(), draw_max);
    draw_manager->LableAndCollect(c_isoE);

    TCanvas* c_EoEH = new TCanvas("c_EoEH", "c_EoEH", 1000, 600);
    c_EoEH->SetLogy();

    DrawParComparison(c_EoEH, h_EoEH_e, h_EoEH_jet_e, h_EoEH_pi, h_EoEH_else, draw_max);
    DrawVerticalLine(c_EoEH, eFinder->get_mEoEH_min(), draw_max);
    draw_manager->LableAndCollect(c_EoEH);

    TCanvas* c_PIDe = new TCanvas("c_PIDe", "c_PIDe", 1000, 600);
    c_PIDe->SetLogy();

    DrawParComparison(c_PIDe, h_PIDe_e, h_PIDe_jet_e, h_PIDe_pi, h_PIDe_else, draw_max);
    draw_manager->LableAndCollect(c_PIDe);

    TCanvas* c_PIDh = new TCanvas("c_PIDh", "c_PIDh", 1000, 600);
    c_PIDh->SetLogy();
    
    DrawParComparison(c_PIDh, h_PIDh_e, h_PIDh_jet_e, h_PIDh_pi, h_PIDh_else, draw_max);
    DrawVerticalLine(c_PIDh, eFinder->get_mPID_veto(), draw_max);
    draw_manager->LableAndCollect(c_PIDh);

    TCanvas* c_EminusPz = new TCanvas("c_EminusPz", "c_EminusPz", 1000, 600);

    DrawTCComparison(c_EminusPz, h_TrackEminusPz, h_CalEminusPz, draw_max);
    DrawVerticalLine(c_EminusPz, 2*Ee, draw_max);
    draw_manager->LableAndCollect(c_EminusPz);

    TCanvas* c_reco_mul = new TCanvas("c_reco_mul", "c_reco_mul", 1000, 600);
    c_reco_mul->SetLogy();

    h_cand_mul->Draw("HIST");
    h_cand_mul->SetLineColor(kGray+2);
    h_cand_mul->GetYaxis()->SetRangeUser(1, h_cand_mul->GetMaximum()*1.5);

    h_cand_mul_eHighPt->Draw("HIST SAME");
    h_cand_mul_eHighPt->SetLineColor(kBlue);

    h_cand_mul_oHighPt->Draw("HIST SAME");
    h_cand_mul_oHighPt->SetLineColor(kOrange+7);

    TLegend* leg_mul = new TLegend(0.6, 0.6, 0.8, 0.88);
    leg_mul->SetBorderSize(0);
    leg_mul->SetFillStyle(0);
    leg_mul->AddEntry(h_cand_mul, "All candidates", "L");
    leg_mul->AddEntry(h_cand_mul_eHighPt, "Scat. e has highest p_{T}", "L");
    leg_mul->AddEntry(h_cand_mul_oHighPt, "Others have highest p_{T}", "L");
    leg_mul->Draw();

    draw_manager->LableAndCollect(c_reco_mul,2);

    TCanvas* c_n_clusters_n_tracks = new TCanvas("c_n_clusters_n_tracks", "c_n_clusters_n_tracks", 1000, 600);
    h_n_clusters_n_tracks->Scale(1.0/h_n_clusters_n_tracks->GetEntries());
    h_n_clusters_n_tracks->Draw("COLZ TEXT");
    draw_manager->LableAndCollect(c_n_clusters_n_tracks,2);

    // TCanvas* c_n_cluster_in_cone = new TCanvas("c_n_cluster_in_cone", "c_n_cluster_in_cone", 1000, 600);
    // h_n_cluster_in_cone->Draw("HIST");
    // h_n_cluster_in_cone->SetLineColor(kGray+2);
    // h_n_cluster_in_cone_found->Draw("HIST SAME");
    // h_n_cluster_in_cone_found->SetLineColor(kRed);
    // draw_manager->LableAndCollect(c_n_cluster_in_cone,2);

    TCanvas* c_pID_eff = new TCanvas("c_pID_eff", "c_pID_eff", 1000, 600);

    TH1D* h_eff_bar = (TH1D*)h_pID_eff->GetCopyPassedHisto();
    TH1D* h_eff_tot = (TH1D*)h_pID_eff->GetCopyTotalHisto();
    h_eff_bar->Divide(h_eff_bar, h_eff_tot, 1.0, 1.0, "B");  // binomial errors
    h_eff_bar->SetFillColor(kP8Green);
    h_eff_bar->SetFillStyle(3003);
    h_eff_bar->SetLineColor(kP8Green);
    h_eff_bar->SetMarkerColor(kP8Green);
    h_eff_bar->Draw("HIST SAME");

    h_eff_bar->GetYaxis()->SetRangeUser(0.0, 1.5);
    gPad->Update();

    h_pID_pur->Draw("SAME");
    h_pID_suc->Draw("SAME");
    h_pID_eff->Draw("SAME");

    TLegend* leg_pID = new TLegend(0.6, 0.72, 0.8, 0.92);
    leg_pID->SetBorderSize(0);
    leg_pID->SetFillStyle(0);
    leg_pID->AddEntry(h_pID_suc, "Purity in general", "LP");
    leg_pID->AddEntry(h_pID_pur, "Purity if pID exists", "LP");
    leg_pID->AddEntry(h_pID_eff, "Purity for e candidates", "LP");
    leg_pID->Draw();

    draw_manager->LableAndCollect(c_pID_eff);

    TCanvas* c_pID_pur_eta;
    DrawPurityCanvas(c_pID_pur_eta, "c_pID_pur_eta", h_cur_pur_eta, h_eVeto_pID_pur_eta, h_e_pID_pur_eta, h_pi_pID_pur_eta, h_K_pID_pur_eta, h_p_pID_pur_eta, draw_manager);

    TCanvas* c_pID_pur_p;
    DrawPurityCanvas(c_pID_pur_p, "c_pID_pur_p", h_cur_pur_p, h_eVeto_pID_pur_p, h_e_pID_pur_p, h_pi_pID_pur_p, h_K_pID_pur_p, h_p_pID_pur_p, draw_manager);

    TCanvas* c_eID_pur_eta;
    DrawPurityCanvas(c_eID_pur_eta, "c_eID_pur_eta", h_cur_pur_eta, h_eVeto_eID_pur_eta, h_e_eID_pur_eta, h_pi_eID_pur_eta, h_K_eID_pur_eta, h_p_eID_pur_eta, draw_manager);

    TCanvas* c_eID_pur_p;
    DrawPurityCanvas(c_eID_pur_p, "c_eID_pur_p", h_cur_pur_p, h_eVeto_eID_pur_p, h_e_eID_pur_p, h_pi_eID_pur_p, h_K_eID_pur_p, h_p_eID_pur_p, draw_manager);

    TCanvas* c_trk_pur_eta;
    DrawPurityCanvas(c_trk_pur_eta, "c_trk_pur_eta", h_cur_pur_eta, h_eVeto_trk_pur_eta, h_e_trk_pur_eta, h_pi_trk_pur_eta, h_K_trk_pur_eta, h_p_trk_pur_eta, draw_manager);

    TCanvas* c_trk_pur_p;
    DrawPurityCanvas(c_trk_pur_p, "c_trk_pur_p", h_cur_pur_p, h_eVeto_trk_pur_p, h_e_trk_pur_p, h_pi_trk_pur_p, h_K_trk_pur_p, h_p_trk_pur_p, draw_manager);

    TCanvas* c_pt_theta = new TCanvas("c_pt_theta", "c_pt_theta", 1000, 600);
    h_pt_theta->Draw("COLZ");
    draw_manager->LableAndCollect(c_pt_theta);

    // TCanvas* c_pid_score = eFinder->MakePlots();
    // draw_manager->LableAndCollect(c_pid_score);

    // Save

    outFile->cd();
    outTree->Write(outTree->GetName(), 2);

    if ( file0 == 0 )
        draw_manager->SaveToTree(outFile);

    return;
}

void DefineHistograms() {

    h_nTPts_e = new TH1D("h_nTPts_e", "Number of Track Points for e; N_{Track Points}; Counts", 14, -0.5, 13.5);
    h_nTPts_jet_e = new TH1D("h_nTPts_jet_e", "Number of Track Points for other e's; N_{Track Points}; Counts", 14, -0.5, 13.5);
    h_nTPts_pi = new TH1D("h_nTPts_pi", "Number of Track Points for #pi; N_{Track Points}; Counts", 14, -0.5, 13.5);
    h_nTPts_else = new TH1D("h_nTPts_else", "Number of Track Points for others; N_{Track Points}; Counts", 14, -0.5, 13.5);

    h_EoP_e = new TH1D("h_EoP_e", "EoP e; E/p; Counts", 100, 0., 2.);
    h_EoP_jet_e = new TH1D("h_EoP_jet_e", "EoP other e's; E/p; Counts", 100, 0., 2.);
    h_EoP_pi = new TH1D("h_EoP_pi", "EoP pi; E/p; Counts", 100, 0., 2.);
    h_EoP_else = new TH1D("h_EoP_else", "EoP; E/p; Counts", 100, 0., 2.);

    h_EoP_gapF = new TH1D("h_EoP_gapF", "EoP gapF; E/p; Counts", 100, 0., 2.);
    h_EoP_gapB = new TH1D("h_EoP_gapB", "EoP gapB; E/p; Counts", 100, 0., 2.);
    h_EoP_gapF_mod = new TH1D("h_EoP_gapF_mod", "EoP gapF; E/p; Counts", 100, 0., 2.);
    h_EoP_gapB_mod = new TH1D("h_EoP_gapB_mod", "EoP gapB; E/p; Counts", 100, 0., 2.);

    h_isoE_e = new TH1D("h_isoE_e", "Isolation Energy; Iso. E; Counts", 110, 0., 1.1);
    h_isoE_jet_e = new TH1D("h_isoE_jet_e", "Isolation Energy other e's; Iso. E; Counts", 110, 0., 1.1);
    h_isoE_pi = new TH1D("h_isoE_pi", "Isolation Energy; Iso. E; Counts", 110, 0., 1.1);
    h_isoE_else = new TH1D("h_isoE_else", "Isolation Energy; Iso. E; Counts", 110, 0., 1.1);

    h_EoEH_e = new TH1D("h_EoEH_e", "E/E+H e; E/E+H; Counts", 110, 0., 1.1);
    h_EoEH_jet_e = new TH1D("h_EoEH_jet_e", "E/E+H other e's; E/E+H; Counts", 110, 0., 1.1);
    h_EoEH_pi = new TH1D("h_EoEH_pi", "E/E+H pi; E/E+H; Counts", 110, 0., 1.1);
    h_EoEH_else = new TH1D("h_EoEH_else", "E/E+H; E/E+H; Counts", 110, 0., 1.1);

    h_EoEH_gapB = new TH1D("h_EoEH_gapB", "E/E+H gapB; E/E+H; Counts", 110, 0., 1.1);
    h_EoEH_gapF = new TH1D("h_EoEH_gapF", "E/E+H gapF; E/E+H; Counts", 110, 0., 1.1);
    h_EoEH_gapF_mod = new TH1D("h_EoEH_gapF_mod", "E/E+H gapF; E/E+H; Counts", 110, 0., 1.1);
    h_EoEH_gapB_mod = new TH1D("h_EoEH_gapB_mod", "E/E+H gapB; E/E+H; Counts", 110, 0., 1.1);

    h_PIDe_e = new TH1D("h_PIDe_e", "PID e; PID; Counts", 100, 0., 1.);
    h_PIDe_jet_e = new TH1D("h_PIDe_jet_e", "PID other e's; PID; Counts", 100, 0., 1.);
    h_PIDe_pi = new TH1D("h_PIDe_pi", "PID pi; PID; Counts", 100, 0., 1.);
    h_PIDe_else = new TH1D("h_PIDe_else", "PID; PID; Counts", 100, 0., 1.);

    h_PIDh_e = new TH1D("h_PIDh_e", "PID e; PID L_{h}/(L_{h}+L_{e}); Counts", 100, 0., 1.);
    h_PIDh_jet_e = new TH1D("h_PIDh_jet_e", "PID other e's; PID L_{h}/(L_{h}+L_{e}); Counts", 100, 0., 1.);
    h_PIDh_pi = new TH1D("h_PIDh_pi", "PID pi; PID L_{h}/(L_{h}+L_{e}); Counts", 100, 0., 1.);
    h_PIDh_else = new TH1D("h_PIDh_else", "PID; PID L_{h}/(L_{h}+L_{e}); Counts", 100, 0., 1.);

    h_pt_theta = new TH2D("h_pt_theta", "p_{T} vs #theta; #theta; p_{T}", 180, 0., 180, 50, 0., 50.);

    h_TrackEminusPz = new TH1D("h_TrackEminusPz", "#Sigma(E - Pz); #Sigma(E - Pz); Counts", 200, 0., 50.);
    h_CalEminusPz = new TH1D("h_CalEminusPz", "#Sigma(E - Pz); #Sigma(E - Pz); Counts", 200, 0., 50.);

    h_n_scat_elec = new TH1D("h_n_scat_elec", "Number of scattered electrons; N_{e}; Counts", 10, -0.5, 9.5);
    h_n_clusters_n_tracks = new TH2D("h_n_clusters_n_tracks", "Number of clusters vs number of tracks; N_{tracks}; N_{clusters}", 5, -0.5, 4.5, 5, -0.5, 4.5);

    h_cand_mul = new TH1D("h_cand_mul", "Scattered electron candidates multiplicity; N_{candidates}; Counts", 10, -0.5, 9.5);
    h_cand_mul_eHighPt = new TH1D("h_cand_mul_eHighPt", "Scattered electron candidates multiplicity (high p_{T,e}); N_{candidates}; Counts", 10, -0.5, 9.5);
    h_cand_mul_oHighPt = new TH1D("h_cand_mul_oHighPt", "Scattered electron candidates multiplicity (high p_{T,others}); N_{candidates}; Counts", 10, -0.5, 9.5);

    h_n_cluster_in_cone = new TH1D("h_n_cluster_in_cone", "Number of clusters in isolation cone; N_{clusters in cone}; Counts", 20, -0.5, 19.5);
    h_n_cluster_in_cone_found = new TH1D("h_n_cluster_in_cone_found", "Number of clusters in isolation cone for found electrons; N_{clusters in cone}; Counts", 20, -0.5, 19.5);

    h_pID_pur = new TEfficiency("h_pID_pur", ";PDG;Purity", 5, -0.5, 4.5);
    TH1* h_total = const_cast<TH1*>(h_pID_pur->GetTotalHistogram());
    h_total->GetXaxis()->SetBinLabel(1, "Not e");
    h_total->GetXaxis()->SetBinLabel(2, "e");
    h_total->GetXaxis()->SetBinLabel(3, "#pi");
    h_total->GetXaxis()->SetBinLabel(4, "K");
    h_total->GetXaxis()->SetBinLabel(5, "p");
    h_total->LabelsOption("h", "X");
    h_pID_pur->SetMarkerColor(kP8Blue);
    h_pID_pur->SetLineColor(kP8Blue);
    h_pID_pur->SetLineWidth(2);
    h_pID_pur->SetMarkerStyle(21);

    h_pID_suc = new TEfficiency("h_pID_suc", ";PDG;Success Rate", 5, -0.5, 4.5);
    h_pID_suc->SetMarkerColor(kP8Red);
    h_pID_suc->SetLineColor(kP8Red);
    h_pID_suc->SetLineWidth(2);
    h_pID_suc->SetMarkerStyle(20);

    h_pID_eff = new TEfficiency("h_pID_eff", ";PDG;Efficiency", 5, -0.5, 4.5);
    h_pID_eff->SetMarkerColor(kP8Green);
    h_pID_eff->SetLineColor(kP8Green);
    h_pID_eff->SetMarkerStyle(29);

    //

    int color[5] = {kP10Blue, kP10Brown, kP10Green, kP10Ash, kP10Red};
    int marker[5] = {20, 21, 22, 23, 29};

    TEfficiency** eID_eta[5] = {&h_e_eID_pur_eta, &h_pi_eID_pur_eta, &h_K_eID_pur_eta, &h_p_eID_pur_eta, &h_eVeto_eID_pur_eta};
    TEfficiency** pID_eta[5] = {&h_e_pID_pur_eta, &h_pi_pID_pur_eta, &h_K_pID_pur_eta, &h_p_pID_pur_eta, &h_eVeto_pID_pur_eta};
    TEfficiency** eID_p[5] = {&h_e_eID_pur_p, &h_pi_eID_pur_p, &h_K_eID_pur_p, &h_p_eID_pur_p, &h_eVeto_eID_pur_p};
    TEfficiency** pID_p[5] = {&h_e_pID_pur_p, &h_pi_pID_pur_p, &h_K_pID_pur_p, &h_p_pID_pur_p, &h_eVeto_pID_pur_p};
    TEfficiency** trk_eta[5] = {&h_e_trk_pur_eta, &h_pi_trk_pur_eta, &h_K_trk_pur_eta, &h_p_trk_pur_eta, &h_eVeto_trk_pur_eta};
    TEfficiency** trk_p[5] = {&h_e_trk_pur_p, &h_pi_trk_pur_p, &h_K_trk_pur_p, &h_p_trk_pur_p, &h_eVeto_trk_pur_p};

    const char* eID_eta_names[5] = {"h_e_eID_pur_eta", "h_pi_eID_pur_eta", "h_K_eID_pur_eta", "h_p_eID_pur_eta", "h_eVeto_eID_pur_eta"};
    const char* pID_eta_names[5] = {"h_e_pID_pur_eta", "h_pi_pID_pur_eta", "h_K_pID_pur_eta", "h_p_pID_pur_eta", "h_eVeto_pID_pur_eta"};
    const char* eID_p_names[5] = {"h_e_eID_pur_p", "h_pi_eID_pur_p", "h_K_eID_pur_p", "h_p_eID_pur_p", "h_eVeto_eID_pur_p"};
    const char* pID_p_names[5] = {"h_e_pID_pur_p", "h_pi_pID_pur_p", "h_K_pID_pur_p", "h_p_pID_pur_p", "h_eVeto_pID_pur_p"};
    const char* trk_eta_names[5] = {"h_e_trk_pur_eta", "h_pi_trk_pur_eta", "h_K_trk_pur_eta", "h_p_trk_pur_eta", "h_eVeto_trk_pur_eta"};
    const char* trk_p_names[5] = {"h_e_trk_pur_p", "h_pi_trk_pur_p", "h_K_trk_pur_p", "h_p_trk_pur_p", "h_eVeto_trk_pur_p"};

    for (int i = 0; i < 5; ++i) {
        SetupPurityEff(*eID_eta[i], eID_eta_names[i], "#eta", 20, -5., 5., color[i], marker[i]);
        SetupPurityEff(*pID_eta[i], pID_eta_names[i], "#eta", 20, -5., 5., color[i], marker[i]);
        SetupPurityEff(*eID_p[i], eID_p_names[i], "p [GeV/c]", 150, 0., 150., color[i], marker[i]);
        SetupPurityEff(*pID_p[i], pID_p_names[i], "p [GeV/c]", 150, 0., 150., color[i], marker[i]);
        SetupPurityEff(*trk_eta[i], trk_eta_names[i], "#eta", 20, -5., 5., color[i], marker[i]);
        SetupPurityEff(*trk_p[i], trk_p_names[i], "p [GeV/c]", 150, 0., 150., color[i], marker[i]);
    }

    SetupPurityEff(h_cur_pur_eta, "h_cur_pur_eta", "#eta", 20, -5., 5., kP10Yellow, 20);
    h_cur_pur_eta->SetLineStyle(1);
    h_cur_pur_eta->SetLineWidth(2);

    SetupPurityEff(h_cur_pur_p, "h_cur_pur_p", "p [GeV/c]", 150, 0., 150., kP10Yellow, 20);
    h_cur_pur_p->SetLineStyle(1);
    h_cur_pur_p->SetLineWidth(2);

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

    TLegend* leg = new TLegend(0.6, 0.6, 0.8, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(ht, "Using E_{Track}", "L");
    leg->AddEntry(hc, "Using E_{Cluster}", "L");
    leg->Draw();

    return;
}

void DrawParComparison(TCanvas* &c, TH1D* &h1, TH1D* &h2, TH1D* &h3, TH1D* &h4, double &draw_max) {

    c->cd();

    h4->Draw("HIST");
    h4->SetLineColor(kGray+2);
    // h4->SetFillColor(kGray);
    // h4->SetFillStyle(3003);
    draw_max = 1.2*std::max({h1->GetMaximum(), h2->GetMaximum(), h3->GetMaximum(), h4->GetMaximum()});
    h4->SetMaximum(draw_max);

    h3->Draw("HIST SAME");
    h3->SetLineColor(kBlue);
    // h3->SetFillColor(kGreen+3);
    // h3->SetFillStyle(3003);  

    h2->Draw("HIST SAME");
    h2->SetLineColor(kViolet);
    // h2->SetFillColor(kBlue);
    // h2->SetFillStyle(3003);

    h1->Draw("HIST SAME");
    h1->SetLineWidth(2);
    h1->SetLineColor(kRed);
    h1->SetFillColor(kRed);
    h1->SetFillStyle(3003);

    double lx = h1->GetName() == TString("h_EoP_e") ? 0.7 : 0.42;
    TLegend* leg = new TLegend(lx, 0.6, lx+0.25, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(h1, "Electrons", "L");
    leg->AddEntry(h2, "Other e's", "L");
    leg->AddEntry(h3, "Pions", "L");
    leg->AddEntry(h4, "Others", "L");
    leg->Draw();

    return;
}

void CreateOutputTree(TString outFileName) {

	outFile = new TFile(outFileName, "RECREATE");
	outTree = new TTree("T_eID", "T_eID");

    outTree->Branch("eID_status", &eID_status);
    outTree->Branch("eRecon_status", &eRecon_status);
    outTree->Branch("mc_PDG", &mc_PDG);
    outTree->Branch("EminusPz", &EminusPz);
    outTree->Branch("EoP", &EoP);

	outTree->Branch("mc_xB", &mc_xB);
	outTree->Branch("mc_Q2", &mc_Q2);
	outTree->Branch("mc_W2", &mc_W2);
	outTree->Branch("mc_y",	 &mc_y);
	outTree->Branch("mc_nu", &mc_nu);

    outTree->Branch("rec_xB", &rec_xB);
	outTree->Branch("rec_Q2", &rec_Q2);
	outTree->Branch("rec_W2", &rec_W2);
	outTree->Branch("rec_y",  &rec_y);
	outTree->Branch("rec_nu", &rec_nu);

    outTree->Branch("vMC_e", &vMC_e);
	outTree->Branch("vTRACK_e", &vTRACK_e);
	outTree->Branch("vCLUSTER_e", &vCLUSTER_e);
    outTree->Branch("vMC_rec", &vMC_rec);
    outTree->Branch("vMC_hfs", &vMC_hfs);
    outTree->Branch("vREC_hfs", &vREC_hfs);

    return;
}

void ResetVariables() {

	eID_status = NO_MC;
    eRecon_status = NO_REC;
    mc_PDG = -999;
    EminusPz = -999;
    EoP = -999;

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

    vMC_e.SetPxPyPzE(0, 0, 0, 0);
	vTRACK_e.SetPxPyPzE(0, 0, 0, 0);
	vCLUSTER_e.SetPxPyPzE(0, 0, 0, 0);
    vMC_rec.SetPxPyPzE(0, 0, 0, 0);

    vMC_hfs.clear();
    vREC_hfs.clear();   

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

void SetupPurityEff(TEfficiency* &h, const char* name, const char* xAxisTitle, int nBins, double xMin, double xMax, int color, int marker)
{
    h = new TEfficiency(name, Form(";%s;Purity", xAxisTitle), nBins, xMin, xMax);
    h->SetMarkerColor(color);
    h->SetLineColor(color);
    h->SetMarkerStyle(marker);
}

void DrawPurityCanvas(TCanvas* &c, const char* canvasName,
    TEfficiency* h_base, TEfficiency* h_veto, TEfficiency* h_e, TEfficiency* h_pi, TEfficiency* h_K, TEfficiency* h_p,
    DrawManager* draw_manager)
{
    c = new TCanvas(canvasName, canvasName, 1000, 600);
    TH1D* h_base_line = nullptr;

    if (h_base) {
        h_base_line = (TH1D*)h_base->GetCopyPassedHisto();
        TH1D* h_base_total = (TH1D*)h_base->GetCopyTotalHisto();
        h_base_line->SetName(Form("%s_base_line", canvasName));
        h_base_line->SetDirectory(nullptr);
        h_base_line->Divide(h_base_line, h_base_total, 1.0, 1.0, "B");
        h_base_line->SetLineColor(h_base->GetLineColor());
        h_base_line->SetLineStyle(1);
        h_base_line->SetLineWidth(2);
        h_base_line->SetFillStyle(3003);
        h_base_line->SetFillColor(h_base->GetLineColor());
        h_base_line->SetMarkerSize(0.0);
        h_base_line->GetYaxis()->SetRangeUser(0.0, 1.5);
        h_base_line->Draw("HIST");
    }

    h_veto->Draw(h_base_line ? "SAME" : "");
    gPad->Update();
    h_veto->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.0, 1.5);

    h_e->Draw("SAME");
    h_pi->Draw("SAME");
    h_K->Draw("SAME");
    h_p->Draw("SAME");

    TLegend* leg = new TLegend(0.6, 0.7, 0.8, 0.92);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    if (h_base_line)
        leg->AddEntry(h_base_line, "baseline", "L");
    leg->AddEntry(h_veto, "electron veto", "LP");
    leg->AddEntry(h_e, "electron", "LP");
    leg->AddEntry(h_pi, "pion", "LP");
    leg->AddEntry(h_K, "kaon", "LP");
    leg->AddEntry(h_p, "proton", "LP");
    leg->Draw();

    draw_manager->LableAndCollect(c);
}

void FillEidPurity(int reco_pid, int mc_pdg, double eta, double momentum)
{
    const bool is_mc_electron = (mc_pdg == 0 || std::abs(mc_pdg) == 11);
    const bool is_mc_pion = (std::abs(mc_pdg) == 211);
    const bool is_mc_kaon = (std::abs(mc_pdg) == 321);
    const bool is_mc_proton = (std::abs(mc_pdg) == 2212);

    if (reco_pid != 11) {
        h_pID_eff->Fill(!is_mc_electron, 0);
        h_eVeto_eID_pur_eta->Fill(!is_mc_electron, eta);
        h_eVeto_eID_pur_p->Fill(!is_mc_electron, momentum);
    }
    if (reco_pid == 11) {
        h_pID_eff->Fill(is_mc_electron, 1);
        h_e_eID_pur_eta->Fill(is_mc_electron, eta);
        h_e_eID_pur_p->Fill(is_mc_electron, momentum);
    }
    if (reco_pid == 211) {
        h_pID_eff->Fill(is_mc_pion, 2);
        h_pi_eID_pur_eta->Fill(is_mc_pion, eta);
        h_pi_eID_pur_p->Fill(is_mc_pion, momentum);
    }
    if (reco_pid == 321) {
        h_pID_eff->Fill(is_mc_kaon, 3);
        h_K_eID_pur_eta->Fill(is_mc_kaon, eta);
        h_K_eID_pur_p->Fill(is_mc_kaon, momentum);
    }
    if (reco_pid == 2212) {
        h_pID_eff->Fill(is_mc_proton, 4);
        h_p_eID_pur_eta->Fill(is_mc_proton, eta);
        h_p_eID_pur_p->Fill(is_mc_proton, momentum);
    }
}

void FillNegTrackPurity(int reco_pid, int mc_pdg, double eta, double momentum)
{
    const bool is_mc_electron = (mc_pdg == 0 || std::abs(mc_pdg) == 11);
    const bool is_mc_pion = (std::abs(mc_pdg) == 211);
    const bool is_mc_kaon = (std::abs(mc_pdg) == 321);
    const bool is_mc_proton = (std::abs(mc_pdg) == 2212);

    if (reco_pid != 11) {
        h_eVeto_trk_pur_eta->Fill(!is_mc_electron, eta);
        h_eVeto_trk_pur_p->Fill(!is_mc_electron, momentum);
    }
    if (reco_pid == 11) {
        h_e_trk_pur_eta->Fill(is_mc_electron, eta);
        h_e_trk_pur_p->Fill(is_mc_electron, momentum);
    }
    if (reco_pid == 211) {
        h_pi_trk_pur_eta->Fill(is_mc_pion, eta);
        h_pi_trk_pur_p->Fill(is_mc_pion, momentum);
    }
    if (reco_pid == 321) {
        h_K_trk_pur_eta->Fill(is_mc_kaon, eta);
        h_K_trk_pur_p->Fill(is_mc_kaon, momentum);
    }
    if (reco_pid == 2212) {
        h_p_trk_pur_eta->Fill(is_mc_proton, eta);
        h_p_trk_pur_p->Fill(is_mc_proton, momentum);
    }
}

void FillPidPurity(int reco_pid, int par_type, double eta, double momentum)
{
    if (reco_pid == 0) {
        return;
    }

    const bool is_mc_electron = (par_type == 0 || std::abs(par_type) == 11);
    const bool is_mc_pion = (std::abs(par_type) == 211);
    const bool is_mc_kaon = (std::abs(par_type) == 321);
    const bool is_mc_proton = (std::abs(par_type) == 2212);

    if (reco_pid != 11) {
        h_pID_pur->Fill(!is_mc_electron, 0);
        h_eVeto_pID_pur_eta->Fill(!is_mc_electron, eta);
        h_eVeto_pID_pur_p->Fill(!is_mc_electron, momentum);
    }
    if (reco_pid == 11) {
        h_pID_pur->Fill(is_mc_electron, 1);
        h_e_pID_pur_eta->Fill(is_mc_electron, eta);
        h_e_pID_pur_p->Fill(is_mc_electron, momentum);
    }
    if (reco_pid == 211) {
        h_pID_pur->Fill(is_mc_pion, 2);
        h_pi_pID_pur_eta->Fill(is_mc_pion, eta);
        h_pi_pID_pur_p->Fill(is_mc_pion, momentum);
    }
    if (reco_pid == 321) {
        h_pID_pur->Fill(is_mc_kaon, 3);
        h_K_pID_pur_eta->Fill(is_mc_kaon, eta);
        h_K_pID_pur_p->Fill(is_mc_kaon, momentum);
    }
    if (reco_pid == 2212) {
        h_pID_pur->Fill(is_mc_proton, 4);
        h_p_pID_pur_eta->Fill(is_mc_proton, eta);
        h_p_pID_pur_p->Fill(is_mc_proton, momentum);
    }
}

void FillPidSuccess(int reco_pid, int par_type)
{
    const bool is_mc_electron = (par_type == 0 || std::abs(par_type) == 11);
    const bool is_mc_pion = (std::abs(par_type) == 211);
    const bool is_mc_kaon = (std::abs(par_type) == 321);
    const bool is_mc_proton = (std::abs(par_type) == 2212);

    if (!is_mc_electron) {
        h_pID_suc->Fill(reco_pid != 11, 0);
    }
    if (is_mc_electron) {
        h_pID_suc->Fill(reco_pid == 11, 1);
    }
    if (is_mc_pion) {
        h_pID_suc->Fill(reco_pid == 211, 2);
    }
    if (is_mc_kaon) {
        h_pID_suc->Fill(reco_pid == 321, 3);
    }
    if (is_mc_proton) {
        h_pID_suc->Fill(reco_pid == 2212, 4);
    }
}