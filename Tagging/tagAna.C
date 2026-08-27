
#include "../GlobalUtil/preLoadLib.hh"
#include "tagAna.h"
#include <cmath>

void tagAna(int Ee, int Eh, int beam_type, int select_region=0, int sr=0, int file0=-1)
{
    std::cout << "** Analysing far-forwards, energy is set to: " << Ee << "x" << Eh << std::endl;

    constexpr double nucleon_mass = 0.9382720813; // GeV
    ion_momentum_per_nucleon =std::sqrt(Eh * Eh - nucleon_mass * nucleon_mass);

    // Standard setup

    AnaManager* ana_manager = new AnaManager("tag");
    ana_manager->Initialize(select_region, sr, file0, beam_type);
    ana_manager->SetBeamEnergy(Ee, Eh);

    std::string type_title[6] = {"e^{3}He", "ep", "#gammap", "ep w. BeamBG", "ep", "ep DVMP"};
    std::string energy_title = beam_type ? Form("%dx%d GeV", Ee, Eh) : Form("%dx%d GeV/A", Ee, Eh);
    draw_manager = new DrawManager(type_title[beam_type], energy_title, ana_manager->campaign);
    draw_manager->SetEPIC();

    if (select_region)
    {
        if ( beam_type )
            draw_manager->SetQ2min(pow(10,sr));
        else
            draw_manager->SetQ2range(pow(10,sr), pow(10,sr+1));
    }
    
    // .. input setup
    auto reader = podio::ROOTReader();
    // reader.openFiles(ana_manager->GetLocalInputNames());
    reader.openFiles(ana_manager->GetInputNames());

    // .. output setup;
    CreateOutputTree(ana_manager->GetOutputName()); 

    // .. FarForward setup
    setup_omd();
    setup_rp();
    setup_zdc();

    h_xq2_pt = new TH1F(Form("h_xq2_pt"), ";p_{T} (GeV/c);Counts", 50, 0, 0.2);
    h_xq2_pt_tag = new TH1F(Form("h_xq2_pt_tag"), ";p_{T} (GeV/c);Counts", 50, 0, 0.2);
    h_xq2_pt_theta = new TH2F(Form("h_xq2_pt_theta"), ";#theta (mrad);p_{T} (GeV/c)", 50, 0, 1, 80, 0, 0.16);
    h_tag_mul[0] = new TH1F(Form("h_tag_mul_p"), ";Tag multiplicity;Counts", 10, 0, 10);
    h_tag_mul[1] = new TH1F(Form("h_tag_mul_n"), ";Tag multiplicity;Counts", 10, 0, 10);
    h_spectator_status = new TH1F("h_spectator_status", ";Spectator status;Events", 7, -0.5, 6.5);
    h_spectator_status->GetXaxis()->SetBinLabel(1, "Unclassified");
    h_spectator_status->GetXaxis()->SetBinLabel(2, "en ancestry");
    h_spectator_status->GetXaxis()->SetBinLabel(3, "ep initial");
    h_spectator_status->GetXaxis()->SetBinLabel(4, "Incomplete");
    h_spectator_status->GetXaxis()->SetBinLabel(5, "Mapping failed");
    h_spectator_status->GetXaxis()->SetBinLabel(6, "en kin clean");
    h_spectator_status->GetXaxis()->SetBinLabel(7, "en kin ambiguous");

    // Analysis loop

    for( size_t ev = 0; ev < reader.getEntries("events"); ev++ )
    // for( size_t ev = 0; ev < 10; ev++ )
    {
        auto raw = reader.readNextEntry("events");
        if(!raw) 
        {
            std::cerr << "readNextEntry returned null at event " << ev << "\n";
            continue;
        }
        podio::Frame event(std::move(raw));
        
        if(ev%100==0) 
            cout << "Analysing event " << ev << "/" << reader.getEntries("events") << std::endl;

        find_spectators(&event);
        h_spectator_status->Fill(spectator_status);
        process_ff(&event);
        outTree->Fill();
            
        // ResetVariables();
    }

    TCanvas* c_pt = new TCanvas("c_pt", "c_pt", 1400, 400);
    c_pt->Divide(3,1);

    c_pt->cd(1);
    h_xq2_pt->Draw("HIST SAME");
    h_xq2_pt->SetLineColor(kBlue);

    h_xq2_pt_tag->Draw("HIST SAME");
    h_xq2_pt_tag->SetLineColor(kRed);

    TLegend* leg = new TLegend(0.6, 0.7, 0.88, 0.88);
    leg->SetTextSize(0.05);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->AddEntry(h_xq2_pt, "MC", "L");
    leg->AddEntry(h_xq2_pt_tag, "Tagged", "L");
    leg->Draw();

    c_pt->cd(2);
    h_xq2_pt_theta->Draw("COLZ");

    c_pt->cd(3);
    const double tag_mul_max = std::max(h_tag_mul[0]->GetMaximum(), h_tag_mul[1]->GetMaximum());
    h_tag_mul[0]->SetMinimum(0.0);
    h_tag_mul[0]->SetMaximum(tag_mul_max > 0.0 ? 1.5 * tag_mul_max : 1.0);
    h_tag_mul[0]->Draw("HIST");
    h_tag_mul[0]->SetLineColor(kRed);

    h_tag_mul[1]->Draw("HIST SAME");
    h_tag_mul[1]->SetLineColor(kBlue);

    TLegend* leg_mul = new TLegend(0.6, 0.7, 0.88, 0.88);
    leg_mul->SetTextSize(0.05);
    leg_mul->SetBorderSize(0);
    leg_mul->SetFillColor(0);
    leg_mul->SetFillStyle(0);
    leg_mul->AddEntry(h_tag_mul[0], "ep DIS", "L");
    leg_mul->AddEntry(h_tag_mul[1], "en DIS", "L");
    leg_mul->Draw();

    // draw_manager->LableAndCollect(c_pt, 2);
    plot_ff();

    draw_manager->SaveToTree(outFile);

    outFile->cd();
    outTree->Write(outTree->GetName(), 2);
    h_spectator_status->Write();

    std::cout << "Spectator status counts: unclassified=" << h_spectator_status->GetBinContent(1)
              << ", en ancestry=" << h_spectator_status->GetBinContent(2)
              << ", ep initial=" << h_spectator_status->GetBinContent(3)
              << ", incomplete=" << h_spectator_status->GetBinContent(4)
              << ", mapping failed=" << h_spectator_status->GetBinContent(5)
              << ", en kin clean=" << h_spectator_status->GetBinContent(6)
              << ", en kin ambiguous=" << h_spectator_status->GetBinContent(7) << '\n';

    // outFile->Close();

    return;
}

void setup_rp()
{
    rpFinder = new FarForward("rp", "ForwardRomanPot");
    rpFinder->define_histograms();

    double rotate = 1./cos(-0.04545);
    double xShift[4] = {-1131.19, -1131.19, -1208.43, -1208.43};
    double zRange[4][2] = {{32535, 32555}, {32555, 32575}, {34235, 34255},{34255, 34275}};

    rpFinder->set_rotation(rotate);
    for ( int i = 0; i < 4; i ++ )
    {
        rpFinder->set_xShift(i, xShift[i]);
        rpFinder->set_zRange(i, zRange[i][0], zRange[i][1]);
    }

    ffFinder.push_back(rpFinder);

    return;
}

void setup_omd()
{
    omdFinder = new FarForward("omd", "ForwardOffMTracker");    
    omdFinder->define_histograms();

    double rotate = 1./cos(-0.047);
    // double xShift[4] = {-780.0, -780.0, -870.0, -870.0};
    // double zRange[4][2] = {{22490, 22512}, {22512, 22528}, {24490, 24510},{24510, 24530}};
    double xShift[4] = {-941.0+20, -941.0+20, -1032.0+20, -1032.0+20};
    double zRange[4][2] = {{25500, 25512}, {25515, 25528}, {26950, 27010}, {27015, 27030}};

    omdFinder->set_rotation(rotate);
    for ( int i = 0; i < 4; i ++ )
    {
        omdFinder->set_xShift(i, xShift[i]);
        omdFinder->set_zRange(i, zRange[i][0], zRange[i][1]);
    }

    ffFinder.push_back(omdFinder);

    return;
}

void setup_zdc()
{
    // zdcFinder = new FarForward("zdc", "HcalFarForwardZDCClusters"); 
    zdcFinder = new FarForward("zdc", "ReconstructedHcalFarForwardZDCNeutrals");    
    zdcFinder->define_histograms();
    return;
}

void process_ff(const podio::Frame* event)
{
    is_tagged = false;
    n_proton_tracks = 0;
    zdc_energy = 0.0;
    fZDCn = false;

    for (int s = 0; s < 2; s ++ )
    {
        if (spec.size() > s)
        {
            for ( int j = 0; j < 4; j ++ )
            {
                spec[s]->det_hits[0][j] = 0;
                spec[s]->det_hits[1][j] = 0;
            }
        }
    }

    for (int i = 0; i < ffFinder.size(); i ++ )
    {
        std::vector<int> index;
        for ( int s = 0; s < spec.size(); s ++ )
            index.push_back(spec[s]->mc_index);
    
        ffFinder[i]->SetEvent(event);
        ffFinder[i]->GetHits(index);

        for ( int s = 0; s < spec.size(); s ++ )
        {
            for ( int j = 0; j < 4; j ++ )
                spec[s]->det_hits[i][j] = ffFinder[i]->GetMCHits(s,j);
        }       

        ffFinder[i]->fill_hit_histograms();
        ffFinder[i]->ClearHits();
    }

    zdcFinder->SetEvent(event);
    zdc_energy = zdcFinder->GetEnergy();
    zdcFinder->fill_energy_histograms();

    // tagging
    n_proton_tracks = rpFinder->get_n_tracks() + omdFinder->get_n_tracks();

    // Reconstructed event-level multiplicity is independent of whether the
    // truth spectators could be identified.
    if (struck_type == 0 || struck_type == 1)
        h_tag_mul[struck_type]->Fill(n_proton_tracks);

    if (n_proton_tracks >= 2)
        is_tagged = true;

    if (zdcFinder->is_ZDC_neutron && zdc_energy > 40.0)
    {
        fZDCn = true;

        if ( is_tagged )
            is_tagged = false;
    }

    // mc tagging
    for ( int s = 0; s < spec.size(); s ++ )
    {
        int n_omd_mc = *std::min_element(spec[s]->det_hits[0], spec[s]->det_hits[0] + 4); // omd first
        int n_rp_mc = *std::min_element(spec[s]->det_hits[1], spec[s]->det_hits[1] + 4);
        
        // cout << "Struck type " << struck_type << " Spec " << s << " MC hits in RP: " << n_rp_mc << ", OMD: " << n_omd_mc << std::endl;

        if ( n_rp_mc + n_omd_mc > 0)
        {
            spec[s]->tagged = true;
            if ( spec[s]->pbg == ID_PROTON)
                h_xq2_pt_tag->Fill(spec[s]->vec.Pt());
        }
            
    }

    for ( int j = 0; j < 4; j ++ )
    {
        if ( spec.size() > 0 )        
        {
            Spec1_omdHits[j] = spec[0]->det_hits[0][j];
            Spec1_rpHits[j] = spec[0]->det_hits[1][j];
        }
        else        
        {
            Spec1_rpHits[j] = 0;
            Spec1_omdHits[j] = 0;
        }

        if ( spec.size() > 1 )        
        {
            Spec2_omdHits[j] = spec[1]->det_hits[0][j];
            Spec2_rpHits[j] = spec[1]->det_hits[1][j];
        }
        else        
        {
            Spec2_rpHits[j] = 0;
            Spec2_omdHits[j] = 0;
        }
    }

    for ( auto s : spec )
    {
        SpecPBG.push_back(s->pbg);
        SpecVec.push_back(s->vec);
        SpecTag.push_back(s->tagged);
    }
    
    return;
}

void plot_ff()
{
    outFile->cd();
    std::vector<TCanvas*> canvases;

    for ( auto ff : ffFinder )
    {
        canvases = ff->draw_histograms();
        for (auto &c : canvases)
            draw_manager->LableAndCollect(c);
    }

    canvases.clear();

    canvases = zdcFinder->draw_histograms();
    for (auto &c : canvases)
        draw_manager->LableAndCollect(c);

    // canvases.clear();

    return;
}

void CreateOutputTree(TString outFileName) {

	outFile = new TFile(outFileName, "RECREATE");
	outTree = new TTree("T_Tag", "T_Tag");

    outTree->Branch("fTagged", &is_tagged);
    outTree->Branch("nTracks", &n_proton_tracks);
    outTree->Branch("E_ZDC", &zdc_energy);
    outTree->Branch("fZDCn", &fZDCn);
    outTree->Branch("SpectatorStatus", &spectator_status);
    outTree->Branch("NFinalProtons", &n_final_protons);

    outTree->Branch("Spec1_omdHits", Spec1_omdHits, "Spec1_omdHits[4]/I");
    outTree->Branch("Spec1_rpHits", Spec1_rpHits, "Spec1_rpHits[4]/I");
    outTree->Branch("Spec2_omdHits", Spec2_omdHits, "Spec2_omdHits[4]/I");
    outTree->Branch("Spec2_rpHits", Spec2_rpHits, "Spec2_rpHits[4]/I");

    outTree->Branch("SpecPBG", &SpecPBG);
    outTree->Branch("SpecVec", &SpecVec);
    outTree->Branch("SpecTag", &SpecTag);

    outTree->Branch("OtherPBG", &OtherPBG);
    outTree->Branch("OtherVec", &OtherVec);

    return;
}

// Some issues with BeAGLE/eic-smear ..
// They are not properly handling the spectator particles in the simulation.
// Made codes with ChatGPT 5.6 Sol Light for some workarounds for now
void find_spectators(const podio::Frame* event)
{
    struck_type = -1;
    spectator_status = SPECTATOR_UNCLASSIFIED;
    n_final_protons = -1;

    for (auto* p : spec)
        delete p;
    spec.clear();

    SpecPBG.clear();
    SpecVec.clear();
    OtherPBG.clear();
    OtherVec.clear();
    SpecTag.clear();

    const auto& mcHeadon = static_cast<const edm4hep::MCParticleCollection&>(*(event->get("MCParticlesHeadOnFrameNoBeamFX")));
    const auto& mc = static_cast<const edm4hep::MCParticleCollection&>(*(event->get("MCParticles")));

    // Fixed BEAGLE layout for this production.
    if (mcHeadon.size() <= 12)
        return;

    const int struckPDG = mcHeadon[12].getPDG();

    if (struckPDG == 2212)
        struck_type = 0; // ep: spectator pair is p+n
    else if (struckPDG == 2112)
        struck_type = 1; // en: spectator pair is p+p
    else
        return;

    std::vector<int> selected;
    bool usedKinematicFallback = false;

    if (struck_type == 1) {
        // en: trace each of the two original He-3 protons to its unique
        // final-state proton daughter. This avoids selecting hard protons.
        int initialProtons = 0;
        for (int i : {5, 6, 7}) {
            if (i >= static_cast<int>(mcHeadon.size()))
                continue;

            const auto initialNucleon = mcHeadon[i];
            if (initialNucleon.getPDG() != ID_PROTON)
                continue;

            ++initialProtons;
            std::vector<int> finalProtonDaughters;
            for (const auto daughter : initialNucleon.getDaughters()) {
                if (daughter.getPDG() == ID_PROTON && daughter.getGeneratorStatus() == 1)
                    finalProtonDaughters.push_back(daughter.getObjectID().index);
            }

            if (finalProtonDaughters.size() != 1)
                continue;
            selected.push_back(finalProtonDaughters[0]);
        }

        if (initialProtons != 2 || selected.size() != 2 || selected[0] == selected[1]) {
            // eic-smear can flatten generator ancestry. Use the kinematic
            // selection only as an explicitly labeled fallback.
            selected.clear();
            struct Candidate { int index; double score; };
            std::vector<Candidate> protons;
            const double pzBeam = std::abs(mcHeadon[2].getMomentum().z);
            const double dpzScale = std::max(0.30 * pzBeam, 1.0);

            for (std::size_t i = 0; i < mcHeadon.size(); ++i) {
                const auto particle = mcHeadon[i];
                if (particle.getPDG() != ID_PROTON || particle.getGeneratorStatus() != 1)
                    continue;
                const auto& p = particle.getMomentum();
                const double pt = std::hypot(p.x, p.y);
                const double score = std::pow((p.z - pzBeam) / dpzScale, 2) + std::pow(pt / 0.5, 2);
                protons.push_back({static_cast<int>(i), score});
            }

            std::sort(protons.begin(), protons.end(), [](const Candidate& a, const Candidate& b) { return a.score < b.score; });
            n_final_protons = static_cast<int>(protons.size());
            if (protons.size() < 2) {
                spectator_status = SPECTATOR_INCOMPLETE;
                return;
            }
            selected.push_back(protons[0].index);
            selected.push_back(protons[1].index);
            usedKinematicFallback = true;
        }
    }

    if (struck_type == 0) {
        // ep: save the initial proton spectator from original He-3
        // nucleons 6,7,8 (zero-based indices 5,6,7).
        std::vector<int> initialProtons;
        for (int i : {5, 6, 7}) {
            if (i >= static_cast<int>(mcHeadon.size()))
                continue;

            const auto particle = mcHeadon[i];
            if (particle.getPDG() != 2212 || particle.getGeneratorStatus() != 14)
                continue;

            initialProtons.push_back(i);
        }

        if (initialProtons.size() != 1) {
            spectator_status = SPECTATOR_INCOMPLETE;
            return;
        }

        const auto particle = mcHeadon[initialProtons[0]];

        auto* spec_info = new spectator_info{};
        spec_info->pbg = particle.getPDG();
        spec_info->mc_index = -1; // No transported MCParticle is available.
        spec_info->tagged = false;
        const TLorentzVector lab = boost_beagle_spectator(particle, mcHeadon);
        spec_info->vec.SetPxPyPzE(lab.Px(), lab.Py(), lab.Pz(), lab.E());
        spec.push_back(spec_info);

        h_xq2_pt->Fill(spec_info->vec.Pt());
        h_xq2_pt_theta->Fill(spec_info->vec.Theta() * 1000.0, spec_info->vec.Pt());

        spectator_status = SPECTATOR_EP_INITIAL;
        return;
    }

    // Validate the same-index mapping before saving either en spectator.
    for (const int index : selected) {
        if (index < 0 || static_cast<std::size_t>(index) >= mcHeadon.size() || static_cast<std::size_t>(index) >= mc.size()) {
            spectator_status = SPECTATOR_MAPPING_FAILED;
            return;
        }

        const auto headOnParticle = mcHeadon[index];
        const auto transportedParticle = mc[index];
        if (headOnParticle.getPDG() != transportedParticle.getPDG() || headOnParticle.getGeneratorStatus() != transportedParticle.getGeneratorStatus()) {
            std::cerr << "Head-on/MC mapping failed at index " << index << '\n';
            spectator_status = SPECTATOR_MAPPING_FAILED;
            return;
        }
    }

    for (const int index : selected) {
        const auto headOnParticle = mcHeadon[index];
        auto* spec_info = new spectator_info{};
        spec_info->pbg = headOnParticle.getPDG();
        spec_info->mc_index = index;
        spec_info->tagged = false;
        spec_info->vec.SetPxPyPzE(headOnParticle.getMomentum().x, headOnParticle.getMomentum().y, headOnParticle.getMomentum().z, headOnParticle.getEnergy());
        spec.push_back(spec_info);

        h_xq2_pt->Fill(spec_info->vec.Pt());
        h_xq2_pt_theta->Fill(spec_info->vec.Theta() * 1000.0, spec_info->vec.Pt());
    }

    if (!usedKinematicFallback)
        spectator_status = SPECTATOR_EN_ANCESTRY;
    else if (n_final_protons == 2)
        spectator_status = SPECTATOR_EN_KINEMATIC_CLEAN;
    else
        spectator_status = SPECTATOR_EN_KINEMATIC_AMBIGUOUS;
}

TLorentzVector boost_beagle_spectator( const edm4hep::MCParticle& spectator, const edm4hep::MCParticleCollection& mcHeadon) {
    
    // Incoming electron in the head-on frame.
    const auto beamElectron = mcHeadon[0];
    TLorentzVector k(beamElectron.getMomentum().x, beamElectron.getMomentum().y, beamElectron.getMomentum().z, beamElectron.getEnergy());

    // Find the primary final-state scattered electron.
    int scatteredIndex = -1;
    double largestElectronEnergy = -1.0;

    for (std::size_t i = 0; i < mcHeadon.size(); ++i) {
        const auto particle = mcHeadon[i];

        if (particle.getPDG() != 11 ||
            particle.getGeneratorStatus() != 1)
            continue;

        if (particle.getEnergy() > largestElectronEnergy) {
            largestElectronEnergy = particle.getEnergy();
            scatteredIndex = static_cast<int>(i);
        }
    }

    if (scatteredIndex < 0)
        return {};

    const auto scatteredElectron = mcHeadon[scatteredIndex];
    TLorentzVector kPrime(scatteredElectron.getMomentum().x, scatteredElectron.getMomentum().y, scatteredElectron.getMomentum().z, scatteredElectron.getEnergy());

    // Virtual photon in the head-on collider frame.
    TLorentzVector qLab = k - kPrime;

    const double nucleonMass = spectator.getMass();
    const double ionEnergyPerNucleon = std::sqrt(ion_momentum_per_nucleon * ion_momentum_per_nucleon + nucleonMass * nucleonMass);
    const double beta = ion_momentum_per_nucleon / ionEnergyPerNucleon;

    // Transform q and incoming electron into the nuclear rest frame.
    TLorentzVector qRest = qLab;
    TLorentzVector kRest = k;

    qRest.Boost(0.0, 0.0, -beta);
    kRest.Boost(0.0, 0.0, -beta);

    // BeAGLE local z-axis is along q in the nuclear rest frame.
    const TVector3 zAxis = qRest.Vect().Unit();

    // BeAGLE x-axis convention: opposite the component of the incoming
    // electron perpendicular to q. This convention reproduced the
    // known en spectator transverse momenta.
    TVector3 xAxis = kRest.Vect() - kRest.Vect().Dot(zAxis) * zAxis;
    xAxis = -xAxis.Unit();

    // Right-handed coordinate system.
    const TVector3 yAxis = zAxis.Cross(xAxis).Unit();

    const auto& p = spectator.getMomentum();

    // Convert BeAGLE local components into nuclear-rest-frame
    // collider coordinates.
    const TVector3 pRest = p.x * xAxis + p.y * yAxis + p.z * zAxis;
    TLorentzVector result(pRest.X(), pRest.Y(), pRest.Z(), spectator.getEnergy());

    // Nuclear rest frame -> head-on collider frame.
    result.Boost(0.0, 0.0, beta);

    return result;
}
