
#include "../GlobalUtil/preLoadLib.hh"
#include "tagAna.h"

void tagAna(int Ee, int Eh, int beam_type, int select_region=0, int sr=0, int file0=-1)
{
    std::cout << "** Analysing far-forwards, energy is set to: " << Ee << "x" << Eh << std::endl;

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
    ana_manager->GetInputNames();
    // reader.openFiles(ana_manager->GetLocalInputNames());
    reader.openFiles(ana_manager->GetInputNames());

    // .. output setup;
    CreateOutputTree(ana_manager->GetOutputName()); 

    // .. FarForward setup
    setup_omd();
    setup_rp();
    setup_zdc();

    // Analysis loop

    for( size_t ev = 0; ev < reader.getEntries("events"); ev++ )
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
        process_ff(&event);
        outTree->Fill();
            
        // ResetVariables();
    }

    plot_ff();

    outFile->cd();
    outTree->Write(outTree->GetName(), 2);

    // outFile->Close();

    return;
}

void find_spectators(const podio::Frame* event)
{
    spec.clear();

    const auto& mcparts = static_cast<const edm4hep::MCParticleCollection&>(*(event->get("MCParticles")));
    for ( auto mcp : mcparts )
    {
        if ( mcp.getGeneratorStatus() == 1 )
        {
            if( mcp.getPDG() == ID_NEUTRON || mcp.getPDG() == ID_PROTON )
            {
                // cout << mcp << endl;

                int count_id = 0;
                int collect_id[2] = {0,0};
                for (auto it = mcp.parents_begin(), end = mcp.parents_end(); it != end; ++it) 
                {
                    if ( count_id < 2)
                        collect_id[count_id] = it->getObjectID().index;

                    count_id ++;
                }

                if( collect_id[0] == 3 && collect_id[1] == 4 ) // spectators from He3
                {
                    spectator_info* spec_info = new spectator_info;
                    spec_info->pbg = mcp.getPDG();
                    spec_info->mc_index = mcp.getObjectID().index;
                    spec_info->vec.SetPxPyPzE(mcp.getMomentum().x, mcp.getMomentum().y, mcp.getMomentum().z, mcp.getEnergy());
                    spec.push_back(spec_info);
                }
            }
        }
    }

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
    double xShift[4] = {-941.0, -941.0, -1032.0, -1032.0};
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
    zdcFinder = new FarForward("zdc", "ReconstructedFarForwardZDCNeutrals");    
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
        int index[2] = {spec[0]->mc_index, spec[1]->mc_index};
        ffFinder[i]->SetEvent(event);
        ffFinder[i]->GetHits(index);

        for ( int s = 0; s < 2; s ++ )
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

    if (n_proton_tracks >= 2)
        is_tagged = true;

    if (zdcFinder->is_ZDC_neutron && zdc_energy > 40.0)
    {
        fZDCn = true;

        if ( is_tagged )
            is_tagged = false;
    }

    // mc tagging
    for ( int s = 0; s < 2; s ++ )
    {
        int n_rp_mc = *std::min_element(spec[s]->det_hits[0], spec[s]->det_hits[0] + 4);
        int n_omd_mc = *std::min_element(spec[s]->det_hits[1], spec[s]->det_hits[1] + 4);
        if ( n_rp_mc + n_omd_mc > 0)
            spec[s]->tagged = true;
    }

    for ( int j = 0; j < 4; j ++ )
    {
        Spec1_rpHits[j] = spec[0]->det_hits[0][j];
        Spec1_omdHits[j] = spec[0]->det_hits[1][j];
        Spec2_rpHits[j] = spec[1]->det_hits[0][j];
        Spec2_omdHits[j] = spec[1]->det_hits[1][j];
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

    draw_manager->SaveToTree(outFile);

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

    outTree->Branch("Spec1_rpHits", Spec1_rpHits, "Spec1_rpHits[4]/I");
    outTree->Branch("Spec1_omdHits", Spec1_omdHits, "Spec1_omdHits[4]/I");
    outTree->Branch("Spec2_rpHits", Spec2_rpHits, "Spec2_rpHits[4]/I");
    outTree->Branch("Spec2_omdHits", Spec2_omdHits, "Spec2_omdHits[4]/I");

    return;
}