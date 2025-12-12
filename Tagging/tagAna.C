#include "tagAna.h"

void tagAna(int Ee, int Eh, int select_region, int sr, int all_file, int analyse_p)
{
    std::cout << "** Analysing far-forwards, energy is set to: " << Ee << "x" << Eh << std::endl;

    // Standard setup

    SetePICStyle();

    AnaManager* ana_manager = new AnaManager("tag");
    ana_manager->Initialize(select_region, sr, all_file, analyse_p);

    // .. input setup
    auto reader = podio::ROOTReader();
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

void setup_rp()
{
    rpFinder = new FarForward("rp", "ForwardRomanPotRecHits");
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
    omdFinder = new FarForward("omd", "ForwardOffMTrackerRecHits");    
    omdFinder->define_histograms();

    double rotate = 1./cos(-0.047);
    // double xShift[4] = {-780.0, -780.0, -870.0, -870.0};
    // double zRange[4][2] = {{22490, 22512}, {22512, 22528}, {24490, 24510},{24510, 24530}};
    double xShift[4] = {-941.0, -941.0, -1032.0, -1032.0};
    double zRange[4][2] = {{25500, 25512}, {25515, 25528}, {26950, 27010}, {27015, 27030}};

    // for ( int i = 0; i < 4; i ++ )
    //     for ( int j = 0; j < 2; j ++ )
    //         zRange[i][j] += 3000; // shift by 1m downstream

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
    zdcFinder = new FarForward("zdc", "HcalFarForwardZDCClusters"); 
    // zdcFinder = new FarForward("zdc", "ReconstructedFarForwardZDCNeutrals");    
    zdcFinder->define_histograms();
    return;
}

void process_ff(const podio::Frame* event)
{
    is_tagged = false;
    n_proton_tracks = 0;
    zdc_energy = 0.0;

    for ( auto ff : ffFinder )
    {
        ff->SetEvent(event);
        ff->GetHits();
        ff->fill_hit_histograms();
    }

    n_proton_tracks = rpFinder->get_n_tracks() + omdFinder->get_n_tracks();

    if (rpFinder->get_n_tracks() + omdFinder->get_n_tracks() >= 2)
        is_tagged = true;

    zdcFinder->SetEvent(event);
    zdc_energy = zdcFinder->GetEnergy();
    zdcFinder->fill_energy_histograms();

    return;
}

void plot_ff()
{
    outFile->cd();
    std::vector<TCanvas*> canvases;

    for ( auto ff : ffFinder )
    {
        canvases = ff->draw_histograms();
        
        for (auto c : canvases)
            c->Write();

        canvases.clear();
    }
        
    canvases = zdcFinder->draw_histograms();
    for (auto c : canvases)
        c->Write();

    canvases.clear();
    
    return;
}

void CreateOutputTree(TString outFileName) {

	outFile = new TFile(outFileName, "RECREATE");
	outTree = new TTree("T_Tag", "T_Tag");

    outTree->Branch("fTagged", &is_tagged);
    outTree->Branch("nTracks", &n_proton_tracks);
    outTree->Branch("E_ZDC", &zdc_energy);

    return;
}