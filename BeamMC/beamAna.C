// Find inclusive scattered electrons

#include "../GlobalUtil/preLoadLib.hh"
#include "beamAna.h"

void beamAna(int Ee, int Eh, int beam_type, int select_region=0, int sr=0, int file0=-1)
{
    std::cout << "** Analysing incoming beam, energy is set to: " << Ee << "x" << Eh << std::endl;

    // Standard setup

    AnaManager* ana_manager = new AnaManager("beam");
    ana_manager->Initialize(select_region, sr, file0, beam_type);
    ana_manager->SetBeamEnergy(Ee, Eh);

    // .. input setup
    auto reader = podio::ROOTReader();
    reader.openFiles(ana_manager->GetInputNames());

    // .. output setup;
    CreateOutputTree(ana_manager->GetOutputName()); 

    // .. Beam class setup
    BeamMC* inFinder = new BeamMC(Ee, Eh);

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
        inFinder->SetEvent(&event);

        if(ev%100==0) 
            cout << "Analysing event " << ev << "/" << reader.getEntries("events") << std::endl;

        inFinder->GetMCinfo(vectE, vectN, N_PDG);
        inFinder->GetSpecInfo(SpecPBG, SpecIndex, SpecVec, OtherPBG, OtherIndex, OtherVec);
        
        outTree->Fill();
        ResetVariables();
    }

    outFile->cd();
    outTree->Write(outTree->GetName(), 2);

    return;
}

void DefineHistograms() {

    return;
}


void CreateOutputTree(TString outFileName) {

	outFile = new TFile(outFileName, "RECREATE");
	outTree = new TTree("T_Beam", "T_Beam");

    outTree->Branch("N_PDG", &N_PDG);
    outTree->Branch("vectE", &vectE);
	outTree->Branch("vectN", &vectN);

    outTree->Branch("SpecPBG", &SpecPBG);
    outTree->Branch("SpecIndex", &SpecIndex);
    outTree->Branch("SpecVec", &SpecVec);

    outTree->Branch("OtherPBG", &OtherPBG);
    outTree->Branch("OtherIndex", &OtherIndex);
    outTree->Branch("OtherVec", &OtherVec);

    return;
}

void ResetVariables() {

	N_PDG = 0;

	vectE.SetPxPyPzE(0, 0, 0, 0);
	vectN.SetPxPyPzE(0, 0, 0, 0);

    SpecPBG.clear();
    SpecVec.clear();
    SpecIndex.clear();

    OtherPBG.clear();
    OtherVec.clear();
    OtherIndex.clear();

	return;
}

