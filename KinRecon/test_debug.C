#include "kRecAna.h"

void test_debug()
{
    const std::string group[3] = {"1to10", "10to100", "100to1000"};
    
    for ( int i = 0; i < 3; i ++ )
    {
        std::cout << "Loading group " << group[i] << std::endl;
        
        TFile* beamFile = new TFile(Form("../data/en_25_10_2/root_files/10x166_%s_beam_combined.root", group[i].c_str()));
        TFile* eidFile  = new TFile(Form("../data/en_25_10_2/root_files/10x166_%s_eIDrecon_combined.root", group[i].c_str()));
        
        if (!beamFile || beamFile->IsZombie()) {
            std::cout << "ERROR: Cannot open beam file" << std::endl;
            continue;
        }
        if (!eidFile || eidFile->IsZombie()) {
            std::cout << "ERROR: Cannot open eid file" << std::endl;
            continue;
        }
        
        std::cout << "Files opened successfully" << std::endl;
        
        TTreeReader beam_reader("T_Beam", beamFile);
        TTreeReader eid_reader("T_eID", eidFile);
        
        std::cout << "Beam entries: " << beam_reader.GetEntries() << std::endl;
        std::cout << "EID entries: " << eid_reader.GetEntries() << std::endl;
        
        beamFile->Close();
        eidFile->Close();
        delete beamFile;
        delete eidFile;
    }
}
