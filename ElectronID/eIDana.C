// Find inclusive scattered electrons

#include "eIDana.h"

void eIDana(int Ee, int Eh, int select_region, int sr, int is_truth_eID, int all_file, int analyse_p)
{
    std::cout << "** Analysing inclusive electrons, energy is set to: " << Ee << "x" << Eh << std::endl;

    // Standard setup

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

    // Analysis loop

    for( size_t ev = 0; ev < reader.getEntries("events"); ev++ )
    {
        const auto event = podio::Frame(reader.readNextEntry("events"));
        eFinder->SetEvent(&event);

        if(ev%100==0) 
        cout << "Analysing event " << ev << "/" << reader.getEntries("events") << std::endl;

        // Generator information (mcID)
        edm4hep::MCParticleCollection e_mc = eFinder->GetMCElectron();

        // Use MC to find reconstructed electron (TruthID)
        auto e_truth = eFinder->GetTruthReconElectron();

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

            // h_eID_status->Fill(eID_status);

            // if ( eID_status != 0 && eID_status != -211 )
            //     std::cout << "EID Status: " << eID_status << " id " << std::endl;
        }	

        // Calculate kinematic variables using MC electron
		TLorentzVector kprime;
		kprime.SetXYZM(e_mc[0].getMomentum().x, e_mc[0].getMomentum().y, e_mc[0].getMomentum().z, MASS_ELECTRON);
		CalculateElectronKinematics(Ee, Eh, kprime, mc_xB, mc_Q2, mc_W2, mc_y, mc_nu);

        outTree->Fill();
        ResetVariables();
    }

    // Save

    outFile->cd();
    outTree->Write(outTree->GetName(), 2);

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

	eID_status = NO_FOUND;
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