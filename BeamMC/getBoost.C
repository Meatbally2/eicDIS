#include "../GlobalUtil/Beam.h"
#include "../GlobalUtil/Boost.h"
#include "../GlobalUtil/Constants.hh"

#include <Math/LorentzVector.h>
using ROOT::Math::PxPyPzEVector;

#include <Math/LorentzRotation.h>
using ROOT::Math::LorentzRotation;

void getBoost() {

    TFile* file = new TFile(Form("18x275_beam.root"));

    TTreeReader reader("T_Beam", file);
    TTreeReaderValue<int> N_PDG(reader, "N_PDG");
    TTreeReaderValue<PxPyPzEVector> vectE(reader, "vectE");
    TTreeReaderValue<PxPyPzEVector> vectN(reader, "vectN");

    Long64_t nEntries = reader.GetEntries();

    PxPyPzEVector avg_e(0,0,0,0);
    PxPyPzEVector avg_n(0,0,0,0);
    for( size_t ev = 0; ev < nEntries; ev++ ) 
    {
        reader.Next();
        // if(ev%100==0) 
        // cout << "Analysing file " << i << " event " << ev << "/" << nEntries << "\t\r" << std::flush;

        int n_pdg = *N_PDG;
        PxPyPzEVector beam_e = *vectE;
        PxPyPzEVector beam_n = *vectN;

        avg_e += beam_e;
        avg_n += beam_n;
    }

    avg_e /= nEntries;
    avg_n /= nEntries;

    std::cout << "Average electron beam: " << Form("(%f, %f, %f, %f)", avg_e.Px(), avg_e.Py(), avg_e.Pz(), avg_e.E()) << std::endl;
    std::cout << "Average nucleon beam: " << Form("(%f, %f, %f, %f)", avg_n.Px(), avg_n.Py(), avg_n.Pz(), avg_n.E()) << std::endl;
    std::cout << "Original lab frame, e, p (mrad) " << avg_e.Theta() * 1000 << " " << avg_n.Theta() * 1000 << std::endl;

    TVector3 ve(avg_e.Px(),avg_e.Py(),avg_e.Pz());
    TVector3 vn(avg_n.Px(),avg_n.Py(),avg_n.Pz());

    // TVector3 ve(0,0,-18);
    // TVector3 vn(0,0,sqrt(275*275-MASS_PROTON*MASS_PROTON));
    
    const PxPyPzEVector ei(
        eicrecon::round_beam_four_momentum(
            ve,
            avg_e.M(),
            {-1*18},
            0.0)
        );

    const PxPyPzEVector ni(
        eicrecon::round_beam_four_momentum(
            vn,
            avg_n.M(),
            {275},
            -0.025)
        );

    std::cout << "initial electron beam: " << Form("(%f, %f, %f, %f)", ei.Px(), ei.Py(), ei.Pz(), ei.E()) << std::endl;
    std::cout << "initial nucleon beam: " << Form("(%f, %f, %f, %f)", ni.Px(), ni.Py(), ni.Pz(), ni.E()) << std::endl;
    std::cout << "Created vector, e, p (mrad) " << ei.Theta() * 1000 << " " << ni.Theta() * 1000 << std::endl;

    LorentzRotation boost = eicrecon::determine_boost(ei, ni);

    PxPyPzEVector boosted_e = boost(avg_e);
    PxPyPzEVector boosted_n = boost(avg_n);
    std::cout << "Boosted electron beam: " << Form("(%f, %f, %f, %f)", boosted_e.Px(), boosted_e.Py(), boosted_e.Pz(), boosted_e.E()) << std::endl;
    std::cout << "Boosted nucleon beam: " << Form("(%f, %f, %f, %f)", boosted_n.Px(), boosted_n.Py(), boosted_n.Pz(), boosted_n.E()) << std::endl;
    std::cout << "Incoming headon, e, p (mrad) " << boosted_e.Theta() * 1000 << " " << boosted_n.Theta() * 1000 << std::endl;

    TFile* outFile = new TFile("18x275_boost.root", "RECREATE");
    TTree* outTree = new TTree("T_Boost", "T_Boost");
    outTree->Branch("boost", &boost);
    outTree->Fill();
    outFile->cd();
    outTree->Write(outTree->GetName(), 2);
    outFile->Close();
    file->Close();

    return;
}