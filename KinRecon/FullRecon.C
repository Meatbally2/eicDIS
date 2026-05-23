#include "../GlobalUtil/DrawManager.cc"
#include "../GlobalUtil/getBoost.h"
#include "../GlobalUtil/luminosityTable.h"

#include "kinematics.cc"
#include "reconMethod.hh"

#include <Math/LorentzVector.h>
using ROOT::Math::PxPyPzEVector;

#include <Math/LorentzRotation.h>
using ROOT::Math::LorentzRotation;

enum {X,Y,Q2};

const int n_group = 3;
// const std::string group[3] = {"1to10", "10to100", "100to1000"};
const std::string group[3] = {"minQ2=1", "minQ2=10", "minQ2=100"};

void load_eID_selection(TH2F* &h_eID, string setting)
{
    std::ifstream eID_table(Form("eHe3_10x166_eID_selectrion.txt"));

    std::string line;
    while ( getline(eID_table, line) )
    {
        double x = 0;
        double q2 = 0;
        int id = 0;

        std::stringstream ss(line);
        ss >> x >> q2 >> id;

        int xbin = h_eID->GetXaxis()->FindBin(x);
        int qbin = h_eID->GetYaxis()->FindBin(q2);

        h_eID->SetBinContent(xbin, qbin, id);
        // std::cout << x << " " << q2 << " " << f2 << " " << f1 << std::endl;
    }
    eID_table.close();

    return;
}

void FullRecon(int beam_type, int Ee, int Eh)
{
    // ECCE binning
    // set_ECCE_bin();

    // ePIC plotting style setup
    std::string type_title[3] = {"e^{3}He", "ep", "#gammap"};
    std::string energy_title = beam_type ? Form("%dx%d GeV", Ee, Eh) : Form("%dx%d GeV/A", Ee, Eh);
    std::string campaign = beam_type ? "25.10.0" : "25.10.2";
    DrawManager* draw_manager = new DrawManager(type_title[beam_type], energy_title, campaign);
    draw_manager->SetEPIC();

    double text_lumi = 10; // fb^-1
    if ( beam_type == 0 && Ee == 10 && Eh == 166 )
        text_lumi = 1.5; // fb^-1
    draw_manager->SetLumi(text_lumi);

    double total_lumi = text_lumi; // fb^-1

    //

    TH2F*h_bin_eID = BookTH2(Form("h_xq2_eID"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    load_eID_selection(h_bin_eID, Form("%dx%d",(int)Ee, (int)Eh));

    TH2F* h_xq_mc = BookTH2(Form("H_XQ_MC"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    TH2F* h_xq_rec = BookTH2(Form("H_XQ_REC"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    TH2F* h_xq_he3 = BookTH2(Form("H_XQ_HE3"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    TH2F* h_xq_stb = BookTH2(Form("H_XQ_STB"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);

    for ( int i = 0; i < n_group; i ++ )
    {
        TH2F* h_xq_mc_tmp[2] = {BookTH2(Form("H_XQ_MC%d_n", i), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature), BookTH2(Form("H_XQ_MC%d_p", i), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature)};
        TH2F* h_xq_rec_tmp[2] = {BookTH2(Form("H_XQ_REC%d_n", i), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature), BookTH2(Form("H_XQ_REC%d_p", i), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature)};
        TH2F* h_xq_he3_tmp[2] = {BookTH2(Form("H_XQ_HE3%d_n", i), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature), BookTH2(Form("H_XQ_HE3%d_p", i), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature)};
        TH2F* h_xq_stb_tmp[2] = {BookTH2(Form("H_XQ_STB%d_n", i), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature), BookTH2(Form("H_XQ_STB%d_p", i), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature)};

        // load files and trees
        // TFile* beamFile = new TFile(Form("../data/en_25_10_2/root_files/10x166_%s_beam_combined.root", group[i].c_str()));
        // TFile* tagFile  = new TFile(Form("../data/en_25_10_2/root_files/10x166_%s_tag_combined.root", group[i].c_str()));
        // TFile* eidFile  = new TFile(Form("../data/en_25_10_2/root_files/10x166_%s_eIDrecon_combined.root", group[i].c_str()));
        
        TFile* beamFile = new TFile(Form("../data/en_25_10_2/root_files/eHe3_%dx%d_%s_beam_combined.root", Ee, Eh, group[i].c_str()));
        TFile* tagFile  = new TFile(Form("../data/en_25_10_2/root_files/eHe3_%dx%d_%s_tag_combined.root", Ee, Eh, group[i].c_str()));
        // TFile* eidFile  = new TFile(Form("../data/en_25_10_2/root_files/eHe3_%dx%d_%s_eIDrecon_combined.root", Ee, Eh, group[i].c_str()));
        TFile* eidFile  = new TFile(Form("../data/en_25_10_2/root_files/eHe3_%dx%d_%s_eid_combined.root", Ee, Eh, group[i].c_str()));

        TTreeReader beam_reader("T_Beam", beamFile);
        TTreeReaderValue<int> N_PDG(beam_reader, "N_PDG");
        TTreeReaderValue<PxPyPzEVector> beam_e(beam_reader, "vectE");
        TTreeReaderValue<PxPyPzEVector> beam_n(beam_reader, "vectN");

        TTreeReader tag_reader("T_Tag", tagFile);
        TTreeReaderValue<bool> fTagged(tag_reader, "fTagged");

        TTreeReader eid_reader("T_eID", eidFile);
        TTreeReaderValue<int>    status(eid_reader, "eID_status");
        TTreeReaderValue<double> mc_xB(eid_reader, "mc_xB");
        TTreeReaderValue<double> mc_Q2(eid_reader, "mc_Q2");
        TTreeReaderValue<double> mc_y(eid_reader, "mc_y");
        TTreeReaderValue<double> mc_W2(eid_reader, "mc_W2");
        TTreeReaderValue<double> mc_nu(eid_reader, "mc_nu");

        TTreeReaderValue<PxPyPzEVector> vMCe(eid_reader, "vMC_e");
        TTreeReaderValue<PxPyPzEVector> vTRe(eid_reader, "vTRACK_e");
        TTreeReaderValue<PxPyPzEVector> vCLe(eid_reader, "vCLUSTER_e");
        TTreeReaderValue<std::vector<PxPyPzEVector>> vMCf(eid_reader, "vMC_hfs");
        TTreeReaderValue<std::vector<PxPyPzEVector>> vTRf(eid_reader, "vREC_hfs");

        int en_num = 0;
        int ep_num = 0;

        Long64_t nEntries = beam_reader.GetEntries();
        for ( int ev = 0; ev < nEntries; ev ++ )
        {
            beam_reader.Next();
            tag_reader.Next();
            eid_reader.Next();
            
            if ( *status < FOUND_E || *N_PDG == 0 )
                continue;

            // if ( *status < FOUND_TRUTH || *N_PDG == 0 )
            //     continue;

            int fill_index = -1;
            if ( *N_PDG == ID_PROTON )
            {
                ep_num ++;
                fill_index = 1;
            }
            else if ( *N_PDG == ID_NEUTRON )
            {
                en_num ++;
                fill_index = 0;
            }
                
            LorentzRotation true_boost = *N_PDG == ID_PROTON ? getBoost( Ee, Eh, MASS_ELECTRON, MASS_PROTON) : getBoost( Ee, Eh, MASS_ELECTRON, MASS_NEUTRON);
            PxPyPzEVector boost_beam_e = true_boost(*beam_e);
            PxPyPzEVector boost_beam_n = true_boost(*beam_n);
            PxPyPzEVector boosted_vMCe = true_boost(*vMCe);
            // PxPyPzEVector boosted_vMCe = *vMCe;

            LorentzRotation boost = getBoost( Ee, Eh, MASS_ELECTRON, MASS_NEUTRON);
            PxPyPzEVector boosted_vTRe;
            PxPyPzEVector boosted_vCLe;

            if ( vTRe->E() > 0 )
                boosted_vTRe = boost(*vTRe);
            else
                continue;

            if ( vCLe->E() > 0 )
                boosted_vCLe = boost(*vCLe);

            double pxsum = 0;
            double pysum = 0;
            double pzsum = 0;
            double Esum  = 0;
            for ( const auto& p : *vTRf )
            {
                PxPyPzEVector boosted_p = boost(p);
                pxsum += boosted_p.Px();
                pysum += boosted_p.Py();
                pzsum += boosted_p.Pz();
                Esum  += boosted_p.E();
            }
            double pt_had = sqrt(pxsum*pxsum + pysum*pysum);
            double sigma_h = Esum - pzsum;

            std::vector<float> MCreco = calc_elec_method(boosted_vMCe.E(), boosted_vMCe.theta(), 0, 0, boost_beam_e.E(), boost_beam_n.E());

            if (*mc_y > 0.01 && *mc_y < 0.95 && *mc_W2 > 4 && *mc_Q2 > 2 && (*fTagged == 1 || beam_type > 0)) 
            // if (*mc_y > 0.01 && *mc_y < 0.95 && *mc_W2 > 4 && *mc_Q2 > 2) 
                h_xq_mc_tmp[fill_index]->Fill(MCreco[X], MCreco[Q2]);

            std::vector<float> reco;
            std::vector<float> Ereco = calc_elec_method(boosted_vTRe.E(), boosted_vTRe.theta(), pt_had, sigma_h, Ee, Eh);

            if ( Ereco[Q2] < 10 && boosted_vCLe.E() > 0)
            {
                Ereco.clear();
                Ereco = calc_elec_method(boosted_vCLe.E(), boosted_vCLe.theta(), pt_had, sigma_h, Ee, Eh);
            }
            
            reco = Ereco;
            if ( Ereco[Y] > 0.045 || vTRf->size() == 0)
                reco = Ereco;
            else
            {
                reco = calc_da_method(boosted_vTRe.E(), boosted_vTRe.theta(), pt_had, sigma_h, Ee, Eh);
                if ( reco[Q2] < 10 && boosted_vCLe.E() > 0)
                {
                    reco.clear();
                    reco = calc_da_method(boosted_vCLe.E(), boosted_vCLe.theta(), pt_had, sigma_h, Ee, Eh);
                }
            }
                
            if (reco[Y] < 0.01 || reco[Y] > 0.95 || reco[Q2] < 2)
                continue;

            double W2 = MASS_NEUTRON*MASS_NEUTRON + reco[Q2] * (1./reco[X] - 1);
            if ( *mc_W2 < 4) 
                continue;

             h_xq_he3_tmp[fill_index]->Fill(reco[X], reco[Q2]);

            if ( beam_type == 0)
                if (*fTagged != 1)
                    continue;
            
            h_xq_rec_tmp[fill_index]->Fill(reco[X], reco[Q2]);

            if ( h_xq_stb_tmp[fill_index]->GetXaxis()->FindBin(MCreco[X]) == h_xq_stb_tmp[fill_index]->GetXaxis()->FindBin(reco[X]) )
                if ( h_xq_stb_tmp[fill_index]->GetYaxis()->FindBin(MCreco[Q2]) == h_xq_stb_tmp[fill_index]->GetYaxis()->FindBin(reco[Q2]) )
                    h_xq_stb_tmp[fill_index]->Fill(reco[X], reco[Q2]);
        }

        std::cout << "En: " << en_num << ", Ep: " << ep_num << std::endl;

        double gen_lumi = get_lumi(0, Ee, Eh, i, en_num, 0);
        h_xq_mc_tmp[0]->Scale( total_lumi / gen_lumi );
        h_xq_rec_tmp[0]->Scale( total_lumi / gen_lumi );
        h_xq_he3_tmp[0]->Scale( total_lumi / gen_lumi );
        h_xq_stb_tmp[0]->Scale( total_lumi / gen_lumi );

        gen_lumi = get_lumi(0, Ee, Eh, i, 0, ep_num);
        h_xq_mc_tmp[1]->Scale( total_lumi / gen_lumi );
        h_xq_rec_tmp[1]->Scale( total_lumi / gen_lumi );
        h_xq_he3_tmp[1]->Scale( total_lumi / gen_lumi );
        h_xq_stb_tmp[1]->Scale( total_lumi / gen_lumi );
        
        for ( int fill = 0; fill < 2; fill ++ )
        {
            h_xq_mc->Add(h_xq_mc_tmp[fill]);
            h_xq_rec->Add(h_xq_rec_tmp[fill]);
            h_xq_he3->Add(h_xq_he3_tmp[fill]);
            h_xq_stb->Add(h_xq_stb_tmp[fill]);
        }
    }

    set_2d_scale(h_xq_mc);
    TCanvas* c_xq2_mc = draw_2d_standard(h_xq_mc, "c_xq2_mc", "mc events", 700, 600, true, true);
    TCanvas* c_xq2_he3 = draw_2d_standard(h_xq_he3, "c_xq2_he3", "he3 events", 700, 600, true, true);
    TCanvas* c_xq2_rec = draw_2d_standard(h_xq_rec, "c_xq2_rec", "rec events", 700, 600, true, true);

    TH2F* h_xq_stb_copy = (TH2F*)h_xq_stb->Clone();
    process_eff_hist(h_xq_stb_copy, h_xq_mc);
    TCanvas* c_xq2_stb = draw_2d_efficiency(h_xq_stb_copy, "c_bin_stb", "bin stb", 1400, 600, false, true);

    draw_manager->LableAndCollect(c_xq2_mc);
    draw_manager->LableAndCollect(c_xq2_he3);
    draw_manager->LableAndCollect(c_xq2_rec);
    draw_manager->LableAndCollect(c_xq2_stb);

    // c_xq2_mc->SaveAs(Form("../data/en_25_10_2/FullRecon_mc_%dx%d_he3.png", (int)Ee, (int)Eh));
    // c_xq2_rec->SaveAs(Form("../data/en_25_10_2/FullRecon_rec_%dx%d_cuts.png", (int)Ee, (int)Eh));
    // c_xq2_stb->SaveAs(Form("../data/en_25_10_2/FullRecon_binStability_%dx%d.png", (int)Ee, (int)Eh));

    TFile*outFile = new TFile(Form("../data/en_25_10_2/FullRecon_%dx%d_new.root", (int)Ee, (int)Eh), "RECREATE");
    // TFile*outFile = new TFile(Form("../data/en_25_10_2/FullRecon_%dx%d_ECCEbins.root", (int)Ee, (int)Eh), "RECREATE");
    outFile->cd();
    draw_manager->SaveToTree(outFile);
    h_xq_mc->Write();
    h_xq_he3->Write();
    h_xq_rec->Write();
    h_xq_stb->Write();
}