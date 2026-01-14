#include "../GlobalUtil/DrawManager.cc"
#include "../GlobalUtil/Constants.hh"
#include "../GlobalUtil/drawHelper.h"
#include "../GlobalUtil/luminosityTable.h"
#include "../GlobalUtil/getBoost.h"
#include <Math/LorentzVector.h>
using ROOT::Math::PxPyPzEVector;

int n_group = 3;
const std::string group[4] = {"1", "10", "100", "1000"};

void eff(int Ee, int Eh)
{
    // ePIC plotting style setup
    std::string type_title = "e^{3}He";
    std::string energy_title = Form("%dx%d GeV", Ee, Eh);
    DrawManager* draw_manager = new DrawManager(type_title, energy_title, "25.10.2");
    draw_manager->SetEPIC();

    double text_lumi = 10; // fb^-1
    if ( Ee == 10 && Eh == 166 )
        text_lumi = 8.65; // fb^-1
    draw_manager->SetLumi(text_lumi);

    // get generated lumi
    double total_lumi = text_lumi*3; // fb^-1
    std::vector<double> gen_lumi = get_lumi(0, Ee, Eh);
    if ( gen_lumi.empty() )
    {
        std::cout << "** Lumi table not found! " << std::endl;
        return;
    }

    // histograms
    TH2F* h_xq2_MCall = BookTH2(Form("h_xq2_MCall"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    TH2F* h_xq2_MCen = BookTH2(Form("h_xq2_MCen"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    TH2F* h_xq2_DTeff = BookTH2(Form("h_xq2_DTeff"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    TH2F* h_xq2_DTpur = BookTH2(Form("h_xq2_DTpur"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    TH2F* h_xq2_ZDCeff = BookTH2(Form("h_xq2_ZDCeff"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    TH2F* h_xq2_ZDCpur = BookTH2(Form("h_xq2_ZDCpur"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    TH2F* h_xq2_FFeff = BookTH2(Form("h_xq2_FFeff"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    TH2F* h_xq2_FFpur = BookTH2(Form("h_xq2_FFpur"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);

    TH2F* h_xq2_pt_theta = new TH2F(Form("h_xq2_pt_theta"), ";#theta (mrad);p_{T} (GeV/c)", 200, 0, 20, 150, 0, 1.5);
    TH2F* h_xq2_pt_theta_acp = new TH2F(Form("h_xq2_pt_theta_acp"), ";#theta (mrad);p_{T} (GeV/c)", 200, 0, 20, 150, 0, 1.5);
    TH1F* h_xq2_pt = new TH1F(Form("h_xq2_pt"), ";p_{T} (GeV/c);Counts", 150, 0, 1.5);
    TH1F* h_xq2_theta = new TH1F(Form("h_xq2_theta"), ";#theta (mrad);Counts", 200, 0, 20);
    TH1F* h_xq2_pt_acp = new TH1F(Form("h_xq2_pt_acp"), ";p_{T} (GeV/c);Counts", 150, 0, 1.5);
    TH1F* h_xq2_theta_acp = new TH1F(Form("h_xq2_theta_acp"), ";#theta (mrad);Counts", 200, 0, 20);
    TH1F* h_xq2_pt_other = new TH1F(Form("h_xq2_pt_other"), ";p_{T} (GeV/c);Counts", 150, 0, 1.5);
    TH1F* h_xq2_theta_other = new TH1F(Form("h_xq2_theta_other"), ";#theta (mrad);Counts", 200, 0, 20);

    for ( int i = 0; i < n_group; i ++ )
    {
        TFile* beamFile = new TFile(Form("../data/en_25_10_2/root_files/eHe3_%dx%d_minQ2=%s_beam_combined.root", Ee, Eh, group[i].c_str()));
        TFile* tagFile  = new TFile(Form("../data/en_25_10_2/root_files/eHe3_%dx%d_minQ2=%s_tag_combined.root", Ee, Eh, group[i].c_str()));
        TFile* eidFile  = new TFile(Form("../data/en_25_10_2/root_files/eHe3_%dx%d_minQ2=%s_eid_combined.root", Ee, Eh, group[i].c_str()));
        
        TTreeReader beam_reader("T_Beam", beamFile);
        TTreeReaderValue<int> N_PDG(beam_reader, "N_PDG");
        TTreeReaderValue<PxPyPzEVector> beam_e(beam_reader, "vectE");
        TTreeReaderValue<PxPyPzEVector> beam_n(beam_reader, "vectN");
        TTreeReaderValue<std::vector<int>> SpecPBG(beam_reader, "SpecPBG");
        TTreeReaderValue<std::vector<PxPyPzEVector>> SpecVec(beam_reader, "SpecVec");
        TTreeReaderValue<std::vector<int>> OtherPBG(beam_reader, "OtherPBG");
        TTreeReaderValue<std::vector<PxPyPzEVector>> OtherVec(beam_reader, "OtherVec");

        TTreeReader tag_reader("T_Tag", tagFile);
        TTreeReaderValue<bool> fTagged(tag_reader, "fTagged");
        TTreeReaderValue<int> nTracks(tag_reader, "nTracks");
        TTreeReaderValue<double> E_ZDC(tag_reader, "E_ZDC");
        TTreeReaderValue<bool> fZDCn(tag_reader, "fZDCn");

        TTreeReader eid_reader("T_eID", eidFile);
        TTreeReaderValue<int>    status(eid_reader, "eID_status");
        TTreeReaderValue<double> mc_xB(eid_reader, "mc_xB");
        TTreeReaderValue<double> mc_Q2(eid_reader, "mc_Q2");
        TTreeReaderValue<double> mc_y(eid_reader, "mc_y");
        TTreeReaderValue<double> mc_W2(eid_reader, "mc_W2");
        TTreeReaderValue<double> mc_nu(eid_reader, "mc_nu");

        Long64_t nEntries = beam_reader.GetEntries();

        TH2F* h_tmp_MCall = BookTH2(Form("h_tmp_MCall_%i", i), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
        TH2F* h_tmp_MCen = BookTH2(Form("h_tmp_MCen_%i", i), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
        TH2F* h_tmp_DTeff = BookTH2(Form("h_tmp_DTeff_%i", i), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
        TH2F* h_tmp_DTpur = BookTH2(Form("h_tmp_DTpur_%i", i), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
        TH2F* h_tmp_ZDCeff = BookTH2(Form("h_tmp_ZDCeff_%i", i), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
        TH2F* h_tmp_ZDCpur = BookTH2(Form("h_tmp_ZDCpur_%i", i), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
        TH2F* h_tmp_FFeff = BookTH2(Form("h_tmp_FFeff_%i", i), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
        TH2F* h_tmp_FFpur = BookTH2(Form("h_tmp_FFpur_%i", i), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);

        TH2F* h_tmp_pt_theta = new TH2F(Form("h_tmp_pt_theta_%i", i), ";#theta (mrad);p_{T} (GeV/c)", 200, 0, 20, 150, 0, 1.5);
        TH2F* h_tmp_pt_theta_acp = new TH2F(Form("h_tmp_pt_theta_acp_%i", i), ";#theta (mrad);p_{T} (GeV/c)", 200, 0, 20, 150, 0, 1.5);
        TH1F* h_tmp_pt = new TH1F(Form("h_tmp_pt_%i", i), ";p_{T} (GeV/c);Counts", 150, 0, 1.5);
        TH1F* h_tmp_theta = new TH1F(Form("h_tmp_theta_%i", i), ";#theta (mrad);Counts", 200, 0, 20);
        TH1F* h_tmp_pt_acp = new TH1F(Form("h_tmp_pt_acp_%i", i), ";p_{T} (GeV/c);Counts", 150, 0, 1.5);
        TH1F* h_tmp_theta_acp = new TH1F(Form("h_tmp_theta_acp_%i", i), ";#theta (mrad);Counts", 200, 0, 20);
        TH1F* h_tmp_pt_other = new TH1F(Form("h_tmp_pt_other_%i", i), ";p_{T} (GeV/c);Counts", 150, 0, 1.5);
        TH1F* h_tmp_theta_other = new TH1F(Form("h_tmp_theta_other_%i", i), ";#theta (mrad);Counts", 200, 0, 20);

        for ( int ev = 0; ev < nEntries; ev ++ )
        {
            beam_reader.Next();
            tag_reader.Next();
            eid_reader.Next();

            if ( *status == NO_MC )
                continue;

            if (*mc_y < 0.01 || *mc_y > 0.95) 
                continue;

            if (*mc_W2 < 4)
                continue;

            if (*mc_Q2 < 2)
                continue;

            LorentzRotation true_boost = *N_PDG == ID_PROTON ? getBoost( Ee, Eh, MASS_ELECTRON, MASS_PROTON) : getBoost( Ee, Eh, MASS_ELECTRON, MASS_NEUTRON);

            h_tmp_MCall->Fill(*mc_xB, *mc_Q2);

            if ( *N_PDG == ID_NEUTRON )
            {
                h_tmp_MCen->Fill(*mc_xB, *mc_Q2);

                if (*nTracks >= 2)
                {
                    h_tmp_DTpur->Fill(*mc_xB, *mc_Q2);

                    if ( *fZDCn != 1 || *E_ZDC < 40.0 )
                        h_tmp_ZDCpur->Fill(*mc_xB, *mc_Q2);
                }

                if ( *fTagged == 1 )
                    h_tmp_FFpur->Fill(*mc_xB, *mc_Q2);

                for ( auto p : *SpecVec )
                {
                    PxPyPzEVector boosted_p = true_boost(p);
                    double pt = sqrt(boosted_p.Px()*boosted_p.Px() + boosted_p.Py()*boosted_p.Py());
                    double theta = boosted_p.Theta()*1000; // mrad
                    h_tmp_pt_theta->Fill(theta, pt);
                    h_tmp_pt->Fill(pt);
                    h_tmp_theta->Fill(theta);
                }

                int index = 0;
                for ( auto p : *OtherVec )
                {
                    if ( (*OtherPBG)[index] == ID_NEUTRON )
                        continue;
                    
                    PxPyPzEVector boosted_p = true_boost(p);
                    double pt = sqrt(boosted_p.Px()*boosted_p.Px() + boosted_p.Py()*boosted_p.Py());
                    double theta = boosted_p.Theta()*1000; // mrad
                    h_tmp_pt_other->Fill(pt);
                    h_tmp_theta_other->Fill(theta);

                    index ++;
                }
            }

            // if ( *fTagged == 1 )
            if (*nTracks >= 2)
            {
                h_tmp_DTeff->Fill(*mc_xB, *mc_Q2);

                if ( *fZDCn != 1 || *E_ZDC < 40.0 )
                    h_tmp_ZDCeff->Fill(*mc_xB, *mc_Q2);

                // if ( *N_PDG == ID_NEUTRON )
                {
                    for ( auto p : *SpecVec )
                    {
                        PxPyPzEVector boosted_p = true_boost(p);
                        double pt = sqrt(boosted_p.Px()*boosted_p.Px() + boosted_p.Py()*boosted_p.Py());
                        double theta = boosted_p.Theta()*1000; // mrad
                        h_tmp_pt_theta_acp->Fill(theta, pt);
                        h_tmp_pt_acp->Fill(pt);
                        h_tmp_theta_acp->Fill(theta);
                    }

                    int index = 0;
                    for ( auto p : *OtherVec )
                    {
                        if ( (*OtherPBG)[index] == ID_NEUTRON )
                        continue;

                        PxPyPzEVector boosted_p = true_boost(p);
                        double pt = sqrt(boosted_p.Px()*boosted_p.Px() + boosted_p.Py()*boosted_p.Py());
                        double theta = boosted_p.Theta()*1000; // mrad
                        h_tmp_pt_other->Fill(pt);
                        h_tmp_theta_other->Fill(theta);

                        index ++;
                    }
                }
            }    

            if ( *fTagged == 1 )
                h_tmp_FFeff->Fill(*mc_xB, *mc_Q2);
        }

        h_tmp_MCall->Scale(total_lumi/gen_lumi[i]);
        h_tmp_MCen->Scale(total_lumi/gen_lumi[i]);
        h_tmp_DTeff->Scale(total_lumi/gen_lumi[i]);
        h_tmp_DTpur->Scale(total_lumi/gen_lumi[i]);
        h_tmp_ZDCeff->Scale(total_lumi/gen_lumi[i]);
        h_tmp_ZDCpur->Scale(total_lumi/gen_lumi[i]);
        h_tmp_FFeff->Scale(total_lumi/gen_lumi[i]);
        h_tmp_FFpur->Scale(total_lumi/gen_lumi[i]);

        h_tmp_pt_theta->Scale(total_lumi/gen_lumi[i]);
        h_tmp_pt->Scale(total_lumi/gen_lumi[i]);
        h_tmp_theta->Scale(total_lumi/gen_lumi[i]);
        h_tmp_pt_acp->Scale(total_lumi/gen_lumi[i]);
        h_tmp_theta_acp->Scale(total_lumi/gen_lumi[i]);
        h_tmp_pt_other->Scale(total_lumi/gen_lumi[i]);
        h_tmp_theta_other->Scale(total_lumi/gen_lumi[i]);

        h_xq2_MCall->Add(h_tmp_MCall);
        h_xq2_MCen->Add(h_tmp_MCen);
        h_xq2_DTeff->Add(h_tmp_DTeff);
        h_xq2_DTpur->Add(h_tmp_DTpur);
        h_xq2_ZDCeff->Add(h_tmp_ZDCeff);
        h_xq2_ZDCpur->Add(h_tmp_ZDCpur);
        h_xq2_FFeff->Add(h_tmp_FFeff);
        h_xq2_FFpur->Add(h_tmp_FFpur);

        h_xq2_pt_theta->Add(h_tmp_pt_theta);
        h_xq2_pt_theta_acp->Add(h_tmp_pt_theta_acp);
        h_xq2_pt->Add(h_tmp_pt);
        h_xq2_theta->Add(h_tmp_theta);
        h_xq2_pt_acp->Add(h_tmp_pt_acp);
        h_xq2_theta_acp->Add(h_tmp_theta_acp);
        h_xq2_pt_other->Add(h_tmp_pt_other);
        h_xq2_theta_other->Add(h_tmp_theta_other);

        beamFile->Close();
        tagFile->Close();
        eidFile->Close();
    }

    TCanvas* c_spec_pt_v_th = new TCanvas("c_spec_pt_v_th", "c_spec_pt_v_th",1600, 1000);
    c_spec_pt_v_th->Divide(2,2);

    c_spec_pt_v_th->cd(1);
    h_xq2_pt_theta->SetStats(0);
    h_xq2_pt_theta->Draw("COLZ");
    h_xq2_pt_theta->GetXaxis()->CenterTitle();
    h_xq2_pt_theta->GetYaxis()->CenterTitle();
    gPad->SetLogz();

    c_spec_pt_v_th->cd(2);
    h_xq2_pt_theta_acp->SetStats(0);
    h_xq2_pt_theta_acp->Draw("COLZ");
    h_xq2_pt_theta_acp->GetXaxis()->CenterTitle();
    h_xq2_pt_theta_acp->GetYaxis()->CenterTitle();
    gPad->SetLogz();

    c_spec_pt_v_th->cd(3);
    h_xq2_pt->SetStats(0);
    h_xq2_pt->Draw("HIST");
    h_xq2_pt_acp->Draw("HIST SAME");
    h_xq2_pt_acp->SetLineColor(kBlue);
    h_xq2_pt_other->Draw("HIST SAME");
    h_xq2_pt_other->SetLineColor(kRed);
    h_xq2_pt->GetXaxis()->CenterTitle();
    h_xq2_pt->GetYaxis()->CenterTitle();
    gPad->SetLogy();

    c_spec_pt_v_th->cd(4);
    h_xq2_theta->SetStats(0);
    h_xq2_theta->Draw("HIST");
    h_xq2_theta_acp->Draw("HIST SAME");
    h_xq2_theta_acp->SetLineColor(kBlue);
    h_xq2_theta_other->Draw("HIST SAME");
    h_xq2_theta_other->SetLineColor(kRed);
    h_xq2_theta->GetXaxis()->CenterTitle();
    h_xq2_theta->GetYaxis()->CenterTitle();
    gPad->SetLogy();

    for ( int i = 0; i < 4; i ++ )
    {
        c_spec_pt_v_th->cd(i+1);
        gPad->SetRightMargin(0.15);
        gPad->SetTickx(1);
        gPad->SetTicky(1);
    }

    // draw_manager->LableAndCollect(c_spec_pt_v_th,1);

    set_2d_scale(h_xq2_MCall);
    TCanvas* c_xq2_MCall = draw_2d_standard(h_xq2_MCall, "c_xq2_MCall", "all MC events", 700, 600, true, true);
    TCanvas* c_xq2_MCen = draw_2d_standard(h_xq2_MCen, "c_xq2_MCen", "all n events", 700, 600, true, true);
    TCanvas* c_xq2_DTeff = draw_2d_standard(h_xq2_DTeff, "c_xq2_DTeff", "tagged eff events", 700, 600, true, true);
    TCanvas* c_xq2_DTpur = draw_2d_standard(h_xq2_DTpur, "c_xq2_DTpur", "tagged purity events", 700, 600, true, true);
    TCanvas* c_xq2_ZDCeff = draw_2d_standard(h_xq2_ZDCeff, "c_xq2_ZDCeff", "ZDC eff events", 700, 600, true, true);
    TCanvas* c_xq2_ZDCpur = draw_2d_standard(h_xq2_ZDCpur, "c_xq2_ZDCpur", "ZDC purity events", 700, 600, true, true);
    TCanvas* c_xq2_FFeff = draw_2d_standard(h_xq2_FFeff, "c_xq2_FFeff", "FF eff events", 700, 600, true, true);
    TCanvas* c_xq2_FFpur = draw_2d_standard(h_xq2_FFpur, "c_xq2_FFpur", "FF purity events", 700, 600, true, true);

    TH2F* h_xq2_MCep = (TH2F*)h_xq2_MCall->Clone();
    h_xq2_MCep->Add(h_xq2_MCen,-1);

    TH2F* h_xq2_MCen_copy = (TH2F*)h_xq2_MCen->Clone();
    process_eff_hist(h_xq2_MCen_copy, h_xq2_MCall);
    TCanvas* c_xq2_fen = draw_2d_efficiency(h_xq2_MCen_copy, "c_fen", "frac en", 1400, 600, false, true);

    TH2F* h_xq2_DTeff_copy = (TH2F*)h_xq2_DTeff->Clone();
    process_eff_hist(h_xq2_DTeff_copy, h_xq2_MCen);
    TCanvas* c_xq2_DT_eff = draw_2d_efficiency(h_xq2_DTeff_copy, "c_DT_eff", "DT eff", 1400, 600, false, true);

    TH2F* h_xq2_DTpur_copy = (TH2F*)h_xq2_DTpur->Clone();
    process_eff_hist(h_xq2_DTpur_copy, h_xq2_DTeff);
    TCanvas* c_xq2_DT_pur = draw_2d_efficiency(h_xq2_DTpur_copy, "c_DT_pur", "DT pur", 1400, 600, false, true);

    TH2F* h_xq2_ZDCeff_copy = (TH2F*)h_xq2_ZDCeff->Clone();
    process_eff_hist(h_xq2_ZDCeff_copy, h_xq2_DTeff);
    TCanvas* c_xq2_ZDC_eff = draw_2d_efficiency(h_xq2_ZDCeff_copy, "c_ZDC_eff", "ZDC eff", 1400, 600, false, true);

    TH2F* h_xq2_ZDCpur_copy = (TH2F*)h_xq2_ZDCpur->Clone();
    process_eff_hist(h_xq2_ZDCpur_copy, h_xq2_ZDCeff);
    TCanvas* c_xq2_ZDC_pur = draw_2d_efficiency(h_xq2_ZDCpur_copy, "c_ZDC_pur", "ZDC pur", 1400, 600, false, true);

    TH2F* h_xq2_FFeff_copy = (TH2F*)h_xq2_FFeff->Clone();
    process_eff_hist(h_xq2_FFeff_copy, h_xq2_MCen);
    TCanvas* c_xq2_FF_eff = draw_2d_efficiency(h_xq2_FFeff_copy, "c_FF_eff", "FF eff", 1400, 600, false, true);

    TH2F* h_xq2_FFpur_copy = (TH2F*)h_xq2_FFpur->Clone();
    process_eff_hist(h_xq2_FFpur_copy, h_xq2_FFeff);
    TCanvas* c_xq2_FF_pur = draw_2d_efficiency(h_xq2_FFpur_copy, "c_FF_pur", "FF pur", 1400, 600, false, true);

    draw_manager->LableAndCollect(c_xq2_MCall);
    draw_manager->LableAndCollect(c_xq2_MCen);
    draw_manager->LableAndCollect(c_xq2_fen);
    draw_manager->LableAndCollect(c_xq2_DTeff);
    draw_manager->LableAndCollect(c_xq2_DTpur);
    draw_manager->LableAndCollect(c_xq2_ZDCeff);
    draw_manager->LableAndCollect(c_xq2_ZDCpur);
    draw_manager->LableAndCollect(c_xq2_DT_eff);
    draw_manager->LableAndCollect(c_xq2_DT_pur);
    draw_manager->LableAndCollect(c_xq2_ZDC_eff);
    draw_manager->LableAndCollect(c_xq2_ZDC_pur);   
    draw_manager->LableAndCollect(c_xq2_FFeff);
    draw_manager->LableAndCollect(c_xq2_FFpur);

    c_xq2_DT_eff->SaveAs(Form("../data/eID/ff_eff_%dx%d_e3He_25.10.2.png", Ee, Eh));
    c_xq2_DT_pur->SaveAs(Form("../data/eID/ff_pur_%dx%d_e3He_25.10.2.png", Ee, Eh));
}