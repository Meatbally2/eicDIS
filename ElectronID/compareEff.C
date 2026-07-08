#include "../GlobalUtil/DrawManager.cc"
#include "../GlobalUtil/drawHelper.h"
#include "../GlobalUtil/Constants.hh"

// void mask_zero_bins(TH2F* hist)
// {
//     for (int ix = 1; ix <= hist->GetNbinsX(); ++ix) {
//         for (int iy = 1; iy <= hist->GetNbinsY(); ++iy) {
//             double value = hist->GetBinContent(ix, iy);
//             if (std::abs(value) < 1e-12)
//                 hist->SetBinContent(ix, iy, TMath::QuietNaN());
//         }
//     }
// }

void compareEff()
{
    gStyle->SetPaintTextFormat(".1f");

    // ePIC plotting style setup
    DrawManager* draw_manager = new DrawManager("ep", "10x100 GeV", "26.04.1");
    draw_manager->SetEPIC("Work in Progress");
    draw_manager->SetDISrange(0.01, 0.95, 4, 2);
    draw_manager->SetLumi(1);

    // histograms
    TFile* file_sig = new TFile(Form("../data/ep_26_04_1/ep_10x100_eID_eff.root"));
    TFile* file_bkg = new TFile(Form("../data/ep_26_04_1/ep_10x100_eID_eff_beamBG.root"));

    TH2F* h_xq2_eff[2];
    TH2F* h_xq2_pur[2];
    TH2F* h_xq2_eID[2];

    for ( int i = 0; i < 2; i ++ )
    {
        TFile* file = i == 0 ? file_sig : file_bkg;
        h_xq2_eff[i] = (TH2F*)file->Get("h_xq2_eff");
        h_xq2_pur[i] = (TH2F*)file->Get("h_xq2_pur");
        h_xq2_eID[i] = (TH2F*)file->Get("h_xq2_eID");
    }

    TH2F* h_eff_diff = BookTH2(Form("h_eff_eff"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    TH2F* h_pur_diff = BookTH2(Form("h_pur_diff"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    TH2F* h_eID_diff = BookTH2(Form("h_eID_diff"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);

    h_eff_diff->SetMinimum(-100);
    h_pur_diff->SetMinimum(-100);
    h_eID_diff->SetMinimum(-100);
    
    for ( int ix = 0; ix < h_eff_diff->GetXaxis()->GetNbins(); ix ++ )
    {
      for ( int iq = 0; iq < h_eff_diff->GetYaxis()->GetNbins(); iq ++ )
      {
        // if ( h_xq2_eff[0]->GetBinContent(ix+1, iq+1) == 0 || h_xq2_eff[1]->GetBinContent(ix+1, iq+1) == 0 )
        //     continue;

        bool x_same = h_xq2_eff[0]->GetXaxis()->GetBinCenter(ix+1) == h_xq2_eff[1]->GetXaxis()->GetBinCenter(ix+1);
        bool q_same = h_xq2_eff[0]->GetYaxis()->GetBinCenter(iq+1) == h_xq2_eff[1]->GetYaxis()->GetBinCenter(iq+1);
        if ( !x_same || !q_same )
        {
          std::cout << "ix: " << ix << ", iq: " << iq << " x_same: " << x_same << ", q_same: " << q_same << std::endl;
          std::cout << "h_xq2_eff[0] x: " << h_xq2_eff[0]->GetXaxis()->GetBinCenter(ix+1) << ", q: " << h_xq2_eff[0]->GetYaxis()->GetBinCenter(iq+1) << std::endl;
          std::cout << "h_xq2_eff[1] x: " << h_xq2_eff[1]->GetXaxis()->GetBinCenter(ix+1) << ", q: " << h_xq2_eff[1]->GetYaxis()->GetBinCenter(iq+1) << std::endl;
          return;
        }

        double eff_diff = (h_xq2_eff[1]->GetBinContent(ix+1, iq+1) - h_xq2_eff[0]->GetBinContent(ix+1, iq+1))*100;
        double pur_diff = (h_xq2_pur[1]->GetBinContent(ix+1, iq+1) - h_xq2_pur[0]->GetBinContent(ix+1, iq+1))*100;
        double eID_diff = (h_xq2_eID[1]->GetBinContent(ix+1, iq+1) - h_xq2_eID[0]->GetBinContent(ix+1, iq+1))*100;

        bool eff_empty = (h_xq2_eff[1]->GetBinContent(ix+1, iq+1) == 0 && h_xq2_eff[0]->GetBinContent(ix+1, iq+1) == 0);
        bool pur_empty = (h_xq2_pur[1]->GetBinContent(ix+1, iq+1) == 0 && h_xq2_pur[0]->GetBinContent(ix+1, iq+1) == 0);
        bool eID_empty = (h_xq2_eID[1]->GetBinContent(ix+1, iq+1) == 0 && h_xq2_eID[0]->GetBinContent(ix+1, iq+1) == 0);

        h_eff_diff->SetBinContent(ix+1, iq+1, eff_empty ? -999 : eff_diff);
        h_pur_diff->SetBinContent(ix+1, iq+1, pur_empty ? -999 : pur_diff);
        h_eID_diff->SetBinContent(ix+1, iq+1, eID_empty ? -999 : eID_diff);

        // double val = h_eID_diff->GetBinContent(ix, iq);
        // if ( val == 0 ) h_eID_diff->SetBinContent(ix, iq, 1e-6);

        // std::cout << "ix: " << ix << ", iq: " << iq << " h_xq2_eff[0]: " << h_xq2_eff[0]->GetBinContent(ix+1, iq+1) << ", h_xq2_eff[1]: " << h_xq2_eff[1]->GetBinContent(ix+1, iq+1) << ", eff_diff: " << eff_diff << std::endl;
      }
    }

    set_2d_scale(h_xq2_eff[0]);

    // mask_zero_bins(h_eff_diff);
    // mask_zero_bins(h_pur_diff);
    // mask_zero_bins(h_eID_diff);

    h_eff_diff->GetZaxis()->SetRangeUser(-20, 2);
    h_pur_diff->GetZaxis()->SetRangeUser(-20, 2);
    h_eID_diff->GetZaxis()->SetRangeUser(-20, 2);

    h_eff_diff->GetZaxis()->SetTitle("C_{sig.+BG}(%) - C_{sig.}(%)");
    h_pur_diff->GetZaxis()->SetTitle("C_{sig.+BG}(%) - C_{sig.}(%)");
    h_eID_diff->GetZaxis()->SetTitle("C_{sig.+BG}(%) - C_{sig.}(%)");

    h_eff_diff->GetZaxis()->SetTitleOffset(1.8);
    h_pur_diff->GetZaxis()->SetTitleOffset(1.8);
    h_eID_diff->GetZaxis()->SetTitleOffset(1.8);

    h_eff_diff->GetZaxis()->SetTitleSize(25);
    h_pur_diff->GetZaxis()->SetTitleSize(25);
    h_eID_diff->GetZaxis()->SetTitleSize(25);

    h_eff_diff->GetZaxis()->SetTitleFont(43);
    h_pur_diff->GetZaxis()->SetTitleFont(43);
    h_eID_diff->GetZaxis()->SetTitleFont(43);

    TCanvas* c_xq2_eff = draw_2d_efficiency(h_eff_diff, "c_xq2_eff", "xq2 eff diff", 1400, 600, false, true, "+.1f");
    draw_manager->LableAndCollect(c_xq2_eff);

    TCanvas* c_xq2_pur = draw_2d_efficiency(h_pur_diff, "c_xq2_pur", "xq2 pur diff", 1400, 600, false, true, "+.1f");
    draw_manager->LableAndCollect(c_xq2_pur);

    TCanvas* c_xq2_eID = draw_2d_efficiency(h_eID_diff, "c_xq2_eID", "xq2 eID diff", 1400, 600, false, true, "+.1f");
    draw_manager->LableAndCollect(c_xq2_eID);

    // for ( int ix = 0; ix < h_eff_diff->GetXaxis()->GetNbins(); ix ++ )
    // {
    //     for ( int iq = 0; iq < h_eff_diff->GetYaxis()->GetNbins(); iq ++ )
    //     {
    //         if ( h_eff_diff->GetBinContent(ix+1, iq+1) != 0 )
    //             std::cout << "ix: " << ix << ", iq: " << iq << ", eff_diff: " << h_eff_diff->GetBinContent(ix+1, iq+1) << std::endl;
    //         if ( h_pur_diff->GetBinContent(ix+1, iq+1) != 0 )
    //             std::cout << "ix: " << ix << ", iq: " << iq << ", pur_diff: " << h_pur_diff->GetBinContent(ix+1, iq+1) << std::endl;
    //         if ( h_eID_diff->GetBinContent(ix+1, iq+1) != 0 )
    //             std::cout << "ix: " << ix << ", iq: " << iq << ", eID_diff: " << h_eID_diff->GetBinContent(ix+1, iq+1) << std::endl;
    //     }
    // }

    TCanvas* c_xq2_eff0 = draw_2d_efficiency(h_xq2_eff[0], "c_xq2_eff0", "xq2 eff 0", 1400, 600, false, true);
    TCanvas* c_xq2_pur0 = draw_2d_efficiency(h_xq2_pur[0], "c_xq2_pur0", "xq2 pur 0", 1400, 600, false, true);
    TCanvas* c_xq2_eID0 = draw_2d_efficiency(h_xq2_eID[0], "c_xq2_eID0", "xq2 eID 0", 1400, 600, false, true);

    TCanvas* c_xq2_eff1 = draw_2d_efficiency(h_xq2_eff[1], "c_xq2_eff1", "xq2 eff 1", 1400, 600, false, true);
    TCanvas* c_xq2_pur1 = draw_2d_efficiency(h_xq2_pur[1], "c_xq2_pur1", "xq2 pur 1", 1400, 600, false, true);
    TCanvas* c_xq2_eID1 = draw_2d_efficiency(h_xq2_eID[1], "c_xq2_eID1", "xq2 eID 1", 1400, 600, false, true);

    c_xq2_eff->Modified();
    c_xq2_eff->Update();
    c_xq2_pur->Modified();
    c_xq2_pur->Update();
    c_xq2_eID->Modified();
    c_xq2_eID->Update();
}