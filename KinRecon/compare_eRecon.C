#include "../GlobalUtil/DrawManager.cc"
#include "../GlobalUtil/Constants.hh"
#include "../GlobalUtil/drawHelper.h"
// #include "GlobalUtil/bin.h"

const std::string group[4] = {"minQ2=1", "minQ2=10", "minQ2=100", "minQ2=1000"};

void compare_eRecon(int beam_type, int Ee, int Eh)
{
    // ePIC plotting style setup
    std::string type_title[5] = {"e^{3}He", "ep", "#gammap", "ep w. BeamBG", "ep"};
    std::string energy_title = beam_type ? Form("%dx%d GeV", Ee, Eh) : Form("%dx%d GeV/A", Ee, Eh);
    std::string campaign = beam_type ? "25.10.0" : "25.10.2";
    if ( beam_type == 4 )
        campaign = "25.10.4";
    DrawManager* draw_manager = new DrawManager(type_title[beam_type], energy_title, campaign);
    draw_manager->SetEPIC();

    double text_lumi = 10; // fb^-1
    if ( beam_type == 0 && Ee == 10 && Eh == 166 )
        text_lumi = 8.65; // fb^-1
    if ( beam_type == 4 && Ee == 10 && Eh == 130 )
        text_lumi = 5.33; // fb^-1
    if ( beam_type == 4 && Ee == 10 && Eh == 250 )
        text_lumi = 9.18; // fb^-1
    draw_manager->SetLumi(text_lumi);

    std::string setting = beam_type == 0 ? Form("eHe3_%dx%d", (int)Ee, (int)Eh) : Form("ep_%dx%d", (int)Ee, (int)Eh);
    TFile* f = new TFile(Form("%s_ReconStats.root", setting.c_str()));
    if( f == NULL )
    {
        printf("File not found!\n");
        return;
    }

    TH2F* h_xq2= BookTH2(Form("h_xq2"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    TH2F* h_xq2_eff = BookTH2(Form("h_xq2_eff"), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    h_xq2->SetStats(0);
    h_xq2_eff->SetStats(0);

    TCanvas* c_recon[5];
    TH2F* h_recon[5];
    string recon_name[5] = {"E", "JB", "DA", "Sig", "ESig"};

    for ( int i = 0; i < 5; i ++ )
    {
        c_recon[i] = (TCanvas*)f->Get(Form("c_%s_Method_eff", recon_name[i].c_str()));
        if ( c_recon[i] == NULL )
            printf("%s CANVAS not found\n", recon_name[i].c_str());

        h_recon[i] = (TH2F*)c_recon[i]->GetPrimitive(Form("H_%s_Method_Efficiency_all", recon_name[i].c_str()));
        if ( h_recon[i] == NULL )
            printf("%s HIST not found\n", recon_name[i].c_str());
    }

    TGraph* g_recon[5];
    for ( int i = 0; i < 5; i ++ )
        g_recon[i] = new TGraph();

    ofstream outfile(Form("%s_eID_selection_test.txt", setting.c_str()));
    
    for ( int ix = 0; ix < h_xq2_eff->GetXaxis()->GetNbins(); ix ++ )
    {
        for ( int iq = 0; iq < h_xq2_eff->GetYaxis()->GetNbins(); iq ++ )
        {
            double max = 0;
            int select = -1;
            for ( int i = 0; i < 5; i ++ )
            {
                if ( h_recon[i]->GetBinContent(ix+1,iq+1) > max )
                {
                    max = h_recon[i]->GetBinContent(ix+1,iq+1);
                    select = i;
                }
            }

            if ( select != -1 )
            {
                double x = h_xq2_eff->GetXaxis()->GetBinCenter(ix+1);
                double q = h_xq2_eff->GetYaxis()->GetBinCenter(iq+1);
                h_xq2_eff->SetBinContent(ix+1,iq+1,max);
                g_recon[select]->AddPoint(x, q);
                outfile << x << " " << q << " " << select << std::endl;
            }
        }
    }

    outfile.close();

    set_2d_scale(h_xq2_eff);
    TCanvas* c_xq_eff = draw_2d_efficiency(h_xq2_eff, "c_xq_eff", "combined eff", 1400, 600, false, true);

    TCanvas* c_selected = draw_2d_standard(h_xq2, "c_selected", "all events", 1400, 600, false, true);
    //  = new TCanvas("c_selected", "c_selected", 1400, 600);
    // c_selected->SetLogx();
    // c_selected->SetLogy();
    h_xq2->Draw();

    TLegend* leg = new TLegend(0.18, 0.46, 0.55, 0.66);
    leg->SetTextSize(0.042);
	leg->SetBorderSize(0);
	leg->SetFillColor(0);
    leg->SetFillStyle(0);

    int color[5] = {2, 4, 1, kOrange+1, kGreen+3};
    string long_name[5] = {"Electron", "Jacquet-Blondel", "Double-Angle", "Sigma", "E-Sigma"};

    for ( int i = 0; i < 5; i ++ )
    {
        if ( g_recon[i]->GetN() < 1 )
            continue;

        g_recon[i]->SetMarkerStyle(20);
        g_recon[i]->SetMarkerSize(2);
        g_recon[i]->SetMarkerColor(color[i]);

        g_recon[i]->Draw("P SAME");
        leg->AddEntry(g_recon[i], Form("%s", long_name[i].c_str()), "P");
    }

    leg->Draw();

    const int n_line = 2;
    double s = 4*16*166;
    double y_line[n_line] = { 0.045, 0.2};
	double x_range[2] = {x_gen_min, 1};

    double q_range[2];
    TGraph* g_yline[n_line+1];
    for ( int i = 0; i < n_line; i ++ )
    {
        q_range[0] = (s - MASS_NEUTRON*MASS_NEUTRON)*x_range[0]*y_line[i];
        q_range[1] = (s - MASS_NEUTRON*MASS_NEUTRON)*x_range[1]*y_line[i];

        g_yline[i] = new TGraph(2, x_range, q_range);
        g_yline[i]->SetLineColor(kBlack);
	    g_yline[i]->SetLineWidth(2);
        g_yline[i]->SetLineStyle(9);

        if ( i == 0)
	        g_yline[i]->Draw("SAME L");
    }

    // Add text labels for y lines
    TLatex* latex = new TLatex();
    latex->SetTextSize(0.035);
    latex->SetTextAlign(12); // left-aligned

    for ( int i = 0; i < n_line; i ++ )
    {
        double label_x = 1.05; // x position for label (adjust as needed)
        double label_q = (s - MASS_NEUTRON*MASS_NEUTRON)*label_x*y_line[i];
        if ( i == 0)
            latex->DrawLatex(label_x, label_q, Form("y = %.3f", y_line[i]));
    }

    c_selected->Modified();
    c_selected->Update();

    q_range[0] = 1;
    // q_range[1] = q_gen_max;
    q_range[1] = pow(10,gPad->GetUymax());
    x_range[0] = 0.025;
    x_range[1] = 0.025;
    g_yline[n_line] = new TGraph(2, x_range, q_range);
    g_yline[n_line]->SetLineColor(kBlack);
    g_yline[n_line]->SetLineWidth(2);
    // g_yline[n_line]->Draw("SAME L");

    draw_manager->LableAndCollect(c_selected);
    // c_selected->SaveAs(Form("../data/en_25_10_2/eID_selection_%dx%d.png", Ee, Eh));
}