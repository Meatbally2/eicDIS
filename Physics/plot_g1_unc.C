#include "calculate_a1.h"
#include "bin.h"
#include "BinHelper.h"

#include "../GlobalUtil/Constants.hh"
#include "../GlobalUtil/drawHelper.h"
#include "../GlobalUtil/DrawManager.cc"

#include "TPolyMarker.h"

// const int n_file = 3;
// std::string setting[3] = {"18x275", "10x100", "5x41"};
// std::string file_type = "p";

const int n_file = 1;
std::string setting[n_file] = {"10x166"};
std::string file_type = "n";

void plot_g1_unc()
{
    std::string type_title = "eHe3";
    std::string energy_title = Form("%dx%d GeV", 10, 166);
    DrawManager* draw_manager = new DrawManager(type_title, energy_title, "25.10.2");
    draw_manager->SetEPIC();
    draw_manager->SetLumi(8.65);

    TH2F* h_xq2 = BookTH2("xq2",  ";x_{B};Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);

    TCanvas* c_selected = new TCanvas("c", "g1", 1600, 1200);
    c_selected->SetLogx();
    c_selected->SetLogy();
    gPad->SetTickx(1);
    gPad->SetTicky(1);
    gPad->SetLeftMargin(0.12);
    gPad->SetBottomMargin(0.12);
    gPad->SetTopMargin(0.043);
    gPad->SetRightMargin(0.03);
    h_xq2->Draw();
    h_xq2->SetStats(0);

    if ( file_type == "n" )
    {
        h_xq2->GetXaxis()->SetRangeUser(4e-4,1);
        h_xq2->GetYaxis()->SetRangeUser(1,1e4);
    }
    else
    {
        h_xq2->GetXaxis()->SetRangeUser(1e-4,1);
        h_xq2->GetYaxis()->SetRangeUser(1,1e4);
    }
    h_xq2->Draw();

    // TScatter* g_g1[3];
    vector<TPolyMarker*> m_ex;
    int draw_color[3] = {kRed-7, kAzure+7, kGreen+3};
    double draw_alpha[3] = {0.45, 0.55, 0.65};
    for ( int i = 0; i < n_file; i ++ )
    {
        // g_g1[i] = new TGraph();
        vector<double> x_col;
        vector<double> q_col;
        // vector<double> g1_col;
        vector<double> g1_err_col;
        vector<double> color_col;

        std::ifstream file(Form("../data/en_25_10_2/%s_xq2_%s_asymm_table.txt", setting[i].c_str(), file_type.c_str()));

        std::string line;
        while ( getline(file, line) )
        {
            double x = 0;
            double q2 = 0;
            double a1 = 0;
            double a1err = 0;
            double g1 = 0;
            double g1err = 0;
            double g1err_all = 0;

            std::stringstream ss(line);
            ss >> x >> q2 >> a1 >> a1err >> g1 >> g1err >> g1err_all;

            // cout << x << " " << q << " " << a1 << " " << a1err << " " << g1 << " " << g1err << " " << g1err_all << " " << bins[ix][iq]->f1 << " " <<  bins[ix][iq]->f1err << endl;
            cout << setting[i] << " " << x << " " << q2 << " " << g1 << " " << g1err << " " << g1err_all << endl;

            x_col.push_back(x);
            q_col.push_back(q2);
            g1_err_col.push_back(g1err/g1);
            color_col.push_back(2);

            double point_x[1] = {x};
            double point_q[1] = {q2};
            auto m = new TPolyMarker(1, point_x, point_q);
            if ( abs(g1err/g1)*100 > 60 )
            {
                m->SetMarkerStyle(27);
                m->SetMarkerColor(draw_color[i]);
                m->SetMarkerSize(3+i);
                m_ex.push_back(m);
            }
            else
            {
                m->SetMarkerStyle(20);
                m->SetMarkerColorAlpha(draw_color[i],draw_alpha[i]);
                m->SetMarkerSize(abs(g1err/g1)*10);
                m->Draw();

                auto m_ref_50 = new TPolyMarker(1, point_x, point_q);
                m_ref_50->SetMarkerStyle(24);
                m_ref_50->SetMarkerColor(kGray+3);
                m_ref_50->SetMarkerSize(5);
                m_ref_50->Draw();

                auto m_ref_30 = new TPolyMarker(1, point_x, point_q);
                m_ref_30->SetMarkerStyle(24);
                m_ref_30->SetMarkerColor(kGray+3);
                m_ref_30->SetMarkerSize(2.5);
                m_ref_30->Draw();

                // auto m_ref_10 = new TPolyMarker(1, point_x, point_q);
                // m_ref_10->SetMarkerStyle(24);
                // m_ref_10->SetMarkerColor(kGray+3);
                // m_ref_10->SetMarkerSize(1);
                // m_ref_10->Draw();
            }

            // std::cout << abs(g1err/g1) << endl;
        }
    }

    for ( auto m : m_ex )
        m->Draw();

    TLegend* leg = new TLegend(0.17, 0.43, 0.38, 0.50);
    leg->SetTextSize(0.02);
	leg->SetBorderSize(0);
	leg->SetFillColor(0);
    leg->SetFillStyle(0);

    TGraph* g_leg[4];
    // double point_size[3] = {10, 30, 50};
    // for (int i = 0; i < 3; i ++ )
    // {
    //     double point_x[1] = {0};
    //     double point_q[1] = {0};
    //     g_leg[i] = new TGraph();
    //     g_leg[i]->SetMarkerStyle(24);
    //     // g_leg[i]->SetMarkerColorAlpha(kGray+2,0.5);
    //     g_leg[i]->SetMarkerColorAlpha(kBlack,1);
    //     g_leg[i]->SetMarkerSize(point_size[i]*0.1);
    //     leg->AddEntry(g_leg[i], Form("%.0f%%", point_size[i]),"P");
    // }
    // g_leg[2]->SetMarkerStyle(24);

    double point_size[2] = {25, 50};
    for (int i = 0; i < 2; i ++ )
    {
        double point_x[1] = {0};
        double point_q[1] = {0};
        g_leg[i] = new TGraph();
        g_leg[i]->SetMarkerStyle(24);
        // g_leg[i]->SetMarkerColorAlpha(kGray+2,0.5);
        g_leg[i]->SetMarkerColorAlpha(kBlack,1);
        g_leg[i]->SetMarkerSize(point_size[i]*0.1);
        leg->AddEntry(g_leg[i], Form("%.0f%%", point_size[i]),"P");
    }

    g_leg[3] = new TGraph();
    g_leg[3]->SetMarkerStyle(27);
    // g_leg[3]->SetMarkerColorAlpha(kGray+2,0.5);
    g_leg[3]->SetMarkerColorAlpha(kBlack,1);
    g_leg[3]->SetMarkerSize(4);
    leg->AddEntry(g_leg[3], Form("> 60%%"),"P");

    if ( file_type == "n" )
        leg->SetHeader("     #Deltag_{1}^{n} / g_{1}^{n}");
    else
        leg->SetHeader("     #Deltag_{1}^{p} / g_{1}^{p}");
    leg->Draw();

    TLatex* text[3];
    double text_y = 1.25;
    double text_x[3] = {2e-4, 5e-3, 1.2e-1};
    
    if ( file_type == "n" )
    {
        text_y = 1.8;
        text_x[0] = 3e-2;
    }
        
    for ( int i = 0; i < n_file; i ++ )
    {
        text[i] = new TLatex(text_x[i], text_y, Form("%s", setting[i].c_str()));
        text[i]->SetTextColorAlpha(draw_color[i], draw_alpha[i]);
        text[i]->SetTextSize(0.03);
        text[i]->Draw();
    }

    draw_manager->LableAndCollect(c_selected);
    c_selected->SaveAs(Form("../data/en_25_10_2/g1_unc_%s.png", file_type.c_str()));
    return;
}