#include "calculate_a1.h"
#include "bin.h"
#include "BinHelper.h"

#include "../GlobalUtil/Constants.hh"
#include "../GlobalUtil/drawHelper.h"
#include "../GlobalUtil/DrawManager.cc"

// std::string setting = "18x275";
// std::string file_type = "p";

std::string setting = "10x166";
std::string file_type = "n";

void plot_a1_unc()
{
    // SetePICStyle();

    TH2F* h_xq2 = BookTH2("xq2",  ";x_{B};Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);

    TCanvas* c_selected = new TCanvas("c", "a1", 1600, 1200);
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

    vector<TPolyMarker*> m_ex;
    int draw_color[3] = {kRed-7, kAzure+7, kGreen+3};
    double draw_alpha[3] = {0.45, 0.55, 0.65};
    double draw_factor = 5;

    vector<double> x_col;
    vector<double> q_col;
    vector<double> a1_stats_col;
    vector<double> a1_sys_col;
    vector<double> a1_norm_col;
    vector<double> color_col;

    std::ifstream file(Form("../data/en_25_10_2/%s_xq2_%s_asymm_table.txt", setting.c_str(), file_type.c_str()));

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

        x_col.push_back(x);
        q_col.push_back(q2);
        a1_stats_col.push_back(a1err/abs(a1));

        double a1_sys = sqrt(1.5*1.5+3*3) * 0.01;
        double a1_norm = 2.9 * 0.01;
        a1_sys_col.push_back(a1_sys);
        a1_norm_col.push_back(a1_norm);

        color_col.push_back(2);

        std::cout << abs(a1err/a1)*100 << std::endl;
        
        double point_x[1] = {x};
        double point_q[1] = {q2};
        auto m = new TPolyMarker(1, point_x, point_q);
        if ( abs(a1err/a1)*100 > 100 )
        {
            m->SetMarkerStyle(27);
            m->SetMarkerColor(draw_color[0]);
            m->SetMarkerSize(3);
            m_ex.push_back(m);
        }
        else
        {
            m->SetMarkerStyle(20);
            m->SetMarkerColorAlpha(draw_color[0],draw_alpha[0]);
            m->SetMarkerSize(abs(a1err/a1)*draw_factor);
            m->Draw();

            auto m_ref_50 = new TPolyMarker(1, point_x, point_q);
            m_ref_50->SetMarkerStyle(24);
            m_ref_50->SetMarkerColor(kGray+3);
            m_ref_50->SetMarkerSize(50*0.01*draw_factor);
            m_ref_50->Draw();
        }

        // auto m_norm = new TPolyMarker(1, point_x, point_q);
        // m_norm->SetMarkerStyle(20);
        // m_norm->SetMarkerColorAlpha(draw_color[1],draw_alpha[1]);
        // m_norm->SetMarkerSize(a1_norm*draw_factor);
        // m_norm->Draw();

        // auto m_sys = new TPolyMarker(1, point_x, point_q);
        // m_sys->SetMarkerStyle(20);
        // m_sys->SetMarkerColorAlpha(draw_color[2],draw_alpha[2]);
        // m_sys->SetMarkerSize(a1_sys*draw_factor);
        // m_sys->Draw();

        auto m_sys_total = new TPolyMarker(1, point_x, point_q);
        m_sys_total->SetMarkerStyle(20);
        m_sys_total->SetMarkerColorAlpha(draw_color[1],draw_alpha[1]);
        m_sys_total->SetMarkerSize(sqrt(a1_sys*a1_sys+a1_norm*a1_norm)*draw_factor);
        m_sys_total->Draw();
    }


    for ( auto m : m_ex )
        m->Draw();

    TLegend* leg = new TLegend(0.2, 0.52, 0.4, 0.85);
    leg->SetTextSize(0.03);
	leg->SetBorderSize(0);
	leg->SetFillColor(0);
    leg->SetFillStyle(0);

    TGraph* g_leg[6];

    double a1_sys = sqrt(1.5*1.5+3*3) * 0.01;
    double a1_norm = 3.5 * 0.01;
    double point_x[1] = {0};
    double point_q[1] = {0};
    g_leg[4] = new TGraph();
    g_leg[4]->SetMarkerStyle(20);
    g_leg[4]->SetMarkerColorAlpha(draw_color[1],draw_alpha[1]);
    g_leg[4]->SetMarkerSize(1);
    // g_leg[4]->SetMarkerSize(sqrt(a1_sys*a1_sys+a1_norm*a1_norm)*draw_factor);
    // leg->AddEntry(g_leg[4], Form("%.0f%% Systematic", sqrt(a1_sys*a1_sys+a1_norm*a1_norm)*100),"P");
    leg->AddEntry(g_leg[4], Form("Systematic"),"P");

    g_leg[5] = new TGraph();
    g_leg[5]->SetMarkerStyle(20);
    g_leg[5]->SetMarkerColorAlpha(draw_color[0],draw_alpha[0]);
    g_leg[5]->SetMarkerSize(1);
    leg->AddEntry(g_leg[5], Form("Statistical"),"P");

    double point_size[2] = {25, 50};
    for (int i = 0; i < 2; i ++ )
    {
        // double point_x[1] = {0};
        // double point_q[1] = {0};
        g_leg[i] = new TGraph();
        g_leg[i]->SetMarkerStyle(24);
        // g_leg[i]->SetMarkerColorAlpha(kGray+2,0.5);
        g_leg[i]->SetMarkerColorAlpha(kBlack,1);
        g_leg[i]->SetMarkerSize(point_size[i]*0.01*draw_factor);
        leg->AddEntry(g_leg[i], Form("%.0f%%", point_size[i]),"P");
    }

    g_leg[3] = new TGraph();
    g_leg[3]->SetMarkerStyle(27);
    // g_leg[3]->SetMarkerColorAlpha(kGray+2,0.5);
    g_leg[3]->SetMarkerColorAlpha(kBlack,1);
    g_leg[3]->SetMarkerSize(4);
    leg->AddEntry(g_leg[3], Form("> 100%%"),"P");

    if ( file_type == "n" )
        leg->SetHeader("    #DeltaA_{1}^{n} / A_{1}^{n}");
    else
        leg->SetHeader("   #DeltaA_{1}^{p} / A_{1}^{p} : 18x275 GeV");
    leg->Draw();

    // TLatex* text[3];
    // double text_y = 1.25;
    // double text_x[3] = {2e-4, 5e-3, 1.2e-1};
    
    // if ( file_type == "n" )
    // {
    //     text_y = 1.8;
    //     text_x[0] = 3e-2;
    // }
        
    // for ( int i = 0; i < n_file; i ++ )
    // {
    //     text[i] = new TLatex(text_x[i], text_y, Form("%s", setting[i].c_str()));
    //     text[i]->SetTextColorAlpha(draw_color[i], draw_alpha[i]);
    //     text[i]->SetTextSize(0.03);
    //     text[i]->Draw();
    // }

    return;
}