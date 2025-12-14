#include "calculate_a1.h"
#include "bin.h"
#include "BinHelper.h"

#include "../GlobalUtil/Constants.hh"
#include "../GlobalUtil/drawHelper.h"
#include "../GlobalUtil/DrawManager.cc"

// std::string type = "p";
// const double e_lepton = 18;
// const double e_nucleon = 275.;
// std::string setting = "18x275";
// const double m_nucleon = MASS_PROTON;


std::string type = "n";
const double e_lepton = 10;
const double e_nucleon = 166.;
std::string setting = "10x166";
const double m_nucleon = MASS_NEUTRON;



const double ePol = 0.7;
const double iPol = 0.7;


struct asymm_data 
{
    std::vector<double> q2;
    std::vector<double> value;
    std::vector<double> error;
    std::vector<double> sys_error; // point to point sys, uncorrelated uncertainty
    
    std::vector<double> q2_low;
    std::vector<double> value_low;
    TGraphErrors* gr_data;
    TGraphErrors* gr_sys;
    TGraphErrors* gr_low;
    std::vector<TLatex*> text_x;

    std::vector<TGraphErrors*> gr_norm;
    std::vector<double> norm_q2;
    std::vector<double> norm_value;
    std::vector<double> norm_error; // correlated normalization uncertainty
};

void format_data(TGraphErrors* &g)
{
    g->SetMarkerStyle(20);
    g->SetMarkerColor(kBlue);
    g->SetLineColor(kBlue);
    return;
}

void format_low(TGraphErrors* &g)
{
    g->SetMarkerStyle(24);
    g->SetMarkerColor(kBlue);
    g->SetLineColor(kBlue);
    return;
}

void draw_asymm()
{
    std::string type_title = "eHe3";
    std::string energy_title = Form("%dx%d GeV", 10, 166);
    DrawManager* draw_manager = new DrawManager(type_title, energy_title, "25.10.2");
    draw_manager->SetEPIC();
    draw_manager->SetLumi(8.65);

    std::string type = m_nucleon == MASS_NEUTRON ? "n" : "p";

    TH2F* h_xq2 = BookTH2("xq2",  ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);

    std::vector<std::vector<KEbin*>> bins = create_bins(h_xq2, e_lepton, e_nucleon);
    if ( type == "p" )
        load_bins_from_text(Form("../data/en_25_10_2/%s_xq2_collection.txt", setting.c_str()), h_xq2, bins);
    else
        load_bins_from_text_simple(Form("../data/en_25_10_2/%s_xq2_%s_collection.txt", setting.c_str(), type.c_str()), h_xq2, bins);

    asymm_data a1_col;

    double x_last = 0;
    double q_last = 0;
    double a1_last = 0;

    for ( int ix = 0; ix < h_xq2->GetXaxis()->GetNbins(); ix ++ )
    {
        for ( int iq = 0; iq < h_xq2->GetYaxis()->GetNbins(); iq ++ )
        {
            bins[ix][iq]->process_bin(m_nucleon, ePol, iPol);

            if ( bins[ix][iq]->y <= 0.01 || bins[ix][iq]->y >= 0.95)
				continue;

            h_xq2->SetBinContent(ix+1, iq+1, bins[ix][iq]->n_events);

            if ( bins[ix][iq]->n_events == 0 )
                continue;
 
            double x = h_xq2->GetXaxis()->GetBinCenter(ix+1);
			double q = h_xq2->GetYaxis()->GetBinCenter(iq+1);

            double a1 = type == "p" ? find_a1p(x, q) : find_a1n(x, q);
            double a1err = bins[ix][iq]->unc_a1;
            double a1sys = a1 * sqrt(1.5*1.5+3*3) * 0.01;
            double a1norm = a1 * 3.5 * 0.01;

            a1 = a1 - 3 * TMath::Log10(x);
        
            if ( a1err < 0.25 )
            {
                a1_col.q2.push_back(q);
                a1_col.value.push_back(a1);
                a1_col.error.push_back(a1err);
                // a1_col.sys_error.push_back(sqrt(a1sys*a1sys+a1err*a1err));
                a1_col.sys_error.push_back(a1sys);
                // a1_col.norm_error.push_back(a1norm);
            }
            else
            {
                a1_col.q2_low.push_back(q);
                a1_col.value_low.push_back(a1);
            }

            if ( x_last == 0 )
            {
                x_last = x;
            }
            else if ( x_last != x )
            {
                int iq = h_xq2->GetYaxis()->FindBin(q_last);
                double last_bin_width = h_xq2->GetYaxis()->GetBinWidth(iq+1);

                TLatex* ta1 = new TLatex(q_last+last_bin_width, a1_last, Form("x = %f", x_last));
                ta1->SetTextAlign(11);
                ta1->SetTextSize(0.02);
                a1_col.text_x.push_back(ta1);

                a1_col.gr_norm.push_back(new TGraphErrors(a1_col.norm_q2.size(), &a1_col.norm_q2[0], &a1_col.norm_value[0], 0, &a1_col.norm_error[0]));
                a1_col.norm_q2.clear();
                a1_col.norm_value.clear();
                a1_col.norm_error.clear();

                x_last = x;
            }
            q_last = q;
            a1_last = a1;

            a1_col.norm_q2.push_back(q);
            a1_col.norm_value.push_back(a1);
            a1_col.norm_error.push_back(a1norm);
        }
    }
 
    int iq = h_xq2->GetYaxis()->FindBin(q_last);
    double last_bin_width = h_xq2->GetYaxis()->GetBinWidth(iq+1);

    TLatex* ta1 = new TLatex(q_last+last_bin_width, a1_last, Form("x = %f", x_last));
    ta1->SetTextAlign(11);
    ta1->SetTextSize(0.02);
    a1_col.text_x.push_back(ta1);

    a1_col.gr_norm.push_back(new TGraphErrors(a1_col.norm_q2.size(), &a1_col.norm_q2[0], &a1_col.norm_value[0], 0, &a1_col.norm_error[0]));
    a1_col.norm_q2.clear();
    a1_col.norm_value.clear();
    a1_col.norm_error.clear();

    // A1
    TCanvas* c_a1 = new TCanvas("c_a1","", 800, 1000);
	c_a1->SetLogx();
    c_a1->SetLeftMargin(0.13);
    c_a1->SetTickx(1);
    c_a1->SetTicky(1);

    a1_col.gr_data = new TGraphErrors(a1_col.q2.size(), &a1_col.q2[0], &a1_col.value[0], 0, &a1_col.error[0]);
    a1_col.gr_data->Draw("AP SAME");

    format_data(a1_col.gr_data);
    a1_col.gr_data->SetTitle("");
    gStyle->SetTitleFontSize(0.08);
    a1_col.gr_data->GetXaxis()->SetTitle(Form("Q^{2} (GeV/c^{2})^{2}"));
    a1_col.gr_data->GetYaxis()->SetTitle(Form("A_{1}^{%s} - 3 #times log_{10}(x)", type.c_str()));
    format_graph(a1_col.gr_data);
    a1_col.gr_data->SetMinimum(0);
	a1_col.gr_data->SetMaximum(13);
    if ( type == "n" )
        a1_col.gr_data->SetMaximum(11);
	a1_col.gr_data->GetXaxis()->SetLimits(1,2e5);

    a1_col.gr_low = new TGraphErrors(a1_col.q2_low.size(), &a1_col.q2_low[0], &a1_col.value_low[0], 0, 0);
    a1_col.gr_low->Draw("P SAME");
    format_low(a1_col.gr_low);

    a1_col.gr_sys = new TGraphErrors(a1_col.q2.size(), &a1_col.q2[0], &a1_col.value[0], 0, &a1_col.sys_error[0]);
    a1_col.gr_sys->Draw("E SAME");
    format_data(a1_col.gr_sys);
    format_graph(a1_col.gr_sys);
    a1_col.gr_sys->SetMarkerColor(kViolet);
    a1_col.gr_sys->SetLineColor(kViolet);

    for ( auto g : a1_col.gr_norm )
    {
        g->Draw("3 SAME");
        format_data(g);
        format_graph(g);
        // g->SetFillStyle(3003);
        g->SetFillColorAlpha(kAzure, 0.5);
    }

    a1_col.gr_data->Draw("P SAME");
    a1_col.gr_low->Draw("P SAME");
    a1_col.gr_sys->Draw("E SAME");

    c_a1->Update();
    for ( auto t : a1_col.text_x )
        t->Draw("SAME");

    TLegend* leg = new TLegend(0.55, 0.55, 0.85, 0.78);
    // leg->SetHeader(Form("%.0f #times %0.f GeV", e_lepton, e_nucleon), "C");
	leg->AddEntry(a1_col.gr_data, "Statistical unc." , "E");
    leg->AddEntry(a1_col.gr_sys, "Point-by-point sys." , "E");
    leg->AddEntry(a1_col.gr_norm[0], "Normalization sys." , "F");
    leg->AddEntry(a1_col.gr_low, Form("Lower statistics"), "P");
	leg->SetTextSize(0.03);
	leg->SetBorderSize(0);
	leg->SetFillColor(0);
	leg->Draw();

    // SetePICStyle();
    // c_a1->Update();

    // set_2d_scale(h_xq2);
    // TCanvas* c_xq = draw_2d_standard(h_xq2, "c_xq", "all events", 650, 600, true, true);

    draw_manager->LableAndCollectSpecial2(c_a1);
    c_a1->SaveAs(Form("../data/en_25_10_2/asymm_a1_%s_%s.png", type.c_str(), setting.c_str()));
    c_a1->SaveAs(Form("../data/en_25_10_2/asymm_a1_%s_%s.pdf", type.c_str(), setting.c_str()));
}