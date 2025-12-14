#include "calculate_a1.h"
#include "bin.h"
#include "BinHelper.h"

#include "../GlobalUtil/Constants.hh"
#include "../GlobalUtil/drawHelper.h"
#include "../GlobalUtil/DrawManager.cc"

// const double e_lepton = 18;
// const double e_nucleon = 110.;
// std::string setting = "18x110";
// double bin_shift = 0.0;

const double e_lepton = 10;
const double e_nucleon = 166.;
std::string setting = "10x166";
double bin_shift = 0.0;

TH2F* h_xq2 = BookTH2("xq2",  ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);

struct asymm_data 
{
    std::vector<double> x;
    std::vector<double> value;
    std::vector<double> error;
    std::vector<double> error_all;
    TGraphErrors* gr_data;
    TGraphErrors* gr_data_all_error;
    // TGraphErrors* gr_err_only;
};

void fill_data(std::string type, asymm_data g1_col[], std::vector<int> qbin)
{
    std::string file_type = type;
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

        for ( int i = 0; i < qbin.size(); i ++ )
        {
            if ( h_xq2->GetYaxis()->FindBin(q2) == qbin[i] )
            {
                int xbin = h_xq2->GetXaxis()->FindBin(x);
                double width = h_xq2->GetXaxis()->GetBinWidth(xbin);
                double x_shift = 0.1*width*i;
                // double x_shift = 0;

                // if ( a1err < 0.25 )
                {
                    g1_col[i].x.push_back(x+x_shift);
                    g1_col[i].value.push_back(g1);
                    // g1_col[i].error.push_back(g1err_all);
                    g1_col[i].error.push_back(g1err);
                }
            }
        }
    }

    return;
}

void format_pad()
{
    gPad->SetTopMargin(0.02);
    gPad->SetBottomMargin(0.2);
	gPad->SetLeftMargin(0.13);
    gPad->SetTickx(1);
    gPad->SetTicky(1);
	gPad->SetLogx();
	gPad->SetGridy();
}

void compare_g1n()
{
    std::string type_title = "eHe3";
    std::string energy_title = Form("%dx%d GeV", 10, 166);
    DrawManager* draw_manager = new DrawManager(type_title, energy_title, "25.10.2");
    draw_manager->SetEPIC();
    draw_manager->SetLumi(8.65);

    std::vector<int> qbin;
    qbin.push_back(h_xq2->GetYaxis()->FindBin(2));
    qbin.push_back(h_xq2->GetYaxis()->FindBin(10));
    // qbin.push_back(h_xq2->GetYaxis()->FindBin(100));
    // qbin.push_back(h_xq2->GetYaxis()->FindBin(900));

    asymm_data g1_n_ff_col[qbin.size()];
    fill_data("n", g1_n_ff_col, qbin);


    // plot

    TLegend* leg_g1;
    leg_g1 = new TLegend(0.65, 0.22, 0.85, 0.41);
    leg_g1->SetTextSize(0.05);
    leg_g1->SetBorderSize(0);
    leg_g1->SetFillColor(0);
    leg_g1->SetFillStyle(0);
    // leg_g1->SetHeader("e^{3}He 10x166 GeV", "");

    TCanvas* c_g1 = new TCanvas("c_g1", "c_g1", 800, 400);
    // c_g1->Divide(1,qbin.size());
    format_pad();

    int color[3] = {kRed, kBlue, kGreen+3};
    for ( int i = 0; i < qbin.size(); i ++ )
    {
        // c_g1->cd(i+1);
        
        g1_n_ff_col[i].gr_data = new TGraphErrors(g1_n_ff_col[i].x.size(), &g1_n_ff_col[i].x[0], &g1_n_ff_col[i].value[0], 0, &g1_n_ff_col[i].error[0]);
        if ( i == 0 )
            g1_n_ff_col[i].gr_data->Draw("ACPE1");
        else
            g1_n_ff_col[i].gr_data->Draw("CPE1 SAME");

        g1_n_ff_col[i].gr_data->SetTitle("");
        format_graph(g1_n_ff_col[i].gr_data);
        g1_n_ff_col[i].gr_data->GetXaxis()->SetTitle("x_{B}");
        g1_n_ff_col[i].gr_data->GetXaxis()->SetTitleOffset(1.1);
        g1_n_ff_col[i].gr_data->GetYaxis()->SetTitle("g_{ 1}^{ n}");
        g1_n_ff_col[i].gr_data->GetYaxis()->SetTitleOffset(1.8);
        g1_n_ff_col[i].gr_data->SetLineColor(color[i]);
        g1_n_ff_col[i].gr_data->SetFillColor(color[i]);
        g1_n_ff_col[i].gr_data->SetFillStyle(3003);

        // leg_g1[i]->SetNColumns(2);
        leg_g1->AddEntry(g1_n_ff_col[i].gr_data, Form("Q^{2} #approx %.0f GeV^{2}", h_xq2->GetYaxis()->GetBinCenter(qbin[i])), "PE");
    }

    leg_g1->Draw();
    
    g1_n_ff_col[0].gr_data->GetXaxis()->SetLimits(1e-4,1);
    // g1_n_ff_col[1].gr_data->GetXaxis()->SetLimits(1e-3,1);
    // g1_n_ff_col[2].gr_data->GetXaxis()->SetLimits(1e-2,1);

    // g1_n_ff_col[0].gr_data->SetMinimum(-0.16);
    // g1_n_ff_col[0].gr_data->SetMaximum(0.09);

    // g1_n_ff_col[1].gr_data->SetMinimum(-0.16);
    // g1_n_ff_col[1].gr_data->SetMaximum(0.04);

    // g1_n_ff_col[2].gr_data->SetMinimum(-0.2);
    // g1_n_ff_col[2].gr_data->SetMaximum(0.31);

    draw_manager->LableAndCollectSpecial(c_g1);
    c_g1->SaveAs(Form("../data/en_25_10_2/compare_g1n_%s_%s.png", setting.c_str(), "n"));
}