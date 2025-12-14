#include "helper.h"

const double ePol = 0.7;
const double iPol = 0.7;
const double m_nucleon = MASS_PROTON;
// const double m_nucleon = MASS_NEUTRON;

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

void calc_data()
{
    std::vector<std::vector<KEbin*>> bins = create_bins(h_xq2, e_lepton, e_nucleon);
    load_function(h_xq2, bins, "p");
    load_function(h_xq2, bins, "n");
    load_function(h_xq2, bins, "he3");
}

void fill_data(std::string type, asymm_data a1_col[], std::vector<int> qbin)
{
    std::string file_type = type;
    if ( type == "p" )
        file_type = "n";

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
                // double x_shift = 0.1*width*i;
                double x_shift = 0;

                if ( type == "he3" )
                    x_shift = 0.1*width*1;
                
                if ( type == "nx" )
                    x_shift = 0.05*width*1;

                if ( type == "p" )
                    a1 = find_a1p(x, q2);

                // if ( type == "n" )
                //     a1 = find_a1n(x, q2);

                // if ( a1err < 0.25 )
                {
                    a1_col[i].x.push_back(x+x_shift);
                    a1_col[i].value.push_back(a1);
                    a1_col[i].error.push_back(a1err);
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

void compare_unc()
{
    std::string type_title = "eHe3";
    std::string energy_title = Form("%dx%d GeV", 10, 166);
    DrawManager* draw_manager = new DrawManager(type_title, energy_title, "25.10.2");
    draw_manager->SetEPIC();
    draw_manager->SetLumi(8.65);

    std::vector<int> qbin;
    qbin.push_back(h_xq2->GetYaxis()->FindBin(2));
    qbin.push_back(h_xq2->GetYaxis()->FindBin(10));
    qbin.push_back(h_xq2->GetYaxis()->FindBin(100));
    // qbin.push_back(h_xq2->GetYaxis()->FindBin(900));

    asymm_data a1_n_ff_col[qbin.size()];
    asymm_data a1_n_ex_col[qbin.size()];
    asymm_data a1_he3_col[qbin.size()];
    asymm_data a1_p_col[qbin.size()];

    fill_data("he3", a1_he3_col, qbin);
    fill_data("n", a1_n_ff_col, qbin);
    fill_data("p", a1_p_col, qbin);
    fill_data("nx", a1_n_ex_col, qbin);

    for ( int i = 0; i < qbin.size(); i ++ )
    {
        for ( int j = 0; j < a1_n_ff_col[i].error.size(); j ++ )
            printf("unc diff: %f %\n", (a1_n_ex_col[i].error[j] - a1_n_ff_col[i].error[j])/a1_n_ex_col[i].error[j] * 100 );
    }

    // plot

    TLegend* leg_a1[qbin.size()];

    TCanvas* c_a1 = new TCanvas("c_a1", "c_a1", 800, 300*qbin.size()+100);
    c_a1->Divide(1,qbin.size());

    for ( int i = 0; i < qbin.size(); i ++ )
    {
        c_a1->cd(i+1);
        format_pad();

        a1_n_ex_col[i].gr_data = new TGraphErrors(a1_n_ex_col[i].x.size(), &a1_n_ex_col[i].x[0], &a1_n_ex_col[i].value[0], 0, &a1_n_ex_col[i].error[0]);
        a1_n_ex_col[i].gr_data->Draw("ACPE1");

        a1_n_ex_col[i].gr_data->SetTitle("");
        format_graph(a1_n_ex_col[i].gr_data);
        a1_n_ex_col[i].gr_data->GetXaxis()->SetTitle("x_{B}");
        a1_n_ex_col[i].gr_data->GetXaxis()->SetTitleOffset(1.1);
        a1_n_ex_col[i].gr_data->GetYaxis()->SetTitle("A_{1}");
        a1_n_ex_col[i].gr_data->GetYaxis()->SetTitleOffset(1.8);
        a1_n_ex_col[i].gr_data->SetLineColor(kRed);
        a1_n_ex_col[i].gr_data->SetFillColor(kRed);
        a1_n_ex_col[i].gr_data->SetFillStyle(3003);

        a1_n_ff_col[i].gr_data = new TGraphErrors(a1_n_ff_col[i].x.size(), &a1_n_ff_col[i].x[0], &a1_n_ff_col[i].value[0], 0, &a1_n_ff_col[i].error[0]);
        a1_n_ff_col[i].gr_data->Draw("CPE1 SAME");
        a1_n_ff_col[i].gr_data->SetLineColor(kBlue);
        a1_n_ff_col[i].gr_data->SetFillColor(kBlue);
        a1_n_ff_col[i].gr_data->SetFillStyle(3003);

        a1_he3_col[i].gr_data = new TGraphErrors(a1_he3_col[i].x.size(), &a1_he3_col[i].x[0], &a1_he3_col[i].value[0], 0, &a1_he3_col[i].error[0]);
        a1_he3_col[i].gr_data->Draw("CPE1 SAME");
        a1_he3_col[i].gr_data->SetLineColor(kOrange+1);
        a1_he3_col[i].gr_data->SetFillColor(kOrange+1);
        a1_he3_col[i].gr_data->SetFillStyle(3003);

        a1_p_col[i].gr_data = new TGraphErrors(a1_p_col[i].x.size(), &a1_p_col[i].x[0], &a1_p_col[i].value[0], 0, 0);
        a1_p_col[i].gr_data->Draw("C SAME");
        a1_p_col[i].gr_data->SetLineColor(kBlack);
        a1_p_col[i].gr_data->SetLineStyle(10);
        a1_p_col[i].gr_data->SetLineWidth(2);
        
        leg_a1[i] = new TLegend(0.13, 0.25, 0.85, 0.35);
        leg_a1[i]->SetTextSize(0.06);
        leg_a1[i]->SetBorderSize(0);
        leg_a1[i]->SetFillColor(0);
        leg_a1[i]->SetFillStyle(0);

        // double q_bin_min = h_xq2->GetYaxis()->GetBinLowEdge(qbin[i]);
        // double q_bin_max = q_bin_min + h_xq2->GetYaxis()->GetBinWidth(qbin[i]);
        // leg_a1[i]->SetHeader(Form("Q^{2} = [%.2f, %.2f) GeV^{2}", q_bin_min, q_bin_max), "C");
        // leg_a1[i]->SetHeader(Form("Q^{2} = %.0f GeV^{2}", h_xq2->GetYaxis()->GetBinCenter(qbin[i])));
        leg_a1[i]->SetNColumns(5);
        leg_a1[i]->AddEntry((TObject*)0, Form("Q^{2} #approx %.0f GeV^{2}", h_xq2->GetYaxis()->GetBinCenter(qbin[i])),"");
        leg_a1[i]->AddEntry(a1_n_ff_col[i].gr_data, "n, double tagged", "LP");
        leg_a1[i]->AddEntry(a1_n_ex_col[i].gr_data, "n, traditional", "LP");
        leg_a1[i]->AddEntry(a1_he3_col[i].gr_data, "^{3}He", "LP");
        leg_a1[i]->AddEntry(a1_p_col[i].gr_data, "p", "LP");
        leg_a1[i]->Draw();
    }
    
    // a1_n_ex_col[0].gr_data->GetXaxis()->SetLimits(1e-4,1e-1);
    // a1_n_ex_col[1].gr_data->GetXaxis()->SetLimits(1e-3,1);
    // a1_n_ex_col[2].gr_data->GetXaxis()->SetLimits(1e-2,1);

    a1_n_ex_col[0].gr_data->SetMinimum(-0.16);
    a1_n_ex_col[0].gr_data->SetMaximum(0.09);

    a1_n_ex_col[1].gr_data->SetMinimum(-0.16);
    a1_n_ex_col[1].gr_data->SetMaximum(0.04);

    a1_n_ex_col[2].gr_data->SetMinimum(-0.5);
    a1_n_ex_col[2].gr_data->SetMaximum(0.5);

    draw_manager->LableAndCollect(c_a1);
    c_a1->SaveAs(Form("../data/en_25_10_2/compare_a1_unc_%s_%dx%d.png", setting.c_str(), (int)e_lepton, (int)e_nucleon));
}