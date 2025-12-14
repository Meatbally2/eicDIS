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

struct asymm_data 
{
    std::vector<double> q2;
    std::vector<double> value;
    std::vector<double> error;
    std::vector<double> q2_low;
    std::vector<double> value_low;
    TGraphErrors* gr_data;
    TGraphErrors* gr_low;
    std::vector<TLatex*> text_x;
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

double find_a1_he3(double x, double q2, double f2p, double f2n, double f2_he3)
{
    double a1p = find_a1p(x, q2);
    double pp = -0.028;

    double a1n = find_a1n(x, q2);
    double pn = 0.86;

    // double a1_he3 = a1n * pn*f2n*(1+0.056/pn) / f2_he3 + 2*(f2p/f2_he3)*pp*a1p*(1-0.014/(2*pp));
    double a1_he3 = pn*(f2n/f2_he3)*a1n + 2*pp*(f2p/f2_he3)*a1p; // simplified model

    // cout << x << " " << q2 << " " << a1p << " " << a1n << " " << a1_he3 << endl; 

    return a1_he3;
}



void print_he3_a1()
{
	TH2F* h_xq2_bin = BookTH2("xq2",  ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);

    TFile* file = new TFile(Form("../data/en_25_10_2/FullRecon_%dx%d.root", (int)e_lepton, (int)e_nucleon));
	TH2F* h_xq2 = (TH2F*) file->Get("H_XQ_HE3");
    // TH2F* h_xq2 = (TH2F*) file->Get("h_xq2_mc_all");

	if ( h_xq2 == NULL )
	{
		std::cout << "histogram not found" << std::endl;
		return;
	}

    // make plots
	set_2d_scale(h_xq2);
	TCanvas* c_xq = draw_2d_standard(h_xq2, "c_xq", "n events", 650, 600, true, true);

	std::vector<std::vector<KEbin*>> bins = create_bins(h_xq2, e_lepton, e_nucleon);

	for ( int ix = 0; ix < h_xq2->GetXaxis()->GetNbins(); ix ++ )
    {
        for ( int iq = 0; iq < h_xq2->GetYaxis()->GetNbins(); iq ++ )
        {
			double x  = h_xq2->GetXaxis()->GetBinCenter(ix+1);
            double q2 = h_xq2->GetYaxis()->GetBinCenter(iq+1);
            double n = h_xq2->GetBinContent(ix+1, iq+1);

            if ( n!= 0 )
            {
                // cout << x << endl;
                bins[ix][iq]->set_bin(x, q2, n, n);
            }    
            	
        }
    }

    load_function(h_xq2_bin, bins, "p");
    load_function(h_xq2_bin, bins, "n");
    load_function(h_xq2_bin, bins, "he3");

    double Q2max = 4*e_nucleon*e_lepton;
	double nu_max = Q2max/(2*m_nucleon);

    asymm_data a1_col;
    asymm_data g1_col;

    double x_last = 0;
    double q_last = 0;
    double a1_last = 0;

    ofstream outfile(Form("../data/en_25_10_2/%s_xq2_he3_asymm_table.txt", setting.c_str()));

    for ( int ix = 0; ix < h_xq2->GetXaxis()->GetNbins(); ix ++ )
    {
        for ( int iq = 0; iq < h_xq2->GetYaxis()->GetNbins(); iq ++ )
        {
            if ( !bins[ix][iq]->is_averaged )
                continue;

            bins[ix][iq]->process_bin(m_nucleon, ePol, iPol);

            if ( bins[ix][iq]->y <= 0.01 || bins[ix][iq]->y >= 0.95)
                continue;

            double x = h_xq2->GetXaxis()->GetBinCenter(ix+1);
            double q2 = h_xq2->GetYaxis()->GetBinCenter(iq+1);

            double func_a1 = find_a1_he3(x, q2, bins[ix][iq]->f2p, bins[ix][iq]->f2n, bins[ix][iq]->f2_he3);
            double a1 = func_a1 - 3 * TMath::Log10(x);
            double err = bins[ix][iq]->unc_a1;

            // cout << x << " " << q2 << " " << bins[ix][iq]->f2p << " " << bins[ix][iq]->f2n << " " << bins[ix][iq]->f2_he3 << " " << func_a1 << endl;

            outfile << x << " " << q2 << " " << " " << func_a1 << " " << err << " " << 0. << " " << 0. << " " << 0. << std::endl;

            // cout << x << endl;

            if ( err < 0.25 )
            {
                a1_col.q2.push_back(q2);
                a1_col.value.push_back(a1);
                a1_col.error.push_back(err);
            }
            else
            {
                a1_col.q2_low.push_back(q2);
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

                x_last = x;
            }
            q_last = q2;
            a1_last = a1;

        }
    }
    
    outfile.close();

    int iq = h_xq2->GetYaxis()->FindBin(q_last);
    double last_bin_width = h_xq2->GetYaxis()->GetBinWidth(iq+1);

    TLatex* ta1 = new TLatex(q_last+last_bin_width, a1_last, Form("x = %f", x_last));
    ta1->SetTextAlign(11);
    ta1->SetTextSize(0.02);
    a1_col.text_x.push_back(ta1);

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
    a1_col.gr_data->GetYaxis()->SetTitle(Form("A_{1}^{^{3}He} - 3 #times log_{10}(x)"));
    format_graph(a1_col.gr_data);
    a1_col.gr_data->SetMinimum(0);
	a1_col.gr_data->SetMaximum(11);
	a1_col.gr_data->GetXaxis()->SetLimits(1,5e4);

    a1_col.gr_low = new TGraphErrors(a1_col.q2_low.size(), &a1_col.q2_low[0], &a1_col.value_low[0], 0, 0);
    a1_col.gr_low->Draw("P SAME");
    format_low(a1_col.gr_low);

    c_a1->Update();
    for ( auto t : a1_col.text_x )
        t->Draw("SAME");

    TLegend* leg = new TLegend(0.57, 0.73, 0.85, 0.83);
	leg->AddEntry(a1_col.gr_data, Form("%.0f #times %0.f GeV", e_lepton, e_nucleon), "PE");
    leg->AddEntry(a1_col.gr_low, Form("Lower statistics"), "P");
	leg->SetTextSize(0.03);
	leg->SetBorderSize(0);
	leg->SetFillColor(0);
	leg->Draw();

    return;
}