
#include "calculate_a1.h"
#include "bin.h"
#include "BinHelper.h"

#include "../GlobalUtil/Constants.hh"
#include "../GlobalUtil/drawHelper.h"
#include "../GlobalUtil/DrawManager.cc"

const double ePol = 0.7;
const double iPol = 0.7;
// const double m_nucleon = 938.272088*1e-3;
const double m_nucleon = MASS_NEUTRON;

// const double e_lepton = 18;
// const double e_nucleon = 110.;
// std::string setting = "18x110";
// double bin_shift = 0.0;

const double e_lepton = 10;
const double e_nucleon = 166.;
std::string setting = "10x166";
double bin_shift = 0.0;

void asymm()
{
	const int n_x_bin = 25;
    const int n_q_bin = 25;

	TFile* file = new TFile(Form("../data/en_25_10_2/FullRecon_%dx%d.root", (int)e_lepton, (int)e_nucleon));
	TH2F* h_xq2 = (TH2F*) file->Get("H_XQ_REC");
	// TH2F* h_xq2 = (TH2F*) file->Get("h_xq2_mc_n");

	if ( h_xq2 == NULL )
	{
		std::cout << "histogram not found" << std::endl;
		return;
	}

	std::vector<std::vector<KEbin*>> bins = create_bins(h_xq2, e_lepton, e_nucleon);

	for ( int ix = 0; ix < h_xq2->GetXaxis()->GetNbins(); ix ++ )
    {
        for ( int iq = 0; iq < h_xq2->GetYaxis()->GetNbins(); iq ++ )
        {
			double x  = h_xq2->GetXaxis()->GetBinCenter(ix+1);
            double q2 = h_xq2->GetYaxis()->GetBinCenter(iq+1);
            double n = h_xq2->GetBinContent(ix+1, iq+1);

            if ( n!= 0 )    
            	bins[ix][iq]->set_bin(x, q2, n, n);
        }
    }

    // make plots
	set_2d_scale(h_xq2);
	TCanvas* c_xq = draw_2d_standard(h_xq2, "c_xq", "n events", 650, 600, true, true);

	double Q2max = 4*e_nucleon*e_lepton;
	double nu_max = Q2max/(2*m_nucleon);

	std::vector<TLatex*> t_a1_x;
	std::vector<double> x_col;
	std::vector<std::vector<double>> q2_col;
	std::vector<std::vector<double>> a1_col;
	std::vector<std::vector<double>> er_col;
	std::vector<std::vector<double>> q2_col_low_acp;
	std::vector<std::vector<double>> a1_col_low_acp;
	std::vector<std::vector<double>> er_col_low_acp;

	ofstream outfile(Form("../data/en_25_10_2/%s_xq2_n_collection.txt", setting.c_str()));

	for ( int ix = 0; ix < h_xq2->GetXaxis()->GetNbins(); ix ++ )
	{
		std::vector<double> q2_set;
		std::vector<double> a1_set;
		std::vector<double> er_set;
		std::vector<double> q2_set_low_acp;
		std::vector<double> a1_set_low_acp;
		std::vector<double> er_set_low_acp;

		double x_first = 0;
		double last_bin_width = 0;
		for ( int iq = 0; iq < h_xq2->GetYaxis()->GetNbins(); iq ++ )
		{
			if ( !bins[ix][iq]->is_averaged )
				continue;

			// cout << h_xq2->GetXaxis()->GetBinCenter(ix+1) << " " << h_xq2->GetYaxis()->GetBinCenter(iq+1) << endl;
			bins[ix][iq]->process_bin(m_nucleon, ePol, iPol);
			// cout <<  " end bin " << endl;

			if ( bins[ix][iq]->y <= 0.01 || bins[ix][iq]->y >= 0.95)
				continue;

			double x = h_xq2->GetXaxis()->GetBinCenter(ix+1);
			double q = h_xq2->GetYaxis()->GetBinCenter(iq+1);
			double print_a1 = find_a1n(x, q);
			double a1 = print_a1 - 3 * TMath::Log10(x);
			double err = bins[ix][iq]->unc_a1;

			last_bin_width = h_xq2->GetYaxis()->GetBinWidth(iq+1);

			if ( err < 0.3 )
			{
				q2_set.push_back(q+bin_shift*last_bin_width);
				a1_set.push_back(a1);
				er_set.push_back(err);
			}
			else
			{
				q2_set_low_acp.push_back(q+bin_shift*last_bin_width);
				a1_set_low_acp.push_back(a1);
				er_set_low_acp.push_back(err);
			}
			
			outfile << x << " " << q << " " << bins[ix][iq]->n_events << " " << print_a1 << " " << err << std::endl;
			
			if ( x_first == 0)
				x_first = x;
		}

		double text_x = 0;
		double text_q = 0;
		if ( !a1_set.empty() )
		{
			q2_col.push_back(q2_set);
			a1_col.push_back(a1_set);
			er_col.push_back(er_set);
			text_x = a1_set.back();
			text_q = q2_set.back()+last_bin_width;
		}

		double text_x_low_acp = 0;
		double text_q_low_acp = 0;
		if ( !a1_set_low_acp.empty() )
		{
			q2_col_low_acp.push_back(q2_set_low_acp);
			a1_col_low_acp.push_back(a1_set_low_acp);
			er_col_low_acp.push_back(er_set_low_acp);
			text_x_low_acp = a1_set_low_acp.back();
			text_q_low_acp = q2_set_low_acp.back()+last_bin_width;
		}

		double print_x = text_q > text_q_low_acp ? text_x : text_x_low_acp;
		double print_q = text_q > text_q_low_acp ? text_q : text_q_low_acp;

		if ( print_q != 0 )
		{
			TLatex* t = new TLatex(print_q, print_x, Form("x = %f", h_xq2->GetXaxis()->GetBinCenter(ix+1)));
			t->SetTextAlign(11);
			t->SetTextSize(0.02);
			t_a1_x.push_back(t);

			x_col.push_back(h_xq2->GetXaxis()->GetBinCenter(ix+1));
		}
	}

	outfile.close();

	std::vector<TGraphErrors*> g_a1_q2;
	for ( int i = 0; i < q2_col.size(); i ++ )
	{
		// a1
		TGraphErrors* g = new TGraphErrors(q2_col[i].size(), &q2_col[i][0], &a1_col[i][0], 0, &er_col[i][0]);
		g->SetMarkerStyle(20);
		g->SetMarkerColor(kBlue);
		g->SetLineColor(kBlue);

		g->SetTitle("");
		// g->SetTitle(Form("x = %f", x_col[i]));
		gStyle->SetTitleFontSize(0.08);
		g->GetXaxis()->SetTitle(Form("Q^{2} (GeV/c^{2})^{2}"));
		g->GetYaxis()->SetTitle(Form("A_{1}^{n} - 3 #times log_{10}(x)"));
		// g->GetYaxis()->SetTitle(Form("A_{1}^{p}"));
		g_a1_q2.push_back(g);
		format_graph(g);
	}

	std::vector<TGraphErrors*> g_a1_q2_low_acp;
	for ( int i = 0; i < q2_col_low_acp.size(); i ++ )
	{
		// a1 - low acceptance
		TGraphErrors* g = new TGraphErrors(q2_col_low_acp[i].size(), &q2_col_low_acp[i][0], &a1_col_low_acp[i][0], 0, 0);
		g->SetMarkerStyle(24);
		g->SetMarkerColor(kBlue);
		g_a1_q2_low_acp.push_back(g);
	}

	int n_col = 4;
    int n_row = g_a1_q2.size() / n_col;
    if ( n_row < (double)g_a1_q2.size() / n_col )
        n_row ++;

	TCanvas* c_a1_q2 = new TCanvas("c_a1_q2", "A1 vs Q2", 1800, 1400);
	c_a1_q2->Divide(n_row, n_col);

	for ( int i = 0; i < g_a1_q2.size(); i ++ )
	{
		c_a1_q2->cd(i+1);
		gPad->SetLogx();
		gPad->SetTickx(1);
  		gPad->SetTicky(1);
  		gPad->SetLeftMargin(0.2);
  		gPad->SetBottomMargin(0.2);
		// g_a1_q2[i]->GetXaxis()->SetLimits(1,1e5);
		// g_a1_q2[i]->GetYaxis()->SetTitleOffset(2);
		g_a1_q2[i]->Draw("AP SAME");

		// g_a1_q2[i]->GetXaxis()->SetLabelSize(15);
		// g_a1_q2[i]->GetXaxis()->SetTitleSize(15);
		// g_a1_q2[i]->GetYaxis()->SetLabelSize(15);
		// g_a1_q2[i]->GetYaxis()->SetTitleSize(15);

		g_a1_q2[i]->SetMinimum(a1_col[i][0]-er_col[i][0]*2);
		g_a1_q2[i]->SetMaximum(a1_col[i][0]+er_col[i][0]*2);

		// c_a1_q2->Update();
		// gPad->Update();
		// TLatex* t = new TLatex(gPad->GetUxmax(), gPad->GetUymax(), Form("x = %f", x_col[i]));
		// std::cout << (gPad->GetUxmin() + gPad->GetUymin())/2. << " " << gPad->GetUxmin() << " " << gPad->GetUymin() << std::endl;
		// TLatex* t = new TLatex((gPad->GetUxmin()+gPad->GetUxmax())/2., gPad->GetUymax()*0.3, Form("x = %f", x_col[i]));
		// t->SetTextAlign(11);
		// t->SetTextSize(0.05);
		// t->Draw("SAME");
	}

	TCanvas* c_compare = new TCanvas("c_comp","", 800, 1000);
	c_compare->SetLogx();
	// c_compare->SetTickx(1);
  	// c_compare->SetTicky(1);
	// c_compare->SetLogy();
	c_compare->SetLeftMargin(0.13);
	g_a1_q2[0]->Draw("AP");
	g_a1_q2[0]->SetMinimum(0);
	g_a1_q2[0]->SetMaximum(12);
	g_a1_q2[0]->GetXaxis()->SetLimits(1,3e4);
	// g_a1_q2[0]->GetYaxis()->SetLimits(-1,1);
	// g_a1_q2[0]->SetMinimum(-1);
	// g_a1_q2[0]->SetMaximum(1);

	c_compare->Update();
	t_a1_x[0]->Draw("SAME");

	for ( int i = 1; i < g_a1_q2.size(); i ++ )
	{
		g_a1_q2[i]->Draw("P SAME");
		// g_a1_q2[i]->GetYaxis()->SetRange(-1,1);

		g_a1_q2[i]->SetMinimum(0);
		g_a1_q2[i]->SetMaximum(14);
		g_a1_q2[i]->GetXaxis()->SetLimits(1,2e5);
	}

	for ( int i = 0; i < g_a1_q2_low_acp.size(); i ++ )
		g_a1_q2_low_acp[i]->Draw("P SAME");

	for ( int i = 0; i < t_a1_x.size(); i ++ )
		t_a1_x[i]->Draw("SAME");

	TLegend* leg = new TLegend(0.6, 0.80, 0.88, 0.88);
	leg->AddEntry(g_a1_q2[0], Form("%.0f #times %0.f GeV", e_lepton, e_nucleon), "PE");
	leg->SetTextSize(0.03);
	leg->SetBorderSize(0);
	leg->SetFillColor(0);
	leg->Draw();

	return;
}