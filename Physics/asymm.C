
#include "calculate_a1.h"
#include "bin.h"
#include "BinHelper.h"

#include "../GlobalUtil/Constants.hh"
#include "../GlobalUtil/drawHelper.h"
#include "../GlobalUtil/DrawManager.cc"

const double ePol = 0.7;
const double iPol = 0.7;

double bin_shift = 0;

void asymm(int beam_type, int Ee, int Eh)
{
    // ePIC plotting style setup
    std::string type_title[5] = {"e^{3}He", "ep", "#gammap", "ep w. BeamBG", "ep"};
    std::string energy_title = beam_type ? Form("%dx%d GeV", Ee, Eh) : Form("%dx%d GeV/A", Ee, Eh);
    std::string campaign = beam_type ? "25.10.2" : "25.10.0";
    if ( beam_type == 4 )
        campaign = "25.10.4";
    DrawManager* draw_manager = new DrawManager(type_title[beam_type], energy_title, campaign);
    draw_manager->SetEPIC();

    draw_manager->SetDISrange(0.01, 0.95, 4, 2);

    double text_lumi = 1.0; // fb^-1
    if ( beam_type == 0 && Ee == 10 && Eh == 166 )
        text_lumi = 1.5; // fb^-1
    if ( beam_type == 4 && Ee == 10 && Eh == 130 )
        text_lumi = 1.0; // fb^-1
    if ( beam_type == 4 && Ee == 10 && Eh == 250 )
        text_lumi = 2.5; // fb^-1
    draw_manager->SetLumi(text_lumi);
     
	const double m_nucleon = beam_type == 0 ? MASS_NEUTRON : MASS_PROTON;

	std::string date = beam_type == 0 ? "en_25_10_2" : "ep_25_10_0";
	if ( beam_type == 4 )
		date = "ep_25_10_4";
	TFile* file = new TFile(Form("../data/%s/FullRecon_%dx%d_new.root", date.c_str(), (int)Ee, (int)Eh));
	if ( !file )
	{
		std::cout << "file not found" << std::endl;
		return;
	}

	TH2F* h_xq2 = (TH2F*) file->Get("H_XQ_REC");
	// TH2F* h_xq2 = (TH2F*) file->Get("h_xq2_mc_n");

	if ( h_xq2 == NULL )
	{
		std::cout << "histogram not found" << std::endl;
		return;
	}

	std::vector<std::vector<KEbin*>> bins = create_bins(h_xq2, Ee, Eh);

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

	double Q2max = 4*Eh*Ee;
	double nu_max = Q2max/(2*m_nucleon);

	std::vector<TLatex*> t_a1_x;
	std::vector<double> x_col;
	std::vector<std::vector<double>> q2_col;
	std::vector<std::vector<double>> a1_col;
	std::vector<std::vector<double>> er_col;
	std::vector<std::vector<double>> q2_col_low_acp;
	std::vector<std::vector<double>> a1_col_low_acp;
	std::vector<std::vector<double>> er_col_low_acp;

	std::string setting = beam_type == 0 ? Form("eHe3_%dx%d", (int)Ee, (int)Eh) : Form("ep_%dx%d", (int)Ee, (int)Eh);
	ofstream outfile(beam_type == 0 ? Form("../data/%s/%s_xq2_n_collection_new.txt", date.c_str(), setting.c_str()) : Form("../data/%s/%s_xq2_p_collection_new.txt", date.c_str(), setting.c_str()));

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

			double x = h_xq2->GetXaxis()->GetBinCenter(ix+1);
			double q = h_xq2->GetYaxis()->GetBinCenter(iq+1);
			double print_a1 = beam_type == 0 ? find_a1n(x, q) : find_a1p(x, q);
			double a1 = print_a1 - 3 * TMath::Log10(x);

			// cout << h_xq2->GetXaxis()->GetBinCenter(ix+1) << " " << h_xq2->GetYaxis()->GetBinCenter(iq+1) << endl;
			if ( beam_type == 0 )
				bins[ix][iq]->process_bin(m_nucleon, ePol, iPol, print_a1);
			else
				bins[ix][iq]->process_bin(m_nucleon, ePol, iPol, print_a1);
			// cout <<  " end bin " << endl;

			if ( bins[ix][iq]->y <= 0.01 || bins[ix][iq]->y >= 0.95)
				continue;

			double err = bins[ix][iq]->unc_a1;
			double a2err = bins[ix][iq]->unc_a2;
			double a1_sys = bins[ix][iq]->sys_a1;
			double a1_norm = bins[ix][iq]->norm_a1;
			double a2_sys = bins[ix][iq]->sys_a2;
			double a2_norm = bins[ix][iq]->norm_a2;

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
			
			outfile << x << " " << q << " " << bins[ix][iq]->n_events << " " << print_a1 << " " << err << " " << a2err << " " << a1_sys << " " << a1_norm << " " << a2_sys << " " << a2_norm << std::endl;
			
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
		if ( beam_type == 0 )
			g->GetYaxis()->SetTitle(Form("A_{1}^{n} - 3 #times log_{10}(x)"));
		else
			g->GetYaxis()->SetTitle(Form("A_{1}^{p} - 3 #times log_{10}(x)"));
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
	leg->AddEntry(g_a1_q2[0], Form("%d #times %d GeV", Ee, Eh), "PE");
	leg->SetTextSize(0.03);
	leg->SetBorderSize(0);
	leg->SetFillColor(0);
	leg->Draw();

	return;
}