#include "bin.h"

const int n_run = 2;
std::string leg_name[n_run] = {"ep", "en"};
int marker[n_run] = {20, 22};
int color[n_run] = {kRed, kBlue};


std::vector<std::vector<KEbin*>> create_bins(TH2F* h_xq2, double eE, double nE)
{
    std::vector<std::vector<KEbin*>> bins;
    for ( int ix = 0; ix < h_xq2->GetXaxis()->GetNbins(); ix ++ )
    {
        std::vector<KEbin*> new_bin_set;
        for ( int iq = 0; iq < h_xq2->GetYaxis()->GetNbins(); iq ++ )
        {
            KEbin* new_bin = new KEbin(h_xq2->GetXaxis()->GetBinCenter(ix+1), h_xq2->GetYaxis()->GetBinCenter(iq+1), eE, nE);
            new_bin_set.emplace_back(std::move(new_bin));
        }
        bins.emplace_back(std::move(new_bin_set));
    }

    return bins;
}

void load_bins_from_text_simple(std::string file, TH2F* h_xq2, std::vector<std::vector<KEbin*>> &bins)
{
    std::ifstream data_file(file);

    std::string line;
    while ( getline(data_file, line) )
    {
        double x = 0;
        double q2 = 0;
        double n = 0;

        std::stringstream ss(line);
        ss >> x >> q2 >> n;

        int xbin = h_xq2->GetXaxis()->FindBin(x);
        int qbin = h_xq2->GetYaxis()->FindBin(q2);

        if ( n != 0 )
            bins[xbin-1][qbin-1]->set_bin(x, q2, n, 0);
            // bins[xbin-1][qbin-1]->set_bin(x, q2, n*(0.4/10), n_raw);


    }

    data_file.close();

    return;
}

void load_bins_from_text(std::string file, TH2F* h_xq2, std::vector<std::vector<KEbin*>> &bins)
{
    std::ifstream data_file(file);

    std::string line;
    while ( getline(data_file, line) )
    {
        double x = 0;
        double q2 = 0;
        double weighted_x = 0;
        double weighted_q2 = 0;
        double y = 0;
        double nu = 0;
        double n = 0;
        double n_raw = 0;

        std::stringstream ss(line);
        ss >> x >> q2 >> weighted_x >> weighted_q2 >> y >> nu >> n >> n_raw;

        int xbin = h_xq2->GetXaxis()->FindBin(x);
        int qbin = h_xq2->GetYaxis()->FindBin(q2);

        if ( n != 0 )
            bins[xbin-1][qbin-1]->set_bin(x, q2, n, n_raw);
            // bins[xbin-1][qbin-1]->set_bin(x, q2, n*(0.4/10), n_raw);


    }

    data_file.close();

    return;
}

void format_graph(TGraph* g)
{
    // g->GetXaxis()->SetTitle("Q^{2} (GeV/c^{2})^{2}");
    g->GetXaxis()->CenterTitle();
    g->GetXaxis()->SetTitleOffset(1.4);
   
    // g->GetYaxis()->SetTitle(ylabel.c_str());
    g->GetYaxis()->CenterTitle();
    g->GetYaxis()->SetTitleOffset(1.4);

    g->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    g->GetXaxis()->SetLabelSize(24);
    g->GetXaxis()->SetTitleSize(28);
    g->GetXaxis()->SetTitleFont(43);

    g->GetYaxis()->SetNdivisions(505);
    g->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    g->GetYaxis()->SetLabelSize(24);
    g->GetYaxis()->SetTitleSize(28);
    g->GetYaxis()->SetTitleFont(43);

    // g->GetYaxis()->ChangeLabel(1, -1, -1, -1, -1, -1, " ");

    // g->SetMarkerStyle(marker[style_id]);
    // g->SetMarkerColor(color[style_id]);
    g->SetMarkerSize(0.2);

    return;
}

void format_graph(TGraph* g, std::string ylabel, int style_id)
{
    // g->GetXaxis()->SetTitle("Q^{2} (GeV/c^{2})^{2}");
    // g->GetXaxis()->CenterTitle();
    // g->GetXaxis()->SetTitleOffset(1.2);
   
    // g->GetYaxis()->SetTitle(ylabel.c_str());
	// g->GetYaxis()->CenterTitle();

	g->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
	g->GetXaxis()->SetLabelSize(20);
	g->GetXaxis()->SetTitleSize(20);
	g->GetXaxis()->SetTitleFont(43);

    g->GetYaxis()->SetNdivisions(505);
	g->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
	g->GetYaxis()->SetLabelSize(20);
	g->GetYaxis()->SetTitleSize(20);
	g->GetYaxis()->SetTitleFont(43);

    g->SetMarkerStyle(marker[style_id]);
    g->SetMarkerColor(color[style_id]);
    // g->SetMarkerSize(0.5);

	return;
}