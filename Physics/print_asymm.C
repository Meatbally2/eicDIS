#include "calculate_a1.h"
#include "bin.h"
#include "BinHelper.h"

#include "../GlobalUtil/Constants.hh"
#include "../GlobalUtil/drawHelper.h"
#include "../GlobalUtil/DrawManager.cc"

// 

// const double e_lepton = 18;
// const double e_nucleon = 275.;
// std::string setting = "18x275";

// const double e_lepton = 10;
// const double e_nucleon = 100.;
// std::string setting = "10x100";

// const double e_lepton = 5;
// const double e_nucleon = 41.;
// std::string setting = "5x41";

// const double m_nucleon = MASS_PROTON;

// 

// const double e_lepton = 18;
// const double e_nucleon = 110.;
// std::string setting = "18x110";

const double e_lepton = 10;
const double e_nucleon = 166.;
std::string setting = "10x166";

const double m_nucleon = MASS_NEUTRON;

// 

const double ePol = 0.7;
const double iPol = 0.7;

void print_asymm()
{
    std::string type = m_nucleon == MASS_NEUTRON ? "n" : "p";

	TH2F* h_xq2 = BookTH2("xq2",  ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);

    std::vector<std::vector<KEbin*>> bins = create_bins(h_xq2, e_lepton, e_nucleon);
    if ( type == "p" )
        load_bins_from_text(Form("../epic_data/%s_xq2_collection.txt", setting.c_str()), h_xq2, bins);
    else
        load_bins_from_text_simple(Form("../data/en_25_10_2/%s_xq2_%s_collection.txt", setting.c_str(), type.c_str()), h_xq2, bins);

    // load structre functions
    std::ifstream fs_file(Form("fs_%s_table.txt", type.c_str()));

    std::string line;
    while ( getline(fs_file, line) )
    {
        double x = 0;
        double q2 = 0;
        double f1 = 0;
        double f2 = 0;
        double a1 = 0;
        double g1 = 0;
        double f2err = 0;
        double f1err = 0;

        std::stringstream ss(line);
        ss >> x >> q2 >> f2 >> f1 >> a1 >> g1 >> f2err >> f1err;

        int xbin = h_xq2->GetXaxis()->FindBin(x);
        int qbin = h_xq2->GetYaxis()->FindBin(q2);

        bins[xbin-1][qbin-1]->save_f2(f2, f2err);
        bins[xbin-1][qbin-1]->save_f1(f1, f1err);

        // std::cout << x << " " << q2 << " " << f2 << " " << f1 << std::endl;
    }
    fs_file.close();

    // caculate a1, g1
    ofstream outfile(Form("../data/en_25_10_2/%s_xq2_%s_asymm_table.txt",setting.c_str(), type.c_str()));

    for ( int ix = 0; ix < h_xq2->GetXaxis()->GetNbins(); ix ++ )
    {
        for ( int iq = 0; iq < h_xq2->GetYaxis()->GetNbins(); iq ++ )
        {
            // if ( !bins[ix][iq]->is_averaged )
			// 	continue;
            
            bins[ix][iq]->process_bin(m_nucleon, ePol, iPol);

            if ( bins[ix][iq]->y <= 0.01 || bins[ix][iq]->y >= 0.95)
				continue;

            h_xq2->SetBinContent(ix+1, iq+1, bins[ix][iq]->n_events);
            // h_xq2_not_scaled->SetBinContent(ix+1, iq+1, bins[ix][iq]->n_raw);

            double x = h_xq2->GetXaxis()->GetBinCenter(ix+1);
			double q = h_xq2->GetYaxis()->GetBinCenter(iq+1);

            double a1;
            double g1;

            if ( type == "p" )
            {
                a1 = find_a1p(x, q);
			    g1 = find_a1p(x, q) * bins[ix][iq]->f1;
            }
            else
            {
                a1 = find_a1n(x, q);
			    g1 = find_a1n(x, q) * bins[ix][iq]->f1;
            }
            
            double a1err = bins[ix][iq]->unc_a1;
            double a1sys = a1 * sqrt(1.5*1.5+3*3) * 0.01;
            double a1norm = a1 * 3.5 * 0.01;
            double a1tot = sqrt(a1err*a1err + a1sys*a1sys + a1norm*a1norm);

            double g1err = a1tot * bins[ix][iq]->f1;
            // double g1err_all = g1err+abs(g1)*sqrt((a1tot/a1)*(a1tot/a1) + (bins[ix][iq]->f1err/bins[ix][iq]->f1)*(bins[ix][iq]->f1err/bins[ix][iq]->f1));
            double g1err_all = abs(g1)*sqrt((a1tot/a1)*(a1tot/a1) + (bins[ix][iq]->f1err/bins[ix][iq]->f1)*(bins[ix][iq]->f1err/bins[ix][iq]->f1));

            if ( bins[ix][iq]->n_events != 0 )
            {
                outfile << x << " " << q << " " << a1 << " " << a1err << " " << g1 << " " << g1err << " " << g1err_all << std::endl;
                // cout << x << " " << q << " " << a1 << " " << a1err << " " << g1 << " " << g1err << " " << g1err_all << " " << bins[ix][iq]->f1 << " " <<  bins[ix][iq]->f1err << endl;
            }
        }
    }
    outfile.close();

    return;
}