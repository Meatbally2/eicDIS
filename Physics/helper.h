#include "calculate_a1.h"
#include "bin.h"
#include "BinHelper.h"

#include "../GlobalUtil/Constants.hh"
#include "../GlobalUtil/drawHelper.h"
#include "../GlobalUtil/DrawManager.cc"

void load_function(TH2F* h_xq2, std::vector<std::vector<KEbin*>> &bins, std::string name)
{
    // load structre functions
    // if ( name == "he3" )
    //     name = "d";

    std::ifstream fs_file(Form("functions/fs_%s_table.txt", name.c_str()));

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

        if ( name == "d" )
            f2 = f2 + bins[xbin-1][qbin-1]->f2p;

        // cout << name << " " << x << " " << q2 << " " << f2 << endl;

        if ( name == "p" )
        {
            bins[xbin-1][qbin-1]->f2p = f2;
            bins[xbin-1][qbin-1]->g1p = g1;
        }
            
        if ( name == "n" )
        {
            bins[xbin-1][qbin-1]->f2n = f2;
            bins[xbin-1][qbin-1]->g1n = g1;
        }
            
        // if ( name == "d" )
        if ( name == "he3" )
            bins[xbin-1][qbin-1]->f2_he3 = f2*3;
    }
    fs_file.close();

    return;
}