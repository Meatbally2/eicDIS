#ifndef HISTMANAGER_HH
#define HISTMANAGER_HH

#include <vector>
#include <string>
#include <TH1F.h>
#include <TH2F.h>

#include "../GlobalUtil/drawHelper.h"

class HistManager{

public:

    HistManager();
    ~HistManager();

    struct HistGroup1D
    {
        TH1F* sum;
        TH1F* ind[4][2]; // 4 groups, n vs p
    };

    struct HistGroup2D
    {
        TH2F* sum;
        TH2F* ind[4][2]; // 4 groups, n vs p
    };

    std::vector<HistGroup1D> hist_group_1d; 
    std::vector<HistGroup2D> hist_group_2d;

    HistGroup1D MakeHistograms(std::string name, std::string title, int n_bins, double x_min, double x_max);
    HistGroup1D BookHistograms(std::string name, std::string title, int n_bins, double x_min, double x_max);
    HistGroup2D MakeHistograms(std::string name, std::string title, double n_bin_x, double x_min, double x_max, double n_bin_y, double y_min, double y_max);
    HistGroup2D BookHistograms(std::string name, std::string title, double n_bin_x, double x_min, double x_max, double n_bin_y, double y_min, double y_max);
    HistGroup2D BookMixedTH2(std::string name, std::string title, int n_bin_x, double x_min, double x_max, bool x_log, int n_bin_y, double y_min, double y_max, bool y_log);
    
    void SumHistograms(int group_index, double n_scale, double p_scale);
};

#endif
