#include "HistManager.hh"

HistManager::HistManager() {
}

HistManager::~HistManager() {
}

HistManager::HistGroup1D HistManager::MakeHistograms(std::string name, std::string title, int n_bins, double x_min, double x_max)
{
    HistGroup1D hist_1d;
    hist_1d.sum = new TH1F(name.c_str(), title.c_str(), n_bins, x_min, x_max);
    hist_1d.sum->SetStats(0);
    hist_1d.sum->GetXaxis()->CenterTitle();
    hist_1d.sum->GetYaxis()->CenterTitle();

    for ( int i = 0; i < 4; i ++ )
    {
        for ( int j = 0; j < 2; j ++ )
        {
            hist_1d.ind[i][j] = new TH1F(Form("%s_%d_%d", name.c_str(), i, j), title.c_str(), n_bins, x_min, x_max);
        }
    }
    hist_group_1d.push_back(hist_1d);
    return hist_1d;
}

HistManager::HistGroup1D HistManager::BookHistograms(std::string name, std::string title, int n_bins, double x_min, double x_max)
{
    HistGroup1D hist_1d;
    hist_1d.sum = BookTH1(name.c_str(), title.c_str(), n_bins, x_min, x_max);

    for ( int i = 0; i < 4; i ++ )
    {
        for ( int j = 0; j < 2; j ++ )
        {
            hist_1d.ind[i][j] = BookTH1(Form("%s_%d_%d", name.c_str(), i, j), title.c_str(), n_bins, x_min, x_max);
        }
    }
    hist_group_1d.push_back(hist_1d);

    return hist_1d;
}

HistManager::HistGroup2D HistManager::MakeHistograms(std::string name, std::string title, double n_bin_x, double x_min, double x_max, double n_bins_y, double y_min, double y_max)
{
    HistGroup2D hist_2d;
    hist_2d.sum = new TH2F(name.c_str(), title.c_str(), n_bin_x, x_min, x_max, n_bins_y, y_min, y_max);
    hist_2d.sum->SetStats(0);
    hist_2d.sum->GetXaxis()->CenterTitle();
    hist_2d.sum->GetYaxis()->CenterTitle();

    for ( int i = 0; i < 4; i ++ )
    {
        for ( int j = 0; j < 2; j ++ )
        {
            hist_2d.ind[i][j] = new TH2F(Form("%s_%d_%d", name.c_str(), i, j), title.c_str(), n_bin_x, x_min, x_max, n_bins_y, y_min, y_max);
        }
    }
    hist_group_2d.push_back(hist_2d);

    return hist_2d;
}

HistManager::HistGroup2D HistManager::BookHistograms(std::string name, std::string title, double n_bin_x, double x_min, double x_max, double n_bins_y, double y_min, double y_max)
{
    HistGroup2D hist_2d;
    hist_2d.sum = BookTH2(name.c_str(), title.c_str(), n_bin_x, x_min, x_max, n_bins_y, y_min, y_max, kLightTemperature);

    for ( int i = 0; i < 4; i ++ )
    {
        for ( int j = 0; j < 2; j ++ )
        {
            hist_2d.ind[i][j] = BookTH2(Form("%s_%d_%d", name.c_str(), i, j), title.c_str(), n_bin_x, x_min, x_max, n_bins_y, y_min, y_max, kLightTemperature);
        }
    }
    hist_group_2d.push_back(hist_2d);

    return hist_2d;
}

void HistManager::SumHistograms(int group_index, double n_scale, double p_scale)
{
    for ( auto& h : hist_group_1d )
    {
        h.ind[group_index][0]->Scale(n_scale);
        h.ind[group_index][1]->Scale(p_scale);

        h.sum->Add(h.ind[group_index][0]);
        h.sum->Add(h.ind[group_index][1]);
    }

    for ( auto& h : hist_group_2d )
    {
        h.ind[group_index][0]->Scale(n_scale);
        h.ind[group_index][1]->Scale(p_scale);

        h.sum->Add(h.ind[group_index][0]);
        h.sum->Add(h.ind[group_index][1]);
    }
}