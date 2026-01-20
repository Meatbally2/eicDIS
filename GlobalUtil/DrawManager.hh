#ifndef DRAWMANAGER_HH
#define DRAWMANAGER_HH


#include "ePIC_style.C"

class DrawManager{

public:

	DrawManager(std::string type_, std::string energy_, std::string campaign_);
    DrawManager(std::string type_, std::string campaign_);
    ~DrawManager();

    void SetEPIC(std::string plot_version_ = "Internal");
    void LableAndCollect(TCanvas* &c, int draw_position = 0);
    void SaveToTree(TFile* &outFile);
    void SetCampaign(std::string campaign_);
    void SetLumi(double L);
    void SetDISrange(double Y_min, double Y_max, double W2_min, double Q2_min);
    void SetQ2range(double Q2_min, double Q2_max);
    void SetQ2min(double Q2_min);
    void SetPlotType(std::string plot_type_);
    // std::pair<double, double> PixelToNDC(int x_px, int y_px);
    // void LableAndCollectSpecial(TCanvas* &c);
    // void LableAndCollectSpecial2(TCanvas* &c);

    std::string type;
    std::string energy;
    std::string campaign;
    std::string plot_type;
    std::string plot_version;
    std::vector<TCanvas*> canvas_list;

    bool setLumi = false;
    bool setDIS = false;
    bool setQ2range = false;
    bool setQ2min = false;
    bool multiple_energies = false;
    double lumi = 0; // in fb-1
    double Ymin = 0;
    double Ymax = 0;
    double W2min = 0;
    double Q2min = 0;
    double Q2max = 0;

private:
};

#endif
