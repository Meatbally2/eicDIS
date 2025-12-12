#ifndef DRAWMANAGER_HH
#define DRAWMANAGER_HH


#include "ePIC_style.c"

class DrawManager{

public:

	DrawManager(std::string type_, std::string energy_, std::string campaign_);
    ~DrawManager();

    void SetEPIC();
    void LableAndCollect(TCanvas* &c);
    void SaveToTree(TFile* &outFile);
    void SetLumi(double L);
    void SetQ2range(double Q2_min, double Q2_max);
    void SetQ2min(double Q2_min);

    std::string type;
    std::string energy;
    std::string campaign;
    std::vector<TCanvas*> canvas_list;

    bool setLumi = false;
    bool setQ2range = false;
    bool setQ2min = false;
    double lumi = 0; // in fb-1
    double Q2min = 0;
    double Q2max = 0;

private:
};

#endif
