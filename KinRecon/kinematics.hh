#ifndef KINEMATICS_HH
#define KINEMATICS_HH

#include "../GlobalUtil/drawHelper.h"
#include "../GlobalUtil/Constants.hh"
#include <cmath>

class Kinematics{

public:

    Kinematics(std::string method_name, int draw_id);
    ~Kinematics();

    void reset();

    TH2F* h_xq;
    TH2F* h_xq_eff;

    TH1F* h_dq[4];
    TH1F* h_dx[4];
    TH1F* h_dy[4];

    TH1F* h_dq_ycut[6];
    TH1F* h_dx_ycut[6];
    TH1F* h_dy_ycut[6];

    TH2F* h_qvq;
    TH2F* h_xvx;
    TH2F* h_yvy;
    
    // void fill_histograms();
    // void draw_histograms();

    float p;
    float eta;
    float phi;

    float xB;
    float Q2;
    float W2;
    float y;
    float nu;

    // outgoing e
    float E;
    float theta;
    float del_h;

    void StyleHist(TH1F* h);
    void StyleHist(TH2F* h);
    void Process_and_QA(double mc_xB, double mc_y, double mc_Q2);
    void save_variables(std::vector<float> reco, double Mp);

    std::string method_name;
    int draw_id;

private:
    
};

#endif