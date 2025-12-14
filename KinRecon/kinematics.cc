#include "kinematics.hh"

Kinematics::Kinematics(std::string method_name_, int draw_id_) 
{
    method_name = method_name_;
    draw_id = draw_id_;

    h_xq = BookTH2(Form("H_%s_Events", method_name.c_str()), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);
    h_xq_eff = BookTH2(Form("H_%s_Efficiency", method_name.c_str()), ";x;Q^{2} (GeV/c^{2})^{2}", n_x_bin, -5, 0, n_q_bin,  0, 5, kLightTemperature);

    std::string QrangeName[4] = {"lowQ", "midQ", "highQ", "ultraQ"};
    std::string Qrange[4] = {"Q^{2} < 10 (GeV/c^{2})^{2}", "10 (GeV/c^{2})^{2} #leq Q^{2} < 100 (GeV/c^{2})^{2}", "100 (GeV/c^{2})^{2} #leq Q^{2} < 1000 (GeV/c^{2})^{2}", "Q^{2} #geq 1000 (GeV/c^{2})^{2}"};
    for ( int i = 0; i < 4; i ++ )
    {
        h_dq[i] = new TH1F(Form("H_%s_dQ_%s", method_name.c_str(), QrangeName[i].c_str()), Form("%s;(Q^{2}_{reco}-Q^{2}_{true})/Q^{2}_{true};Counts", Qrange[i].c_str()), 100, -1, 1);
        h_dx[i] = new TH1F(Form("H_%s_dX_%s", method_name.c_str(), QrangeName[i].c_str()), Form("%s;(x_{reco}-x_{true})/x_{true};Counts", Qrange[i].c_str()), 100, -1, 1);
        h_dy[i] = new TH1F(Form("H_%s_dY_%s", method_name.c_str(), QrangeName[i].c_str()), Form("%s;(y_{reco}-y_{true})/y_{true};Counts", Qrange[i].c_str()), 100, -1, 1);

        StyleHist(h_dq[i]);
        StyleHist(h_dx[i]);
        StyleHist(h_dy[i]);
    }

    h_qvq = BookTH2(Form("H_%s_QVQ", method_name.c_str()), Form("%s;Q_{true}^{2} (GeV/c^{2})^{2};Q_{reco}^{2} (GeV/c^{2})^{2}", method_name.c_str()), 200, 0, 4, 200,  0, 4, kLightTemperature);
    h_xvx = BookTH2(Form("H_%s_XVX", method_name.c_str()), Form("%s;x_{true};x_{reco}", method_name.c_str()), 100, -4, 0, 100,  -4, 0, kLightTemperature);
    // h_yvy = BookTH2(Form("H_%s_YVY", method_name.c_str()), Form("%s;y_{true};y_{reco}", method_name.c_str()), 100, -2, 0, 100,  -2, 0, kLightTemperature);
    h_yvy = new TH2F(Form("H_%s_YVY", method_name.c_str()), Form("%s;y_{true};y_{reco}", method_name.c_str()), 100, 0, 1, 100,  0, 1);

    StyleHist(h_qvq);
    StyleHist(h_xvx);
    StyleHist(h_yvy);

    std::string YrangeName[6] = {"lowY_lowX", "midY_lowX", "highY_lowX, lowY_highX", "midY_highX", "highY_highX"};
    std::string Yrange[6] = {"0.01 < y #leq 0.035, x #leq 0.025", "0.035 < y #leq 0.2, x #leq 0.025", "0.2 < y < 0.95, x #leq 0.025", "0.01 < y #leq 0.035, x > 0.025", "0.035 < y #leq 0.2, x > 0.025", "0.2 < y < 0.95, x > 0.025"};
    for ( int i = 0; i < 6; i ++ )
    {
        h_dq_ycut[i] = new TH1F(Form("H_%s_dQ_%s", method_name.c_str(), YrangeName[i].c_str()), Form("%s;(Q^{2}_{reco}-Q^{2}_{true})/Q^{2}_{true};Counts", Yrange[i].c_str()), 100, -0.5, 0.5);
        h_dx_ycut[i] = new TH1F(Form("H_%s_dX_%s", method_name.c_str(), YrangeName[i].c_str()), Form("%s;(x_{reco}-x_{true})/x_{true};Counts", Yrange[i].c_str()), 100, -0.5, 0.5);
        h_dy_ycut[i] = new TH1F(Form("H_%s_dY_%s", method_name.c_str(), YrangeName[i].c_str()), Form("%s;(y_{reco}-y_{true})/y_{true};Counts", Yrange[i].c_str()), 100, -0.5, 0.5);

        StyleHist(h_dq_ycut[i]);
        StyleHist(h_dx_ycut[i]);
        StyleHist(h_dy_ycut[i]);
    }
}

Kinematics::~Kinematics() {
}

void Kinematics::StyleHist(TH2F* h) 
{
    h->SetStats(0);
    h->SetLineColor(draw_id);

    h->GetXaxis()->CenterTitle();
    h->GetXaxis()->SetTitleOffset(1.3);
   
    h->GetYaxis()->SetTitleOffset(1.6);
    h->GetYaxis()->CenterTitle();

    h->GetXaxis()->SetNdivisions(5);
    h->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    h->GetXaxis()->SetLabelSize(20);
    h->GetXaxis()->SetTitleSize(20);
    h->GetXaxis()->SetTitleFont(43);

    h->GetYaxis()->SetNdivisions(5);
    h->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    h->GetYaxis()->SetLabelSize(20);
    h->GetYaxis()->SetTitleSize(20);
    h->GetYaxis()->SetTitleFont(43);
}

void Kinematics::StyleHist(TH1F* h) 
{
    h->SetStats(0);
    h->SetLineColor(draw_id);

    h->GetXaxis()->CenterTitle();
    h->GetXaxis()->SetTitleOffset(1.3);
   
    h->GetYaxis()->SetTitleOffset(1.5);
    h->GetYaxis()->CenterTitle();

    h->GetXaxis()->SetNdivisions(5);
    h->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    h->GetXaxis()->SetLabelSize(20);
    h->GetXaxis()->SetTitleSize(20);
    h->GetXaxis()->SetTitleFont(43);

    h->GetYaxis()->SetNdivisions(5);
    h->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    h->GetYaxis()->SetLabelSize(20);
    h->GetYaxis()->SetTitleSize(20);
    h->GetYaxis()->SetTitleFont(43);
}

void Kinematics::Process_and_QA(double mc_xB, double mc_y, double mc_Q2) 
{
    if ( xB == -999 )
        return;

    h_xq->Fill(xB, Q2);

    if ( h_xq_eff->GetXaxis()->FindBin(mc_xB) == h_xq_eff->GetXaxis()->FindBin(xB) )
        if ( h_xq_eff->GetYaxis()->FindBin(mc_Q2) == h_xq_eff->GetYaxis()->FindBin(Q2) )
            h_xq_eff->Fill(xB, Q2);

    int Qrange = 0;
    if ( mc_Q2 > 1000 )
        Qrange = 3;
    else if ( mc_Q2 > 100 )
        Qrange = 2;
    else if ( mc_Q2 > 10 )
        Qrange = 1;

    h_dq[Qrange]->Fill((Q2-mc_Q2)/mc_Q2);
    h_dx[Qrange]->Fill((xB-mc_xB)/mc_xB);
    h_dy[Qrange]->Fill((y-mc_y)/mc_y);

    h_qvq->Fill(mc_Q2, Q2);
    h_xvx->Fill(mc_xB, xB);
    h_yvy->Fill(mc_y, y);

    int Yrange = xB > 0.025 ? 3 : 0;
    if ( y > 0.2 )
        Yrange += 2;
    else if ( y > 0.035 )
        Yrange += 1;

    h_dq_ycut[Yrange]->Fill((Q2-mc_Q2)/mc_Q2);
    h_dx_ycut[Yrange]->Fill((xB-mc_xB)/mc_xB);
    h_dy_ycut[Yrange]->Fill((y-mc_y)/mc_y);

    // if ( y > 0.95  )
    //     cout << method_name << " Q2 " << mc_Q2 << " " << Q2 << " xB " << mc_xB << " " << xB << " y " << mc_y << " " << y << endl;

    return;
}

void Kinematics::save_variables(std::vector<float> reco, double Mp) 
{
    xB = reco[0];
    y = reco[1];
    Q2 = reco[2];
    W2 = Mp*Mp + Q2 * (1./xB - 1);
    nu = Q2/(2*Mp*xB);
}

void Kinematics::reset()
{
    p = -999;
    eta = -999;
    phi = -999;

    xB = -999;
    Q2 = -999;
    W2 = -999;
    y = -999;
    nu = -999;

    E = -999;
    theta = -999;
    del_h = -999;
}