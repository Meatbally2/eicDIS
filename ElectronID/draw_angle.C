void draw_angle(TCanvas* &c, TH2F* &h, double angle, double Ee, double En)
{
    double theta_deg = 180-angle; // formula from yellow report assume 0 deg is forward scattering for electron
    double theta = theta_deg * M_PI / 180.0;

    // double Ee = 10.0; 
    // double En = 166.0;

    // double Ee = 18.0; 
    // double En = 275.0;
    
    TAxis* xa = h->GetXaxis();
    int nx = xa->GetNbins();
    double last_Q2 = 0;

    TGraph* g_theta = new TGraph();
    for (int i = 1; i <= nx; ++i) {
        double low = xa->GetBinLowEdge(i);
        double high = xa->GetBinUpEdge(i);
        double center = xa->GetBinCenter(i);

        double xB[3] = {low, center, high};
        for ( int j = 0; j < 3; ++j) 
        {
            double Eprime = 2*Ee*En*xB[j] / (En*(1+cos(theta))*xB[j] + Ee*(1-cos(theta)));
            double Q2 = 2*Ee*Eprime*(1 - cos(theta));
            g_theta->AddPoint(xB[j], Q2);
            // std::cout << "xB: " << xB[j] << ", Q2: " << Q2 << std::endl;

            last_Q2 = Q2;
        }
    }

    g_theta->SetLineColor(kRed);
    g_theta->SetLineWidth(2);
    g_theta->SetLineStyle(7);

    c->cd();
    g_theta->Draw("L SAME");

    TLatex *lt = new TLatex(1.05, last_Q2, Form("#theta_{e} = %.0f#circ", angle));
    // lt->SetNDC(kTRUE);
    lt->SetTextFont(42);
    lt->SetTextSize(0.03);
    lt->SetTextColor(kRed);
    lt->Draw();


    return;
}