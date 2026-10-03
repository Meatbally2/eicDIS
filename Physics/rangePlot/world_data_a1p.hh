// masses in GeV
#define MASS_ELECTRON 0.000511
#define MASS_PROTON   0.93827
#define MASS_NEUTRON  0.93957
#define MASS_HELIUM3  2.8094

struct WorldA1pDataSet
{
    std::string name;
    std::vector<double> x;
    std::vector<double> q2;
    int plot_color;
    int plot_marker;
};

std::vector<WorldA1pDataSet> a1_data;

struct WorldA1pValueSet
{
    std::string name;
    std::vector<double> x;
    std::vector<double> q2;
    std::vector<double> a1;
    std::vector<double> a1_err;
    int plot_color;
    int plot_marker;
};

std::vector<WorldA1pValueSet> a1_value;

// SLAC E80, p, A_LL, 1976, https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.37.1261
void slac_e80()
{
    std::vector<double> q2 = {1.680, 2.735, 1.418};
    std::vector<double> ALL = {0.191, 0.215, 0.141};
    std::vector<double> ALL_err = {0.057, 0.089, 0.058};

    std::vector<double> W = {2.059, 2.519, 2.560};
    std::vector<double> x;
    for ( int i = 0; i < 3; i ++ )
    {
        x.push_back(q2[i] / (W[i]*W[i] - MASS_PROTON*MASS_PROTON + q2[i]));
        // printf("%f %f\n", x[i], ALL[i]);
    }

    a1_data.push_back({"SLAC E80", x, q2, kBlack, 20});
}

// SLAC E130, p, A_LL, 1978, https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.41.70
void slac_e130()
{
    std::vector<double> q2 = {4.09, 1.68, 1.68, 2.74, 2.74, 1.02, 1.42, 2.95, 1.70};
    std::vector<double> ALL = {0.213, 0.131, 0.188, 0.062, 0.148, 0.109, 0.177};
    std::vector<double> ALL_err = {0.057, 0.039, 0.066, 0.031, 0.073, 0.081, 0.057};

    std::vector<double> W = {2.0, 3.0, 3.0, 3.0, 3.0, 4.9, 5.0, 5.0, 10.0};
    std::vector<double> x;
    for ( int i = 0; i < 9; i ++ )
    {
        x.push_back(q2[i] / (W[i]*W[i] - MASS_PROTON*MASS_PROTON + q2[i]));
        // printf("%f %f\n", x[i], ALL[i]);
    }

    a1_data.push_back({"SLAC E130", x, q2, kGreen + 4, 21});
}

// CERN EMC, p, A1 = ALL / D, 1988, https://cds.cern.ch/record/183730/files/198801416.pdf
void cern_emc()
{
    std::vector<double> x = {0.015, 0.025, 0.035, 0.050, 0.078, 0.124, 0.175, 0.248, 0.344, 0.466};
    std::vector<double> q2 = {3.5, 4.5, 6.0, 8.0, 10.3, 12.9, 15.2, 18.0, 22.5, 29.5};
    std::vector<double> A1 = {0.021, 0.087, 0.013, 0.094, 0.139, 0.169, 0.360, 0.469, 0.517, 0.657};
    std::vector<double> A1_stat_err = {0.035, 0.043, 0.054, 0.048, 0.049, 0.063, 0.087, 0.088, 0.141, 0.175};
    std::vector<double> A1_syst_err = {0.017, 0.022, 0.024, 0.028, 0.037, 0.045, 0.057, 0.065, 0.068, 0.065};
    double A1_norm_err = 0.096; // 9.6% normalization error

    a1_data.push_back({"CERN EMC", x, q2, kBlue, 22});
}

// CERN SMC, p, A1, A2, 1997, https://www.sciencedirect.com/science/article/pii/S0370269397011064?via%3Dihub
// They took a1_data in 1994 and 1996, included here should be a combination of both.
void cern_smc()
{
    std::vector<double> x = {0.005, 0.008, 0.014, 0.025, 0.035, 0.049, 0.077, 0.122, 0.173, 0.242, 0.342, 0.480};
    std::vector<double> q2 = {1.3, 2.1, 3.6, 5.7, 7.8, 10.4, 14.9, 21.3, 27.8, 35.6, 45.9, 58.0};
    std::vector<double> A1 = {0.017, 0.047, 0.035, 0.058, 0.067, 0.115, 0.176, 0.267, 0.318, 0.400, 0.568, 0.658};
    std::vector<double> A1_stat_err = {0.018, 0.016, 0.014, 0.018, 0.022, 0.019, 0.019, 0.025, 0.035, 0.036, 0.058, 0.079};
    std::vector<double> A1_syst_err = {0.003, 0.004, 0.003, 0.005, 0.005, 0.008, 0.013, 0.018, 0.021, 0.028, 0.042, 0.055};

    a1_data.push_back({"CERN SMC", x, q2, kRed + 1, 23});
}

// SLAC E143, p, A_LL, some region has A_LT, 1998, https://journals.aps.org/prd/pdf/10.1103/PhysRevD.58.112003
void slac_e143()
{
    std::vector<double> x = {
        0.031, 0.035, 0.039, 0.044, 0.049, 0.056, 0.063, 0.071, 0.079, 0.090,
        0.101, 0.113, 0.128, 0.144, 0.162, 0.182, 0.205, 0.230, 0.259, 0.292,
        0.329, 0.370, 0.416, 0.468, 0.526, 0.592, 0.666, 0.749
    };
    std::vector<double> q2 = {
        1.27, 1.40, 1.52, 1.65, 1.78, 1.91, 2.04, 2.19, 2.41, 2.55,
        2.85, 3.13, 3.41, 3.71, 4.03, 4.34, 4.15, 4.37, 5.26, 5.53,
        6.01, 6.29, 6.56, 6.79, 7.72, 7.97, 9.26, 9.52
    };

    std::vector<double> A1 = {
        0.063, 0.122, 0.083, 0.110, 0.123, 0.130, 0.123, 0.152, 0.180, 0.168,
        0.214, 0.217, 0.232, 0.238, 0.272, 0.294, 0.314, 0.325, 0.437, 0.432,
        0.441, 0.471, 0.616, 0.654, 0.713, 0.649, 0.612, 0.914
    };
    std::vector<double> A1_stat_err = {
        0.034, 0.025, 0.023, 0.021, 0.020, 0.019, 0.018, 0.018, 0.018, 0.016,
        0.016, 0.017, 0.017, 0.018, 0.019, 0.020, 0.020, 0.022, 0.026, 0.030,
        0.036, 0.041, 0.049, 0.059, 0.090, 0.118, 0.182, 0.273
    };
    std::vector<double> A1_syst_err = {
        0.009, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.009, 0.009, 0.011,
        0.011, 0.011, 0.011, 0.012, 0.012, 0.013, 0.014, 0.015, 0.016, 0.017,
        0.021, 0.023, 0.022, 0.026, 0.027, 0.030, 0.033, 0.041
    };   

    a1_data.push_back({"SLAC E143", x, q2, kBlue, 33});
}

// SLAC E155, p, g1/F1, 2000, https://pdf.sciencedirectassets.com/271623/1-s2.0-S0370269300X03036/1-s2.0-S0370269300010145/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEFoaCXVzLWVhc3QtMSJHMEUCIBdIQv1wWtTWTgLYfSQNdYxncPH6d37m5j%2F5K73Jr5dFAiEA%2Brs3pvk4WwjHUUD6kTJStxFEW8zQig%2B92YTKVWF%2BRbgqsgUIIxAFGgwwNTkwMDM1NDY4NjUiDJzgqkoPiuvxxrD%2BZCqPBcNBeOHA4TBPDQqUQt5XbJuARt4NIK3sWfPYqzK5%2FjveN9X%2Fllnj8Gtk6jPjcJKx3q2xpIBbmhzQPPrZGZYA3x8HBn8MJTei1tldTWBED43ReM%2Bd1MYrb0wm%2F2H6HsPCFq6f5%2FJSO4L%2FC40hacRGZQyaviPCn3gDtDOZRJviUqExNZApzczpYiaWTcKAz96w5%2B1T8Co%2BTgx5YDML%2FM9kanKi3sDryrBRWv1Ll%2Bp4cnHtYOzQBI6Dq8iipxEew1lb%2BQOJjTBQqspW7TuSOeMmJ4YRtkiGJBWHXVo%2BEYhSWGEhGMhhwQGHDx9sIUrgmu1MRfvgY4%2BSo0Ps6F3Veaa7X%2FMzzY7BnsJeLQzouztHCsdjzp%2BT0yXE97Vy2SArGowkcZ%2FGcM5JSILwpLcaLygF225vrK1BMeBM8mE04Y7e3V1JzcZ9CJXXjObvhNaOaJfO%2FsE9Y5D7DuxEbo3jIIQB1BIw2jXkT0NbXWaXfT57rhWwKvHtTofs1hJ79R1LUSrqwjzVfCItt6lxMyDfMTUualOJeXLaK3MbSuEN0v8AzKPemqS6OlRlpsuwanfTQtjyaIkIGY2lrc83Ww6dSgn%2FEUXNBLObU%2BGKPKhvMJBX6a8ZR6nC6cb0TF07T5jJ2vJdjA%2FOYP%2BbfHq0yx%2BX%2Bvyfs5GJ8Lxj4bLyT0EtegnG8wqQ2fnmkZ1AYCAbvjV1tRT91VUZwM7Iv1V7ngWxwzVltuPiAO4rjtWOB9Y%2FZh7pMKMo%2BuyA8wZ2iPgB69%2BEMnc3imn8ErARbMfsdwT3%2BPQeO7iD%2B%2BEHgj3fEjlidvPkfcFJB%2Br2ccFAL6J8ocdIeA1bgLDAdEhfDtL2uT5sd7dx5OdCt7z0K0vZQRpBEK81WFsw5NCZzwY6sQEJC3Vqj%2B7%2BKjAuekv2Gl7gnzfeEsnEJs2aCYDuSMbDBZxy0uVkAPC05J0n93T68HLrMW%2BnDcNlqe3BEHs47SNf17gtvbKu9N%2BZQh0vFKefwNe7QiXs%2FnDPyJHG7UzVFWbgB5jcmJba13TUvuuPXwG8tD4h1zEQZW7XKsh2iT61zK5F0Xdym1hu5uegxDHwXZdpcVOvrOyt1kGDk%2BNoikefQE3qJF%2BoLkjte5BTM%2BWOOpo%3D&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20260420T191617Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTY3NAZWWPV%2F20260420%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=62799506d0195d28d2afe851fd8e1b8e39ffcda2a604c93bc848ed8e0ac0f5c9&hash=b5a1b3cd6ff0dedde972fdfcc7ab599fd0a72051b45fb79b4d8b2faa9610e053&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S0370269300010145&tid=spdf-3ec56b1b-426f-41c4-a546-13157dfdbeaa&sid=3223362b49a54048ea6889c7bb7a31e04f95gxrqa&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&rh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=10145c0b045656065652&rr=9ef67b85f88db7a0&cc=us
// g2 is available from E155X
void slac_e155()
{
    std::vector<double> x = {
        0.015, 0.025, 0.035, 0.050, 0.050, 0.080, 0.080, 0.125, 0.125, 0.125,
        0.175, 0.175, 0.175, 0.250, 0.250, 0.250, 0.350, 0.350, 0.350,
        0.500, 0.500, 0.500, 0.750, 0.750
    };
    std::vector<double> q2 = {
         1.22,  1.59,  2.05,  2.58,  4.01,  3.24,  5.36,  4.03,  7.17, 10.99,
         4.62,  8.90, 13.19,  5.06, 10.64, 17.21,  5.51, 12.60, 22.73,
         5.77, 14.02, 26.86, 15.70, 34.72
    };
 
    // Proton g1^p/F1^p
    std::vector<double> g1_over_F1 = {
        0.048, 0.057, 0.070, 0.111, 0.222, 0.155, 0.150, 0.186, 0.209, 0.307,
        0.273, 0.247, 0.305, 0.358, 0.353, 0.396, 0.424, 0.466, 0.500,
        0.561, 0.561, 0.507, 0.622, 0.559
    };
    std::vector<double> g1_over_F1_stat_err = {
        0.009, 0.008, 0.008, 0.009, 0.088, 0.009, 0.011, 0.012, 0.007, 0.051,
        0.023, 0.012, 0.022, 0.023, 0.011, 0.014, 0.049, 0.020, 0.025,
        0.058, 0.024, 0.042, 0.091, 0.405
    };
    std::vector<double> g1_over_F1_syst_err = {
        0.004, 0.006, 0.007, 0.009, 0.009, 0.013, 0.013, 0.018, 0.018, 0.018,
        0.023, 0.023, 0.023, 0.030, 0.030, 0.030, 0.039, 0.039, 0.038,
        0.048, 0.048, 0.048, 0.050, 0.050
    };

    a1_data.push_back({"SLAC E155", x, q2, kRed + 1, 33});
}

void slac_e155x()
{
    std::vector<double> x  = {0.021, 0.026, 0.038, 0.061, 0.098, 0.155, 0.245, 0.380, 0.580, 0.780};
    std::vector<double> q2 = {0.80,  0.90,  1.10,  1.40,  2.30,  3.70,  5.00,  7.10,  8.40,  8.20};

    a1_data.push_back({"SLAC E155x", x, q2, kOrange + 1, 26});
}


// DESY HERMES, p, 1998, https://arxiv.org/pdf/hep-ex/9807015
void desy_hermes()
{
    std::vector<double> x = {
        0.023, 0.028, 0.033, 0.040, 0.047, 0.056, 0.067, 0.080, 0.095, 0.114,
        0.136, 0.162, 0.193, 0.230, 0.274, 0.327, 0.389, 0.464, 0.550, 0.660
    };
    std::vector<double> q2 = {
        0.92, 1.01, 1.11, 1.24, 1.39, 1.56, 1.73, 1.90, 2.09, 2.26,
        2.44, 2.63, 2.81, 3.02, 3.35, 3.76, 4.25, 4.80, 5.51, 7.36
    };
 
    // g1^p/F1^p
    std::vector<double> g1_over_F1 = {
        0.064, 0.080, 0.069, 0.090, 0.112, 0.122, 0.103, 0.163, 0.163, 0.198,
        0.249, 0.244, 0.307, 0.336, 0.354, 0.468, 0.494, 0.520, 0.784, 0.615
    };
    std::vector<double> g1_over_F1_stat_err = {
        0.013, 0.012, 0.011, 0.011, 0.011, 0.012, 0.012, 0.013, 0.014, 0.015,
        0.016, 0.017, 0.019, 0.022, 0.025, 0.031, 0.042, 0.053, 0.075, 0.108
    };
    std::vector<double> g1_over_F1_syst_err = {
        0.004, 0.006, 0.005, 0.006, 0.007, 0.008, 0.007, 0.011, 0.011, 0.013,
        0.016, 0.016, 0.020, 0.022, 0.023, 0.031, 0.032, 0.034, 0.051, 0.040
    };
    std::vector<double> g1_over_F1_R_err = {
        0.002, 0.003, 0.003, 0.004, 0.004, 0.004, 0.004, 0.006, 0.006, 0.007,
        0.009, 0.008, 0.010, 0.010, 0.011, 0.014, 0.015, 0.016, 0.024, 0.020
    };

    a1_data.push_back({"DESY HERMES", x, q2, kMagenta + 1, 29});
}

// JLab EG-CLAS, p, A_LL, 2003, https://journals.aps.org/prl/supplemental/10.1103/PhysRevLett.120.062501/EG4_ND3_tables.pdf

// JLab RSS, d, A1, A2, 2007 https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.98.132003
void jlab_rss()
{
    // Figure 2 (approx): A1 (open triangles) and A2 (filled circles) vs W
    std::vector<double> W = {
        1.10, 1.13, 1.16, 1.20, 1.23, 1.27, 1.30, 1.33, 1.36,
        1.39, 1.42, 1.45, 1.48, 1.51, 1.54, 1.57, 1.60, 1.63,
        1.66, 1.69, 1.72, 1.75, 1.78, 1.81, 1.84, 1.87, 1.90
    };

    std::vector<double> x;
    std::vector<double> q2;
    for ( int i = 0; i < W.size(); i ++ )
    {
        q2.push_back(1.3);
        x.push_back(q2[i] / (W[i]*W[i] - MASS_PROTON*MASS_PROTON + q2[i]));
        // printf("%f %f\n", x[i], ALL[i]);
    }

    std::vector<double> A1 = {
        0.18,  0.02, -0.30, -0.22, -0.12,  0.05,  0.22,  0.40,  0.52,
        0.43,  0.45,  0.53,  0.70,  0.73,  0.72,  0.68,  0.62,  0.54,
        0.50,  0.51,  0.58,  0.56,  0.50,  0.41,  0.32,  0.33,  0.36
    };

    std::vector<double> A2 = {
        -0.10, -0.02, -0.04, -0.10, -0.05,  0.03,  0.10,  0.26,  0.18,
        0.16,  0.19,  0.24,  0.30,  0.24,  0.18,  0.17,  0.19,  0.20,
        0.30,  0.26,  0.23,  0.18,  0.14,  0.13,  0.10,  0.09,  0.20
    };

    // Approximate statistical errors read from vertical bars in Fig. 2
    std::vector<double> A1_stat_err = {
        0.20, 0.14, 0.02, 0.03, 0.03, 0.03, 0.04, 0.04, 0.04,
        0.04, 0.04, 0.04, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03,
        0.03, 0.03, 0.03, 0.03, 0.03, 0.04, 0.05, 0.05, 0.06
    };

    std::vector<double> A2_stat_err = {
        0.25, 0.14, 0.05, 0.04, 0.04, 0.04, 0.04, 0.05, 0.05,
        0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
        0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.12
    };

    // Approximate systematic errors from Fig. 2 error bands
    // (upper band for A1, lower band for A2 per caption)
    std::vector<double> A1_syst_err = {
        0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
        0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
        0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05
    };

    std::vector<double> A2_syst_err = {
        0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03,
        0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03,
        0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03
    };

    a1_data.push_back({"JLab RSS", x, q2, kCyan + 2, 47});
}

// CERN COMPASS, p, A_LL, 2016 https://arxiv.org/pdf/1503.08935, https://arxiv.org/pdf/1001.4654
void cern_compass()
{
    std::vector<double> x = {
        0.0035, 0.0036, 0.0038, 0.0044, 0.0045, 0.0046, 0.0055, 0.0055, 0.0056,
        0.0069, 0.0069, 0.0071, 0.0089, 0.0089, 0.0090, 0.0116, 0.0117, 0.0120,
        0.0164, 0.0165, 0.0168, 0.0239, 0.0240, 0.0246, 0.0341, 0.0343, 0.0347,
        0.0473, 0.0480, 0.0492, 0.0740, 0.0754, 0.0800, 0.1190, 0.1210, 0.1250,
        0.1710, 0.1720, 0.1750, 0.2220, 0.2220, 0.2240, 0.2890, 0.2900, 0.2960,
        0.4030, 0.4050, 0.4130, 0.5610, 0.5670, 0.5750,
        0.0046, 0.0055, 0.0070, 0.0090, 0.0147, 0.0247, 0.0346, 0.0487, 0.0765,
        0.1220, 0.1720, 0.2220, 0.2900, 0.4050, 0.5680
    };

    std::vector<double> q2 = {
        1.03, 1.10, 1.22, 1.07, 1.24, 1.44, 1.11, 1.36, 1.68,
        1.14, 1.50, 2.02, 1.17, 1.62, 2.41, 1.21, 1.75, 2.92,
        1.26, 1.92, 3.74, 1.55, 2.49, 5.16, 2.18, 3.50, 7.07,
        2.65, 5.00, 10.4, 4.91, 10.7, 19.7, 8.23, 17.8, 31.7,
        12.9, 26.9, 43.8, 16.1, 32.1, 52.4, 21.7, 42.1, 66.3,
        28.4, 53.1, 85.1, 29.8, 50.4, 96.1,
        1.10, 1.20, 1.37, 1.59, 2.14, 3.24, 4.36, 6.05, 9.42,
        14.9, 20.9, 26.7, 34.6, 47.1, 62.1
    };

    const std::vector<double> A1p = {
        0.059, -0.004,  0.002,  0.006,  0.021,  0.023,  0.009,  0.026,  0.022,
        0.033,  0.041,  0.006,  0.007,  0.029,  0.015,  0.044,  0.040,  0.044,
        0.087,  0.100,  0.063,  0.072,  0.079,  0.079,  0.103,  0.099,  0.083,
        0.128,  0.136,  0.103,  0.147,  0.203,  0.129,  0.291,  0.263,  0.243,
        0.299,  0.316,  0.344,  0.405,  0.340,  0.268,  0.397,  0.374,  0.392,
        0.396,  0.400,  0.631,  0.420,  0.750,  0.870,
        0.006, 0.019, 0.035, 0.033, 0.047, 0.076, 0.115, 0.130, 0.172,
        0.218, 0.286, 0.446, 0.453, 0.594, 0.855
    };

    const std::vector<double> A1p_stat_err = {
        0.029, 0.027, 0.032, 0.021, 0.020, 0.022, 0.024, 0.020, 0.020,
        0.020, 0.015, 0.014, 0.027, 0.018, 0.014, 0.026, 0.017, 0.011,
        0.034, 0.020, 0.011, 0.030, 0.025, 0.011, 0.035, 0.041, 0.015,
        0.040, 0.026, 0.015, 0.031, 0.020, 0.023, 0.038, 0.028, 0.034,
        0.045, 0.045, 0.050, 0.060, 0.066, 0.060, 0.057, 0.077, 0.062,
        0.086, 0.120, 0.088, 0.230, 0.230, 0.150,
        0.017, 0.013, 0.009, 0.010, 0.006, 0.009, 0.012, 0.012, 0.013,
        0.017, 0.024, 0.032, 0.033, 0.049, 0.094
    };

    const std::vector<double> A1p_syst_err = {
        0.014, 0.012, 0.012, 0.008, 0.008, 0.011, 0.011, 0.008, 0.008,
        0.009, 0.007, 0.007, 0.013, 0.007, 0.006, 0.013, 0.011, 0.005,
        0.015, 0.011, 0.006, 0.016, 0.011, 0.008, 0.016, 0.018, 0.013,
        0.023, 0.016, 0.011, 0.016, 0.017, 0.021, 0.024, 0.021, 0.021,
        0.027, 0.036, 0.029, 0.043, 0.035, 0.045, 0.035, 0.050, 0.035,
        0.051, 0.060, 0.054, 0.100, 0.120, 0.090,
        0.008, 0.006, 0.005, 0.005, 0.004, 0.006, 0.009, 0.010, 0.012,
        0.015, 0.020, 0.030, 0.032, 0.043, 0.068
    };

    a1_data.push_back({"CERN COMPASS", x, q2, kCyan + 2, 24});

    std::vector<double> A1p_err;
    for (size_t i = 0; i < A1p_stat_err.size(); ++i)
        A1p_err.push_back(sqrt(A1p_stat_err[i]*A1p_stat_err[i] + A1p_syst_err[i]*A1p_syst_err[i]));
    a1_value.push_back({"CERN COMPASS", x, q2, A1p, A1p_err, kCyan + 2, 24});
}

// CERN COMPASS, p, A_LL, 2015 https://pdf.sciencedirectassets.com/271623/1-s2.0-S0370269318X00057/1-s2.0-S0370269318302405/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEHUaCXVzLWVhc3QtMSJIMEYCIQD5si1Ok6ue%2B7pmw9PIq2824o1xz%2FYse4%2FFfoKNXW%2BE%2FQIhANvn3jIrrfWHXg48TaTUhI%2FJ44YbfRA91JzUIgJHAxcKKrMFCD4QBRoMMDU5MDAzNTQ2ODY1IgxfxZTK1DHPvbKczIkqkAX6xt8VKj7Gu9hUh1ImKLWoKjyLn0hdAeDBlxYnlYBSv54umZxIayMtAuTi2YXLfbKdTCxUDyEQmYEMCPi2waDMf4imEIXokh7ZBKpaOlEMfoLHavfYBuC%2BYuJhQzGgzkIWVVb8Hu0tAhLScWq%2FINgkxeRE%2Fv%2BSpjX665anX2yP3JHCwvQD1fsRQQR7rH%2FYZFbQ9ctvkd9xDhLxV6yNY%2FWrdbZJ1ISCPtrQifMcr5HZnWVEqcltP3HyiQu6SUGvm%2FlQSB1v3My5Pa6dNI%2BwXysYoIHJ%2Fj970UW%2B5DaZdCjyim20WcR95%2BmwqcspBts8Eqp2mMBrjzlKHRhhmPNhvHwNyQNX1edEwSTNDEeTudMuxCB738oJnKTKtY%2FfTPFuvc%2FU3Z24IpLNJHHe%2FFCsBx7ykIdhvYGufEr6kaxzXb76kVobjEq%2F6VhcB52F5JNeoBLipi1h54ps4UbuRdS5aRhIFACECs7PiQVLr4sv9hiibSSc6Ajz0pQJC2n5qQltoPgNLDIZDwTEd6Df94jyyZqTjgmaIIVcVt6vyk3S%2FbvMRGKH9PhTrAvv6wZwHgZ7I7F%2BVHXDFDHLTL0mUe0fRQjHyCQYOMqEoU%2BhnxOPOq5rXR%2B3C%2BqaVmBQAVhQOciG9pZ%2BGonlBghuCXXE%2FRNaUVLIWLRvKi2nt0%2BN0MTedcjfvEF%2F6gsmc1mGFQ1fohug1Tes%2FDZgg3jOEBi4CbXWbFo1RwyBnVnh%2BiMB4cjYS1GBfiarq0kApf8ljvLTLIEn7tD%2F3PR%2BbR7jzK9%2Fo68rLM20%2FNONW5hR4ePEzeKxUsFyNOWDgZzgZ2INI40cb63E8bzgBLdMCEQipQn8BjRGj%2BgBUIcgK2FjwK1fF1lzhpXTCjCky5%2FPBjqwAVvpRkyiw5QnV5VszW7EyebEzlC1K%2FmnS1QUL%2BQgxiBMUpang2etSv9ohhBaN39QZSyqHxpAIH1CUGrUHlM6uN751IiNxuZw%2BTsEG%2FfWUQuBH38r6gxHjnoDm7OTlM9SGEx5WMD6kWMQzUrtWjqH0dgcykSOyEUS3JGEZJo2tJgJ5qfvbyyttDEfxtuTQZiM7QYu6qt1AIlKxOZsZxXtqTIz5dfjy1MQqlaPNxlB3cyH&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20260421T214433Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTYSV245EZU%2F20260421%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=304dc5c06a4503ff510122d9db6a71bd9995db5df6f411569065ea189418d40c&hash=7527defbba4306593c863ac51db3e5b29de8f9d1de0202014bdca34d0a4ef62e&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S0370269318302405&tid=spdf-7996abe6-d99d-407a-8440-281e7f0e79bb&sid=8dda62a319065340ca688ca-e362b7849f14gxrqa&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&rh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=0f165d0607085d5551&rr=9eff92175af22dae&cc=us
void cern_compass_lowq()
{
    const std::vector<double> x = {
        0.000052, 0.000081, 0.00013, 0.00020, 0.00032,
        0.00050, 0.00079, 0.0013, 0.0020, 0.0031,
        0.0049, 0.0077, 0.012, 0.019, 0.028
    };

    const std::vector<double> q2 = {
        0.0076, 0.013, 0.022, 0.037, 0.059,
        0.094, 0.15, 0.24, 0.38, 0.57,
        0.67, 0.71, 0.76, 0.82, 0.91
    };

    const std::vector<double> a1 = {
        0.0073, 0.0064, 0.0079, 0.0065, 0.0052,
        0.0117, 0.0147, 0.0171, 0.0073, 0.0167,
        0.0273, 0.0386, 0.0540, 0.0330, 0.0200
    };
        
    a1_data.push_back({"CERN COMPASS low Q^{2}", x, q2, kBlue, 20});
}

// JLab SANE, p, A1, A2


void load_world_a1p_data()
{
    jlab_rss();

    cern_compass();
    cern_compass_lowq();
    cern_smc();
    cern_emc();

    desy_hermes();

    slac_e155();
    slac_e143();
    slac_e130();
    slac_e80();
}

void load_world_a2p_data()
{
    jlab_rss();

    slac_e155();
    slac_e143();
}