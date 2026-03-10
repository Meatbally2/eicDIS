#pragma once

enum{EHE3,EP,PI_BG,BEAM_BG,EP_PYTHIA6,EP_DVMP};

std::vector<double> get_lumi(int beam_type, int Ee, int Eh)
{
    std::vector<double> lumi;

    if ( beam_type == EHE3 )
    {
        if ( Ee == 10 && Eh == 166 )
        {
            double cs[3][2] = {{0.198440424611563, 0.205327493968226}, {4.04371412707044E-02, 4.41976212963417E-02}, {1.36416909784756E-03, 1.69583242740138E-03}};
            double ev[3][2] = {{333675, 666325}, {333694, 666306}, {333365, 666640}};

            for ( int i = 0; i < 3; i ++ )
            {
                double gen_lumi = ev[i][0]/(cs[i][0]*(1e-34/1e-43)) + ev[i][1]/(cs[i][1]*(1e-34/1e-43));
                lumi.push_back(gen_lumi);
            }
        }

        if ( Ee == 18 && Eh == 110 )
        {
            double cs[3][2] = {{0.205033149243416, 0.21151100919885}, {4.27967567829114E-02, 4.63590143556042E-02}, {1.5447864536253E-03, 1.8963099176877E-03}};
            double ev[3][2] = {{333574, 666426}, {333199, 333199}, {332978, 667024}};

            for ( int i = 0; i < 3; i ++ )
            {
                double gen_lumi = ev[i][0]/(cs[i][0]*(1e-34/1e-43)) + ev[i][1]/(cs[i][1]*(1e-34/1e-43));
                lumi.push_back(gen_lumi);
            }
        }

        if ( Ee == 10 && Eh == 110 )
        {
            double cs[3][2] = {{0.184023165609017, 0.192267385649895}, {3.5058173624293300E-02, 3.92002645972039E-02}, {9.9564457630334600E-04, 1.27882734303755E-03}};
            double ev[3][2] = {{332731, 667269}, {333163, 666837}, {333358, 666645}};

            for ( int i = 0; i < 3; i ++ )
            {
                double gen_lumi = ev[i][0]/(cs[i][0]*(1e-34/1e-43)) + ev[i][1]/(cs[i][1]*(1e-34/1e-43));
                lumi.push_back(gen_lumi);
            }
        }

        if ( Ee == 5 && Eh == 41 )
        {
            double cs[3][2] = {{0.128613950877546, 0.144569529424099}, {1.55974793402807E-02, 1.91516550407188E-02}, {7.69468056138372E-04, 1706.45007082295}};
            double ev[3][2] = {{333382, 666817}, {333446, 666695}, {333837, 666194}};

            for ( int i = 0; i < 3; i ++ )
            {
                double gen_lumi = ev[i][0]/(cs[i][0]*(1e-34/1e-43)) + ev[i][1]/(cs[i][1]*(1e-34/1e-43));
                lumi.push_back(gen_lumi);
            }
        }
    }

    if ( beam_type == EP )
    {
        if ( Ee == 18 && Eh == 275 )
        {
            double gen_lumi[4] = {6.73335E-03, 7.21215E-02, 1.48384E+00, 6.29541E+01}; // fb^-1
            for ( int i = 0; i < 4; i ++ )
                lumi.push_back(gen_lumi[i]);
        }
        if ( Ee == 10 && Eh == 100 )
        {
            double gen_lumi[4] = {8.99265E-03, 1.25042E-01, 3.72397E+00, 7.32992E+02}; // fb^-1
            for ( int i = 0; i < 4; i ++ )
                lumi.push_back(gen_lumi[i]);
        }
        if ( Ee == 5 && Eh == 41 )
        {
            double gen_lumi[3] = {1.24919E-02, 2.45752E-01, 1.74313E+01}; // fb^-1
            for ( int i = 0; i < 3; i ++ )
                lumi.push_back(gen_lumi[i]);
        }
    }

    if ( beam_type == EP_PYTHIA6 )
    {
         if ( Ee == 10 && Eh == 130 )
        {
            // double cs[4] = {0.63928175476836346, 4.7191206939697819E-002, 1.5146539778048273E-003, 45.369017792412620};
            double cs[4] = {0.63928175476836346, 4.7191206939697819E-002, 1.5146539778048273E-003, 8.1678590974353973E-006}; 

            for ( int i = 0; i < 4; i ++ )
            {
                // double gen_lumi = 500000./(cs[i]*(1e-34/1e-43));
                // if ( i == 3)
                //     gen_lumi = 50000./(cs[i]*(1e-34/1e-36));
                //     // gen_lumi = 50./(cs[i]*(1e-34/1e-43));
                // lumi.push_back(gen_lumi);

                cs[i] = cs[i]*1e-6/1e-15; // convert from microbarn to fb
                double gen_lumi = 500000./cs[i];
                if ( i == 3)
                    // gen_lumi = 50./(cs[i]);
                    gen_lumi = 50000./cs[i];
                    // gen_lumi = 50000./(cs[i]*1e-6);
                lumi.push_back(gen_lumi);

                std::cout << "Process " << i << ", cross section: " << cs[i] << " fb, gen lumi: " << gen_lumi << " fb^-1" << std::endl;
            }
        }

        if ( Ee == 10 && Eh == 250 )
        {
            // double cs[4] = {0.70873205053123534, 5.9758065033598498E-002, 2.2867629538197093E-003, 22.527455461166682};
            double cs[4] = {0.70873205053123534, 5.9758065033598498E-002, 2.2867629538197093E-003, 2.7215279714732617E-005};

            for ( int i = 0; i < 4; i ++ )
            {
                // double gen_lumi = 500000./(cs[i]*(1e-34/1e-43));
                // if ( i == 3)
                //     gen_lumi = 50000./(cs[i]*(1e-34/1e-37));
                    // gen_lumi = 50./(cs[i]*(1e-34/1e-43));

                cs[i] = cs[i]*1e-6/1e-15; // convert from microbarn to fb
                double gen_lumi = 500000./cs[i];
                if ( i == 3)
                    // gen_lumi = 50./cs[i];
                     gen_lumi = 50000./cs[i];
                    // gen_lumi = 50000./(cs[i]*1e-6);
                lumi.push_back(gen_lumi);
                
                std::cout << "Process " << i << ", cross section: " << cs[i] << " fb, gen lumi: " << gen_lumi << " fb^-1" << std::endl;
            }
        }
    }

    if ( beam_type == PI_BG )
    {
        double gen_lumi = 9.26e-4; // fb^-1
        lumi.push_back(gen_lumi);
    }

    return lumi;
}

