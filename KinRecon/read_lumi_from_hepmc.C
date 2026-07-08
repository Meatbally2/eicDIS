
void read_lumi_from_hepmc()
{
    std::string address = "root://dtn-eic.jlab.org/";

    // std::string path = "volatile/eic/EPIC/EVGEN/SIDIS/pythia6-eic/1.0.0/10x100/q2_0to1";
    // std::string type = "pythia_ep_noradcor_10x100_q2_0.000000001_1.0";

    std::string path = "volatile/eic/EPIC/EVGEN/SIDIS/pythia6-eic/1.2.0/ep_noradcor/10x100/q2_100to1000";
    std::string type = "pythia6-eic_1.2.0_ep_noradcor_10x100_q2_100_1000";

    std::string extention = "ab.hepmc3.tree.root";


    double avg_cs = 0;
    double total_lumi = 0;

    for ( int i = 0; i < 20; i ++ )
    {
        std::string file_name = Form("%s/%s/%s_run%03d.%s", address.c_str(), path.c_str(), type.c_str(), i+1, extention.c_str());
        // std::string file_name = "root://dtn-eic.jlab.org:1094//volatile/eic/EPIC/EVGEN/SIDIS/pythia6-eic/1.0.0/10x100/q2_0to1/pythia_ep_noradcor_10x100_q2_0.000000001_1.0_run1.ab.hepmc3.tree.root";
        TFile* hepmc3_file = TFile::Open(Form("%s", file_name.c_str()));
        if ( hepmc3_file == NULL )
        {
            printf("%s not found\n", file_name.c_str());
            continue;
        }

        TTreeReader hepmc3_reader("hepmc3_tree", hepmc3_file);
        TTreeReaderArray<string> att_string(hepmc3_reader, "hepmc3_event.attribute_string");
        TTreeReaderArray<int> att_id(hepmc3_reader, "hepmc3_event.attribute_id");
        TTreeReaderArray<string> att_name(hepmc3_reader, "hepmc3_event.attribute_name");

        Int_t nEvents = hepmc3_reader.GetEntries();

        double run_cs = 0;
        hepmc3_reader.SetEntry(nEvents-1);
        for ( int i = 0; i < att_string.GetSize(); i++ )
            if ( att_id[i] == 0 && att_name[i] == "GenCrossSection")
                run_cs = std::stod(att_string[i]); 
                avg_cs += run_cs;

        std::cout << "gen cs: " << run_cs << " fb" << std::endl;

        double run_lumi = nEvents / run_cs; 
        std::cout << "gen lumi: " << run_lumi << " fb^-1" << std::endl;
        total_lumi += run_lumi;
    }

    avg_cs /= 20.;
    std::cout << "average gen cs: " << avg_cs << " fb" << std::endl;
    std::cout << "total gen lumi: " << total_lumi << " fb^-1" << std::endl;
}