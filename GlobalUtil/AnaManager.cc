#include "AnaManager.hh"

// fixed based on simulation campaign
const std::string n_group[3] = {"1to10", "10to100", "100to1000"};
const std::string p_group[4] = {"minQ2=1", "minQ2=10", "minQ2=100", "minQ2=1000"};

const std::string address =  "root://dtn-rucio.jlab.org:1094//volatile/eic/EPIC";

AnaManager::AnaManager(std::string ana_name_) : ana_name(ana_name_) {
}

AnaManager::~AnaManager() {
}

void AnaManager::Initialize(bool is_select_region_, int region_index_, int starting_file, bool is_analyse_protons_) 
{
    is_select_region = is_select_region_;
    region_index = region_index_;
    starting_file = starting_file_;
    is_analyse_protons = is_analyse_protons_;

    return;
}

std::string AnaManager::GetOutputName()
{
    std::string outname;

    if ( is_analyse_protons )
        outname = is_select_region ? Form("18x275_%s_%s.root", p_group[region_index].c_str(), ana_name.c_str()) : Form("18x275_%s.root", ana_name.c_str());
    else
        outname = is_select_region ? Form("10x166_%s_%s.root", n_group[region_index].c_str(), ana_name.c_str()) : Form("10x166_%s.root", ana_name.c_str());

    return outname;
}

vector<std::string> AnaManager::GetInputNames()
{
    std::vector<std::string> inFiles;

    int n_set = is_analyse_protons ? 4 : 3;
    
    int total_file = 0;
    
    for ( int r = 0; r < n_set; r ++ )
    {
        if ( is_select_region )
            if ( r != region_index )
                continue;
        
        std::string file_name;
        if ( is_analyse_protons )
            file_name = Form("../data/ep_25_10_0/18x275minQ2=%.0f_filelist.txt", pow(10,r));
        else
            file_name = Form("../data/en_25_08_0/10x166minQ2=%.0f_filelist.txt", pow(10,r));

        std::ifstream data_file(file_name);

        int line_c = 0;
        std::string line;
        while ( getline(data_file, line) )
        {
            if ( starting_file < 0 )
                if ( line_c >= 10 )
                    break;

            if ( starting_file >= 0 )
                if ( line_c < starting_file )
                {
                    line_c ++;
                    continue;
                } 
                else if ( line_c >= starting_file + 2000 )
                    break;

            std::string fname;
            std::stringstream ss(line);
            ss >> fname;
            fname.erase(0, 5);
            inFiles.push_back(address+fname);

            line_c ++;
            total_file ++;
        }

        data_file.close();
    }

    cout << "total of " << total_file << " files are found" << endl << endl;

    return inFiles;
}