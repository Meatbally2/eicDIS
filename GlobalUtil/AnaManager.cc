#include "AnaManager.hh"

// fixed based on simulation campaign
const std::string n_group[3] = {"1to10", "10to100", "100to1000"};
const std::string p_group[4] = {"minQ2=1", "minQ2=10", "minQ2=100", "minQ2=1000"};

// const std::string address =  "root://dtn-rucio.jlab.org:1094//volatile/eic/EPIC";
const std::string address =  "root://dtn-rucio.jlab.org:1094/";

const int step = 2000; // number of files to process in one batch

AnaManager::AnaManager(std::string ana_name_) : ana_name(ana_name_) {
}

AnaManager::~AnaManager() {
}

void AnaManager::Initialize(bool is_select_region_, int region_index_, int starting_file_index_, bool is_analyse_protons_) 
{
    is_select_region = is_select_region_;
    region_index = region_index_;
    file_index = starting_file_index_;
    starting_file = starting_file_index_*step;
    is_analyse_protons = is_analyse_protons_;

    return;
}

void AnaManager::InitializeForLocal(std::string type_) 
{
    file_type = type_;
    is_analyse_protons = true;
    is_select_region = false;
    region_index = -1;
    starting_file = -1; 

    return;
}

std::string AnaManager::GetOutputName()
{
    std::string outname;

    if ( starting_file >= 0 )
        ana_name += Form("_f%d", file_index);

    if ( is_analyse_protons )
        outname = is_select_region ? Form("tmp/18x275_%s_%s.root", p_group[region_index].c_str(), ana_name.c_str()) : Form("tmp/18x275_%s.root", ana_name.c_str());
    else
        outname = is_select_region ? Form("tmp/10x166_%s_%s.root", n_group[region_index].c_str(), ana_name.c_str()) : Form("tmp/10x166_%s.root", ana_name.c_str());

    return outname;
}

vector<std::string> AnaManager::GetLocalInputNames()
{
    std::vector<std::string> inFiles;

    for ( int r = 0; r < 10; r ++ )
    {
        std::string file_name = Form("../data/BG_Study/18x275_%s/eicrecon_%d_to_%d.root", file_type.c_str(), r*100, r*100+99);
        std::cout << "File " << r << " : " << file_name << std::endl;
        inFiles.push_back(file_name);
    }     

    return inFiles;
}

vector<std::string> AnaManager::GetLowQInputNames()
{
    std::vector<std::string> inFiles;

    for ( int r = 0; r <= 1709; r ++ )
    {
        // if ( starting_file < 0 )
        //     if ( r > 10 )
        //         break;

        if ( r == 1519 )
            continue;

        std::string fname = Form("epic:/RECO/25.05.0/epic_craterlake/SIDIS/pythia6-eic/1.0.0/18x275/q2_0to1/pythia_ep_noradcor_18x275_q2_0.000000001_1.0_run9.ab.%04d.eicrecon.edm4eic.root", r);
        fname.erase(0, 5);
        std::cout << "File " << r << " : " << fname << std::endl;
        inFiles.push_back(address+fname);
    }

    return inFiles;
}

vector<std::string> AnaManager::GetInputNames()
{
    std::vector<std::string> inFiles;

    int n_set = is_analyse_protons ? 4 : 3;
    
    int total_file = 0;

    std::string campaign = is_analyse_protons ? "25.10.0" : "25.10.2";

    std::string file_name = "../data/" + campaign + "_manifest.txt";
    // std::string file_name = "../data/test.txt";
    
    std::string prefix = "/volatile/eic/EPIC/RECO/" + campaign + "/epic_craterlake/DIS/";
    // std::string phys_group = "/epic_craterlake/DIS/";

    for ( int r = 0; r < n_set; r ++ )
    {
        if ( is_select_region )
            if ( r != region_index )
                continue;

        // std::cout << " r " << r << std::endl;
        
        // std::string file_name;
        // if ( is_analyse_protons )
        //     file_name = Form("../data/ep_25_05_0/18x275minQ2=%.0f_filelist.txt", pow(10,r));
        // else
        //     file_name = Form("../data/en_25_05_0/10x166minQ2=%.0f_filelist.txt", pow(10,r));

        std::string target;
        if ( is_analyse_protons )
        {
            std::string gen_group = "NC/"; 
            std::string beam_group = "18x275/";
            std::string sample_group = Form("minQ2=%.0f",pow(10,r));
            target = gen_group + beam_group + sample_group;
        }
        else
        {
            std::string gen_group = "BeAGLE1.03.02-1.2/"; 
            std::string beam_group = "eHe3/10x166/";
            std::string sample_group;
            if ( r == 0 )
                sample_group = Form("q2_2to10/");
            else if ( r == 1 )
                sample_group = Form("q2_10to100/");
            else
                sample_group = Form("q2_100to10000/");

            target = gen_group + beam_group + sample_group;
            
            // std::cout << "sample_group: " << sample_group << std::endl;
            // std::cout << "target: " << target << std::endl;            
        }

        // std::cout << "prefix: " << prefix << std::endl;
        // std::cout << "prefix size: " << prefix.size() << " target size: " << target.size() << std::endl; 

   
        std::ifstream data_file(file_name);

        int line_c = 0;
        std::string line;
        while ( getline(data_file, line) )
        {
            std::string fname;
            std::stringstream ss(line);
            ss >> fname;
            // fname.erase(0, 5);

            int compare = line.compare(prefix.size(), target.size(), target);
            // std::cout << "checking: " << line.substr(prefix.size(), target.size()) << std::endl;
            // std::cout << "compare: " << compare << std::endl;

            if ( compare != 0 )
                continue;

            if ( starting_file < 0 )
                if ( line_c >= 3 )
                    break;

            if ( starting_file >= 0 )
            {
                if ( line_c < starting_file )
                {
                    line_c ++;
                    continue;
                } 
                else if ( line_c >= starting_file + step )
                    break;
            }
                
            // inFiles.push_back(address+fname);
            inFiles.push_back(address+line);

            // std::cout << "File " << total_file << ": " << fname << std::endl;
            //  std::cout << "File " << total_file << ": " << address+line << std::endl;

            line_c ++;
            total_file ++;
        }

        data_file.close();
        
    }

    cout << "total of " << total_file << " files are found" << endl << endl;

    // for ( const auto& line : inFiles )
    //     std::cout << line << std::endl;

    return inFiles;
}