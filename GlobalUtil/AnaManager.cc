#include "AnaManager.hh"

// fixed based on simulation campaign
const std::string n_group[3] = {"1to10", "10to100", "100to1000"};
const std::string p_group[4] = {"minQ2=1", "minQ2=10", "minQ2=100", "minQ2=1000"};

const std::vector<std::string> kDefaultXrdHosts = {
    "root://epicxrd1.sdcc.bnl.gov:1095",
    "root://xcache.jlab.org:1095",
    "root://se01.af.uchicago.edu:1095"
};

const int step = 200; // number of files to process in one batch
// const int step = 2;

namespace {
std::string BuildTargetForRegion(int beam_type, int Ee, int Eh, int r);
std::vector<std::string> ResolveRemoteFilesForTarget(const std::string& campaign, int Ee, int Eh, int beam_type, const std::string& target, std::string* selectedHost, std::string* selectedPrefix);

} // namespace

AnaManager::AnaManager(std::string ana_name_) : ana_name(ana_name_) {
}

AnaManager::~AnaManager() {
}

void AnaManager::UpdateCampaign()
{
    campaign = "26.07.1";

    if (beam_type == EHE3 && Ee == 9 && Eh == 166)
        campaign = "26.07.2";

    if (beam_type == PI_BG)
        campaign = "26.03.0";
}

void AnaManager::Initialize(bool is_select_region_, int region_index_, int starting_file_index_, int beam_type_) 
{
    is_select_region = is_select_region_;
    region_index = region_index_;
    file_index = starting_file_index_;
    starting_file = starting_file_index_*step;
    // starting_file = 1000;
    beam_type = beam_type_;

    UpdateCampaign();
    return;
}

void AnaManager::SetBeamEnergy(int Ee_, int Eh_)
{
    Ee = Ee_;
    Eh = Eh_;
    UpdateCampaign();
    return;
}

void AnaManager::InitializeForLocal(std::string type_) 
{
    file_type = type_;
    beam_type = true;
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

    std::string prefix[6] = {"eHe3", "ep_pythia8", "piBG", "beamBG", "ep_pythia6", "ep"};
    outname = is_select_region ? Form("tmp/%s_%dx%d_%s_%s.root", prefix[beam_type].c_str(), Ee, Eh, p_group[region_index].c_str(), ana_name.c_str()) : Form("tmp/%s_%dx%d_%s.root", prefix[beam_type].c_str(), Ee, Eh, ana_name.c_str());

    if ( beam_type == PI_BG )
        outname = is_select_region ? Form("tmp/%s_%dx%d_%s_%s.root", prefix[beam_type].c_str(), Ee, Eh, "minQ2=0", ana_name.c_str()) : Form("tmp/%s_%dx%d_%s.root", prefix[beam_type].c_str(), Ee, Eh, ana_name.c_str());

    return outname;
}

vector<std::string> AnaManager::GetInputNames()
{
    std::cout << "** Gathering input files for analysis..." << std::endl;

    std::vector<std::string> inFiles;

    int n_set = beam_type ? 4 : 3;
    if ( beam_type == PI_BG )
        n_set = 1;
    
    int total_file = 0;
    
    std::cout << "n_set: " << n_set << std::endl;

    for ( int r = 0; r < n_set; r ++ )
    {
        if ( is_select_region )
            if ( r != region_index )
                continue;

        const std::string target = BuildTargetForRegion(beam_type, Ee, Eh, r);
        if (target.empty()) {
            std::cerr << "o.O Error: invalid beam type specified." << std::endl;
            return {};
        }

        std::cout << "searching: " << campaign << " .. " << target << std::endl;

        std::string selectedHost;
        std::string usedPrefix;
        std::vector<std::string> remoteFiles = ResolveRemoteFilesForTarget(campaign, Ee, Eh, beam_type, target, &selectedHost, &usedPrefix);

        if (remoteFiles.empty()) {
            std::cout << "No remote ROOT files found for target: " << target << " (tried both /volatile and /eic namespaces)." << std::endl;
            continue;
        }

        std::cout << "Resolved host: " << selectedHost << std::endl;
        std::cout << "Resolved prefix: " << usedPrefix << std::endl;

        int line_c = 0;
        for (const auto& file : remoteFiles)
        {
            if ( file_index == -1 )
                if ( line_c >= 1 )
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

            inFiles.push_back(file);
            line_c ++;
            total_file ++;
        }

        std::cout << "collected files in region set: " << line_c << std::endl;
    }

    cout << "total of " << total_file << " files are found" << endl;

    // Validate input files
    std::vector<std::string> validFiles = ValidateFiles(inFiles);
    if (validFiles.empty()) {
        std::cerr << "No valid input files found." << std::endl;
        return {};
    }

    // std::cout << "Number of valid input files: " << validFiles.size() << " .. passing them to PODIO" << std::endl;

    // for ( const auto& line : inFiles )
    //     std::cout << line << std::endl;

    return validFiles;
}

vector<std::string> AnaManager::ValidateFiles(std::vector<std::string>& fileNames)
{
    // Provided by Sakib R. :D

    std::cout << "** Validating input files..." << std::endl;

    std::vector<std::string> validFiles;

    int count = 0;
    for (const auto& fileName : fileNames) {

        // std::cerr << "Validating file " << count++ << ": " << fileName << std::endl;

        TFile* file = TFile::Open(fileName.c_str());
        
        if (!file || file->IsZombie()) {
            std::cerr << "!!! Skipping corrupted or inaccessible file: " << fileName << std::endl;
            if (file) delete file;
            continue;
        }

        // Some files may open but contain no keys (recovered/empty). Treat as invalid.
        TList* keys = file->GetListOfKeys();
        Long64_t nkeys = keys ? keys->GetSize() : 0;
        if (nkeys == 0) {
            std::cerr << "!!! Skipping file with no keys (empty/recovered): " << fileName << std::endl;
            delete file;
            continue;
        }
        
        if (file->TestBit(TFile::kRecovered)) {
            std::cerr << "!!! Skipping recovered file (possibly corrupted): " << fileName << std::endl;
            delete file;
            continue;
        }
        
        delete file;
        validFiles.push_back(fileName);
    }
    
    return validFiles;
}

namespace {

std::vector<std::string> ResolveRemoteFilesForTarget(const std::string& campaign, int Ee, int Eh, int beam_type, const std::string& target, std::string* selectedHost, std::string* selectedPrefix)
{
    std::vector<std::string> hosts;

    // Highest precedence: explicit ordered list from EICDIS_XRD_HOSTS.
    const char* envHosts = gSystem->Getenv("EICDIS_XRD_HOSTS");
    if (envHosts && envHosts[0] != '\0') {
        std::string token;
        std::string rawHosts(envHosts);
        for (size_t i = 0; i <= rawHosts.size(); ++i) {
            char c = (i < rawHosts.size()) ? rawHosts[i] : ',';
            if (c == ',') {
                if (!token.empty()) {
                    while (!token.empty() && token.back() == '/') {
                        token.pop_back();
                    }
                    bool exists = false;
                    for (const auto& host : hosts) {
                        if (host == token) {
                            exists = true;
                            break;
                        }
                    }
                    if (!exists) {
                        hosts.push_back(token);
                    }
                    token.clear();
                }
            } else if (c != ' ' && c != '\t' && c != '\n' && c != '\r') {
                token.push_back(c);
            }
        }
    }

    // Optional single preferred host from EICDIS_XRD_HOST.
    const char* envHost = gSystem->Getenv("EICDIS_XRD_HOST");
    if (envHost && envHost[0] != '\0') {
        std::string preferred(envHost);
        while (!preferred.empty() && preferred.back() == '/') {
            preferred.pop_back();
        }
        if (!preferred.empty()) {
            std::vector<std::string> reordered;
            reordered.push_back(preferred);
            for (const auto& h : hosts) {
                if (h != preferred) {
                    reordered.push_back(h);
                }
            }
            hosts = reordered;
        }
    }

    // Fallback defaults.
    for (const auto& h : kDefaultXrdHosts) {
        bool exists = false;
        for (const auto& host : hosts) {
            if (host == h) {
                exists = true;
                break;
            }
        }
        if (!exists) {
            hosts.push_back(h);
        }
    }

    const std::string volatilePrefix = "/volatile/eic/EPIC/RECO/" + campaign + "/epic_craterlake/";
    const std::string eicPrefix = "/eic/EPIC/RECO/" + campaign + "/epic_craterlake/";
    // Newer campaigns are generally under /eic; keep /volatile as fallback.
    const std::vector<std::string> prefixes = {eicPrefix, volatilePrefix};

    for (const auto& prefix : prefixes) {
        const std::string remoteDir = prefix + target;
        for (const auto& host : hosts) {
            const std::string command = "xrdfs " + host + " ls -R " + remoteDir + " 2>/dev/null || true";
            TString cmdOut = gSystem->GetFromPipe(command.c_str());
            if (cmdOut.Length() == 0) {
                continue;
            }

            std::vector<std::string> files;
            TObjArray* lines = cmdOut.Tokenize("\n");
            if (!lines) {
                continue;
            }

            for (int i = 0; i < lines->GetEntriesFast(); ++i) {
                TObject* obj = lines->At(i);
                if (!obj) {
                    continue;
                }

                std::string line = obj->GetName();
                if (line.empty()) {
                    continue;
                }

                const std::string rootSuffix = ".root";
                if (line.size() < rootSuffix.size() || line.compare(line.size() - rootSuffix.size(), rootSuffix.size(), rootSuffix) != 0) {
                    continue;
                }
                files.push_back(host + "/" + line);
            }

            delete lines;
            if (!files.empty()) {
                if (selectedHost) {
                    *selectedHost = host;
                }
                if (selectedPrefix) {
                    *selectedPrefix = prefix;
                }
                return files;
            }
        }
    }

    if (selectedHost) {
        *selectedHost = "";
    }
    if (selectedPrefix) {
        *selectedPrefix = "";
    }
    return {};
}

std::string BuildTargetForRegion(int beam_type, int Ee, int Eh, int r)
{
    if (beam_type == AnaManager::EHE3)
    {
        std::string gen_group = Eh == 166 ? "DIS/BeAGLE1.03.02-3.1/" : "DIS/BeAGLE1.03.02-1.2/";
        std::string beam_group = Form("eHe3/%dx%d/", Ee, Eh);
        std::string sample_group;
        if (r == 0)
        {
            sample_group = Eh == 166 ? "q2_1to10/" : "q2_2to10/";
        }     
        else if (r == 1)
            sample_group = "q2_10to100/";
        else
            sample_group = "q2_100to10000/";

        return gen_group + beam_group + sample_group;
    }

    if (beam_type == AnaManager::EP || beam_type == AnaManager::EP_CC)
    {
        std::string gen_group = beam_type == AnaManager::EP ? "DIS/NC/" : "DIS/CC/";
        std::string beam_group = Form("%dx%d/", Ee, Eh);
        std::string sample_group = Form("minQ2=%.0f/", pow(10, r));
        return gen_group + beam_group + sample_group;
    }

    if (beam_type == AnaManager::PI_BG)
    {
        std::string gen_group = "SIDIS/pythia6-eic/1.0.0/";
        std::string beam_group = Form("%dx%d/", Ee, Eh);
        std::string sample_group = "q2_0to1/";
        return gen_group + beam_group + sample_group;
    }

    if (beam_type == AnaManager::BEAM_BG)
    {
        std::string gen_group = "Bkg_Exact1S_2us/GoldCt/10um/DIS/NC/";
        std::string beam_group = Form("%dx%d/", Ee, Eh);
        std::string sample_group = Form("minQ2=%.0f/", pow(10, r));
        return gen_group + beam_group + sample_group;
    }

    if (beam_type == AnaManager::EP_PYTHIA6)
    {
        std::string gen_group = "DIS/pythia6.428-1.0/NC/noRad/ep/";
        std::string beam_group = Form("%dx%d/", Ee, Eh);
        std::string sample_group = Form("q2_%dto%d/", (int)pow(10, r), (int)pow(10, r + 1));

        if (Ee == 10 && Eh == 100)
        {
            gen_group = "SIDIS/pythia6-eic/1.2.0/ep_noradcor/";

            if (r == 3)
                sample_group = Form("q2_%dto%d/", (int)pow(10, r), (int)pow(10, r + 2));
        }

        return gen_group + beam_group + sample_group;
    }

    return "";
}

} // namespace