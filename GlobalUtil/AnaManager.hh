#ifndef ANAMANAGER_HH
#define ANAMANAGER_HH

class AnaManager{

public:

	AnaManager(std::string ana_name_);
    ~AnaManager();

    void Initialize(bool is_select_region_, int region_index_, int starting_file_index_, bool is_analyse_protons_);
    void InitializeForLocal(std::string type_);
    void SetBeamEnergy(int Ee_, int Eh_);

    std::string GetOutputName();
    vector<std::string> GetInputNames();
    vector<std::string> GetLocalInputNames();
    vector<std::string> GetLowQInputNames();
    vector<std::string> ValidateFiles(std::vector<std::string>& fileNames);

    std::string campaign;
    
private:
    int Ee;
    int Eh;
    bool is_select_region;
    int region_index;
    int file_index;
    int starting_file;
    bool is_analyse_protons;
    std::string ana_name;
    std::string file_type;
};

#endif
