#ifndef ANAMANAGER_HH
#define ANAMANAGER_HH

class AnaManager{

public:

	AnaManager(std::string ana_name_);
    ~AnaManager();

    void Initialize(bool is_select_region_, int region_index_, int starting_file_index_, int beam_type_);
    void InitializeForLocal(std::string type_);
    void SetBeamEnergy(int Ee_, int Eh_);

    std::string GetOutputName();
    vector<std::string> GetInputNames();
    vector<std::string> GetLocalInputNames();
    vector<std::string> GetLowQInputNames();
    vector<std::string> ValidateFiles(std::vector<std::string>& fileNames);

    std::string campaign;

    enum{EHE3,EP,PI_BG,BEAM_BG,EP_PYTHIA6,EP_DVMP};
    
private:
    int Ee;
    int Eh;
    bool is_select_region;
    int region_index;
    int file_index;
    int starting_file;
    int beam_type;
    std::string ana_name;
    std::string file_type;
};

#endif
