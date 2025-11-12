#ifndef ANAMANAGER_HH
#define ANAMANAGER_HH

class AnaManager{

public:

	AnaManager(std::string ana_name_);
    ~AnaManager();

    void Initialize(bool is_select_region_, int region_index_, int starting_file_index_, bool is_analyse_protons_);
    void InitializeForLocal(std::string type_);

    std::string GetOutputName();
    vector<std::string> GetInputNames();
    vector<std::string> GetLocalInputNames();
    vector<std::string> GetLowQInputNames();

private:
    bool is_select_region;
    int region_index;
    int starting_file;
    bool is_analyse_protons;
    std::string ana_name;
    std::string file_type;
};

#endif
