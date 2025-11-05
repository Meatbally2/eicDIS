#ifndef ANAMANAGER_HH
#define ANAMANAGER_HH

class AnaManager{

public:

	AnaManager(std::string ana_name_);
    ~AnaManager();

    void Initialize(bool is_select_region_, int region_index_, int is_analyse_all_files_, bool is_analyse_protons_);

    std::string GetOutputName();
    vector<std::string> GetInputNames();

private:
    bool is_select_region;
    int region_index;
    int starting_file;
    bool is_analyse_protons;
    std::string ana_name;
};

#endif
