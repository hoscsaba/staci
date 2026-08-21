// data_io_rapidxml.h
#ifndef DATA_IO_RAPIDXML_H
#define DATA_IO_RAPIDXML_H

#include <string>
#include <fstream>
#include <vector>
#include <memory>
#include "Csomopont.h"
#include "Agelem.h"
#include "rapidxml.hpp"
#include "rapidxml_utils.hpp"
#include "rapidxml_print.hpp"

class data_io {
public:
    data_io(const char* xml_file);
    void load_system(std::vector<std::unique_ptr<Csomopont> >& cspok,
                     std::vector<std::unique_ptr<Agelem> >& agelemek);
    void load_ini_values(std::vector<Csomopont*>& cspok, std::vector<Agelem*>& agelemek);
    void save_results(double FolyMenny, double sum_of_inflow, double sum_of_demand,
                      const std::vector<Csomopont*>& cspok,
                      const std::vector<Agelem*>& agelemek,
                      bool conv_reached, int staci_debug_level);
    void save_mod_prop(const std::vector<Csomopont*>& cspok,
                       const std::vector<Agelem*>& agelemek,
                       const std::string& eID, const std::string& pID,
                       bool is_property_general);
    void save_mod_prop_all_elements(const std::vector<Csomopont*>& cspok,
                                    const std::vector<Agelem*>& agelemek,
                                    const std::string& pID);
    void save_transport(int mode,
                        const std::vector<Csomopont*>& cspok,
                        const std::vector<Agelem*>& agelemek);
    std::string read_setting(const std::string& which);

private:
    const char* xml_file;
    bool debug;
    int csp_db, ag_db;
    std::vector<std::string> edge_type;
    int edge_type_number;
    std::vector<int> edge_type_occur;
    void curve_reader(const std::string& id, rapidxml::xml_node<>* node,
                      std::vector<double>& px, std::vector<double>& py);
    double string_to_double(const std::string& s, const std::string& elem_name,
                            const std::string& tag_name, double def_value);

    // helpers
    rapidxml::xml_node<>* get_child(rapidxml::xml_node<>* parent, const char* name, int idx = 0);
    void replace_value(rapidxml::xml_node<>* node, const char* field,
                       double val, const std::string& msg);
};

#endif // DATA_IO_RAPIDXML_H
