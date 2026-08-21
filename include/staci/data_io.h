#include <string>
#include <fstream>
#include <memory>
#include <vector>
#include "Csomopont.h"
#include "Agelem.h"
#include "xmlParser.h"

class EpanetReader;
class EpanetDocument;

class data_io
{

public:
    /// Konstruktor
    data_io(const char *xml_fnev, bool epanet_extended = false);
    ~data_io();
    void load_system(vector<Csomopont *> &cspok, vector<Agelem *> &agelemek);
    void load_ini_values(vector<Csomopont *> &cspok, vector<Agelem *> &agelemek);
    void save_results(double FolyMenny, double sum_of_inflow, double sum_of_demand, vector<Csomopont *> cspok, vector<Agelem *> agelemek, bool conv_reached, int staci_debug_level);
    void save_mod_prop(vector<Csomopont *> cspok, vector<Agelem *> agelemek, string eID, string pID, bool is_property_general);
    void save_mod_prop_all_elements(vector<Csomopont *> cspok, vector<Agelem *> agelemek, string pID);
    void save_transport(int mode, vector<Csomopont *> cspok, vector<Agelem *> agelemek);
    string read_setting(string which);
    const EpanetDocument *get_epanet_document() const;

private:
    /// debug info a kepernyore
    const char *xml_fnev;
    bool is_epanet_input;
    std::unique_ptr<EpanetReader> epanet_reader;
    /// debug info a kepernyore
    bool debug;
    /// csomopontok es agak szama
    int csp_db, ag_db;
    /// agelem tipusok
    vector<string> edge_type;
    int edge_type_number;
    vector<int> edge_type_occur;
    /// gorbek kiolvasasa
    void curve_reader(const string id, const XMLNode node, vector<double> &, vector<double> &);
    double string_to_double( const string &s , const string &elem_name, const string &tag_name, const double &def_value);
    void replace_value(XMLNode &root_node, string fieldname, double val, string msg);
//    void replace_value(XMLNode &root_node, string fieldname, string val, string msg);
    void replace_value2(XMLNode &root_node, string fieldname1, int i, string fieldname2, double val, string msg);
  //  void replace_value2(XMLNode &root_node, string fieldname1, int i, string fieldname2, string val, string msg);
    void replace_value4(XMLNode &root_node, string fieldname1, int i,
                             string fieldname2, string fieldname3, string fieldname4,  double val, string msg);
    void warn_epanet_write_unsupported(const string &operation) const;

};
