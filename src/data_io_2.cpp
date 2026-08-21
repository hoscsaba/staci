// data_io_rapidxml.cpp
#include "data_io_rapidxml.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cstring>
#include <algorithm>

using namespace rapidxml;
using namespace std;

// helper to find child by name and index
xml_node<>* data_io::get_child(xml_node<>* parent, const char* name, int idx) {
    xml_node<>* node = parent->first_node(name);
    while (node && idx--) node = node->next_sibling(name);
    return node;
}

data_io::data_io(const char* xml_file_)
    : xml_file(xml_file_), debug(false) {
    // initialize edge types
    edge_type = {"press","pipe","pump","valve","pool","channel","channel1","channel2","overflow","checkvalve"};
    edge_type_number = static_cast<int>(edge_type.size());
    edge_type_occur.assign(edge_type_number, 0);
    // read debug flag
    file<> fdoc(xml_file);
    xml_document<> doc;
    doc.parse<0>(fdoc.data());
    xml_node<>* root = doc.first_node("staci");
    xml_node<>* settings = get_child(root, "settings");
    if (settings) {
        xml_node<>* dbg = get_child(settings, "cpp_xml_debug");
        if (dbg && strcmp(dbg->value(), "true") == 0)
            debug = true;
    }
}

void data_io::load_system(vector<unique_ptr<Csomopont> >& cspok,
                          vector<unique_ptr<Agelem> >& agelemek) {
    file<> fdoc(xml_file);
    xml_document<> doc;
    doc.parse<0>(fdoc.data());
    xml_node<>* root = doc.first_node("staci");
    xml_node<>* nodes = get_child(root, "nodes");
    xml_node<>* edges = get_child(root, "edges");

    // count
    csp_db = 0;
    for (xml_node<>* n = nodes->first_node("node"); n; n = n->next_sibling("node")) {
        csp_db++;
    }
    ag_db = 0;
    for (xml_node<>* e = edges->first_node("edge"); e; e = e->next_sibling("edge")) {
        ag_db++;
    }
    if (debug) cout << "Nodes: " << csp_db << "  Edges: " << ag_db << endl;

    // load nodes
    for (xml_node<>* n = nodes->first_node("node"); n; n = n->next_sibling("node")) {
        xml_node<>* endn = get_child(n, "endnode");
        if (endn && strcmp(endn->value(), "false") == 0) {
            string id = get_child(n, "id")->value();
            double height = string_to_double(get_child(n, "height")->value(), id, "<height>", 0);
            double demand = string_to_double(get_child(n, "demand")->value(), id, "<demand>", 0);
            double cl_input = string_to_double(get_child(n, "cl_input")->value(), id, "<cl_input>", 0);
            double pressure = string_to_double(get_child(n, "pressure")->value(), id, "<pressure>", 0);
            double density = string_to_double(get_child(n, "density")->value(), id, "<density>", 0);
            double tt = string_to_double(get_child(n, "travel_time")->value(), id, "<travel_time>", 0);
            cspok.push_back(make_unique<Csomopont>(id, height, demand, cl_input, pressure, density, tt));
        }
    }

    // load edges
    for (xml_node<>* e = edges->first_node("edge"); e; e = e->next_sibling("edge")) {
        xml_node<>* spec = get_child(e, "edge_spec");
        int type_idx = -1;
        for (int j = 0; j < edge_type_number; ++j) {
            if (spec->first_node(edge_type[j].c_str())) { type_idx = j; break; }
        }
        string id = get_child(e, "id")->value();
        string from = get_child(e, "node_from")->value();
        string to   = get_child(e, "node_to")->value();
        string s_area = get_child(e, "aref")->value();
        replace(s_area.begin(), s_area.end(), ',', '.');
        double aref = atof(s_area.c_str());
        double density = atof(get_child(e, "density")->value());
        double mfr = atof(get_child(e, "mass_flow_rate")->value());
        double tt   = atof(get_child(e, "travel_time")->value());
        switch (type_idx) {
            case 0: { // press
                double pres = atof(get_child(spec, "press")->first_node("pressure")->value());
                agelemek.push_back(make_unique<KonstNyomas>(id, aref, from, density, pres, mfr, tt));
                break;
            }
            case 1: { // pipe
                double L = atof(get_child(spec, "pipe")->first_node("length")->value());
                double D = atof(get_child(spec, "pipe")->first_node("diameter")->value());
                double rough = atof(get_child(spec, "pipe")->first_node("roughness")->value());
                double cl_k = atof(get_child(spec, "pipe")->first_node("cl_k")->value());
                double cl_w = atof(get_child(spec, "pipe")->first_node("cl_w")->value());
                agelemek.push_back(make_unique<Cso>(id, from, to, density, L, D, rough, cl_k, cl_w, mfr));
                break;
            }
            // ... implement other cases similarly ...
            default:
                cerr << "Unknown edge spec for id=" << id << endl;
        }
    }
}

// ... implement other methods (curve_reader, save_results, etc.) similarly using rapidxml ...


double data_io::string_to_double(const string& s, const string&, const string&, double def) {
    string tmp = s;
    replace(tmp.begin(), tmp.end(), ',', '.');
    try { return stod(tmp); } catch (...) { return def; }
}
