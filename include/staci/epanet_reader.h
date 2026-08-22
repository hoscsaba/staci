#ifndef STACI_EPANET_READER_H
#define STACI_EPANET_READER_H

#include <cstddef>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include "epanet_document.h"

class Agelem;
class Csomopont;

class EpanetReader {
public:
    explicit EpanetReader(const std::string &filename, bool extended_period = false);

    void load_system(std::vector<std::unique_ptr<Csomopont> > &nodes,
                     std::vector<std::unique_ptr<Agelem> > &edges);
    std::string read_setting(const std::string &name) const;
    const EpanetDocument &document() const { return document_; }

private:
    struct Record {
        std::size_t line_number;
        std::string text;
        std::vector<std::string> fields;
    };

    struct Warning {
        std::string section;
        std::string element;
        std::size_t line_number;
        std::string message;
    };

    std::string filename_;
    EpanetDocument document_;
    std::map<std::string, std::vector<Record> > sections_;
    std::map<std::string, std::vector<double> > patterns_;
    std::map<std::string, std::vector<std::pair<double, double> > > curves_;
    std::map<std::string, std::string> settings_;
    std::vector<Warning> warnings_;

    std::string flow_units_;
    std::string headloss_model_;
    std::string default_pattern_;
    std::string quality_mode_;
    std::string quality_chemical_name_;
    std::string quality_units_;
    std::string quality_trace_node_;
    double specific_gravity_;
    double demand_multiplier_;
    bool us_customary_units_;
    bool extended_period_;

    void parse_file();
    void parse_options();
    void parse_patterns();
    void parse_curves();
    void configure_settings();
    void warn_for_unsupported_sections();
    void print_report(std::size_t node_count,
                      std::size_t edge_count,
                      std::size_t demand_component_count,
                      std::size_t initial_quality_count,
                      std::size_t reservoir_head_pattern_count) const;

    const std::vector<Record> &records(const std::string &section) const;
    bool has_records(const std::string &section) const;
    void add_warning(const std::string &section,
                     const std::string &element,
                     std::size_t line_number,
                     const std::string &message);

    double number(const Record &record,
                  std::size_t field,
                  const std::string &label,
                  double default_value);
    double flow_to_m3_per_hour(double value) const;
    double length_to_metres(double value) const;
    double pipe_diameter_to_metres(double value) const;
    double tank_diameter_to_metres(double value) const;
    double pump_power_to_watts(double value) const;
    double first_pattern_multiplier(const std::string &pattern_id) const;
};

#endif
