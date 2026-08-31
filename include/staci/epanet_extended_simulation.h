#ifndef STACI_EPANET_EXTENDED_SIMULATION_H
#define STACI_EPANET_EXTENDED_SIMULATION_H

#include <string>

class Staci;

class EpanetExtendedSimulation {
public:
    EpanetExtendedSimulation(const std::string &input_filename,
                             const std::string &output_prefix);
    void run(Staci &system);

private:
    std::string input_filename_;
    std::string output_prefix_;
};

class EpanetSteadyQualitySimulation {
public:
    EpanetSteadyQualitySimulation(const std::string &input_filename,
                                  const std::string &output_prefix,
                                  const std::string &quality_mode,
                                  bool compute_sensitivity);
    void run(Staci &system);

private:
    std::string input_filename_;
    std::string output_prefix_;
    std::string quality_mode_;
    bool compute_sensitivity_;
};

#endif
