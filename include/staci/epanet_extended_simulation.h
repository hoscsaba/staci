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

#endif
