#ifndef EPANET_POWER_PUMP_H
#define EPANET_POWER_PUMP_H

#include "Agelem.h"

class EpanetPowerPump : public Agelem {
public:
    EpanetPowerPump(const string &name,
                    const string &node_from,
                    const string &node_to,
                    double density,
                    double hydraulic_power_watts,
                    double initial_mass_flow_rate);

    double f(vector<double> x) override;
    vector<double> df(vector<double> x) override;
    void Ini(int mode, double value) override;
    double Get_dprop(string property) override;
    string Info() override;
    string GetType() const override { return "EpanetPowerPump"; }

private:
    double hydraulic_power;
    double safe_mass_flow() const;
};

#endif
