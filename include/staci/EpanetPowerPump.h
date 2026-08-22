#ifndef EPANET_POWER_PUMP_H
#define EPANET_POWER_PUMP_H

#include "Agelem.h"
#include "EpanetPump.h"

class EpanetPowerPump : public Agelem, public EpanetPumpConfigurable {
public:
    EpanetPowerPump(const string &name,
                    const string &node_from,
                    const string &node_to,
                    double density,
                    double hydraulic_power_watts,
                    double initial_mass_flow_rate);

    double f(const vector<double> &state) override;
    vector<double> df(const vector<double> &state) override;
    void Ini(int mode, double value) override;
    double Get_dprop(const string &property) override;
    void Set_dprop(const string &property, double value) override;
    string Info() override;
    string_view GetType() const noexcept override { return "EpanetPowerPump"; }
    void SetEpanetPumpMetadata(const EpanetPumpMetadata &value) override;
    const EpanetPumpMetadata &GetEpanetPumpMetadata() const override { return metadata; }
    void SetOperatingSpeed(double value) override;
    double GetOperatingSpeed() const override { return operating_speed; }

private:
    double hydraulic_power;
    double operating_speed;
    EpanetPumpMetadata metadata;
    double safe_mass_flow() const;
    double effective_power() const;
};

#endif
