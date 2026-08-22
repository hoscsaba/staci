#ifndef STACI_EPANET_PUMP_H
#define STACI_EPANET_PUMP_H

#include <string>
#include <utility>
#include <vector>

struct EpanetPumpMetadata {
    std::string definition;
    std::string head_curve_id;
    std::vector<std::pair<double, double> > head_curve_points;
    double base_speed = 1.0;
    bool initial_setting_specified = false;
    double initial_setting = 1.0;
    std::string speed_pattern_id;
    std::vector<double> speed_pattern_values;

    bool energy_price_specified = false;
    bool energy_price_global = false;
    double energy_price = 0.0;
    std::string energy_pattern_id;
    bool energy_pattern_global = false;
    std::vector<double> energy_pattern_values;
    std::string efficiency_curve_id;
    bool efficiency_curve_global = false;
    std::vector<std::pair<double, double> > efficiency_curve_points;
    bool demand_charge_specified = false;
    double demand_charge = 0.0;
};

class EpanetPumpConfigurable {
public:
    virtual ~EpanetPumpConfigurable() = default;
    virtual void SetEpanetPumpMetadata(const EpanetPumpMetadata &metadata) = 0;
    virtual const EpanetPumpMetadata &GetEpanetPumpMetadata() const = 0;
    virtual void SetOperatingSpeed(double speed) = 0;
    virtual double GetOperatingSpeed() const = 0;
};

#endif
