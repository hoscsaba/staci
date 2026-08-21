#include "EpanetPowerPump.h"

#include <cmath>
#include <sstream>
#include <stdexcept>

EpanetPowerPump::EpanetPowerPump(const string &name,
                                 const string &node_from,
                                 const string &node_to,
                                 double density,
                                 double hydraulic_power_watts,
                                 double initial_mass_flow_rate)
    : Agelem(name, 1.0, initial_mass_flow_rate, density),
      hydraulic_power(hydraulic_power_watts) {
    csp_db = 2;
    cspe_nev = node_from;
    cspv_nev = node_to;
}

double EpanetPowerPump::safe_mass_flow() const {
    const double minimum_mass_flow = 1.0e-6;
    if (std::fabs(mp) >= minimum_mass_flow)
        return mp;
    return mp < 0.0 ? -minimum_mass_flow : minimum_mass_flow;
}

double EpanetPowerPump::f(const vector<double> &x) {
    if (!enabled)
        return mp;
    const double flow = safe_mass_flow();
    const double pump_head = hydraulic_power / (g * flow);
    return (x[1] - x[0]) + (x[3] - x[2]) - pump_head;
}

vector<double> EpanetPowerPump::df(const vector<double> &) {
    if (!enabled)
        return vector<double>{0.0, 0.0, 1.0, 0.0};
    const double flow = safe_mass_flow();
    vector<double> result;
    result.push_back(-1.0);
    result.push_back(1.0);
    result.push_back(hydraulic_power / (g * flow * flow));
    result.push_back(0.0);
    return result;
}

void EpanetPowerPump::Ini(int mode, double value) {
    if (mode != 0)
        mp = enabled ? value : 0.0;
}

double EpanetPowerPump::Get_dprop(const string &property) {
    if (property == "Aref")
        return Aref;
    if (property == "mass_flow_rate")
        return mp;
    if (property == "power")
        return hydraulic_power;
    if (property == "headloss" || property == "headloss_per_unit_length")
        return std::fabs(hydraulic_power / (g * safe_mass_flow()));
    if (property == "concentration" || property == "konc_atlag")
        return konc_atlag;
    return 0.0;
}

void EpanetPowerPump::Set_dprop(const string &property, double value) {
    if (property == "power")
        hydraulic_power = value;
    else if (property == "concentration" || property == "konc_atlag")
        konc_atlag = value;
    else
        throw std::invalid_argument("Unsupported EpanetPowerPump property: " + property);
}

string EpanetPowerPump::Info() {
    std::ostringstream out;
    out << Agelem::Info();
    out << "\n       nodes : " << cspe_nev << " -> " << cspv_nev;
    out << "\n       power : " << hydraulic_power << " W\n";
    return out.str();
}
