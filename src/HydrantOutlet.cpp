#include "HydrantOutlet.h"
#include <cmath>
#include <stdexcept>

HydrantOutlet::HydrantOutlet(const string &id, const string &node,
                             double density, double area, double coefficient)
    : Agelem(id, area, 0.0, density) {
  if (!std::isfinite(area) || area <= 0 || !std::isfinite(coefficient) ||
      coefficient <= 0 || !std::isfinite(density) || density <= 0)
    throw std::invalid_argument(
        "Hydrant area, density and total K must be positive and finite");
  resistance_ = coefficient / (2 * g * density * density * area * area);
  if (!std::isfinite(resistance_) || resistance_ <= 0)
    throw std::invalid_argument("Hydrant resistance is outside numeric range");
  csp_db = 1;
  cspe_nev = node;
  cspv_nev = "<atmosphere>";
  enabled = false;
}

bool HydrantOutlet::flowing(const vector<double> &x) const {
  return enabled && (x.at(0) > 0 || mp > 0);
}

double HydrantOutlet::f(const vector<double> &x) {
  // Semismooth active set: zero outflow at nonpositive pressure, no inflow.
  if (!flowing(x))
    return mp;
  return resistance_ * mp * std::fabs(mp) - x.at(0);
}

vector<double> HydrantOutlet::df(const vector<double> &x) {
  if (!flowing(x))
    return {0, 0, 1, 0};
  return {-1, 0, 2 * resistance_ * std::fabs(mp), 0};
}

void HydrantOutlet::Ini(int mode, double value) {
  mp = enabled && mode ? value : 0;
}

double HydrantOutlet::discharge(double head) const {
  return head > 0 ? std::sqrt(head / resistance_) / ro : 0;
}

double HydrantOutlet::Get_dprop(const string &property) {
  if (property == "mass_flow_rate")
    return mp;
  if (property == "headloss")
    return resistance_ * mp * std::fabs(mp);
  if (property == "length" || property == "headloss_per_unit_length")
    return 0;
  throw std::invalid_argument("Unsupported hydrant property: " + property);
}

void HydrantOutlet::Set_dprop(const string &property, double) {
  throw std::invalid_argument("Unsupported hydrant property: " + property);
}
