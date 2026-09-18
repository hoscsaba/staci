#ifndef STACI_HYDRANT_OUTLET_H
#define STACI_HYDRANT_OUTLET_H
#include "Agelem.h"

// Atmospheric outlet at the junction elevation. K is TOTAL head resistance:
// h = K (Q/A)^2 / (2g), including outlet kinetic head. Outflow is positive.
class HydrantOutlet : public Agelem {
public:
  HydrantOutlet(const string &id, const string &node, double density,
                double area_m2, double total_loss_coefficient);
  string_view GetType() const noexcept override { return "HydrantOutlet"; }
  double f(const vector<double> &x) override;
  vector<double> df(const vector<double> &x) override;
  void Ini(int mode, double value) override;
  double Get_dprop(const string &property) override;
  void Set_dprop(const string &property, double value) override;
  double discharge(double pressure_head_m) const;

private:
  double resistance_;
  bool flowing(const vector<double> &x) const;
};
#endif
