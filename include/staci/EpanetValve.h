#pragma once
#include <utility>
#include "JelleggorbesFojtas.h"

// EPANET control valves share the existing status/control interface.
class EpanetValve : public JelleggorbesFojtas {
public:
  enum class Kind { PRV, PSV, PBV, FCV, GPV };
  EpanetValve(const string&, const string&, const string&, double, double,
              Kind, double, double, double, const string&,
              const vector<pair<double,double>>&);
  bool Is_enabled() const override { return enabled && state_!=EpanetTcvStatus::Closed; }
  double f(const vector<double>&) override;
  vector<double> df(const vector<double>&) override;
  void report_operating_status();
  std::string operating_violation() const;
  bool update_status(const vector<double>&);
  void SetEpanetTcvSetting(double) override;
  void SetEpanetTcvStatus(EpanetTcvStatus) override;
  double GetEpanetTcvSetting() const noexcept override { return raw_setting_; }
  double GetEpanetTcvMinorLoss() const noexcept override { return minor_; }
  EpanetTcvStatus GetEpanetTcvStatus() const noexcept override { return enabled ? state_ : EpanetTcvStatus::Closed; }
  bool CanExportAsEpanetTcv() const noexcept override { return false; }
  double Get_dprop(const string&) override;
  void Set_dprop(const string&, double) override;
  string_view GetType() const noexcept override { return "EpanetValve"; }
  const char* valve_type() const noexcept;
  bool fixed_status() const noexcept { return fixed_; }
  double setting_si() const noexcept { return setting_; }
  const vector<pair<double,double>>& curve_points() const noexcept { return curve_; }
  const string& curve_id() const noexcept { return curve_id_; }
private:
  Kind kind_;
  double raw_setting_, factor_, setting_, minor_, resistance_;
  EpanetTcvStatus state_ = EpanetTcvStatus::Active;
  bool fixed_ = false;
  string curve_id_;
  vector<pair<double,double>> curve_;
  pair<double,double> loss() const;
};
