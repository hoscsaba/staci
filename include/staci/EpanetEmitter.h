#pragma once
#include "Agelem.h"
#include "diagnostics.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>
// A one-node outlet to atmosphere. Positive flow leaves the junction.
class EpanetEmitter final : public Agelem {
 double coefficient_, exponent_;
 std::pair<double,double> loss() const {
  const double q=mp/ro, power=1.0/exponent_;
  const double gradient=power*std::pow(std::max(std::abs(q)/coefficient_,1e-30),power-1)/coefficient_/ro;
  const double minimum=1.07639104167e-9;
  if(gradient<minimum) return {minimum*mp,minimum};
  return {std::copysign(std::pow(std::abs(q)/coefficient_,power),q),gradient};
 }
public:
 EpanetEmitter(const string& id,const string& node,double density,double coefficient,double exponent)
 : Agelem(id,1,0,density),coefficient_(coefficient),exponent_(exponent) {
  if(!std::isfinite(density)||density<=0 || !std::isfinite(coefficient)||coefficient<=0 || !std::isfinite(exponent)||exponent<=0) throw std::invalid_argument("Emitter +id+ requires positive finite density, coefficient and exponent.");
  csp_db=1;cspe_nev=node;cspv_nev="<atmosphere>";FolyTerf=0;
 }
 double f(const vector<double>& x) override {return enabled ? x[0]-loss().first : mp;}
 vector<double> df(const vector<double>&) override {return enabled ? vector<double>{1,0,-loss().second,0} : vector<double>{0,0,1,0};}
 void ReportOperatingStatus() const {
  if(enabled && mp < -1e-5) diagnostics::warning("EPANET.EMITTER_BACKFLOW","Emitter '"+nev+"' at junction '"+cspe_nev+"' has negative pressure and atmospheric inflow. This reproduces EPANET's signed emitter equation; check whether this backflow is physical for the intended model.");
 }
 void InitializePressure(double pressure) { mp=ro*coefficient_*std::copysign(std::pow(std::abs(pressure),exponent_),pressure); }
 void Ini(int mode,double value) override {mp=mode?value:1;}
 string_view GetType() const noexcept override {return "EpanetEmitter";}
 double Get_dprop(const string& p) override {
  if(p=="mass_flow_rate") return mp;
  if(p=="coefficient") return coefficient_;
  if(p=="exponent") return exponent_;
  throw std::invalid_argument("Emitter '"+nev+"': unsupported property '"+p+"'.");
 }
 void Set_dprop(const string& p,double v) override {
  if(!std::isfinite(v)||v<=0) throw std::invalid_argument("Emitter '"+nev+"' requires a positive finite "+p);
  if(p=="coefficient") coefficient_=v;
  else if(p=="exponent") exponent_=v;
  else throw std::invalid_argument("Emitter '"+nev+"': unsupported property '"+p+"'.");
 }
};
