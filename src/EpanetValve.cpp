#include "EpanetValve.h"
#include "diagnostics.h"
#include <sstream>
#include <algorithm>
#include <cmath>
#include <stdexcept>

EpanetValve::EpanetValve(const string& id,const string& a,const string& b,
 double rho,double area,Kind kind,double raw,double factor,double minor,
 const string& curve_id,const vector<pair<double,double>>& curve)
 : JelleggorbesFojtas(id,a,b,rho,area,{0,100},{0,0},50,1), kind_(kind),
 raw_setting_(raw),factor_(factor),setting_(raw*factor),minor_(minor),
 resistance_(minor*0.05093871858091765/(area*area*rho*rho)),curve_id_(curve_id),curve_(curve) { state_=EpanetTcvStatus::Open; }
const char* EpanetValve::valve_type() const noexcept {
 return kind_==Kind::PSV ? "PSV" : kind_==Kind::PBV ? "PBV" : kind_==Kind::PRV ? "PRV" : kind_==Kind::FCV ? "FCV" : "GPV";
}
void EpanetValve::SetEpanetTcvSetting(double value) {
 if(kind_==Kind::GPV) throw std::invalid_argument("GPV '"+nev+"' uses a curve ID; numeric settings are not supported.");
 if(!std::isfinite(value)||(value<0 && kind_!=Kind::PRV && kind_!=Kind::PSV)) throw std::invalid_argument("Valve '"+nev+"' requires a finite setting; negative values are supported only for PRV/PSV pressure settings.");
 raw_setting_=value;setting_=value*factor_;fixed_=false;state_=EpanetTcvStatus::Active;enabled=true;
}
void EpanetValve::SetEpanetTcvStatus(EpanetTcvStatus s) {
 if(kind_==Kind::GPV && s==EpanetTcvStatus::Active) throw std::invalid_argument("GPV '"+nev+"' supports OPEN/CLOSED status; ACTIVE regulation requires a PRV or FCV.");
 state_=s;fixed_=s!=EpanetTcvStatus::Active;enabled=s!=EpanetTcvStatus::Closed;
}
pair<double,double> EpanetValve::loss() const {
 if(kind_!=Kind::GPV) {
  const double gradient=2*resistance_*std::abs(mp), minimum=1.07639104167e-9;
  if(resistance_==0) return {1.07639104167e-8*mp,1.07639104167e-8};
  if(gradient<minimum) return {minimum*mp,minimum};
  return {resistance_*mp*std::abs(mp),gradient};
 }
 double q=std::abs(mp)/ro;size_t i=1;
 while(i+1<curve_.size() && q>curve_[i].first) ++i;
 double slope=(curve_[i].second-curve_[i-1].second)/(curve_[i].first-curve_[i-1].first);
 double h=curve_[i-1].second+slope*(q-curve_[i-1].first);
 return {(mp<0?-1:1)*h,std::max(slope/ro,1e-12)};
}
double EpanetValve::f(const vector<double>& x) {
 if(!enabled || state_==EpanetTcvStatus::Closed) return x[1]+x[3]-x[0]-x[2]+1.07639104167e6*mp;
 if(state_==EpanetTcvStatus::Active && kind_==Kind::PSV) return x[0]-setting_-1e-9*mp;
 if(!fixed_ && kind_==Kind::PBV && setting_>0 && resistance_*mp*mp<=setting_) return x[1]+x[3]-x[0]-x[2]+setting_+1e-9*mp;
 if(state_==EpanetTcvStatus::Active && kind_==Kind::PRV) return x[1]-setting_+1e-9*mp;
 if(state_==EpanetTcvStatus::Active && kind_==Kind::FCV) return mp-setting_*ro+1e-9*(x[1]+x[3]-x[0]-x[2]);
 return x[1]+x[3]-x[0]-x[2]+loss().first;
}
vector<double> EpanetValve::df(const vector<double>&) {
 if(!enabled || state_==EpanetTcvStatus::Closed) return {-1,1,1.07639104167e6,0};
 if(state_==EpanetTcvStatus::Active && kind_==Kind::PSV) return {1,0,-1e-9,0};
 if(!fixed_ && kind_==Kind::PBV && setting_>0 && resistance_*mp*mp<=setting_) return {-1,1,1e-9,0};
 if(state_==EpanetTcvStatus::Active && kind_==Kind::PRV) return {0,1,1e-9,0};
 if(state_==EpanetTcvStatus::Active && kind_==Kind::FCV) return {-1e-9,1e-9,1,0};
 return {-1,1,loss().second,0};
}
bool EpanetValve::update_status(const vector<double>& x) {
 if(!enabled || fixed_ || kind_==Kind::GPV || kind_==Kind::PBV) return false;
 const double h1=x[0]+x[2],h2=x[1]+x[3],tol=0.0001524,qtol=2.8316846592e-8*ro;
 auto s=state_;
 if(kind_==Kind::PRV) {
  const double target=x[3]+setting_,hml=resistance_*mp*mp;
  if(s==EpanetTcvStatus::Active) {
   if(mp < -qtol) s=EpanetTcvStatus::Closed;
   else if(h1-hml < target-tol) s=EpanetTcvStatus::Open;
  } else if(s==EpanetTcvStatus::Open) {
   if(mp < -qtol) s=EpanetTcvStatus::Closed;
   else if(h2 >= target+tol) s=EpanetTcvStatus::Active;
  } else {
   if(h1>=target+tol && h2<target-tol) s=EpanetTcvStatus::Active;
   else if(h1<target-tol && h1>h2+tol) s=EpanetTcvStatus::Open;
  }
 } else if(kind_==Kind::PSV) {
  const double target=x[2]+setting_,hml=resistance_*mp*mp;
  if(s==EpanetTcvStatus::Active) {
   if(mp < -qtol) s=EpanetTcvStatus::Closed;
   else if(h2+hml > target+tol) s=EpanetTcvStatus::Open;
  } else if(s==EpanetTcvStatus::Open) {
   if(mp < -qtol) s=EpanetTcvStatus::Closed;
   else if(h1 < target-tol) s=EpanetTcvStatus::Active;
  } else {
   if(h2>target+tol && h1>h2+tol) s=EpanetTcvStatus::Open;
   else if(h1>=target+tol && h1>h2+tol) s=EpanetTcvStatus::Active;
  }
 } else {
  if(h1-h2 < -tol || mp < -qtol) s=EpanetTcvStatus::Open;
  else if(s==EpanetTcvStatus::Open && mp>=setting_*ro) s=EpanetTcvStatus::Active;
 }
 bool changed=s!=state_;state_=s;return changed;
}
double EpanetValve::Get_dprop(const string& p) {
 if(p=="status") return !Is_enabled()?0:state_==EpanetTcvStatus::Open?1:2;
 if(p=="diameter") return std::sqrt(4*Aref/3.14159265358979323846);
 if(p=="tcv_setting" || p=="setting") return raw_setting_;
 if(p=="minor_loss" || p=="tcv_minor_loss") return minor_;
 if(p=="mass_flow_rate") return mp;
 return JelleggorbesFojtas::Get_dprop(p);
}
void EpanetValve::Set_dprop(const string& p,double v) {
 if(p=="tcv_setting" || p=="setting") SetEpanetTcvSetting(v);
 else if(p=="minor_loss" || p=="tcv_minor_loss" || p=="diameter") {
  if(!std::isfinite(v) || v<0 || (p=="diameter" && v==0))
   throw std::invalid_argument("Valve '"+nev+"' has an invalid "+p+" value.");
  if(p=="diameter") Aref=3.14159265358979323846*v*v/4;
  else minor_=v;
  resistance_=minor_*0.05093871858091765/(Aref*Aref*ro*ro);
 }
 else if(p=="status") SetEpanetTcvStatus(v==0?EpanetTcvStatus::Closed:v==1?EpanetTcvStatus::Open:EpanetTcvStatus::Active);
 else JelleggorbesFojtas::Set_dprop(p,v);
}

void EpanetValve::report_operating_status() {
 if(Is_enabled() && kind_==Kind::FCV && !fixed_ && state_==EpanetTcvStatus::Open && Get_Q()<setting_-1e-8) {
  std::ostringstream message;
  message << "FCV '" << nev << "' cannot maintain its flow setting (" << setting_
          << " m3/s); available pressure gives " << Get_Q()
          << " m3/s. The valve is fully open, as in EPANET. Check source heads, downstream pressure and the requested flow.";
  diagnostics::warning("EPANET.FCV_UNATTAINABLE",message.str());
 }
}

// The tiny hydraulic coupling terms regularize the Newton matrix. They must
// not turn a physically closed valve or an infeasible flow setpoint into a
// water supply at enormous pressure differences.
std::string EpanetValve::operating_violation() const {
 const double q = mp / ro;
 std::ostringstream message;
 if ((!enabled || state_ == EpanetTcvStatus::Closed) && std::abs(q) > 1e-6) {
  message << valve_type() << " '" << nev << "' is closed but the computed state requires "
          << q << " m3/s through it. Closed-link numerical regularization cannot supply a disconnected demand. "
          << "Check sources and demands on both sides of the valve, pump shutdowns and flow directions.";
 } else if (enabled && kind_ == Kind::FCV && state_ == EpanetTcvStatus::Active &&
            std::abs(q - setting_) > std::max(1e-6, std::abs(setting_) * 1e-4)) {
  message << "FCV '" << nev << "' is active at " << setting_
          << " m3/s but continuity requires " << q << " m3/s in the computed state. "
          << "The flow setting conflicts with the network supply/demand balance. "
          << "Check alternate sources, downstream demand and the valve setting; tightening tolerances cannot remove this conflict.";
 }
 return message.str();
}
