#include "Agelem.h"

class JelleggorbesFojtas : public Agelem
{
private:
  vector<double> e, zeta, p;
  double allas, veszt;
  int fokszam;
  double headloss;
public:
  JelleggorbesFojtas(const string &nev,
   const string &a_cspe_nev,
   const string &a_cspv_nev,
   const double a_ro,
   const double Aref,
   const vector<double> &e,
   const vector<double> &zeta,
   const double allas,
   const double a_mp);
  void Update_zeta();
  ~JelleggorbesFojtas() override;
  string Info() override;
  double f(const vector<double> &state) override;
  vector<double> df(const vector<double> &state) override;
  void Ini(int mode, double value) override;
  void Set_dprop(const string &property, double value) override;
  double Get_dprop(const string &property) override;
  string_view GetType() const noexcept override
  {
    return "JelleggorbesFojtas";
  }
};

