#include "Agelem.h"

class JelleggorbesFojtas : public Agelem
{
private:
  vector<double> e, zeta, p;
  double allas, veszt;
  int fokszam;
  double headloss;
public:
  JelleggorbesFojtas(const string nev,
   const string a_cspe_nev,
   const string a_cspv_nev,
   const double a_ro,
   const double Aref,
   const vector <double>e,
   const vector <double> zeta,
   const double allas,
   const double a_mp);
  void Update_zeta();
  ~JelleggorbesFojtas();
  string Info();
  double f(vector<double>);
  vector<double> df(vector<double>);
  void Ini(int mode, double value);
  void Set_dprop(string mit, double mire);
  double Get_dprop(string mit);
  string GetType()
  {
    return "JelleggorbesFojtas";
  }
  // VIGYAZAT! VALAMIERT VAN EGY Get_Tipus is, ami "Jelleggorbes fojtas"-t ad vissza!!!
};

