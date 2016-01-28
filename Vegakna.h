#include "Agelem.h"

class Vegakna:public Agelem{
  private:
    /// Fenékszint [m]-ben
    double Hf;
    /// Vízszint a fenéktõl [m]-ben
    double H;
    /// A teljes nyomásszint a geodetikus alapszintre számolva (ro*g*(Hf+H))
    double p;
  public:
    Vegakna(const string nev, 
        const string cspnev, 
        const double a_ro, 
        const double Aref, 
        const double Hf, 
        const double H, 
        const double a_mp, 
        const double a_tt);
    ~Vegakna();
    string Info();
    double f(vector<double>);
    vector<double> df(vector<double>);
    void Ini(int mode, double value);
    void Set_dprop(string mit, double mire);
    string GetType() {
        return "Vegakna";
    }
    double Get_dprop(string mit);
  };
