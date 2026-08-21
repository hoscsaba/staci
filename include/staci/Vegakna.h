#include "Agelem.h"

class Vegakna:public Agelem{
  private:
    /// Fenekszint [m]-ben
    double Hf;
    /// Vizszint a fenektol [m]-ben
    double H;
    /// A teljes nyomasszint a geodetikus alapszintre sz�molva (ro*g*(Hf+H))
    double p;
  public:
    Vegakna(const string &nev,
        const string &cspnev,
        const double a_ro, 
        const double Aref, 
        const double Hf, 
        const double H, 
        const double a_mp, 
        const double a_tt);
    ~Vegakna() override;
    string Info() override;
    double f(const vector<double> &state) override;
    vector<double> df(const vector<double> &state) override;
    void Ini(int mode, double value) override;
    void Set_dprop(const string &property, double value) override;
    string_view GetType() const noexcept override {
        return "Vegakna";
    }
    double Get_dprop(const string &property) override;
  };
