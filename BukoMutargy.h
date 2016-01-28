#include "Agelem.h"

class BukoMutargy:public Agelem
  {
  private:
    /// Fen�kszint
    double Hf;
    /// Nyitva van?
    bool is_opened;
    /// Buk�si sz�less�g
    double width;
    /// Buk�si magass�g
    double overflow_height;
    /// �t�ml�si t�nyez�
    double discharge_coeff;
    /// Szelep ellen�ll�st�nyez�
    double valve_coeff;
    
  public:
    /// Konstruktor
    BukoMutargy(const string nev, const string a_cspe_nev, const string a_cspv_nev, const double a_ro, const double Aref,
                const double Hf, const bool is_opened, const double width, const double overflow_height,
                const double discharge_coeff, const double valve_coeff, const double a_mp);
    /// Destruktor
    ~BukoMutargy();
    /// Inform�ci�
    string Info();
    /// �gegyenlet �rt�ke
    double f(vector<double>);
    /// �gegyenlet linariz�ltja
    vector<double> df(vector<double>);
    /// Inicializ�ci�
    void Ini(int mode, double value);
    /// Keresztmetszeti jellemzok szamitasa
    string GetType() {
        return "BukoMutargy";
    }
    void Set_dprop(string mit, double mire);
  };


