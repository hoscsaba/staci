#include "Agelem.h"

class BukoMutargy:public Agelem
  {
  private:
    /// Fenékszint
    double Hf;
    /// Nyitva van?
    bool is_opened;
    /// Bukási szélesség
    double width;
    /// Bukási magasság
    double overflow_height;
    /// Átömlési tényezõ
    double discharge_coeff;
    /// Szelep ellenállástényezõ
    double valve_coeff;
    
  public:
    /// Konstruktor
    BukoMutargy(const string nev, const string a_cspe_nev, const string a_cspv_nev, const double a_ro, const double Aref,
                const double Hf, const bool is_opened, const double width, const double overflow_height,
                const double discharge_coeff, const double valve_coeff, const double a_mp);
    /// Destruktor
    ~BukoMutargy();
    /// Információ
    string Info();
    /// Ágegyenlet értéke
    double f(vector<double>);
    /// Ágegyenlet linarizáltja
    vector<double> df(vector<double>);
    /// Inicializáció
    void Ini(int mode, double value);
    /// Keresztmetszeti jellemzok szamitasa
    string GetType() {
        return "BukoMutargy";
    }
    void Set_dprop(string mit, double mire);
  };


