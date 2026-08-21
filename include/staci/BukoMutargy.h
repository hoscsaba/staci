#include "Agelem.h"

class BukoMutargy: public Agelem
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
    BukoMutargy(const string &nev, const string &a_cspe_nev, const string &a_cspv_nev, const double a_ro, const double Aref,
                const double Hf, const bool is_opened, const double width, const double overflow_height,
                const double discharge_coeff, const double valve_coeff, const double a_mp);
    /// Destruktor
    ~BukoMutargy() override;
    /// Információ
    string Info() override;
    /// Ágegyenlet értéke
    double f(const vector<double> &state) override;
    /// Ágegyenlet linarizáltja
    vector<double> df(const vector<double> &state) override;
    /// Inicializáció
    void Ini(int mode, double value) override;
    /// Keresztmetszeti jellemzok szamitasa
    string_view GetType() const noexcept override {
        return "BukoMutargy";
    }
    void Set_dprop(const string &property, double value) override;
    double Get_dprop(const string &property) override;
};
