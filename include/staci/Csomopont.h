#ifndef CSOMOPONT_H
#define CSOMOPONT_H

#include <string>
#include <vector>

struct EpanetDemandComponent
{
    double base_demand_m3s = 0.0;
    std::string pattern_id;
    std::vector<double> pattern_values;
    std::string category;
    bool primary = false;
};

struct EpanetInitialQuality
{
    bool specified = false;
    double source_value = 0.0;
    std::string mode = "NONE";
    std::string chemical_name;
    std::string units;
    std::string trace_node;
};

class Csomopont
{
public:
    /// Bemen� �gak nyilv�ntart�sa
    std::vector<int> ag_be;
    /// Kimen� �gak nyilv�ntart�sa
    std::vector<int> ag_ki;

private:
    /// A csomopont neve
    std::string nev;
    /// HEAD = p[Pa]/ro/g
    double p_head;
    /// Suruseg
    double ro;
    /// Geodetikus magassag
    double h;
    /// Fogyasztas
    double fogy;
    /// Betaplalt koncentracio
    double cl_be;
    /// Csomoponti atlagkoncentracio
    double konc_atlag;
    /// Csomoponti atloagos vizkor
    double tt; // travel time
    /// Node rank
    int rank;
    /// WR Graph node index for igraph
    int index;
    /// custom data
    double user1, user2;
    /// Structured EPANET junction demands, retained separately from the
    /// aggregate demand used by the steady-state hydraulic equations.
    std::vector<EpanetDemandComponent> epanet_demand_components;
    double epanet_demand_multiplier;
    /// Original per-node EPANET initial-quality value and its interpretation.
    EpanetInitialQuality epanet_initial_quality;


public:
    /// Konstruktor
    Csomopont(const std::string nev,
              const double h,
              const double fogy,
              const double cl_be,
              const double pressure,
              const double ro,
              const double tt);
    /// Copy Konstruktor
    Csomopont(const Csomopont &csp);
    /// Destruktor
    ~Csomopont();
    // Info
    std::string Info(bool check_if_lonely);
    /// Nyomas beallitasa
    void Set_p(double x)
    {
        p_head = x;
    }
    /// Fogyasztas erteke
    double Get_fogy()
    {
        return fogy;
    }
    /// Fogyaszt�s �rt�ke
    void Set_fogy(double a_fogy)
    {
        fogy = a_fogy;
    }
    /// Az elem neve
    std::string Get_nev()
    {
        return nev;
    }
    /// Get head (p_head)
    double Get_p()
    {
        return p_head;
    }
    /// Geodetikus magassag erteke
    double Get_h()
    {
        return h;
    }
    /// Getting rank
    int Get_rank()
    {
        return rank;
    }
    /// Increasing the rank
    void Inc_rank()
    {
        rank++;
    }
    /// WR Setting index
    void Set_index(int a_index)
    {
        index = a_index;
    }
    int Get_index()
    {
        return index;
    }
    /// Setter for user1
    void Set_user1(double val) {
        user1 = val;
    }
    /// Getter for user2
    void Set_user2(double val) {
        user2 = val;
    }
    /// Setter for user1
    double Get_user1() {
        return user1;
    }
    /// Setter for user2
    double Get_user2() {
        return user2;
    }

    void SetEpanetDemandComponents(
        const std::vector<EpanetDemandComponent> &components,
        double demand_multiplier);
    const std::vector<EpanetDemandComponent> &GetEpanetDemandComponents() const;
    double GetEpanetDemandMultiplier() const;

    void SetEpanetInitialQuality(const EpanetInitialQuality &quality);
    const EpanetInitialQuality &GetEpanetInitialQuality() const;


    /// Inicializ�ci�
    void Ini(int mode, double value);
    /// Double �rt�kek be�ll�t�sa
    void Set_dprop(std::string mit, double value);
    /// Double �rt�kek lek�r�se
    double Get_dprop(std::string mit);


};
#endif
