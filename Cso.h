#include "Agelem.h"

class Cso : public Agelem
{
private:
    double L, D, lambda;
    double erdesseg; // csofal abszolut erdesseg
    int f_count;
    double cl_k, cl_w;
    double headloss; // Megoszló, teljes hossz mentén és m-enként
    int friction_model_type; // 0 - Darcy-Weisbach (DW az adatfajlban), 1 - Hazen-Williams (HW)
public:
    Cso(const string nev, const string a_cspe_nev, const string a_cspv_nev, const double a_ro,
        const double L, const double D, const double lambda,
        const double cl_k, const double cl_w, const double mp);
    ~Cso();
    string Info();
    double f(vector<double>);
    vector<double> df(vector<double>);
    void Ini(int mode, double value);
    double surlodas();
    double Get_dprop(string mit);
    double Get_dfdmu(string mit);
    void Set_dprop(string mit, double mire);
    string GetType()
    {
        return "Cso";
    }    
    double ComputeHeadloss();
    double ComputeHeadlossDerivative();
    void Set_friction_model(string friction_model);
};
