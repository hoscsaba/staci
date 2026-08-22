#include "Agelem.h"

class Cso : public Agelem, public ParameterSensitivity, public FrictionModelConfigurable
{
private:
    double L, D, lambda;
    double minor_loss;
    double erdesseg; // csofal abszolut erdesseg
    int f_count;
    double cl_k, cl_w;
    double headloss; // Megoszló, teljes hossz mentén és m-enként
    int friction_model_type; // 0 - Darcy-Weisbach (DW az adatfajlban), 1 - Hazen-Williams (HW)
    bool check_valve;
    bool CheckValveClosed(const vector<double> &state) const;
public:
    Cso(const string &nev, const string &a_cspe_nev, const string &a_cspv_nev, const double a_ro,
        const double L, const double D, const double lambda,
        const double cl_k, const double cl_w, const double mp);
    ~Cso() override;
    string Info() override;
    double f(const vector<double> &state) override;
    vector<double> df(const vector<double> &state) override;
    void Ini(int mode, double value) override;
    double surlodas();
    double Get_dprop(const string &property) override;
    double Get_dfdmu(const string &property) override;
    void Set_dprop(const string &property, double value) override;
    string_view GetType() const noexcept override
    {
        return "Cso";
    }    
    double ComputeHeadloss();
    double ComputeHeadlossDerivative();
    void Set_friction_model(const string &friction_model) override;
    void SetCheckValve(bool value) { check_valve = value; }
    bool IsCheckValve() const { return check_valve; }
};
