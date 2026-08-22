#ifndef KONSTNYOMAS_H
#define KONSTNYOMAS_H

#include "Agelem.h"

struct EpanetHeadPattern
{
    bool specified = false;
    double base_head_m = 0.0;
    string pattern_id;
    vector<double> pattern_values;
};

class KonstNyomas: public Agelem
{
private:
    double p;
    EpanetHeadPattern epanet_head_pattern;
public:
    KonstNyomas(const string &nev,
                double Aref,
                const string &cspnev,
                double a_ro,
                double p,
                double mp,
                double tt);

    ~KonstNyomas() override;
    string Info() override;
    double f(const vector<double> &state) override;
    vector<double> df(const vector<double> &state) override;
    void Ini(int mode, double value) override;
    void Set_dprop(const string &property, double value) override;
    string_view GetType() const noexcept override
    {
        return "KonstNyomas";
    }
    double Get_dprop(const string &property) override;
    void SetEpanetHeadPattern(const EpanetHeadPattern &head_pattern);
    const EpanetHeadPattern &GetEpanetHeadPattern() const;
};

#endif
