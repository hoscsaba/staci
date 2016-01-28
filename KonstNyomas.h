#include "Agelem.h"

class KonstNyomas: public Agelem
{
private:
    double p;
public:
    KonstNyomas(const string nev,
                const double Aref,
                const string cspnev,
                const double a_ro,
                const double p,
                const double mp,
                const double tt);

    ~KonstNyomas();
    string Info();
    double f(vector<double>);
    vector<double> df(vector<double>);
    void Ini(int mode, double value);
    void Set_dprop(string mit, double mire);
    string GetType()
    {
        return "KonstNyomas";
    }
};
