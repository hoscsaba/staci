using namespace std;

#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <math.h>
#include "Agelem.h"
#include "Vegakna.h"

Vegakna::Vegakna(const string a_nev, const string a_csp_nev, const double a_ro, const double Aref,
                 const double a_Hf, const double a_H, const double a_mp, const double a_tt) :
        Agelem(a_nev, Aref, a_mp, a_ro, a_tt) {
    //Kotelezo adatok minden Agelemnel:
    tipus = "Vegakna";
    csp_db = 1;
    cspe_nev = a_csp_nev;
    cspv_nev = "<nincs>";
    // Fenekszint es viszint beallitasa
    Hf = a_Hf;
    H = a_H;
    // Nyomas
    p = (Hf + H) * ro * g;
    //cout<<"\n a_tt="<<a_tt<<", tt_start="<<tt_start<<", tt_end="<<tt_end<<"\n";
    //cin.get();
}

//--------------------------------------------------------------
Vegakna::~Vegakna() {
}

//--------------------------------------------------------------
string Vegakna::Info() {
    ostringstream strstrm;
    strstrm << Agelem::Info();
    strstrm << "\n        kapcsolodas : " << cspe_nev << "(index:" << cspe_index
            << ")\n";
    strstrm << "        adatok : fenek magassag [m] : " << Hf << endl;
    strstrm << "                 vizszint [m]       : " << H << endl;
    return strstrm.str();
}

/*//--------------------------------------------------------------
 double Vegakna::f(vector<double> x) {
 double ere=x.at(1)+x.at(3)-p/ro/g;
 return ere;
 }*/

//--------------------------------------------------------------
double Vegakna::f(vector<double> x) {
    //    double pe=x[0]*ro*g;
    //    double pv=x[1]*ro*g;
    //    double he=x[2];
    //    double hv=x[3];
    double L = 1;
    double D = 0.5;
    double A = D * D * 3.14 / 4;
    double c = 0.02 * L / D / 2. / 9.81 / ro / ro / A / A;
    double ere = x[0] + x[2] + 0. * c * mp * fabs(mp) - (Hf + H);

    //cout<<endl<<" --> POOL: Hf+H="<<(Hf + H)<<" ?= he+ze="<<x[0]<<" + "<<x[2]<<" = "<<x[0] + x[2] << ", Q="<<(mp/1000*3600);

    return ere;
}

//--------------------------------------------------------------
vector<double> Vegakna::df(vector<double> x) {
    vector<double> ere;
    ere.push_back(1.0);
    ere.push_back(0.0);
    double L = 1;
    double D = 0.5;
    double A = D * D * 3.14 / 4;
    double c = 0.02 * L / D / 2. / 9.81 / ro / ro / A / A;
    ere.push_back(0 * 2 * c * fabs(mp));
    //ere.push_back(-p / ro / g);
    ere.push_back(-(Hf + H));

    return ere;
}
/*
 //--------------------------------------------------------------
 vector<double> Vegakna::df(vector<double> x) {
 vector<double> ere;
 ere.push_back(0.0);
 ere.push_back(1.0);
 ere.push_back(0);
 ere.push_back(-p/ro/g);

 return ere;
 }
 */
//--------------------------------------------------------------
void Vegakna::Ini(int mode, double value) {
    if (mode == 0)
        mp = 1;
    else
        mp = value;
}

//--------------------------------------------------------------
void Vegakna::Set_dprop(string mit, double mire) {

    if (mit == "bottom_level")
        Hf = mire;
    else if (mit == "water_level")
        H = mire;
    else if ((mit == "concentration") || (mit == "konc_atlag"))
        konc_atlag = mire;
    else
        cout << endl
             << "HIBA! Vegakna::Set_dprop(mit), ismeretlen bemenet: mit="
             << mit << endl << endl;
}

//--------------------------------------------------------------
double Vegakna::Get_dprop(string mit) {

    double out = 0.0;
    if (mit == "bottom_level")
        out = Hf;
    else if (mit == "mass_flow_rate")
        out = mp;
    else if (mit == "water_level")
        out = H;
    else if (mit == "headloss")
        out = 0.0;
    else if (mit == "headloss_per_unit_length")
        out = 0.0;
    else if ((mit == "concentration") || (mit == "konc_atlag"))
        out = konc_atlag;
    else {
        cout << endl
             << "HIBA! Vegakna::Get_dprop(mit), ismeretlen bemenet: mit="
             << mit << endl << endl;
        out = 0.0;
    }
    return out;
}


