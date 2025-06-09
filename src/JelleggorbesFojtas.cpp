using namespace std;

#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include "Agelem.h"
#include "JelleggorbesFojtas.h"

JelleggorbesFojtas::JelleggorbesFojtas(const string a_nev,
                                       const string a_cspe_nev,
                                       const string a_cspv_nev,
                                       const double a_ro,
                                       const double a_Aref,
                                       vector<double> a_e,
                                       vector<double> a_zeta,
                                       double a_allas,
                                       const double a_mp) :
        Agelem(a_nev, a_Aref, a_mp, a_ro) {
    //Kotelezo adatok minden Agelemnel:
    tipus = "Jelleggorbes fojtas";
    csp_db = 2;
    cspe_nev = a_cspe_nev;
    cspv_nev = a_cspv_nev;
    // jelleggorbe adatok
    e = a_e;
    zeta = a_zeta;
    allas = a_allas;
    Update_zeta();

}

//--------------------------------------------------------------
void JelleggorbesFojtas::Update_zeta() {
    double VESZT_MIN = 0.0;

    vector<double> allasv(1), vesztv(1);
    allasv.at(0) = allas;
    vesztv = interp(e, zeta, allasv);
    veszt = vesztv.at(0);
    if (fabs(veszt) < VESZT_MIN)
        veszt = copysign(VESZT_MIN,veszt);

//double      copysign( double x, double y );
//Composes a floating point value with the magnitude of x and the sign of y.

// ostringstream strstrm;
// strstrm << endl << nev<<": az aktualis " << fixed << allas
// << "%-os allasnal a vesztesegtenyezo: zeta=" << veszt << endl;
// cout << setprecision(3);
// cout<<strstrm.str();
}


//--------------------------------------------------------------
JelleggorbesFojtas::~JelleggorbesFojtas() {
}

//--------------------------------------------------------------
string JelleggorbesFojtas::Info() {
    ostringstream strstrm;
    strstrm << Agelem::Info();
    strstrm << endl << "  kapcsolodas : " << cspe_nev << "(index:" << cspe_index
            << ") --> " << cspv_nev << "(index:" << cspv_index << ")\n";
    cout << setprecision(3);
    vector<double>::iterator it;
    strstrm << "       adatok : e [%]      = ";
    for (it = e.begin(); it != e.end(); it++)
        strstrm << *(it) << "  ";
    strstrm << "\n";
    strstrm << "                zeta [-]   = ";
    for (it = zeta.begin(); it != zeta.end(); it++)
        strstrm << *(it) << "  ";
    strstrm << "\n";
    strstrm << endl << "       Az aktualis " << fixed << allas
            << "%-os allasnal a vesztesegtenyezo: zeta=" << veszt << endl;
    if (veszt<0){
        strstrm<< endl<<"       Mivel ez az értéek negatív Kv értéket jelent!"<< endl;
    }

    return strstrm.str();
}

//--------------------------------------------------------------
double JelleggorbesFojtas::f(vector<double> x) {
    double ere;
    double pe = x[0] * ro * g;
    double pv = x[1] * ro * g;
    double he = x[2];
    double hv = x[3];

    double loss_coeff = veszt;
    // if veszt<0, it is a Kv value:
    if (veszt<0){
        // Kv = Q(m3/h) / sqrt(dp (bar))
        // Kv = Q(m3/s)*3600 / sqrt(dp *1e5 (Pa))
        // Kv = mp(kg/s)/ro(kg/m3)*3600 / sqrt(dh(m)*ro(kg/m3)*g(m/s^2) *1e5 (Pa))

        loss_coeff = 1.e5/pow(ro,3.)/9.81/pow(veszt,2.)*pow(3600.,2.);
    }

    ere = (pv - pe) / ro / g + (hv - he) + loss_coeff * mp * fabs(mp);

    headloss = (pe - pv) / ro / g + (he - hv);

    return ere;
}

//--------------------------------------------------------------
vector<double> JelleggorbesFojtas::df(vector<double> x) {
    vector<double> ere;

 double loss_coeff = veszt;
    // if veszt<0, it is a Kv value:
    if (veszt<0){
        // Kv = Q(m3/h) / sqrt(dp (bar))
        // Kv = Q(m3/s)*3600 / sqrt(dp *1e5 (Pa))
        // Kv = mp(kg/s)/ro(kg/m3)*3600 / sqrt(dh(m)*ro(kg/m3)*g(m/s^2) *1e5 (Pa))

        loss_coeff = 1.e5/pow(ro,3.)/9.81/pow(veszt,2.)*pow(3600.,2.);
    }

    ere.push_back(-1.0);
    ere.push_back(+1.0);
    ere.push_back(2. * loss_coeff * fabs(mp));
    ere.push_back(0.0);

    return ere;
}

//--------------------------------------------------------------
void JelleggorbesFojtas::Ini(int mode, double value) {
    if (mode == 0)
        mp = 1;
    else
        mp = value;
}

//--------------------------------------------------------------
void JelleggorbesFojtas::Set_dprop(string mit, double mire) {
    if (mit == "position") {
        allas = mire;
        Update_zeta();
    } else if ((mit == "concentration") || (mit == "konc_atlag")) {
        konc_atlag = mire;
    } else {
        cout << endl << endl << "HIBA! JelleggorbesFojtas::Set_dprop(mit), ismeretlen bemenet: mit="
             << mit << endl << endl;
    }
}

//--------------------------------------------------------------
double JelleggorbesFojtas::Get_dprop(string mit) {

    double out = 0.0;
    if (mit == "position")
        out = allas;
    else if (mit == "veszt")
        out = veszt;
    else if (mit == "adzeta")
        out = veszt;
    else if (mit == "mass_flow_rate")
        out = mp;
    else if (mit == "headloss")
        out = headloss;
    else if (mit == "headloss_per_unit_length")
        out = headloss;
    else if ((mit == "concentration") || (mit == "konc_atlag"))
        out = konc_atlag;
    else if ((mit == "length") || (mit == "L"))
        out = 0.5;
    else {
        cout << endl << "HIBA! Jelleggorbes::Get_dprop(mit), ismeretlen bemenet: mit="
             << mit << endl << endl;
        out = 0.0;
    }
    return out;
}