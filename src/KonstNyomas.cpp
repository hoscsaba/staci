#include <cmath>
#include <vector>
#include <iostream>
#include <sstream>
#include "Agelem.h"
#include "KonstNyomas.h"

using namespace std;

KonstNyomas::KonstNyomas(
    const string &a_nev,
    double a_Aref,
    const string &a_csp_nev,
    double a_ro,
    double a_p,
    double a_mp,
    double a_tt) :

    Agelem(a_nev, a_Aref, a_mp, a_ro, a_tt) {
    //Kotelezo adatok minden Agelemnel:
    csp_db = 1;

    /*cspe_nev="<nincs>";
    cspv_nev=a_csp_nev;
    */
    cspv_nev = "<nincs>";
    cspe_nev = a_csp_nev;
    // nyomas
    p = a_p;
}

//--------------------------------------------------------------
KonstNyomas::~KonstNyomas() {
}

//--------------------------------------------------------------
void KonstNyomas::SetEpanetHeadPattern(const EpanetHeadPattern &head_pattern) {
    epanet_head_pattern = head_pattern;
}

//--------------------------------------------------------------
const EpanetHeadPattern &KonstNyomas::GetEpanetHeadPattern() const {
    return epanet_head_pattern;
}

//--------------------------------------------------------------
string KonstNyomas::Info() {
    ostringstream strstrm;
    strstrm << Agelem::Info();
    strstrm << endl << "  kapcsolodas : " << cspe_nev << "(index:" << cspe_index << ")\n";
    return strstrm.str();
}

//--------------------------------------------------------------
double KonstNyomas::f(const vector<double> &x) {
    double ere = x[0] - p / ro / g;
    return ere;
}

//--------------------------------------------------------------
vector<double> KonstNyomas::df(const vector<double> &x) {
    vector<double> ere;
    ere.push_back(1);
    ere.push_back(0);
    ere.push_back(0);
    ere.push_back(-p / ro / g);

    return ere;
}

//--------------------------------------------------------------
void KonstNyomas::Ini(int mode, double value) {
    if (mode == 0)
        mp = 0.01;
    else
        mp = value;
}

//--------------------------------------------------------------
void KonstNyomas::Set_dprop(const string &mit, double mire) {

    if (mit == "head") {
        p = mire * ro * g;
    } else if ((mit == "concentration") || (mit == "konc_atlag")) {
        konc_atlag = mire;
    } else {
        cout << endl << "HIBA! KonstNyomas::Set_dprop(mit), ismeretlen bemenet: mit="
             << mit << endl << endl;
    }
}

double KonstNyomas::Get_dprop(const string &mit) {
    double out = 0.0;
    if (mit == "head")
        out = p / ro / g;
    else if (mit == "headloss")
        out = 0.;
    else if (mit == "headloss_per_unit_length")
        out = 0.;
    else if ((mit == "concentration") || (mit == "konc_atlag"))
        out = konc_atlag;
    else {
        cout << endl
             << "HIBA! Cso::Get_dprop(mit), ismeretlen bemenet: mit=" << mit << endl
             << endl;
        cout << endl << "Name of KonstNyomas: " << nev << endl;
        cin.get();
    }
    return out;
}
