using namespace std;

#include <cmath>
#include <vector>
#include <iostream>
#include <sstream>
#include "Agelem.h"
#include "KonstNyomas.h"

KonstNyomas::KonstNyomas(
        const string a_nev,
        const double a_Aref,
        const string a_csp_nev,
        const double a_ro,
        const double a_p,
        const double a_mp,
        const double a_tt) :

        Agelem(a_nev, a_Aref, a_mp, a_ro, a_tt) {
    //Kotelezo adatok minden Agelemnel:
    tipus = "Konstans nyomas";
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
string KonstNyomas::Info() {
    ostringstream strstrm;
    strstrm << Agelem::Info();
    strstrm << endl << "  kapcsolodas : " << cspv_nev << "(index:" << cspe_index << ")\n";
    return strstrm.str();
}

//--------------------------------------------------------------
double KonstNyomas::f(vector<double> x) {
    double ere = x[0] - p / ro / g;
    return ere;
}

//--------------------------------------------------------------
vector<double> KonstNyomas::df(vector<double> x) {
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
void KonstNyomas::Set_dprop(string mit, double mire) {

    if ((mit == "concentration") || (mit == "konc_atlag")) {
        konc_atlag = mire;
    } else {
        cout << endl << "HIBA! KonstNyomas::Set_dprop(mit), ismeretlen bemenet: mit="
             << mit << endl << endl;
    }
}