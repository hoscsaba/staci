using namespace std;
#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <cmath>
#include "Agelem.h"
#include "VisszacsapoSzelep.h"

VisszacsapoSzelep::VisszacsapoSzelep(const string a_nev,
                                     const string a_cspe_nev,
                                     const string a_cspv_nev,
                                     const double a_ro,
                                     const double Aref,
                                     const double a_dzeta_e,
                                     const double a_dzeta_v,
                                     const double a_mp) :
    Agelem(a_nev, Aref, a_mp, a_ro)
{
    //Kotelezo adatok minden Agelemnel:
    tipus = "VisszacsapoSzelep";
    csp_db = 2;
    cspe_nev = a_cspe_nev;
    cspv_nev = a_cspv_nev;
    // Visszafolyás irányban az ellenállás
    dzeta_e = a_dzeta_e;
    dzeta_v = a_dzeta_v;
}

//--------------------------------------------------------------
VisszacsapoSzelep::~VisszacsapoSzelep()
{
}

//--------------------------------------------------------------
string VisszacsapoSzelep::Info()
{
    ostringstream strstrm;
    strstrm << Agelem::Info();
    strstrm << "\n     tipusa : " << tipus;
    strstrm << "\nkapcsolodas : " << cspv_nev << "(index:" << cspv_index << ")\n";
    strstrm << "\n\t átfolyási tényezõ folyásirányban : " << dzeta_e;
    strstrm << "\n\t átfolyási tényezõ zárás iráynyban: " << dzeta_v << "\n";
    return strstrm.str();
}

//--------------------------------------------------------------
double VisszacsapoSzelep::f(vector<double> x)
{
    double ere, veszt;
    double pe = x[0] * ro * g;
    double pv = x[1] * ro * g;
    double he = x[2];
    double hv = x[3];

    if (mp > 0)
        veszt = dzeta_e;
    else
        veszt = dzeta_v;

    ere = (pv - pe) / ro / g + (hv - he) + veszt * mp * fabs(mp);

    /*if (nev == "CHECKVALVE676")
    {
        cout << Info();
        cout << endl << nev << ", f=" << ere << ", mp=" << mp;
        cout << endl << "he+ze=" << x[0] << "+" << he << "=" << x[0] + he;
        cout << endl << "hv+zv=" << x[1] << "+" << hv << "=" << x[1] + hv;
        cout << endl << "veszt=" << veszt;
    }*/
    return ere;
}

//--------------------------------------------------------------
vector<double> VisszacsapoSzelep::df(vector<double> x)
{
    vector<double> ere;
    if (mp > 0)
    {
        ere.push_back(-1.0);
        ere.push_back(1.0);
        ere.push_back(dzeta_e * 2 * mp);
        ere.push_back(0);
    }
    else
    {
        ere.push_back(-1.0);
        ere.push_back(1.0);
        ere.push_back(-dzeta_v * 2 * mp);
        ere.push_back(0);

    }

    return ere;
}

//--------------------------------------------------------------
void VisszacsapoSzelep::Ini(int mode, double value)
{
    if (mode == 0)
        mp = 1;
    else
        mp = value;
}

//--------------------------------------------------------------
void VisszacsapoSzelep::Set_dprop(string mit, double mire)
{
    //    if (mit=="diameter")
    //      D=mire;
    //    else
    //      {
    cout << endl << "HIBA! VisszacsapoSzelep::Set_dprop(mit), ismeretlen bemenet: mit="
         << mit << endl << endl;
    //      }
}

double VisszacsapoSzelep::Get_dprop(string mit) {
    double out = 0.0;
    // if (mit == "Aref")
    //   out = Aref;
    // else if (mit == "lambda")
    //   out = lambda;
    // else {
    cout << endl
         << "HIBA! Cso::Get_dprop(mit), ismeretlen bemenet: mit=" << mit << endl
         << endl;
    cout << endl << "Name of VisszacsapoSzelep: " << nev << endl;
    cin.get();
    // }
    return out;
}