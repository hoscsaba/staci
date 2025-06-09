#include <iomanip>
#include <string>
#include <vector>
#include "Agelem.h"
#include "Szivattyu.h"
#include "nr.h"

using namespace std;

Szivattyu::Szivattyu(const string a_nev, const string a_cspe_nev,
                     const string a_cspv_nev, const double a_ro, const double Aref, vector<double> a_q,
                     vector<double> a_H, const double a_mp) : Agelem(a_nev, Aref, a_mp, a_ro) {
    //Kotelezo adatok minden Agelemnel:
    tipus = "Szivattyu";
    csp_db = 2;
    cspe_nev = a_cspe_nev;
    cspv_nev = a_cspv_nev;
    // jelleggorbe adatok
    q = a_q;
    H = a_H;
    fokszam = 3; // sum a_i x^i, tehat ha fokszam=2, csak egyenesrol van szo.

    // A terfogataram m3/h-ban erkezik:
    for (unsigned int i = 0; i < q.size(); i++) q.at(i) /= 3600.;

    //cout<<"\n Gorbeillesztes a jellegorbe pontokra:\n";
    Mat_DP A(fokszam, fokszam);
    Vec_DP b(fokszam);
    Vec_INT indx(fokszam);
    DP d;
    //    for (int i=0; i<q.size(); i++) cout<<endl<<"++++ q("<<i<<")="<<q.at(i)<<", H("<<i<<")="<<H.at(i);
    // b eloallitasa:
    for (int i = 0; i < fokszam; i++) {
        b[i] = 0.0;
        for (unsigned int k = 0; k < q.size(); k++)
            b[i] += H.at(k) * pow(q.at(k), i);
        //        cout<<"\nb["<<i<<"]="<<b[i]<<endl<<"A["<<i<<",:]=";
        for (int j = 0; j < fokszam; j++) {
            A[i][j] = 0.0;
            for (unsigned int k = 0; k < q.size(); k++)
                A[i][j] += pow(q.at(k), i) * pow(q.at(k), j);
            //            cout<<"\t"<<A[i][j];
        }
    }
    //cout<<endl;
    NR::ludcmp(A, indx, d);
    NR::lubksb(A, indx, b);

    // A polinom adatainak visszatoltese:
    for (int i = 0; i < fokszam; i++)
        p.push_back(b[i]);

    mer_szorzo = 10;

    if (p.at(1) > 0)
        if (debug_level > 0)
            cout << endl << "\tWarning! PUMP: " << nev << ":  dH/dQ(0)=" << p.at(1) << " > 0 !!!";

    //cout<<endl<<"Az illesztett polinom egyenlete:"<<endl;
    //for (int i=0; i<fokszam; i++) cout<<"\t"<<p[i];
    //cout<<endl;

    /*      cout << endl << endl;
          vector<double> xx;
          vector<double> yy;
          double Qmax = q.at(q.size() - 1);
          cout << endl << nev << ": Qmax=" << Qmax * 3600;
          cin.get();
          double xmin = -0.2*Qmax;
          double xmax = 1.2 * Qmax;
          double dx = (xmax - xmin) / 50;
          double x, y;
          for (unsigned int i = 0; i < 50; i++)
          {
              x = xmin + i * dx;
              y = 0;
              cout << endl << x * 3600 << "\t" << PumpCharCurve(x);
          }
          cout << endl;
          cin.get();
      */
}

//--------------------------------------------------------------
Szivattyu::~Szivattyu() {}

//--------------------------------------------------------------
string Szivattyu::Info() {
    ostringstream strstrm;
    strstrm << Agelem::Info();
    strstrm << endl << "       tipusa : " << tipus;
    strstrm << endl << "  kapcsolodas : " << cspe_nev << "(index:" << cspe_index << ") --> " << cspv_nev << "(index:"
            << cspv_index << ")\n";
    cout << setprecision(3);
    vector<double>::iterator it;
    strstrm << "       adatok : Q [m3/h]= ";
    for (unsigned int i = 0; i < q.size(); i++) strstrm << q.at(i) * 3600 << "  ";
    strstrm << "\n";
    strstrm << "                H [m]   = ";
    for (unsigned int i = 0; i < H.size(); i++) strstrm << H.at(i) << "  ";
    strstrm << "\n";
    strstrm << "       illesztett polinom: H= sum(i=0,fokszam) p_i * q^i" << endl << "\t";
    strstrm << scientific << setprecision(3);
    for (unsigned int i = 0; i < p.size(); i++) strstrm << "p[" << i << "]=" << p.at(i) << "  ";
    strstrm << "\n";
    cin.get();
    return strstrm.str();
}

//! Pump branch equation
/*!
The function evaluates the curve fit by the constructor \sa Szivattyu()
\param x vector<double>
\param x(0)=pstart/ro/g
\param x(1)=pend/ro/g
\param x(2)=zstart
\param x(3)=zend
\return (double) function error (should be zero)
\sa PumpCharCurve()
*/double Szivattyu::f(vector<double> x) {
    double ere;
    double pe = x[0] * ro * g;
    double pv = x[1] * ro * g;
    double he = x[2];
    double hv = x[3];

    ere = (pv - pe) / ro / g - PumpCharCurve(mp / ro) + (hv - he);

    return ere;
}

//! Pump performance curve evaluation
/*!
The function evaluates the curve fit by the constructor \sa Szivattyu()
\param qq (double) flow rate in m^3/s
\return (double) head in m
\sa Szivattyu() and f()
*/
double Szivattyu::PumpCharCurve(double qq) {

    double He = 0.0;
    double qmax = q.at(q.size() - 1);
    double Hmin = H.at(H.size() - 1);
    if (qq < 0) {
        He = -mer_szorzo * p[0] / qmax * qq + p[0];
        if (debug_level > 0)
            cout << endl << "\tWarning! PUMP: " << nev << ": Q=" << (qq * 3600) << " m^3/h < 0, H=" << He
                 << " m, p[0]=H(0)=" << p[0];
    } else {
        if (qq < qmax)
            for (int i = 0; i < fokszam; i++) He += p[i] * pow(qq, i);
        else {
            He = -mer_szorzo * p[0] / qmax * (qq - qmax) + Hmin;
            if (debug_level > 0)
                cout << endl << "\tWarning! PUMP: " << nev << ": Q=" << (qq * 3600)
                     << " m^3/h > Qmax, extrapolating on the performace curve gives H=" << He << " m";
        }
    }
    return He;
}


//--------------------------------------------------------------
vector<double> Szivattyu::df(vector<double> x) {
    vector<double> ere;
    ere.push_back(-1.0);
    ere.push_back(+1.0);

    // Szep megoldas: analitikus derivalt (HCs. 2014.07.30.)
    double der = 0.0;
    double qmax = q.at(q.size() - 1);
    if ((mp < 0) || (mp > qmax * ro))
        der = -mer_szorzo * p[0] / qmax;
    else
        for (int i = 1; i < fokszam; i++) der += p[i] * i * pow(mp / ro, i - 1);
    der /= -ro;

    //--------------------------------
    // HCs. 2014.07.30.
    // EZ FONTOS, NEGATIV MEREDEKSEGU SZIV: KARAKTERISZTIKA ESETEN
    // 0-VAL FELULIRJUK
    if (der < 0)
        der = 0;
    //--------------------------------

    ere.push_back(der);

    ere.push_back(0.0);

    return ere;
}

//--------------------------------------------------------------
void Szivattyu::Ini(int mode, double value) {
    //if (mode==0)
    //mp=fabs(q.at(1)-q.at(q.size()-1))/2*ro;
    //else mp=value;
    if (mode != 0)
        mp = value;
    //mp=fabs(q.at(1)-q.at(q.size()-1))/2*ro;
    //else
}

//--------------------------------------------------------------
void Szivattyu::Set_dprop(string mit, double mire) {
    if ((mit == "concentration") || (mit == "konc_atlag")) {
        konc_atlag = mire;
    } else {
        cout << endl << "HIBA! Szivattyu::Set_dprop(mit), ismeretlen bemenet: mit="
             << mit << endl << endl;
    }
}

//--------------------------------------------------------------
double Szivattyu::Get_dprop(string mit) {

    double out = 0.0;
    if (mit == "Aref")
        out = Aref;
    else if (mit == "mass_flow_rate")
        out = mp;
    else if ((mit == "concentration") || (mit == "konc_atlag"))
        out = konc_atlag;
    else if (mit == "headloss")
        out = abs(PumpCharCurve(mp / ro));
    else if (mit == "headloss_per_unit_length")
        out = abs(PumpCharCurve(mp / ro));
    else {
        cout << endl << "HIBA! Szivattyu::Get_dprop(mit), ismeretlen bemenet: mit="
             << mit << endl << endl;
        out = 0.0;
    }
    return out;
}
