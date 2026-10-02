#include <iomanip>
#include <string>
#include <vector>
#include "Agelem.h"
#include "Szivattyu.h"
#include "nr.h"

using namespace std;

Szivattyu::Szivattyu(const string &a_nev, const string &a_cspe_nev,
                     const string &a_cspv_nev, double a_ro, double Aref,
                     const vector<double> &a_q, const vector<double> &a_H,
                     double a_mp) : Agelem(a_nev, Aref, a_mp, a_ro),
                                    operating_speed(1.0) {
    //Kotelezo adatok minden Agelemnel:
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
\param x State vector containing the upstream and downstream pressure heads
and elevations in positions 0 through 3.
\return (double) function error (should be zero)
\sa PumpCharCurve()
*/bool Szivattyu::reverse_closed(const vector<double>& x) {
    if(metadata.definition!="HEAD" || metadata.head_curve_points.empty()) return false;
    return mp < -2.8316846592e-8*ro || (mp<=0 && x[1]+x[3]-x[0]-x[2] > PumpCharCurve(0)+0.0001524);
}

double Szivattyu::f(const vector<double> &x) {
    hydraulic_closed_=reverse_closed(x);
    if (!enabled || operating_speed <= 0.0 || reverse_closed(x))
        return mp;
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
    if (operating_speed <= 0.0)
        return 0.0;
    return operating_speed * operating_speed *
           BasePumpCharCurve(qq / operating_speed);
}

double Szivattyu::BasePumpCharCurve(double qq) {
    if (metadata.definition == "HEAD" && !metadata.head_curve_points.empty()) {
        auto points = metadata.head_curve_points;
        if (points.size() == 1) {
            const double q1 = points[0].first, h1 = points[0].second;
            points = {{0.0, 1.33334*h1}, {q1, h1}, {2*q1, 0.0}};
        }
        if (points.size() == 3 && points[0].first == 0.0) {
            const double h0 = points[0].second;
            const double exponent = std::log((h0-points[2].second)/(h0-points[1].second)) /
                std::log(points[2].first/points[1].first);
            const double resistance = (h0-points[1].second)/std::pow(points[1].first, exponent);
            return h0 - resistance * std::pow(std::max(0.0, qq), exponent);
        }
        std::size_t upper = 1;
        while (upper+1 < points.size() && qq > points[upper].first) ++upper;
        const auto a = points[upper-1], b = points[upper];
        return a.second + (qq-a.first)*(b.second-a.second)/(b.first-a.first);
    }

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
vector<double> Szivattyu::df(const vector<double> &x) {
    vector<double> ere;
    if (!enabled || operating_speed <= 0.0 || reverse_closed(x)) {
        ere.push_back(0.0);
        ere.push_back(0.0);
        ere.push_back(1.0);
        ere.push_back(0.0);
        return ere;
    }
    ere.push_back(-1.0);
    ere.push_back(+1.0);

    // Szep megoldas: analitikus derivalt (HCs. 2014.07.30.)
    double der = -operating_speed *
                 BasePumpCharCurveDerivative(mp / ro / operating_speed) / ro;

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
        mp = enabled && operating_speed > 0.0 ? value : 0.0;
    //mp=fabs(q.at(1)-q.at(q.size()-1))/2*ro;
    //else
}

//--------------------------------------------------------------
void Szivattyu::Set_dprop(const string &mit, double mire) {
    if ((mit == "concentration") || (mit == "konc_atlag")) {
        konc_atlag = mire;
    } else if (mit == "speed") {
        SetOperatingSpeed(mire);
        metadata.base_speed = mire;
    } else if (mit == "status") {
        Set_enabled(mire != 0.0);
    } else {
        cout << endl << "HIBA! Szivattyu::Set_dprop(mit), ismeretlen bemenet: mit="
             << mit << endl << endl;
    }
}

//--------------------------------------------------------------
double Szivattyu::Get_dprop(const string &mit) {

    double out = 0.0;
    if (mit == "Aref")
        out = Aref;
    else if (mit == "mass_flow_rate")
        out = mp;
    else if (mit == "speed")
        out = operating_speed;
    else if (mit == "base_speed")
        out = metadata.base_speed;
    else if (mit == "speed_pattern_length")
        out = static_cast<double>(metadata.speed_pattern_values.size());
    else if (mit == "efficiency_curve_points")
        out = static_cast<double>(metadata.efficiency_curve_points.size());
    else if (mit == "head_curve_points")
        out = static_cast<double>(metadata.head_curve_points.size());
    else if (mit == "status")
        out = enabled && operating_speed > 0.0 && !hydraulic_closed_ ? 1.0 : 0.0;
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

//--------------------------------------------------------------
double Szivattyu::BasePumpCharCurveDerivative(double qq) const {
    if (metadata.definition == "HEAD" && !metadata.head_curve_points.empty()) {
        auto points = metadata.head_curve_points;
        if (points.size() == 1) {
            const double q1 = points[0].first, h1 = points[0].second;
            points = {{0.0, 1.33334*h1}, {q1, h1}, {2*q1, 0.0}};
        }
        if (points.size() == 3 && points[0].first == 0.0) {
            const double h0 = points[0].second;
            const double exponent = std::log((h0-points[2].second)/(h0-points[1].second)) /
                std::log(points[2].first/points[1].first);
            return -(h0-points[1].second)/std::pow(points[1].first, exponent) * exponent *
                std::pow(std::max(1e-12, qq), exponent-1);
        }
        std::size_t upper = 1;
        while (upper+1 < points.size() && qq > points[upper].first) ++upper;
        return (points[upper].second-points[upper-1].second)/(points[upper].first-points[upper-1].first);
    }
    const double qmax = q.at(q.size() - 1);
    if (qq < 0.0 || qq > qmax)
        return -mer_szorzo * p[0] / qmax;
    double derivative = 0.0;
    for (int i = 1; i < fokszam; ++i)
        derivative += p[i] * i * pow(qq, i - 1);
    return derivative;
}

void Szivattyu::SetEpanetPumpMetadata(const EpanetPumpMetadata &value) {
    metadata = value;
    SetOperatingSpeed(metadata.base_speed);
}

void Szivattyu::SetOperatingSpeed(double value) {
    if (value < 0.0)
        throw invalid_argument("EPANET pump speed cannot be negative.");
    operating_speed = value;
}
