using namespace std;

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <cmath>
#include <ctime>
#include <stdlib.h>
#include "Agelem.h"


Agelem::Agelem(const string a_nev, const double a_Aref, const double a_mp, const double a_ro) {
    // tipus
    string tipus;
    // tomegaram es suruseg
    mp = a_mp;

    ro = a_ro;

    // azonosito
    nev = a_nev;
    Aref = a_Aref;
    // Csomopontok vektora
    cspe_index = -1, cspv_index = -1;
    string cspe_nev = "nincs_cspe_nev", cspv_nev = "nincs_cspv_nev";
    FolyTerf = 0;
    debug_level = 0;

    tt_start = 0.;
    tt_end = 0.;

    if (fabs(a_ro) < 1.0e-3)
        error("Agelem constructor", "Density (ro) is zero up to machine precision!");

    force_more_iter = false;
    update_diameter = false;

}

//--------------------------------------------------------------
Agelem::Agelem(const string a_nev, const double a_Aref, const double a_mp, const double a_ro, const double a_tt) {
    // tipus
    string tipus;
    // tomegaram es suruseg
    mp = a_mp;

    ro = a_ro;

    // azonosito
    nev = a_nev;
    Aref = a_Aref;
    // Csomopontok vektora
    cspe_index = -1, cspv_index = -1;
    string cspe_nev = "nincs_cspe_nev", cspv_nev = "nincs_cspv_nev";
    FolyTerf = 0;
    debug_level = 0;

    tt_start = a_tt * 3600; // Az adatfajlbol oraban olvassuk ki, oraban
    tt_end = a_tt * 3600; // Az adatfajlbol oraban olvassuk ki, oraban

    if (fabs(a_ro) < 1.0e-3)
        error("Agelem constructor", "Density (ro) is zero up to machine precision!");

}

//--------------------------------------------------------------
Agelem::~Agelem() {
}

//--------------------------------------------------------------
string Agelem::Info() {
    ostringstream strstrm;
    strstrm << "\n Agelem neve  : " << nev;
    strstrm << "\n        tipusa: " << tipus;
    strstrm << "\n        ro    : " << ro << " [kg/m^3]";
    strstrm << "\n        Aref  : " << Aref << " [m^2]";
    strstrm << "\n        mp    : " << mp / ro * 3600 << " [m3/h]";
    return strstrm.str();
}

//--------------------------------------------------------------
void Agelem::add_csp(const int a_cspe_index, const int a_cspv_index) {
    cspe_index = a_cspe_index;
    cspv_index = a_cspv_index;
}

//--------------------------------------------------------------
vector<double> Agelem::Get_res(string mit) {
    vector<double> x;
    return x;
}

//--------------------------------------------------------------
vector<double> Agelem::interp(vector<double> x, vector<double> y, vector<double> xg) {
    vector<double> yg;
    double xp, xn, yp, yn;

    double xmin = 1e100, xmax = -1e100;
    for (unsigned int i = 0; i < x.size(); i++) {
        if (x.at(i) < xmin)
            xmin = x.at(i);
        if (x.at(i) > xmax)
            xmax = x.at(i);
    }
    if (abs(xmin) < 1e-8)
        xmin = 0.0;

    for (unsigned int i = 0; i < xg.size(); i++) {
        if (xg.at(i) < xmin || xg.at(i) > xmax) {
            cout << endl << endl
                 << "!!!element/basic/interp interp hiba!!! Interpolacio kiserlet a tartmanyon kivulre:";
            cout << endl << "\t xmin=" << xmin << " <? " << xg.at(i) << " <? xmax=" << xmax
                 << endl;
        }
        //      if (xg.at(i)<xmin*0.999)
        //          xg.at(i)=xmin;
        //      if (xg.at(i)>xmax*1.001)
        //          xg.at(i)=xmax;
    }

    for (unsigned int i = 0; i < xg.size(); i++) {
        unsigned int j = 0;
        bool megvan = false;
        while ((!megvan) && (j < x.size() - 1)) {
            double ize = (xg.at(i) - x.at(j)) * (xg.at(i) - x.at(j + 1));
            //cout<<endl<<"xg="<<xg.at(i)<<"  x.at("<<j<<")="<<x.at(j)<<"  x.at("<<j+1<<")="<<x.at(j+1)<<"  y.at("<<j<<")="<<y.at(j)<<"  y.at("<<j+1<<")="<<y.at(j+1)<<" ize="<<ize;
            if (ize < 1e-10)
                megvan = true;
            else
                j++;
        }
        if (j == x.size() - 1)
            j--;
        xp = x.at(j);
        xn = x.at(j + 1);
        yp = y.at(j);
        yn = y.at(j + 1);
        yg.push_back(yp + (xg.at(i) - xp) * (yn - yp) / (xn - xp));
        //cout<<endl<<"==>xg="<<xg.at(i)<<"  x.at(j)="<<x.at(j)<<"  x.at(j+1)="<<x.at(j+1)<<"  y.at(j)="<<y.at(j)<<"  y.at(j+1)="<<y.at(j+1)<<" yi="<<yg.at(i);
        //int int1; cin>>int1;
    }
    /*
     cout<<endl<<endl<<"Csatorna::interp, az eredeti vektorok:";
     for (int i=0; i<x.size(); i++) cout<<endl<<scientific<<"\t x ="<<x.at(i)<<" y ="<<y.at(i);
     cout<<endl<<"Az interpolalt ertekek ";
     for (int i=0; i<xg.size(); i++) cout<<endl<<scientific<<"\t xg="<<xg.at(i)<<" yg="<<yg.at(i);
     cout<<endl<<"Csatorna::interp: kesz (kerek egy egesz erteket...)"<<endl<<endl;
     int int1; cin>>int1;
     */
    return yg;
}

//--------------------------------------------------------------
void Agelem::set_up_grid(double a_konc, vector<double> a_vel, double a_cL) {
    //konc.clear(); vel.clear();
    // Vizminoseg adatok:
    for (unsigned int i = 0; i < a_vel.size(); i++)
        konc.push_back(a_konc);
    for (unsigned int i = 0; i < a_vel.size(); i++)
        vel.push_back(a_vel.at(i));
    cL = a_cL;
    cT = cL / fabs(mean(vel));
    cdt = cT / vel.size();
}

//--------------------------------------------------------------
string Agelem::show_grid(double ido) {
    ostringstream strstrm;
    strstrm << endl << "A " << nev << " agelem seb.- es konc.eloszlasa  t=" << ido
            << "s-ban, cL=" << cL << "m, cT=" << cT << "s, cdt=" << cdt << "s\n";
    strstrm << "\n\tv= ";
    for (unsigned int i = 0; i < vel.size(); i++)
        strstrm << scientific << showpos << setprecision(2) << vel.at(i) << " ";
    strstrm << "\n\tc= ";
    for (unsigned int i = 0; i < vel.size(); i++)
        strstrm << konc.at(i) << " ";
    strstrm << "\n";
    return strstrm.str();
}

//--------------------------------------------------------------
double Agelem::mean(vector<double> x) {
    double mean = 0;
    for (unsigned int i = 0; i < x.size(); i++)
        mean += x.at(i);
    return mean / x.size();
}

//--------------------------------------------------------------
void Agelem::error(string fv, string msg) {
    ostringstream strstrm;
    strstrm.str("");
    strstrm << "\n\n******** ERROR *********";
    strstrm << "\n\telement name: " << nev;
    strstrm << "\n\tmethod      : " << fv;
    strstrm << "\n\tmessage     : " << msg << "\n\n";
    logfile_write(strstrm.str(), 0);
    cout << strstrm.str();
    exit(0);
}

//--------------------------------------------------------------
void Agelem::logfile_write(string msg, int msg_debug_level) {
    if (debug_level >= msg_debug_level) {
        ofstream outfile(out_file.c_str(), ios::app);
        outfile << msg;
        outfile.close();
    }
}

void Agelem::SetLogFile()
{
    out_file = nev + ".out";
    ofstream outputFile;
    outputFile.open(out_file.c_str());
    outputFile.close();
};