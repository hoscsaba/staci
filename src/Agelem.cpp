#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <cmath>
#include <ctime>
#include <limits>
#include <stdexcept>
#include <stdlib.h>
#include "Agelem.h"

using namespace std;

Agelem::Agelem(const string &a_nev, double a_Aref, double a_mp, double a_ro)
    : Agelem(a_nev, a_Aref, a_mp, a_ro, 0.0) {}

//--------------------------------------------------------------
Agelem::Agelem(const string &a_nev, double a_Aref, double a_mp, double a_ro, double a_tt)
    : mp(a_mp),
      ro(a_ro),
      Aref(a_Aref),
      nev(a_nev),
      cspe_index(-1),
      cspv_index(-1),
      cspe_nev("nincs_cspe_nev"),
      cspv_nev("nincs_cspv_nev"),
      csp_db(0),
      FolyTerf(0.0),
      head_loss(0.0),
      out_file(),
      debug_level(0),
      tt_start(a_tt * 3600.0),
      tt_end(a_tt * 3600.0),
      user1(0.0),
      user2(0.0),
      enabled(true),
      konc(),
      vel(),
      konc_atlag(0.0),
      ido(0.0),
      cL(1.0),
      cT(0.0),
      cdt(0.0),
      force_more_iter(false),
      update_diameter(false) {
    if (fabs(a_ro) < 1.0e-3)
        error("Agelem constructor", "Density (ro) is zero up to machine precision!");
}

//--------------------------------------------------------------
string Agelem::Info() {
    ostringstream strstrm;
    strstrm << "\n Agelem neve  : " << nev;
    strstrm << "\n        type  : " << GetType();
    strstrm << "\n        ro    : " << ro << " [kg/m^3]";
    strstrm << "\n        Aref  : " << Aref << " [m^2]";
    strstrm << "\n        mp    : " << mp / ro * 3600 << " [m3/h]";
    return strstrm.str();
}

//--------------------------------------------------------------
void Agelem::add_csp(int a_cspe_index, int a_cspv_index) {
    cspe_index = a_cspe_index;
    cspv_index = a_cspv_index;
}

//--------------------------------------------------------------
vector<double> Agelem::interp(const vector<double> &x, const vector<double> &y,
                              const vector<double> &xg) {
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
void Agelem::set_up_grid(double a_konc, const vector<double> &a_vel, double a_cL) {
    if (a_vel.empty())
        throw std::invalid_argument("Agelem::set_up_grid requires at least one velocity value.");
    konc.clear();
    vel.clear();
    // Vizminoseg adatok:
    for (unsigned int i = 0; i < a_vel.size(); i++)
        konc.push_back(a_konc);
    for (unsigned int i = 0; i < a_vel.size(); i++)
        vel.push_back(a_vel.at(i));
    cL = a_cL;
    const double mean_velocity = fabs(mean(vel));
    if (mean_velocity <= 1.0e-12) {
        cT = std::numeric_limits<double>::infinity();
        cdt = std::numeric_limits<double>::infinity();
    } else {
        cT = cL / mean_velocity;
        cdt = cT / vel.size();
    }
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
double Agelem::mean(const vector<double> &x) {
    double mean = 0;
    for (unsigned int i = 0; i < x.size(); i++)
        mean += x.at(i);
    return mean / x.size();
}

//--------------------------------------------------------------
void Agelem::error(const string &fv, const string &msg) {
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
void Agelem::logfile_write(const string &msg, int msg_debug_level) {
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
