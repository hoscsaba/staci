using namespace std;
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <cmath>
#include <ctime>
#include <stdlib.h>
#include "Csomopont.h"

Csomopont::Csomopont(const string a_nev,
                     const double a_h,
                     const double a_fogy,
                     const double a_cl_be,
                     const double pressure,
                     const double a_ro,
                     const double a_tt)
{
    // A fogyasztas m3/h-ban erkezik...
    ro = a_ro;
    fogy = a_fogy / 3600 * ro;
    h = a_h;
    p_head = 1 * ro * 9.81;
    // azonosito
    nev = a_nev;
    // kapcsoldo elemek elojeles azonositoja
    vector<int> elemek;
    //double PI=3.14159265;
    cl_be = a_cl_be; // gr/s
    p_head = pressure;

    konc_atlag = 0.0;
    tt = a_tt*3600.; // travel time, órában jön, de s-ban fogunk számolni


    // if (a_nev=="NODE_3998554"){
    //     cout<<"\n\n "<<a_nev<<": tt="<<tt<<"\n\n";
    //     cin.get();
    // }
}

//--------------------------------------------------------------
Csomopont::~Csomopont()
{
}

//--------------------------------------------------------------
string Csomopont::Info()
{
    ostringstream strstrm;
    strstrm << "\n Csomopont neve: " << nev;
    strstrm << "\n       magassag: " << h << " m";
    strstrm << "\n           head: " << p_head << " m (=p[Pa]/ro/g)";
    strstrm << "\n       pressure: " << p_head *ro * 9.81 << " Pa";
    strstrm << "\n        suruseg: " << ro << " kg/m3";
    strstrm << "\n        elvetel: " << fogy << " kg/s = " << fogy * 3600 / ro << " m3/h";
    strstrm << "\n         vizkor: " << tt << " s = " << tt/ 60. << " min";
    strstrm << "\n    erkezo agak: ";
    for (vector<int>::iterator it = ag_be.begin(); it != ag_be.end(); it++)
        strstrm << *it << " ";
    strstrm << "\n    indulo agak: ";
    for (vector<int>::iterator it = ag_ki.begin(); it != ag_ki.end(); it++)
        strstrm << *it << " ";
    strstrm << endl;
    strstrm << "   beadagolt klor: " << cl_be << " gr/m3\n";

    if ((ag_be.size() + ag_ki.size()) == 0)
    {
        strstrm << "\n!!! PANIC !!! Lonely node: " << nev
        << " !!!\n";
        cout << strstrm.str();
        exit(-1);
    }

    return strstrm.str();
}

//--------------------------------------------------------------
void Csomopont::Ini(int mode, double value)
{
    if (mode == 0)
        p_head = 300. - h;
    else
        p_head = value - h;
}

//--------------------------------------------------------------
void Csomopont::Set_dprop(string mit, double value)
{
    if (mit == "konc_be")
        cl_be = value;
    if (mit == "demand")
        // m3/h-ban kell megadni, de belul kg/s-ban taroljuk
        fogy = value / 3600 * ro;
    if (mit == "konc_atlag")
        konc_atlag = value;
    if (mit == "tt")
        tt = value;
}

//--------------------------------------------------------------
double Csomopont::Get_dprop(string mit)
{

    bool megvan = false;
    double outdata = 0.0;

    if (mit == "cl_be")
    {
        megvan = true;
        outdata = cl_be;
    }

    if (mit == "konc_atlag")
    {
        megvan = true;
        outdata = konc_atlag;
    }

    if (mit == "demand")
    {
        megvan = true;
        // kg/s-ban taroljuk, de m3/h-ban adjuk vissza
        outdata = fogy*3600/ro;
    }

    if (mit == "head")
    {
        megvan = true;
        outdata = p_head/ro/9.81;
    }
    if (mit == "pressure")
    {
        megvan = true;
        outdata = p_head;
    }

    if (mit == "ro")
    {
        megvan = true;
        outdata = ro;
    }

    if (mit == "tt")
    {
        megvan = true;
        outdata = tt;
    }

    if (!megvan)
    {
        cout << endl << endl << "Csomopont::Get_dprop() hibas arg.:" << mit;
        cout << ", helyes ertekek: konc_be|konc_atlag|consumption|pressure" << endl << endl;
    }
    return outdata;
}