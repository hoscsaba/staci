//! The pipe class.
/*!
Fully filled pipe only!
*/

using namespace std;
#include <string>
#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>
#include "Agelem.h"
#include "Cso.h"

Cso::Cso(const string a_nev, const string a_cspe_nev, const string a_cspv_nev, const double a_ro,
         const double a_L, const double a_D, const double a_erdesseg,
         const double a_cl_k, const double a_cl_w, const double a_mp) :
    Agelem(a_nev, a_D *a_D *pi / 4, a_mp, a_ro)
{
    //Kotelezo adatok minden Agelemnel:
    tipus = "Cso";
    csp_db = 2;
    cspe_nev = a_cspe_nev;
    cspv_nev = a_cspv_nev;
    // hossz, atmero, lambda
    L = a_L;
    D = a_D;
    erdesseg = a_erdesseg;
    f_count = 0;
    cl_k = a_cl_k;
    cl_w = a_cl_w;
    FolyTerf = D * D * pi / 4 * L;
    lambda = 0.02;

}

//--------------------------------------------------------------
//! Sets the friction model
/*!
DW (friction_model_type = 0) -> Darcy Wiesenbach
HW (friction_model_type = 1) -> Hazen-Williams
*/
void Cso::Set_friction_model(string a_fric_type)
{
    if (a_fric_type == "DW")
        friction_model_type = 0;
    else
    {
        if (a_fric_type == "HW")
            friction_model_type = 1;
        else
        {
            cout << endl << "HIBA! Cso::ComputeHeadloss, ismeretlen surlodasi modell (DW|HW) : " << a_fric_type << endl << endl;
        }
    }
}

//--------------------------------------------------------------
Cso::~Cso()
{
}

//--------------------------------------------------------------
string Cso::Info()
{
    ostringstream strstrm;
    strstrm << Agelem::Info();
    strstrm << "\n  kapcsolodas : " << cspe_nev << "(index:" << cspe_index << ") --> "
            << cspv_nev << "(index:" << cspv_index << ")\n";
    strstrm << "       adatok : L=" << L << "[m], D=" << D << "[m], erdesseg=" << erdesseg
            << "[mm], lambda=" << lambda << "[-]";
    strstrm << endl << "                klor lebomlasi allando         : cl_k="
            << cl_k;
    strstrm << endl << "                klor lebomlasi allando a falnal: cl_w="
            << cl_w << endl;

    return strstrm.str();
}

//--------------------------------------------------------------
double Cso::f(vector<double> x)
{
    double ere, tag1;
    double pe = x[0] * ro * g;
    double pv = x[1] * ro * g;
    double he = x[2];
    double hv = x[3];

    tag1 = ro * g * (hv - he);
    ere = pv - pe + tag1 + ComputeHeadloss();
    f_count++;
    lambda = surlodas();
    /* if (fabs(ere / ro / g) > 1e3){
         cout << endl << nev << ", f=" << ere << ", mp=" << mp << ", lambda=" << lambda << ", v=" << (mp / ro / (D * D * pi / 4));
         cout<<endl<<"he+ze="<<x[0]<<"+"<<he<<"="<<x[0]+he;
         cout<<endl<<"hv+zv="<<x[1]<<"+"<<hv<<"="<<x[1]+hv;
     }*/
    /*cin.get();*/
    return ere / ro / g;
}

//--------------------------------------------------------------
vector<double> Cso::df(vector<double> x)
{
    //double pe=x[0]*ro*g;
    //double pv=x[1]*ro*g;
    double he = x[2];
    double hv = x[3];
    vector<double> ere;
    ere.push_back(-ro * g);
    ere.push_back(+ro * g);
    ere.push_back(ComputeHeadlossDerivative());
    ere.push_back(-ro * g * (hv - he));

    for (unsigned int i = 0; i < ere.size(); i++)
        ere.at(i) /= ro * g;
    return ere;
}

//--------------------------------------------------------------
void Cso::Ini(int mode, double value)
{
    if (mode == 0)
        Set_mp(1.);
    else
        Set_mp(value);
}

//--------------------------------------------------------------
//! Sets the friction coefficient lambda
/*!
Based on the actual friction model (DW, HW) computes the friction coefficient.

For DW (friction_model_type = 0) -> Darcy Wiesenbach model parameter 'erdesseg' is pipe surface roughness in mm

For HW (friction_model_type = 1) -> Hazen-Williams model parameter 'erdesseg' is the Hazen-Williamd constant. If the user-supplied value is less tha 10, it is overwritten to 10.

For any of these models, if parameter erdesseg is negative, it is assumed that lambda=-erdesseg
*/
double Cso::surlodas()
{

    double v_min = 0.1;
    double v = mp / ro / (D * D * pi / 4);
    if (fabs(v) < v_min)
        v = v_min;
    double nu = 1e-6;
    double lambda_min = 0.001;
    double lambda_max = 1.0;
    double dp;

    if (friction_model_type == 0) // Darcy-Wiesenbach
    {
        if (erdesseg <= 0)
            lambda = -erdesseg;
        else
        {
            if (f_count >= 0)
            {
                double Re = fabs(v) * D / nu;

                double hiba = 1.0e10, ize = 0.0, lambda_uj = 0.0;
                unsigned int i = 0;
                while ((hiba > 1e-6) && (i < 10))
                {
                    ize = -2.0 * log10(erdesseg / 1000 / D / 3.71 + 2.51 / Re / sqrt(lambda));
                    lambda_uj = 1 / ize / ize;
                    hiba = fabs((lambda - lambda_uj) / lambda);
                    lambda = lambda_uj;
                }
                if (i > 8)
                    cout << endl << endl << "WARNING: " << nev << endl
                         << "\t\t pipe " << nev << " friction factor coefficient iteration #" << i;
            }
            else
                lambda = 0.02;
        }
    }

    if (friction_model_type == 1) // Hazen-Williams, C_factor around 100
    {
        if (erdesseg <= 0)
            lambda = -erdesseg;
        else
        {
            if (f_count >= 0)
            {
                // v=0.849*C_factor*Rh^0.63*s^0.54
                // s= h'/L
                /*double C_factor=100;*/
                //double Rh = D / 4; // A/K=(D^2*pi/4)/(D*pi)=D/4
                double C_factor = erdesseg;

                if (C_factor < 10.)
                {
                    cout << endl << "\tWARNING: " << " pipe " << nev
                         << " friction factor is set to Hazen-Williams but the friction factor is too small (C_HW=" << C_factor << ").";
                    cout << " -> OVERRIDING by C_HW=10.";
                    C_factor = 10.;
                    erdesseg = 10.;
                }

                dp = L / pow(C_factor, 1.85) / pow(D, 4.87) * 7.88 / pow(0.85, 1.85) * pow(fabs(v * Aref), 0.85) * (v * Aref) * ro * g;

                lambda = fabs(dp / (L / D * ro / 2 * v * fabs(v)));

            }
            else
                lambda = 0.02;
        }
    }

    //cout << endl <<nev<< ": lambda=" << lambda << endl;
    //    cin.get();

    if (fabs(lambda) < lambda_min)
    {
        //cout<<endl<<"\t WARNING: "<<nev<<": v="<<v<<"m/s, lambda="<<lambda<<"<"<<lambda_min<<", overriding by "<<lambda_min;
        lambda = lambda_min;
    }

    if (fabs(lambda) > lambda_max)
    {
        //cout<<endl<<"\t WARNING: "<<nev<<": v="<<v<<"m/s, lambda="<<lambda<<">"<<lambda_max<<", overriding by "<<lambda_max;
        lambda = lambda_max;
    }

    return lambda;
}

//--------------------------------------------------------------
double Cso::Get_dprop(string mit)
{

    double out = 0.0;
    if (mit == "Aref")
        out = Aref;
    else if (mit == "lambda")
        out = lambda;
    else if (mit == "diameter")
        out = D;
    else if (mit == "L")
        out = L;
    else if (mit == "Rh")
        out = D / 2.;
    else if (mit == "cl_k")
        out = cl_k;
    else if (mit == "cl_w")
        out = cl_w;
    else if (mit == "erdesseg")
        out = erdesseg;
    else if (mit == "headloss")
        out = fabs(ComputeHeadloss() / ro / g);
    else if (mit == "headloss_per_unit_length")
        out = fabs(ComputeHeadloss() / ro / g / L);
    else if (mit == "mass_flow_rate")
        out = mp;
    else
    {
        cout << endl << "HIBA! Cso::Get_dprop(mit), ismeretlen bemenet: mit="
             << mit << endl << endl;
        out = 0.0;
    }
    return out;
}

//--------------------------------------------------------------
double Cso::Get_dfdmu(string mit)
{

    double out = 0.0;
    if (mit == "diameter")
        out = -5.*surlodas() * L / pow(D, 6) * 8 / ro / pow(pi, 2) * mp * abs(mp);
    //else if (mit == "lambda")
    //    out = lambda;
    else
    {
        cout << endl << "HIBA! Cso::Get_dprop(mit), ismeretlen bemenet: mit="
             << mit << endl << endl;
        out = 0.0;
    }
    return out/ro/g;
}

//--------------------------------------------------------------
void Cso::Set_dprop(string mit, double mire)
{
    if (mit == "diameter")
        D = mire;
    else
    {
        cout << endl << "HIBA! Cso::Set_dprop(mit), ismeretlen bemenet: mit="
             << mit << endl << endl;
    }
}


//--------------------------------------------------------------
//! Computes the head loss in Pa
/*!
dp'=lambda*L/D*ro/2*v*fabs(v)
*/
double Cso::ComputeHeadloss()
{
    double headloss = 0.0;
    double v = mp / ro / Aref;

    headloss = surlodas() * L / D * ro / 2.*v * fabs(v);

    return headloss;
}

//--------------------------------------------------------------
//! Computes the head loss derivative w.r.t. mass flow rate
/*!
dp'=lambda*L/D*ro/2*v*fabs(v)
d dp'/dmp=lambda*L/D*ro/2*1/(ro*A)^2*abs(v)
*/

double Cso::ComputeHeadlossDerivative()
{
    double der;
    der = surlodas() * L / pow(D, 5) * 8 / ro / pow(pi, 2) * 2 * abs(mp);
    return der;
}
