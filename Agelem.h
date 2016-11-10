#ifndef AGELEM_H
#define AGELEM_H
#include <iostream>
#include <vector>
using namespace std;

class Agelem
{
protected:
    /// Tomegaram, [kg/s]
    double mp;
    /// Suruseg, [kg/m3]
    double ro;
    /// Referencia keresztmetszet, [m2], a sebesseg kiszamitasahoz hasznaljuk
    double Aref;
    /// Pi
    static const double pi = 3.1416;
    /// g
    static const double g = 9.81;
    /// Az agelem neve
    string nev;
    /// Az agelem tipusa, pl. Cso, Csatorna, stb.
    string tipus;
    /// Az elejen es vegen levo csompont indexe
    int cspe_index, cspv_index;
    /// Az elejen es vegen levo csompont neve
    string cspe_nev, cspv_nev;
    int csp_db;
    /// Folyadekterfogat
    double FolyTerf;
    /// Nyomaseses
    double head_loss;
    /// Standard error message
    void error(string fv, string msg);
    string out_file;
    int debug_level;
    /// vizkor az elejen
    double tt_start;
    /// vizkor a vegen
    double tt_end;
    

public:
    void logfile_write(string msg, int msg_debug_level);
    /// Konstruktor
    Agelem(const string nev, const double Aref, const double mp, const double ro);
    /// Konstruktor with travel time
    Agelem(const string nev, const double Aref, const double mp, const double ro, const double tt);
    /// Destruktor
    virtual ~Agelem();
    /// Informacio
    virtual string Info();
    /// Csomopont beallitasa
    virtual void add_csp(const int cspe_index, const int cspv_index);
    /// Az agegyenlet erteke, nullara rendezve, v.o.m.-ben
    virtual double f(vector<double>) = 0;
    /// Jacobi: df/dhe, df/dhv, df/dmp, konstans tag
    virtual vector<double> df(vector<double>) = 0;
    /// Inicializacio, mode=0 -> automatikus, mode=1 -> value beirasa
    virtual void Ini(int mode, double value) = 0;
    /// Get double property, Cso es Csatorna akarja elulirja
    virtual double Get_dprop(string mit)
    {
        return 0.0;
    }
    /// Get double property, Cso es Csatorna akarja elulirja
    double Get_Aref()
    {
        return Aref;
    }
    /// Get equation derivative w.r.t. parameter
    virtual double Get_dfdmu(string mit)
    {
        return 0.0;
    }
    virtual void Set_dprop(string mit, double mire) {

    };

    virtual void Set_friction_model(string friction_model) {};

    /// Tomegaram visszaadasa
    double Get_mp()
    {
        return mp;
    }
    /// Tomegaram beallitasa
    void Set_mp(double x)
    {
        mp = x;
        //cout<<endl<<"mp set to "<<mp<<endl;
    }
    /// Terfogataram visszaadasa
    double Get_Q()
    {
        return mp / ro;
    }
    /// Sebesseg visszaadasa
    double virtual Get_v()
    {
        /*cout << "\n" << nev << "\t mp=" << mp << "\t ro=" << ro << "\t Aref=" << Aref << "\t v="
             << (mp / ro / Aref);*/
        return mp / ro / Aref;
    }
    /// Nev visszaadasa
    string Get_nev()
    {
        return nev;
    }
    /// Tipus visszaadasa
    string Get_Tipus()
    {
        return tipus;
    }
    /// cspe_nev visszaadasa
    string Get_Cspe_Nev()
    {
        return cspe_nev;
    }
    /// cspv_nev visszaadasa
    string Get_Cspv_Nev()
    {
        return cspv_nev;
    }
    /// cspe_index visszaadasa
    int Get_Cspe_Index()
    {
        return cspe_index;
    }
    /// cspv_nev visszaadasa
    int Get_Cspv_Index()
    {
        return cspv_index;
    }
    /// csp darabszam visszaadasa
    int Get_Csp_db()
    {
        return csp_db;
    }
    /// vizkor a cso elejen
    double Get_tt_start()
    {
        return tt_start;
    }
    /// vizkor a cso vegen
    double Get_tt_end()
    {
        return tt_end;
    }
    /// vizkor a cso elejen - beallitas
    void Set_tt_start(double tmp)
    {
        tt_start = tmp;
    }
    /// vizkor a cso vegen - beallitas
    void Set_tt_end(double tmp)
    {
        tt_end = tmp;
    }
    virtual double Get_FolyTerf()
    {
        return FolyTerf;
    }
    virtual void Set_FolyTerf()
    {
    }
    virtual void build_res()
    {
    }
    virtual string GetType() = 0;
    /// Eredmenyvektor visszaadasa (csak Csatorna eseten)
    virtual vector<double> Get_res(string which);
    /// Matlab-szeru linearis interpolacio
    vector<double> interp(vector<double> x, vector<double> y, vector<double> xg);
    /// Atlag szamitasa
    double mean(vector<double> x);
    /// Vizminoseghez numerikus halo
    void set_up_grid(double a_konc, vector<double> a_vel, double a_cL);
    /// Vizminoseghez numerikus halo kiirasa
    string show_grid(double ido);
    /// Vizminoseghez koncentracio- es sebessegeloszlas
    vector<double> konc, vel;
    /// Atlagos koncentracio
    double konc_atlag;
    /// Vizminoseg szamitashoz belso ido nyilvantartasa
    double ido;
    /// Hossz, Csatorna es Cso eseten maga a hossz, az osszes tobbi elemnel 1[m]
    double cL;
    /// Az agelem vegigjarasahoz szukseges ido es idolepes
    double cT, cdt;
    /// Hozz�f�r�s nyom�ses�shez
    double Get_head_loss()
    {
        return head_loss;
    }
    /// Nyom�ses�s �r�sa
    void Set_head_loss(double a_head_loss)
    {
        head_loss = a_head_loss;
    }
    double Get_ro()
    {
        return ro;
    }
    void Set_cdt(double a_cdt)
    {
        cdt = a_cdt;
    }
    void virtual SetLogFile(string fnev)
    {
        out_file = fnev;
    }
    void Set_Aref(double a_Aref)
    {
        Aref = a_Aref;
        //cout << endl << "Agelem " << nev << ": referencia keresztmetszet beallitasa: " << a_Aref << "-> Aref=" << Aref << endl;

    }
};
#endif
