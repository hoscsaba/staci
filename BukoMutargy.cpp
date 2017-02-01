using namespace std;
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <cmath>
#include <ctime>
//#include "Agelem.h"
#include "BukoMutargy.h"

BukoMutargy::BukoMutargy(const string a_nev, const string a_cspe_nev,
                         const string a_cspv_nev, const double a_ro, const double Aref, const double a_Hf,
                         const bool a_is_opened, const double a_width,
                         const double a_overflow_height, const double a_discharge_coeff,
                         const double a_valve_coeff, const double a_mp) :
  Agelem(a_nev, Aref, a_mp, a_ro)
{
  //Kotelezo adatok minden Agelemnel:
  tipus = "BukoMutargy";
  csp_db = 2;
  cspe_nev = a_cspe_nev;
  cspv_nev = a_cspv_nev;
  //A buko mutagy adatai:
  is_opened = a_is_opened;
  Hf = a_Hf;
  width = a_width;
  overflow_height = a_overflow_height;
  discharge_coeff = a_discharge_coeff;
  valve_coeff = a_valve_coeff;
}

//--------------------------------------------------------------
BukoMutargy::~BukoMutargy()
{
}

//--------------------------------------------------------------
string BukoMutargy::Info()
{
  ostringstream strstrm;
  strstrm << Agelem::Info();
  strstrm << endl << "  kapcsolodas : " << cspe_nev << "(index:" << cspe_index
          << ") ==> " << cspv_nev << "(index:" << cspv_index << ")";
  strstrm << "       adatok : nyitva a szelep?   : ";
  if (is_opened)
    strstrm << " igen";
  else
    strstrm << " nem";
  strstrm << "\n                fen?kszint         : " << Hf;
  strstrm << "\n                sz?less?g          : " << width;
  strstrm << "\n                buk?magass?g       : " << overflow_height;
  strstrm << "\n                ?tbuk?si t?nyez?   : " << discharge_coeff;
  strstrm << "\n                szelep ?tf. t?nyez?: " << valve_coeff << endl;
  return strstrm.str();
}

//--------------------------------------------------------------
double BukoMutargy::f(vector<double> x)
{
  double ere;
  double pe = x[0] * ro * g;
  double pv = x[1] * ro * g;
  double he = x[2];
  double hv = x[3];

  if (is_opened)
  {
    ere = (pv - pe) / ro / g + (hv - he) + valve_coeff * mp * fabs(mp);
  }
  else
  {
    double hve = pe / ro / g + he - Hf;
    double hvu = pv / ro / g + hv - Hf;
    double Q = mp / ro;
    double ve = Q / width / hve;
    double vu = Q / width / hvu;
    //      cout<<endl<<endl<<"hve="<<hve<<"hvu="<<hvu<<", mp="<<mp<<", Hb="<<overflow_height;
    //int int1; cin>>int1;

    if (hve > overflow_height)
    {
      // Norm?l ir?ny, de nincs buk?s
      if (hvu > overflow_height)
      {
        ere = (hvu + vu * vu / 2 / g) - (hve + ve * ve / 2 / g);
      }
      else
      {
        // Norm?l ir?ny, van buk?s
        double H = hve - overflow_height;
        ere = mp / ro - width * H * discharge_coeff * sqrt(2 * g * H);
      }
    }
    else
    {
      // Ellent?tes ir?ny, buk?s
      if (hvu > overflow_height)
      {
        double H = hvu - overflow_height;
        ere = -mp / ro - width * H * discharge_coeff * sqrt(2 * g * H);
      }
      // beszoptuk, mindk?t oldalon a buk?si szint alatt vagyunk
      else
      {
        ere = mp / ro;
      }
    }
  }

  return ere;
}

//--------------------------------------------------------------
vector<double> BukoMutargy::df(vector<double> x)
{
  vector<double> ere, temp = x;
  double f1 = f(x), f2, dtemp;

  // df/dhe
  x.at(0) *= 1.01;
  f2 = f(x);
  ere.push_back((f2 - f1) / 0.01 / x.at(0));

  // df/dhv
  x = temp;
  x.at(1) *= 1.01;
  f2 = f(x);
  ere.push_back((f2 - f1) / 0.01 / x.at(1));

  // df/dmp
  dtemp = mp;
  mp *= 1.01;
  f2 = f(x);
  ere.push_back((f2 - f1) / 0.01 / mp);
  mp = dtemp;

  // Konstans tag:
  ere.push_back(0.0);

  // Az?rt ellen?rizni kell a buk?s lehet?s?g?t:
  double pe = x[0] * ro * g;
  double pv = x[1] * ro * g;
  double he = x[2];
  double hv = x[3];
  double hve = pe / ro / g + he - Hf;
  double hvu = pv / ro / g + hv - Hf;

  if ((hve > overflow_height) && (hvu < overflow_height))
    ere.at(1) = 0.0;
  if ((hvu > overflow_height) && (hve < overflow_height))
    ere.at(0) = 0.0;

  return ere;
}

//--------------------------------------------------------------
void BukoMutargy::Ini(int mode, double value)
{
  if (mode == 0)
    mp = 1;
  else
    mp = value;
}

//--------------------------------------------------------------
void BukoMutargy::Set_dprop(string mit, double mire) {
//    if (mit=="diameter")
//      D=mire;
//    else
//      {
  cout << endl << "HIBA! BukoMutargy::Set_dprop(mit), ismeretlen bemenet: mit="
       << mit << endl << endl;
//      }
}

double BukoMutargy::Get_dprop(string mit) {
  double out = 0.0;
  // if (mit == "Aref")
  //   out = Aref;
  // else if (mit == "lambda")
  //   out = lambda;
  // else {
  cout << endl
       << "HIBA! Cso::Get_dprop(mit), ismeretlen bemenet: mit=" << mit << endl
       << endl;
  cout << endl << "Name of BukoMutargy: " << nev << endl;
  exit(-1);
  // }
  return out;
}