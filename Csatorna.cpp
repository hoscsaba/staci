//using namespace std;
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include "Agelem.h"
#include "Csatorna.h"

//Csatorna::Csatorna(const string a_nev, const string a_cspe_nev,
//      const string a_cspv_nev, const double Aref, const double a_L,
//      const double a_ze, const double a_zv, const double a_erd,
//      const int a_int_steps, const int a_debugl, const double a_width,
//      const double a_Hmax, const double a_cl_k, const double a_cl_w) :
//  Agelem(a_nev, Aref) {
//  //Kotelezo adatok minden Agelemnel:
//  tipus="Csatorna";
//  csp_db=2;
//  cspe_nev=a_cspe_nev;
//  cspv_nev=a_cspv_nev;
//  // specifikus adatok
//  L =a_L;
//  ze =a_ze;
//  zv =a_zv;
//  lejtes =(ze-zv)/L;
//  erdesseg =a_erd;
//  geo_tipus=0;
//  // Numerikus eredmenyek tarolasara:
//  db=a_int_steps;
//  // Debug f�jl be�ll�t�sa.
//  df_name=nev+".out";
//  debug_level=a_debugl;
//  if (debug_level>4)
//      debug_level=4;
//  f_count=0;
//  df_count=0;
//  Hmax=a_Hmax;
//  B=a_width;
//  res_ready=false;
//  cl_k=a_cl_k;
//  cl_w=a_cl_w;
//}

//--------------------------------------------------------------
/// K�r km.
Csatorna::Csatorna(const string a_nev, const string a_cspe_nev,
		const string a_cspv_nev, const double a_ro, const double Aref, const double a_L,
		const double a_ze, const double a_zv, const double a_erd,
		const int a_int_steps, const int a_debugl, const double a_dia,
		const double a_cl_k, const double a_cl_w, const bool a_is_reversed, const double a_mp) :
		Agelem(a_nev, a_dia * a_dia * pi / 4, a_mp, a_ro) {
	//Kotelezo adatok minden Agelemnel:
	tipus = "Csatorna";
	csp_db = 2;
	cspe_nev = a_cspe_nev;
	cspv_nev = a_cspv_nev;
	// specifikus adatok
	L = a_L;
	ze = a_ze;
	zv = a_zv;
	lejtes = (ze - zv) / L;
	erdesseg = a_erd;
	geo_tipus = 1;
	// Numerikus eredmenyek tarolasara:
	db = a_int_steps;

	// Set debug level
	debug_level = a_debugl;
	if (debug_level > 4)
		debug_level = 4;

	f_count = 0;
	df_count = 0;
	dia = a_dia;
	Hmax = dia;
	D_fake = dia;

	res_ready = false;
	write_res = false;
	cl_k = a_cl_k;
	cl_w = a_cl_w;
	is_reversed = a_is_reversed;

	if (fabs(erdesseg) == 0.0)
		error("Konstruktor",
				"A csatorna erdessege zerus, ez nem megengedett!!! Kerem javitsa az adatot!");

	if (fabs(lejtes) < (0.01 / 100.))
		warning("Konstruktor",
				"A csatorna lejtese kisebb, mint 0.1mm/m (0.01%), ez biztosan helyes adat?");

	double (Csatorna::*pt2fun)(double, double, double) = NULL;
}

//--------------------------------------------------------------
Csatorna::~Csatorna() {
}

/// Info
/**
 * @return info string
 */
string Csatorna::Info() {
	ostringstream strstrm;
	strstrm << Agelem::Info();
	cout << setprecision(3);
	strstrm << endl << "  kapcsolodas : " << cspe_nev << "(index:" << cspe_index
			<< ") ==> " << cspv_nev << "(index:" << cspv_index << ")";
	strstrm << endl << "       adatok : hossz [m]         : " << L;
	strstrm << endl << "                z eleje [m]       : " << ze;
	strstrm << endl << "                z vege [m]        : " << zv;
	strstrm << endl << "                lejtes [%]        : " << (lejtes * 100);
	strstrm << endl << "                Hmax [m]          : " << Hmax;
	strstrm << endl << "                erdesseg [mm]             : " << erdesseg;
	surlodas();
	strstrm << endl << "                lambda [-]                : " << lambda;
	strstrm << endl << "                Manning surl.teny.(n) [-] : " << surl;
	strstrm << endl << "                geometria tipus           : ";
	switch (geo_tipus) {
	case 0:
		strstrm << B << "m szelessegu, " << Hmax
		<< "m magassagu teglalap keresztmetszet";
		break;
	case 1:
		strstrm << dia << "m atmeroju kor keresztmetszet";
		break;
	default:
		error("Info()", " Ismeretlen geometria tipus!");
	}
	strstrm << endl << "                numerikus megoldo beallitasai (db)     : "
			<< db;
	strstrm << endl << "                klor lebomlasi allando (cl_k)          : "
			<< cl_k;
	strstrm << endl << "                klor lebomlasi allando a falnal (cl_w) : "
			<< cl_w;
	strstrm << endl << "                debug szint                            : "
			<< debug_level << endl;
	strstrm << endl << "                logfile                                : "
			<< out_file << endl;


	if (debug_level > 0) {
		time_t ido = time(0);
		ostringstream strstrm1;
		strstrm1 << "Staci\nBME Hidrodinamikai Rendszerek Tanszek\n" << ctime(&ido)
                		 << "\n";
		ofstream outfile(out_file.c_str(), ios::trunc);
		outfile << strstrm1.str();
		outfile.close();
		logfile_write(strstrm.str(), 1);
	}
	return strstrm.str();
}

/// �gegyenlet �s deriv�ltak sz�m�t�sa
/**
 * @param x[0] = pe/ro/g - nyom�s [v.o.m.] az �gelem elej�n
 * @param x[1] = pv/ro/g - nyom�s [v.o.m.] az �gelem v�g�n
 * @param x[2] = he/ro/g - �gelem elej�n a csom�pont nullszintje [m]
 * @param x[3] = hv/ro/g - �gelem v�g�n a csom�pont nullszintje [m]
 * @return f �rt�ke
 *
 * �gelem elej�n az abszol�t nyom�sszint: he(aknafen�k) + pe/ro/g = x[0]+x[1]
 * �gelem elej�n az abszol�t nyom�sszint: hv(aknafen�k) + pv/ro/g = x[0]+x[1]
 *
 * Az elj�r�s kisz�m�tja a deriv�ltakat is, �gy @see df csak visszaadja az �rt�ket.
 */

double Csatorna::f(vector<double> x) {
	pe = x[0] * ro * g;
	pv = x[1] * ro * g;
	he = x[2];
	hv = x[3];
	double ye = pe / ro / g + he - ze;
	double yv = pv / ro / g + hv - zv;
	ostringstream strstrm;

	headloss = fabs(ye - yv);

	strstrm << endl << endl << " evaluating f (#" << f_count << ")" << endl
			<< "-----------------------------------------";
	strstrm << endl << "  Data:" << scientific << showpos << " pe=" << pe / ro / g << "m, he="
			<< he << "m, pv=" << pv / ro / g << "m, hv=" << hv << "m, mp=" << mp;
	strstrm << endl << "        ye=" << ye << "m, ye+ze=" << (ze + ye) << "m, yv=" << yv
			<< ", yv+zv=" << (zv + yv) << "m";
	logfile_write(strstrm.str(), 1);

	//double (Csatorna::*pt2fun)(double, double, double) = NULL;

	jac.clear();
	num_eval_jac.clear();
	for (unsigned int ii = 0; ii < 3; ii++) {
		jac.push_back(0.0);
		num_eval_jac.push_back(true);
	}

	which_case(ye, yv);

	//***************
	// SZAMITASOK
	//***************

	logfile_write("\n\n\t  * evaluating f...", 3);
	double ff = (*this.*pt2fun)(ye, yv, mp);

	// Jacobi elemeinek szamitasa
	double dx, df;

	// dfdye szamitasa:
	logfile_write("\n\n\t  * evaluating df/dye...", 3);
	if (num_eval_jac[0]) {
		dx = -0.001 * ye;
		df = (*this.*pt2fun)(ye + dx, yv, mp);
		jac[0] = (df - ff) / dx;
	}

	// dfdyv szamitasa:
	logfile_write("\n\n\t  * evaluating df/dyv...", 3);
	if (num_eval_jac[1]) {
		dx = -0.001 * yv;
		df = (*this.*pt2fun)(ye, yv + dx, mp);
		jac[1] = (df - ff) / dx;
	}

	// dfdmp szamitasa:
	logfile_write("\n\n\t  * evaluating df/dmp...", 3);
	if (num_eval_jac[2]) {
		double mdot_old = mp;
		dx = -0.001 * mp;
		mp = mp + dx;
		which_case(ye, yv);
		df = (*this.*pt2fun)(ye, yv, mp);
		jac[2] = (df - ff) / dx;
		mp = mdot_old;
	}

	// konstans tag, csak linearizalas eseten van jelentosege
	jac.push_back(0.);

	strstrm.str("");
	strstrm << "\n\n\n\t      f  = " << scientific << ff;
	strstrm << "\n\t df/dye  = " << jac[0];
	strstrm << "\n\t df/dyv  = " << jac[1];
	strstrm << "\n\t df/dmp  = " << jac[2];
	logfile_write(strstrm.str(), 1);

	f_count++;

	return ff;
}


void Csatorna::which_case(const double ye, const double yv) {
	ostringstream strstrm;


	strstrm.str("");
	if ((ye > dia) && (yv > dia)) {
		// telt eset
		eset = "teltszelveny";
		strstrm << "\n\n\t " << eset << ": ye=" << ye << "m > D and yv=" << yv << "m >D (D=" << dia << "m)";
		pt2fun = &Csatorna::f_telt;
	} else {
		if (ye < 0) {
			if (yv + zv < ze) {
				// 0.a. eset
				eset = "0.a.";
				strstrm << "\n\n\t case " << eset << ":  ye=" << ye << "m <0 and yv=" << yv << "m < ze=" << ze << "m";
				pt2fun = &Csatorna::f_0;
			} else {
				if (mp > 0) {
					// 0.b.i. eset
					eset = "0.b.i.";
					strstrm << "\n\n\t case " << eset << ":  ye=" << ye << "m <0, yv=" << yv << "m >0, mp=" << mp
							<< " kg/s >0";

					pt2fun = &Csatorna::f_0;
				} else {
					// 0.b.ii. eset
					eset = "0.b.ii.";
					strstrm << "\n\n\t case " << eset << ":  ye=" << ye << "m <0, yv=" << yv << "m >0, mp=" << mp
							<< " kg/s <0";

					pt2fun = &Csatorna::f_1;
					//eval_jac[0]= false;
				}
			}
		} else {
			if (mp < 0) {
				// 1. eset
				eset = "1.";
				strstrm << "\n\n\t case " << eset << ":  ye=" << ye << "m >0, mp=" << mp << " kg/s <0";

				pt2fun = &Csatorna::f_1;
			} else {
				vector<double> yn = normal_szint(mp / ro);
				double yc = kritikus_szint(mp / ro);

				switch (yn.size()) {
				case 1:
					if (yn[0] > yc) {
						// 2.a.i. eset
						eset = "2.a.i.";
						strstrm << "\n\n\t case " << eset << ":  ye=" << ye << " m >0, mp=" << mp
								<< " kg/s >0, yn>yc";
						strstrm << "\n\t\t(1 db normalszint: yn=" << yn[0] << " m > yc=" << yc << " m)";

						pt2fun = &Csatorna::f_2c;
					} else {
						// 2.a.ii. eset
						eset = "2.a.ii.";
						strstrm << "\n\n\t case " << eset << ":  ye=" << ye << " m >0, mp=" << mp
								<< " kg/s >0, yn<yc";
						strstrm << "\n\t\t(1 db normalszint: yn=" << yn[0] << " m < yc=" << yc << " m)";

						pt2fun = &Csatorna::f_2aii;
					}
					break;

				case 2:

					double yn1, yn2;
					// yn1<yn2, rendez�s sz�ks�ges lehet.
					if (yn[0] > yn[1]) {
						yn1 = yn[1];
						yn2 = yn[0];
					} else {
						yn1 = yn[0];
						yn2 = yn[1];
					}

					if (yc < yn1) {
						eset = "2.b.i.";
						strstrm << "\n\n\t case " << eset << ":  ye=" << ye << "m >0, mp=" << mp
								<< "kg/s >0, yc<yn1";
						strstrm << "\n\t\t(2 db normalszint: yc=" << yc << "m < yn1=" << yn1 << "m, yn2=" << yn2
								<< "m)";

						pt2fun = &Csatorna::f_2c;
					} else {
						if (yc < yn2) {
							eset = "2.b.ii.";
							strstrm << "\n\n\t case " << eset << ":  ye=" << ye << "m >0, mp=" << mp
									<< "kg/s >0, yn1<yc<yn2";
							strstrm << endl << "\n\t\t(2 db normalszint: yn1=" << yn1 << "m < yc=" << yc
									<< "m < yn2=" << yn2 << "m";

							pt2fun = &Csatorna::f_2aii;
						} else {
							eset = "2.b.iii.";
							strstrm << "\n\n\t case " << eset << ":  ye=" << ye << "m >0, mp=" << mp
									<< "kg/s >0, yn2<yc";
							strstrm << "\n\t\t(2 db normalszint: yn1=" << yn1 << "m < yn2=" << yn2 << "m < yc="
									<< yc << "m";

							//pt2fun = &Csatorna::f_2c;
							pt2fun = &Csatorna::f_2aii;
						}
					}
					break;

				case 0:
					eset = "2.c.";
					strstrm << "\n\n\t case " << eset << ":  ye=" << ye << "m >0, mp=" << mp << "kg/s >0, yn nincs";
					strstrm << "\n\t\t(0 db normalszint, yc=" << yc << "m)";

					pt2fun = &Csatorna::f_2c;

					break;

				default: {
					strstrm.str("Bajjjj van: yn.size()=");
					strstrm << yn.size();
					error("f()", strstrm.str());
				}
				}
			}
		}
	}
	logfile_write(strstrm.str(), 1);
}

/// �gegyenlet deriv�ltjai (sz�m�t�s: @see f)
/**
 * @return df �rt�ke
 */
vector<double> Csatorna::df(vector<double> x) {
	return jac;
}

/// �gelem inicializ�l�sa
/**
 * @param mode - (0) automatikus, (1) mp=value
 * @param value csak mode=1 eset�n
 */

void Csatorna::Ini(int mode, double value) {
	if (mode == 0)
		mp = 1.e-3;
	else
		mp = value;
}

/// Keresztmetszeti jellemz�k sz�m�t�sa
/**
 * @param y v�zszint
 * @param &A nedves�tett ter�let
 * @param &B nedves�tett ker�let
 * @param &Rh hidraulikai sug�r (A/B)
 */

void Csatorna::keresztmetszet(const double yy, double &A, double &B, double &Rh) {
	double K = 0.0;
	double y = yy;
	if (y < 0) {
		//ostringstream msg;
		//msg.str("");
		//msg<<"Negat�v szint, y="<<y;
		//warning("keresztmetszet()", msg.str());
		y = dia / 1000.;
	} else {
		switch (geo_tipus) {
		case 0:   // teglalap
		{
			A = B * y;
			K = B + 2 * y;
			break;
		}
		case 1:   // kor
		{
			if (y < dia) {
				double teta, R = dia / 2;
				teta = acos(1 - y / R);
				K = 2 * R * teta;
				A = dia * dia / 4. * (teta - sin(teta) * cos(teta));
				B = dia * sin(teta);
			} else {
				double h = dia / 100;
				A = dia * dia * pi / 4.;
				B = h;
				K = dia * pi;
			}
			break;
		}
		default:
			ostringstream strstrm;
			strstrm.str("");
			strstrm << "Ismeretlen geometria tipus:" << geo_tipus;
			error("keresztmetszet()", strstrm.str());
		}

		Rh = A / K;
	}
}

/// Kritikus szint sz�m�t�sa k�r keresztmetszet eset�n
/**
 * yc(0) -> -inf, yc(D)->1, �gy a kett� k�z�tt pontosan egy gy�khely van.
 * Megold�s: intervallumfelez�ssel
 *
 * @param Q t�rfogat�ram
 * @return kritikus szint
 */

double Csatorna::kritikus_szint(const double Q) {
	double fun = 1.0, yj = dia, yb = dia, A, B, Rh, Q2 = Q * Q;
	ostringstream strstrm;

	strstrm << scientific << showpos;
	strstrm << "\n\n\tKritikus szint szamitasa: Q=" << Q << " m3/s";

	// Cs�kken� yb mellett keres�nk egy negat�v �rt�ket
	unsigned int iter = 0;
	strstrm << "\n\t\tEl�sz�m�t�s: csokkeno yb menten negativ f ertek keresese:";
	while (fun > 0) {
		yb /= 2.;
		keresztmetszet(yb, A, B, Rh);
		fun = 1. - Q2 * B / (A * A * A) / g;
		iter++;

		strstrm << "\n\t\t\t iter=" << iter << ", y=" << yb << " m, f=" << fun;
		if (iter > iter_max) {
			logfile_write(strstrm.str(), 0);
			error("kritikus_szint()", "HIBA (1) - Nem tal�lok negat�v �rt�ket!");
		}
	}
	logfile_write(strstrm.str(), 4);

	// H�rm�dszer
	fun = 1e10;
	iter = 0;
	double yk = 0.0;

	strstrm.str("");
	strstrm << "\n\n\t\tSz�m�t�s:";

	keresztmetszet(yb, A, B, Rh);
	double fb = 1. - Q2 * B / (A * A * A) / g;
	keresztmetszet(yj, A, B, Rh);
	double fj = 1. - Q2 * B / (A * A * A) / g;

	while ((fabs(fun) > 1.e-10) && (fabs(yj - yb) > 1.e-10)) {
		strstrm << "\n\t\t\titer=" << iter << ": yb=" << yb << "(fb=" << fb << "),  yj=" << yj
				<< "(fj=" << fj << ")";

		yk = (yj + yb) / 2;
		//yk=yb-fb*(yj-yb)/(fj-fb);

		keresztmetszet(yk, A, B, Rh);
		fun = 1. - Q2 * B / (A * A * A) / g;

		if (fun > 0) {
			yj = yk;
			fj = fun;
		} else {
			yb = yk;
			fb = fun;
		}

		iter++;

		if (iter > iter_max) {
			logfile_write(strstrm.str(), 0);
			error("kritkus_szint()", "HIBA (2) - T�l sok iter�ci�!");
		}
	}
	logfile_write(strstrm.str(), 4);

	strstrm.str("");
	strstrm << "\n\n\t\tKritikus szint: yc = " << yk << " m";
	logfile_write(strstrm.str(), 2);

	return yk;
}

/// Norm�lszint sz�m�t�sa k�r keresztmetszet eset�n
/**
 * yn(0) -> -inf
 * ha Q<Q_D -> 1db gy�khely
 * ha Q_D<Q<Q_cs -> 2db gy�khely
 * ha Q>Q_cs -> 0db gy�khely
 * Megold�s: h�rm�dszer
 *
 * @param Q t�rfogat�ram
 * @return norm�lszint(ek)
 */

vector<double> Csatorna::normal_szint(const double Q) {
	vector<double> yn;

	double y_cs = 0.938181 * dia;
	double A, B, Rh, Q2 = Q * Q;

	keresztmetszet(y_cs, A, B, Rh);
	double C = pow(Rh, 1.0 / 6.0) / surlodas();
	double Q_cs = sqrt(lejtes * A * A * C * C * Rh);

	keresztmetszet(dia, A, B, Rh);
	C = pow(Rh, 1.0 / 6.0) / surlodas();
	double Q_D = sqrt(lejtes * A * A * C * C * Rh);

	ostringstream strstrm;
	strstrm << scientific << showpos;
	strstrm << "\n\n\tNormalszint szamitasa: Q=" << Q << " m3/s, Q_cs=" << Q_cs
			<< "m3/s, Q_D=" << Q_D << "m3/s";
	logfile_write(strstrm.str(), 4);

	if (Q < Q_cs) {
		strstrm.str("\n\n\t\t1. normalszint szamitasa:");
		unsigned int iter = 0;
		double fun = 1.e10;
		double y = dia / 2;
		while (fun > 0) {
			keresztmetszet(y, A, B, Rh);
			C = pow(Rh, 1.0 / 6.0) / surlodas();
			fun = lejtes - Q2 / A / A / C / C / Rh;
			strstrm << "\n\t\t -> " << " y=" << y << ", f=" << fun;

			y /= 2.;
			iter++;
			if (iter > iter_max) {
				logfile_write(strstrm.str(), 0);
				error("normal_szint()",
						"HIBA (1) - Nem talalok negativ erteket!");
			}
		}
		logfile_write(strstrm.str(), 4);

		// 1. norm�lszint sz�m�t�sa: 0<y<y_cs
		// H�rm�dszer
		fun = 1e10;
		iter = 0;
		double yk = 0.0;

		strstrm.str("\n\n\t\tSzamitas:");

		double yb = y, yj = y_cs;
		keresztmetszet(yb, A, B, Rh);
		C = pow(Rh, 1.0 / 6.0) / surlodas();
		double fb = lejtes - Q2 / A / A / C / C / Rh;
		keresztmetszet(yj, A, B, Rh);
		C = pow(Rh, 1.0 / 6.0) / surlodas();
		double fj = lejtes - Q2 / A / A / C / C / Rh;

		while ((fabs(fun) > 1e-10) && fabs(yb - yj) > 1e-10) {
			//yk=yb-fb*(yj-yb)/(fj-fb);
			yk = (yb + yj) / 2.;

			keresztmetszet(yk, A, B, Rh);
			C = pow(Rh, 1.0 / 6.0) / surlodas();
			fun = lejtes - Q2 / A / A / C / C / Rh;

			strstrm << "\n\t\t\titer=" << iter << ": yb=" << yb << "(fb=" << fb << ")  yk="
					<< yk << "(fk=" << fun << "), yj=" << yj << "(fj=" << fj << ")";
			//cout<<"\n\t\t\titer="<<iter<<": yb="<<yb<<"(fb="<<fb<<")  yk=" <<yk <<"(fk="<<fun <<"), yj="<<yj<<"(fj="<<fj<<")";;
			if (fun > 0) {
				yj = yk;
				fj = fun;
			} else {
				yb = yk;
				fb = fun;
			}
			iter++;

			if (iter > iter_max) {
				logfile_write(strstrm.str(), 0);
				error("normal_szint()",
						"HIBA (2) - T�l sok iter�ci�, 1. norm�lszint sz�m�t�sa!");
			}
		}
		yn.push_back(yk);
		logfile_write(strstrm.str(), 4);

		if (Q > Q_D) {
			// 2. norm�lszint sz�m�t�sa: 0<y<y_cs
			strstrm.str("");
			strstrm << "\n\n\t2. norm�lszint sz�m�t�sa:";

			double yb = y_cs, yj = dia;
			keresztmetszet(yb, A, B, Rh);
			C = pow(Rh, 1.0 / 6.0) / surlodas();
			double fb = lejtes - Q2 / A / A / C / C / Rh;
			keresztmetszet(yj, A, B, Rh);
			C = pow(Rh, 1.0 / 6.0) / surlodas();
			double fj = lejtes - Q2 / A / A / C / C / Rh;

			double fun = 1e10, yk;
			unsigned int iter = 0;
			while ((fabs(fun) > 1e-10) && (fabs(yb - yj) > 1e-10)) {
				//yk=yb-fb*(yj-yb)/(fj-fb);
				yk = (yb + yj) / 2.;
				keresztmetszet(yk, A, B, Rh);
				C = pow(Rh, 1.0 / 6.0) / surlodas();
				fun = lejtes - Q2 / A / A / C / C / Rh;

				strstrm << "\n\t\t\titer=" << iter << ": yb=" << yb << "(fb=" << fb
						<< ")  yk=" << yk << "(fk=" << fun << "), yj=" << yj << "(fj="
						<< fj << ")  |yb-yj|=" << fabs(yb - yj);
				if (fun > 0) {
					yj = yk;
					fj = fun;
				} else {
					yb = yk;
					fb = fun;
				}
				iter++;

				if (iter > iter_max) {
					logfile_write(strstrm.str(), 0);
					error("normal_szint()",
							"HIBA (3) - T�l sok iter�ci� h�rm�dszer k�zben, 2. norm�lszint sz�m�t�sa!");
				}
				//cout<<endl<<iter<<": yk="<<yk<<"  |yb-yj|="<<fabs(yb-yj)<<"  fun="<<fabs(fun);

			}
			yn.push_back(yk);
			logfile_write(strstrm.str(), 4);
		}
	}

	strstrm.str("");
	strstrm << "\n\n\t\tNormalszint(ek):";
	switch (yn.size()) {
	case 0:
		strstrm << " Nincs";
		logfile_write(strstrm.str(), 2);
		break;
	case 1:
		strstrm << " yn1 = " << yn[0] << " m";
		logfile_write(strstrm.str(), 2);
		break;
	case 2:
		strstrm << " yn1 = " << yn[0] << " m";
		strstrm << " yn2 = " << yn[1] << " m";
		logfile_write(strstrm.str(), 2);
		break;
	default:
		strstrm.str("");
		strstrm << " yn.size()=" << yn.size() << " ???";
		error("normal_szint", strstrm.str());
	}

	return yn;
}

/// Folyad�kfelsz�nt le�r� dy/dx=f(x,y) KDE jobboldala
/**
 * @param x x koordin�ta
 * @param y y koordin�ta
 * @return dy/dx �rt�ke
 */

double Csatorna::nyf_ode(const double x, const double y, const double mp) {
	double Q, A, B, C, Rh, ere = 0.0;
	Q = mp / ro;

	if (y < dia) {
		keresztmetszet(y, A, B, Rh);
		C = pow(Rh, 1. / 6.) / surlodas();
		ere = (lejtes - fabs(Q) * Q / A / A / C / C / Rh) / (1 - Q * Q / A / A / A / g * B);
	} else {
		A = dia * dia * pi / 4.;
		Rh = dia / 4.0;
		C = pow(Rh, 1. / 6.) / surlodas();
		ere = lejtes - fabs(Q) * Q / A / A / C / C / Rh;
	}

	return ere;
}

/// ODE solver
/**
 * @param y0  depth at x=0
 * @param dx0 initial dx step
 * @param x0  initial x value (defines direction of integration)
 * @return y result of integration
 */

double Csatorna::ode_megoldo(double y0, double dx0, double x0, double mp) {
	ostringstream strstrm;
	double A, B, Rh, y = y0, x = x0, dx = dx0, sumdx = 0, yy = y0;
	int i = 0;

	// clear vectors and save initial point
	if (write_res) {
		xres.clear();
		yf.clear();
		yres.clear();
		vres.clear();
		xres.push_back(x);
		if (dx0 > 0)
			yf.push_back(ze - lejtes * x);
		else
			yf.push_back(zv + lejtes * (L - x));
		yres.push_back(yy);

		keresztmetszet(yy, A, B, Rh);
		vres.push_back(mp / ro / A);
	}

	strstrm.str("");
	strstrm << endl << endl << "\t\t\tSolving open-surface DE in channel " << nev << " ...\n";
	strstrm << scientific << setprecision(4) << showpos;
	strstrm << " x=" << x << " y=" << y << " (mp=" << mp << " kg/s)";
	logfile_write(strstrm.str(), 2);

	bool last_step = false;

	while (!last_step) {
		strstrm.str("");
		double y1, y2, f_at_y, f_at_y2;
		double hiba_max = 1e-6, dx_min = 1e-10, hiba = 2.0 * hiba_max, dxuj;

		while (hiba > hiba_max) {
			// Is this the last step?
			if (x + dx > L) {
				dx = L - x;
				last_step = true;
			}
			if (x + dx < 0) {
				dx = -x;
				last_step = true;
			}

			// one full step...        .
			f_at_y = nyf_ode(x, y, mp);
			y1 = y + dx * f_at_y;

			// two falf steps...
			y2 = y + dx / 2 * f_at_y;
			f_at_y2 = nyf_ode(x + dx / 2, y2, mp);
			y2 = y2 + dx / 2 * f_at_y2;

			// compare results
			hiba = fabs(y1 - y2);

			// Info
			strstrm.str("");
			strstrm << endl << "\t\t  x=" << scientific << setprecision(4) << showpos << x
					<< " dx=" << dx << " y=" << y << ", f(x,y)=" << f_at_y;
			strstrm << ", y1=" << y1 << ", f(x+dx/2)=" << f_at_y2 << ", y2=" << y2
					<< ", hiba=" << hiba << ", hiba_max=" << hiba_max;
			logfile_write(strstrm.str(), 3);
			strstrm.str("");

			// Lepes beallitasa
			if ((hiba > hiba_max) || (y1 < 0.)) {
				dxuj = dx / 2.;
				logfile_write(" dx -> dx/2", 3);
				if (fabs(dxuj) < dx_min) {
					strstrm << endl << "!!! Feladom, a " << nev
							<< " csatornaban az eloirt dx_min=" << dx_min
							<< " lepeskozzel sem tudom elerni a megadott "
							<< hiba_max << " hibahatart!!!\n A szamitas folytatodik, de pontatlan lehet!!!" << endl;

					error("ode_megoldo()", strstrm.str());
				}
			} else {
				logfile_write(" OK", 3);
				if (hiba < hiba_max / 10) {
					dxuj = dx * 2.;
					logfile_write(" dx -> dx*2", 3);
				} else
					dxuj = dx;

				// Lepes lezarasa
				y = y2;
				x += dx;
				sumdx += fabs(dx);
				i++;
			}
			dx = dxuj;
		}

		// Collect results
		if (write_res) {
			xres.push_back(x);
			if (dx0 > 0)
				yf.push_back(ze - lejtes * x);
			else
				yf.push_back(zv + lejtes * (L - x));
			yres.push_back(y);
			keresztmetszet(y, A, B, Rh);
			vres.push_back(mp / ro / A);
		}
		logfile_write(strstrm.str(), 4);
	}
	strstrm.str("");
	strstrm << "\n\t\t\t => solution: x=" << x << " y=" << y;
	logfile_write(strstrm.str(), 2);

	return y;
}

/// Logfile-ba �r�s
/**
 * @param msg �zenet
 * @param msg_debug_level csak akkor �r�dik be msg a nev.out f�jlba, ha msg_debug_level>=debug_level
 */

/*void Csatorna::logfile_write(string msg, int msg_debug_level) {
    if (debug_level>=msg_debug_level) {
        ofstream outfile(df_name.c_str(), ios::app);
        outfile<<msg;
        outfile.close();
    }
}*/

/// Eredm�nyvektorok �p�t�se
/**
 * xres: x koordinata
 * yf  : fenek y koordinata
 * yres: nyomasvonal y koordinata
 * Hres: vizszint y koordinata
 * vres: sebesseg
 * !!! Folyasiranyban kell rajzolni!
 */

void Csatorna::build_res() {
	write_res = true;
	res_ready = true;
	double ye = pe / ro / g + he - ze;
	double yv = pv / ro / g + hv - zv;

	xres.clear();
	yf.clear();
	yres.clear();
	Hres.clear();
	vres.clear();

	// Empty channel hadnled separately
	bool megvolt = false;
	if (!strcmp(eset.c_str(), "teltszelveny")) {
		megvolt = true;
		xres.push_back(0);
		xres.push_back(L);
		yf.push_back(ze);
		yf.push_back(zv);
		yres.push_back(ye);
		yres.push_back(yv);
	}

	if (!strcmp(eset.c_str(), "0.a.") || !strcmp(eset.c_str(), "0.b.i.")) {
		// f_0
		megvolt = true;
		xres.push_back(0);
		xres.push_back(L);
		yf.push_back(ze);
		yf.push_back(zv);
		yres.push_back(ye);
		yres.push_back(yv);
	}

	if (!strcmp(eset.c_str(), "0.b.ii.") ) {
		// f_1
		megvolt = true;
		double ff = f_1(ye, yv, mp);
		// Az eredmenyvektorokat folyasiranyban rajzoljuk fel!
		is_reversed = true;
	}

	if (!strcmp(eset.c_str(), "1.")) {
		// f_1
		megvolt = true;
		double ff = f_1(ye, yv, mp);
		// Az eredmenyvektorokat folyasiranyban rajzoljuk fel!
		is_reversed = false;
	}

	if (!strcmp(eset.c_str(), "2.a.ii.") || !strcmp(eset.c_str(), "2.b.ii.")) {
		// f_2aii
		megvolt = true;
		double ff = f_2aii(ye, yv, mp);
	}

	if (!strcmp(eset.c_str(), "2.a.i.") || !strcmp(eset.c_str(), "2.b.i.")
			|| !strcmp(eset.c_str(), "2.b.iiii.") || !strcmp(eset.c_str(),
					"2.c.")) {
		// f_2c
		megvolt = true;
		double ff = f_2c(ye, yv, mp);
	}

	if (!megvolt) {
		string msg = "Ismeretlen eset: " + eset;
		error("build_res()", msg);
	}

	// Hmax korrekcio
	for (unsigned int i = 0; i < yres.size(); i++) {

		if (yres.at(i) > dia)
			Hres.push_back(dia);
		else if (yres.at(i) < 0)
			Hres.push_back(0.0);
		else
			Hres.push_back(yres.at(i));
	}

	// Az eredm�nyek konzisztenci�ja miatt mindig utolag szamitjuk yres-bol.
	double A, B, Rh;
	for (unsigned int i = 0; i < yres.size(); i++) {
		keresztmetszet(Hres.at(i), A, B, Rh);
		vres.push_back(mp / ro / A);
		yf.at(i) -= zv;
		//      if (i==0) {
		//          //veszt.at(i)=0.0;
		//      } else {
		//          C2=pow(Rh, 1./3.)/surlodas()/surlodas();
		//          //veszt.at(i)=veszt.at(i-1)+vere.at(i)*vere.at(i)/C2/Rh *(xere.at(i) -xere.at(i-1));
		//      }
	}

	// Ha meg volt forditva, vissza kell forgatni a cuccost.
	if (is_reversed) {
		vector<double> tmp_xres(xres.size()), tmp_yf(xres.size()),
				tmp_yres(xres.size()), tmp_Hres(xres.size()),
				tmp_vres(xres.size());
		int max = xres.size();
		for (unsigned int i = 0; i < xres.size(); i++) {
			tmp_xres.at(i) = L - xres.at(max - i - 1);
			tmp_yf.at(i) = yf.at(max - i - 1);
			tmp_yres.at(i) = yres.at(max - i - 1);
			tmp_Hres.at(i) = Hres.at(max - i - 1);
			tmp_vres.at(i) = vres.at(max - i - 1);
		}
		xres.clear();
		yf.clear();
		yres.clear();
		Hres.clear();
		vres.clear();

		xres = tmp_xres;
		yf = tmp_yf;
		yres = tmp_yres;
		Hres = tmp_Hres;
		vres = tmp_vres;
	}

	// Info
	ostringstream strstrm;
	strstrm << "\n\n\tBuilding solution vectors...";
	strstrm << "\n\n\t\t ye+ze = " << (pe / ro / g + he) << "m, ye=" << ye << "m";
	strstrm << "\n\t\t yv+zv = " << (pv / ro / g + hv) << "m, yv=" << yv << "m";
	strstrm << "\n\t\t mp    = " << mp << "kg/s";
	strstrm << "\n\t\t case  : " << eset;
	strstrm << "\n\t\t length of solution vectors : " << yres.size() << "\n";

	strstrm << "\n\n Solution vectors:\n\n";
	strstrm
	<< "\n\t  i    x      \t   yf     \t  yf+D   \t   y   \t     p (mwc) \t  yf+p    \t   v    \t     A";
	strstrm << scientific << setprecision(3);
	double sumh;

	double AA, BB, RRh;
	for (unsigned int i = 0; i < yres.size(); i++) {
		keresztmetszet(Hres.at(i), AA, BB, RRh);

		// !!!!!!!!!!!!!!!!!!
		// At kell allitani Aref-et, mert kulonben a sebesseg a teljes keresztmetszettel lesz szamolva es az hulyeseg.
		if (i == 0)
			Set_Aref(AA);
		// !!!!!!!!!!!!!!!!!!

		sumh = yf.at(i) + yres.at(i);
		strstrm << "\n\t  " << i << "\t" << xres.at(i) << "\t" << yf.at(i) << "\t";
		strstrm << (yf.at(i) + dia) << "\t" << Hres.at(i) << "\t";
		strstrm << yres.at(i) << "\t" << sumh << "\t" << vres.at(i) << "\t" << AA;
	}
	logfile_write(strstrm.str(), 2);
}

//--------------------------------------------------------------
vector<double> Csatorna::Get_res(string mit) {
	vector<double> empty;

	if (!res_ready)
		build_res();

	if (mit == "xf")
		return xres;
	else if (mit == "x")
		return xres;
	else if (mit == "yf")
		return yf;
	else if (mit == "y")
		return Hres;
	else if (mit == "p")
		return yres;
	else if (mit == "v")
		return vres;
	else {
		ostringstream strstrm;
		strstrm.str("");
		strstrm << "Nincs ilyen valtozo: " << mit;
		error("Get_res()", strstrm.str());
		return empty;
	}
}

//--------------------------------------------------------------
double Csatorna::surlodas() {
	// surl: Manning-allando
	if (erdesseg <= 0) {
		lambda = -erdesseg;
		surl = -erdesseg;
	} else {
		if (f_count >= 0) {
			Hmax = dia;
			double ize = 2.0 * log(14.8 * (Hmax / 2) / (erdesseg / 1000));
			lambda = 1 / ize / ize;
		} else
			lambda = 0.02;
		surl = pow(Hmax / 2., 1. / 6.) * sqrt(lambda / 8.0 / g);
	}
	//cout<<endl<<nev<<": surl="<<surl<<" (lambda="<<lambda<<")";

	return surl;
}

//--------------------------------------------------------------
double Csatorna::Get_dprop(string mit) {
	double out = 0.0;
	if (mit == "Aref")
		out = Aref;
	else if (mit == "lambda")
		out = lambda;
	else if (mit == "L")
		out = L;
	else if (mit == "Rh") {
		double A, B, Rh;
		keresztmetszet(Hmax, A, B, Rh);
		out = Rh;
	} else if (mit == "cl_k")
		out = cl_k;
	else if (mit == "cl_w")
		out = cl_w;
	else if (mit == "headloss")
		out = headloss;
	else if (mit == "headloss_per_unit_length")
		out = headloss / L;
	else if (mit == "ze")
		out = ze;
	else if (mit == "zv")
		out = zv;
	else if (mit == "lejtes")
		out = lejtes;
	else if ((mit == "concentration") || (mit == "konc_atlag"))
		out = konc_atlag;
	else {
		ostringstream msg;
		msg.str("");
		msg << "Ismeretlen bemenet: mit=" << mit;
		error("Get_dprop()", msg.str());
	}
	return out;
}

//--------------------------------------------------------------
void Csatorna::Set_FolyTerf() {
	double A, B, Rh, szum = 0;
	if (eset == "-2") {
		szum = 0;
	} else {
		for (unsigned int i = 1; i < Hres.size(); i++) {
			keresztmetszet(Hres.at(i), A, B, Rh);
			szum += A * fabs(xres.at(i) - xres.at(i - 1));
			//cout<<endl<<"\t\t i="<<i<<", x(i)="<<xere.at(i)<<", x(i-1)="<<xere.at(i
			//-1) <<", dV="<<A*fabs(xere.at(i)-xere.at(i-1));
		}
	}
	//cout<<endl<<nev<<"\t V= "<<szum<<" m^3";
	FolyTerf = szum;
}

/// �gegyenlet (0.a.) ye<0, yv<0 �s (0.b.i.) ye<0, yv>0, Q>0 esetben
/**
 * D/100 �tm�r�j� sz�v�sz�lon teltszelv�ny� �raml�s ye-yv nyom�sk�l�nbs�g hat�s�ra
 * @param ye eleje v�zszint
 * @param yv v�ge v�zszint
 * @param mp t�meg�ram
 * @return �gegyenlet �rt�ke
 */
double Csatorna::f_telt(double ye, double yv, double mp) {
	// Matlab kod:
	//  Df=D/100;
	//  Af=Df^2*pi/4;
	//  ff=yea-yva-0.02*L/Df/2/9.81*Q*abs(Q)/Af^2;
	//  eval_jac=[1 1 1];

	ostringstream strstrm;

	double v = mp / ro / (D_fake * D_fake * pi / 4);
	double dh = 0.02 * L / D_fake / 2 / 9.81 * v * fabs(v);

	double f = ze + ye - (zv + yv) - dh;

	// If the pipe was empty in the previous step, we reduce the diameter.
	if (D_fake > dia) {
		strstrm.str("\n\t\t  f_telt: fake pipe, D_fake=");
		strstrm << D_fake << "m > D=" << dia << "m";
		if (fabs(f) < 1.e-3) {
			D_fake /= 2.;
			v = mp / ro / (D_fake * D_fake * pi / 4);
			dh = 0.02 * L / D_fake / 2 / 9.81 * v * fabs(v);
			f = ze + ye - (zv + yv) - dh;
			strstrm << "\n\t\t\t  D=" << D_fake << ", f=" << f;
		}
	} else
		D_fake = dia;
	// end of d_fake update

	strstrm << "\n\t\t\t f = (ze+ye)-(zv+yv)-dp = " << f;
	logfile_write(strstrm.str(), 3);

	jac[0] = +1.0;
	num_eval_jac[0] = false;
	jac[1] = -1.0;
	num_eval_jac[1] = false;

	return f;
}

/// �gegyenlet (0.a.) ye<0, yv<ze
/**
 * D/100 �tm�r�j� sz�v�sz�lon teltszelv�ny� �raml�s ye-yv nyom�sk�l�nbs�g hat�s�ra
 * @param ye eleje v�zszint
 * @param yv v�ge v�zszint
 * @param mp t�meg�ram
 * @return �gegyenlet �rt�ke
 */
double Csatorna::f_0(double ye, double yv, double mp) {

	//  Matlab k�d:
	//  Df=D/100;
	//  Af=Df^2*pi/4;
	//  ff=yea-yva-0.02*L/Df/2/9.81*Q*abs(Q)/Af^2;
	//  eval_jac=[1 1 1];

	ostringstream strstrm;
	strstrm.str("\n\t\t  f_0: fake cso, D=");
	strstrm << D_fake;

	double v = mp / ro / (D_fake * D_fake * pi / 4);
	double dh = 0.02 * L / D_fake / 2 / 9.81 * v * fabs(v);

	double f = ze + ye - (zv + yv) - dh;

	if ((fabs(f) < 1.e-3) && (D_fake > dia / 1000.)) {
		D_fake /= 2.;
		v = mp / ro / (D_fake * D_fake * pi / 4);
		dh = 0.02 * L / D_fake / 2 / 9.81 * v * fabs(v);
		f = ze + ye - (zv + yv) - dh;
		strstrm << "\n\t\t\t  D=" << D_fake << ", f=" << f;
		//cout<<"\n\t\t\t  "<<nev<<": D="<<D_fake<<", f="<<f;
	}

	strstrm << "\n\t\t\t f = (ze+ye)-(zv+yv)-dp = " << f;
	logfile_write(strstrm.str(), 2);

	jac[0] = +1.0;
	num_eval_jac[0] = false;
	jac[1] = -1.0;
	num_eval_jac[1] = false;

	return f;
}

/// �gegyenlet negat�v t�meg�ram eset�n
/**
 * @param ye eleje v�zszint
 * @param yv v�ge v�zszint
 * @param mp t�meg�ram
 * @return �gegyenlet �rt�ke
 */
double Csatorna::f_1(double ye, double yv, double mp) {
	//  Matlab k�d
	//  eval_jac=[1 1 1];
	//  ykrit=kritikus_szint(Q,D,g,n);
	//  if ye<ykrit
	//      ye=1.01*ykrit;
	//      eval_jac(1)=0;
	//  end
	//  [x,y]=ode15s(@nyf_ode,[0,L],ye,options,Q);
	//  ff=yv-y(end);

	ostringstream strstrm;
	strstrm.str("\n\t\t  f_1: ellentetes aramlas, integralas iranya: -->");
	logfile_write(strstrm.str(), 2);

	double ykrit = kritikus_szint(mp / ro);

	strstrm.str("");
	strstrm << "\n\t\t\t y_krit=" << ykrit; //<<", ye_min="<<ye_min;

	if (ye < ykrit) {
		ye = ykrit * 1.01;
		jac[0] = 0.0;
		num_eval_jac[0] = false;
		strstrm
		<< "\n\t\t\t o !!! Nem lehetseges a kritikus szint alatti bearamlas";
		strstrm << "\n\t\t\t          -> ye=1.01*ykrit=" << ye << " �s eval_jac[0]=0";
	}
	logfile_write(strstrm.str(), 2);

	// Integralas
	double f = yv - ode_megoldo(ye, L / db, 0., mp);
	jac[1] = 1.0;
	num_eval_jac[1] = false;

	strstrm.str("");
	strstrm << "\n\t\t\t f = yv-yv_int = " << f;
	logfile_write(strstrm.str(), 2);

	return f;
}

/// Open-channel flow equation for positive flow rate, y_normal > y_crit
/**
 * @param ye water depth @ x=0
 * @param yv water depth @ x=L
 * @param mp mass flow rate
 * @return error of branch equation
 */
double Csatorna::f_2c(double ye, double yv, double mp) {
	ostringstream strstrm;
	strstrm.str("\n\t\t  f_2c: normal depth not found, integrateing backwards: <--");
	logfile_write(strstrm.str(), 2);

	double ykrit = kritikus_szint(mp / ro);

	strstrm.str("");
	strstrm << "\n\t\t\t y_crit=" << ykrit;

	if (yv < ykrit) {
		yv = ykrit * 1.01;
		jac[1] = 0.0;
		num_eval_jac[1] = false;
		strstrm
		<< "\n\t\t\t o !!! y(x=L) < y_crit, which is not possible.";
		strstrm << "\n\t\t\t          -> setting y(L) = 1.01 * y_crit=" << yv << " and eval_jac[1]=0";
	}
	logfile_write(strstrm.str(), 2);

	// Integrate backwards
	double f = ye - ode_megoldo(yv, -L / db, L, mp);
	jac[0] = 1.0;
	num_eval_jac[0] = false;

	strstrm.str("");
	strstrm << "\n\t\t\t f = y(0)-y(0)_integrated = " << f;
	logfile_write(strstrm.str(), 2);

	return f;
}

/// �gegyenlet pozit�v t�meg�ram eset�n, yn<yk
/**
 * @param ye eleje v�zszint
 * @param yv v�ge v�zszint
 * @param mp t�meg�ram
 * @return �gegyenlet �rt�ke
 */
double Csatorna::f_2aii(double ye, double yv, double mp) {
	//  Matlab k�d:
	//  eval_jac=[1 1 1];
	//  [x,y]=ode15s(@nyf_ode,[0,L],ye,options,Q);
	//
	//  ykrit=kritikus_szint(Q,D,g,n);
	//  if ye<ykrit
	//      ff=ye-normal_szint(Q,D,lejt,n);
	//      eval_jac(2)=0;
	//      fprintf('\n\t\t f_2aii: (b) normal   ye=%+5.3e, yv=%+5.3e, Q=%+5.3e, ff=%+5.3e, yv=%+5.3e',ye,yv,Q,ff,yv);
	//  else
	//      ff=yv-y(end);
	//      fprintf('\n\t\t f_2aii: (c) szubkritikus ye=%+5.3e, yv=%+5.3e, Q=%+5.3e, ff=%+5.3e, yv=%+5.3e',ye,yv,Q,ff,yv);
	//  end

	ostringstream strstrm;
	strstrm.str("\n\t\t  f_2aii: yn<ykrit, integralas iranya: -->");
	logfile_write(strstrm.str(), 2);

	double f, ykrit = kritikus_szint(mp / ro);

	if (ye < ykrit) {
		vector<double> ynv = normal_szint(mp / ro);
		double yn;
		if (ynv.size() == 1)
			yn = ynv[0];
		else {
			if (ynv[0] > ynv[1])
				yn = ynv[1];
			else
				yn = ynv[0];
		}
		f = ye - yn;
		jac[0] = 1.0;
		num_eval_jac[0] = false;
		if (yv < yn)
			jac[1] = 0.0;
		else
			jac[1] = 0.0;
		num_eval_jac[1] = false;
		strstrm.str("");
		strstrm << "\n\t\t\t o !!! ye=" << ye << " < ykrit=" << ykrit << ",  jac[1]=0";
		strstrm << "\n\t\t\t f = ye-yn = " << f;

		// build_res
		if (write_res) {
			if (yv > yn) {
				//cout<<"\n\n"<<nev
				//      <<": 2.a.ii. : vizugras folyadekszint szamitasa...\n";
				double dx = -L / 100., x = L, y1 = yv, f_at_y;

				while ((x > 0) && (y1 > ykrit)) {
					xres.push_back(x);
					yf.push_back((L - x) / L * (ze - zv) + zv);
					yres.push_back(y1);
					//cout<<"\n\t x="<<x<<", yf="<<x/L*(ze-zv)<<", y="<<y1;
					f_at_y = nyf_ode(x, y1, mp);
					y1 += dx * f_at_y;
					x += dx;
				}
				xres.push_back(x);
				yf.push_back((L - x) / L * (ze - zv) + zv);
				yres.push_back(yn);

				xres.push_back(0);
				yf.push_back(ze);
				yres.push_back(yn);
			} else {
				xres.push_back(0);
				xres.push_back(L);
				yres.push_back(yn);
				yres.push_back(yn);
				yf.push_back(ze);
				yf.push_back(zv);
			}
		}
	} else {
		f = yv - ode_megoldo(ye, L / db, 0., mp);
		strstrm.str("");
		strstrm << "\n\t\t\t f = yv-yv_int = " << f;
	}

	logfile_write(strstrm.str(), 2);

	return f;
}

/// Figyelmeztetes k�perny�re �s eredm�nyf�jlba
/**
 * @param fv h�v� f�ggv�ny
 * @param msg �zenet
 */
void Csatorna::warning(string fv, string msg) {
	ostringstream strstrm;
	strstrm.str("");
	strstrm << "\n\n******** FIGYELMEZTETES *********";
	strstrm << "\n\telem neve: " << nev;
	strstrm << "\n\tfuggveny : " << fv;
	strstrm << "\n\tuzenet   : " << msg << "\n\n";
	logfile_write(strstrm.str(), 0);
	cout << strstrm.str();
}

/// Hibauzenet k�perny�re �s eredm�nyf�jlba
/**
 * @param fv h�v� f�ggv�ny
 * @param msg �zenet
 */
void Csatorna::error(string fv, string msg) {
	ostringstream strstrm;
	strstrm.str("");
	strstrm << "\n\n******** HIBA *********";
	strstrm << "\n\telem neve: " << nev;
	strstrm << "\n\tfuggveny : " << fv;
	strstrm << "\n\tuzenet   : " << msg << "\n\n";
	logfile_write(strstrm.str(), 0);
	cout << strstrm.str();
	exit(0);
}

//--------------------------------------------------------------
void Csatorna::Set_dprop(string mit, double mire) {
	//    if (mit=="diameter")
	//      D=mire;
	//    else
	//      {
	cout << endl << "HIBA! Csatorna::Set_dprop(mit), ismeretlen bemenet: mit="
			<< mit << endl << endl;
	//      }
}

//--------------------------------------------------------------
void Csatorna::SetLogFile(string fnev) {
	out_file = nev + ".out";
	ofstream outputFile;
	outputFile.open(out_file.c_str());
	outputFile.close();
}
