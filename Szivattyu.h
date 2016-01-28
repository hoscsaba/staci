#include "Agelem.h"

class Szivattyu: public Agelem {
private:
	vector<double> q, H, p;
	double a, b, c;
	int fokszam;
	double mer_szorzo;
public:
	Szivattyu(const string nev, const string a_cspe_nev,
			const string a_cspv_nev, const double a_ro, const double Aref, const vector<double> q,
			const vector<double> H, const double a_mp);
	~Szivattyu();
	string Info();
	double f(vector<double>);
	vector<double> df(vector<double>);
	void Ini(int mode, double value);
	void Set_dprop(string mit, double mire);
	string GetType() {
		return "Szivattyu";
	}
	double PumpCharCurve(double q);
	double Get_PumpHeadAt(double q) {
		return PumpCharCurve(q);
	}
};
