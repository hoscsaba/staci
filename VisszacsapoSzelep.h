#include "Agelem.h"

class VisszacsapoSzelep : public Agelem {
private:
	/// Visszafolyás esetén az ellenállásteényezõ
	double dzeta_e, dzeta_v;

public:
	VisszacsapoSzelep(const string nev, const string cspenev,
	                  const string cspvnev, const double a_ro, const double Aref, const double veszt_e,
	                  const double veszt_v, const double a_mp);
	~VisszacsapoSzelep();
	string Info();
	double f(vector<double>);
	vector<double> df(vector<double>);
	void Ini(int mode, double value);
	void Set_dprop(string mit, double mire);
	string GetType() {
		return "VisszacsapoSzelep";
	}
	double Get_dprop(string mit);
};
