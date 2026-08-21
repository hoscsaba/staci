#include "Agelem.h"

class VisszacsapoSzelep : public Agelem {
private:
	/// Visszafolyás esetén az ellenállásteényezõ
	double dzeta_e, dzeta_v;

public:
	VisszacsapoSzelep(const string &nev, const string &cspenev,
	                  const string &cspvnev, const double a_ro, const double Aref, const double veszt_e,
	                  const double veszt_v, const double a_mp);
	~VisszacsapoSzelep() override;
	string Info() override;
	double f(const vector<double> &state) override;
	vector<double> df(const vector<double> &state) override;
	void Ini(int mode, double value) override;
	void Set_dprop(const string &property, double value) override;
	string_view GetType() const noexcept override {
		return "VisszacsapoSzelep";
	}
	double Get_dprop(const string &property) override;
};
