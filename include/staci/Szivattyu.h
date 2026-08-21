#include "Agelem.h"

class Szivattyu: public Agelem {
private:
	vector<double> q, H, p;
	double a, b, c;
	int fokszam;
	double mer_szorzo;
public:
	Szivattyu(const string &nev, const string &a_cspe_nev,
			const string &a_cspv_nev, double a_ro, double Aref, const vector<double> &q,
			const vector<double> &H, double a_mp);
	~Szivattyu() override;
	string Info() override;
	double f(const vector<double> &state) override;
	vector<double> df(const vector<double> &state) override;
	void Ini(int mode, double value) override;
	void Set_dprop(const string &property, double value) override;
	double Get_dprop(const string &property) override;
	string_view GetType() const noexcept override {
		return "Szivattyu";
	}
	double PumpCharCurve(double q);
	double Get_PumpHeadAt(double q) {
		return PumpCharCurve(q);
	}
	vector<double> GetCurveFlowM3PerHour() const {
		vector<double> result = q;
		for (double &value : result) value *= 3600.0;
		return result;
	}
	const vector<double> &GetCurveHead() const {
		return H;
	}
};
