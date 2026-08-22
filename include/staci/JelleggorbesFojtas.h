#include "Agelem.h"

enum class EpanetTcvStatus
{
  Active,
  Open,
  Closed
};

class JelleggorbesFojtas : public Agelem
{
private:
  vector<double> e, zeta, p;
  double allas, veszt;
  int fokszam;
  double headloss;
  bool epanet_tcv = false;
  double epanet_tcv_setting = 0.0;
  double epanet_tcv_minor_loss = 0.0;
  EpanetTcvStatus epanet_tcv_status = EpanetTcvStatus::Active;
  void SetConstantEpanetLoss(double loss_coefficient);
public:
  JelleggorbesFojtas(const string &nev,
   const string &a_cspe_nev,
   const string &a_cspv_nev,
   const double a_ro,
   const double Aref,
   const vector<double> &e,
   const vector<double> &zeta,
   const double allas,
   const double a_mp);
  void Update_zeta();
  ~JelleggorbesFojtas() override;
  string Info() override;
  double f(const vector<double> &state) override;
  vector<double> df(const vector<double> &state) override;
  void Ini(int mode, double value) override;
  void Set_dprop(const string &property, double value) override;
  double Get_dprop(const string &property) override;
  void SetEpanetTcvMetadata(double setting, double minor_loss,
                            EpanetTcvStatus status = EpanetTcvStatus::Active);
  void SetEpanetTcvSetting(double setting);
  void SetEpanetTcvStatus(EpanetTcvStatus status);
  bool HasEpanetTcvMetadata() const noexcept { return epanet_tcv; }
  bool CanExportAsEpanetTcv() const noexcept;
  double GetEpanetTcvSetting() const noexcept;
  double GetEpanetTcvMinorLoss() const noexcept;
  EpanetTcvStatus GetEpanetTcvStatus() const noexcept;
  string_view GetType() const noexcept override
  {
    return "JelleggorbesFojtas";
  }
};

