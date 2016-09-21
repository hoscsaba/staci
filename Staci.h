#include <vector>
#include "AnyOption.h"
#include "Csomopont.h"
#include "Agelem.h"
#include "nr.h"

class Staci
{
public:
    vector<Csomopont *> cspok;
    vector<Agelem *> agelemek;
    void export_connected_nodes();
    Staci(int argc, char *argv[]);
    Staci(string spr_filename);
    ~Staci();
    string Get_out_file()
    {
        return out_file.c_str();
    }
    int Get_debug_level()
    {
        return debug_level;
    }
    void Set_debug_level(int a_debug_level)
    {
        debug_level = a_debug_level;
    }
    void build_system();
    void build_system_old();
    void ini();
    void ini(const Staci *IniStaci);
    void list_system();
    string list_results();
    void save_results(bool conv_reached);
    void save_mod_prop();
    bool solve_system();
    void solve_residence_time();
    void residence_time_step(string& max_ID, double& max_VAL, double& mean_VAL);
    bool solve_system_old();
    void set_up_transport();
    void solve_transport(int mode);
    double get_oldest();
    void transport_step(double dt);
    double teta(double c, const int i);
    void export_to_aisee();
    void post_process_res_file();
    void logfile_write(string msg, int msg_debug_level);
    void set_res_file(string xml_fnev)
    {
        res_file = xml_fnev;
    }
    void set_ini_file(string xml_fnev)
    {
        ini_file = xml_fnev;
    }
    void set_out_file(string fnev)
    {
        out_file = fnev;
    }
    void set_progress_file(string fnev)
    {
        progress_file = fnev;
    }
    string get_def_file()
    {
        return def_file;
    }
    int get_mode()
    {
        return mode;
    }
    double tt_length, cl_length; //, mp_init, p_init;
    void save_transport(int mode);
    void copy_file(const string in_f_nev, const string out_f_nev);
    void set_van_ini(const bool van_e)
    {
        van_ini = van_e;
    }
    void Set_FolyTerf();
    double m_get_dprop();
    void m_set_dprop();
    double get_dprop(string ID, string prop);
    void set_dprop(string ID, string prop, double val);

    // ADATMODOSITASHOZ
    string new_def_file, element_ID, property_ID;
    double newValue;

    void list_all_elements();

    //void Print_Jacobian(Mat_DP jac);
    //void Print_Jacobian(vector<vector<double> > jac);
    void Print_Jacobian();
    void Print_dfdmu();
    void Print_dxdmu();
    void Compute_dxdmu();

    void set_do_save_file(const bool save_it)
    {
        do_save_file = save_it;
    }

private:
    bool van_ini;
    int debug_level, iter_max, mode;
    double p_init, mp_init;
    double e_mp_max, e_p_max, relax, relax_mul;
    double m_relax, m_relax_mul, m_RELAX_MIN, m_RELAX_MAX;
    double FolyMenny;
    string friction_model;
    stringstream m_ss;
    string def_file, out_file, ini_file, aisee_file, res_file;
    AnyOption *opt;

    double transp_dt, ido; //, transp_dt_ki;
    string progress_file;

    void SetInitialParameters();

    void progress_file_write(double percent);
    void get_command_line_options(int argc, char *argv[]);
    string trim(string s, const string drop);
    void Compute_Head_Losses();

    void build_vectors(Vec_DP &x, Vec_DP &f, bool create_sparse_pattern);
    void build_vectors_frozen_Jacobian(Vec_DP &x, Vec_DP &f);

    void nr_solver(Vec_DP x, Vec_DP f);
    void pbcg_solver(Vec_I_DP x, Vec_I_DP f);
    bool umfpack_solver(Vec_I_DP x, Vec_I_DP f);

    string iter_info(Vec_DP x, Vec_DP f, int iter, double e_mp, double e_p);
    void
    compute_error(Vec_DP f, double &e_mp, double &e_p, double &e_mp_r,
                  double &e_p_r, bool &konv_ok);
    void update_relax(double e_mp, double e_p, double &e_mp_r, double &e_p_r);

    vector<vector<double> > m_jac; // Jacobia
    vector<double>  m_dfdmu; // derivative of eqs. w.r.t. parameter
    vector<double>  m_dxdmu; // derivative of eqs. w.r.t. parameter
    int m_nnz;  /*!< Number of nonzero entries of the Jacobian. */

    void Compute_dfdmu();
    

    // UMFPACK
    /* nemnulla elemek */
    vector<vector<bool> > m_is_element_empty;
    /* Ti[k] is row index of entry k, as matrix is scanned columnwise */
    vector<int> m_Ti;
    /* Tj[k] is column index of entry k, as matrix is scanned columnwise */
    vector<int> m_Tj;
    /* value of entry k, as matrix is scanned columnwise */
    vector<double> m_Tx;

    string convert_to_hr_min(double seconds);

    bool do_save_file;

};
