/**
* @file Staci.h
* @brief Header file of Staci class
* @author Csaba Hos
* @date 09/20/2017
*/


#include <vector>
#include <memory>
#include "AnyOption.h"
#include "Csomopont.h"
#include "Agelem.h"
#include "nr.h"
#include "epanet_document.h"
#include <string>

struct UmfpackSymbolicDeleter {
    void operator()(void *symbolic) const noexcept;
};

struct UmfpackNumericDeleter {
    void operator()(void *numeric) const noexcept;
};

class Staci
{
public:
    // Non-owning compatibility views. The pointed-to objects are owned by
    // owned_cspok and owned_agelemek for the complete lifetime of this Staci.
    vector<Csomopont *> cspok;
    vector<Agelem *> agelemek;
    void export_connected_nodes();
    Staci(int argc, char *argv[]);
    Staci(string spr_filename);
    ~Staci();
    Staci(const Staci &) = delete;
    Staci &operator=(const Staci &) = delete;
    Staci(Staci &&) noexcept = default;
    Staci &operator=(Staci &&) noexcept = default;
    void print_header();
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
    void save_mod_prop(bool is_general_property);
    void save_mod_prop_all_elements(string property_ID);
    void save_modified_network();
    void export_epanet(const string &filename);
    bool solve_system();
    void solve_residence_time();
    void compute_demand_sensitivity();
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
    void Compute_rank();// WR Computing rank of a node

    // ADATMODOSITASHOZ
    string new_def_file, element_ID, property_ID;
    double newValue;


    void list_all_elements();

    vector<vector<double> > SM_MassFlowRates; // Sensitivity Matrix
    vector<vector<double> > SM_Pressures; // Sensitivity Matrix
    vector<string> SM_row_name; // name of the row element ID
    vector<double> SM_row_sum_MassFlowRates; // row-wise sum of SM elements
    vector<double> SM_row_sum_Pressures; // row-wise sum of SM elements
    vector<string> SM_col_name; // Name of the column element ID
    vector<double> SM_col_sum_MassFlowRates; // column-wise sum of SM elements
    vector<double> SM_col_sum_Pressures; // column-wise sum of SM elements
    vector<double> SM_col_ss_MassFlowRates; // column-wise square sum of SM elements
    vector<double> SM_col_ss_Pressures; // column-wise square sum of SM elements
    vector<double> SM_row_ss_MassFlowRates; // row-wise square sum of SM elements
    vector<double> SM_row_ss_Pressures; // row-wise square sum of SM elements

    double SM_sum_sum_MassFlowRates; // abs sum of the whole SM
    double SM_sum_sum_Pressures; // abs sum of the whole SM
    double SM_ss_con_MassFlowRates; // abs sum of the  whole SM with condition
    double SM_ss_con_Pressures; // abs sum of the  whole SM with condition
    double SM_arp_MassFlowRates; // average relative perturbation that cause 10% change
    double SM_arp_Pressures;

    //void Print_Jacobian(Mat_DP jac);
    //void Print_Jacobian(vector<vector<double> > jac);
    void Print_Jacobian();
    void Print_dfdmu();
    void Print_dxdmu();
    void Compute_dxdmu();
    void Compute_Sensitivity_Matrix(string parameter, int scale);
    void Print_matrix(vector<vector<double> > M);
    void Save_matrix(vector<vector<double> > M, string filename);
    void Save_matrix(vector<vector<double> > M, ofstream &file);
    void Save_Sensitivity();
    void set_do_save_file(const bool save_it)
    {
        do_save_file = save_it;
    }

    double get_sum_of_consumption();
    double get_sum_of_pos_consumption();
    double get_sum_of_neg_consumption();

    double GetMinPipeDiameter(int &idx);
    double GetMaxPipeDiameter(int &idx);
    double GetMinPipeLength(int &idx);
    double GetMaxPipeLength(int &idx);
    double GetSumPipeLength();
    double GetSumPipeVolume();
    double GetMinConsumption(int &idx);
    double GetMaxConsumption(int &idx);
    double GetMinGeoHeight(int &idx);
    double GetMaxGeoHeight(int &idx);
    double GetMaxAbsFlowRate(int &idx);
    double GetMinAbsFlowRate(int &idx);
    double GetMaxPressure(int &idx);
    double GetMinPressure(int &idx);

    void Statistics();
    void Add_edge(const int a, const int b) { // WR: collecting edge vector
        edge.push_back(a);
        edge.push_back(b);
    }
    vector<int> Get_edge() {
        return edge;
    }
    void ProgressBar(int i, int n);// WR Progress bar for long calcs, n: no. of steps, i: current step
    void Avr_absmax_stddev(vector<double> x, double &a, double &m, double &s);// WR Calculates average, max, standard deviation of a vector x
    vector<vector<double> > CSVRead(ifstream &file);// WR Reading doubles from file, separeted with ','
    bool perform_demand_sensitivity_analysis;

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
    std::unique_ptr<AnyOption> opt;
    std::vector<std::unique_ptr<Csomopont> > owned_cspok;
    std::vector<std::unique_ptr<Agelem> > owned_agelemek;
    EpanetDocument epanet_document;
    bool has_epanet_document;

    double transp_dt, ido; //, transp_dt_ki;
    string progress_file;

    void SetInitialParameters();
    void rebuild_network_views();

    void progress_file_write(double percent);
    void get_command_line_options(int argc, char *argv[]);
    string trim(string s, const string drop);
    void Compute_Head_Losses();

    void build_vectors(Vec_DP &x, Vec_DP &f, bool create_sparse_pattern);
    void build_vectors_frozen_Jacobian(Vec_DP &x, Vec_DP &f);

    void pbcg_solver(Vec_I_DP x, Vec_I_DP f);
    bool umfpack_solver(const Vec_DP &x, const Vec_DP &f);
    bool solve_sparse_system(const vector<double> &rhs, vector<double> &solution);
    void build_sparse_pattern();
    int sparse_position(int column, int row) const;
    double sparse_value(int row, int column) const;
    void invalidate_numeric_factorization();

    string iter_info(const Vec_DP &x, const Vec_DP &f, int iter, double e_mp, double e_p);
    void
    compute_error(const Vec_DP &f, double &e_mp, double &e_p, double &e_mp_r,
                  double &e_p_r, bool &konv_ok);
    void update_relax(double e_mp, double e_p, double &e_mp_r, double &e_p_r);

    vector<double>  m_dfdmu; // derivative of eqs. w.r.t. parameter
    vector<double>  m_dxdmu; // derivative of eqs. w.r.t. parameter
    int m_nnz;  /*!< Number of nonzero entries of the Jacobian. */

    void Compute_dfdmu();


    // UMFPACK compressed-column matrix. Its sparsity pattern depends only on
    // network topology and is therefore shared by all Newton and EPS steps.
    vector<int> m_Ap;
    vector<int> m_Ai;
    vector<double> m_Ax;
    vector<int> m_edge_flow_position;
    vector<int> m_edge_start_head_position;
    vector<int> m_edge_end_head_position;
    bool m_sparse_pattern_valid = false;
    bool m_numeric_factorization_valid = false;
    std::unique_ptr<void, UmfpackSymbolicDeleter> m_umfpack_symbolic;
    std::unique_ptr<void, UmfpackNumericDeleter> m_umfpack_numeric;

    string convert_to_hr_min(double seconds);

    // If set to false, the result file will not be saved!
    // Only set false for sensitivity analysis
    bool do_save_file;

    vector<double> col_abs_sum(vector<vector<double> > M);
    vector<double> col_sqr_sum(vector<vector<double> > M);
    vector<double> row_abs_sum(vector<vector<double> > M);
    vector<double> row_sqr_sum(vector<vector<double> > M);
    double abs_sum_sum(vector<vector<double> > M);
    double Cond_sum_sum(vector<vector<double> > M, string par_name);
    double Avr_rel_pert(vector<vector<double> > M, string par_name);
    void rescale_vec(vector<double>& vec);
    void print_worst_iter(const Vec_DP &x, const Vec_DP &f, const int a_debug_level);

    vector<int> edge;

    void Save_Sensitivity_Matrix(string fname);
    double GetAbsMaxCoeff(vector< vector<double> > M);
};
