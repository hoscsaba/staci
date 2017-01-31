// #include <stdio.h>
#include <vector>
#include <sstream>
#include "Staci.h"
#include "xmlParser.h"
#include <Eigen/Dense>

#include <ga/ga.h>
#include <ga/GASimpleGA.h>
#include <ga/GAListGenome.h>

#include <igraph/igraph.h>

float Objective(GAGenome &);

using namespace std;
using namespace Eigen;

#define cout STD_COUT

// These global variables are coming from the datafile
string fname;
int popsize, ngen;
double pcross, pmut;
string weight_type;
string weight_type_mod;

Staci *wds;
int n_n, n_p, n_comm;
int fcount;
double f_ce;
vector<int> pipe_idx;
VectorXd p, W;
MatrixXd A;
double sumW;

int n_c;
double Qmod;

bool info;

vector< vector<double> > jac;

void LoadMatrices(MatrixXd &A, VectorXd &W, VectorXd &p, string weight_type);

void LoadSystem(string fname);

void save_state();
void save_state(const igraph_vector_t *v);

void logfile_write(string msg, int debug_level);

int Mutator(GAGenome & c, float prob_mut);
int Q_Mutator(GAGenome & c, float prob_mut);
int D_Mutator(GAGenome & c, float prob_mut);

void Initializer(GAGenome & c);
void Q_Initializer(GAGenome & c);
void D_Initializer(GAGenome & c);

float Objective(GAGenome &);
float Q_Objective(GAGenome &);
float D_Objective(GAGenome &);

void Optimize();
void Q_Optimize();
void D_Optimize();

void Load_Settings();

void Save_Membership(GAGenome & c);
void Save_Membership2();

double GetAbsMaxCoeff(vector< vector<double> > M);

int global_debug_level;

vector<int> best;
double best_Q;
double obj_offset;

stringstream strstrm;
string logfilename;

string obj_type;

void print(igraph_t *g);
int print_vector(const igraph_vector_t *v);
int print_matrix(const igraph_matrix_t *m);
int igraph_community_eigenvector();
void copy_to_best(const igraph_vector_t *v);

void PerformSensitivityAnalysis(bool is_edge_prop, string par, string fname);

vector< vector<double> > SM;
vector< vector<double> > SM_MFR ;
vector< vector<double> > SM_PR ;

struct val_and_ID {
    double val;
    string ID;
    bool operator>(const val_and_ID& rhs) { return val > rhs.val; }
    bool operator>=(const val_and_ID& rhs) { return val >= rhs.val; }
    bool operator<(const val_and_ID& rhs) { return val < rhs.val; }
    bool operator<=(const val_and_ID& rhs) { return val <= rhs.val; }
};

bool comparison_function1(const val_and_ID& lhs, const val_and_ID& rhs ) { return lhs.val > rhs.val; }

struct val_and_ID_and_comm {
    double val;
    string ID;
    int comm;
    bool operator>(const val_and_ID_and_comm& rhs) { return val > rhs.val; }
    bool operator>=(const val_and_ID_and_comm& rhs) { return val >= rhs.val; }
    bool operator<(const val_and_ID_and_comm& rhs) { return val < rhs.val; }
    bool operator<=(const val_and_ID_and_comm& rhs) { return val <= rhs.val; }
};

bool comparison_function2(const val_and_ID_and_comm& lhs, const val_and_ID_and_comm& rhs ) { return lhs.val > rhs.val; }


int main(int argc, char **argv) {

    best_Q = -100.;
    obj_offset = 1.0;
    info = false;
    fcount = 0;

    // Clear logfile
    ofstream ofs(logfilename.c_str(), std::ios::out | std::ios::trunc);
    ofs.close();

    // Load settings
    Load_Settings();

    // Load system
    LoadSystem(fname);

    // Sensitivity matrix
    if ((obj_type == "A-optimality") || (obj_type == "D-optimality")) {
        if (weight_type == "friction_coeff") {
            PerformSensitivityAnalysis(true /*is_edge_prop*/, "friction_coeff", "sensitivity_matrix_friction_coeff.csv");
        }
        else if (weight_type == "diameter") {
            PerformSensitivityAnalysis(true/*is_edge_prop*/, "diameter", "sensitivity_matrix_diameter.csv");
        }
        else if (weight_type == "demand") {
            PerformSensitivityAnalysis(false /*is_edge_prop*/, "demand", "sensitivity_matrix_demand.csv");
        }
        else {
            cout << endl << endl << "ERROR: illegal weight_type: >" << weight_type << "<" << endl;
            exit(-1);
        }
    }

    for (int i = 0; i < wds->SM_Pressures.size(); i++) {
        vector<double> tmp = wds->SM_Pressures.at(i);
        SM.push_back(tmp);
        tmp.clear();
    }

    double maxSM = GetAbsMaxCoeff(SM);
    for (int i = 0; i < SM.size(); i++)
        for (int j = 0; j < SM.at(0).size(); j++)
            SM.at(i).at(j) /= maxSM;

    // Number of funeval for complete enumeration
    f_ce = pow((double) n_comm, (double) n_n);

    // Load matrices
    LoadMatrices(A, W, p, weight_type);

    // let it go...
    Optimize();

    return 0;
}

void Save_Membership(GAGenome & c) {
    GA1DArrayAlleleGenome<int> &genome = (GA1DArrayAlleleGenome<int> &) c;

    stringstream strstrm;
    strstrm.str("");
    strstrm << endl << "Saving membership to membership.txt...";

    FILE * pFile;

    pFile = fopen ("membership.txt", "w");
    fprintf(pFile, "n_comm  : %d :\n", n_comm);
    fprintf(pFile, "n_nodes : %d\n", n_n);
    int node_count = 0;
    vector<int> n_nonempty(n_comm, 0);
    for (int i = 0; i < n_comm; i++) {
        fprintf(pFile, "#%d; ", i);
        strstrm << endl << "comm. #" << i << ": ";
        for (int j = 0 ; j < n_n; j++) {
            if (genome.gene(j) == i) {
                fprintf (pFile, "%s; ", wds->cspok.at(j)->Get_nev().c_str());
                strstrm << wds->cspok.at(j)->Get_nev().c_str() << " ";
                node_count++;
                n_nonempty.at(i) = 1;
            }
        }
        fprintf (pFile, "\n");
    }
    int sum_nonempty = 0;
    for (unsigned int i = 0; i < n_comm; i++)
        sum_nonempty += n_nonempty.at(i);

    fprintf(pFile, "non empty comm: %d\n", sum_nonempty);
    fprintf(pFile, "    empty comm: %d\n", n_c - sum_nonempty);

    fclose (pFile);

    if (n_n != node_count)
        strstrm << " \n\n???? ERROR n_n=" << n_n << " ?= node_count=" << node_count << " WRF????\n\n";
    else
        strstrm << " done.\n";
    logfile_write(strstrm.str(), 0);

}

void Save_Membership2() {
    // GA1DArrayAlleleGenome<int> &genome = (GA1DArrayAlleleGenome<int> &) c;

    stringstream strstrm;
    strstrm.str("");
    strstrm << endl << "Saving membership to membership.txt...";

    FILE * pFile;

    pFile = fopen ("membership.txt", "w");
    fprintf(pFile, "n_comm  : %d :\n", n_comm);
    fprintf(pFile, "n_nodes : %d\n", n_n);
    int node_count = 0;
    vector<int> n_nonempty(n_comm, 0);
    for (int i = 0; i < n_comm; i++) {
        fprintf(pFile, "#%d; ", i);
        strstrm << endl << "comm. #" << i << ": ";
        for (int j = 0 ; j < n_n; j++) {
            if (best.at(j) == i) {
                fprintf (pFile, "%s; ", wds->cspok.at(j)->Get_nev().c_str());
                strstrm << wds->cspok.at(j)->Get_nev().c_str() << " ";
                node_count++;
                n_nonempty.at(i) = 1;
            }
        }
        fprintf (pFile, "\n");
    }
    int sum_nonempty = 0;
    for (unsigned int i = 0; i < n_comm; i++)
        sum_nonempty += n_nonempty.at(i);

    fprintf(pFile, "non empty comm: %d\n", sum_nonempty);
    fprintf(pFile, "    empty comm: %d\n", n_c - sum_nonempty);

    fclose (pFile);

    if (n_n != node_count)
        strstrm << " \n\n???? ERROR n_n=" << n_n << " ?= node_count=" << node_count << " WRF????\n\n";
    else
        strstrm << " done.\n";
    logfile_write(strstrm.str(), 0);

}


double GetAbsMaxCoeff(vector< vector<double> > M) {
    double maxval = -1.;
    for (unsigned int i = 0; i < M.size(); i++) {
        for (unsigned int j = 0; j < M.at(i).size(); j++) {
            if (maxval < abs(M.at(i).at(j)))
                maxval = abs(M.at(i).at(j));
        }
    }
    return maxval;
}

void PerformSensitivityAnalysis(bool is_edge_prop, string par, string fname) {
    wds->Compute_Sensitivity_Matrix(par, 0);

    int debug = false;

    if (debug) {
        cout << endl << "Entering PerformSensitivityAnalysis()..." << flush;
        cin.get();
    }

    SM_MFR = wds->SM_MassFlowRates;
    SM_PR  = wds->SM_Pressures;
    int n_edges = wds->agelemek.size();
    int n_nodes = wds->cspok.size();
    // cout << endl << "SM_PR:";
    // for (unsigned int i = 0; i < SM_PR.size(); i++) {
    //     cout << endl;
    //     for (unsigned int j = 0; j < SM_PR.at(i).size(); j++) {
    //         printf("%+5.3e ", SM_PR.at(i).at(j));
    //     }
    // }

    double max_SM_MFR = GetAbsMaxCoeff(SM_MFR);
    double max_SM_PR  = GetAbsMaxCoeff(SM_PR);

    if (is_edge_prop) {
        for (int i = 0 ; i < n_edges; i++) {
            for (int j = 0 ; j < n_edges; j++)
                SM_MFR.at(i).at(j) /= max_SM_MFR;
            for (int j = 0 ; j < n_nodes; j++)
                SM_PR.at(i).at(j) /= max_SM_PR;
        }
    }
    else {
        for (int i = 0 ; i < n_nodes; i++) {

            for (int j = 0 ; j < n_edges; j++)
                SM_MFR.at(i).at(j) /= max_SM_MFR;
            for (int j = 0 ; j < n_nodes; j++)
                SM_PR.at(i).at(j) /= max_SM_PR;
        }
    }

    FILE * pFile;

    pFile = fopen (fname.c_str(), "w");
    fprintf(pFile, "parameter: %s\n", par.c_str());
    fprintf(pFile, "n_edges : %d\n", n_edges);
    fprintf(pFile, "n_nodes : %d\n", n_nodes);
    fprintf(pFile, "n_cols  : %d ( = n_edges + n_nodes )\n", n_edges + n_nodes);
    fprintf(pFile, "n_rows  : %d ( = n_edges )\n", n_edges);
    fprintf(pFile, "column names; ");
    for (int j = 0 ; j < n_edges; j++)
        fprintf (pFile, "%s; ", wds->agelemek.at(j)->Get_nev().c_str());
    for (int j = 0 ; j < n_nodes; j++)
        fprintf (pFile, "%s; ", wds->cspok.at(j)->Get_nev().c_str());
    fprintf (pFile, "\n");

    if (is_edge_prop) {
        for (int i = 0 ; i < n_edges; i++) {
            fprintf(pFile, "%s @ %s; ", par.c_str(), wds->agelemek.at(i)->Get_nev().c_str());
            for (int j = 0 ; j < n_edges; j++)
                fprintf (pFile, "%+7.5e; ", SM_MFR.at(i).at(j) );
            for (int j = 0 ; j < n_nodes; j++)
                fprintf (pFile, "%+7.5e; ", SM_PR.at(i).at(j));
            fprintf (pFile, "\n");
        }
    }
    else {
        for (int i = 0 ; i < n_nodes; i++) {
            fprintf(pFile, "%s @ %s; ", par.c_str(), wds->cspok.at(i)->Get_nev().c_str());
            for (int j = 0 ; j < n_edges; j++)
                fprintf (pFile, "%+7.5e; ", SM_MFR.at(i).at(j));
            for (int j = 0 ; j < n_nodes; j++)
                fprintf (pFile, "%+7.5e; ", SM_PR.at(i).at(j));
            fprintf (pFile, "\n");
        }
    }

    // Add total sensitivities
    vector<val_and_ID> v_edges;
    vector<val_and_ID> v_nodes;
    val_and_ID s;

    if (is_edge_prop) {
        for (int i = 0 ; i < n_edges; i++) {
            double tmp = 0.;
            for (int j = 0 ; j < n_edges; j++)
                tmp += SM_MFR.at(j).at(i) * SM_MFR.at(j).at(i);
            s.val =  sqrt(tmp);
            s.ID = wds->agelemek.at(i)->Get_nev().c_str();
            v_edges.push_back(s);
        }
        for (int i = 0 ; i < n_nodes; i++) {
            double tmp = 0.;
            for (int j = 0 ; j < n_edges; j++) {
                tmp += SM_PR.at(j).at(i) * SM_PR.at(j).at(i);
            }
            s.val = sqrt(tmp);
            s.ID = wds->cspok.at(i)->Get_nev().c_str();
            v_nodes.push_back(s);
        }
    }
    else {
        for (int i = 0 ; i < n_edges; i++) {
            double tmp = 0.;
            for (int j = 0 ; j < n_nodes; j++)
                tmp += SM_MFR.at(j).at(i) * SM_MFR.at(j).at(i) ;
            s.val =  sqrt(tmp);
            s.ID = wds->agelemek.at(i)->Get_nev().c_str();
            v_edges.push_back(s);
        }
        for (int i = 0 ; i < n_nodes; i++) {
            double tmp = 0.;
            for (int j = 0 ; j < n_nodes; j++)
                tmp += SM_PR.at(j).at(i) * SM_PR.at(j).at(i) ;
            s.val = sqrt(tmp);
            s.ID = wds->cspok.at(i)->Get_nev().c_str();
            v_nodes.push_back(s);
        }
    }

    // fprintf (pFile, "\nColumn-wise sum of absolute values:");
    // for (int i = 0; i < n_edges; i++)
    //     fprintf(pFile, "\n %10s: %7.5e", v_edges.at(i).ID.c_str(), v_edges.at(i).val);
    // fprintf(pFile, "\n");
    // for (int i = 0; i < n_nodes; i++)
    //     fprintf(pFile, "\n %10s: %7.5e", v_nodes.at(i).ID.c_str(), v_nodes.at(i).val);

    sort(v_edges.begin(), v_edges.end(), comparison_function1);
    sort(v_nodes.begin(), v_nodes.end(), comparison_function1);

    fprintf (pFile, "\nColumn-wise L2 norm of sensitivities, after sorting:");
    for (int i = 0; i < n_edges; i++)
        fprintf(pFile, "\n %10s; %7.5e", v_edges.at(i).ID.c_str(), v_edges.at(i).val);
    fprintf(pFile, "\n");
    for (int i = 0; i < n_nodes; i++)
        fprintf(pFile, "\n %10s; %7.5e", v_nodes.at(i).ID.c_str(), v_nodes.at(i).val);

    fclose (pFile);

    if (debug) {
        cout << endl << "Leaving PerformSensitivityAnalysis()." << endl << flush;
        cin.get();
    }
}

void Optimize() {
    cout << endl << "RUNNING OPTMIZATION..." << endl;

    if (obj_type == "modularity") {
        Q_Optimize();
    }
    else if ((obj_type == "D-optimality") || (obj_type == "A-optimality")) {
        D_Optimize();
    }
    else {
        cout << endl << "ERROR: !!bad obj_type value: " << obj_type << " !!!" << endl << endl;
        exit(-1);
    }
}

void Q_Optimize() {
    GAAlleleSet<int> alleles;
    for (int i = 0; i < n_comm; i++)
        alleles.add(i);

    GA1DArrayAlleleGenome<int> genome(n_n, alleles, Objective);

    genome.initializer(Initializer);
    genome.mutator(Mutator);

    GASteadyStateGA ga(genome);

    ga.set(gaNpopulationSize, popsize);
    ga.set(gaNpCrossover, pcross);
    ga.set(gaNpMutation, pmut);
    ga.set(gaNnGenerations, ngen);
    ga.set(gaNscoreFrequency, 100);
    ga.set(gaNflushFrequency, 100);
    ga.set(gaNselectScores, GAStatistics::Maximum | GAStatistics::Minimum | GAStatistics::Mean);
    ga.set(gaNscoreFilename, "bog.dat");

    ga.evolve();
    cout << endl << "done." << endl;

    // genome = ga.statistics().bestIndividual();
    //double best_obj = Objective(genome);

    // strstrm.str("");
    // strstrm << endl << "BEST SOLUTION FOUND:" << endl;

    // for (unsigned int i = 0; i < n_n; i++)
        // strstrm << genome.gene(i) << " ";
    // strstrm << endl;
    // logfile_write(strstrm.str(), 0);

    // info = true;
    // double best_obj = Objective(genome);

    strstrm << endl << "Some more info:";
    strstrm << endl << "\t datafile   : " << fname;
    strstrm << endl << "\t # of nodes : " << n_n;
    strstrm << endl << "\t # of pipes : " << n_p;
    strstrm << endl << "\t D          : " << (wds->GetMinPipeDiameter()) << " ... " << (wds->GetMaxPipeDiameter())
            << " m";
    strstrm << endl << "\t L          : " << (wds->GetMinPipeLength()) << " ... " << (wds->GetMaxPipeLength()) << " m";
    strstrm << endl << "\t sum. L     : " << (wds->GetSumPipeLength()) << " m ";
    strstrm << endl << "\t max. cons. : " << (wds->GetMaxConsumption()) << " kg/s = "
            << 3.6 * (wds->GetMaxConsumption())
            << " m3/h" << endl;

    logfile_write(strstrm.str(), 0);

    // Save_Membership(genome);
    Save_Membership2();


    // Add total sensitivities
    bool is_edge_prop;

    if (weight_type_mod == "friction_coeff") {
        is_edge_prop = true;
        PerformSensitivityAnalysis(true /*is_edge_prop*/, "friction_coeff", "sensitivity_matrix_friction_coeff.csv");
    }
    else if (weight_type_mod == "diameter") {
        is_edge_prop = true;
        PerformSensitivityAnalysis(true/*is_edge_prop*/, "diameter", "sensitivity_matrix_diameter.csv");
    }
    else if (weight_type_mod == "demand") {
        is_edge_prop = false;
        PerformSensitivityAnalysis(false /*is_edge_prop*/, "demand", "sensitivity_matrix_demand.csv");
    }
    else {
        cout << endl << endl << "ERROR: illegal weight_type: >" << weight_type << "<" << endl;
        exit(-1);
    }


    // vector<val_and_ID_and_comm> v_edges;
    vector<val_and_ID_and_comm> v_nodes;
    val_and_ID_and_comm s;
    // int n_edges = wds->agelemek.size();
    int n_nodes = wds->cspok.size();

    for (int i = 0 ; i < n_nodes; i++) {
        double tmp = 0.;
        for (int j = 0 ; j < SM_PR.size(); j++) {
            tmp += SM_PR.at(j).at(i) * SM_PR.at(j).at(i);
        }
        s.val = sqrt(tmp);
        s.ID = wds->cspok.at(i)->Get_nev().c_str();
        s.comm = best.at(i);//genome.gene(i);
        v_nodes.push_back(s);
    }

    sort(v_nodes.begin(), v_nodes.end(), comparison_function2);

    strstrm.str("");
    strstrm << "\n\nColumn-wise L2 norm of sensitivities, after sorting:";
    for (int i_comm = 0; i_comm < n_comm; i_comm++) {
        strstrm << "\n\nCommunity #" << i_comm;
        for (int i = 0; i < n_nodes; i++) {
            if (v_nodes.at(i).comm == i_comm)
                strstrm << endl << "\t " << v_nodes.at(i).ID.c_str() << ";\t" << v_nodes.at(i).val;
        }
    }
    logfile_write(strstrm.str(), 0);
}

void D_Optimize() {

    vector<int> best(n_comm, 0);

    GAAlleleSet<int> alleles;
    for (int i = 0; i < n_n; i++)
        alleles.add(i);

    GA1DArrayAlleleGenome<int> genome(n_comm, alleles, Objective);

    genome.initializer(Initializer);
    genome.mutator(Mutator);

    GASteadyStateGA ga(genome);
    if (obj_type == "A-optimality") {
        ga.minimize();
        best_Q = +100.;
    }

    ga.set(gaNpopulationSize, popsize);
    ga.set(gaNpCrossover, pcross);
    ga.set(gaNpMutation, pmut);
    ga.set(gaNnGenerations, ngen);
    ga.set(gaNscoreFrequency, 100);
    ga.set(gaNflushFrequency, 100);
    ga.set(gaNselectScores, GAStatistics::Maximum | GAStatistics::Minimum | GAStatistics::Mean);
    ga.set(gaNscoreFilename, "bog.dat");
    // ga.elitism();

    ga.evolve();
    cout << endl << "done." << endl;

    // genome = ga.statistics().bestIndividual();

    // strstrm.str("");
    // strstrm << endl << "BEST SOLUTION FOUND:" << endl;


    vector<int> best_sorted;
    for (unsigned int i = 0; i < n_comm; i++)
        best_sorted.push_back(best.at(i));
        // best_sorted.push_back(genome.gene(i));
    sort(best_sorted.begin(), best_sorted.begin() + best_sorted.size());

    for (unsigned int i = 0; i < n_comm; i++)
        strstrm << "\t" << wds->cspok.at(best_sorted.at(i))->Get_nev();
    // strstrm << "\tNode #" << genome.gene(i) << " is " << wds->cspok.at(genome.gene(i))->Get_nev() << endl;
    logfile_write(strstrm.str(), 0);

    // info = true;
    // GA1DArrayAlleleGenome<int> tmp = (GA1DArrayAlleleGenome<int> &) c;
    // double best_obj = Objective(genome);

    strstrm.str("");
    strstrm << endl << "Some more info:";
    strstrm << endl << "\t datafile   : " << fname;
    strstrm << endl << "\t # of nodes : " << n_n;
    strstrm << endl << "\t # of pipes : " << n_p;
    strstrm << endl << "\t D          : " << (wds->GetMinPipeDiameter()) << " ... " << (wds->GetMaxPipeDiameter()) << " m";
    strstrm << endl << "\t L          : " << (wds->GetMinPipeLength()) << " ... " << (wds->GetMaxPipeLength()) << " m";
    strstrm << endl << "\t sum. L     : " << (wds->GetSumPipeLength()) << " m ";
    strstrm << endl << "\t max. cons. : " << (wds->GetMaxConsumption()) << " kg/s = " << 3.6 * (wds->GetMaxConsumption()) << " m3/h" << endl;

    logfile_write(strstrm.str(), 0);

// Save_Membership(genome);

}

void Initializer(GAGenome & c) {
    if (obj_type == "modularity") {
        Q_Initializer(c);
    }
    else if ((obj_type == "D-optimality") || (obj_type == "A-optimality")) {
        D_Initializer(c);
    }
    else {
        cout << endl << "ERROR: !!bad obj_type value: " << obj_type << " !!!" << endl << endl;
        exit(-1);
    }
}

void Q_Initializer(GAGenome & c) {
    GA1DArrayAlleleGenome<int> &child = (GA1DArrayAlleleGenome<int> &) c;

    bool ini_info = false;

    if (n_comm == 2) {
        igraph_community_eigenvector();
        for (int i = 0; i < n_n; i++)
            child.gene(i, best.at(i));
        if (ini_info) {
            strstrm.str();
            strstrm << endl << endl << "igraph_Initializer with c_comm=" << n_comm << "):" << endl;
            for (int i = 0; i < child.length(); i++)
                strstrm << child.gene(i) << " ";
            logfile_write(strstrm.str(), 0);
        }
    }
    else {
        int wmax = n_comm - 1;
        int wmin = 0;
        for (int i = 0; i < n_n; i++)
            child.gene(i, wmin + (rand() % (int)(wmax - wmin + 1)));
        if (ini_info) {
            strstrm.str();
            strstrm << endl << endl << "Random initializer with c_comm=" << n_comm << "):" << endl;
            for (int i = 0; i < n_n; i++)
                strstrm << child.gene(i) << " ";
            logfile_write(strstrm.str(), 0);
        }
    }
}

void D_Initializer(GAGenome & c) {
    GA1DArrayAlleleGenome<int> &child = (GA1DArrayAlleleGenome<int> &) c;

    int wmax = n_n - 1;
    int wmin = 0;
    int tmp;
    vector<int> already_used;
    for (int i = 0; i < n_comm; i++) {
        if (already_used.size() > 0) {
            // cout << endl << "already_used: ";
            // for (int k = 0; k < already_used.size(); k++)
            // cout << already_used.at(k) << " ";
            bool is_ok = false;
            // for (int k = 0; k < 100; k++) {
            while (!is_ok) {
                tmp = wmin + (rand() % (int)(wmax - wmin + 1));
                // cout << endl << "\t tmp=" << tmp;
                is_ok = true;
                for (unsigned int j = 0; j < already_used.size(); j++) {
                    if (already_used.at(j) == tmp) {
                        is_ok = false;
                        // cout << " = already_used(" << j << ")";
                        // cin.get();
                    }
                }
                // cin.get();
            }
        }
        else {
            tmp = wmin + (rand() % (int)(wmax - wmin + 1));
        }

        already_used.push_back(tmp);
        child.gene(i, tmp);
        // cout<<endl<<" --- "<<tmp;
    }

    // strstrm.str();
    // strstrm << endl << "Random initializer with c_comm=" << n_comm << ": ";
    // for (int i = 0; i < child.size(); i++)
    //     strstrm << child.gene(i) << " ";
    // // cin.get();

    logfile_write(strstrm.str(), 0);

}


int Mutator(GAGenome & c, float prob_mut) {
    if (obj_type == "modularity") {
        Q_Mutator(c, prob_mut);
    }
    else if ((obj_type == "D-optimality") || (obj_type == "A-optimality")) {
        D_Mutator(c, prob_mut);
    }
    else {
        cout << endl << "ERROR: !!bad obj_type value: " << obj_type << " !!!" << endl << endl;
        exit(-1);

    }
}

int Q_Mutator(GAGenome & c, float prob_mut) {
    GA1DArrayAlleleGenome<int> &child = (GA1DArrayAlleleGenome<int> &) c;

    // Copy of child
    vector<int> copy_of_child(child.size());
    for (unsigned int i = 0; i < child.size(); i++)
        copy_of_child.at(i) = child.gene(i);

    bool mut_info = false;
    info = mut_info;
    int mut_round = 0;
    bool next_round = true;
    double obj_before = 0., obj_after = 1.;
    int nMut = 0;

    if (prob_mut <= 0.0) return (0);

    if (mut_info) {
        cout << endl << endl << "MUTATION:" << endl << endl << "Initial gene:" << flush;
        // cin.get();
    }

    int wmax = 10;
    int wmin = 0;
    int walk_steps_max = wmin + (rand() % (int)(wmax - wmin + 1));
    // cout<<endl<<"walk_steps_max:"<<walk_steps_max;

    obj_before = Objective(child);
    while (next_round && (mut_round < 10)) {
        if (mut_info)
            cout << endl << endl << "ROUND " << mut_round;

        int n_c = 0;
        vector<bool> cut_idx(n_n, false);
        vector<int> cut_pair(n_n, 0);
        for (unsigned int i = 0; i < n_p; i++) {
            int idx_n1 = wds->agelemek.at(pipe_idx.at(i))->Get_Cspe_Index();
            int idx_n2 = wds->agelemek.at(pipe_idx.at(i))->Get_Cspv_Index();

            if (child.gene(idx_n1) != child.gene(idx_n2)) {
                n_c++;
                cut_idx.at(idx_n1) = true;
                cut_idx.at(idx_n2) = true;
                cut_pair.at(idx_n1) = idx_n2;
                cut_pair.at(idx_n2) = idx_n1;
                if (mut_info) {
                    cout << endl << "\tCut: " << idx_n1 << "(" << child.gene(idx_n1);
                    cout << ") --x-- " << idx_n2 << "(" << child.gene(idx_n2) << ")";
                }
            }
        }
        if (mut_info)
            cout << endl;

        for (unsigned int i = 0; i < n_n; i++) {
            if (cut_idx.at(i)) {
                if (GAFlipCoin(0.5)) {
                    int idx_from = i;
                    int idx_to   = cut_pair.at(idx_from);
                    int other_node_community = child.gene(idx_from);
                    cut_idx.at(i) = false;
                    cut_idx.at(cut_pair.at(i)) = false;

                    child.gene(idx_to, other_node_community);
                    if (mut_info) {
                        cout << endl << "edge from node # " << idx_from << " to node # " << idx_to << " : 2. node moved to comm." << other_node_community;
                    }

                    for (int walk_steps = 0; walk_steps < walk_steps_max; walk_steps++) {
                        // Now try to find a pipe connected to node_to NOT terminating at node_from
                        bool found = false;
                        for (int i = 0; i < pipe_idx.size(); i++) {
                            int idx_n1 = wds->agelemek.at(pipe_idx.at(i))->Get_Cspe_Index();
                            int idx_n2 = wds->agelemek.at(pipe_idx.at(i))->Get_Cspv_Index();
                            if ((idx_n1 == idx_to) && (idx_n2 != idx_from)) {
                                if (mut_info) {
                                    cout << endl << "\t" << wds->agelemek.at(pipe_idx.at(i))->Get_nev() << ": " << idx_n1 << " <-- " << idx_n2;
                                }
                                found = true;
                                idx_from = idx_n1;
                                idx_to = idx_n2;

                                break;
                            }
                            if ((idx_n2 == idx_to) && (idx_n1 != idx_from)) {
                                if (mut_info) {
                                    cout << endl << "\t" << wds->agelemek.at(pipe_idx.at(i))->Get_nev() << ": " << idx_n1 << " --> " << idx_n2;
                                }
                                found = true;
                                idx_from = idx_n2;
                                idx_to = idx_n1;

                                break;
                            }
                        }
                        if (found) {
                            child.gene(idx_to, other_node_community);
                            if (mut_info) {
                                cout << endl << "edge from node # " << idx_from << " to node # " << idx_to << " : 2. node moved to comm." << other_node_community;
                            }
                            nMut++;
                        }
                        else {
                            if (mut_info) {
                                cout << endl << "No appropriate edge was found, stopping.";
                                break;
                            }
                        }
                    } if (mut_info) {
                        cin.get();
                    }
                }
            }
        }
        if (mut_info)
            cout << endl;

        obj_after = Objective(child);

        if (obj_after > obj_before) {
            next_round = true;
            mut_round++;
            obj_before = obj_after;
            for (unsigned int i = 0; i < child.size(); i++)
                copy_of_child.at(i) = child.gene(i);
        }
        else {
            next_round = false;
            for (unsigned int i = 0; i < child.size(); i++)
                child.gene(i, copy_of_child.at(i));
        }
    }
    info = false;

    // if no improvement was made, switch to "normal" mutation
    if ((mut_round == 0) || (nMut == 0)) {
        nMut = 0;
        for (int i = 0; i < child.length(); i++) {
            if (GAFlipCoin(prob_mut)) {
                child.gene(i, child.alleleset().allele());
                nMut++;
            }
        }
        if (mut_info)
            cout << endl << endl << "Done with FlipCoin mutation, nMut=" << nMut << endl;
    }

    if (mut_info) {
        cout << endl << "Leaving Q_Mutator()..." << endl;
        cin.get();
    }
    return ((int)nMut);
}

int D_Mutator(GAGenome & c, float prob_mut) {
    GA1DArrayAlleleGenome<int> &child = (GA1DArrayAlleleGenome<int> &) c;

    bool mut_info = false;
    int nMut = 0;

    if (prob_mut <= 0.0) return (0);

    if (mut_info) {
        cout << endl << endl << "MUTATION:" << endl << endl << "Initial gene: ";
        for (int i = 0; i < n_comm; i++)
            cout << " " << child.gene(i);
        // cin.get();
    }

    for (int i = 0; i < child.length(); i++) {
        if (GAFlipCoin(prob_mut)) {
            child.gene(i, child.alleleset().allele());
            nMut++;
        }
    }

    if (mut_info) {
        cout <<  endl << "After mutation:" << flush;
        for (int i = 0; i < n_comm; i++)
            cout << " " << child.gene(i);
        cout << endl << "Leaving M_Mutator()..." << endl;
        // cin.get();
    }

    return ((int)nMut);
}


void save_state() {
    for (int i = 0; i < n_n; i++) {
        wds->cspok.at(i)->Set_dprop("concentration", best.at(i));
    }

    for (int i = 0; i < wds->agelemek.size(); i++) {
        int idx_n1 = wds->agelemek.at(i)->Get_Cspe_Index();
        int idx_n2 = wds->agelemek.at(i)->Get_Cspv_Index();
        if (idx_n1 < 0)
            idx_n1 = idx_n2;
        if (idx_n2 < 0)
            idx_n2 = idx_n1;
        double comm1 = (double) best.at(idx_n1);
        double comm2 = (double) best.at(idx_n2);
        double comm = (comm1 + comm2) / 2.;
        // cout<<endl<<"edge #"<<i<<": "<<wds->agelemek.at(i)->Get_nev()<<", community # "<<comm;
        wds->agelemek.at(i)->Set_dprop("concentration", comm);
    }

    wds->save_mod_prop_all_elements("concentration");
}

float Objective(GAGenome & c) {
    if (obj_type == "modularity") {
        Q_Objective(c);
    }
    else if ((obj_type == "D-optimality") || (obj_type == "A-optimality")) {
        D_Objective(c);
    }
    else {
        cout << endl << "ERROR: !!bad obj_type value: " << obj_type << " !!!" << endl << endl;
        exit(-1);

    }
}

float D_Objective(GAGenome & c) {
    GA1DArrayAlleleGenome<int> &genome = (GA1DArrayAlleleGenome<int> &) c;

    fcount++;
    bool obj_info = false;

    vector<int> tmp(n_comm);
    for (int i = 0; i < n_comm; i++)
        tmp.at(i) = genome.gene(i);

    bool is_same_gene = false;
    for (int i = 0; i < tmp.size(); i++) {
        for (int j = 0; j < tmp.size(); j++) {
            if ((genome.gene(i) == genome.gene(j)) && (i != j))
                is_same_gene = true;
        }
    }


    if (obj_info) {
        cout << endl << "gene:";
        for (int i = 0; i < n_comm; i++)
            cout << " " << genome.gene(i);
        cout << endl;


        cout << endl << "SM:" << endl;
        for (int i = 0; i < SM.size(); i++) {
            for (int j = 0; j < SM.at(i).size(); j++) {
                printf("%+5.3f ", SM.at(i).at(j));
            }
            cout << endl;
        }
    }

    // vector< vector<double> > SM_tr(SM.at(0).size(),SM.size());
    MatrixXd jac(SM.size(), n_comm);
    for (int i = 0; i < n_comm; i++) {
        int col = genome.gene(i);
        for (int j = 0; j < SM.size(); j++) {
            jac(j, i) = SM.at(j).at(col);
        }
    }
    if (obj_info) {
        cout << endl << "jac:" << endl;
        cout << endl << jac;

        cout << endl << "jac.transpose()*jac:" << endl;
        cout << endl << jac.transpose()*jac;
    }

    double Q;
    if (obj_type == "A-optimality") {
        double TINY_NUM = 1.e-8;
        if (!is_same_gene) {
            // MINIMIZE tr(Cov)
            MatrixXd Cov(n_comm, n_comm);
            MatrixXd Cur = jac.transpose() * jac;
            if (Cur.determinant() > TINY_NUM) {
                Cov = Cur.inverse();

                Q = 0.;
                for (unsigned int i = 0; i < n_comm; i++)
                    Q += sqrt(Cov(i, i));
                // Q = sqrt(Cov.trace());
                Q /= n_comm;
                if (obj_info) {
                    cout << endl << endl << "jac.transpose() * jac" << endl << jac.transpose() * jac;
                    cout << endl << "Cov=" << endl << Cov;
                    cout << endl << "Tr(J^T*J)=" << Cur.trace() << ", Q=" << Q << endl;
                }
            }
            else {
                Q = 1.e5;
                if (obj_info) {
                    cout << endl << "Cur=" << endl << Cur;
                    cout << endl << "det(Cur)=" << Cur.determinant() << " < " << TINY_NUM << endl;
                }
            }
            // cout << "Q=" << Q;
        }
        else {
            if (obj_info) {
                cout << endl << "Same gene found" << endl;
            }
            Q = 1.e5;
            // cout << endl << "SAME GENE!" << endl;
        }
        // cin.get();
    }
    else if (obj_type == "D-optimality") {
        // MINIMIZE det(Cov) = MAXIMIZE det(Cov^-1)
        if (!is_same_gene) {
            Q = (jac.transpose() * jac).determinant();
            Q = pow(Q, 1. / 2. / (double)n_comm);
        }
        else
            Q = 0.;
        // EigenSolver<MatrixXd> es(jac.transpose()*jac,false);

        // cout<<endl<<"Eigenvalues:"<<es.eigenvalues();

        // vector<double> eigs(es.eigenvalues().size());
        // for (int i=0; i<es.eigenvalues().size(); i++)
        //     eigs.at(i)=abs(es.eigenvalues()[i].real());

        // std::sort(eigs.begin(), eigs.begin()+eigs.size());

        // Q=1.;
        // for (unsigned int i=0; i<n_comm; i++)
        //     Q*=eigs.at(eigs.size()-1-i);

        // Q = pow(Q,1./(double)n_comm);
    }
    else {
        cout << endl << endl << "???WTF???" << endl << endl;
        exit(-1);
    }

    if (obj_info)
        cin.get();


    bool add_info = false;
    if ((obj_type == "D-optimality") && (Q > best_Q))
        add_info = true;
    if ((obj_type == "A-optimality") && (Q < best_Q))
        add_info = true;
    if ((add_info) || (info)) {
        best = tmp;
        best_Q = Q;

        //save_state();

        strstrm.str("");
        strstrm << endl << "Gene   : ";
        for (int i = 0; i < n_comm; i++)
            strstrm << best.at(i) << " ";

        if (fcount < 1000)
            strstrm << endl << "fcount : " << fcount << " (" << ((double) fcount) / f_ce * 100. << " % of compl. enum.)";
        else if (fcount < 1e6)
            strstrm << endl << "fcount : " << ((double) fcount) / 1000. << "k (" << ((double) fcount) / f_ce * 100. << " % of compl. enum.)";
        else
            strstrm << endl << "fcount : " << ((double) fcount) / 1000000. << "M (" << ((double) fcount) / f_ce * 100.
                    << " % of compl. enum.)";
        strstrm << endl << "obj    : " << Q << endl << "----------------------------------------------------";
        logfile_write(strstrm.str(), 0);
    }

    if (isnan(Q))
        Q = 0.;

    return Q;

}

float Q_Objective(GAGenome & c) {
    GA1DArrayAlleleGenome<int> &genome = (GA1DArrayAlleleGenome<int> &) c;

    fcount++;

    // Extract genome to more convenient form
    vector<int> tmp(n_n);
    vector<int> empty_comm(n_comm, 0);
    for (int i = 0; i < n_n; i++) {
        tmp.at(i) = genome.gene(i);
        for (int j = 0; j < n_comm; j++) {
            if (j == tmp.at(i))
                empty_comm.at(j) = 1;
        }
    }

    int sum_nonempty = 0;
    for (int j = 0; j < n_comm; j++)
        sum_nonempty += empty_comm.at(j);

    // Compute number of cuts
    n_c = 0;
    for (unsigned int i = 0; i < n_p; i++) {
        int idx_n1 = wds->agelemek.at(pipe_idx.at(i))->Get_Cspe_Index();
        int idx_n2 = wds->agelemek.at(pipe_idx.at(i))->Get_Cspv_Index();

        if (genome.gene(idx_n1) != genome.gene(idx_n2))
            n_c++;
    }

    // Modularity with weight W, p = abs(Apn.transpose()) * w
    Qmod = 0.;
    double Qtmp;
    for (int m = 0; m < n_comm; m++) {
        Qtmp = 0.;
        for (int i = 0; i < n_n; i++) {
            if (genome.gene(i) == m) {
                Qtmp += p(i);
            }
        }
        Qmod += (Qtmp / 2. / sumW) * (Qtmp / 2. / sumW);
    }

    double Q = 1. - ((double) n_c / (double) n_p) - Qmod;

    if ((Q > best_Q) || (info)) {
        best = tmp;
        best_Q = Q;

        save_state();

        strstrm.str("");
        strstrm << endl << "Gene:" << endl;
        for (int i = 0; i < n_n; i++)
            strstrm << best.at(i) << " ";

        if (fcount < 1000)
            strstrm << endl << "fcount  : " << fcount << " (" << ((double) fcount) / f_ce * 100.
                    << " % of compl. enum.)";
        if (fcount < 1e6)
            strstrm << endl << "fcount  : " << ((double) fcount) / 1000. << "k (" << ((double) fcount) / f_ce * 100.
                    << " % of compl. enum.)";
        else
            strstrm << endl << "fcount  : " << ((double) fcount) / 1000000. << "M (" << ((double) fcount) / f_ce * 100.
                    << " % of compl. enum.)";
        strstrm << endl << "n_c/np  : " << ((double) n_c / (double) n_p) << " (=n_c/np with n_c = " << n_c << ", np = "
                << n_p
                << ")";
        strstrm << endl << "empty c.: " << n_comm-sum_nonempty;
        strstrm << endl << "Qmod    : " << Qmod;
        strstrm << endl << "Q       : " << Q << endl << "----------------------------------------------------";
        logfile_write(strstrm.str(), 0);

    }

    return (obj_offset + Q);

}

void LoadMatrices(MatrixXd & A, VectorXd & W, VectorXd & p, string weight_type) {

    A = MatrixXd::Zero(n_n, n_n);
    MatrixXd Apn = MatrixXd::Zero(n_p, n_n);
    MatrixXd abs_Apn = MatrixXd::Zero(n_p, n_n);
    VectorXd k = VectorXd::Zero(n_n);
    p = VectorXd::Zero(n_n);
    W = VectorXd::Zero(n_p);

    // Topological incidence matrix
    for (unsigned int j = 0; j < n_p; j++) {
        int idx_n1 = wds->agelemek.at(pipe_idx.at(j))->Get_Cspe_Index();
        int idx_n2 = wds->agelemek.at(pipe_idx.at(j))->Get_Cspv_Index();
        Apn(j, idx_n1)++;
        Apn(j, idx_n2)--;
        abs_Apn(j, idx_n1)++;
        abs_Apn(j, idx_n2)++;
    }

    A = Apn.transpose() * Apn;
    for (unsigned i = 0; i < n_n; i++) {
        A(i, i) = 0;
        for (unsigned int j = 0; j < n_n; j++) {
            if (A(i, j) != 0) {
                A(i, j) = 1.;
                // k(i)++;
            }
        }
    }

    if (0 == strcmp(weight_type.c_str(), "topology")) {
        for (unsigned int i = 0; i < W.size(); i++)
            W(i) = 1;
    }

    if (0 == strcmp(weight_type.c_str(), "dp")) {
        double dp, dp_max = -1., weight_min = 1.e-5;
        for (unsigned int i = 0; i < W.size(); i++) {
            // dp = wds->agelemek.at(pipe_idx.at(i))->Get_dprop("mass_flow_rate");
            dp = wds->agelemek.at(pipe_idx.at(i))->Get_dprop("headloss");
            // dp = wds->agelemek.at(pipe_idx.at(i))->Get_dprop("length");
            W(i) = abs(dp);
            if (W(i) > dp_max)
                dp_max = W(i);
            //cout<<endl<<wds->agelemek.at(pipe_idx.at(i))->Get_nev()<<" dp="<<dp;
        }
        for (unsigned int i = 0; i < W.size(); i++) {
            W(i) /= dp_max;
            if (W(i) < weight_min)
                W(i) = weight_min;
            // cout<<endl<<wds->agelemek.at(pipe_idx.at(i))->Get_nev()<<" weight="<<W(i);
        }
    }

    if (0 == strcmp(weight_type.c_str(), "sensitivity")) {
        // Add total sensitivities
        bool is_edge_prop;

        if (weight_type_mod == "friction_coeff") {
            is_edge_prop = true;
            PerformSensitivityAnalysis(true /*is_edge_prop*/, "friction_coeff", "sensitivity_matrix_friction_coeff.csv");
        }
        else if (weight_type_mod == "diameter") {
            is_edge_prop = true;
            PerformSensitivityAnalysis(true/*is_edge_prop*/, "diameter", "sensitivity_matrix_diameter.csv");
        }
        else if (weight_type_mod == "demand") {
            is_edge_prop = false;
            PerformSensitivityAnalysis(false /*is_edge_prop*/, "demand", "sensitivity_matrix_demand.csv");
        }
        else {
            cout << endl << endl << "ERROR: illegal weight_type: >" << weight_type << "<" << endl;
            exit(-1);
        }
        vector<double> node_weights;
        int n_nodes = wds->cspok.size();
        for (int i = 0 ; i < n_nodes; i++) {
            double tmp = 0.;
            for (int j = 0 ; j < SM_PR.size(); j++) {
                tmp += SM_PR.at(j).at(i) * SM_PR.at(j).at(i);
            }
            node_weights.push_back(sqrt(tmp));
        }
        for (unsigned int i = 0; i < W.size(); i++) {
            int idx_n1 = wds->agelemek.at(i)->Get_Cspe_Index();
            int idx_n2 = wds->agelemek.at(i)->Get_Cspv_Index();
            double w1, w2;
            if (idx_n1 < 0)
                w1 = 0;
            else
                w1 = node_weights.at(idx_n1);
            if (idx_n2 < 0)
                w2 = 0;
            else
                w2 = node_weights.at(idx_n2);
            W(i) = (w1 + w2) / 2.;
            // cout << endl << wds->agelemek.at(i)->Get_nev() << ": ";
            // cout << "idx_n1:" << idx_n1 << ", w1=" << w1 << " ---> ";
            // cout << "idx_n2:" << idx_n2 << ", w2=" << w2;
        }
    }
    // cin.get();
    // Final computations
    p = abs_Apn.transpose() * W;
    sumW = W.sum();

}

void LoadSystem(string fname) {


    wds = new Staci(fname.c_str());
    wds->build_system();
    wds->ini();
    wds->solve_system();
    wds->set_res_file(wds->get_def_file());
    wds->save_results(true);

    n_n = wds->cspok.size();

    for (unsigned int i = 0; i < wds->agelemek.size(); i++) {
        bool add_this = true;
        string tipus = wds->agelemek.at(i)->GetType();
        // cout<<endl<<"i = "<<i<<", tipus : "<<tipus<<", ID : "<<wds->agelemek.at(i)->Get_nev();

        if (0 == strcmp(tipus.c_str(), "Vegakna"))
            add_this = false;
        if (0 == strcmp(tipus.c_str(), "KonstNyomas"))
            add_this = false;

        if (add_this) {
            pipe_idx.push_back(i);
            //cout<<" added.";
        }
    }
    n_p = pipe_idx.size();

}

void logfile_write(string msg, int debug_level) {
    if (global_debug_level >= debug_level) {
        ofstream logfile;
        logfile.open(logfilename.c_str(), ios::app);
        logfile << msg;
        logfile.close();
        cout << msg;
    }
}

void Load_Settings() {

    XMLNode xMainNode = XMLNode::openFileHelper("staci_split_settings.xml", "settings");
    global_debug_level = atoi(xMainNode.getChildNode("global_debug_level").getText());
    n_comm = atoi(xMainNode.getChildNode("n_comm").getText());
    weight_type = xMainNode.getChildNode("weight_type").getText();
    weight_type_mod = xMainNode.getChildNode("weight_type_mod").getText();
    fname = xMainNode.getChildNode("fname").getText();
    logfilename = xMainNode.getChildNode("logfilename").getText();
    obj_type = xMainNode.getChildNode("obj_type").getText();

    popsize = atoi(xMainNode.getChildNode("popsize").getText());
    ngen = atoi(xMainNode.getChildNode("ngen").getText());
    pmut = atof(xMainNode.getChildNode("pmut").getText());
    pcross = atof(xMainNode.getChildNode("pcross").getText());

    stringstream msg;
    msg.str("");
    msg << "\n================================ ";
    msg << "\n        STACI SPLIT";
    msg << "\n================================\n";
    msg << "Settings:" << endl << endl;
    msg << "\t global_debug_level     : " << global_debug_level << endl;
    msg << "\t n_comm                 : " << n_comm << endl;
    msg << "\t obj_type               : " << obj_type << endl;
    msg << "\t weight_type            : " << weight_type << endl;
    msg << "\t weight_type_mod        : " << weight_type_mod << " (used to allocate nodes in communities)" << endl;

    bool is_weight_ok = false;
    if (obj_type == "modularity") {
        if ( (weight_type == "topology") || (weight_type == "dp") || (weight_type == "sensitivity")) {
            // msg << " (ok)" << endl;
            is_weight_ok = true;
        }
    }

    if ((obj_type == "A-optimality") || (obj_type == "D-optimality")) {
        if (  (weight_type == "friction_coeff") || (weight_type == "demand") || (weight_type == "diameter")) {
            // msg << " (ok)" << endl;
            is_weight_ok = true;
        }
    }

    if (!is_weight_ok) {
        if (obj_type == "modularity") {
            cout << endl << "ERROR in Load_Settings() unknown weight_type >" << weight_type << "<";
            cout << endl << "obj_type = " << obj_type << " -> possible weigths: topology, dp";
            cout << endl << "Exiting..." << endl;
            exit(-1);
        }
        else if ((obj_type == "A-optimality") || (obj_type == "D-optimality")) {
            cout << endl << "ERROR in Load_Settings() unknown weight_type >" << weight_type << "<";
            cout << endl << "obj_type = " << obj_type << " -> possible weigths: friction_coeff, diameter, demand";
            cout << "Exiting..." << endl;
            exit(-1);
        }
        else {
            cout << endl << "ERROR in Load_Settings() unknown obj_type >" << obj_type << "<";
            cout << endl << "obj_type = " << obj_type << " -> possible weigths: modularity, A-optimality, D-optimality";
            cout << endl << "Exiting..." << endl;
            exit(-1);
        }
    }





    msg << "\t fname                  : " << fname << endl;
    msg << "\t logfilename            : " << logfilename << endl << endl;
    msg << "\t popsize : " << popsize << endl;
    msg << "\t ngen    : " << ngen << endl;
    msg << "\t pmut    : " << pmut << endl;
    msg << "\t pcross  : " << pcross << endl << endl;

    logfile_write(msg.str(), 0);
}

void print(igraph_t *g) {
    igraph_vector_t el;
    long int i, j, n;
    char ch = igraph_is_directed(g) ? '>' : '-';

    igraph_vector_init(&el, 0);
    igraph_get_edgelist(g, &el, 0);
    n = igraph_ecount(g);

    for (i = 0, j = 0; i < n; i++, j += 2) {
        printf("%ld --%c %ld: %ld\n",
               (long)VECTOR(el)[j], ch, (long)VECTOR(el)[j + 1], (long)EAN(g, "weight", i));
    }
    printf("\n");

    igraph_vector_destroy(&el);
}

int print_vector(const igraph_vector_t *v) {
    long int i, n = igraph_vector_size(v);
    for (i = 0; i < n; i++) {
        printf("%.2g", (double)VECTOR(*v)[i]);
        if (i != n - 1) { printf(" "); }
    }
    printf("\n");
    return 0;
}

int print_matrix(const igraph_matrix_t *m) {
    long int i, j, nrow = igraph_matrix_nrow(m), ncol = igraph_matrix_ncol(m);
    for (i = 0; i < nrow; i++) {
        for (j = 0; j < ncol; j++) {
            printf("%.2g", (double)MATRIX(*m, i, j));
            if (j != ncol - 1) { printf(" "); }
        }
        printf("\n");
    }
    return 0;
}

int igraph_community_eigenvector() {

    igraph_t g;
    igraph_matrix_t mat;
    long int i, j;

    igraph_matrix_init(&mat, n_n, n_n);
    for (i = 0; i < n_n; i++) {
        for (j = 0; j < n_n; j++) {
            MATRIX(mat, i, j) = (double) A(i, j);
        }
    }

    igraph_i_set_attribute_table(&igraph_cattribute_table);

    // igraph_weighted_adjacency(&g, &mat, IGRAPH_ADJ_UNDIRECTED,0,1);
    igraph_adjacency(&g, &mat, IGRAPH_ADJ_UNDIRECTED);

    //igraph_t g;
    igraph_matrix_t merges;
    igraph_vector_t membership;
    igraph_vector_t x;
    igraph_arpack_options_t options;

    igraph_matrix_init(&merges, 0, 0);
    igraph_vector_init(&membership, 0);
    igraph_vector_init(&x, 0);
    igraph_arpack_options_init(&options);

// weights
    igraph_vector_t modularity, weights;
    igraph_vector_init(&modularity, 0);
    igraph_vector_init(&weights, 0);
    igraph_vector_resize(&weights, igraph_ecount(&g));

    // cout<<endl<<"igraph_ecount(&g)="<<igraph_ecount(&g)<<", W="<<W.size();
    // cin.get();
    // long int n=igraph_vector_size(weights);
    for (i = 0; i < n_p; i++)
        VECTOR(weights)[i] = W(i);

//============== EIGENVECTORS =========================

    igraph_community_leading_eigenvector(&g, &weights, &merges,
                                         &membership, n_c + 1,
                                         &options, /*modularity=*/ 0,
                                         /*start=*/ 0, /*eigenvalues=*/ 0,
                                         /*eigenvectors=*/ 0, /*history=*/ 0,
                                         /*callback=*/ 0,
                                         /*callback_extra=*/ 0);

// ============== GREEDY =========================

    // igraph_community_fastgreedy(&g, &weights, &merges, &modularity, &membership);

// ============== MULTILEVEL =========================

    // igraph_matrix_t memberships;
    // igraph_matrix_init(&memberships,0,0);
    // igraph_community_multilevel(&g, &weights, &membership, &memberships, &modularity);

// ============== OPTIMAL =========================
// igraph_real_t modularity2;
//     igraph_community_optimal_modularity(&g, &modularity2, &membership, &weights);

// ============== RANDOM WALK =========================

// igraph_community_walktrap(&g, &weights, 50, &merges, &modularity,  &membership);

// ============== EDGE BETWEENNESS =========================
// igraph_vector_t edges, eb;
//  igraph_vector_init(&edges, 0);
//  igraph_vector_init(&eb, 0);
// igraph_community_edge_betweenness(&g, &edges, &eb, 0 /*merges */,
//                    0 /*bridges */, /*modularity=*/ 0,
//                    &membership,
//                    IGRAPH_UNDIRECTED,
//                    &weights);

//============== END =========================

    // print_matrix(&merges);
    // print_vector(&membership);
    // save_state(&membership);
    copy_to_best(&membership);

    igraph_vector_destroy(&x);
    igraph_vector_destroy(&membership);
    igraph_matrix_destroy(&merges);
    igraph_destroy(&g);

    return 0;
}

void copy_to_best(const igraph_vector_t *v) {
    long int length = igraph_vector_size(v);
    best.clear();
    for (int i = 0; i < length; i++)
        best.push_back( (int) VECTOR(*v)[i]);
}

void save_state(const igraph_vector_t *v) {
    for (int i = 0; i < n_n; i++) {
        wds->cspok.at(i)->Set_dprop("concentration", (double)VECTOR(*v)[i]);
    }

    for (int i = 0; i < wds->agelemek.size(); i++) {
        int idx_n1 = wds->agelemek.at(i)->Get_Cspe_Index();
        int idx_n2 = wds->agelemek.at(i)->Get_Cspv_Index();
        if (idx_n1 < 0)
            idx_n1 = idx_n2;
        if (idx_n2 < 0)
            idx_n2 = idx_n1;
        double comm1 = (double)VECTOR(*v)[idx_n1];
        double comm2 = (double)VECTOR(*v)[idx_n2];
        double comm = (comm1 + comm2) / 2.;
        // cout<<endl<<"edge #"<<i<<": "<<wds->agelemek.at(i)->Get_nev()<<", community # "<<comm;
        wds->agelemek.at(i)->Set_dprop("concentration", comm);
    }

    wds->save_mod_prop_all_elements("concentration");
}
