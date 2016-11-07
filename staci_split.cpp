#include <stdio.h>
#include <ga/ga.h>
#include <vector>
#include "Staci.h"
#include <iomanip>
#include <boost/tokenizer.hpp>
#include <numeric>
#include <list>
#include "xmlParser.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <ga/GASimpleGA.h>
#include <ga/GAListGenome.h>

float Objective(GAGenome &);

float Objective2(GAGenome &);

using namespace std;
using namespace boost;
using namespace Eigen;

#define cout STD_COUT

// These global variables are coming from the datafile
string fname;
int popsize, ngen;
double pcross, pmut;

Staci *wds;
int n_n, n_p, n_comm;
int fcount;
double f_ce;
vector<int> pipe_idx;
VectorXd p;
double sumW;

int n_c;
double Qmod;

bool info;

void LoadMatrices(VectorXd &p, string weight_type);

void LoadSystem(string fname);

void save_state();

void logfile_write(string msg, int debug_level);

string Load_Settings();

int global_debug_level;

vector<int> best;
double best_Q;
double obj_offset;

stringstream strstrm;
string logfilename;

int main(int argc, char **argv) {

    best_Q = -100.;
    obj_offset = 1.0;
    info = false;
    fcount = 0;

    string msg = Load_Settings();
    logfile_write(msg,0);

    // Clear logfile
    ofstream ofs(logfilename.c_str(), std::ios::out | std::ios::trunc);
    ofs.close();

    // Load system
    LoadSystem(fname);
    f_ce = pow((double) n_comm, (double) n_n);

    // Load matrices
    LoadMatrices(p, "uniform");

    // Optimize
    cout << endl << "RUNNING OPTMIZATION...";
    cout << flush;

    int *aset = new int[n_comm];
    for (int i = 0; i < n_comm; i++)
        aset[i] = i;
    GAAlleleSet<int> alleles(n_comm, aset);

    GA1DArrayAlleleGenome<int> genome(n_n - 1, alleles, Objective);

    GASteadyStateGA ga(genome);

    ga.set(gaNpopulationSize, popsize);    // population size
    ga.set(gaNpCrossover, pcross);        // probability of crossover
    ga.set(gaNpMutation, pmut);        // probability of mutation
    ga.set(gaNnGenerations, ngen);        // number of generations
    ga.set(gaNscoreFrequency, 100);        // how often to record scores
    ga.set(gaNflushFrequency, 100);    // how often to dump scores to file
    ga.set(gaNselectScores,        // which scores should we track?
           GAStatistics::Maximum | GAStatistics::Minimum | GAStatistics::Mean);
    ga.set(gaNscoreFilename, "bog.dat");

    ga.evolve();

    genome = ga.statistics().bestIndividual();
    double best_obj = Objective(genome);

    cout << endl << "done." << endl;
    strstrm.str("");
    strstrm << endl << "BEST SOLUTION FOUND:" << endl;

    strstrm << "0 ";
    for (unsigned int i = 0; i < n_n - 1; i++)
        strstrm << genome.gene(i) << " ";
    strstrm << endl;

    best_obj = Objective(genome);

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

// Write Settings to logfile
    logfile_write(strstrm.str(), 0);

    return 0;
}

void save_state() {
    for (int i = 0; i < n_n; i++){
        wds->cspok.at(i)->Set_dprop("concentration", best.at(i));
        // cout<<endl<<"node #"<<i<<": "<<wds->cspok.at(i)->Get_nev()<<", comm. #"<<best.at(i);
    }

    for (int i = 0; i < wds->agelemek.size(); i++) {
        int idx_n1 = wds->agelemek.at(i)->Get_Cspe_Index();
        int idx_n2 = wds->agelemek.at(i)->Get_Cspv_Index();
        if (idx_n1<0)
            idx_n1=idx_n2;
        if (idx_n2<0)
            idx_n2=idx_n1;
        double comm1 = (double) best.at(idx_n1);
        double comm2 = (double) best.at(idx_n2);
        double comm = (comm1+comm2)/2.;
        // cout<<endl<<"edge #"<<i<<": "<<wds->agelemek.at(i)->Get_nev()<<", community # "<<comm;
        wds->agelemek.at(i)->Set_dprop("concentration", comm);
    }

    wds->save_mod_prop_all_elements("concentration");
}

float
Objective(GAGenome &c) {
    GA1DArrayAlleleGenome<int> &genome = (GA1DArrayAlleleGenome<int> &) c;

    fcount++;

    // Extract genome to more convenient form
    vector<int> tmp(n_n);
    tmp.at(0) = 0;
    for (int i = 0; i < n_n - 1; i++)
        tmp.at(i + 1) = genome.gene(i);

    // Compute number of cuts
    n_c = 0;
    for (unsigned int i = 0; i < n_p; i++) {
        int idx_n1 = wds->agelemek.at(pipe_idx.at(i))->Get_Cspe_Index();
        int idx_n2 = wds->agelemek.at(pipe_idx.at(i))->Get_Cspv_Index();

        if (tmp.at(idx_n1) != tmp.at(idx_n2)) {
            n_c++;
            /*if (info) {
                cout << endl << "Name: " << wds->agelemek.at(i)->Get_nev();
                cout << " node from: #" << idx_n1 << " comm. #" << tmp.at(idx_n1);
                cout << " -------> node to: #" << idx_n2 << " comm. #" << tmp.at(idx_n2);
                cout << " !! PIPE IS CUT !! -> n_c=" << n_c;
            }*/
        }
    }

    // Modularity with weight W, p = abs(Apn.transpose()) * w
    Qmod = 0.;
    double Qtmp;
    for (unsigned int m = 0; m < n_comm; m++) {
        Qtmp = 0.;
        for (unsigned int i = 0; i < n_n; i++) {
            if (tmp.at(i) == m) {
                Qtmp += p(i);
            }
        }
        Qmod += (Qtmp / 2 / sumW) * (Qtmp / 2 / sumW);
    }

    double Q = 1 - ((double) n_c / (double) n_p) - Qmod;

    if (Q > best_Q) {
        best = tmp;
        best_Q = Q;

        save_state();

        strstrm.str("");
        if (fcount < 1000)
            strstrm << endl << "fcount : " << fcount << " (" << ((double) fcount) / f_ce * 100.
                    << " % of compl. enum.)";
        if (fcount < 1e6)
            strstrm << endl << "fcount : " << ((double) fcount) / 1000. << "k (" << ((double) fcount) / f_ce * 100.
                    << " % of compl. enum.)";
        else
            strstrm << endl << "fcount : " << ((double) fcount) / 1000000. << "M (" << ((double) fcount) / f_ce * 100.
                    << " % of compl. enum.)";
        strstrm << endl << "n_c/np : " << ((double) n_c / (double) n_p) << " (=n_c/np with n_c = " << n_c << ", np = "
                << n_p
                << ")";
        strstrm << endl << "Qmod   : " << Qmod;
        strstrm << endl << "Q      : " << Q << endl << "----------------------------------------------------";
        logfile_write(strstrm.str(), 0);
        
    }

    return (obj_offset + Q);

}

void LoadMatrices(VectorXd &p, string weight_type) {

    MatrixXd A = MatrixXd::Zero(n_n, n_n);
    MatrixXd Apn = MatrixXd::Zero(n_p, n_n);
    MatrixXd abs_Apn = MatrixXd::Zero(n_p, n_n);
    VectorXd k = VectorXd::Zero(n_n);
    p = VectorXd::Zero(n_n);
    VectorXd W = VectorXd::Zero(n_p);

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
                k(i)++;
            }
        }
    }

    // Now set the weight
    // Node rank = topology
    if (0 == strcmp(weight_type.c_str(), "uniform")) {
        for (unsigned int i = 0; i < W.size(); i++)
            W(i) = 1;
    }

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
        string tipus = wds->agelemek.at(i)->Get_Tipus();
        if (0 == strcmp(tipus.c_str(), "Cso"))
            pipe_idx.push_back(i);
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

string Load_Settings() {

    XMLNode xMainNode = XMLNode::openFileHelper("staci_split_settings.xml", "settings");
    global_debug_level = atoi(xMainNode.getChildNode("global_debug_level").getText());
    n_comm = atoi(xMainNode.getChildNode("n_comm").getText());
    fname = xMainNode.getChildNode("fname").getText();
    logfilename = xMainNode.getChildNode("logfilename").getText();

    //string fname = "jwrpm_test2.spr";
    //string fname = "jwrpm_test_wiki.spr";
    //string fname = "VIZ-SOPTVR-T-81-input_mod.spr";
    //string fname = "BUK-AMEVAV-1-80-input_mod.spr";
    //string fname = "LOV-SOPRKV-1-input_mod.spr";
    //string fname = "VIZ-SOPTVR-E-input_mod.spr";

    popsize = atoi(xMainNode.getChildNode("popsize").getText());
    ngen = atoi(xMainNode.getChildNode("ngen").getText());
    pmut = atof(xMainNode.getChildNode("pmut").getText());
    pcross = atof(xMainNode.getChildNode("pcross").getText());

    stringstream msg;
    msg.str("");
    msg << "Settings:" << endl << endl;
    msg << "\t global_debug_level     : " << global_debug_level << endl;
    msg << "\t n_comm                 : " << n_comm << endl << endl;
    msg << "\t fname                  : " << fname << endl;
    msg << "\t logfilename            : " << logfilename << endl;
    msg << "\t popsize : " << popsize << endl;
    msg << "\t ngen    : " << ngen << endl;
    msg << "\t pmut    : " << pmut << endl;
    msg << "\t pcross  : " << pcross << endl << endl;

    return msg.str();
}


