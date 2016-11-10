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

void LoadMatrices(MatrixXd &A, VectorXd &W, VectorXd &p, string weight_type);

void LoadSystem(string fname);

void save_state();
void save_state(const igraph_vector_t *v);

void logfile_write(string msg, int debug_level);

// int QMutator(GAGenome & c, float prob_mut);
int QMutator(GAGenome & c, float prob_mut);
void igraph_Initializer(GAGenome & c);

void Optimize();

void Load_Settings();

int global_debug_level;

vector<int> best;
double best_Q;
double obj_offset;

stringstream strstrm;
string logfilename;

void print(igraph_t *g);
int print_vector(const igraph_vector_t *v);
int print_matrix(const igraph_matrix_t *m);
int igraph_community_eigenvector();
void copy_to_best(const igraph_vector_t *v);

void PerformSensitivityAnalysis();

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
    PerformSensitivityAnalysis();

    // Number of funeval for complete enumeration
    f_ce = pow((double) n_comm, (double) n_n);

    // Load matrices
    LoadMatrices(A, W, p, weight_type);

    // let it go...
    Optimize();

    return 0;
}

void PerformSensitivityAnalysis() {
    wds->Compute_Sensitivity_Matrix("friction_coeff", 1);

    vector< vector<double> > SM_MFR = wds->SM_MassFlowRates;
    vector< vector<double> > SM_PR  = wds->SM_Pressures;

    FILE * pFile;

    int n_edges = wds->agelemek.size();
    int n_nodes = wds->cspok.size();

    pFile = fopen ("sensitivity_matrix.txt", "w");
    fprintf(pFile, "parameter: friction_coeff\n");
    fprintf(pFile, "n_edges : %d\n", n_edges);
    fprintf(pFile, "n_nodes : %d\n", n_nodes);
    fprintf(pFile, "n_cols  : %d ( = n_edges + n_nodes )\n", n_edges+n_nodes);
    fprintf(pFile, "n_rows  : %d ( = n_edges )\n", n_edges);
    for (int i = 0 ; i < n_edges; i++) {
        fprintf(pFile,"friction_coeff @ %s, ",wds->agelemek.at(i)->Get_nev().c_str());
        for (int j = 0 ; j < n_edges; j++) 
            fprintf (pFile, "%+7.5e, ", SM_MFR.at(i).at(j));
        for (int j = 0 ; j < n_nodes; j++) 
            fprintf (pFile, "%+7.5e, ", SM_PR.at(i).at(j));
        fprintf (pFile, "\n");
    }
    fclose (pFile);
}

void Optimize() {

    cout << endl << "RUNNING OPTMIZATION..." << endl;

    int *aset = new int[n_comm];
    for (int i = 0; i < n_comm; i++)
        aset[i] = i;
    GAAlleleSet<int> alleles(n_comm, aset);

    GA1DArrayAlleleGenome<int> genome(n_n, alleles, Objective);
    genome.mutator(QMutator);
    genome.initializer(igraph_Initializer);

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

    genome = ga.statistics().bestIndividual();
    //double best_obj = Objective(genome);

    strstrm.str("");
    strstrm << endl << "BEST SOLUTION FOUND:" << endl;

    for (unsigned int i = 0; i < n_n; i++)
        strstrm << genome.gene(i) << " ";
    strstrm << endl;

    info = true;
    double best_obj = Objective(genome);

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

}

void igraph_Initializer(GAGenome & c) {
    GA1DArrayAlleleGenome<int> &child = (GA1DArrayAlleleGenome<int> &) c;

    if (n_comm == 2) {
        igraph_community_eigenvector();
        for (int i = 0; i < n_n; i++)
            child.gene(i, best.at(i));
        strstrm.str();
        strstrm << endl << endl << "igraph_Initializer with c_comm=" << n_comm << "):" << endl;
        for (int i = 0; i < child.length(); i++)
            strstrm << child.gene(i) << " ";
        logfile_write(strstrm.str(), 0);
    }
    else {
        int wmax = n_comm - 1;
        int wmin = 0;
        for (int i = 0; i < n_n; i++)
            child.gene(i, wmin + (rand() % (int)(wmax - wmin + 1)));

        strstrm.str();
        strstrm << endl << endl << "Random initializer with c_comm=" << n_comm << "):" << endl;
        for (int i = 0; i < n_n; i++)
            strstrm << child.gene(i) << " ";
        logfile_write(strstrm.str(), 0);
    }
}

// int QMutator(GAGenome & c, float prob_mut){
//     GA1DArrayAlleleGenome<int> &child = (GA1DArrayAlleleGenome<int> &) c;

//     // Copy of child
//     vector<int> copy_of_child(child.size());
//     for (unsigned int i=0; i<child.size(); i++)
//         copy_of_child.at(i)=child.gene(i);

//     bool mut_info = false;
//     info = mut_info;
//     int mut_round=0;
//     bool next_round=true;
//     double obj_before=0., obj_after=1.;
//     int nMut = 0;

//     if(prob_mut <= 0.0) return(0);

//     if (mut_info){
//         cout<<endl<<endl<<"MUTATION:"<<endl<<endl<<"Initial gene:"<<flush;
//         // cin.get();
//     }

//     obj_before = Objective(child);
//     while (next_round && (mut_round<10)){
//         if (mut_info)
//             cout<<endl<<endl<<"ROUND "<<mut_round;

//         int n_c = 0;
//         vector<bool> cut_idx(n_n,false);
//         vector<int> cut_pair(n_n,0);
//         for (unsigned int i = 0; i < n_p; i++) {
//             int idx_n1 = wds->agelemek.at(pipe_idx.at(i))->Get_Cspe_Index();
//             int idx_n2 = wds->agelemek.at(pipe_idx.at(i))->Get_Cspv_Index();

//             if (child.gene(idx_n1) != child.gene(idx_n2)){
//                 n_c++;
//                 cut_idx.at(idx_n1)=true;
//                 cut_idx.at(idx_n2)=true;
//                 cut_pair.at(idx_n1)=idx_n2;
//                 cut_pair.at(idx_n2)=idx_n1;
//                 if (mut_info){
//                     cout<<endl<<"\tCut: "<<idx_n1 <<"("<<child.gene(idx_n1);
//                     cout<<") --x-- "<<idx_n2<<"("<<child.gene(idx_n2)<<")";
//                 }
//             }
//         }
//         if (mut_info)
//             cout<<endl;

//         for(unsigned int i=0; i<n_n; i++){
//             if (cut_idx.at(i)){
//                 if (GAFlipCoin(0.5)){
//                     int other_node_community = child.gene(cut_pair.at(i));
//                     if (mut_info){
//                         cout<<endl<<" --> node ="<<i<<" moved to comm."<<other_node_community;
//                     }
//                     cut_idx.at(i)=false;
//                     cut_idx.at(cut_pair.at(i))=false;

//                     child.gene(i,other_node_community);
//                     nMut++;
//                 }
//             }
//         }
//         if (mut_info)
//             cout<<endl;

//         obj_after = Objective(child);

//         if (obj_after>obj_before){
//             next_round=true;
//             mut_round++;
//             obj_before=obj_after;
//             for (unsigned int i=0; i<child.size(); i++)
//                 copy_of_child.at(i)=child.gene(i);
//         }
//         else{
//             next_round=false;
//             for (unsigned int i=0; i<child.size(); i++)
//                 child.gene(i,copy_of_child.at(i));
//         }
//     }
//     info = false;

//     // if no improvement was made, switch to "normal" mutation
//     if ((mut_round==0)||(nMut==0)) {
//         nMut = 0;
//         for (int i=0; i<child.length(); i++){
//             if (GAFlipCoin(prob_mut)){
//                 child.gene(i, child.alleleset().allele());
//                 nMut++;
//             }
//         }
//         if (mut_info)
//             cout<<endl<<endl<<"Done with FlipCoin mutation, nMut="<<nMut<<endl;
//     }

//     if (mut_info){
//         cout<<endl<<"Leaving QMutator()..."<<endl;
//         cin.get();
//     }
//     return((int)nMut);
// }

int QMutator(GAGenome & c, float prob_mut) {
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
        cout << endl << "Leaving QMutator()..." << endl;
        cin.get();
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

float Objective(GAGenome &c) {
    GA1DArrayAlleleGenome<int> &genome = (GA1DArrayAlleleGenome<int> &) c;

    fcount++;

    // Extract genome to more convenient form
    vector<int> tmp(n_n);
    //tmp.at(0) = 0;
    for (int i = 0; i < n_n; i++)
        tmp.at(i) = genome.gene(i);

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

void LoadMatrices(MatrixXd &A, VectorXd &W, VectorXd &p, string weight_type) {

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
        double dp, dp_max = -1., weight_min = 0.0001;
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
    fname = xMainNode.getChildNode("fname").getText();
    logfilename = xMainNode.getChildNode("logfilename").getText();

    popsize = atoi(xMainNode.getChildNode("popsize").getText());
    ngen = atoi(xMainNode.getChildNode("ngen").getText());
    pmut = atof(xMainNode.getChildNode("pmut").getText());
    pcross = atof(xMainNode.getChildNode("pcross").getText());

    stringstream msg;
    msg.str("");
    msg << "Settings:" << endl << endl;
    msg << "\t global_debug_level     : " << global_debug_level << endl;
    msg << "\t n_comm                 : " << n_comm << endl;
    msg << "\t weight_type            : " << weight_type;

    if ( (weight_type == "topology") || (weight_type == "dp") ) {
        msg << " (ok)" << endl;
    }
    else {
        cout << endl << "!!!! UNKNOWN WEIGHT_TYPE!!!" << endl << "Exiting..." << endl;
        exit(-1);
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
    for (int i = 0; i < n_n; i++)
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
