#include <stdio.h>
#include <ga/ga.h>
#include <vector>
#include "Staci.h"
#include <iomanip>
#include <boost/tokenizer.hpp>
#include <numeric>
#include <list>
#include "xmlParser.h"

using namespace std;
using namespace boost;

#define cout STD_COUT

int global_debug_level;
int Staci_debug_level;

string dir_name;
string fname_prefix;
string logfilename;
string best_logfilename;
string sollwert_dfile;

void Load_WDS_Systems();

void Load_Sollwert_Datafiles();

void Update_Reservoirs(unsigned int);

void logfile_write(string msg, int debug_level);

int Num_of_Periods;
int Start_of_Periods;
vector<Staci *> wds;

vector<string> Sollwert_Node_Staci_ID;
vector<int> Sollwert_Node_Staci_Idx;
vector<vector<double> > Sollwert_Node_Values;
vector<vector<double> > Staci_Node_Values;

vector<string> Sollwert_Pool_Staci_ID;
vector<int> Sollwert_Pool_Staci_Idx;
vector<vector<double> > Sollwert_Pool_Values;
vector<vector<double> > Staci_Pool_Values;

vector<double> Pool_Surfs;

vector<string> pipe_name;
vector<bool> pipe_is_active;
unsigned int num_of_active_pipes;
vector<double> pipe_origD;

void Set_Up_Active_Pipes();

void Spoil_Active_Pipes();

bool bool_Spoil_Active_Pipes;
double Dmin;

double dt;

int Find_Pressure_Index(string PressureName);

int Find_Pool_Index(string PressureName);

void Set_Initial_Pool_Levels();

double get_A(string PoolName);

double Compute_Error();

float Objective(GAGenome &);

bool do_PerformSensitivityAnalysis;

void PerformSensitivityAnalysis();

vector<vector<double> > pipe_Sensitivity_MassFlowRates;
vector<vector<double> > pipe_Sensitivity_Pressures;

std::vector<std::string> csv_read_row(std::istream &in, char delimiter);

std::vector<std::string> csv_read_row(std::string &in, char delimiter);

int Obj_Eval;

string Load_Settings();

void PrintBestDataFile(GAGenome &);

bool last_computation_OK;
int popsize;
int ngen;
float pmut;
float pcross;
float best_obj;

int main(int argc, char **argv) {

    Obj_Eval = 0;
    best_obj = 1.e10;
    last_computation_OK = false;
    do_PerformSensitivityAnalysis = false;

    // Load parameters, including logfile name
    string msg_set = Load_Settings();

    // Clear logfile and best_individual.txt
    ofstream ofs(logfilename.c_str(), std::ios::out | std::ios::trunc);
    ofs.close();
    ofstream ofs1(best_logfilename.c_str(), std::ios::out | std::ios::trunc);
    ofs1.close();

    // Write Settings to logfile
    logfile_write(msg_set, 0);

    // Load WDSs and data files
    Load_WDS_Systems();
    Load_Sollwert_Datafiles();
    Set_Initial_Pool_Levels();

    // Initial computation
    for (unsigned i = 0; i < (wds.size()); i++) {
        wds.at(i)->Set_debug_level(Staci_debug_level);
        bool success = wds.at(i)->solve_system();
        if (success) {
            wds.at(i)->Compute_Sensitivity_Matrix("diameter", 1);

            /*for (unsigned jj=0; jj<wds.at(i)->SM_row_name.size(); jj++){
                cout<<endl<<" SM row "<<jj<<": "<<wds.at(i)->SM_row_name.at(jj);
                cout<<"-> m:"<<wds.at(i)->SM_row_sum_MassFlowRates.at(jj);
                cout<<", p:"<<wds.at(i)->SM_row_sum_Pressures.at(jj);
                cin.get();
            }*/
            if (i > Start_of_Periods)
                Update_Reservoirs(i);
        } else {
            stringstream msg;
            msg.str("");
            msg << "\n\n ERROR: wds.at(" << i << ") could not be soled!";
            logfile_write(msg.str(), 0);
            exit(-1);
        }
    }

    // Get difference
    logfile_write("Evaluating initial state.\n", 0);
    double err_pres = Compute_Error();

    // Optimization

    Set_Up_Active_Pipes();
    if (bool_Spoil_Active_Pipes)
        Spoil_Active_Pipes();

    // See if we've been given a seed to use (for testing purposes).  When you
    // specify a random seed, the evolution will be exactly the same each time
    // you use that seed number.

    unsigned int seed = 0;
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i++], "seed") == 0) {
            seed = atoi(argv[i]);
        }
    }

    // Declare variables for the GA parameters and set them to some default values.



    // Create a phenotype for two variables.  The number of bits you can use to
    // represent any number is limited by the type of computer you are using.  In
    // this case, we use 16 bits to represent a floating point number whose value
    // can range from -5 to 5, inclusive.  The bounds on x1 and x2 can be applied
    // here and/or in the objective function.
    GABin2DecPhenotype map;
    double LOWER_BOUND, UPPER_BOUND;
    stringstream msg;
    logfile_write("Setting up LOWER_BOUND, UPPER_BOUND\n", 0);
    msg.str("");
    msg << "\n\n:";
    double dia;
    for (unsigned int i = 0; i < pipe_name.size(); i++) {
        if (pipe_is_active.at(i)) {
            dia = wds.at(0)->get_dprop(pipe_name.at(i), "diameter");
            LOWER_BOUND = 0.5 * dia;
            UPPER_BOUND = 2.0 * dia;
            map.add(16, LOWER_BOUND, UPPER_BOUND);
            msg << endl << "\t" << pipe_name.at(i) << " D=" << dia << ", LB=" << LOWER_BOUND << ", UB=" << UPPER_BOUND;
        }
    }
    msg << endl;
    logfile_write(msg.str(), 1);

    // Create the template genome using the phenotype map we just made.

    GABin2DecGenome genome(map, Objective);

    // Now create the GA using the genome and run it.  We'll use sigma truncation
    // scaling so that we can handle negative objective scores.

    logfile_write("Starting optimization...\n", 0);

    GASimpleGA ga(genome);
    GASigmaTruncationScaling scaling;
    ga.minimize();
    ga.populationSize(popsize);
    ga.nGenerations(ngen);
    ga.pMutation(pmut);
    ga.pCrossover(pcross);
    ga.scaling(scaling);

    ga.scoreFrequency(10);
    ga.flushFrequency(10);
    ga.scoreFilename("bog.dat");
    ga.evolve(seed);

    do_PerformSensitivityAnalysis = true;
    ga.statistics().write("bog_stats.dat");
    genome = ga.statistics().bestIndividual();
    Objective(genome);
    PrintBestDataFile(genome);

    do_PerformSensitivityAnalysis = false;

    logfile_write("Optimization finished.\n",0);
    msg.str("");
    msg << "Best solution objective is " << Objective(genome) << endl << endl;
    logfile_write(msg.str(),0);

    return 0;
}

void PrintBestDataFile(GAGenome &x) {
    GABin2DecGenome &genome = (GABin2DecGenome &) x;

    ofstream resfile;
    resfile.open("best.dat", ios::out);
    resfile << scientific << setprecision(5);

    for (unsigned int j = 0; j < genome.nPhenotypes(); j++)
        if (pipe_is_active.at(j))
            resfile << pipe_name.at(j) << " ; ";
    resfile << endl;
    for (unsigned int j = 0; j < genome.nPhenotypes(); j++)
        resfile << genome.phenotype(j) << " ; ";
    resfile << endl;
    for (unsigned int j = 0; j < genome.nPhenotypes(); j++)
        resfile << genome.phenotype(j) / pipe_origD.at(j) << " ; ";
    resfile << endl;

    for (unsigned int i = 0; i < Sollwert_Pool_Staci_ID.size(); i++)
        resfile << "Sollwert " << Sollwert_Pool_Staci_ID.at(i) << " ; " << Sollwert_Pool_Staci_ID.at(i) << " ; ";

    for (unsigned int i = 0; i < Sollwert_Node_Staci_ID.size(); i++)
        resfile << "Sollwert " << Sollwert_Node_Staci_ID.at(i) << " ; " << Sollwert_Node_Staci_ID.at(i) << " ; ";

    for (unsigned int j = 0; j < Num_of_Periods; j++) {
        resfile << endl;
        for (unsigned int i = 0; i < Sollwert_Pool_Staci_ID.size(); i++)
            resfile << Sollwert_Pool_Values.at(i).at(j) << " ; " << Staci_Pool_Values.at(i).at(j) << " ; ";
        for (unsigned int i = 0; i < Sollwert_Node_Staci_ID.size(); i++)
            resfile << Sollwert_Node_Values.at(i).at(j) << " ; " << Staci_Node_Values.at(i).at(j) << " ; ";
    }
    resfile.close();
}

float
Objective(GAGenome &x) {
    GABin2DecGenome &genome = (GABin2DecGenome &) x;

    bool success = false;

    for (unsigned int i = 0; i < wds.size(); i++) {
        int k = 0;
        if (i > 0) {
            wds.at(i)->ini(wds.at(i - 1));
        } else
            wds.at(i)->ini();
        for (unsigned int j = 0; j < pipe_name.size(); j++)
            if (pipe_is_active.at(j)) {
                wds.at(i)->set_dprop(pipe_name.at(j), "diameter", genome.phenotype(k));
                k++;
            }
        last_computation_OK = wds.at(i)->solve_system();
        if (last_computation_OK) {
            if (i > 0)
                Update_Reservoirs(i);
        } else {
            break;
        }
    }

    if (do_PerformSensitivityAnalysis)
        PerformSensitivityAnalysis();

    double err;
    if (last_computation_OK) {
        err = Compute_Error();
        if ((Obj_Eval % 100) == 0) {
            stringstream str;
            str.str("");
            str << "\n Objective evaluation #" << Obj_Eval << "\tError: " << err;
            logfile_write(str.str(), 1);
        }
    } else {
        stringstream msg;
        msg.str("");
        msg << "\n\n WARNING: problem could not be solved!";
        logfile_write(msg.str(), 0);
        err = 1.e5;
    }

    // Save best into the best_logfile
    if (err < best_obj) {
        best_obj = err;
        ofstream best_logfile;
        best_logfile.open(best_logfilename.c_str(), ios::app);
        best_logfile << Obj_Eval << "; " << err << scientific << setprecision(5);
        for (unsigned int i = 0; i < genome.nPhenotypes(); i++)
            best_logfile << "; " << genome.phenotype(i);
        best_logfile << endl;
        best_logfile.close();
        PrintBestDataFile(genome);
    }

    Obj_Eval++;

    return err;
}

void PerformSensitivityAnalysis() {

    vector<double> tmp1;
    vector<double> tmp2;
    for (unsigned i = 0; i < wds.size(); i++) {
        tmp1.clear();
        tmp2.clear();
        wds.at(i)->Compute_Sensitivity_Matrix("diameter", 1);
        for (unsigned int j = 0; j < wds.at(i)->agelemek.size(); j++) {
            string name1 = wds.at(i)->agelemek.at(j)->Get_nev();
            for (unsigned int k = 0; k < pipe_name.size(); k++) {
                if (pipe_is_active.at(k)) {
                    string name2 = pipe_name.at(k);
                    if (0 == strcmp(name1.c_str(), name2.c_str())) {
                        tmp1.push_back(wds.at(i)->SM_row_sum_MassFlowRates.at(j));
                        tmp2.push_back(wds.at(i)->SM_row_sum_Pressures.at(j));
                        //cout<<endl<<"\t i="<<name1<<" ? "<<name2;
                        //cout<<" S_MFR: "<<wds.at(i)->SM_row_sum_MassFlowRates.at(j);
                        //cout<<" S_Prs: "<<wds.at(i)->SM_row_sum_Pressures.at(j);
                        break;
                    }
                }
            }
        }
        pipe_Sensitivity_MassFlowRates.push_back(tmp1);
        pipe_Sensitivity_Pressures.push_back(tmp2);
    }

    stringstream msg;
    msg.str("");
    msg << "\n Size of pipe_Sensitivity_MassFlowRates is " << pipe_Sensitivity_MassFlowRates.size() << " x "
         << pipe_Sensitivity_MassFlowRates.at(0).size();
    msg << "\n Size of pipe_Sensitivity_Pressures     is " << pipe_Sensitivity_Pressures.size() << " x "
         << pipe_Sensitivity_Pressures.at(0).size();
    logfile_write(msg.str(),1);

    // Now save to data files
    // First, Mass Flow Rate sensitivity
    ofstream resfile;
    resfile.open("sensitivity_MassFlowRates.dat", ios::out);
    resfile << scientific << setprecision(5);


    for (unsigned int j = 0; j < pipe_name.size(); j++)
        if (pipe_is_active.at(j))
            resfile << pipe_name.at(j) << " ; ";
    resfile << endl;

    vector<double> sum1(num_of_active_pipes, 0.0);

    for (unsigned int i = 0; i < wds.size(); i++) {
        int k = 0;
        for (unsigned int j = 0; j < pipe_name.size(); j++) {
            if (pipe_is_active.at(j)) {
                resfile << pipe_Sensitivity_MassFlowRates.at(i).at(k) << ";";
                sum1.at(k) += fabs(pipe_Sensitivity_MassFlowRates.at(i).at(k));
                k++;
            }
        }
        resfile << endl;
    }

    for (unsigned int j = 0; j < pipe_name.size(); j++)
        if (pipe_is_active.at(j))
            resfile << (sum1.at(j)/wds.size()) << ";";
    resfile.close();

    // Second, Pressure sensitivity
    resfile.open("sensitivity_Pressures.dat", ios::out);
    resfile << scientific << setprecision(5);

    for (unsigned int j = 0; j < pipe_name.size(); j++)
        if (pipe_is_active.at(j))
            resfile << pipe_name.at(j) << " ; ";
    resfile << endl;

    vector<double> sum2(num_of_active_pipes, 0.0);

    for (unsigned int i = 0; i < wds.size(); i++) {
        int k = 0;
        for (unsigned int j = 0; j < pipe_name.size(); j++) {
            if (pipe_is_active.at(j)) {
                resfile << pipe_Sensitivity_Pressures.at(i).at(k) << ";";
                sum2.at(k) += fabs(pipe_Sensitivity_Pressures.at(i).at(k));
                k++;
            }
        }
        resfile << endl;
    }

    for (unsigned int j = 0; j < pipe_name.size(); j++)
        if (pipe_is_active.at(j))
            resfile << (sum2.at(j)/wds.size()) << ";";
    resfile.close();

}

void Set_Up_Active_Pipes() {

    num_of_active_pipes = 0;

    for (unsigned int i = 0; i < wds.at(0)->agelemek.size(); i++)
        if (strcmp(wds.at(0)->agelemek.at(i)->Get_Tipus().c_str(), "Cso") == 0) {
            pipe_origD.push_back(wds.at(0)->agelemek.at(i)->Get_dprop("diameter"));
            pipe_name.push_back(wds.at(0)->agelemek.at(i)->Get_nev());
            if (wds.at(0)->agelemek.at(i)->Get_dprop("diameter") > Dmin) {
                pipe_is_active.push_back(true);
                num_of_active_pipes++;
            } else
                pipe_is_active.push_back(false);
        }

    stringstream msg;
    msg.str("");
    msg << "Setting active pipes (number: " << num_of_active_pipes << ")...\n";
    logfile_write(msg.str(), 0);
    msg.str("");
    for (unsigned int i = 0; i < pipe_name.size(); i++)
        if (pipe_is_active.at(i))
            msg << pipe_name.at(i) << " ";
    logfile_write(msg.str(), 1);
}

void Spoil_Active_Pipes() {

    stringstream msg;
    msg.str("");
    msg << endl << endl << "Spoiling active pipes....";
    for (unsigned int i = 0; i < wds.at(0)->agelemek.size(); i++) {
        string name1 = wds.at(0)->agelemek.at(i)->Get_nev();
        if (strcmp(wds.at(0)->agelemek.at(i)->Get_Tipus().c_str(), "Cso") == 0) {
            for (unsigned j = 0; j < pipe_name.size(); j++) {
                string name2 = pipe_name.at(j);
                if ((pipe_is_active.at(j)) && (0 == strcmp(name1.c_str(), name2.c_str()))) {
                    float mul = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
                    double Dold = wds.at(0)->agelemek.at(i)->Get_dprop("diameter");
                    double Dnew = (0.5 + mul) * Dold;
                    for (unsigned k = 0; k < wds.size(); k++)
                        wds.at(k)->agelemek.at(i)->Set_dprop("diameter", Dnew);
                    msg << endl << "\tSpoiling pipe " << name1 << " D = " << Dold << " -> " << Dnew << " m";
                }
            }
        }
    }
    logfile_write(msg.str(), 1);
}

double Compute_Error() {
    stringstream str_msg;
    str_msg.str("");
    str_msg << "\n\nComputing errors...\n";
    logfile_write(str_msg.str(), 2);
    str_msg.str("");

    // Extract data needed...
    double p_Sollwert, p_Staci;
    for (unsigned int i = 0; i < Sollwert_Node_Staci_ID.size(); i++) {
        for (unsigned int j = 0; j < Num_of_Periods; j++) {
            p_Sollwert = Sollwert_Node_Values.at(i).at(j);

            p_Staci = wds.at(j)->cspok.at(Sollwert_Node_Staci_Idx.at(i))->Get_p();

            Staci_Node_Values.at(i).at(j) = p_Staci;

            char msg[100];
            sprintf(msg, "\n\t %s (Staci name: %s): p(%2d)= FVM: %5.2f, Staci: %5.2f bar",
                    Sollwert_Node_Staci_ID.at(i).c_str(),
                    wds.at(j)->cspok.at(Sollwert_Node_Staci_Idx.at(i))->Get_nev().c_str(),
                    j, p_Sollwert, p_Staci);
            str_msg << msg;
        }
        str_msg << endl;
    }

    double H_FVM, H_Staci, H_err = 0.;
    for (unsigned int i = 0; i < Sollwert_Pool_Staci_ID.size(); i++) {
        for (unsigned int j = 0; j < Num_of_Periods; j++) {
            H_FVM = Sollwert_Pool_Values.at(i).at(j);
            H_Staci = wds.at(j)->agelemek.at(Sollwert_Pool_Staci_Idx.at(i))->Get_dprop("water_level");
            Staci_Pool_Values.at(i).at(j) = H_Staci;

            char msg[100];
            sprintf(msg, "\n\t %s (Staci name: %s): H(%2d)= FVM: %5.3f, Staci: %5.3f m",
                    Sollwert_Pool_Staci_ID.at(i).c_str(),
                    wds.at(j)->agelemek.at(Sollwert_Pool_Staci_Idx.at(i))->Get_nev().c_str(),
                    j, H_FVM, H_Staci);
            str_msg << msg;

            H_err += fabs(H_FVM - H_Staci);
        }
        str_msg << endl;
    }

    // Compute modified pressure error:
    double p_err = 0.0;
    for (unsigned int i = 0; i < Sollwert_Node_Staci_ID.size(); i++) {
        double sum = std::accumulate(Staci_Node_Values.at(i).begin(), Staci_Node_Values.at(i).end(), 0.0);
        double Staci_mean = sum / Staci_Node_Values.at(i).size();

        sum = std::accumulate(Sollwert_Node_Values.at(i).begin(), Sollwert_Node_Values.at(i).end(), 0.0);
        double FVM_mean = sum / Sollwert_Node_Values.at(i).size();
        for (unsigned int j = 0; j < Num_of_Periods; j++)
            p_err += fabs(
                    (Staci_Node_Values.at(i).at(j) - Staci_mean) - (Sollwert_Node_Values.at(i).at(j) - FVM_mean));

    }

    double weight_p_err = 0.01;
    double err = weight_p_err * p_err + (1.0 - weight_p_err) * H_err;
    str_msg.str("");
    str_msg << "\n  p_err = " << p_err << ", weigth= " << weight_p_err;
    str_msg << "\n  H_err = " << H_err << ", weigth= " << (1.0 - weight_p_err);
    str_msg << "\n  err   = " << err << endl;

    logfile_write(str_msg.str(), 2);

    return err;
}

void Load_WDS_Systems() {
    stringstream fname;
    stringstream msg;

    logfile_write("Building WDS systems....\n", 0);

    for (unsigned int i = 0; i < Num_of_Periods; i++) {
        fname.str("");
        fname << dir_name << fname_prefix << (Start_of_Periods + i) << ".spr";

        msg.str("");
        msg << "\n\t" << fname.str();
        logfile_write(msg.str(), 2);

        wds.push_back(new Staci(fname.str()));
        wds.at(i)->build_system();
        wds.at(i)->ini();
    }
    logfile_write("\n", 2);
}

void Update_Reservoirs(unsigned int i) {

    logfile_write("\n\nUpdating Pools...", 2);

    for (unsigned int j = 0; j < wds.at(i)->agelemek.size(); j++) {

        string type = wds.at(i)->agelemek.at(j)->Get_Tipus();
        if (strcmp(type.c_str(), "Vegakna") == 0) {
            string PoolName = wds.at(i)->agelemek.at(j)->Get_nev();
            double mp_prev = wds.at(i - 1)->agelemek.at(j)->Get_dprop("mass_flow_rate");
            double Q = mp_prev / 1000. * 3600.;
            double H_prev = wds.at(i - 1)->agelemek.at(j)->Get_dprop("water_level");
            double A = get_A(PoolName);
            double H_next = H_prev + Q * dt / A;

            wds.at(i)->agelemek.at(j)->Set_dprop("water_level", H_next);

            char msg[100];
            sprintf(msg, "\n\t%s: t=%2d, Q=%+5.0fm3/h, A=%4.0fm^2, H(%2d)=%4.2fm -> H(%2d)=%4.2fm",
                    (wds.at(i)->agelemek.at(j)->Get_nev()).c_str(), i - 1, Q, A, i - 1, H_prev, i, H_next);
            string str_msg(msg);
            logfile_write(str_msg, 2);
        }
    }
}

void Load_Sollwert_Datafiles() {
    stringstream tmp;
    vector<vector<string> > lines;

    tmp.str("");
    tmp << dir_name << sollwert_dfile;

    std::ifstream in(tmp.str().c_str());

    if (in.fail()) {
        cout << endl << "Load_Sollwert_Datafiles() -> CANNOT FIND DATAFILE " << tmp.str() << "!!!!" << endl;
        exit(-1);
    } else {
        while (in.good())
            lines.push_back(csv_read_row(in, ';'));
    }
    in.close();

    stringstream tmpstr;
    vector<double> tmp_double_vec;

    tmpstr.str("");
    tmpstr << "\n Reading sollwert file " << sollwert_dfile << " ... ";

    for (unsigned int i = 0; i < lines.size(); i++) {
        if (lines.at(i).size() < (Start_of_Periods + Num_of_Periods + 2 + 1)) {
            tmpstr << endl << endl << "!!! ERROR !!!" << endl;
            tmpstr << "\t in Sollwert file " << sollwert_dfile << ", line " << i << endl;
            tmpstr << "\t The line contains " << lines.at(i).size()
                   << " data (separated by ';'), but should be minimum "
                   << Start_of_Periods + Num_of_Periods + 2 << ", " << endl;
            tmpstr << " because Start_of_Periods is " << Start_of_Periods;
            tmpstr << " and Num_of_Periods is " << Num_of_Periods << ". " << endl << endl;
            logfile_write(tmpstr.str(), 0);
            exit(-1);
        }

        string name = lines.at(i).at(0);
        string type = lines.at(i).at(1);
        std::transform(type.begin(), type.end(), type.begin(), ::tolower);

        tmp_double_vec.clear();

        // POOL:
        if (0 == strcmp(type.c_str(), "pool")) {
            Sollwert_Pool_Staci_ID.push_back(name);
            Sollwert_Pool_Staci_Idx.push_back(Find_Pool_Index(name));
            for (unsigned int j = Start_of_Periods; j < Start_of_Periods + Num_of_Periods; j++)
                tmp_double_vec.push_back(atof(lines.at(i).at(2 + j).c_str()));
            Sollwert_Pool_Values.push_back(tmp_double_vec);

            // Add surface
            Pool_Surfs.push_back(get_A(name));

            // Info
            tmpstr << endl << "\t Found pool: " << name << " H(0) = " << tmp_double_vec.at(0) << " m , A = "
                   << get_A(name) << " m2";

            // Make space for Staci values:
            vector<double> aa(tmp_double_vec.size(), 0.);
            Staci_Pool_Values.push_back(aa);
        }

        // NODE:
        if (0 == strcmp(type.c_str(), "node")) {
            Sollwert_Node_Staci_ID.push_back(name);
            Sollwert_Node_Staci_Idx.push_back(Find_Pressure_Index(name));
            for (unsigned int j = Start_of_Periods; j < Start_of_Periods + Num_of_Periods; j++)
                tmp_double_vec.push_back(atof(lines.at(i).at(2 + j).c_str()));
            Sollwert_Node_Values.push_back(tmp_double_vec);

            // Info
            tmpstr << endl << "\t Found node: " << lines.at(i).at(0);

            // Make space for Staci values:
            vector<double> aa(tmp_double_vec.size(), 0.);
            Staci_Node_Values.push_back(aa);
        }
    }
    logfile_write(tmpstr.str(), 1);
}

double get_A(string PoolName) {

    bool found = false;
    double A;

    for (unsigned int i = 0; i < wds.at(0)->agelemek.size(); i++) {
        string name = wds.at(0)->agelemek.at(i)->Get_nev();
        if (strcmp(name.c_str(), PoolName.c_str()) == 0) {
            A = wds.at(0)->agelemek.at(i)->Get_Aref();
            found = true;

            stringstream msg;
            msg.str("");
            msg << "\n\t get_A(): " << PoolName << " was found.";
            logfile_write(msg.str(), 3);

            break;
        }
    }

    if (!found) {
        stringstream msg;
        msg.str("");
        msg << "\n\n get_A(): could not find " << PoolName << "\n\n";
        logfile_write(msg.str(), 0);
        exit(-1);
    }

    return A;
}

int Find_Pool_Index(string PoolName) {
    bool found = false;
    int idx = -1;
    for (unsigned int j = 0; j < wds.at(0)->agelemek.size(); j++) {
        string type = wds.at(0)->agelemek.at(j)->Get_Tipus();
        if (strcmp(type.c_str(), "Vegakna") == 0) {
            string Name = wds.at(0)->agelemek.at(j)->Get_nev();
            if (strcmp(Name.c_str(), PoolName.c_str()) == 0) {
                found = true;
                idx = j;
                break;
            }
        }
    }
    if (!found) {
        stringstream msg;
        msg.str("");
        msg << "\nFind_Pool_Index -> ERROR: " << PoolName << " was not found!!!\n\n";
        logfile_write(msg.str(), 0);
        exit(-1);
    }
    return idx;
}

int Find_Pressure_Index(string PressureName) {
    bool found = false;
    int idx = -1;
    for (unsigned int j = 0; j < wds.at(0)->cspok.size(); j++) {
        string Name = wds.at(0)->cspok.at(j)->Get_nev();
        //cout<<endl<<Name;
        //cin.get();
        if (strcmp(Name.c_str(), PressureName.c_str()) == 0) {
            found = true;
            idx = j;
            //cin.get();
            break;
        }
    }
    if (!found) {
        stringstream msg;
        msg.str("");
        msg << "\nFind_Pressure_Index -> ERROR: " << PressureName << " was not found!!!\n\n";
        logfile_write(msg.str(), 0);
        exit(-1);
    }
    return idx;
}


void Set_Initial_Pool_Levels() {

    stringstream str;
    str.str("");
    logfile_write("Setting up initial pool levels...\n", 0);

    for (unsigned j = 0; j < Sollwert_Pool_Staci_ID.size(); j++) {
        str.str("");
        str << "\n\t Pool name: " << Sollwert_Pool_Staci_ID.at(j);

        double H0 = Sollwert_Pool_Values.at(j).at(0);
        int idx = Sollwert_Pool_Staci_Idx.at(j);
        wds.at(0)->agelemek.at(idx)->Set_dprop("water_level", H0);

        str << ", i.e. agelemek.at(" << idx << ") aka " << wds.at(0)->agelemek.at(idx)->Get_nev() << " -> H0 = "
            << H0;
        logfile_write(str.str(), 1);
    }
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

vector<string> csv_read_row(string &line, char delimiter) {
    stringstream ss(line);
    return csv_read_row(ss, delimiter);
}

vector<string> csv_read_row(istream &in, char delimiter) {
    stringstream ss;
    bool inquotes = false;
    vector<string> row;//relying on RVO
    while (in.good()) {
        char c = in.get();
        if (!inquotes && c == '"') //beginquotechar
        {
            inquotes = true;
        } else if (inquotes && c == '"') //quotechar
        {
            if (in.peek() == '"')//2 consecutive quotes resolve to 1
            {
                ss << (char) in.get();
            } else //endquotechar
            {
                inquotes = false;
            }
        } else if (!inquotes && c == delimiter) //end of field
        {
            row.push_back(ss.str());
            ss.str("");
        } else if (!inquotes && (c == '\r' || c == '\n')) {
            if (in.peek() == '\n') { in.get(); }
            row.push_back(ss.str());
            return row;
        } else {
            ss << c;
        }
    }
}

string Load_Settings() {

    XMLNode xMainNode = XMLNode::openFileHelper("staci_calibrate_settings.xml", "settings");

    global_debug_level = atoi(xMainNode.getChildNode("global_debug_level").getText());
    Staci_debug_level = atoi(xMainNode.getChildNode("Staci_debug_level").getText());
    dir_name = xMainNode.getChildNode("dir_name").getText();
    fname_prefix = xMainNode.getChildNode("fname_prefix").getText();
    logfilename = xMainNode.getChildNode("logfilename").getText();
    best_logfilename = xMainNode.getChildNode("best_logfilename").getText();
    sollwert_dfile = xMainNode.getChildNode("sollwert_dfile").getText();
    Start_of_Periods = atoi(xMainNode.getChildNode("Start_of_Periods").getText());
    Num_of_Periods = atoi(xMainNode.getChildNode("Num_of_Periods").getText());
    string Text_Spoil_Active_Pipes = xMainNode.getChildNode("Spoil_Active_Pipes").getText();
    dt = atof(xMainNode.getChildNode("dt").getText());
    Dmin = atof(xMainNode.getChildNode("Dmin").getText());
    popsize = atoi(xMainNode.getChildNode("popsize").getText());
    ngen = atoi(xMainNode.getChildNode("ngen").getText());
    pmut = atof(xMainNode.getChildNode("pmut").getText());
    pcross = atof(xMainNode.getChildNode("pcross").getText());

    if (0 == strcmp(Text_Spoil_Active_Pipes.c_str(), "yes"))
        bool_Spoil_Active_Pipes = true;
    else
        bool_Spoil_Active_Pipes = false;


    stringstream msg;
    msg.str("");
    msg << "Settings:" << endl;
    msg << "\t global_debug_level     : " << global_debug_level << endl;
    msg << "\t Staci_debug_level      : " << Staci_debug_level << endl << endl;
    msg << "\t dir_name               : " << dir_name << endl;
    msg << "\t fname_prefix           : " << fname_prefix << endl;
    msg << "\t logfilename            : " << logfilename << endl;
    msg << "\t best_logfilename       : " << best_logfilename << endl;
    msg << "\t sollwert_dfile         : " << sollwert_dfile << endl;
    msg << "\t Start_of_Periods       : " << Start_of_Periods << endl;
    msg << "\t Num_of_Periods         : " << Num_of_Periods << endl;
    msg << "\t length of periods (dt) : " << dt << " hour " << endl;
    msg << "\t spoil active pipes     : " << bool_Spoil_Active_Pipes << endl;
    msg << "\t Dmin                   : " << Dmin << " m" << endl << endl;
    msg << "\t popsize : " << popsize << endl;
    msg << "\t ngen    : " << ngen << endl;
    msg << "\t pmut    : " << pmut << endl;
    msg << "\t pcross  : " << pcross << endl << endl;

    return msg.str();
}
