#include <stdio.h>
#include <vector>
#include "Staci.h"

using namespace std;

vector<Staci *> wds;
int Staci_debug_level;

#define cout STD_COUT
/*
int global_debug_level;


string dir_name;
string fname_prefix;
string logfilename;
string best_logfilename;
string Eredmenyek_FVM_dfile;
string FVM_dfile;

void Load_WDS_Systems();

void Load_FVM_Datafiles();

void Update_Reservoirs(unsigned int);

void logfile_write(string msg, int debug_level);

int Num_of_Periods;
//vector<string> pipe_names;


vector<string> FVM_Pressure_StaciID;
vector<int> FVM_Pressure_Staci_Idx;
vector<string> FVM_Pressure_Info;
vector<vector<double> > FVM_Pressure_Values;
vector<vector<double> > Staci_Pressure_Values;

vector<string> FVM_Pool_StaciID;
vector<int> FVM_Pool_Staci_Idx;
vector<string> FVM_Pool_Info;
vector<vector<double> > FVM_Pool_tmp_Values;
vector<vector<double> > FVM_Pool_Values;
vector<vector<double> > Staci_Pool_Values;

vector<string> FVM_Pool_IDs;
vector<string> FVM_Pool_Names;
vector<int> FVM_Pool_Name_Idx; // this staci element refers to the actual FVM Pool

vector<double> PoolSurfs;

vector<string> pipe_name;
vector<bool> pipe_is_active;
unsigned int num_of_active_pipes;
vector<double> pipe_origD;

void Set_Up_Active_Pipes();

double Dmin;

int Find_Pressure_Index(string PressureName);

int Find_Pool_Index(string PressureName);

void Define_Pools();

double get_A(string PoolName);

double Compute_Error();

float Objective(GAGenome &);

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
float best_obj;*/

int main(int argc, char **argv) {

    Staci_debug_level = 0;

    // Load WDSs and data files
    vector <string> orig_dfile_names;
    orig_dfile_names.push_back("Anytown_DC1.spr");


    // Initial computation
    for (unsigned i = 0; i < orig_dfile_names.size(); i++) {
        wds.push_back(new Staci(orig_dfile_names.at(i)));
        wds.at(i)->build_system();
        wds.at(i)->ini();
        //wds.at(i)->Set_debug_level(Staci_debug_level);
    }


    return 0;
}
/*
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
        resfile << genome.phenotype(j)/pipe_origD.at(j) << " ; ";
    resfile << endl;


    for (unsigned int i = 0; i < FVM_Pool_Names.size(); i++)
        resfile << "FVM " << FVM_Pool_Names.at(i) << " ; " << FVM_Pool_Names.at(i) << " ; ";
    for (unsigned int i = 0; i < FVM_Pressure_StaciID.size(); i++)
        resfile << "FVM " << FVM_Pressure_StaciID.at(i) << " ; " << FVM_Pressure_StaciID.at(i) << " ; ";

    for (unsigned int j = 0; j < Num_of_Periods; j++) {
        resfile << endl;
        for (unsigned int i = 0; i < FVM_Pool_Names.size(); i++)
            resfile << FVM_Pool_Values.at(i).at(j) << " ; " << Staci_Pool_Values.at(i).at(j) << " ; ";
        for (unsigned int i = 0; i < FVM_Pressure_StaciID.size(); i++)
            resfile << FVM_Pressure_Values.at(i).at(j) << " ; " << Staci_Pressure_Values.at(i).at(j) << " ; ";
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
            Update_Reservoirs(i);
        } else {
            break;
        }
    }

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
    msg << "\n Active pipes (number: " << num_of_active_pipes << "):\n";
    for (unsigned int i = 0; i < pipe_name.size(); i++)
        if (pipe_is_active.at(i))
            msg << pipe_name.at(i) << " ";
    logfile_write(msg.str(), 0);
}

double Compute_Error() {
    stringstream str_msg;
    str_msg.str("");
    str_msg << "\n\nComputing errors...\n";
    logfile_write(str_msg.str(), 2);
    str_msg.str("");

    // Extract data needed...
    double p_FVM, p_Staci;
    for (unsigned int i = 0; i < FVM_Pressure_StaciID.size(); i++) {
        for (unsigned int j = 0; j < Num_of_Periods; j++) {
            p_FVM = FVM_Pressure_Values.at(i).at(j);
            p_Staci = wds.at(j)->cspok.at(FVM_Pressure_Staci_Idx.at(i))->Get_p();
            Staci_Pressure_Values.at(i).at(j) = p_Staci;

            char msg[100];
            sprintf(msg, "\n\t %s (Staci name: %s): p(%2d)= FVM: %5.2f, Staci: %5.2f bar",
                    FVM_Pressure_StaciID.at(i).c_str(),
                    wds.at(j)->cspok.at(FVM_Pressure_Staci_Idx.at(i))->Get_nev().c_str(),
                    j, p_FVM, p_Staci);
            str_msg << msg;
        }
        str_msg << endl;
    }

    double H_FVM, H_Staci, H_err = 0.;
    for (unsigned int i = 0; i < FVM_Pool_StaciID.size(); i++) {
        for (unsigned int j = 0; j < Num_of_Periods; j++) {
            H_FVM = FVM_Pool_Values.at(i).at(j);
            H_Staci = wds.at(j)->agelemek.at(FVM_Pool_Staci_Idx.at(i))->Get_dprop("water_level");
            Staci_Pool_Values.at(i).at(j) = H_Staci;

            char msg[100];
            sprintf(msg, "\n\t %s (Staci name: %s): H(%2d)= FVM: %5.3f, Staci: %5.3f m",
                    FVM_Pool_StaciID.at(i).c_str(),
                    wds.at(j)->agelemek.at(FVM_Pool_Staci_Idx.at(i))->Get_nev().c_str(),
                    j, H_FVM, H_Staci);
            str_msg << msg;

            H_err += fabs(H_FVM - H_Staci);
        }
        str_msg << endl;
    }
    logfile_write(str_msg.str(), 2);

    // Compute modified pressure error:
    double p_err = 0.0;
    for (unsigned int i = 0; i < FVM_Pressure_StaciID.size(); i++) {
        double sum = std::accumulate(Staci_Pressure_Values.at(i).begin(), Staci_Pressure_Values.at(i).end(), 0.0);
        double Staci_mean = sum / Staci_Pressure_Values.at(i).size();

        sum = std::accumulate(FVM_Pressure_Values.at(i).begin(), FVM_Pressure_Values.at(i).end(), 0.0);
        double FVM_mean = sum / FVM_Pressure_Values.at(i).size();
        for (unsigned int j = 0; j < Num_of_Periods; j++)
            p_err += fabs(
                    (Staci_Pressure_Values.at(i).at(j) - Staci_mean) - (FVM_Pressure_Values.at(i).at(j) - FVM_mean));

    }

    double weight_p_err = 1.0;
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

    logfile_write("Building WDS systems....", 0);

    for (unsigned int i = 0; i < Num_of_Periods; i++) {
        fname.str("");
        fname << dir_name << fname_prefix << (i + 1) << ".spr";

        msg.str("");
        msg << "\n\t" << fname.str();
        logfile_write(msg.str(), 0);

        wds.push_back(new Staci(fname.str()));
        wds.at(i)->build_system();
        wds.at(i)->ini();
    }
    logfile_write("\n", 0);
}

void Update_Reservoirs(unsigned int i) {

    logfile_write("\n\nUpdating Pools...", 2);

    for (unsigned int j = 0; j < wds.at(i)->agelemek.size(); j++) {
        string type = wds.at(i)->agelemek.at(j)->Get_Tipus();
        if (strcmp(type.c_str(), "Vegakna") == 0) {
            string PoolName = wds.at(i)->agelemek.at(j)->Get_nev();
            double mp_prev = wds.at(i)->agelemek.at(j)->Get_dprop("mass_flow_rate");
            double Q = mp_prev / 1000. * 3600.;
            double H_prev = wds.at(i)->agelemek.at(j)->Get_dprop("water_level");
            double A = get_A(PoolName);
            double H_next = H_prev + Q * 0.5 / A;

            if (i < Num_of_Periods - 1)
                wds.at(i + 1)->agelemek.at(j)->Set_dprop("water_level", H_next);

            char msg[100];
            sprintf(msg, "\n\t%s: t=%2d, Q=%+5.0fm3/h, A=%4.0fm^2, H(%2d)=%4.2fm -> H(%2d)=%4.2fm",
                    (wds.at(i)->agelemek.at(j)->Get_nev()).c_str(), i, Q, A, i, H_prev, i + 1, H_next);
            string str_msg(msg);
            logfile_write(str_msg, 2);
        }
    }
}

void Load_FVM_Datafiles() {
    stringstream tmp, msg;
    vector<vector<string> > lines;

    // POOL DATA FILE
    tmp.str("");
    tmp << dir_name << FVM_dfile;

    std::ifstream in(tmp.str().c_str());

    if (in.fail()) {
        cout << endl << "Load_FVM_Datafiles() -> CANNOT FIND DATAFILE " << tmp.str() << "!!!!" << endl;
        exit(-1);
    } else {
        while (in.good())
            lines.push_back(csv_read_row(in, ';'));
    }
    in.close();

    string tmpstr;
    for (unsigned int i = 0; i < lines.size(); i++) {
        char first_char = lines.at(i).at(0).c_str()[0];
        char second_char = lines.at(i).at(1).c_str()[0];
        //cout << "\n line " << i << "/" << lines.size() << " first_char=" << first_char << ", second_char="
        //     << second_char;
        if ((first_char == '#') && (second_char != '-')) {
            FVM_Pool_StaciID.push_back(lines.at(i).at(0));
            tmpstr = lines.at(i).at(1);
            tmpstr.append(", ");
            tmpstr.append(lines.at(i).at(2));
            FVM_Pool_Info.push_back(tmpstr);
            vector<double> vals;
            for (unsigned int j = 3; j < lines.at(i).size(); j++)
                vals.push_back(atof(lines.at(i).at(j).c_str()));
            FVM_Pool_tmp_Values.push_back(vals);
            //cout << " added " << vals.size() << " data";
            // Make space for Staci values
            vector<double> tmpvec(vals.size());
            Staci_Pool_Values.push_back(tmpvec);
        }
    }
    msg.str("");
    msg << "\nReading datafile " << tmp.str() << " ...";
    for (unsigned int i = 0; i < FVM_Pool_StaciID.size(); i++) {
        msg << "\n " << FVM_Pool_StaciID.at(i) << ", number of data: "
            << FVM_Pool_tmp_Values.at(i).size() << ", info: " << FVM_Pool_Info.at(i);
        if (global_debug_level > 1) {
            msg << endl << "\t";
            for (unsigned int j = 0; j < FVM_Pool_tmp_Values.at(i).size(); j++)
                msg << FVM_Pool_tmp_Values.at(i).at(j) << ", ";
        }
    }
    logfile_write(msg.str(), 0);

    // PRESSURE DATA FILE
    lines.clear();
    tmp.str("");
    tmp << dir_name << Eredmenyek_FVM_dfile;

    std::ifstream in1(tmp.str().c_str());
    //in1(tmp.str());

    if (in1.fail()) {
        cout << endl << "Load_FVM_Datafiles() -> CANNOT FIND DATAFILE " << tmp.str() << "!!!!" << endl;
        exit(-1);
    } else
        while (in1.good())
            lines.push_back(csv_read_row(in1, ';'));
    in1.close();


    //string tmpstr;
    for (unsigned int i = 1; i < lines.size(); i++) {
        string name = lines.at(i).at(1);
        if (0 != (strcmp(name.c_str(), "x"))) {
            //cout<<endl<<"name:"<<name;
            //cin.get();
            //int idx = Find_Pressure_Index(name);
            //cout<< ", idx="<<idx;//<<", staci name: "<<wds.at(0)->cspok.at(idx)->Get_nev();
            //cin.get();
            FVM_Pressure_StaciID.push_back(name);
            FVM_Pressure_Staci_Idx.push_back(Find_Pressure_Index(name));
            tmpstr = lines.at(i).at(2);
            tmpstr.append(", ");
            tmpstr.append(lines.at(i).at(3));
            FVM_Pressure_Info.push_back(tmpstr);
            vector<double> vals;
            for (unsigned int j = 4; j < lines.at(i).size(); j++)
                vals.push_back(atof(lines.at(i).at(j).c_str()));
            FVM_Pressure_Values.push_back(vals);
            // Make space for Staci values
            vector<double> tmpvec(vals.size());
            Staci_Pressure_Values.push_back(tmpvec);
        }
    }

    msg.str("");
    msg << "\n\nReading datafile " << tmp.str() << " ...";
    for (unsigned int i = 0; i < FVM_Pressure_StaciID.size(); i++) {
        msg << "\n " << FVM_Pressure_StaciID.at(i) << ", number of data: "
            << FVM_Pressure_Values.at(i).size() << ", info: " << FVM_Pressure_Info.at(i);
        if (global_debug_level > 1) {
            msg << endl << "\t";
            for (unsigned int j = 0; j < FVM_Pressure_Values.at(i).size(); j++)
                msg << FVM_Pressure_Values.at(i).at(j) << ", ";
        }
    }
    logfile_write(msg.str(), 0);
}

double get_A(string PoolName) {

    bool found = false;
    double A;

    for (unsigned int i = 0; i < FVM_Pool_IDs.size(); i++)
        if (strcmp(FVM_Pool_IDs.at(i).c_str(), PoolName.c_str()) == 0) {
            A = PoolSurfs.at(i);
            found = true;
            break;

            stringstream msg;
            msg.str("");
            msg << "\n\t get_A(): " << PoolName << " was found.";
            logfile_write(msg.str(), 2);
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


void Define_Pools() {

    FVM_Pool_IDs.push_back("POOL288");
    FVM_Pool_Names.push_back("#med_cinkota");
    FVM_Pool_Staci_Idx.push_back(Find_Pool_Index("POOL288"));
    PoolSurfs.push_back(1670.);

    FVM_Pool_IDs.push_back("POOL315");
    FVM_Pool_Names.push_back("#med_gilice");
    FVM_Pool_Staci_Idx.push_back(Find_Pool_Index("POOL315"));
    PoolSurfs.push_back(2783);

    FVM_Pool_IDs.push_back("POOL327");
    FVM_Pool_Names.push_back("#med_kob_regi");
    FVM_Pool_Staci_Idx.push_back(Find_Pool_Index("POOL327"));
    PoolSurfs.push_back(2750. / 2.);

    FVM_Pool_IDs.push_back("POOL324");
    FVM_Pool_Names.push_back("#med_kob_uj");
    FVM_Pool_Staci_Idx.push_back(Find_Pool_Index("POOL324"));
    PoolSurfs.push_back(3333.);

    FVM_Pool_IDs.push_back("POOL276");
    FVM_Pool_Names.push_back("#med_krisztina");
    FVM_Pool_Staci_Idx.push_back(Find_Pool_Index("POOL276"));
    PoolSurfs.push_back(3767.);

    FVM_Pool_IDs.push_back("POOL303");
    FVM_Pool_Names.push_back("#med_rakossz");
    FVM_Pool_Staci_Idx.push_back(Find_Pool_Index("POOL303"));
    PoolSurfs.push_back(1667.);

    FVM_Pool_IDs.push_back("POOL336");
    FVM_Pool_Names.push_back("#med_sanc");
    FVM_Pool_Staci_Idx.push_back(Find_Pool_Index("POOL336"));
    PoolSurfs.push_back(1000.);

    stringstream str;
    str.str("");
    str << "\n\n Setting up pools...";
    logfile_write(str.str(), 0);
    for (unsigned j = 0; j < FVM_Pool_Names.size(); j++) {
        str.str("");
        str << "\n\t Pool name: " << FVM_Pool_Names.at(j) << ", \tID: " << FVM_Pool_IDs.at(j);
        bool found = false;
        for (unsigned int i = 0; i < FVM_Pool_tmp_Values.size(); i++)
            if (0 == strcmp(FVM_Pool_StaciID.at(i).c_str(), FVM_Pool_Names.at(j).c_str())) {
                //cout << endl << FVM_Pool_StaciID.at(i) << ": H(0)=" << FVM_Pool_tmp_Values.at(i).at(0);
                FVM_Pool_Values.push_back(FVM_Pool_tmp_Values.at(i));
                found = true;
                break;
            }
        if (!found) {
            cout << "\n WTF??????";
            exit(-1);
        }
        double H0 = FVM_Pool_Values.at(j).at(0);
        int idx = FVM_Pool_Staci_Idx.at(j);
        wds.at(0)->agelemek.at(idx)->Set_dprop("water_level", H0);
        str << ", i.e. agelemek.at(" << idx << ") aka " << wds.at(0)->agelemek.at(idx)->Get_nev() << " -> H0 = " << H0;

        logfile_write(str.str(), 0);
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
    Eredmenyek_FVM_dfile = xMainNode.getChildNode("Eredmenyek_FVM_dfile").getText();
    FVM_dfile = xMainNode.getChildNode("FVM_dfile").getText();
    Num_of_Periods = atoi(xMainNode.getChildNode("Num_of_Periods").getText());
    Dmin = atof(xMainNode.getChildNode("Dmin").getText());
    popsize = atoi(xMainNode.getChildNode("popsize").getText());
    ngen = atoi(xMainNode.getChildNode("ngen").getText());
    pmut = atof(xMainNode.getChildNode("pmut").getText());
    pcross = atof(xMainNode.getChildNode("pcross").getText());

    stringstream msg;
    msg.str("");
    msg << "Settings:" << endl;
    msg << "\t global_debug_level   : " << global_debug_level << endl;
    msg << "\t Staci_debug_level    : " << Staci_debug_level << endl << endl;
    msg << "\t dir_name             : " << dir_name << endl;
    msg << "\t fname_prefix         : " << fname_prefix << endl;
    msg << "\t logfilename          : " << logfilename << endl;
    msg << "\t best_logfilename     : " << best_logfilename << endl;
    msg << "\t Eredmenyek_FVM_dfile : " << FVM_dfile << endl;
    msg << "\t Num_of_Periods       : " << Num_of_Periods << endl;
    msg << "\t Dmin                 : " << Dmin << endl << endl;
    msg << "\t popsize : " << popsize << endl;
    msg << "\t ngen    : " << ngen << endl;
    msg << "\t pmut    : " << pmut << endl;
    msg << "\t pcross  : " << pcross << endl << endl;

    return msg.str();
}*/
/*

vector<int> data = {5, 16, 4, 7};
vector<int> index(data.size(), 0);
for (int i = 0 ; i != index.size() ; i++) {
index[i] = i;
}
sort(index.begin(), index.end(),
[&](const int& a, const int& b) {
return (data[a] < data[b]);
}
);
for (int i = 0 ; i != index.size() ; i++) {
cout << index[i] << endl;
}*/
