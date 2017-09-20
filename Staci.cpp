#include "Staci.h"
#include <string.h>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <algorithm>
#include "StaciException.h"
#include "data_io.h"

// Wiley:
//#include "/usr/include/suitesparse/umfpack.h"
// Mac:
#include "/usr/local/include/umfpack.h"
// #include <suitesparse/umfpack.h>

// TODO: Clear up debug level behaviour in the doc!
// TODO: log file cleanup!

struct val_and_ID {
  double val, x;
  string ID;
  bool operator>(const val_and_ID& rhs) { return val > rhs.val; }
  bool operator>=(const val_and_ID& rhs) { return val >= rhs.val; }
  bool operator<(const val_and_ID& rhs) { return val < rhs.val; }
  bool operator<=(const val_and_ID& rhs) { return val <= rhs.val; }
};

bool comparison_function1(const val_and_ID& lhs, const val_and_ID& rhs ) { return lhs.val > rhs.val; }

Staci::Staci(int argc, char *argv[]) {
  // Initiate with command line arguments
  opt = new AnyOption();
  get_command_line_options(argc, argv);
  SetInitialParameters();
}

Staci::Staci(string spr_filename) {
  mode = 0;
  // Just to be able to delete in the destructor
  opt = new AnyOption();

  // Initiate directly from file
  def_file = spr_filename;
  SetInitialParameters();
}

void Staci::SetInitialParameters() {
  van_ini = false;

  m_nnz = 0; /*!< Number of nonzero entries of the Jacobian. */
  // az adatfile olvasasa
  if (mode >= 0) {

    data_io datta_io(def_file.c_str());

    datta_io.load_system(cspok, agelemek);

    // beallitasok
    // cout<<endl<<endl<<"Beallitasok:";
    debug_level = atoi(datta_io.read_setting("debug_level").c_str());
    // cout<<endl<<"\tdebug_level  => "<<datta_io.read_setting("debug_level");

    out_file = datta_io.read_setting("out_file");
    // cout<<endl<<"\tout_file     => "<<out_file;

    iter_max = atoi(datta_io.read_setting("iter_max").c_str());
    // cout<<endl<<"\titer_max     => "<<iter_max;

    e_p_max = atof(datta_io.read_setting("e_p_max").c_str());
    // cout<<endl<<"\te_p_max      => "<<e_p_max;

    e_mp_max = atof(datta_io.read_setting("e_mp_max").c_str());
    // cout<<endl<<"\te_mp_max     => "<<e_mp_max;

    m_relax = atof(datta_io.read_setting("relax").c_str());
    //     cout<<endl<<"\trelax        => "<<relax;

    m_relax_mul = atof(datta_io.read_setting("relax_mul").c_str());
    //      cout<<endl<<"\trelax szorzo => "<<relax_mul;

    mp_init = atof(datta_io.read_setting("mp_init").c_str());
    //      cout<<endl<<"\tmp ini       => "<<mp_init;

    p_init = atof(datta_io.read_setting("p_init").c_str());
    //      cout<<endl<<"\tp ini        => "<<p_init;

    friction_model = datta_io.read_setting("friction_model").c_str();

    if (friction_model == "nincs ilyen node!") {
      friction_model = "DW";
      cout << endl << endl << "******** WARNING! ********";
      cout << endl << "<friction_model>DW|HW</friction_model> node not found!";
      cout << endl
           << "Setting DW friction model for all pipes." << endl
           << endl;
    }
  }

  // Set up the log & progress file name
  set_out_file(def_file + ".ros");
  set_progress_file(def_file + ".rps");

  // Delete if previous exists...
  ostringstream msg;
  msg << "\nTrying to delete previous logfile " << out_file.c_str() << "... ";
  if (remove(out_file.c_str()) != 0)
    msg << " file not found, cannot delete it.";
  else
    msg << "file successfully deleted.";

  // Add header
  m_ss.str("");
  m_ss << endl << "======================================";
  m_ss << endl << " STACI v.2.0";
  m_ss << endl << " (c) BME Dept. of Hydrodynamic Systems";
  m_ss << endl << " (c) C. Hos (cshos@hds.bme.hu)";
  m_ss << endl << " info: www.hds.bme.hu\\staci_web";
  m_ss << endl << "======================================" << endl;
  time_t ido = time(0);
  m_ss << endl << " date: " << ctime(&ido);
  m_ss << endl << " input file: " << def_file << endl;
  if (debug_level > 0) {
    // cout << m_ss.str();
    logfile_write(m_ss.str(), 1);
  }
  m_ss.str("");

  // If set to false, the result file will not be saved!
  // Only set false for sensitivity analysis
  do_save_file = true;
}

//--------------------------------------------------------------
Staci::~Staci() { delete opt; }

//--------------------------------------------------------------
void Staci::get_command_line_options(int argc, char *argv[]) {
  opt->setVerbose();         /* print warnings about unknown options */
  opt->autoUsagePrint(true); /* print usage for bad options */

  /// This function sets members 'mode' and 'def_file'.
  /// Possible argumentums:
  /// -h --help
  /// -s --stac              (mode=0)
  /// -i --ini
  /// -t --travel_time       (mode=1)
  /// -c --conc_transp       (mode=2)
  /// -m --mod_prop          (mode=3)
  ///    -e --element_ID
  ///    -p --property_ID
  ///    -n --newValue
  ///    -o --outfile
  /// -l --list_all_elements (mode=4)
  /// -g --get_data          (mode=5)
  ///    -e --element_ID
  ///    -p --property_ID
  /// -r sensitivity         (mode=6)
  /// -e export community structure (mode=7)

  opt->addUsage("");
  opt->addUsage("staci hasznalata: ");
  opt->addUsage("");
  opt->addUsage("\t help:");
  opt->addUsage(
    "\t\t -h  --help                          Help nyomtatasa a kepernyore");
  opt->addUsage(" ");
  opt->addUsage("\t stacioner halozatszamitas: ");
  opt->addUsage(
    "\t\t -s  (--stac) <halofile>.xml           Definicios file, kotelezo");
  opt->addUsage(
    "\t\t -i  (--ini) <resfile>.xml             Inicializacios file, nem "
    "kotelezo");
  opt->addUsage(" ");
  opt->addUsage("\t tartozkodasi ido szamitasa: ");
  opt->addUsage(
    "\t\t -t  (--travel_time) <halofile>.xml    Tartozkodasi ido szamitas a "
    "halofile.xml feladaton");
  opt->addUsage(" ");
  opt->addUsage("\t koncentracioeloszlas szamitasa: ");
  opt->addUsage(
    "\t\t -c  (--conc_transp) <halofile>.xml    Koncentracio eloszlas "
    "szamitasa a halofile.xml feladaton");
  opt->addUsage(" ");
  opt->addUsage("\t parameter megvaltoztatasa: ");
  opt->addUsage(
    "\t\t -m  (--mod_prop) <halofile_regi>.spr -e (--element_ID) <ID_ag/csp> "
    "-p (--property_ID) <ID_adat> -n (--newValue) <uj_ertek> -o (--outfile) "
    "<halofajl_uj>.spr    Adatmodositas: halofajl_regi.xml -> "
    "halofajl_uj.xml");
  opt->addUsage(" ");
  opt->addUsage("\t minden elem listazasa a kepernyora: ");
  opt->addUsage("\t\t -l  (--list_all_elements) <halofile_regi>.spr");
  opt->addUsage(" ");
  opt->addUsage("\t adat kiolvasasa: ");
  opt->addUsage(
    "\t\t -g  (--get_data) <halofile_regi>.spr -e (--element_ID) <ID_ag/csp> "
    "-p (--property_ID) <ID_adat>");
  opt->addUsage(" ");
  opt->addUsage("\t erzekenysegvizsgalat: ");
  opt->addUsage(
    "\t\t -r  (--sensitivity) <halofile_regi>.spr -e (--element_ID) "
    "<ID_ag/csp> -p (--property_ID) <ID_adat>");
  opt->addUsage(" ");
  opt->addUsage("\t check for islands: ");
  opt->addUsage(
    "\t\t -x  (--export_to_connectivity_check) <halofile_regi>.spr");

  opt->addUsage(" ");
  opt->addUsage(" ");

  opt->setOption("stac", 's');
  opt->setOption("ini", 'i');
  opt->setOption("travel_time", 't');
  opt->setOption("conc_transp", 'c');
  opt->setOption("mod_prop", 'm');
  opt->setOption("element_ID", 'e');
  opt->setOption("property_ID", 'p');
  opt->setOption("newValue", 'n');
  opt->setOption("outfile", 'o');
  opt->setOption("list_all_elements", 'l');
  opt->setOption("get_data", 'g');
  opt->setOption("sensitivity", 'r');
  opt->setOption("export_for_connectivity_check", 'x');

  opt->processCommandArgs(argc, argv);

  if (!opt->hasOptions()) {
    mode = -1;
    opt->printUsage();
  } else {
    if (opt->getFlag("help") || opt->getFlag('h')) opt->printUsage();
    if (opt->getValue('s') != NULL || opt->getValue("stac") != NULL) {
      mode = 0;
      def_file = opt->getValue('s');
      if (debug_level > 0) {
        string msg = "\n Steady-state hydraulic simulation, data file: " +
                     def_file + "\n";
        if (debug_level > 0) {
          // cout << msg;
          logfile_write(msg, 1);
        }
      }
    }

    if (opt->getValue('i') != NULL || opt->getValue("ini") != NULL) {
      van_ini = true;
      mode = 0;
      ini_file = opt->getValue('i');
      if (debug_level > 0) {
        string msg = "\n Stacioner halozatszamitas, inicializacios file : " +
                     ini_file + "\n";
        if (debug_level > 0) {
          // cout << msg;
          logfile_write(msg, 1);
        }
      }
    }

    if (opt->getValue('t') != NULL || opt->getValue("travel_time") != NULL) {
      mode = 1;
      def_file = opt->getValue('t');
      cout << "Tartozkodasi ido szamitas, a feladatot tartalmazo file: "
           << def_file << endl;
    }

    if (opt->getValue('c') != NULL || opt->getValue("conc_transp") != NULL) {
      mode = 2;
      def_file = opt->getValue('c');
      cout << "Koncentracioeloszlas szamitas, a feladatot tartalmazo file: "
           << def_file << endl;
    }

    if (opt->getValue('r') != NULL || opt->getValue("sens") != NULL) {
      mode = 0;
      def_file = opt->getValue('r');
      cout << endl
           << "Steady-state hydraulic simulation and sensitivity analysis, "
           "data file: "
           << def_file << endl;
    }

    // ADATMODOSITAS
    //-----------------------------
    if (opt->getValue('m') != NULL || opt->getValue("mod_prop") != NULL) {
      mode = 3;
      def_file = opt->getValue('m');
      cout << endl << "Adatmodositas" << endl;
      /*cout << "\n regi adatfajl: "<<def_file;*/
    }
    if (opt->getValue('e') != NULL || opt->getValue("element_ID") != NULL) {
      element_ID = opt->getValue('e');
      /* cout << "\n element_ID: "<< element_ID ;*/
    }
    if (opt->getValue('p') != NULL || opt->getValue("property_ID") != NULL) {
      property_ID = opt->getValue('p');
      /* cout << "\n property_ID: "<< property_ID ;*/
    }
    if (opt->getValue('n') != NULL || opt->getValue("newValue") != NULL) {
      newValue = atof(opt->getValue('n'));
      /* cout << "\n newValue: "<< newValue;*/
    }
    if (opt->getValue('o') != NULL || opt->getValue("outfile") != NULL) {
      new_def_file = opt->getValue('o');
      /* cout << "\n new_def_file: "<< new_def_file<<endl;*/
    }

    // Minden elem listazasa
    //-----------------------------
    if (opt->getValue('l') != NULL ||
        opt->getValue("list_all_elements") != NULL) {
      mode = 4;
      def_file = opt->getValue('l');
      /*cout << endl << "Adatkiolvasas"<< endl ;*/
    }
    // element_ID es property_ID fent

    // Adat kiolvasasa
    //-----------------------------
    if (opt->getValue('g') != NULL || opt->getValue("get_data") != NULL) {
      mode = 5;
      def_file = opt->getValue('g');
      /*cout << endl << "Adatkiolvasas"<< endl ;*/
    }
    // element_ID es property_ID fent

    // Erzekenysegvizsgalat
    //-----------------------------
    if (opt->getValue('r') != NULL || opt->getValue("sensitivity") != NULL) {
      mode = 6;
      def_file = opt->getValue('r');
    }
    // element_ID es property_ID fentAdatkiolvasas"<< endl ;*/

    // Export for connectivity check
    //-----------------------------
    if (opt->getValue('x') != NULL || opt->getValue("export_for_connectivity_check") != NULL) {
      mode = 7;
      def_file = opt->getValue('x');
    }
    // element_ID es property_ID fentAdatkiolvasas"<< endl ;*/
  }
}

//--------------------------------------------------------------
double Staci::m_get_dprop() {
  bool elem_megvan = false;
  bool prop_megvan = false;
  double outdata = 0.0;

  for (unsigned int i = 0; i < agelemek.size(); i++) {
    if ((agelemek[i]->Get_nev()) == element_ID) {
      /*cout<<endl<<"\t Megvan az agelem: "<<element_ID<<endl;*/
      elem_megvan = true;
      if ((property_ID == "diameter") || (property_ID == "mass_flow_rate") ||
          (property_ID == "bottom_level") || (property_ID == "water_level") ||
          (property_ID == "position")) {
        // prop_megvan = true;
        outdata = agelemek[i]->Get_dprop(property_ID);
        prop_megvan = true;
      }
    }
  }

  for (unsigned int i = 0; i < cspok.size(); i++) {
    if ((cspok[i]->Get_nev()) == element_ID) {
      /*cout<<endl<<"\t Megvan az agelem: "<<element_ID<<endl;*/
      elem_megvan = true;
      if ((property_ID == "pressure") || (property_ID == "demand") ||
          (property_ID == "head")) {
        // prop_megvan = true;
        outdata = cspok[i]->Get_dprop(property_ID);
        prop_megvan = true;
      }
    }
  }

  if (!elem_megvan) {
    cout << endl
         << endl
         << "ERROR!!! Staci::m_get_dprop(): node/edge not found: " << element_ID
         << endl
         << endl;
    exit(-1);
  } else {
    if (!prop_megvan) {
      cout << endl
           << endl
           << "ERROR!!! Staci::m_get_dprop(): node/edge: " << element_ID
           << ", property not found: " << property_ID << endl;
      cout << "  edges: diameter|mass_flow_rate|bottom_level|water_level|position" << endl;
      cout << "  nodes: pressure|head|demand" << endl << endl;
      exit(-1);
    } else
      return outdata;
  }
}

//--------------------------------------------------------------
void Staci::m_set_dprop() {
  bool elem_megvan = false;
  bool prop_megvan = false;
  // double outdata;

  for (unsigned int i = 0; i < agelemek.size(); i++) {
    if ((agelemek[i]->Get_nev()) == element_ID) {
      /*cout<<endl<<"\t Megvan az agelem: "<<element_ID<<endl;*/
      elem_megvan = true;
      if (property_ID == "diameter") {
        agelemek[i]->Set_dprop(property_ID, newValue);
        prop_megvan = true;
      }

      if (property_ID == "bottom_level") {
        agelemek[i]->Set_dprop(property_ID, newValue);
        prop_megvan = true;
      }

      if (property_ID == "water_level") {
        agelemek[i]->Set_dprop(property_ID, newValue);
        prop_megvan = true;
      }

      if (property_ID == "position") {
        agelemek[i]->Set_dprop(property_ID, newValue);
        prop_megvan = true;
      }
    }
  }

  for (unsigned int i = 0; i < cspok.size(); i++) {
    // cout<<endl<<"\t Csomopont: "<<element_ID;
    if ((cspok[i]->Get_nev()) == element_ID) {
      // cout<<endl<<"\t Megvan a csomopont: "<<element_ID<<endl;

      elem_megvan = true;
      if (property_ID == "demand") {
        cspok[i]->Set_dprop(property_ID, newValue);
        prop_megvan = true;
      }
    }
  }

  if (!elem_megvan) {
    cout << endl
         << endl
         << "HIBA!!! Staci::m_set_dprop(): Nincs ilyen elem: " << element_ID
         << endl
         << endl;
    exit(-1);
  }
  if (!prop_megvan) {
    cout << endl
         << endl
         << "HIBA!!! Staci::m_set_dprop(): Elem: " << element_ID
         << ", nincs ilyen adat: " << property_ID << endl
         << endl;
    exit(-1);
  }
}

//--------------------------------------------------------------
double Staci::get_sum_of_consumption() {
  double sum = 0.0;
  for (unsigned int i = 0; i < cspok.size(); i++) {
    sum += fabs(cspok.at(i)->Get_fogy());
  }
  return sum;
}

//--------------------------------------------------------------
double Staci::get_sum_of_pos_consumption() {
  double sum = 0.0, cons;
  for (unsigned int i = 0; i < cspok.size(); i++) {
    cons = cspok.at(i)->Get_fogy();
    if (cons > 0.) sum += cons;
  }
  return sum;
}

//--------------------------------------------------------------
double Staci::get_sum_of_neg_consumption() {
  double sum = 0.0, cons;
  for (unsigned int i = 0; i < cspok.size(); i++) {
    cons = cspok.at(i)->Get_fogy();
    if (cons < 0.) sum += fabs(cons);
  }
  return sum;
}

//--------------------------------------------------------------
double Staci::get_dprop(string in_element_ID, string in_property_ID) {
  bool elem_megvan = false;
  bool prop_megvan = false;
  double outdata = 0.0;

  for (unsigned int i = 0; i < agelemek.size(); i++) {
    if ((agelemek[i]->Get_nev()) == in_element_ID) {
      /*cout<<endl<<"\t Megvan az agelem: "<<in_element_ID<<endl;*/
      elem_megvan = true;
      if ((in_property_ID == "diameter") ||
          (in_property_ID == "mass_flow_rate") ||
          (in_property_ID == "bottom_level") ||
          (in_property_ID == "water_level") || (in_property_ID == "position")) {
        // prop_megvan = true;
        outdata = agelemek[i]->Get_dprop(in_property_ID);
        prop_megvan = true;
      }
    }
  }

  for (unsigned int i = 0; i < cspok.size(); i++) {
    if ((cspok[i]->Get_nev()) == in_element_ID) {
      /*cout<<endl<<"\t Megvan az agelem: "<<in_element_ID<<endl;*/
      elem_megvan = true;
      if ((in_property_ID == "pressure") || (in_property_ID == "demand") ||
          (in_property_ID == "head")) {
        // prop_megvan = true;
        outdata = cspok[i]->Get_dprop(in_property_ID);
        prop_megvan = true;
      }
    }
  }

  if (!elem_megvan) {
    cout << endl
         << endl
         << "HIBA!!! Staci::m_get_dprop(): Nincs ilyen elem: " << in_element_ID
         << endl
         << endl;
    exit(-1);
  } else {
    if (!prop_megvan) {
      cout << endl
           << endl
           << "HIBA!!! Staci::m_get_dprop(): Elem: " << in_element_ID
           << ", nincs ilyen adat: " << in_property_ID << endl
           << endl;
      exit(-1);
    } else
      return outdata;
  }
}

//--------------------------------------------------------------
void Staci::set_dprop(string in_element_ID, string in_property_ID,
                      double in_newValue) {
  bool elem_megvan = false;
  bool prop_megvan = false;

  for (unsigned int i = 0; i < agelemek.size(); i++) {
    if ((agelemek[i]->Get_nev()) == in_element_ID) {
      /*cout<<endl<<"\t Megvan az agelem: "<<in_element_ID<<endl;*/
      elem_megvan = true;
      if (in_property_ID == "diameter") {
        agelemek[i]->Set_dprop(in_property_ID, in_newValue);
        prop_megvan = true;
      }

      if (in_property_ID == "bottom_level") {
        agelemek[i]->Set_dprop(in_property_ID, in_newValue);
        prop_megvan = true;
      }

      if (in_property_ID == "water_level") {
        agelemek[i]->Set_dprop(in_property_ID, in_newValue);
        prop_megvan = true;
      }

      if (in_property_ID == "position") {
        agelemek[i]->Set_dprop(in_property_ID, in_newValue);
        prop_megvan = true;
      }
    }
  }

  for (unsigned int i = 0; i < cspok.size(); i++) {
    // cout<<endl<<"\t Csomopont: "<<in_element_ID;
    if ((cspok[i]->Get_nev()) == in_element_ID) {
      // cout<<endl<<"\t Megvan a csomopont: "<<in_element_ID<<endl;

      elem_megvan = true;
      if (in_property_ID == "demand") {
        cspok[i]->Set_dprop(in_property_ID, in_newValue);
        prop_megvan = true;
      }
    }
  }

  if (!elem_megvan) {
    cout << endl
         << endl
         << "HIBA!!! Staci::m_set_dprop(): Nincs ilyen elem: " << in_element_ID
         << endl
         << endl;
    exit(-1);
  }
  if (!prop_megvan) {
    cout << endl
         << endl
         << "HIBA!!! Staci::m_set_dprop(): Elem: " << in_element_ID
         << ", nincs ilyen adat: " << in_property_ID << endl
         << endl;
    exit(-1);
  }
}

//--------------------------------------------------------------
void Staci::build_system() {
  ostringstream msg1;

  // for (int i = 0; i < agelemek.size(); i++)
  // agelemek.at(i)->Info();
  // cin.get();

  bool stop = false;

  logfile_write("\n Azonos ID-k keresese....", 3);
  // ELEMEK
  string nev1, nev2;
  // int szam = 0;
  for (unsigned int i = 0; i < agelemek.size(); i++) {
    /*printf("\n i=%d",i); cin.get();*/
    // szam = 0;
    nev1 = agelemek.at(i)->Get_nev();
    /*cout<<endl<<nev1;
    cin.get();*/
    for (unsigned int j = 0; j < agelemek.size(); j++) {
      nev2 = agelemek.at(j)->Get_nev();
      // cout<<"\n\t"<<nev1 <<"?=?"<<nev2;
      if (i != j) {
        if (nev1 == nev2) {
          ostringstream msg;
          msg << "\n !!!ERROR!!! edge #" << i << ": " << nev1 << " and edge #" << j << ": " << nev2 << " with same ID!" << endl;
          // cout << msg.str();
          logfile_write(msg.str(), 1);
          stop = true;
        }
      }
    }

    for (unsigned int j = 0; j < cspok.size(); j++) {
      nev2 = cspok.at(j)->Get_nev();
      // cout<<"\n\t"<<nev1 <<"?=?"<<nev2;
      if (i != j) {
        if (nev1 == nev2) {
          ostringstream msg;
          msg << "\n !!!ERROR!!! edge #" << i << ": " << nev1 << " and node #" << j << ": " << nev2 << " with same ID!" << endl;
          // cout << msg.str();
          logfile_write(msg.str(), 1);
          stop = true;
        }
      }
    }
  }

  if (stop)
    exit(-1);
  else
    logfile_write("\t ok.", 3);

  logfile_write("\n\n Rendszer epitese...", 3);
  bool e_megvan = false;
  bool v_megvan = false;
  unsigned int j = 0;
  int cspe = -1, cspv = -1;
  ostringstream strstrm;

  // az vege csp. nem mindig kell...
  for (unsigned int i = 0; i < agelemek.size(); i++) {
    /*cout << "\n\t" << agelemek[i]->Get_nev() << ":\tcspe:";*/

    e_megvan = false;
    j = 0;
    while ((j < cspok.size()) && (!e_megvan)) {
      // log
      strstrm.str("");
      strstrm << "\n\t" << agelemek[i]->Get_nev().c_str()
              << " cspe: " << agelemek[i]->Get_Cspe_Nev().c_str() << " =? "
              << cspok[j]->Get_nev().c_str();
      logfile_write(strstrm.str(), 5);
      // cout<<strstrm.str();
      // cout << strstrm.str();
      if ((agelemek[i]->Get_Cspe_Nev()).compare(cspok[j]->Get_nev()) == 0) {
        e_megvan = true;
        cspe = j;
        cspok.at(j)->Set_index(j); // WR indexing the nodes
        cspok[j]->ag_ki.push_back(i);
        // log
        //                logfile_write(" OK", 3);
        //                 strstrm.str("");
        //                 strstrm << "\n\t" << agelemek[i]->Get_nev() << "
        //                 cspe: "
        //                         << agelemek[i]->Get_Cspe_Nev() << " OK ";
        //                 logfile_write(strstrm.str(), 4);
        //                cout<<strstrm.str();
      }
      j++;
    }
    if (!e_megvan) {
      strstrm.str("");
      strstrm << "\n!!! Nincs meg a " << agelemek[i]->Get_nev().c_str()
              << " agelem eleji csomopont: " << agelemek[i]->Get_Cspe_Nev()
              << " !!!";
      logfile_write(strstrm.str(), 1);
      // cout << strstrm.str();
      StaciException csphiba(strstrm.str());
      throw csphiba;
    } else {
      //        strstrm.str("");
      //        strstrm<<"\n\t"<<agelemek[i]->Get_nev()<<" cspe: "<<cspe;
      //        cout<<strstrm.str();
    }

    // cout<<"\tcspv: ";
    if (agelemek[i]->Get_Csp_db() == 2) {
      v_megvan = false;
      j = 0;
      while ((j < cspok.size()) && (!v_megvan)) {
        // log
        //                      strstrm.str("");
        //                      strstrm << "\n\t"<<agelemek[i]->Get_nev()<<"
        //                      cspv: "
        //                          <<agelemek[i]->Get_Cspv_Nev()<<" =?
        //                          "<<cspok[j]->Get_nev();
        //                      logfile_write(strstrm.str(), 3);
        //                      cout<<strstrm.str();
        if ((agelemek[i]->Get_Cspv_Nev()).compare(cspok[j]->Get_nev()) == 0) {
          v_megvan = true;
          cspv = j;
          cspok.at(j)->Set_index(j); // WR: indexing the nodes
          cspok[j]->ag_be.push_back(i);
          // log
          //                                logfile_write(" OK", 3);
          //                                strstrm.str("");
          //                                strstrm<<"\n\t"<<agelemek[i]->Get_nev()<<"
          //                                cspv: "
          //                                   <<agelemek[i]->Get_Cspv_Nev()<<"
          //                                   OK ";
          //                                logfile_write(strstrm.str(), 3);
          ////                                cout<<strstrm.str();
        }
        j++;
      }
      if (!v_megvan) {
        strstrm.str("");
        strstrm << "\n!!! Nincs meg a " << agelemek[i]->Get_nev().c_str()
                << " agelem vegi csomopont!";
        logfile_write(strstrm.str(), 1);
        // cout << strstrm.str();
      } else {
        //      strstrm.str("");
        //      strstrm<<"\n\t"<<agelemek[i]->Get_nev()<<" cspv: "<<cspv;
        //      cout<<strstrm.str();
      }
    } else {
      strstrm.str("");
      strstrm << "\n\t" << agelemek[i]->Get_nev().c_str() << " cspv: - "
              << cspv;
      //        cout<<strstrm.str();
    }

    if (agelemek[i]->Get_Csp_db() == 2) {
      agelemek[i]->add_csp(cspe, cspv);
      Add_edge(cspe, cspv); // WR: collecting the edge vector for igraph
    }
    else
      agelemek[i]->add_csp(cspe, -1);
  }

  // Surlodas beallitasa
  for (unsigned int i = 0; i < agelemek.size(); i++) {
    if (agelemek[i]->Get_Tipus() == "Cso") {
      /*cout << endl << "friction_model:" << friction_model;
      cin.get();*/
      agelemek[i]->Set_friction_model(friction_model);
    }
  }


  logfile_write("\t ok.", 3);


  // for (unsigned int i = 0; i < cspok.size(); i++) {
  //   double d = cspok.at(i)->Get_dprop("demand");
  //   cspok.at(i)->Set_dprop("demand", d * 0.0001);
  // }
}

//--------------------------------------------------------------
void Staci::build_system_old() {
  cout << "\nRendszer epitese...";
  logfile_write("\n Rendszer epitese\n---------------------------", 2);
  bool e_megvan = false;
  bool v_megvan = false;
  unsigned int j = 0;
  int cspe = 0, cspv = 0;
  ostringstream strstrm;
  // az 1. csp nem mindig kell...
  for (unsigned int i = 0; i < agelemek.size(); i++) {
    // cout<<"\n\t"<<agelemek[i]->Get_nev()<<":\tcspe:";
    if (agelemek[i]->Get_Csp_db() == 2) {
      e_megvan = false;
      j = 0;
      while ((j < cspok.size()) && (!e_megvan)) {
        // log
        strstrm.str("");
        strstrm << "\n\t" << agelemek[i]->Get_nev().c_str()
                << " cspe: " << agelemek[i]->Get_Cspe_Nev() << " =? "
                << cspok[j]->Get_nev();
        logfile_write(strstrm.str(), 3);
        if ((agelemek[i]->Get_Cspe_Nev()).compare(cspok[j]->Get_nev()) == 0) {
          e_megvan = true;
          cspe = j;
          cspok[j]->ag_ki.push_back(i);
          // log
          logfile_write(" OK", 3);
          strstrm.str("");
          strstrm << "\n\t" << agelemek[i]->Get_nev().c_str()
                  << " cspe: " << agelemek[i]->Get_Cspe_Nev() << " OK ";
          logfile_write(strstrm.str(), 2);
        }
        j++;
      }
      if (!e_megvan) {
        strstrm.str("");
        strstrm << "\n!!! Nincs meg a " << agelemek[i]->Get_nev().c_str()
                << " agelem eleji csomopont!";
        logfile_write(strstrm.str(), 1);
        cout << " nincs meg!!!";
        StaciException csphiba(strstrm.str());
        throw csphiba;
      } else {
        strstrm.str("");
        strstrm << "\n\t" << agelemek[i]->Get_nev().c_str()
                << " cspe: " << cspe;
        // cout<<" ok";
      }
    } else {
      strstrm.str("");
      strstrm << "\n\t" << agelemek[i]->Get_nev().c_str() << " cspe: - "
              << cspe;
      // cout<<" - ";
    }

    // cout<<"\tcspv: ";
    v_megvan = false;
    j = 0;
    while ((j < cspok.size()) && (!v_megvan)) {
      // log
      strstrm.str("");
      strstrm << "\n\t" << agelemek[i]->Get_nev().c_str()
              << " cspv: " << agelemek[i]->Get_Cspv_Nev() << " =? "
              << cspok[j]->Get_nev();
      logfile_write(strstrm.str(), 3);
      if ((agelemek[i]->Get_Cspv_Nev()).compare(cspok[j]->Get_nev()) == 0) {
        v_megvan = true;
        cspv = j;
        cspok[j]->ag_be.push_back(i);
        // log
        logfile_write(" OK", 3);
        strstrm.str("");
        strstrm << "\n\t" << agelemek[i]->Get_nev().c_str()
                << " cspv: " << agelemek[i]->Get_Cspv_Nev() << " OK ";
        logfile_write(strstrm.str(), 2);
      }
      j++;
    }
    if (!v_megvan) {
      strstrm.str("");
      strstrm << "\n!!! Nincs meg a " << agelemek[i]->Get_nev().c_str()
              << " agelem vegi csomopont!";
      logfile_write(strstrm.str(), 1);
      // cout<<" nincs meg!!!";
    } else {
      strstrm.str("");
      strstrm << "\n\t" << agelemek[i]->Get_nev().c_str() << " cspv: " << cspv;
      // cout<<" ok";
    }
    if (agelemek[i]->Get_Csp_db() == 2)
      agelemek[i]->add_csp(cspe, cspv);
    else
      // agelemek[i]->add_csp(-1, cspv);
      agelemek[i]->add_csp(cspv, -1);
  }

  /*for (unsigned int i=0; i<cspok.size(); i++)
    cout<<cspok.at(i)->Info(false);
    for (unsigned int i=0; i<agelemek.size(); i++)
    cout<<agelemek.at(i)->Info();
  */
}

//--------------------------------------------------------------
void Staci::logfile_write(string msg, int msg_debug_level) {
  if (debug_level >= msg_debug_level) {
    ofstream outfile(out_file.c_str(), ios::app);
    outfile << msg;
    outfile.close();
    cout << msg;
  }
}

//--------------------------------------------------------------
void Staci::list_system() {
  logfile_write("\n\n Csomopontok:\n--------------------------", 3);
  for (unsigned int i = 0; i < cspok.size(); i++)
    logfile_write(cspok[i]->Info(true), 3);
  logfile_write("\n\n Agelemek:\n--------------------------", 3);
  for (unsigned int i = 0; i < agelemek.size(); i++)
    logfile_write(agelemek[i]->Info(), 3);
}

//--------------------------------------------------------------
string Staci::list_results() {
  ostringstream strstrm;

  // string pot;
  // const
  // int string_hossz = 10;

  strstrm << scientific << setprecision(3) << showpos;
  strstrm << endl << endl << "EREDMENYEK:";
  for (unsigned int i = 0; i < agelemek.size(); i++)
    strstrm << endl
            << "\t" << agelemek[i]->Get_nev()
            << ":\tmp=" << agelemek[i]->Get_mp() << " kg/s"
            << "\tQ=" << (3600 * (agelemek[i]->Get_Q())) << " m3/h"
            << "\tv=" << agelemek[i]->Get_v() << " m/s";
  strstrm
      << endl
      << "\t"
      << "-----------------------------------------------------------------";
  for (unsigned int i = 0; i < cspok.size(); i++)
    strstrm << endl
            << "\t" << cspok[i]->Get_nev()
            << ":\t p=" << cspok[i]->Get_p() * 1000 * 9.81 / 1e5 << " bar"
            << "\tH=" << cspok[i]->Get_p() << " m"
            << ",   H+magassag=" << cspok[i]->Get_p() + cspok[i]->Get_h()
            << " m";
  strstrm << endl << endl;
  return strstrm.str();
}

//--------------------------------------------------------------
void Staci::save_results(bool conv_reached) {
  if (do_save_file) {
    data_io datta_io(res_file.c_str());
    datta_io.save_results(FolyMenny, get_sum_of_neg_consumption()*3.600, get_sum_of_pos_consumption()*3.600, cspok, agelemek, conv_reached, debug_level);
  } else {
    cout << "\n WARNING!\n!!! do_save_file = false, results will not be saved "
         "to the .spr file !!!!\n\n";
  }
}

//--------------------------------------------------------------
void Staci::save_mod_prop(bool is_general_property) {
  data_io datta_io(res_file.c_str());
  datta_io.save_mod_prop(cspok, agelemek, element_ID, property_ID,
                         is_general_property);
}

//--------------------------------------------------------------
void Staci::save_mod_prop_all_elements(string property_ID) {
  data_io datta_io(res_file.c_str());
  datta_io.save_mod_prop_all_elements(cspok, agelemek, property_ID);
}

//--------------------------------------------------------------
void Staci::build_vectors(Vec_DP &x, Vec_DP &f, bool create_sparse_pattern) {
  unsigned long int N = cspok.size() + agelemek.size();
  Vec_DP col(N);
  m_jac.clear();
  m_jac.reserve(N);
  vector<double> jac_row(N, 0.0);

  for (int i = 0; i < N; i++) m_jac.push_back(jac_row);

  if (create_sparse_pattern) {
    m_is_element_empty.clear();
    vector<bool> tmp(N, true);
    for (int i = 0; i < N; i++) m_is_element_empty.push_back(tmp);
  }

  // A matrix mintazatanak mente

  // aktualis x kiszedes az elemekbol
  //---------------------------------------------------
  for (unsigned int i = 0; i < agelemek.size(); i++) {
    x[i] = agelemek[i]->Get_mp();
    //      cout<<"\n x="<<x[i];
  }
  for (unsigned int i = 0; i < cspok.size(); i++) {
    x[agelemek.size() + i] = cspok[i]->Get_p();
    //      cout<<"\n x="<<x[agelemek.size() + i];
  }
  // f es Jacobi kiertekelese az aktualis adatokkal
  int Q_indx, pe_indx, pv_indx;
  vector<double> pevhev(4, 0.0);
  vector<double> jv;

  /*for (int i = 0; i < N; i++)
      for (int j = 0; j < N; j++)
          m_jac[i][j] = 0.0;*/
  for (unsigned int i = 0; i < agelemek.size(); i++) {
    if (agelemek[i]->Get_Csp_db() == 1) {
      pevhev[0] = cspok[agelemek[i]->Get_Cspe_Index()]->Get_p();
      pevhev[1] = 0.0;
      pevhev[2] = cspok[agelemek[i]->Get_Cspe_Index()]->Get_h();
      pevhev[3] = 0.0;
    } else {
      pevhev[0] = cspok[agelemek[i]->Get_Cspe_Index()]->Get_p();
      pevhev[1] = cspok[agelemek[i]->Get_Cspv_Index()]->Get_p();
      pevhev[2] = cspok[agelemek[i]->Get_Cspe_Index()]->Get_h();
      pevhev[3] = cspok[agelemek[i]->Get_Cspv_Index()]->Get_h();
    }
    f[i] = agelemek[i]->f(pevhev);
    jv = agelemek[i]->df(pevhev);
    // cout<<endl<<jv;
    // pevhev.clear();

    Q_indx = i;
    pe_indx = agelemek[i]->Get_Cspe_Index();
    pv_indx = agelemek[i]->Get_Cspv_Index();

    m_jac[i][agelemek.size() + pe_indx] = jv.at(0);
    if (create_sparse_pattern) {
      m_is_element_empty[i][agelemek.size() + pe_indx] = false;
      m_nnz++;
    }

    if (agelemek[i]->Get_Csp_db() == 2) {
      m_jac[i][agelemek.size() + pv_indx] = jv.at(1);
      if (create_sparse_pattern) {
        m_is_element_empty[i][agelemek.size() + pv_indx] = false;
        m_nnz++;
      }
    }

    m_jac[i][Q_indx] = jv.at(2);
    if (create_sparse_pattern) {
      m_is_element_empty[i][Q_indx] = false;
      m_nnz++;
    }

    col[i] = jv.at(3);
    jv.clear();
  }

  for (unsigned int i = 0; i < cspok.size(); i++) {
    f[agelemek.size() + i] = -cspok[i]->Get_fogy();
    for (unsigned int j = 0; j < cspok[i]->ag_be.size(); j++) {
      f[agelemek.size() + i] += agelemek[cspok[i]->ag_be.at(j)]->Get_mp();
      m_jac[agelemek.size() + i][cspok[i]->ag_be.at(j)] = +1.0;
      if (create_sparse_pattern) {
        m_is_element_empty[agelemek.size() + i][cspok[i]->ag_be.at(j)] = false;
        m_nnz++;
      }
    }
    for (unsigned int j = 0; j < cspok[i]->ag_ki.size(); j++) {
      f[agelemek.size() + i] -= agelemek[cspok[i]->ag_ki.at(j)]->Get_mp();
      m_jac[agelemek.size() + i][cspok[i]->ag_ki.at(j)] = -1.0;
      if (create_sparse_pattern) {
        m_is_element_empty[agelemek.size() + i][cspok[i]->ag_ki.at(j)] = false;
        m_nnz++;
      }
    }
    col[agelemek.size() + i] = -cspok[i]->Get_fogy();
  }

  //    int NN = x.size();
  //    m_jac.clear();
  //    vector<double> tmp;
  //    for (int i = 0; i < NN; i++) {
  //        tmp.clear();
  //        for (int j = 0; j < NN; j++)
  //            //tmp.push_back(jac[i][j]);
  //            tmp.push_back(sin((i + 1) * (j + 1)));
  //        m_jac.push_back(tmp);
  //    }
  //
  //    ofstream myfile;
  //    myfile.open("matrix.csv");
  //    for (int i = 0; i < NN; i++) {
  //        for (int j = 0; j < NN; j++) {
  //            myfile << jac[i][j];
  //            if (j < NN - 1)
  //                myfile << ";";
  //        }
  //        myfile << endl;
  //    }
  //    myfile.close();
  //    cout << "\n matrix.csv kiirva.\n";
  //    cin.get();
  //
  //    //  for (int i = 0; i < NN; i++) {
  //    //      cout << endl;
  //    //      for (int j = 0; j < NN; j++)
  //    //          cout << " " << m_jac[i][j];
  //    //  }
  //    //  cout << endl;

  if (debug_level > 6) Print_Jacobian();

  // if (create_sparse_pattern) {
  //   ostringstream msg1;
  //   msg1 << endl
  //        << " Number of nonzero Jacobian entries: " << m_nnz << " out of "
  //        << (N * N);
  //   msg1 << " (" << ((((double)m_nnz) / N / N) * 100) << "%)" << endl;
  //   // msg1 << endl << " Jacobian size check : m_jac.capacity()=" <<
  //   // m_jac.capacity();
  //   // msg1 << endl << " Jacobian size check : m_jac[0].capacity()=" <<
  //   // m_jac[0].capacity()<<endl;
  //   logfile_write(msg1.str(), 1);
  // }
}

//--------------------------------------------------------------
void Staci::build_vectors_frozen_Jacobian(Vec_DP &x, Vec_DP &f) {
  // int N = cspok.size() + agelemek.size();
  vector<double> pevhev;

  for (unsigned int i = 0; i < agelemek.size(); i++)
    x[i] = agelemek[i]->Get_mp();
  for (unsigned int i = 0; i < cspok.size(); i++)
    x[agelemek.size() + i] = cspok[i]->Get_p();

  for (unsigned int i = 0; i < agelemek.size(); i++) {
    if (agelemek[i]->Get_Csp_db() == 1) {
      pevhev.push_back(cspok[agelemek[i]->Get_Cspe_Index()]->Get_p());
      pevhev.push_back(0);
      pevhev.push_back(cspok[agelemek[i]->Get_Cspe_Index()]->Get_h());
      pevhev.push_back(0);
    } else {
      pevhev.push_back(cspok[agelemek[i]->Get_Cspe_Index()]->Get_p());
      pevhev.push_back(cspok[agelemek[i]->Get_Cspv_Index()]->Get_p());
      pevhev.push_back(cspok[agelemek[i]->Get_Cspe_Index()]->Get_h());
      pevhev.push_back(cspok[agelemek[i]->Get_Cspv_Index()]->Get_h());
    }
    f[i] = agelemek[i]->f(pevhev);
    // jv = agelemek[i]->df(pevhev);
    pevhev.clear();
  }

  for (unsigned int i = 0; i < cspok.size(); i++) {
    f[agelemek.size() + i] = -cspok[i]->Get_fogy();
    for (unsigned int j = 0; j < cspok[i]->ag_be.size(); j++)
      f[agelemek.size() + i] += agelemek[cspok[i]->ag_be.at(j)]->Get_mp();

    for (unsigned int j = 0; j < cspok[i]->ag_ki.size(); j++)
      f[agelemek.size() + i] -= agelemek[cspok[i]->ag_ki.at(j)]->Get_mp();

    // col[agelemek.size() + i] = -cspok[i]->Get_fogy();
  }
}

//--------------------------------------------------------------
bool Staci::solve_system() {
  const int N = cspok.size() + agelemek.size();
  Mat_DP invjac(N, N);
  // Mat_DP jac(N, N), invjac(N, N);
  Vec_DP col(N), b(N), x(N), dx(N), f(N), xu(N);
  Vec_INT indx(N);
  int iter = 0;
  double e_mp = 1e10, e_p = 1e10, e_mp_r = 1e10, e_p_r = 1e10;
  bool konv_ok = false;

  // Progress file torlese
  if (debug_level > 0) {
    ofstream pfile(progress_file.c_str(), ios::trunc);
    pfile << fixed << setprecision(1) << 0.0 << "\n";
    pfile.close();
  }

  m_ss.str("");
  m_ss << "\n\nSolving system...\n====================================" << endl;

  if (debug_level > 0) {
    logfile_write(m_ss.str(), 1);
    // cout << m_ss.str();
  }

  // Solver
  build_vectors(x, f, true);

  // Iteracio!!!
  bool comp_ok = true;
  // cout<<endl<<"debug_level:"<<debug_level<<endl; cin.get();
  while ((iter < iter_max + 1) && (!konv_ok)) {
    if (debug_level > 0) progress_file_write((double)iter / iter_max * 100.0);

    if ((e_mp > 0.1 || e_p > 0.1) || (iter % 5 == 0))
      build_vectors(x, f, false);
    else {
      build_vectors_frozen_Jacobian(x, f);
    }

    compute_error(f, e_mp, e_p, e_mp_r, e_p_r, konv_ok);

    logfile_write(iter_info(x, f, iter, e_mp, e_p), 1);

    print_worst_iter(x, f, 3);

    update_relax(e_mp, e_p, e_mp_r, e_p_r);

    comp_ok = umfpack_solver(x, f);

    if ((e_mp < e_mp_max) && (e_p < e_p_max))
      konv_ok = true;

    m_ss.str("");
    for (unsigned int i = 0; i < agelemek.size(); i++)
      m_ss << endl << "\t" << agelemek.at(i)->Get_nev() << ": \tmp=" << x[i] << ", \tf=" << f[i];
    for (unsigned int i = 0; i < cspok.size(); i++)
      m_ss << endl << "\t" << cspok.at(i)->Get_nev() << ": \tp =" << x[agelemek.size() + i] << ", \tf=" << f[agelemek.size() + i];
    logfile_write(m_ss.str(), 4);

    if ((!comp_ok) && (iter > 100)) {
      cout << endl
           << endl
           << "WARNING: Staci::solve_system() -> umfpack_solver did not "
           "provide a solution, switching back to nr_solver!!\n\n";
      nr_solver(x, f);
    } else
      iter++;
  }

  if (!konv_ok)
    print_worst_iter(x, f, 1);

  if (konv_ok) {
    for (unsigned int i = 0; i < agelemek.size(); i++) {
      if (0 == strcmp( (agelemek.at(i)->GetType()).c_str(), "Csatorna")) {
        agelemek.at(i)->build_res();
      }
    }
  }

  Set_FolyTerf();

  return konv_ok;
}

void Staci::print_worst_iter(const Vec_DP x, const Vec_DP f , const int a_debug_level) {

  vector<val_and_ID> v_edges;
  vector<val_and_ID> v_nodes;
  val_and_ID s;

  for (int i = 0 ; i < agelemek.size(); i++) {
    s.val =  fabs(f[i]);
    s.x =  x[i];
    s.ID = agelemek.at(i)->Get_nev().c_str();
    v_edges.push_back(s);
  }

  for (int i = 0 ; i < cspok.size(); i++) {
    s.val =  fabs(f[i + agelemek.size()]);
    s.x =  x[i + agelemek.size()];
    s.ID = cspok.at(i)->Get_nev().c_str();
    v_nodes.push_back(s);
  }

  sort(v_edges.begin(), v_edges.end(), comparison_function1);
  sort(v_nodes.begin(), v_nodes.end(), comparison_function1);

  int N_to_print = 10;
  if (N_to_print > agelemek.size())
    N_to_print = agelemek.size();
  if (N_to_print > cspok.size())
    N_to_print = cspok.size();

  m_ss << "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n";
  m_ss << "The " << N_to_print << " worst edges:";
  for (int i = 0; i < N_to_print; i++)
    m_ss << "\n\t " << v_edges.at(i).ID.c_str() << ": mp=" << v_edges.at(i).x << ", f=" << v_edges.at(i).val;
  m_ss << "\nThe " << N_to_print << " worst nodes:";
  for (int i = 0; i < N_to_print; i++)
    m_ss << "\n\t " << v_nodes.at(i).ID.c_str() << ": p =" << v_nodes.at(i).x << ", f=" << v_nodes.at(i).val;
  m_ss << "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";

  // if (!konv_ok) m_ss << "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n";
  // for (unsigned int i = 0; i < agelemek.size(); i++)
  //   m_ss << "\n\t" << agelemek[i]->Get_nev().c_str() << "("
  //        << agelemek[i]->GetType().c_str() << "): mp=" << x[i]
  //        << ", f=" << f[i];

  // for (unsigned int i = 0; i < cspok.size(); i++)
  //   m_ss << "\n\t" << cspok[i]->Get_nev().c_str()
  //        << ": p=" << x[agelemek.size() + i]
  //        << ", f=" << f[agelemek.size() + i];
  // if (!konv_ok)
  //   m_ss << "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\nBajjj van.\n\n";

  logfile_write(m_ss.str(), a_debug_level);

}

//--------------------------------------------------------------
bool Staci::solve_system_old() {
  typedef enum { Newton_Raphson, Linear } megoldo;
  const int N = cspok.size() + agelemek.size();
  Mat_DP jac(N, N), invjac(N, N);
  Vec_DP col(N), b(N), x(N), dx(N), f(N), xu(N);
  Vec_INT indx(N);
  DP d;
  int iter = 0;
  double e_mp = 1e10, e_p = 1e10, e_mp_r = 1e10, e_p_r = 1e10;
  double RELAX_MAX = 1.0;
  bool konv_ok = false;
  ostringstream strstrm;
  // double RELAX=0.1;
  // double szorzo=1.1;
  relax = 1;
  int akt_megoldo;
  akt_megoldo = (megoldo)Newton_Raphson;

  cout << scientific << setprecision(3) << showpos;

  // Inicializalas
  // for (int i=0; i<agelemek.size(); i++)agelemek[i]->Ini(0,0);
  // for (int i=0; i<cspok.size(); i++)cspok[i]->Ini(0,0);

  strstrm.str("");
  strstrm << "\n\nSzamitas\n====================================" << endl;
  strstrm << "iter_max=" << iter_max << endl;
  strstrm << "e_mp_max=" << e_mp_max << endl;
  strstrm << "e_p_max =" << e_p_max << endl;
  strstrm << "megoldo :" << akt_megoldo << endl;
  logfile_write(strstrm.str(), 1);
  cout << endl << "SZAMITAS..." << endl;

  // Iteracio!!!
  while ((iter < iter_max) && (!konv_ok)) {
    cout << endl << "debug_level:" << debug_level << endl;
    cin.get();
    if (debug_level > 0) progress_file_write((double)iter / iter_max * 100.0);

    strstrm.str("");
    strstrm.setf(ios::dec);
    strstrm.unsetf(ios::showpos);
    strstrm << endl
            << "  " << iter << "./" << iter_max << " iteracio:  RELAX=" << relax
            << "   ";
    logfile_write(strstrm.str(), 1);

    // aktualis x kiszedes az elemekbol
    //-------------------------------------------------------------------------------------------
    for (unsigned int i = 0; i < agelemek.size(); i++)
      x[i] = agelemek[i]->Get_mp();
    for (unsigned int i = 0; i < cspok.size(); i++)
      x[agelemek.size() + i] = cspok[i]->Get_p();

    // f es Jacobi kiertekelese az aktualis adatokkal
    int Q_indx, pe_indx, pv_indx;
    vector<double> pevhev(4, 0.0);
    vector<double> jv;

    for (int i = 0; i < N; i++)
      for (int j = 0; j < N; j++) jac[i][j] = 0.0;

    for (unsigned int i = 0; i < agelemek.size(); i++) {
      if (agelemek[i]->Get_Csp_db() == 1) {
        /*pevhev.push_back(0);
          pevhev.push_back(cspok[agelemek[i]->Get_Cspv_Index()]->Get_p());
          pevhev.push_back(0);
          pevhev.push_back(cspok[agelemek[i]->Get_Cspv_Index()]->Get_h());*/
        pevhev[0] = cspok[agelemek[i]->Get_Cspe_Index()]->Get_p();
        pevhev[1] = 0.0;
        pevhev[2] = cspok[agelemek[i]->Get_Cspe_Index()]->Get_h();
        pevhev[3] = 0.0;
      } else {
        pevhev[0] = cspok[agelemek[i]->Get_Cspe_Index()]->Get_p();
        pevhev[1] = cspok[agelemek[i]->Get_Cspv_Index()]->Get_p();
        pevhev[2] = cspok[agelemek[i]->Get_Cspe_Index()]->Get_h();
        pevhev[3] = cspok[agelemek[i]->Get_Cspv_Index()]->Get_h();
      }
      // cout<<endl<<agelemek[i]->Get_nev();
      f[i] = agelemek[i]->f(pevhev);
      jv = agelemek[i]->df(pevhev);
      // pevhev.clear();

      Q_indx = i;
      pe_indx = agelemek[i]->Get_Cspe_Index();
      pv_indx = agelemek[i]->Get_Cspv_Index();

      jac[i][agelemek.size() + pe_indx] = jv.at(0);
      if (agelemek[i]->Get_Csp_db() == 2)
        jac[i][agelemek.size() + pv_indx] = jv.at(1);

      jac[i][Q_indx] = jv.at(2);
      col[i] = jv.at(3);
      jv.clear();
    }

    for (unsigned int i = 0; i < cspok.size(); i++) {
      f[agelemek.size() + i] = -cspok[i]->Get_fogy();
      for (unsigned int j = 0; j < cspok[i]->ag_be.size(); j++) {
        f[agelemek.size() + i] += agelemek[cspok[i]->ag_be.at(j)]->Get_mp();
        jac[agelemek.size() + i][cspok[i]->ag_be.at(j)] = +1.0;
      }
      for (unsigned int j = 0; j < cspok[i]->ag_ki.size(); j++) {
        f[agelemek.size() + i] -= agelemek[cspok[i]->ag_ki.at(j)]->Get_mp();
        jac[agelemek.size() + i][cspok[i]->ag_ki.at(j)] = -1.0;
      }
      col[agelemek.size() + i] = -cspok[i]->Get_fogy();
    }

    // Hibaszamitas
    //-------------------------------------------------------------------------------------------
    e_mp_r = e_mp;
    e_p_r = e_p;
    e_mp = 0.0;
    e_p = 0.0;
    for (unsigned int i = 0; i < agelemek.size(); i++) e_p += f[i] * f[i];
    for (unsigned int i = 0; i < cspok.size(); i++)
      e_mp += f[agelemek.size() + i] * f[agelemek.size() + i];

    e_mp = pow(e_mp, 0.5);
    e_p = pow(e_p, 0.5);
    /*cout.setf(ios::dec);
    cout.unsetf(ios::showpos);
    cout << endl << "  " << iter << "./" << iter_max;
    cout.setf(ios::scientific);
    cout.setf(ios::showpos);
    cout << " iteracio:\te_mp=" << e_mp << " (max:" << e_mp_max << ")   e_p=" <<
    e_p << " (max:" << e_p_max << ")";*/

    logfile_write(iter_info(x, f, iter, e_mp, e_p), 2);
    cout << iter_info(x, f, iter, e_mp, e_p);

    if ((e_p < e_p_max) && (e_mp < e_mp_max)) konv_ok = true;

    // Javult a megoldas?
    // Ha igen noveljuk a relaxiacios parametert, ha nem, csokkentjuk
    if ((e_mp < e_mp_r) || (e_p < e_p_r))
      relax = relax * relax_mul;
    else
      relax = relax / relax_mul;
    if (relax < 1e-3) relax = 1e-3;
    if (relax > RELAX_MAX) relax = RELAX_MAX;
    cout << "  RELAX=" << relax;

    strstrm.str("");
    strstrm << scientific << setprecision(3) << showpos;
    strstrm << endl << "\nJACOBI:" << endl << "          ";
    string strstrm_nev;
    const unsigned int MAX_NEV_HOSSZ = 12;
    for (unsigned int i = 0; i < agelemek.size(); i++) {
      strstrm_nev = agelemek[i]->Get_nev();
      while (strstrm_nev.size() < MAX_NEV_HOSSZ) strstrm_nev.append(" ");
      strstrm << "\tmp," << strstrm_nev;
    }
    for (unsigned int i = 0; i < cspok.size(); i++) {
      strstrm_nev = cspok[i]->Get_nev();
      while (strstrm_nev.size() < MAX_NEV_HOSSZ) strstrm_nev.append(" ");
      strstrm << "\t p," << strstrm_nev;
    }
    strstrm << endl;

    for (unsigned int i = 0; i < agelemek.size(); i++) {
      strstrm_nev = agelemek[i]->Get_nev();
      while (strstrm_nev.size() < MAX_NEV_HOSSZ) strstrm_nev.append(" ");
      strstrm << strstrm_nev;
      for (unsigned int j = 0; j < agelemek.size() + cspok.size(); j++)
        strstrm << "\t" << jac[i][j];
      strstrm << endl;
    }
    for (unsigned int i = 0; i < cspok.size(); i++) {
      strstrm_nev = cspok[i]->Get_nev();
      while (strstrm_nev.size() < MAX_NEV_HOSSZ) strstrm_nev.append(" ");
      strstrm << strstrm_nev;
      for (unsigned int j = 0; j < agelemek.size() + cspok.size(); j++)
        strstrm << "\t" << jac[agelemek.size() + i][j];
      strstrm << endl;
    }
    logfile_write(strstrm.str(), 3);

    switch (akt_megoldo) {
    case Linear:
      // Linearizalt egyenletek megoldasa:
      //-------------------------------------------------------------------------------------------
      NR::ludcmp(jac, indx, d);
      NR::lubksb(jac, indx, col);
      xu = col;
      for (int i = 0; i < N; i++) {
        xu[i] = x[i] + relax * xu[i];
        dx[i] = xu[i] - x[i];
      }
      break;

    case Newton_Raphson:

      // Jacobi inverzenek szamitasa:
      //-------------------------------
      NR::ludcmp(jac, indx, d);

      for (int j = 0; j < N; j++) {
        for (int i = 0; i < N; i++) col[i] = 0.0;
        col[j] = 1.0;
        NR::lubksb(jac, indx, col);
        for (int i = 0; i < N; i++) invjac[i][j] = col[i];
      }

      // Kiiras
      strstrm.str("");
      strstrm << scientific << setprecision(3) << showpos;
      strstrm << endl << "\nINVERZ JACOBI:" << endl << "          ";
      for (unsigned int i = 0; i < agelemek.size(); i++) {
        strstrm_nev = agelemek[i]->Get_nev();
        while (strstrm_nev.size() < MAX_NEV_HOSSZ) strstrm_nev.append(" ");
        strstrm << "\tmp," << strstrm_nev;
      }
      for (unsigned int i = 0; i < cspok.size(); i++) {
        strstrm_nev = cspok[i]->Get_nev();
        while (strstrm_nev.size() < MAX_NEV_HOSSZ) strstrm_nev.append(" ");
        strstrm << "\t p," << strstrm_nev;
      }
      strstrm << endl;

      for (unsigned int i = 0; i < agelemek.size(); i++) {
        strstrm_nev = agelemek[i]->Get_nev();
        while (strstrm_nev.size() < MAX_NEV_HOSSZ) strstrm_nev.append(" ");
        strstrm << strstrm_nev;
        for (unsigned int j = 0; j < agelemek.size() + cspok.size(); j++)
          strstrm << "\t" << invjac[i][j];
        strstrm << endl;
      }
      for (unsigned int i = 0; i < cspok.size(); i++) {
        strstrm_nev = cspok[i]->Get_nev();
        while (strstrm_nev.size() < MAX_NEV_HOSSZ) strstrm_nev.append(" ");
        strstrm << strstrm_nev;
        for (unsigned int j = 0; j < agelemek.size() + cspok.size(); j++)
          strstrm << "\t" << invjac[agelemek.size() + i][j];
        strstrm << endl;
      }
      logfile_write(strstrm.str(), 3);

      // Visszaszorzas
      //--------------------------------

      for (int i = 0; i < N; i++) {
        dx[i] = 0;
        for (int j = 0; j < N; j++) dx[i] += invjac[i][j] * f[j];
        xu[i] = x[i] - relax * dx[i];
      }
      break;

    default:
      cout << endl
           << endl
           << "HIBA: nem ismerek " << akt_megoldo << " megoldot!!!" << endl;
      break;
    }

    // Visszairas
    //-------------------------------------------------------------------------------------------
    for (unsigned int i = 0; i < agelemek.size(); i++)
      agelemek[i]->Set_mp(xu[i]);
    for (unsigned int i = 0; i < cspok.size(); i++)
      cspok[i]->Set_p(xu[agelemek.size() + i]);

    // Log
    strstrm << scientific << setprecision(3) << showpos;
    strstrm.str("");
    strstrm << endl
            << endl
            << "Elem neve      \t    xr             f(xr)            dx        "
            "     xuj";
    // A nevek legyenek azonos hosszuak...
    // string strstrm_nev;
    // const int MAX_NEV_HOSSZ=15;
    for (unsigned int i = 0; i < agelemek.size(); i++) {
      strstrm_nev = agelemek[i]->Get_nev();
      while (strstrm_nev.size() < MAX_NEV_HOSSZ) strstrm_nev.append(" ");
      strstrm << endl
              << strstrm_nev << "\t" << x[i] << "\t" << f[i] << "\t" << dx[i]
              << "\t" << xu[i];
    }
    for (unsigned int i = 0; i < cspok.size(); i++) {
      strstrm_nev = cspok[i]->Get_nev();
      while (strstrm_nev.size() < MAX_NEV_HOSSZ) strstrm_nev.append(" ");
      strstrm << endl
              << strstrm_nev << "\t" << x[agelemek.size() + i] << "\t"
              << f[agelemek.size() + i] << "\t" << dx[agelemek.size() + i]
              << "\t" << xu[agelemek.size() + i];
    }
    strstrm << endl << endl;
    logfile_write(strstrm.str(), 2);

    strstrm.str("");
    strstrm << "e_mp=" << e_mp << " (max:" << e_mp_max << ")   e_p=" << e_p
            << " (max:" << e_p_max << ")";
    logfile_write(strstrm.str(), 1);
    logfile_write("\n\n", 2);

    iter++;
  }

  strstrm.str("");
  if (konv_ok) {
    strstrm << endl << endl << "Normal befejezes...  :)" << endl;
    strstrm << list_results();
    logfile_write(strstrm.str(), 1);
    cout << strstrm.str();

    //      for (unsigned int i=0; i<agelemek.size(); i++) {
    //
    //          //agelemek.at(i)->build_res();
    //
    //          //agelemek.at(i)->Set_FolyTerf();
    //
    //      }
  } else
    strstrm << endl
            << endl
            << "HIBA: Maximalis lepesszamot elertem, de a megadott hibahatar "
            "felett vagyok..."
            << endl;
  logfile_write(strstrm.str(), 1);
  // int int1; cin>>int1;
  return konv_ok;
}

//--------------------------------------------------------------
void Staci::export_connected_nodes() {
  ofstream nodelistfile("nodelist.txt", ios::trunc);
  for (unsigned int i = 0; i < cspok.size(); i++)
    nodelistfile << cspok[i]->Get_nev() << "\n";
  nodelistfile.close();

  ofstream nodefile("connected_nodes.txt", ios::trunc);
  for (unsigned int i = 0; i < agelemek.size(); i++) {
    if (agelemek[i]->Get_Csp_db() == 2) {
      string innen = cspok[agelemek[i]->Get_Cspe_Index()]->Get_nev();
      string ide = cspok[agelemek[i]->Get_Cspv_Index()]->Get_nev();
      string nev = agelemek[i]->Get_nev();
      nodefile << "\n" << innen << "," << ide;
    }
  }
  nodefile.close();
}

//--------------------------------------------------------------
void Staci::ini() {
  ostringstream strstrm;
  if (!van_ini) {
    strstrm << endl << endl << " Automatic inicialization...";

    strstrm << endl << "\tnodal pressures:\t" << p_init << " mwc ";
    for (unsigned int i = 0; i < cspok.size(); i++)
      cspok.at(i)->Ini(1., p_init);

    strstrm << endl << "\tedge mass flow rates:\t" << mp_init << " kg/s ...";
    for (unsigned int i = 0; i < agelemek.size(); i++)
      agelemek.at(i)->Ini(1., mp_init);

  } else {
    strstrm << endl
            << endl
            << "loading initial values from file " << ini_file << endl;
    data_io datta_io(ini_file.c_str());
    datta_io.load_ini_values(cspok, agelemek);
  }

  // Flush info into logfile and screen
  if (debug_level > 0) {
    logfile_write(strstrm.str(), 1);
    // cout << strstrm.str();
  }

  // Set log file
  // for (unsigned int i = 0; i < agelemek.size(); i++) {
  //        cout << endl << "\t SetLogFile: " << (i+1) << "/" <<
  //        agelemek.size() << ":" << agelemek.at(i)->Get_nev()<<" logfile:
  //        "<<out_file.c_str();
  //            cout<<agelemek.at(i)->Info();
  // agelemek.at(i)->SetLogFile(out_file.c_str());
  // }
}

//--------------------------------------------------------------
void Staci::ini(const Staci * IniStaci) {
  for (unsigned int i = 0; i < cspok.size(); i++)
    cspok.at(i)->Set_p(IniStaci->cspok.at(i)->Get_p());

  for (unsigned int i = 0; i < agelemek.size(); i++)
    agelemek.at(i)->Set_mp(IniStaci->agelemek.at(i)->Get_mp());

  // Set log file
  // for (unsigned int i = 0; i < agelemek.size(); i++) {
  //        cout << endl << "\t SetLogFile: " << (i+1) << "/" <<
  //        agelemek.size() << ":" << agelemek.at(i)->Get_nev()<<" logfile:
  //        "<<out_file.c_str();
  //            cout<<agelemek.at(i)->Info();
  // agelemek.at(i)->SetLogFile(out_file.c_str());
  // }
}

//--------------------------------------------------------------
void Staci::post_process_res_file() {
  ostringstream strstrm;
  string s;
  ifstream ifs(res_file.c_str());
  int i = 0;
  while (!ifs.eof()) {
    getline(ifs, s);
    if (i == 0)
      strstrm << "<?xml version=\"1.0\" encoding=\"ISO-8859-2\"?>";
    else
      strstrm << trim(s, "\t");
    i++;
  }
  ifs.close();

  ofstream ofs(res_file.c_str());
  ofs << strstrm.str();
  ofs.close();
}

//--------------------------------------------------------------
string Staci::trim(string s, const string drop) {
  string r = s.erase(s.find_last_not_of(drop) + 1);
  return r.erase(0, r.find_first_not_of(drop));
}

//--------------------------------------------------------------
void Staci::set_up_transport() {
  unsigned int i;
  van_ini = true;
  data_io datta_io(def_file.c_str());
  datta_io.load_ini_values(cspok, agelemek);
  solve_system();

  vector<double> v;
  double hossz;
  for (i = 0; i < agelemek.size(); i++) {
    v.clear();
    // cout<<endl<<"Agelem: "<<agelemek.at(i)->Get_nev();
    if (agelemek.at(i)->Get_Tipus() == "Csatorna") {
      v = agelemek.at(i)->Get_res("v");
      hossz = agelemek.at(i)->Get_dprop("L");
      agelemek.at(i)->set_up_grid(0.0, v, hossz);
    } else {
      if (agelemek.at(i)->Get_Tipus() == "Cso") {
        for (int j = 0; j < 10; j++) {
          v.push_back(agelemek.at(i)->Get_v());
        }
        hossz = agelemek.at(i)->Get_dprop("L");
        agelemek.at(i)->set_up_grid(0.0, v, hossz);
      } else {
        v.push_back(agelemek.at(i)->Get_v());
        v.push_back(agelemek.at(i)->Get_v());
        hossz = 10.0;
        agelemek.at(i)->set_up_grid(0.0, v, hossz);
        agelemek.at(i)->Set_cdt(10 * 60);
      }
    }

    // cout<<agelemek.at(i)->Info();
    // cout<<agelemek.at(i)->show_grid(ido);
    // int int1;cin>>int1;
  }
  transp_dt = 1e100;
  for (i = 0; i < agelemek.size(); i++) {
    bool kell_eloszlas = false;

    if (agelemek.at(i)->Get_Tipus() == "Cso") {
      if (fabs(agelemek.at(i)->Get_dprop("erdesseg")) > 1e-5)
        kell_eloszlas = true;
    }
    if (agelemek.at(i)->Get_Tipus() == "Csatorna") kell_eloszlas = true;

    if (kell_eloszlas) {
      // cout << ", dt=" << agelemek.at(i)->cdt;
      if (transp_dt > agelemek.at(i)->cdt) {
        transp_dt = agelemek.at(i)->cdt;
        // cout << "*";
      }
    }
  }
  // double transp_dt_min=1; // s
  // if (transp_dt<transp_dt_min)
  //    transp_dt=transp_dt_min;

  tt_length = atof(datta_io.read_setting("tt_length").c_str());

  ostringstream strstrm;
  strstrm << endl
          << endl
          << "transp_dt=" << transp_dt << "s, tt_length=" << tt_length << "h"
          << endl;
  strstrm << endl << "Grid setup complete." << endl;
  logfile_write(strstrm.str(), 2);
  cout << strstrm.str();
}

//--------------------------------------------------------------
void Staci::solve_transport(int mode) {
  ostringstream strstrm;
  strstrm << "Starting computation..." << endl;
  logfile_write(strstrm.str(), 2);
  cout << strstrm.str();

  cout << fixed << setprecision(2) << setw(5) << setfill(' ') << right;
  ido = 0.;
  double counter = 0, szazalek;
  while (ido < tt_length * 3600) {
    szazalek = round(ido / (tt_length * 3600) * 100);
    // szazalek=round(10*szazalek)/10;
    if ((debug_level) > 0) progress_file_write(szazalek);
    if (counter < szazalek) {
      cout << "\n\t\t" << szazalek << " %, time = " << (ido / 3600.);
      cout << "h, oldest fluid particle: " << (get_oldest() / 3600.) << "h";
      counter += 1;
    }
    transport_step(transp_dt);
    ido += transp_dt;
  }
}

//--------------------------------------------------------------
void Staci::save_transport(int mode) {
  res_file = def_file;
  cout << endl << "Kimeneti fajl neve:" << res_file << endl;
  data_io datta_io(res_file.c_str());
  datta_io.save_transport(mode, cspok, agelemek);
}

//--------------------------------------------------------------
double Staci::get_oldest() {
  double tmax = 0;
  double tmp;
  for (unsigned int i = 0; i < cspok.size(); i++) {
    tmp = cspok.at(i)->Get_dprop("konc_atlag");
    if (tmp > tmax) tmax = tmp;
  }
  return tmax;
}

//--------------------------------------------------------------
void Staci::transport_step(double dt) {
  bool transp_debug = false;
  vector<int> temp;
  double seb;

  // 1. lepes: csomopontok atlagos koncentraciojanak szamitasa
  double szaml, nevezo, c, m;
  for (unsigned int i = 0; i < cspok.size(); i++) {
    // Ha meg van adva beadagolt klor, a csomopontban pontosan ennyi
    if (fabs(cspok.at(i)->Get_dprop("cl_be")) > 1e-10) {
      cspok.at(i)->Set_dprop("konc_atlag", cspok.at(i)->Get_dprop("cl_be"));
    }
    // Ha nincs megadva, kevered�st kell sz�molni
    else {
      szaml = 0.0;
      nevezo = 0.0;

      if (transp_debug) cout << endl << "\t bemeno agak: (";
      temp.clear();
      temp = cspok.at(i)->ag_be;
      if (transp_debug) cout << temp.size() << "db)";

      if (temp.size() > 0) {
        for (unsigned int j = 0; j < temp.size(); j++) {
          // pozitiv atlagsebessegu agak kellenek csak
          // agelemek.at(temp.at(j))->show_grid(ido);
          seb = agelemek.at(temp.at(j))->mean(agelemek.at(temp.at(j))->vel);
          if (transp_debug)
            cout << endl
                 << "\t\t" << agelemek.at(temp.at(j))->Get_nev()
                 << ", vatlag=" << fixed << setprecision(3) << seb;
          int utolso = agelemek.at(temp.at(j))->vel.size() - 1;
          if (agelemek.at(temp.at(j))->vel.at(utolso) > 0) {
            c = agelemek.at(temp.at(j))->konc.at(utolso);
            m = agelemek.at(temp.at(j))->Get_mp();
            if (transp_debug) cout << " mp= " << m << " c=" << c;
            szaml += fabs(m) * c;
            nevezo += fabs(m);
          }
        }
      }

      if (transp_debug) cout << endl << "\t kimeno agak: (";
      temp.clear();
      temp = cspok.at(i)->ag_ki;
      if (transp_debug) cout << temp.size() << "db)";

      if (temp.size() > 0) {
        for (unsigned int j = 0; j < temp.size(); j++) {
          // negativ sebessegu agak kellenek csak
          seb = agelemek.at(temp.at(j))->mean(agelemek.at(temp.at(j))->vel);
          if (transp_debug)
            cout << endl
                 << "\t\t" << agelemek.at(temp.at(j))->Get_nev()
                 << ", vatlag=" << fixed << setprecision(3) << seb;

          if (agelemek.at(temp.at(j))->vel.at(0) < 0) {
            c = agelemek.at(temp.at(j))->konc.at(0);
            m = agelemek.at(temp.at(j))->Get_mp();
            if (transp_debug) cout << " mp= " << m << " c=" << c;
            szaml += fabs(m) * c;
            nevezo += fabs(m);
          }
        }
      }
      double c, fogy, cl_be;
      fogy = cspok.at(i)->Get_dprop("fogy");
      cl_be = cspok.at(i)->Get_dprop("cl_be");
      // Ha "kifele" megy a viz, attol a csomopontban nem higul
      if (fogy > 0) {
        fogy = 0.0;
        cl_be = 0.0;
      }

      if (mode == 1) {
        // if (fabs(nevezo - fogy) > 1.e-10)
        c = (szaml) / (nevezo - fogy);
        // else
        //{
        // TODO: logfile_write-ba!
        // cout << endl << "ERROR!!!! Staci::transport_step():
        // fabs(nevezo-fogy)=" << fabs(nevezo - fogy) << "<1.e-10" << endl <<
        // endl;
        // exit(-1);
        //}
      } else {
        // if (nevezo + fabs(fogy) > 1.e-10)
        c = (szaml + cl_be * fabs(fogy)) / (nevezo + fabs(fogy));
        // else
        //{
        // TODO: logfile_write-ba!
        //    cout << endl << "ERROR!!!! Staci::transport_step(): nevezo +
        //    fabs(fogy)" << (nevezo + fabs(fogy)) << "<1.e-10" << endl << endl;
        // exit(-1);
        //}
      }

      cspok.at(i)->Set_dprop("konc_atlag", c);
      if ((cspok.at(i)->Get_nev() == "NODE_1346650") ||
          (cspok.at(i)->Get_nev() == "NODE_1337960")) {
        cout << "\n\t" << cspok.at(i)->Get_nev() << ": szaml=" << szaml
             << ", cl_be*fabs(fogy)=" << (cl_be * fabs(fogy))
             << ", nevezo=" << nevezo;
        cout << "\n\t\t atlagos koncentracio: " << c << endl;
      }
    }
  }

  // cout << endl << "Node update complete.\n";
  // cin.get();

  // 2. lepes: agelemek leptetese
  for (unsigned int i = 0; i < agelemek.size(); i++) {
    // cout << endl << "\t updating " << agelemek.at(i)->Get_nev() << "...";

    // Ha nulla az �rdess�g, akkor csak �sszek�t� cs� �s szkippelj�k
    bool kell_eloszlas = false;

    if (agelemek.at(i)->Get_Tipus() == "Cso") {
      if ((fabs(agelemek.at(i)->Get_dprop("erdesseg"))) > 1e-5)
        kell_eloszlas = true;
    }
    if (agelemek.at(i)->Get_Tipus() == "Csatorna") kell_eloszlas = true;

    //      cout<<endl<<agelemek.at(i)->Get_nev()<<",
    //      kell_eloszlas="<<kell_eloszlas;
    //      int int1;
    //      cin>>int1;

    if (kell_eloszlas) {
      if ((agelemek.at(i)->mean(agelemek.at(i)->vel)) > 0) {
        double konc_eleje, konc, konc_e, vel_e, dx;
        dx = agelemek.at(i)->cL / agelemek.at(i)->vel.size();
        if (agelemek.at(i)->Get_Csp_db() == 1) {
          konc_eleje = 0.0;
          agelemek.at(i)->konc.at(0) = 0.0;
          agelemek.at(i)->konc.at(1) = 0.0;
        } else {
          konc_eleje = cspok.at(agelemek.at(i)->Get_Cspe_Index())
                       ->Get_dprop("konc_atlag");
          for (unsigned int j = agelemek.at(i)->vel.size() - 1; j > 0; j--) {
            konc_e = agelemek.at(i)->konc.at(j - 1);
            vel_e = agelemek.at(i)->vel.at(j - 1);
            konc = agelemek.at(i)->konc.at(j);
            agelemek.at(i)->konc.at(j) =
              konc - dt / dx * vel_e * (konc - konc_e) + dt * teta(konc_e, i);

            /*if (agelemek.at(i)->Get_nev()=="PIPE116") {
              cout<<endl<<"i="<<i<<", konc="<<konc<<", konc_e="
              <<konc_e<<", vel_e="<<vel_e<<", teta="<<teta;
              }*/

            if (agelemek.at(i)->konc.at(j) < 0)
              agelemek.at(i)->konc.at(j) = 0.0;
          }
          agelemek.at(i)->konc.at(0) = konc_eleje;
        }
      } else {
        double konc_vege, konc, konc_u, vel_u, dx;
        dx = agelemek.at(i)->cL / agelemek.at(i)->vel.size();
        if (agelemek.at(i)->Get_Csp_db() == 1) {
          konc_vege = cspok.at((unsigned long)agelemek.at(i)->Get_Cspv_Index())
                      ->Get_dprop("konc_atlag");
          agelemek.at(i)->konc.at(0) = konc_vege;
          agelemek.at(i)->konc.at(1) = konc_vege;
        } else {
          konc_vege = cspok.at((unsigned long)agelemek.at(i)->Get_Cspv_Index())
                      ->Get_dprop("konc_atlag");
          for (unsigned int j = 0; j < agelemek.at(i)->vel.size() - 1; j++) {
            konc_u = agelemek.at(i)->konc.at(j + 1);
            vel_u = agelemek.at(i)->vel.at(j + 1);
            konc = agelemek.at(i)->konc.at(j);
            agelemek.at(i)->konc.at(j) =
              konc - dt / dx * vel_u * (konc_u - konc) + dt * teta(konc_u, i);
            if (agelemek.at(i)->konc.at(j) < 0)
              agelemek.at(i)->konc.at(j) = 0.0;
          }
          agelemek.at(i)->konc.at(agelemek.at(i)->konc.size() - 1) = konc_vege;
        }
      }
    } else {
      double konc =
        cspok.at(agelemek.at(i)->Get_Cspe_Index())->Get_dprop("konc_atlag");

      if (agelemek.at(i)->Get_Csp_db() == 2) {
        double konc_vege =
          cspok.at(agelemek.at(i)->Get_Cspv_Index())->Get_dprop("konc_atlag");
        konc = (konc + konc_vege) / 2;
      }

      // cout << endl << "dt=" << dt << ", " << agelemek.at(i)->show_grid(0.0)
      // << "\n END OF GRID\n";

      for (unsigned int j = 0; j < agelemek.at(i)->konc.size(); j++) {
        // cout << "\t agelemek.at(i)->konc.size()=" <<
        // agelemek.at(i)->konc.size();
        // cin.get();
        agelemek.at(i)->konc.at(j) = konc;
      }
    }
    if (agelemek.at(i)->Get_nev() == "PIPE_602407") {
      cout << endl << "dt=" << dt << ", " << agelemek.at(i)->show_grid(0.0);
      /*cin.get();*/
    }
  }
  // cout << endl << "Edge update complete.\n";
  // cin.get();
}

//--------------------------------------------------------------
double Staci::teta(double konc, const int i) {
  if (mode == 1)
    return 1.0;
  else {
    double Rh, cl_k, cl_w;
    if (agelemek.at(i)->Get_Tipus() == "Csatorna") {
      Rh = agelemek.at(i)->Get_dprop("Rh");
      cl_k = agelemek.at(i)->Get_dprop("cl_k");
      cl_w = agelemek.at(i)->Get_dprop("cl_w");
    } else {
      if (agelemek.at(i)->Get_Tipus() == "Cso") {
        Rh = agelemek.at(i)->Get_dprop("Rh");
        cl_k = agelemek.at(i)->Get_dprop("cl_k");
        cl_w = agelemek.at(i)->Get_dprop("cl_w");
      } else {
        Rh = 1;
        cl_k = 0.0;
        cl_w = 0.0;
      }
    }
    double Sh = 3.65;
    double kf = Sh;
    double tag1 = -cl_k * konc;
    double tag2 = -kf / Rh * (konc - cl_w);
    return (tag1 + 0 * tag2);
  }
}

//--------------------------------------------------------------
void Staci::copy_file(const string in_f_nev, const string out_f_nev) {
  ifstream infile(in_f_nev.c_str(), ifstream::in);
  ofstream outfile(out_f_nev.c_str(), ifstream::trunc);
  outfile << infile.rdbuf();
}

//--------------------------------------------------------------
void Staci::progress_file_write(double percent) {
  ofstream pfile(progress_file.c_str(), ios::trunc);
  pfile << fixed << setprecision(1) << percent << "\n";
  pfile.close();
}

//--------------------------------------------------------------
void Staci::Set_FolyTerf() {
  double osszeg = 0;
  for (unsigned int i = 0; i < agelemek.size(); i++) {
    osszeg += agelemek.at(i)->Get_FolyTerf();
    // cout << endl << "\t" << agelemek.at(i)->Get_nev() << ":\t  " <<
    // agelemek.at(i)->Get_FolyTerf() << ",\t szum:" << osszeg;
  }
  FolyMenny = osszeg;

}

//--------------------------------------------------------------
void Staci::nr_solver(Vec_DP x, Vec_DP f) {
  int N = x.size();
  Vec_DP col(N), b(N), dx(N), xu(N);
  Vec_INT indx(N);
  Mat_DP invjac(N, N), jac(N, N);
  DP d;

  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++) jac[i][j] = m_jac[i][j];

  NR::ludcmp(jac, indx, d);

  for (int j = 0; j < N; j++) {
    for (int i = 0; i < N; i++) col[i] = 0.0;
    col[j] = 1.0;
    NR::lubksb(jac, indx, col);
    for (int i = 0; i < N; i++) invjac[i][j] = col[i];
  }

  for (int i = 0; i < N; i++) {
    dx[i] = 0;
    for (int j = 0; j < N; j++) dx[i] += invjac[i][j] * f[j];
    xu[i] = x[i] - m_relax * dx[i];
  }

  // Visszairas
  //------------------------
  for (unsigned int i = 0; i < agelemek.size(); i++) agelemek[i]->Set_mp(xu[i]);
  for (unsigned int i = 0; i < cspok.size(); i++)
    cspok[i]->Set_p(xu[agelemek.size() + i]);
}

//--------------------------------------------------------------
bool Staci::umfpack_solver(Vec_I_DP xr, Vec_I_DP f) {
  // Build sparse matrix
  /* Ti[k] is row index of entry k, as matrix is scanned columnwise */

  //    cout<<endl<<"\n m_nnz="<<m_nnz<<endl;

  vector<int> vTi;
  vTi.reserve(m_nnz);
  /* Tj[k] is column index of entry k, as matrix is scanned columnwise */
  vector<int> vTj;
  vTj.reserve(m_nnz);
  /* value of entry k, as matrix is scanned columnwise */
  vector<double> vTx;
  vTx.reserve(m_nnz);

  int n = agelemek.size() + cspok.size();
  int nz = 0;
  for (int col = 0; col < n; col++)
    for (int row = 0; row < n; row++)
      if (!m_is_element_empty[row][col]) {
        // vTi.push_back(row);
        // vTj.push_back(col);
        // vTx.push_back(m_jac[row][col]);
        vTi[nz] = row;
        vTj[nz] = col;
        vTx[nz] = m_jac[row][col];
        nz++;
      }

  int Ap[n + 1];
  int Ai[nz];
  double Ax[nz];
  int status;
  int Ti[nz];
  int Tj[nz];
  double Tx[nz];
  double dx[n];
  for (int i = 0; i < n; i++) dx[i] = 0.;

  void *Symbolic, *Numeric;

  for (int i = 0; i < nz; i++) {
    Ti[i] = vTi[i];
    Tj[i] = vTj[i];
    Tx[i] = vTx[i];
  }

  double b[n];
  for (int i = 0; i < n; i++) b[i] = f[i];

  /* convert matrix from triplet form to compressed-column form */
  status = umfpack_di_triplet_to_col(n, n, nz, Ti, Tj, Tx, Ap, Ai, Ax, NULL);

  /* symbolic analysis */
  status = umfpack_di_symbolic(n, n, Ap, Ai, Ax, &Symbolic, NULL, NULL);

  /* LU factorization */
  umfpack_di_numeric(Ap, Ai, Ax, Symbolic, &Numeric, NULL, NULL);

  umfpack_di_free_symbolic(&Symbolic);

  /* solve system */
  umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax, dx, b, Numeric, NULL, NULL);

  umfpack_di_free_numeric(&Numeric);

  bool success = true;
  for (int i = 0; i < n; i++)
    if (isnan(dx[i])) {
      //            cout << "\n\n!!!!\nStaci.cpp, umfpack_solver() -> x[" << i
      //                    << "]=NaN!!!\n\n";
      success = false;
      break;
    }

  // Visszairas
  //------------------------
  if (success) {
    for (unsigned int i = 0; i < agelemek.size(); i++)
      agelemek[i]->Set_mp(xr[i] - m_relax * dx[i]);
    for (unsigned int i = 0; i < cspok.size(); i++)
      cspok[i]->Set_p(xr[agelemek.size() + i] -
                      m_relax * dx[agelemek.size() + i]);
  }

  return success;
}

//--------------------------------------------------------------
void Staci::Compute_Head_Losses() {
  double pe, pv, he, hv, dh;
  // cout << endl << endl << "Nyomasesesek szamitasa...";
  for (unsigned int i = 0; i < agelemek.size(); i++) {
    if (agelemek[i]->Get_Csp_db() == 1) {
      agelemek[i]->Set_head_loss(0.0);
    } else {
      pe = cspok[agelemek[i]->Get_Cspe_Index()]->Get_p();
      pv = cspok[agelemek[i]->Get_Cspv_Index()]->Get_p();
      he = cspok[agelemek[i]->Get_Cspe_Index()]->Get_h();
      hv = cspok[agelemek[i]->Get_Cspv_Index()]->Get_h();
      dh = fabs(pe + he - pv - hv);
      agelemek[i]->Set_head_loss(dh);
    }
  }
  // cout << "\t OK";
}

//--------------------------------------------------------------

void Staci::Print_Jacobian()
//                void Staci::Print_Jacobian(vector<vector<double> > jac)

{
  ostringstream strstrm;
  string strstrm_nev;
  const unsigned int MAX_NEV_HOSSZ = 15;

  strstrm.str("");
  strstrm << scientific << setprecision(3) << showpos;
  strstrm << endl << "\nJACOBI:" << endl << "          ";

  for (unsigned int i = 0; i < agelemek.size(); i++) {
    strstrm_nev = agelemek[i]->Get_nev();
    while (strstrm_nev.size() < MAX_NEV_HOSSZ) strstrm_nev.append(" ");
    strstrm << "\tmp," << strstrm_nev;
  }
  for (unsigned int i = 0; i < cspok.size(); i++) {
    strstrm_nev = cspok[i]->Get_nev();
    while (strstrm_nev.size() < MAX_NEV_HOSSZ) strstrm_nev.append(" ");
    strstrm << "\t p," << strstrm_nev;
  }
  strstrm << endl;

  for (unsigned int i = 0; i < agelemek.size(); i++) {
    strstrm_nev = agelemek[i]->Get_nev();
    while (strstrm_nev.size() < MAX_NEV_HOSSZ) strstrm_nev.append(" ");
    strstrm << strstrm_nev;
    for (unsigned int j = 0; j < agelemek.size() + cspok.size(); j++)
      strstrm << "; " << m_jac[i][j];
    strstrm << endl;
  }
  for (unsigned int i = 0; i < cspok.size(); i++) {
    strstrm_nev = cspok[i]->Get_nev();
    while (strstrm_nev.size() < MAX_NEV_HOSSZ) strstrm_nev.append(" ");
    strstrm << strstrm_nev;
    for (unsigned int j = 0; j < agelemek.size() + cspok.size(); j++)
      strstrm << ";" << m_jac[agelemek.size() + i][j];
    strstrm << endl;
  }

  ofstream JacFile;
  JacFile.open("dfdx.txt");
  JacFile << strstrm.str();
  JacFile.close();
  //  logfile_write(strstrm.str(), 3);
}

//--------------------------------------------------------------
string Staci::iter_info(Vec_DP x, Vec_DP f, int iter, double e_mp, double e_p) {
  m_ss.str("");

  int num = 0;
  for (unsigned int i = 0; i < agelemek.size(); i++)
    if (agelemek.at(i)->force_more_iter)
      num++;

  m_ss.setf(ios::dec);
  m_ss.unsetf(ios::showpos);
  m_ss << endl << " iter. # " << iter << "./" << iter_max;
  m_ss << setprecision(2) << scientific;  // << number << std::endl;
  m_ss << " e_mp=" << e_mp << ",  e_p=" << e_p;
  m_ss << setprecision(2) << fixed;
  m_ss << "  relax.=" << m_relax << ", relax.mul.=" << m_relax_mul;
  m_ss << " (" << num << ")";

  return m_ss.str();
}

//--------------------------------------------------------------
void Staci::compute_error(Vec_DP f, double & e_mp, double & e_p, double & e_mp_r,
                          double & e_p_r, bool & konv_ok) {
  e_mp_r = e_mp;
  e_p_r = e_p;
  e_mp = 0.0;
  e_p = 0.0;
  for (unsigned int i = 0; i < agelemek.size(); i++) e_p += f[i] * f[i];
  for (unsigned int i = 0; i < cspok.size(); i++)
    e_mp += f[agelemek.size() + i] * f[agelemek.size() + i];

  e_mp = pow(e_mp, 0.5);
  e_p = pow(e_p, 0.5);

  if ((e_p < e_p_max) && (e_mp < e_mp_max)) konv_ok = true;
}

//--------------------------------------------------------------
void Staci::update_relax(double e_mp, double e_p, double & e_mp_r,
                         double & e_p_r) {
  m_RELAX_MIN = 0.01;
  m_RELAX_MAX = 1.0;
  double hiba, hiba_r;
  if (e_p > e_mp) {
    hiba = e_p;
    hiba_r = e_p_r;
  } else {
    hiba = e_mp;
    hiba_r = e_mp_r;
  }
  if (hiba < hiba_r)
    m_relax = m_relax * m_relax_mul;
  else
    m_relax = m_relax / 10.;
  if (m_relax < m_RELAX_MIN) m_relax = m_RELAX_MIN;
  if (m_relax > m_RELAX_MAX) m_relax = m_RELAX_MAX;

  e_mp_r = e_mp;
  e_p_r = e_p;
}

//--------------------------------------------------------------
void Staci::list_all_elements() {
  stringstream str;
  for (unsigned int i = 0; i < agelemek.size(); i++)
    str << agelemek.at(i)->Get_Tipus() << ";\t" << agelemek.at(i)->Get_nev()
        << endl;
  for (unsigned int i = 0; i < cspok.size(); i++)
    str << "Node;\t" << cspok.at(i)->Get_nev() << endl;

  cout << endl
       << endl
       << "Element list is in file element_list.txt." << endl
       << endl;

  ofstream OutFile;
  OutFile.open("element_list.txt");
  OutFile << str.str();
  OutFile.close();
}

//--------------------------------------------------------------
void Staci::Compute_dfdmu() {
  bool megvan = false;
  m_dfdmu.clear();

  // element_ID - elem neve
  // property_ID - property neve
  // Jacobi sorok: agelemek, majd csomopontok
  // Jacobi oszlopok: mp, majd p

  for (int i = 0; i < agelemek.size(); i++) {
    if (((property_ID == "diameter") || (property_ID == "friction_coeff"))  &&
        (element_ID == agelemek.at(i)->Get_nev())) {
      double dfdmu = agelemek.at(i)->Get_dfdmu(property_ID);
      // cout<<endl<<"Hali! element:"<<agelemek.at(i)->Get_nev()<<", dfdmu="<<dfdmu;
      // cin.get();
      m_dfdmu.push_back(dfdmu);
      megvan = true;
    } else
      m_dfdmu.push_back(0.0);
  }

  for (int i = 0; i < cspok.size(); i++) {
    if ((property_ID == "demand") && (element_ID == cspok.at(i)->Get_nev())) {
      m_dfdmu.push_back(-1.0);
      megvan = true;
    } else
      m_dfdmu.push_back(0.0);
  }

  if (!megvan) {
    stringstream strstrm;
    strstrm.str("");
    strstrm << "\nStaci::Compute_dfdmu(): !!! element_ID: " << element_ID
            << ", property_ID: " << property_ID
            << " not found (property_ID = demand | diameter)!";
    logfile_write(strstrm.str(), 1);
    cout << strstrm.str();
    StaciException hiba(strstrm.str());
    throw hiba;
  }
}

//--------------------------------------------------------------
void Staci::Compute_dxdmu() {
  Compute_dfdmu();

  // Build sparse matrix
  /* Ti[k] is row index of entry k, as matrix is scanned columnwise */

  vector<int> vTi;
  vTi.reserve(m_nnz);
  /* Tj[k] is column index of entry k, as matrix is scanned columnwise */
  vector<int> vTj;
  vTj.reserve(m_nnz);
  /* value of entry k, as matrix is scanned columnwise */
  vector<double> vTx;
  vTx.reserve(m_nnz);

  int n = agelemek.size() + cspok.size();
  int nz = 0;
  for (int col = 0; col < n; col++)
    for (int row = 0; row < n; row++)
      if (!m_is_element_empty[row][col]) {
        // vTi.push_back(row);
        // vTj.push_back(col);
        // vTx.push_back(m_jac[row][col]);
        vTi[nz] = row;
        vTj[nz] = col;
        vTx[nz] = m_jac[row][col];
        nz++;
      }

  int Ap[n + 1];
  int Ai[nz];
  double Ax[nz];
  int status;
  int Ti[nz];
  int Tj[nz];
  double Tx[nz];
  double dx[n];
  for (int i = 0; i < n; i++) dx[i] = 0.;

  void *Symbolic, *Numeric;

  for (int i = 0; i < nz; i++) {
    Ti[i] = vTi[i];
    Tj[i] = vTj[i];
    Tx[i] = vTx[i];
  }

  double b[n];
  for (int i = 0; i < n; i++) b[i] = -m_dfdmu[i];

  /* convert matrix from triplet form to compressed-column form */
  status = umfpack_di_triplet_to_col(n, n, nz, Ti, Tj, Tx, Ap, Ai, Ax, NULL);

  /* symbolic analysis */
  status = umfpack_di_symbolic(n, n, Ap, Ai, Ax, &Symbolic, NULL, NULL);

  /* LU factorization */
  umfpack_di_numeric(Ap, Ai, Ax, Symbolic, &Numeric, NULL, NULL);

  umfpack_di_free_symbolic(&Symbolic);

  /* solve system */
  umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax, dx, b, Numeric, NULL, NULL);

  umfpack_di_free_numeric(&Numeric);

  // bool success = true;
  for (int i = 0; i < n; i++)
    if (isnan(dx[i])) {
      cout << "\n\n!!!!\nStaci.cpp, dxdmu() -> x[" << i << "]=NaN!!!\n\n";
      break;
    }

  // Visszairas
  //------------------------
  m_dxdmu.clear();

  for (unsigned int i = 0; i < n; i++) m_dxdmu.push_back(dx[i]);
}

//--------------------------------------------------------------
void Staci::Print_dfdmu() {
  // Compute_dfdmu();

  ostringstream strstrm;
  string strstrm_nev;
  // const unsigned int MAX_NEV_HOSSZ = 15;

  strstrm.str("");
  strstrm << scientific << setprecision(3) << showpos;
  strstrm << "Parameter: " << element_ID << ", " << property_ID;

  for (unsigned int i = 0; i < agelemek.size(); i++)
    strstrm << "\n"
            << i << "; (" << agelemek.at(i)->Get_nev() << "); " << m_dfdmu[i];

  for (unsigned int i = 0; i < cspok.size(); i++)
    strstrm << "\n"
            << i << "; (" << cspok.at(i)->Get_nev() << "); "
            << m_dfdmu[agelemek.size() + i];

  ofstream OutFile;
  OutFile.open("dfdmu.txt");
  OutFile << strstrm.str();
  OutFile.close();
}

//--------------------------------------------------------------
void Staci::Save_Sensitivity() {
  /*! \todo Sensitivity is flushed to <concentration>
     *
     *  A <tag> must be created to be used for sensitivity information in the
   * data file.
     */

  for (unsigned int i = 0; i < agelemek.size(); i++) {
    element_ID = agelemek.at(i)->Get_nev();
    property_ID = "concentration";
    double sens = m_dxdmu[i];
    agelemek.at(i)->Set_dprop("konc_atlag", fabs(sens));
    save_mod_prop(true);
  }

  for (unsigned int i = 0; i < cspok.size(); i++) {
    element_ID = cspok.at(i)->Get_nev();
    property_ID = "concentration";
    double sens = m_dxdmu[agelemek.size() + i];
    cspok.at(i)->Set_dprop("konc_atlag", fabs(sens));
    save_mod_prop(true);
  }
}

//--------------------------------------------------------------
void Staci::Print_dxdmu() {
  ostringstream strstrm;
  strstrm.str("");
  strstrm << scientific << setprecision(3) << showpos;
  strstrm << "Parameter: " << element_ID << ", " << property_ID;
  strstrm << endl << "No., expl. , dx/dmu, x , fabs((dx/dmu)/x) ";

  for (unsigned int i = 0; i < agelemek.size(); i++)
    strstrm << "\n"
            << i << ";(mp @ " << agelemek.at(i)->Get_nev()
            << "); s = " << m_dxdmu[i]
            << "; mp(kg/s) = " << agelemek.at(i)->Get_mp()
            << ", s/mp = " << fabs(m_dxdmu[i] / agelemek.at(i)->Get_mp());

  for (unsigned int i = 0; i < cspok.size(); i++)
    strstrm << "\n"
            << i << ";(p @ " << cspok.at(i)->Get_nev() << "); "
            << m_dxdmu[agelemek.size() + i]
            << "; p(vom)=" << cspok.at(i)->Get_p() << ", s/p = "
            << fabs(m_dxdmu[agelemek.size() + i] / cspok.at(i)->Get_h());

  ofstream OutFile;
  stringstream fname;
  fname.str("");
  fname << "dxdmu_" << element_ID << "_" << property_ID << ".txt";
  OutFile.open(fname.str().c_str());
  OutFile << strstrm.str();
  OutFile.close();
}

//--------------------------------------------------------------
/**
       * Builds sensitivity matrix
       * @param <string> parameter : "diameter" | "friction_coeff" | "demand"
       * @param <int> scale : if set to 1, the matrix is rescaled with the largest element
       * @return <->
       *
       * access results with wds->SM_MassFlowRates and/or wds->SM_Pressures
       *
       */
void Staci::Compute_Sensitivity_Matrix(string parameter, int scale) {

  property_ID = parameter;
  int is_edge_prop = -1;  // 0->edge prop, 1-> node prop

  if (strcmp(property_ID.c_str(), "diameter") == 0) is_edge_prop = 0;

  if (strcmp(property_ID.c_str(), "friction_coeff") == 0) is_edge_prop = 0;

  if (strcmp(property_ID.c_str(), "demand") == 0) is_edge_prop = 1;

  if (is_edge_prop == -1) {
    cout << "\n\n ERROR!\n\nCompute_Sensitivity_Matrix(string parameter) -> "
         "unknown parameter: "
         << parameter;
    cout << "\n possible values: diameter|friction_coeff|demand" << endl
         << endl;
    exit(-1);
  }

  SM_row_name.clear();
  SM_col_name.clear();
  SM_MassFlowRates.clear();
  SM_Pressures.clear();

  // Edge property was selected
  if (is_edge_prop == 0) {
    for (unsigned int i = 0; i < agelemek.size(); i++) {
      element_ID = agelemek.at(i)->Get_nev();
      SM_row_name.push_back(agelemek.at(i)->Get_nev());
      Compute_dxdmu();
      // cout<<endl<<"Staci, "<<element_ID<<", "<<property_ID;

      vector<double> tmp1;
      // cout<<endl<<"SM_MassFlowRates:";
      for (unsigned int j = 0; j < agelemek.size(); j++) {
        // printf(" %5.2e",m_dxdmu.at(j));
        tmp1.push_back(m_dxdmu.at(j));
      }
      SM_MassFlowRates.push_back(tmp1);

      vector<double> tmp2;
      // cout<<endl<<"SM_Pressures:";
      for (unsigned int j = 0; j < cspok.size(); j++) {
        // printf(" %7.5e",m_dxdmu.at(agelemek.size() + j));
        tmp2.push_back(m_dxdmu.at(agelemek.size() + j));
      }
      SM_Pressures.push_back(tmp2);
    }
  }

  // Node property was selected
  if (is_edge_prop == 1) {
    for (unsigned int i = 0; i < cspok.size(); i++) {
      element_ID = cspok.at(i)->Get_nev();
      SM_row_name.push_back(cspok.at(i)->Get_nev());
      Compute_dxdmu();

      // cout<<endl<<"!!!!!!!!!!!!"<<endl<<property_ID<<" @ "<<element_ID<<endl;
      // for (unsigned int kkk=0; kkk<agelemek.size(); kkk++)
      //   printf("%+5.3e ",m_dxdmu.at(kkk));
      // cout<<endl;
      // for (unsigned int kkk=0; kkk<cspok.size(); kkk++)
      //   printf("%+5.3e ",m_dxdmu.at(agelemek.size()+kkk));
      // cin.get();

      vector<double> tmp1;
      for (unsigned int j = 0; j < agelemek.size(); j++)
        tmp1.push_back(m_dxdmu.at(j));
      SM_MassFlowRates.push_back(tmp1);

      vector<double> tmp2;
      for (unsigned int j = 0; j < cspok.size(); j++)
        tmp2.push_back(m_dxdmu.at(agelemek.size() + j));
      SM_Pressures.push_back(tmp2);
      if (debug_level < 0 )
        ProgressBar(i, cspok.size());
    }
  }

  // Set up column names
  for (unsigned int j = 0; j < agelemek.size(); j++)
    SM_col_name.push_back(agelemek.at(j)->Get_nev());

  for (unsigned int j = 0; j < cspok.size(); j++)
    SM_col_name.push_back(cspok.at(j)->Get_nev());

  // Compute row-wise and column-wise summitation
  SM_col_sum_MassFlowRates = col_abs_sum(SM_MassFlowRates);
  SM_col_sum_Pressures = col_abs_sum(SM_Pressures);
  SM_col_ss_MassFlowRates = col_sqr_sum(SM_MassFlowRates);
  SM_col_ss_Pressures = col_sqr_sum(SM_Pressures);
  SM_row_sum_MassFlowRates = row_abs_sum(SM_MassFlowRates);
  SM_row_sum_Pressures = row_abs_sum(SM_Pressures);
  SM_row_ss_MassFlowRates = col_sqr_sum(SM_MassFlowRates);
  SM_row_ss_Pressures = col_sqr_sum(SM_Pressures);

  SM_sum_sum_MassFlowRates = abs_sum_sum(SM_MassFlowRates);
  SM_sum_sum_Pressures = abs_sum_sum(SM_Pressures);
  SM_ss_con_MassFlowRates = Cond_sum_sum(SM_MassFlowRates, "MassFlowRates");
  SM_ss_con_Pressures = Cond_sum_sum(SM_Pressures, "Pressures");
  SM_arp_MassFlowRates = Avr_rel_pert(SM_MassFlowRates, "MassFlowRates");
  SM_arp_Pressures = Avr_rel_pert(SM_Pressures , "Pressures");

  // if scale == 1, every vector is rescaled by the maximum
  if (scale == 1) {
    rescale_vec(SM_row_sum_MassFlowRates);
    rescale_vec(SM_row_sum_Pressures);
    rescale_vec(SM_col_sum_MassFlowRates);
    rescale_vec(SM_col_sum_Pressures);
  }
}

//--------------------------------------------------------------

void Staci::rescale_vec(vector<double> &vec) {
  double maxval = -1.e5;
  for (unsigned int j = 0; j < vec.size(); j++)
    if (vec.at(j) > maxval) maxval = vec.at(j);

  for (unsigned int j = 0; j < vec.size(); j++) vec.at(j) = vec.at(j) / maxval;
}

//--------------------------------------------------------------

vector<double> Staci::col_abs_sum(vector<vector<double> > M) {
  vector<double> out(M.at(0).size(), 0.);
  for (unsigned int j = 0; j < M.at(0).size(); j++) {
    double sum = 0.;
    for (unsigned int i = 0; i < M.size(); i++) sum += fabs(M.at(i).at(j));
    out.at(j) = sum;
  }
  return out;
}

//--------------------------------------------------------------
vector<double> Staci::col_sqr_sum(vector<vector<double> > M) {
  vector<double> out(M.at(0).size(), 0.);
  for (unsigned int j = 0; j < M.at(0).size(); j++) {
    double sum = 0.;
    for (unsigned int i = 0; i < M.size(); i++) sum += M.at(i).at(j) * M.at(i).at(j);
    out.at(j) = sum;
  }
  return out;
}

//--------------------------------------------------------------
vector<double> Staci::row_abs_sum(vector<vector<double> > M) {
  vector<double> out(M.size(), 0.);
  for (unsigned int i = 0; i < M.size(); i++) {
    double sum = 0.;
    for (unsigned int j = 0; j < M.at(0).size(); j++) sum += fabs(M.at(i).at(j));
    out.at(i) = sum;
  }
  return out;
}

//--------------------------------------------------------------
vector<double> Staci::row_sqr_sum(vector<vector<double> > M) {
  vector<double> out(M.size(), 0.);
  for (unsigned int i = 0; i < M.size(); i++) {
    double sum = 0.;
    for (unsigned int j = 0; j < M.at(0).size(); j++) sum += M.at(i).at(j) * M.at(i).at(j);
    out.at(i) = sum;
  }
  return out;
}

//--------------------------------------------------------------
double Staci::abs_sum_sum(vector<vector<double> > M) {
  double out = 0;
  for (unsigned int i = 0; i < M.size(); i++) {
    for (unsigned int j = 0; j < M.at(i).size(); j++) {
      out += abs(M.at(i).at(j));
    }
  }
  return out;
}

//-------------------------------------------------------------
// WR conditionally sum sum of the matrix
double Staci::Cond_sum_sum(vector<vector<double> > M, string par_name) {
  double out = 0;
  int idx = 0, n = M.size();
  double maxCons = GetMaxConsumption(idx);
  double treshold, perturb;

  if (par_name == "MassFlowRates") {
    treshold = .1 * GetMaxAbsFlowRate(idx);
    perturb = .1 * maxCons;
  } else {
    if (par_name != "Pressures") {
      cout << " !!!!! WARNING !!!!!! " << endl << " The variable par_name not known, it must be either MassFlowRate or Pressure " << endl << " Switching to Pressure" << endl;
    }
    treshold = 1;
    perturb = .5 * maxCons;
  }

  for (unsigned int i = 0; i < n; i++) {
    for (unsigned int j = 0; j < M.at(i).size(); j++) {
      if (abs(perturb * M.at(i).at(j)) > treshold) {
        out++;
      }
    }
  }
  out = out / (n * M.at(0).size());
  return out;
}

//--------------------------------------------------------------
// WR the average relative perturbation that will cause 1m change in pressure
double Staci::Avr_rel_pert(vector<vector<double> > M, string par_name) {
  int idx = 0, n = M.size(), m = M.at(0).size();
  double maxCons = GetMaxConsumption(idx);
  double treshold, perturb = 0;

  if (par_name == "MassFlowRates") {
    treshold = .1 * GetMaxAbsFlowRate(idx);
  } else {
    if (par_name != "Pressures") {
      cout << " !!!!! WARNING !!!!!! " << endl << " The variable par_name not known, it must be either MassFlowRate or Pressure " << endl << " Switching to Pressure" << endl;
    }
    treshold = 1;
  }

  for (unsigned int i = 0; i < n; i++) {
    for (unsigned int j = 0; j < M.at(i).size(); j++) {
      if (abs(M.at(i).at(j)) > 1e-5) {
        perturb += treshold / abs(M.at(i).at(j));
      }
    }
  }
  perturb = perturb / n / m;
  return perturb;
}

//--------------------------------------------------------------
// WR Saving matrix to a new file with name of "filename"
void Staci::Save_matrix(vector<vector<double> > M, string filename) {
  int n = M.size();
  ofstream file;
  file.open(filename.c_str());
  for (unsigned int i = 0; i < n; i++) {
    for (unsigned int j = 0; j < M.at(i).size(); j++) {
      file << M.at(i).at(j) << ',';
    }
    file << '\n';
  }
}
//--------------------------------------------------------------
// WR Saving matrix to az existing file using the file handler "&file"
void Staci::Save_matrix(vector<vector<double> > M, ofstream &file) {
  int n = M.size();
  for (unsigned int i = 0; i < n; i++) {
    for (unsigned int j = 0; j < M.at(i).size(); j++) {
      file << M.at(i).at(j) << ',';
    }
    file << '\n';
  }
}


//--------------------------------------------------------------
// WR Printing matrix to terminal
void Staci::Print_matrix(vector<vector<double> > M) {
  int n = M.size();
  for (unsigned int i = 0; i < n; i++) {
    for (unsigned int j = 0; j < M.at(i).size(); j++) {
      printf("%5.1e  ", M.at(i).at(j));
    }
    printf("\n");
  }
}

//--------------------------------------------------------------
void Staci::solve_residence_time() {
  string max_ID;
  double max_VAL, mean_VAL, VAL_prev = -61, d_VAL = 1000;

  m_ss.str("");
  m_ss
      << "\n\nComputing residence time...\n===================================="
      << endl;
  if (debug_level > 0) {
    logfile_write(m_ss.str(), 1);
  }
  int step = 1;
  int step_max = 100000;

  while ((step < step_max) && (d_VAL > 0.1)) {
    residence_time_step(max_ID, max_VAL, mean_VAL);
    if (debug_level > 0) {
      if ((step % 20) == 0)
        cout << endl
             << "\t step #" << step
             << " max. res. time: " << convert_to_hr_min(max_VAL) << " ("
             << max_ID << "), mean: " << convert_to_hr_min(mean_VAL);
    }
    d_VAL = fabs(mean_VAL - VAL_prev);
    VAL_prev = mean_VAL;
    step++;
  }
  if (debug_level > 0) {
    m_ss.str("");
    m_ss << endl
         << "\t step #" << step
         << " max. res. time: " << convert_to_hr_min(max_VAL) << " (" << max_ID
         << "), mean: " << convert_to_hr_min(mean_VAL);
    logfile_write(m_ss.str(), 1);
    // cout << m_ss.str();
  }
}

//--------------------------------------------------------------
void Staci::residence_time_step(string & max_ID, double & max_VAL,
                                double & mean_VAL) {
  bool transp_debug = false;
  double mv = -1.;
  double sum = 0;
  double TINY_MASS_FLOW_RATE = 1.e-3;
  double TINY_VEL = 0.000001;  // 1m -> 278 hours

  double MAX_TIME = 168 * 3600.;

  // 1. lepes: csomoponti atlagkor meghatarozasa
  double szaml, nevezo, c, m;
  vector<int> temp;
  // int length_be, length_ki;
  for (unsigned int i = 0; i < cspok.size(); i++) {
    // Csak akkor piszkaljuk a csomoponti eletkort, ha nincs betap.
    // Legalabb 1cm3/h betap legyen
    if (cspok.at(i)->Get_fogy() < (-1.e-6 * 1000. / 3600.)) {
      if (transp_debug)
        cout << endl
             << cspok.at(i)->Get_nev() << ": \n\tbetap miatt adott vizkor: "
             << cspok.at(i)->Get_dprop("tt") / 60. << "min";
      // cin.get();
    } else {
      // Tomegarammal sulyozott atlag
      szaml = 0.0;
      nevezo = 0.0;

      // Bemeno agak
      if (transp_debug)
        cout << endl << cspok.at(i)->Get_nev() << ": \n\tbemeno agak: (";
      temp.clear();
      temp = cspok.at(i)->ag_be;
      // length_be = temp.size();
      if (transp_debug) cout << temp.size() << "db)";

      if (temp.size() > 0) {
        for (unsigned int j = 0; j < temp.size(); j++) {
          if (agelemek.at(temp.at(j))->Get_v() >= 0) {
            c = agelemek.at(temp.at(j))->Get_tt_end();
            double ce = agelemek.at(temp.at(j))->Get_tt_start();
            m = agelemek.at(temp.at(j))->Get_mp();
            if (transp_debug)
              cout << " mp= " << m << " t=" << (c / 60)
                   << " min (tt_start = " << (ce / 60) << "min)";

            // 0 tomegaram eseten 0-val szorozzuk a szamlalot
            if (fabs(m) < TINY_MASS_FLOW_RATE) m = TINY_MASS_FLOW_RATE;
            szaml += fabs(m) * c;
            nevezo += fabs(m);
          }
        }
      }

      // Kimeno agak
      if (transp_debug) cout << endl << "\n\tkimeno agak: (";
      temp.clear();
      temp = cspok.at(i)->ag_ki;
      // length_ki = temp.size();
      if (transp_debug) cout << temp.size() << "db)";

      if (temp.size() > 0) {
        for (unsigned int j = 0; j < temp.size(); j++) {
          if (agelemek.at(temp.at(j))->Get_v() <= 0) {
            c = agelemek.at(temp.at(j))->Get_tt_start();
            double ce = agelemek.at(temp.at(j))->Get_tt_end();
            m = agelemek.at(temp.at(j))->Get_mp();
            if (transp_debug)
              cout << " mp= " << m << " t=" << (c / 60)
                   << " min (tt_end = " << (ce / 60) << "min)";

            // 0 tomegaram eseten 0-val szorozzuk a szamlalot
            if (fabs(m) < TINY_MASS_FLOW_RATE) m = TINY_MASS_FLOW_RATE;
            szaml += fabs(m) * c;
            nevezo += fabs(m);
          }
        }
      }
      // A virtualics csomopontokban siman lehet nevezo=0
      if (fabs(nevezo) < TINY_MASS_FLOW_RATE) nevezo = TINY_MASS_FLOW_RATE;

      cspok.at(i)->Set_dprop("tt", szaml / nevezo);
      sum = sum + szaml / nevezo;

      if (transp_debug)
        cout << "\n\t atlagos vizkor= " << cspok.at(i)->Get_dprop("tt") / 60
             << "min";
    }
  }

  for (unsigned int i = 0; i < cspok.size(); i++) {
    if ((cspok.at(i)->Get_dprop("tt")) > MAX_TIME)
      cspok.at(i)->Set_dprop("tt", MAX_TIME);
  }

  mean_VAL = sum / (double(cspok.size()));

  // 2. lepes: csoelemek leptetese

  double v, L, tt_s, tt_e;
  for (int i = 0; i < agelemek.size(); i++) {
    string type = agelemek.at(i)->GetType();
    int cspe_id = agelemek.at(i)->Get_Cspe_Index();
    int cspv_id = agelemek.at(i)->Get_Cspv_Index();

    if ((type.compare("Cso") == 0) || (type.compare("Csatorna") == 0)) {
      L = agelemek.at(i)->Get_dprop("L");
      v = agelemek.at(i)->Get_v();
      if (v >= 0) {
        // cout<<"\n v>=0!!!";
        if (v < TINY_VEL) {
          tt_s = cspok.at(cspe_id)->Get_dprop("tt");
          tt_e = cspok.at(cspv_id)->Get_dprop("tt");
          v = TINY_VEL;
          if (tt_s > tt_e)
            tt_e = tt_s + L / v;
          else
            tt_s = tt_e + L / v;

          if (tt_s > MAX_TIME) tt_s = MAX_TIME;

          if (tt_e > MAX_TIME) tt_e = MAX_TIME;
        } else {
          int cspe_id = agelemek.at(i)->Get_Cspe_Index();
          tt_s = cspok.at(cspe_id)->Get_dprop("tt");
          tt_e = tt_s + L / v;
          if (tt_e > MAX_TIME) tt_e = MAX_TIME;
        }
        // Get maximum value
        if (tt_e > mv) {
          max_VAL = tt_e;
          mv = tt_e;
          max_ID = agelemek.at(i)->Get_nev();
        }
        // }

        if (tt_s < MAX_TIME)
          agelemek.at(i)->Set_tt_start(tt_s);
        else
          agelemek.at(i)->Set_tt_start(MAX_TIME);
        if (tt_e < MAX_TIME)
          agelemek.at(i)->Set_tt_end(tt_e);
        else
          agelemek.at(i)->Set_tt_end(MAX_TIME);

        if (transp_debug)
          cout << endl
               << agelemek.at(i)->Get_nev() << ": v=" << v << "m/s, L=" << L
               << "m, tt_start=" << (tt_s / 60.)
               << "min, tt_end=" << (tt_e / 60.) << "min";
      } else {
        if (v > -TINY_VEL) {
          tt_s = cspok.at(cspe_id)->Get_dprop("tt");
          tt_e = cspok.at(cspv_id)->Get_dprop("tt");
          v = -TINY_VEL;  // Ki volt kommentelve
          if (tt_s > tt_e)
            tt_e = tt_s - L / v;
          else
            tt_s = tt_e - L / v;
          if (tt_s > MAX_TIME) tt_s = MAX_TIME;

          if (tt_e > MAX_TIME) tt_e = MAX_TIME;
        } else {
          int cspv_id = agelemek.at(i)->Get_Cspv_Index();
          tt_e = cspok.at(cspv_id)->Get_dprop("tt");
          tt_s = tt_e - L / v;
          if (tt_s > MAX_TIME) tt_s = MAX_TIME;

          if (tt_s > mv) {
            max_VAL = tt_s;
            mv = tt_s;
            max_ID = agelemek.at(i)->Get_nev();
          }
        }

        if (tt_s < MAX_TIME)
          agelemek.at(i)->Set_tt_start(tt_s);
        else
          agelemek.at(i)->Set_tt_start(MAX_TIME);
        if (tt_e < MAX_TIME)
          agelemek.at(i)->Set_tt_end(tt_e);
        else
          agelemek.at(i)->Set_tt_end(MAX_TIME);
      }

      // if (agelemek.at(i)->Get_nev()=="PIPE21"){
      //     cout<<endl<<"PIPE21: v="<<v<<", L/v="<<(L/v/3600)<<",
      //     tt_e="<<tt_e/3600<<", tt_s="<<tt_s/3600;
      // }
      //             if (agelemek.at(i)->Get_nev()=="PIPE24"){
      //     cout<<endl<<"PIPE24: v="<<v<<", L/v="<<(L/v/3600)<<",
      //     tt_e="<<tt_e/3600<<", tt_s="<<tt_s/3600;
      // }

      // cin.get();
    } else {
      double v = agelemek.at(i)->Get_v();
      double mp = agelemek.at(i)->Get_mp();

      // cout << endl << "!!!!! "<<agelemek.at(i)->Get_nev()<< " type:
      // "<<type<<endl;

      if ((type.compare("KonstNyomas") == 0) ||
          (type.compare("Vegakna") == 0)) {
        if (mp < 0) {
          // Ez esetben a rendszerbe befele aramlik a kozeg, nem bantjuk a
          // vizkort
          tt_s = agelemek.at(i)->Get_tt_start();
          // tt_e = agelemek.at(i)->Get_tt_end();

          if (transp_debug) {
            cout << endl
                 << agelemek.at(i)->Get_nev()
                 << ": \n\t flow TO the system FROM the element (v=" << v
                 << "m/s) -> prescribed age: tt_start = " << tt_s / 60.
                 << "min";
            //<< ", tt_end = " << tte / 60. << "min";
            // cin.get();
          }
        } else {
          if (transp_debug) {
            tt_s = cspok.at(cspe_id)->Get_dprop("tt");
            agelemek.at(i)->Set_tt_start(tt_s);
            cout << endl
                 << agelemek.at(i)->Get_nev()
                 << ": \n\t flow FROM the system TO the element (v=" << v
                 << "m/s) -> age: tt_start = " << tt_s / 60. << "min";
            // cin.get();
          }
        }
      } else {
        double tt;

        if (fabs(v) < TINY_VEL) {
          double tt1 = cspok.at(cspe_id)->Get_dprop("tt");
          double tt2 = cspok.at(cspv_id)->Get_dprop("tt");
          if (tt1 > tt2)
            tt = tt1;
          else
            tt = tt2;
        } else {
          if (v > 0)
            tt = cspok.at(cspe_id)->Get_dprop("tt");
          else {
            if (cspv_id > -1)
              tt = cspok.at(cspv_id)->Get_dprop("tt");
            else
              tt = cspok.at(cspe_id)->Get_dprop("tt");
          }
        }

        if (tt > MAX_TIME) tt = MAX_TIME;
        agelemek.at(i)->Set_tt_start(tt);
        agelemek.at(i)->Set_tt_end(tt);
        if (transp_debug) {
          cout << endl
               << agelemek.at(i)->Get_nev() << ": \n\tvizkor: " << tt / 60.
               << "min";
          // cin.get();
        }
      }
    }
  }
}

string Staci::convert_to_hr_min(double input_seconds) {
  stringstream str;
  int seconds = round(input_seconds);
  int days = (int)seconds / 60 / 60 / 24;
  int hours = (int)(seconds / 60 / 60) % 24;
  int minutes = (int)(seconds / 60) % 60;

  str << days << " d, " << hours << " h, " << minutes << " m";
  return str.str();
}

double Staci::GetMinPipeDiameter(int &idx) {
  double D, Dmin = 1.e5;
  for (unsigned int i = 0; i < agelemek.size(); i++) {
    if (0 == strcmp((agelemek.at(i)->GetType()).c_str(), "Cso") || 0 == strcmp((agelemek.at(i)->GetType()).c_str(), "Csatorna")) {
      D = agelemek.at(i)->Get_dprop("diameter");
      if (D < Dmin) {
        Dmin = D;
        idx = i;
      }
    }
  }
  return Dmin;
}

double Staci::GetMaxPipeDiameter(int &idx) {
  double D, Dmax = 0.;
  for (unsigned int i = 0; i < agelemek.size(); i++) {
    if (0 == strcmp((agelemek.at(i)->GetType()).c_str(), "Cso") || 0 == strcmp((agelemek.at(i)->GetType()).c_str(), "Csatorna")) {
      D = agelemek.at(i)->Get_dprop("diameter");
      if (D > Dmax) {
        Dmax = D;
        idx = i;
      }
    }
  }
  return Dmax;
}

double Staci::GetMinPipeLength(int &idx) {
  double L, Lmin = 1.e5;
  for (unsigned int i = 0; i < agelemek.size(); i++) {
    if (0 == strcmp((agelemek.at(i)->GetType()).c_str(), "Cso") || 0 == strcmp((agelemek.at(i)->GetType()).c_str(), "Csatorna")) {
      L = agelemek.at(i)->Get_dprop("length");
      if (L < Lmin) {
        Lmin = L;
        idx = i;
      }
    }
  }
  return Lmin;
}

double Staci::GetMaxPipeLength(int &idx) {
  double L, Lmax = 0.;
  for (unsigned int i = 0; i < agelemek.size(); i++) {
    if (0 == strcmp((agelemek.at(i)->GetType()).c_str(), "Cso") || 0 == strcmp((agelemek.at(i)->GetType()).c_str(), "Csatorna")) {
      L = agelemek.at(i)->Get_dprop("length");
      if (L > Lmax) {
        Lmax = L;
        idx = i;
      }
    }
  }
  return Lmax;
}

double Staci::GetSumPipeLength() {
  double Lsum = 0.;
  for (unsigned int i = 0; i < agelemek.size(); i++) {
    if (0 == strcmp((agelemek.at(i)->GetType()).c_str(), "Cso") || 0 == strcmp((agelemek.at(i)->GetType()).c_str(), "Csatorna"))
      Lsum += agelemek.at(i)->Get_dprop("length");
  }
  return Lsum;
}

double Staci::GetSumPipeVolume() {
  double L, D, VolSum = 0.;
  for (unsigned int i = 0; i < agelemek.size(); i++) {
    if (0 == strcmp((agelemek.at(i)->GetType()).c_str(), "Cso") || 0 == strcmp((agelemek.at(i)->GetType()).c_str(), "Csatorna")) {
      L = agelemek.at(i)->Get_dprop("length");
      D = agelemek.at(i)->Get_dprop("diameter");
      VolSum += L * D * D * 3.1416 / 4.;
    }
  }
  return VolSum;
}


double Staci::GetMaxConsumption(int &idx) {
  double c, cmax = -1.e5;
  for (unsigned int i = 0; i < cspok.size(); i++) {
    c = cspok.at(i)->Get_fogy();
    if (c > cmax) {
      cmax = c;
      idx = i;
    }
  }
  return cmax;
}

double Staci::GetMinConsumption(int &idx) {
  double c, cmax = 1.e5;
  for (unsigned int i = 0; i < cspok.size(); i++) {
    c = cspok.at(i)->Get_fogy();
    if (c < cmax) {
      cmax = c;
      idx = i;
    }
  }
  return cmax;
}

double Staci::GetMinGeoHeight(int &idx) {
  double c, cmax = 1.e5;
  for (unsigned int i = 0; i < cspok.size(); i++) {
    c = cspok.at(i)->Get_dprop("height");
    if (c < cmax) {
      cmax = c;
      idx = i;
    }
  }
  return cmax;
}

double Staci::GetMaxGeoHeight(int &idx) {
  double c, cmax = -1.e5;
  for (unsigned int i = 0; i < cspok.size(); i++) {
    c = cspok.at(i)->Get_dprop("height");
    if (c > cmax) {
      cmax = c;
      idx = i;
    }
  }
  return cmax;
}

double Staci::GetMaxAbsFlowRate(int &idx) {
  double c, cmax = -1.e5;
  for (unsigned int i = 0; i < agelemek.size(); i++) {
    c = agelemek.at(i)->Get_Q();
    if (abs(c) > cmax) {
      cmax = c;
      idx = i;
    }
  }
  return abs(cmax);
}

double Staci::GetMinAbsFlowRate(int &idx) {
  double c, cmax = 1.e10;
  for (unsigned int i = 0; i < agelemek.size(); i++) {
    c = agelemek.at(i)->Get_Q();
    if (abs(c) < cmax) {
      cmax = c;
      idx = i;
    }
  }
  return abs(cmax);
}

double Staci::GetMaxPressure(int &idx) {
  double c, cmax = -1.e5;
  for (unsigned int i = 0; i < cspok.size(); i++) {
    c = cspok.at(i)->Get_p();
    if (c > cmax) {
      cmax = c;
      idx = i;
    }
  }
  return cmax;
}

double Staci::GetMinPressure(int &idx) {
  double c, cmax = 1.e10;
  for (unsigned int i = 0; i < cspok.size(); i++) {
    c = cspok.at(i)->Get_p();
    if (c < cmax) {
      cmax = c;
      idx = i;
    }
  }
  return cmax;
}


void Staci::Statistics() {
  int idx;
  m_ss.str("");
  m_ss << endl << endl << "Statistics:";
  m_ss << endl << "======================";

  m_ss << endl << "\t # of nodes : " << agelemek.size();
  m_ss << endl << "\t # of edges : " << cspok.size() << endl;

  m_ss << endl << "\t diameter     : ";
  m_ss << "min : " << GetMinPipeDiameter(idx) << "m (ID: " << agelemek.at(idx)->Get_nev() << ")  ";
  m_ss << "max : " << GetMaxPipeDiameter(idx) << "m (ID: " << agelemek.at(idx)->Get_nev() << ")";

  m_ss << endl << "\t length       : ";
  m_ss << "min : " << GetMinPipeLength(idx) << "m (ID: " << agelemek.at(idx)->Get_nev() << ")  ";
  m_ss << "max : " << GetMaxPipeLength(idx) << "m (ID: " << agelemek.at(idx)->Get_nev() << ")";

  m_ss << endl << "\t sum. length  : " << GetSumPipeLength() << " m ";
  m_ss << endl << "\t pipe volume  : " << GetSumPipeVolume() << " m^3 ";
  m_ss << endl << "\t fluid volume : " << FolyMenny << " m^3 (" << FolyMenny / GetSumPipeVolume() * 100. << "% of pipe volume)" << endl;

  m_ss << endl << "\t max. inflow          : " << -3.6 * GetMinConsumption(idx) << " m3/h (ID: " << cspok.at(idx)->Get_nev() << ")  ";
  m_ss << endl << "\t max. consumption     : " << 3.6 * GetMaxConsumption(idx) << " m3/h (ID: " << cspok.at(idx)->Get_nev() << ")  ";
  m_ss << endl << "\t sum of inflow        : " << 3.6 * get_sum_of_neg_consumption() << " m3/h";
  m_ss << endl << "\t sum of consumption   : " << 3.6 * get_sum_of_pos_consumption() << " m3/h";
  m_ss << endl << "\t consumption - inflow : " << 3.6 * get_sum_of_consumption() << " m3/h" << endl;

  m_ss << endl << "\t node geodetic height : ";
  m_ss << "min : " << GetMinGeoHeight(idx) << "m (ID: " << cspok.at(idx)->Get_nev() << ")  ";
  m_ss << "max : " << GetMaxGeoHeight(idx) << "m (ID: " << cspok.at(idx)->Get_nev() << ")" << endl;

  int N = cspok.size() + agelemek.size();
  m_ss << endl << " \t number of nonzero Jacobian entries: " << m_nnz << " out of " << (N * N);
  m_ss << " (" << ((((double)m_nnz) / N / N) * 100) << "%)" << endl;


  logfile_write(m_ss.str(), 1);
}

// WR Computing rank of a node
void Staci::Compute_rank()
{
  for (int i = 0; i < agelemek.size(); i++) {

    int cspe_idx = agelemek.at(i)->Get_Cspe_Index();
    cspok.at(cspe_idx)->Inc_rank();

    if (agelemek.at(i)->Get_Csp_db() == 2) {
      int cspv_idx = agelemek.at(i)->Get_Cspv_Index();
      cspok.at(cspv_idx)->Inc_rank();
    }
  }
}

// WR Progress bar for long calcs, n: no. of steps, i: current step
void Staci::ProgressBar(int i, int n) {
  if (i + 1 == n) i++;
  int nn;
  if (n > 100) nn = round(n / 100.);
  nn = 1;
  if (remainder(i, nn) == 0) {
    cout << "\x1b[A" << "[*] Wait for it ... [";
    for (int j = 0; j < round((double)i / n * 20); j++)
      cout << "#";
    for (int j = round((double)i / n * 20); j < 20; j++)
      cout << " ";
    cout << "]" << " ... ";
    printf("%3.0f / 100\n", round((double)i / n * 100));
  }
}

// WR Calculates average, max, standard deviation of a vector x
void Staci::Avr_absmax_stddev(vector<double> x, double &a, double &m, double &s) {
  a = 0; m = 0; s = 0;
  int n = x.size();
  for (int i = 0; i < n; i++) {
    a += x.at(i);
    if (abs(x.at(i)) > m) m = abs(x.at(i));
  }
  a /= n;
  for (int i = 0; i < n; i++)
    s += pow(x.at(i) - a, 2);
  s = pow(s / (n - 1), 0.5);
}

// WR Reading doubles from file, separeted with ','
vector<vector<double> > Staci::CSVRead(ifstream &file) {
  string temp, line;
  int k = 0;
  vector<vector<double> > v;
  if (file.is_open()) {
    while (getline(file, line)) {
      for (string::iterator j = line.begin(); j != line.end(); j++) {
        if (*j != ',')
          temp += *j;
        else {
          //v.at(k).push_back(atof(temp.c_str()));
          v.at(k).push_back(stod(temp, 0)); // not working WHY? TODO
          k++;
          temp = "";
        }
      }
      k = 0;
    }
  } else
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
  cout << "Staci::CSVRead(), File is not open when calling CSVRead function!!!" << endl;
  return v;
}

void Staci::compute_demand_sensitivity() {
  Compute_Sensitivity_Matrix("demand", 0);
  Save_Sensitivity_Matrix("SM.txt");

  for (unsigned int j = 0; j < agelemek.size(); j++)
    agelemek.at(j)->Set_user1(SM_col_sum_MassFlowRates.at(j));
  
  for (unsigned int j = 0; j < cspok.size(); j++)
    cspok.at(j)->Set_user1(SM_col_sum_Pressures.at(j));

}

void Staci::Save_Sensitivity_Matrix(string fname) {
  FILE * pFile;

  int n_edges = agelemek.size();
  int n_nodes = cspok.size();
  double max_SM_MassFlowRates = GetAbsMaxCoeff(SM_MassFlowRates);
  double max_SM_Pressures  = GetAbsMaxCoeff(SM_Pressures);

  string par = "demand";
  pFile = fopen (fname.c_str(), "w");
  fprintf(pFile, "parameter: %s\n", par.c_str());
  fprintf(pFile, "n_edges : %d\n", n_edges);
  fprintf(pFile, "n_nodes : %d\n", n_nodes);
  fprintf(pFile, "n_cols  : %d ( = n_edges + n_nodes )\n", n_edges + n_nodes);
  fprintf(pFile, "n_rows  : %d ( = n_edges )\n", n_edges);
  fprintf(pFile, "Note that these are scaled sensitivities:");
  fprintf(pFile, " max_SM_MassFlowRates = %5.3e, max_SM_Pressures = %5.3e\n\n", max_SM_MassFlowRates, max_SM_Pressures);
  fprintf(pFile, "column names; ");
  for (int j = 0 ; j < n_edges; j++)
    fprintf (pFile, "%s; ", agelemek.at(j)->Get_nev().c_str());
  for (int j = 0 ; j < n_nodes; j++)
    fprintf (pFile, "%s; ", cspok.at(j)->Get_nev().c_str());
  fprintf (pFile, "\n");

  for (int i = 0 ; i < n_nodes; i++) {
    fprintf(pFile, "%s @ %s; ", par.c_str(), cspok.at(i)->Get_nev().c_str());
    for (int j = 0 ; j < n_edges; j++)
      fprintf (pFile, "%+7.5e; ", SM_MassFlowRates.at(i).at(j));
    for (int j = 0 ; j < n_nodes; j++)
      fprintf (pFile, "%+7.5e; ", SM_Pressures.at(i).at(j));
    fprintf (pFile, "\n");
  }

  // Add total sensitivities
  vector<pair<double, string> > v_edges;
  vector<pair<double, string> > v_nodes;

  // cout << "\n n_edges = " << n_edges << ", n_nodes=" << n_nodes << ", sum: " << (n_edges + n_nodes);
  // cout << "\n size of SM_MassFlowRates : " << SM_MassFlowRates.size() << " x " << SM_MassFlowRates.at(0).size();
  // cout << "\n size of SM_Pressures  : " << SM_Pressures.size() << " x " << SM_Pressures.at(0).size();
  // cin.get();


  for (int i = 0 ; i < n_edges; i++) {
    double tmp = 0.;
    for (int j = 0 ; j < n_nodes; j++)
      tmp += SM_MassFlowRates.at(j).at(i) * SM_MassFlowRates.at(j).at(i) ;
    v_edges.push_back(pair<double, string>( sqrt(tmp), agelemek.at(i)->Get_nev() ));
  }
  for (int i = 0 ; i < n_nodes; i++) {
    double tmp = 0.;
    for (int j = 0 ; j < n_nodes; j++)
      tmp += SM_Pressures.at(j).at(i) * SM_Pressures.at(j).at(i) ;
    v_nodes.push_back(pair<double, string>( sqrt(tmp), cspok.at(i)->Get_nev() ));
  }


  // for (int k = 0 ; k < v_edges.size(); k++)
  //     cout << "\n ID = " << v_edges.at(k).ID << ", val=" <<  v_edges.at(k).val;
  // for (int k = 0 ; k < v_nodes.size(); k++)
  //     cout << "\n ID = " << v_nodes.at(k).ID << ", val=" <<  v_nodes.at(k).val;
  // cin.get();

  // fprintf (pFile, "\nColumn-wise sum of absolute values:");
  // for (int i = 0; i < n_edges; i++)
  //     fprintf(pFile, "\n %10s: %7.5e", v_edges.at(i).ID.c_str(), v_edges.at(i).val);
  // fprintf(pFile, "\n");
  // for (int i = 0; i < n_nodes; i++)
  //     fprintf(pFile, "\n %10s: %7.5e", v_nodes.at(i).ID.c_str(), v_nodes.at(i).val);

  sort(v_edges.begin(), v_edges.end());
  sort(v_nodes.begin(), v_nodes.end());

  fprintf (pFile, "\nColumn-wise L2 norm of sensitivities, after sorting:");
  for (int i = 0; i < n_edges; i++)
    fprintf(pFile, "\n %10s; %7.5e", v_edges.at(i).second.c_str(), v_edges.at(i).first);
  fprintf(pFile, "\n");
  for (int i = 0; i < n_nodes; i++)
    fprintf(pFile, "\n %10s; %7.5e", v_nodes.at(i).second.c_str(), v_nodes.at(i).first);

  fclose (pFile);
}

double Staci::GetAbsMaxCoeff(vector< vector<double> > M) {
  double maxval = -1.;
  for (unsigned int i = 0; i < M.size(); i++) {
    for (unsigned int j = 0; j < M.at(i).size(); j++) {
      if (maxval < abs(M.at(i).at(j)))
        maxval = abs(M.at(i).at(j));
    }
  }
  return maxval;
}
