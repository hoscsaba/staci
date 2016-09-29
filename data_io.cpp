using namespace std;

#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <vector>
#include <math.h>
#include "KonstNyomas.h"
#include "Cso.h"
#include "Szivattyu.h"
#include "JelleggorbesFojtas.h"
#include "Vegakna.h"
#include "Csatorna.h"
#include "BukoMutargy.h"
#include "VisszacsapoSzelep.h"
#include "data_io.h"
#include "xmlParser.h"

// TODO
// 2014.05.07. HCs. Ha a felhasznalo veletlenul kitorli pl. a klor adatot csompontnal, elhasal az egesz! Hibavisszajelzest!


data_io::data_io(const char *a_xml_fnev) {
    xml_fnev = a_xml_fnev;

    // Itt kenytelen vagyok kulon definialni a rednelkezesre allo elemeket
    // Jobb lenne valahol mashol letrehozni oket, de sajna csak itt megy.
    edge_type.push_back("press");
    edge_type.push_back("pipe");
    edge_type.push_back("pump");
    edge_type.push_back("valve");
    edge_type.push_back("pool");
    edge_type.push_back("channel"); // t�glalap km.
    edge_type.push_back("channel1"); // k�r km.
    edge_type.push_back("channel2"); // tetsz�leges km.
    edge_type.push_back("overflow"); // buk� m�t�rgy
    edge_type.push_back("checkvalve"); // visszacsap� szelep
    edge_type_number = edge_type.size();

    for (int i = 0; i < edge_type_number; i++)
        edge_type_occur.push_back(0);

    XMLNode xMainNode = XMLNode::openFileHelper(xml_fnev, "staci");

    string cpp_xml_debug = xMainNode.getChildNode("settings").getChildNode("cpp_xml_debug").getText();
    if (strcmp(cpp_xml_debug.c_str(), "true") == 0)
        debug = true;
    else
        debug = false;

}

//--------------------------------------------------------------------------------
void data_io::load_system(vector<Csomopont *> &cspok, vector<Agelem *> &agelemek) {
    XMLNode xMainNode = XMLNode::openFileHelper(xml_fnev, "staci");

    XMLNode Node_settings = xMainNode.getChildNode("settings");
    XMLNode Node_nodes = xMainNode.getChildNode("nodes");
    XMLNode Node_edges = xMainNode.getChildNode("edges");

    // Csomopontok es agak szamanak es nevenek kiolvasasa...
    csp_db = Node_nodes.nChildNode("node");
    ag_db = Node_edges.nChildNode("edge");
    if (debug)
        cout << endl << "Adatfajl: " << xml_fnev << endl << "\tcsomopontok szama: "
             << csp_db << endl << "\tagak szama:        " << ag_db;

    // Csomopontok reszletes kiolvasasa...
    if (debug)
        cout << endl << endl << "A csomopontok alapadatai reszletesen: " << endl
             << "-----------------------------------------" << endl;

    string id, is_endnode;
    double height, demand, cl_be, pressure, ro, tt;
    int j = 0;
    for (int i = 0; i < csp_db; i++) {
        is_endnode = Node_nodes.getChildNode("node", i).getChildNode("endnode").getText();
        //cout << endl << " NODE #" << i << "/" << csp_db << ": endnode=" << is_endnode << "   " << (is_endnode == "false");

        if (is_endnode == "false") {
            id = Node_nodes.getChildNode("node", i).getChildNode("id").getText();

            height = string_to_double(Node_nodes.getChildNode("node", i).getChildNode("height").getText(), id,
                                      "<height>", 0.0);

            demand = string_to_double(Node_nodes.getChildNode("node", i).getChildNode("demand").getText(), id,
                                      "<demand>", 0.0);

            cl_be = string_to_double(Node_nodes.getChildNode("node", i).getChildNode("cl_input").getText(), id,
                                     "<cl_input>", 0.0);

            pressure = string_to_double(Node_nodes.getChildNode("node", i).getChildNode("pressure").getText(), id,
                                        "<pressure>", 0.0);

            ro = string_to_double(Node_nodes.getChildNode("node", i).getChildNode("density").getText(), id,
                                  "<pressure>", 0.0);

            tt = string_to_double(Node_nodes.getChildNode("node", i).getChildNode("travel_time").getText(), id,
                                  "<travel_time>", 0.0);

            cspok.push_back(new Csomopont(id, height, demand, cl_be, pressure, ro, tt));
//            if (debug)
//                cout << cspok.at(j)->Info();
            j++;
        }
    }
    csp_db = j--;

    // Agak reszletes kiolvasasa...
    if (debug)
        cout << endl << endl << "Az agak: " << endl
             << "-----------------------------------------" << endl;
    // Eloszor csak az info kedveert megszamoljuk a elemszamot
    int db;
    bool megvan;
    int ezaz = 0;
    if (debug) {
        for (int i = 0; i < ag_db; i++) {
            megvan = false;
            db = 0;
            for (int j = 0; j < edge_type_number; j++) {
                db = Node_edges.getChildNode("edge", i).getChildNode("edge_spec", 0).nChildNode(
                        edge_type.at(j).c_str());
                if (db == 1) {
                    edge_type_occur.at(j) += 1;
                    megvan = true;
                    ezaz = j;
                }
            }
            //cout<<endl<<"\t"<<(i+1)<<". edge: ";
            //if (megvan)
            //  cout<<edge_type.at(ezaz);
            //else
            //  cout<<"????";
        }
    }

    if (debug)
        cout << endl << endl << "Az agak reszletesen:" << endl
             << "-----------------------------------------" << endl;
    string node_from, node_to;
    double density, patlag = 0, pmax = 0, aref, mass_flow_rate, travel_time;
    int darab = 0;
    for (int i = 0; i < ag_db; i++) {
        id = Node_edges.getChildNode("edge", i).getChildNode("id").getText();
        /*cout<<endl<<" id:"<<id;
        cin.get();
        */
        aref = atof(Node_edges.getChildNode("edge", i).getChildNode("aref").getText());
        if (aref < 1e-3 * 1e-3) {
            printf("\nWARNING!!! Reference area of %s is %g, overriding with 1. m^2.", id.c_str(), aref);
            aref = 1.;
        }
        /*cout<<endl<<" aref:"<<aref;
        cin.get();
        */
        node_from = Node_edges.getChildNode("edge", i).getChildNode("node_from").getText();
        /*cout<<endl<<" node_from:"<<node_from;
        cin.get();
        */
        node_to = Node_edges.getChildNode("edge", i).getChildNode("node_to").getText();
        /*cout<<endl<<" node_to:"<<node_to;
        cin.get();
        */
        density = atof(Node_edges.getChildNode("edge", i).getChildNode("density").getText());
        /*cout<<endl<<" density:"<<density;
        cin.get();
        */
        mass_flow_rate = atof(Node_edges.getChildNode("edge", i).getChildNode("mass_flow_rate").getText());
        /*cout<<endl<<" mass_flow_rate:"<<mass_flow_rate;
        cin.get();
        */
        travel_time = atof(Node_edges.getChildNode("edge", i).getChildNode("travel_time").getText());
        if (debug)
            cout << endl << (i + 1) << ". edge: " << id << ", node_from=" << node_from
                 << ", node_to=" << node_to << ", density=" << density << ", aref="
                 << aref;

        for (int j = 0; j < edge_type_number; j++) {
            db = Node_edges.getChildNode("edge", i).getChildNode("edge_spec").nChildNode(edge_type.at(j).c_str());
            if (db == 1) {
                XMLNode elem = Node_edges.getChildNode("edge", i).getChildNode("edge_spec").getChildNode(
                        edge_type.at(j).c_str());
                switch (j) {
                    case 0: //press = Konst_Nyomas
                    {
                        double pres = atof(elem.getChildNode("pressure").getText());
                        if (debug)
                            cout << ", pressure=" << pres << "Pa, szumma_patlag=" << patlag
                                                                                     / 1000 / 9.81 << "vom";
                        patlag += pres;
                        if (pres > pmax)
                            pmax = pres;
                        darab++;
                        agelemek.push_back(
                                new KonstNyomas(id, aref, node_from, density, pres, mass_flow_rate, travel_time));
                        if (debug)
                            cout << " OK";
                        break;
                    }

                    case 1: //pipe = Cso
                    {
                        double L = atof(elem.getChildNode("length").getText());
                        double D = atof(elem.getChildNode("diameter").getText());
                        double erdesseg = atof(elem.getChildNode("roughness").getText());
                        double cl_k = atof(elem.getChildNode("cl_k").getText());
                        double cl_w = atof(elem.getChildNode("cl_w").getText());
                        //double b=atof(elem.getChildNode("distr_consumption").getText());
                        if (debug)
                            cout << ", length=" << L << "m, diameter=" << D
                                 << "m, roughness=" << erdesseg << "mm";
                        agelemek.push_back(
                                new Cso(id, node_from, node_to, density, L, D, erdesseg, cl_k, cl_w, mass_flow_rate));
                        if (debug)
                            cout << " OK";
                        break;
                    }

                    case 2: // pump = Szivattyu
                    {
                        int jgpsz = elem.getChildNode("curve").getChildNode("points").nChildNode("point_x");
                        vector<double> Q(jgpsz), H(jgpsz);
                        curve_reader(id, elem.getChildNode("curve"), Q, H);
                        agelemek.push_back(new Szivattyu(id, node_from, node_to, density, aref, Q, H, mass_flow_rate));
                        if (debug)
                            cout << " OK";
                        break;
                    }

                    case 3: // Valve = JelleggorbesFojtas
                    {
                        double allas = atof(elem.getChildNode("position").getText());

                        if (debug)
                            cout << ", actual setting=" << allas;
                        int jgpsz = elem.getChildNode("curve").getChildNode("points").nChildNode("point_x");
                        vector<double> e(jgpsz), zeta(jgpsz);
                        curve_reader(id, elem.getChildNode("curve"), e, zeta);
                        if (aref < 1e-5) {
                            double Aref_min = 3.1416e-04; // 2cm-es cso belvilaga
                            aref = Aref_min;
                            cout << "Warning! element " << id << " Aref=" << aref << ", overwriting with " << Aref_min
                                 << endl;
                        }
                        if (allas > e.at(e.size() - 1))
                            cout << "Warning! element " << id << " actual setting=" << allas << " > e(end) "
                                 << e.at(e.size() - 1) << endl;
                        if (allas < e.at(0))
                            cout << "Warning! element " << id << " actual setting=" << allas << " < e(0) " << e.at(0)
                                 << endl;

                        agelemek.push_back(new JelleggorbesFojtas(id, node_from, node_to, density, aref, e, zeta, allas,
                                                                  mass_flow_rate));
                        if (debug)
                            cout << " OK";
                        break;
                    }

                    case 4: //pool = Vegakna
                    {
                        double Hb = atof(elem.getChildNode("bottom_level").getText());
                        double Hw = atof(elem.getChildNode("water_level").getText());
                        patlag += Hw + Hb;
                        if ((Hw + Hb) > pmax)
                            pmax = Hw + Hb;
                        darab++;
                        if (debug)
                            cout << ", bottom_level=" << Hb << "m" << ", water_level=" << Hw
                                 << "m, szumma_patlag=" << patlag;

                        agelemek.push_back(
                                new Vegakna(id, node_from, density, aref, Hb, Hw, mass_flow_rate, travel_time));
                        if (debug)
                            cout << " OK";
                        break;
                    }

                        //              case 5: //channel0 = T�glalap km. csatorna
                        //              {
                        //                  double L = atof(elem.getChildNode("length").getText());
                        //                  double ze = atof(elem.getChildNode("start_height").getText());
                        //                  double zv = atof(elem.getChildNode("end_height").getText());
                        //                  //      double inc   = atof(elem.getChildNode("inclination").getText())/100;
                        //                  double erdesseg=atof(elem.getChildNode("roughness").getText());
                        //                  int int_steps= atoi(elem.getChildNode("integral_steps").getText());
                        //                  int debugl = atoi(elem.getChildNode("debug_level").getText());
                        //                  double cl_k=atof(elem.getChildNode("cl_k").getText());
                        //                  double cl_w=atof(elem.getChildNode("cl_w").getText());
                        //                  if (debug)
                        //                      cout<<", length="<<L<<"m"<<", start_height="<<ze
                        //                              <<"m, end_height="<<zv<<"m, roughness="
                        //                              <<erdesseg <<", integral_steps="<<int_steps
                        //                              <<", debug_level=" <<debug;
                        //                  double width = atof(elem.getChildNode("width").getText());
                        //                  double Hmax = atof(elem.getChildNode("max_height").getText());
                        //
                        //                  if (debug)
                        //                      cout<<endl<<"\t\t teglalap geometria: B="<<width
                        //                              <<"m, Hmax="<<Hmax;
                        //                  agelemek.push_back(new Csatorna(id,node_from,node_to,aref,L,ze,zv,erdesseg,int_steps,debugl,width,Hmax,cl_k,cl_w));
                        //                  break;
                        //              }

                    case 6: // kor rm.
                    {
                        double L = atof(elem.getChildNode("length").getText());
                        double ze = atof(elem.getChildNode("start_height").getText());
                        double zv = atof(elem.getChildNode("end_height").getText());
                        bool is_reversed = false;
                        if (ze < zv) {
                            cout
                                    << "\n\n\tNegativ lejtes -> CSOMOPONT CSERE!!!\n\t elem neve: "
                                    << id;
                            cout << "\n\t elotte: node_from:" << node_from
                                 << ", node_to:" << node_to << ", ze=" << ze << ", zv="
                                 << zv;

                            string s_tmp = node_from;
                            node_from = node_to;
                            node_to = s_tmp;

                            double d_tmp = ze;
                            ze = zv;
                            zv = d_tmp;
                            cout << "\n\t utana : node_from:" << node_from
                                 << ", node_to:" << node_to << ", ze=" << ze << ", zv="
                                 << zv;
                            is_reversed = true;
                        }
                        double erdesseg = atof(elem.getChildNode("roughness").getText());
                        int int_steps = atoi(elem.getChildNode("integral_steps").getText());
                        int debugl = atoi(elem.getChildNode("debug_level").getText());
                        double cl_k = atof(elem.getChildNode("cl_k").getText());
                        double cl_w = atof(elem.getChildNode("cl_w").getText());
                        if (debug)
                            cout << ", length=" << L << "m" << ", start_height=" << ze
                                 << "m, end_height=" << zv << "m, roughness="
                                 << erdesseg << ", integral_steps=" << int_steps
                                 << ", debug_level=" << debug;
                        double dia = atof(elem.getChildNode("diameter").getText());

                        if (debug)
                            cout << endl << "\t\t kor geometria: D=" << dia;
                        agelemek.push_back(
                                new Csatorna(id, node_from, node_to, density, dia * dia * M_PI / 4., L, ze, zv,
                                             erdesseg, int_steps, debugl, dia, cl_k, cl_w, is_reversed,
                                             mass_flow_rate));
                        break;
                    }

                        //              case 7: // Felhasznalo altal definialt tipus
                        //              {
                        //                  double L = atof(elem.getChildNode("length").getText());
                        //                  double ze = atof(elem.getChildNode("start_height").getText());
                        //                  double zv = atof(elem.getChildNode("end_height").getText());
                        //                  //                        double inc   = atof(elem.getChildNode("inclination").getText())/100;
                        //                  double erdesseg=atof(elem.getChildNode("roughness").getText());
                        //                  int int_steps= atoi(elem.getChildNode("integral_steps").getText());
                        //                  int debugl = atoi(elem.getChildNode("debug_level").getText());
                        //                  double cl_k=atof(elem.getChildNode("cl_k").getText());
                        //                  double cl_w=atof(elem.getChildNode("cl_w").getText());
                        //                  if (debug)
                        //                      cout<<", length="<<L<<"m"<<", start_height="<<ze
                        //                              <<"m, end_height="<<zv<<"m, roughness="
                        //                              <<erdesseg <<", integral_steps="<<int_steps
                        //                              <<", debug_level=" <<debug;
                        //                  int jgpszb, jgpsza, jgpszk;
                        //                  for (int cdb=0; cdb<elem.nChildNode("curve"); cdb++) {
                        //                      if (elem.getChildNode("curve",cdb).getChildNode("id").getText()=="curve_b")
                        //                          jgpszb=elem.getChildNode("curve").getChildNode("points").nChildNode("point_x");
                        //                      if (elem.getChildNode("curve",cdb).getChildNode("id").getText()=="curve_a")
                        //                          jgpsza=elem.getChildNode("curve").getChildNode("points").nChildNode("point_x");
                        //                      if (elem.getChildNode("curve",cdb).getChildNode("id").getText()=="curve_k")
                        //                          jgpszk=elem.getChildNode("curve").getChildNode("points").nChildNode("point_x");
                        //                  }
                        //                  vector<double> yb(jgpszb), b(jgpszb), ya(jgpsza), a(jgpsza),
                        //                          yk(jgpszk), k(jgpszk);
                        //
                        //                  for (int cdb=0; cdb<elem.nChildNode("curve"); cdb++) {
                        //                      if (elem.getChildNode("curve",cdb).getChildNode("id").getText()=="curve_b")
                        //                          curve_reader(id, elem.getChildNode("curve"), yb, b);
                        //                      if (elem.getChildNode("curve",cdb).getChildNode("id").getText()=="curve_a")
                        //                          curve_reader(id, elem.getChildNode("curve"), ya, a);
                        //                      if (elem.getChildNode("curve",cdb).getChildNode("id").getText()=="curve_k")
                        //                          curve_reader(id, elem.getChildNode("curve"), yk, k);
                        //                  }
                        //                  agelemek.push_back(new Csatorna(id,node_from,node_to,aref,L,ze,zv,erdesseg,int_steps,debugl,yb,b,ya,a,yk,k,cl_k,cl_w));
                        //                  break;
                        //              }

                    case 8: // Buko mutargy
                    {
                        string iso = elem.getChildNode("is_opened").getText();
                        double bh = atof(elem.getChildNode("bottom_height").getText());
                        bool is_opened = false;
                        if (iso == "true")
                            is_opened = true;
                        double wi = atof(elem.getChildNode("width").getText());
                        double oh = atof(elem.getChildNode("overflow_height").getText());
                        double dc = atof(elem.getChildNode("discharge_coeff").getText());
                        double vc = atof(elem.getChildNode("valve_coeff").getText());
                        if (debug)
                            cout << endl << "is_opened:" << iso << ", width=" << wi
                                 << ", overflow_height=" << oh
                                 << ", discharge_coeff=" << dc << ", valve_coeff="
                                 << vc;
                        agelemek.push_back(
                                new BukoMutargy(id, node_from, node_to, density, aref, bh, is_opened, wi, oh, dc, vc,
                                                mass_flow_rate));
                        if (debug)
                            cout << "  OK";
                        break;
                    }
                    case 9: // visszacsap� szelep
                    {
                        double lce = atof(elem.getChildNode("loss_coeff_f").getText()) / 1e5;
                        double lcv = atof(elem.getChildNode("loss_coeff_b").getText()) / 1e5;
                        if (debug) {
                            cout << endl << "loss_coeff_f(orward):" << lce;
                            cout << endl << "loss_coeff_b(ack):" << lcv;
                        }
                        agelemek.push_back(
                                new VisszacsapoSzelep(id, node_from, node_to, density, aref, lce, lcv, mass_flow_rate));
                        if (debug)
                            cout << "  OK";
                        break;
                    }
                    default: {
                        cout << endl << endl
                             << "xml adatfajl feldolgozasa: HIBA!!! Sz@r van, nem talalok ilyen tipusu agelemet => "
                             << i << ". edge, id:" << id << endl;
                    }
                }
            }
        }
/*        printf("\n Size of agelemek is %lu, last one added is %s.",agelemek.size(),agelemek.at(agelemek.size()-1)->Get_nev().c_str());
        cout<<agelemek.at(agelemek.size()-1)->Info();
        cin.get();*/
    }
    if (debug) {
        cout << endl << endl << endl << "Number of edge types:";
        cout << endl << "\tnode:\t" << csp_db;
        for (int j = 0; j < edge_type_number; j++)
            cout << endl << "\t" << edge_type.at(j) << ":\t" << edge_type_occur.at(j);
    }

    //patlag = patlag / darab;

}

//--------------------------------------------------------------------------------
void data_io::curve_reader(const string id, const XMLNode elem,
                           vector<double> &px, vector<double> &py) {

    string curve_id = elem.getChildNode("id").getText();
    string x_val = elem.getChildNode("x_val").getText();
    string y_val = elem.getChildNode("y_val").getText();
    string x_dim = elem.getChildNode("x_dim").getText();
    string y_dim = elem.getChildNode("y_dim").getText();
    int xjgpsz = elem.getChildNode("points").nChildNode("point_x");
    int yjgpsz = elem.getChildNode("points").nChildNode("point_y");

    string tmp;

    if (debug)
        cout << endl << "\t curve_id: " << curve_id << "  x_jgsz=" << xjgpsz
             << ", y_jgpsz=" << yjgpsz;

    if (xjgpsz == yjgpsz) {
        if (debug)
            cout << ", OK, adatok beolvasasa indulhat" << endl;
        double x, y;
        for (int i = 0; i < xjgpsz; i++) {
            tmp = elem.getChildNode("points").getChildNode("point_x", i).getText();
            replace(tmp.begin(), tmp.end(), ',', '.');
            x = atof(tmp.c_str());

            tmp = elem.getChildNode("points").getChildNode("point_y", i).getText();
            replace(tmp.begin(), tmp.end(), ',', '.');
            y = atof(tmp.c_str());

            px.at(i) = x;
            py.at(i) = y;
            if (debug)
                cout << "\t\t " << x_val << "[" << x_dim << "]=" << x << "  " << y_val << "["
                     << y_dim << "]=" << y << endl;
        }
    } else
        cout << endl << "xml adatfajl feldolgozasa: HIBA! jgpsz_x=" << xjgpsz
             << " <-> jgpsz_y=" << yjgpsz << endl << "Atugrom a(z) " << id
             << " elemet..." << endl << endl;

}

//--------------------------------------------------------------------------------
void data_io::save_results(double FolyMenny, vector<Csomopont *> cspok,
                           vector<Agelem *> agelemek, bool conv_reached) {


    XMLNode xMainNode = XMLNode::openFileHelper(xml_fnev, "staci");

    cout << endl << endl << "Output file name:" << xml_fnev << ", saving results...";

    XMLNode Node_settings = xMainNode.getChildNode("settings");
    XMLNode Node_nodes = xMainNode.getChildNode("nodes");
    XMLNode Node_edges = xMainNode.getChildNode("edges");

    /*debug = false;
    string debug_string = Node_settings.getChildNode("cpp_xml_debug").getText();

    if (strcmp(debug_string,"true"))
        debug=true;

    if (debug)
        cout << endl << endl << "cpp_xml_debug=" << debug_string
             << " => adatfajl irasa kozbeni debug bekapcsolasa..." << endl;*/

    Node_settings.getChildNode("solution_exists").deleteText();
    Node_settings.getChildNode("fluid_volume").deleteText();
    if (conv_reached) {
        Node_settings.getChildNode("solution_exists").addText("true");
        ostringstream os;
        os << scientific << setprecision(5) << FolyMenny;
        Node_settings.getChildNode("fluid_volume").addText(os.str().c_str());
    } else
        Node_settings.getChildNode("solution_exists").addText("false");

    // Csomopontok es agak szamanak es nevenek kiolvasasa...
    csp_db = Node_nodes.nChildNode("node");
    ag_db = Node_edges.nChildNode("edge");

    string id, is_endnode;
    bool megvan;
    double p;
    for (int i = 0; i < csp_db; i++) {
        id = Node_nodes.getChildNode("node", i).getChildNode("id").getText();
        megvan = false;
        if (debug)
            cout << endl << "\tcsp id: " << id << "  ";
        is_endnode = Node_nodes.getChildNode("node", i).getChildNode("endnode").getText();
        if (is_endnode == "false") {
            for (unsigned int j = 0; j < cspok.size(); j++) {
                if (id == (cspok.at(j)->Get_nev())) {
                    Node_nodes.getChildNode("node", i).getChildNode("pressure").deleteText();
                    p = cspok.at(j)->Get_dprop("ro") * 9.81 * cspok.at(j)->Get_p();
                    ostringstream os;
                    os << scientific << setprecision(5) << p;
                    Node_nodes.getChildNode("node", i).getChildNode("pressure").addText(os.str().c_str());
                    megvan = true;
                    /*if (debug)
                        cout << " => pressure=" << p << "Pa = " << p / 1000 / 9.81 << "m";*/

                    Node_nodes.getChildNode("node", i).getChildNode("head").deleteText();
                    p = cspok.at(j)->Get_p();
                    os.str("");
                    os << scientific << setprecision(5) << p;
                    Node_nodes.getChildNode("node", i).getChildNode("head").addText(os.str().c_str());

                    Node_nodes.getChildNode("node", i).getChildNode("travel_time").deleteText();
                    p = cspok.at(j)->Get_dprop("tt");
                    os.str("");
                    os << scientific << setprecision(5) << p / 3600.;
                    Node_nodes.getChildNode("node", i).getChildNode("travel_time").addText(os.str().c_str());

                    //cout<<endl<<cspok.at(j)->Get_nev()<<" tt="<<p;
                    //cin.get();
                }
            }
            if (debug && !megvan)
                cout << "nincs meg, atugrom...";
        } else {
            if (debug)
                cout << " ez endnode, nem erdekes...";
        }
    }

    if (debug) {
        cout << endl << endl << "Csomopontok eredmenyeinek kiirasa kesz." << endl;
    }

    double mp, q, v, dh, dhpL, tt, aref;
    for (int i = 0; i < ag_db; i++) {
        id = Node_edges.getChildNode("edge", i).getChildNode("id").getText();
        megvan = false;
        if (debug)
            cout << endl << "\tag id : " << id << "  ";
        for (unsigned int j = 0; j < agelemek.size(); j++) {
            //cout << endl << "\t agelemek.at(" << j << ")->Get_nev()=" << agelemek.at(j)->Get_nev();
            if (id == (agelemek.at(j)->Get_nev())) {
                mp = agelemek.at(i)->Get_mp();
                Node_edges.getChildNode("edge", i).getChildNode("mass_flow_rate").deleteText();
                ostringstream os;
                os << scientific << setprecision(5) << mp;
                Node_edges.getChildNode("edge", i).getChildNode("mass_flow_rate").addText(os.str().c_str());

                q = agelemek.at(i)->Get_Q();
                Node_edges.getChildNode("edge", i).getChildNode("volume_flow_rate").deleteText();
                os.str("");
                os << scientific << setprecision(5) << q * 3600;
                Node_edges.getChildNode("edge", i).getChildNode("volume_flow_rate").addText(os.str().c_str());

                v = agelemek.at(i)->Get_v();
                Node_edges.getChildNode("edge", i).getChildNode("velocity").deleteText();
                os.str("");
                os << scientific << setprecision(5) << v;
                Node_edges.getChildNode("edge", i).getChildNode("velocity").addText(os.str().c_str());

                dh = agelemek.at(i)->Get_dprop("headloss");
                Node_edges.getChildNode("edge", i).getChildNode("headloss").deleteText();
                os.str("");
                os << scientific << setprecision(5) << dh;
                Node_edges.getChildNode("edge", i).getChildNode("headloss").addText(os.str().c_str());

                dhpL = agelemek.at(i)->Get_dprop("headloss_per_unit_length");
                Node_edges.getChildNode("edge", i).getChildNode("head_loss_per_unit_length").deleteText();
                os.str("");
                os << scientific << setprecision(5) << dhpL;
                Node_edges.getChildNode("edge", i).getChildNode("head_loss_per_unit_length").addText(os.str().c_str());

                tt = (agelemek.at(i)->Get_tt_start() + agelemek.at(i)->Get_tt_end()) / 2. / 3600.;
                Node_edges.getChildNode("edge", i).getChildNode("travel_time").deleteText();
                os.str("");
                os << scientific << setprecision(5) << tt;
                Node_edges.getChildNode("edge", i).getChildNode("travel_time").addText(os.str().c_str());

                aref = agelemek.at(i)->Get_Aref();
                Node_edges.getChildNode("edge", i).getChildNode("aref").deleteText();
                os.str("");
                os << scientific << setprecision(5) << aref;
                Node_edges.getChildNode("edge", i).getChildNode("aref").addText(os.str().c_str());


                megvan = true;
                if (debug)
                    cout << " \t=> mass_flow_rate=" << mp
                         << " kg/s, volume_flow_rate=" << q
                         << " m^3/s, velocity=" << v << " m/s"
                         << " m/s, head_loss=" << dh << " m";

                // Jelleggorbes fojtas - interpolalt dzeta ertek
                if (agelemek.at(j)->Get_Tipus() == "Jelleggorbes fojtas") {
                    XMLNode akt_node;

                    if (Node_edges.getChildNode("edge", i).getChildNode("edge_spec").nChildNode("valve") > 0)
                        akt_node = Node_edges.getChildNode("edge", i).getChildNode("edge_spec").getChildNode(
                                "valve").getChildNode("adzeta");

                    akt_node.deleteText();
                    ostringstream os;
                    os << scientific << setprecision(5) << agelemek.at(j)->Get_dprop("veszt");
                    akt_node.addText(os.str().c_str());


                }

                // Cso es Csatorna eseten kiirjuk a referncia keresztmetszet sz�m�tott �rt�k�t
                if (agelemek.at(j)->Get_Tipus() == "Csatorna" || agelemek.at(j)->Get_Tipus() == "Cso") {
                    XMLNode akt_node;
                    if (Node_edges.getChildNode("edge", i).getChildNode("edge_spec").nChildNode("channel1") > 0)
                        akt_node = Node_edges.getChildNode("edge", i).getChildNode("edge_spec").getChildNode(
                                "channel1");
                    if (Node_edges.getChildNode("edge", i).getChildNode("edge_spec").nChildNode("pipe") > 0)
                        akt_node = Node_edges.getChildNode("edge", i).getChildNode("edge_spec").getChildNode("pipe");

                    Node_edges.getChildNode("edge", i).getChildNode("aref").deleteText();

                    ostringstream os;
                    os << scientific << setprecision(5) << agelemek.at(j)->Get_dprop("Aref");
                    Node_edges.getChildNode("edge", i).getChildNode("aref").addText(os.str().c_str());
                }

                // Csatorna eseten kiirjuk: node_from, node_to, lejtes
                if (agelemek.at(j)->Get_Tipus() == "Csatorna") {
                    XMLNode akt_node;
                    if (Node_edges.getChildNode("edge", i).getChildNode("edge_spec").nChildNode("channel1") > 0)
                        akt_node = Node_edges.getChildNode("edge", i).getChildNode("edge_spec").getChildNode(
                                "channel1");

                    Node_edges.getChildNode("edge", i).getChildNode("node_from").deleteText();
                    Node_edges.getChildNode("edge", i).getChildNode("node_from").addText(
                            (agelemek.at(j)->Get_Cspe_Nev()).c_str());

                    Node_edges.getChildNode("edge", i).getChildNode("node_to").deleteText();
                    Node_edges.getChildNode("edge", i).getChildNode("node_to").addText(
                            (agelemek.at(j)->Get_Cspv_Nev()).c_str());
                }

                // Speci resz: Cso es Csatorna eseten kiirjuk a csosurlodasi tenyezo erteket
                if (agelemek.at(j)->Get_Tipus() == "Csatorna" || agelemek.at(j)->Get_Tipus() == "Cso") {
                    XMLNode akt_node;
                    if (Node_edges.getChildNode("edge", i).getChildNode("edge_spec").nChildNode("channel") > 0)
                        akt_node = Node_edges.getChildNode("edge", i).getChildNode("edge_spec").getChildNode("channel");
                    if (Node_edges.getChildNode("edge", i).getChildNode("edge_spec").nChildNode("channel0") > 0)
                        akt_node = Node_edges.getChildNode("edge", i).getChildNode("edge_spec").getChildNode(
                                "channel0");
                    if (Node_edges.getChildNode("edge", i).getChildNode("edge_spec").nChildNode("channel1") > 0)
                        akt_node = Node_edges.getChildNode("edge", i).getChildNode("edge_spec").getChildNode(
                                "channel1");
                    if (Node_edges.getChildNode("edge", i).getChildNode("edge_spec").nChildNode("channel2") > 0)
                        akt_node = Node_edges.getChildNode("edge", i).getChildNode("edge_spec").getChildNode(
                                "channel2");
                    if (Node_edges.getChildNode("edge", i).getChildNode("edge_spec").nChildNode("pipe") > 0)
                        akt_node = Node_edges.getChildNode("edge", i).getChildNode("edge_spec").getChildNode("pipe");

                    akt_node.getChildNode("friction_coeff").deleteText();
                    ostringstream os;
                    os << scientific << setprecision(5) << agelemek.at(j)->Get_dprop("lambda");
                    akt_node.getChildNode("friction_coeff").addText(os.str().c_str());
                }

                if (agelemek.at(j)->Get_Tipus() == "Csatorna") {
                    XMLNode akt_node;
                    if (Node_edges.getChildNode("edge", i).getChildNode("edge_spec").nChildNode("channel") > 0)
                        akt_node = Node_edges.getChildNode("edge", i).getChildNode("edge_spec").getChildNode("channel");
                    if (Node_edges.getChildNode("edge", i).getChildNode("edge_spec").nChildNode("channel0") > 0)
                        akt_node = Node_edges.getChildNode("edge", i).getChildNode("edge_spec").getChildNode(
                                "channel0");
                    if (Node_edges.getChildNode("edge", i).getChildNode("edge_spec").nChildNode("channel1") > 0)
                        akt_node = Node_edges.getChildNode("edge", i).getChildNode("edge_spec").getChildNode(
                                "channel1");
                    if (Node_edges.getChildNode("edge", i).getChildNode("edge_spec").nChildNode("channel2") > 0)
                        akt_node = Node_edges.getChildNode("edge", i).getChildNode("edge_spec").getChildNode(
                                "channel2");

                    akt_node.getChildNode("start_height").deleteText();
                    ostringstream os;
                    os.str("");
                    os << scientific << setprecision(5) << agelemek.at(j)->Get_dprop("ze");
                    akt_node.getChildNode("start_height").addText(os.str().c_str());

                    akt_node.getChildNode("end_height").deleteText();
                    os.str("");
                    os << scientific << setprecision(5) << agelemek.at(j)->Get_dprop("zv");
                    akt_node.getChildNode("end_height").addText(os.str().c_str());

                    akt_node.getChildNode("inclination").deleteText();
                    os.str("");
                    os << scientific << setprecision(5) << (100 * (agelemek.at(j)->Get_dprop("lejtes")));
                    akt_node.getChildNode("inclination").addText(os.str().c_str());
                }


                // Speci resz: Csatorna eseten kiirjuk a sebesseg- es szintelszlast is
                if (agelemek.at(j)->Get_Tipus() == "Csatorna") {
                    XMLNode akt_node;
                    if (Node_edges.getChildNode("edge", i).getChildNode("edge_spec").nChildNode("channel") > 0)
                        akt_node = Node_edges.getChildNode("edge", i).getChildNode("edge_spec").getChildNode("channel");
                    if (Node_edges.getChildNode("edge", i).getChildNode("edge_spec").nChildNode("channel0") > 0)
                        akt_node = Node_edges.getChildNode("edge", i).getChildNode("edge_spec").getChildNode(
                                "channel0");
                    if (Node_edges.getChildNode("edge", i).getChildNode("edge_spec").nChildNode("channel1") > 0)
                        akt_node = Node_edges.getChildNode("edge", i).getChildNode("edge_spec").getChildNode(
                                "channel1");
                    if (Node_edges.getChildNode("edge", i).getChildNode("edge_spec").nChildNode("channel2") > 0)
                        akt_node = Node_edges.getChildNode("edge", i).getChildNode("edge_spec").getChildNode(
                                "channel2");


                    for (int nn = 0; nn < akt_node.nChildNode("curve"); nn++) {
                        XMLNode xnode = akt_node.getChildNode("curve", nn);
                        string melyik = xnode.getChildNode("id").getText();
                        xnode.getChildNode("points").deleteNodeContent();
                        xnode.addChild(XMLNode::parseString("<points> </points>"));
                        ostringstream os;
                        os.str("");
                        os << scientific << setprecision(5);

                        if (melyik == "curve_yf") {
                            // Fen�kg�rbe
                            vector<double> xf = agelemek.at(j)->Get_res("xf");
                            vector<double> yf = agelemek.at(j)->Get_res("yf");

                            for (unsigned int pdb = 0; pdb < xf.size(); pdb++)
                                os << "<point_x>" << xf.at(pdb)
                                   << "</point_x><point_y>" << yf.at(pdb)
                                   << "</point_y>";
                            if (debug)
                                cout << endl << " \t\t fenek kesz";
                        }

                        if (melyik == "curve_y") {
                            // V�zfelszin g�rbe
                            vector<double> x = agelemek.at(j)->Get_res("x");
                            vector<double> y = agelemek.at(j)->Get_res("y");
                            for (unsigned int pdb = 0; pdb < x.size(); pdb++)
                                os << "<point_x>" << x.at(pdb)
                                   << "</point_x><point_y>" << y.at(pdb)
                                   << "</point_y>";
                            if (debug)
                                cout << endl << " \t\t vizfelszin kesz";
                        }

                        if (melyik == "curve_p") {
                            // Nyom�sg�rbe
                            vector<double> x = agelemek.at(j)->Get_res("x");
                            vector<double> y = agelemek.at(j)->Get_res("p");
                            for (unsigned int pdb = 0; pdb < x.size(); pdb++)
                                os << "<point_x>" << x.at(pdb)
                                   << "</point_x><point_y>" << y.at(pdb)
                                   << "</point_y>";
                            if (debug)
                                cout << endl << " \t\t vizfelszin kesz";
                        }

                        if (melyik == "curve_v") {
                            // Sebess�geloszl�s
                            vector<double> x = agelemek.at(j)->Get_res("x");
                            vector<double> v = agelemek.at(j)->Get_res("v");
                            for (unsigned int pdb = 0; pdb < x.size(); pdb++)
                                os << "<point_x>" << x.at(pdb)
                                   << "</point_x><point_y>" << v.at(pdb)
                                   << "</point_y>";
                            if (debug)
                                cout << endl << " \t\t sebesseg kesz";
                        }
                        xnode.getChildNode("points").addChild(XMLNode::parseString(os.str().c_str()));
                    }
                }
            }
        }
        if (debug && !megvan)
            cout << "nincs meg, atugrom...";
    }
    xMainNode.writeToFile(xml_fnev);
    cout << "  ok.\n\n";
}

//--------------------------------------------------------------------------------
void data_io::save_mod_prop(vector<Csomopont *> cspok, vector<Agelem *> agelemek, string eID, string pID) {

    debug = false;

    XMLNode xMainNode = XMLNode::openFileHelper(xml_fnev, "staci");

    cout << endl << "Kimeneti fajl neve:" << xml_fnev << ", modositott adatok mentese... (eID:" << eID << ", pID="
         << pID << ")";

    XMLNode Node_settings = xMainNode.getChildNode("settings");
    XMLNode Node_nodes = xMainNode.getChildNode("nodes");
    XMLNode Node_edges = xMainNode.getChildNode("edges");

    if (debug)
        cout << endl << endl << "cpp_xml_debug=" << Node_settings.getChildNode("cpp_xml_debug").getText()
             << " => adatfajl irasa kozbeni debug bekapcsolasa..." << endl;

    // Csomopontok es agak szamanak es nevenek kiolvasasa...
    csp_db = Node_nodes.nChildNode("node");
    ag_db = Node_edges.nChildNode("edge");

    string id;
    bool megvan;

    for (int i = 0; i < ag_db; i++) {
        id = Node_edges.getChildNode("edge", i).getChildNode("id").getText();
        megvan = false;
        if (debug) {
            cout << endl << "\tag id : " << id << ",  eID=" << eID;

        }
        for (unsigned int j = 0; j < agelemek.size(); j++) {
            if ((id == (agelemek.at(j)->Get_nev())) && (id == eID)) {
                // cout << endl << "\tag id : " << id << ",  eID="<<eID;
                //  cout << endl << "\tpID : " << pID;
                //  cout << endl << "\tval : " << agelemek.at(j)->Get_dprop(pID);
                //  cout << endl << "\ttip : " <<agelemek.at(j)->Get_Tipus();
                // Cso
                if (agelemek.at(j)->Get_Tipus() == "Cso") {
                    XMLNode akt_node;

                    if (Node_edges.getChildNode("edge", i).getChildNode("edge_spec").nChildNode("pipe") > 0)
                        akt_node = Node_edges.getChildNode("edge", i).getChildNode("edge_spec").getChildNode("pipe");

                    akt_node.getChildNode(pID.c_str()).deleteText();
                    ostringstream os;
                    os << scientific << setprecision(5) << agelemek.at(j)->Get_dprop(pID);
                    akt_node.getChildNode(pID.c_str()).addText(os.str().c_str());
                }

                // Vegakna
                if (agelemek.at(j)->Get_Tipus() == "Vegakna") {
                    XMLNode akt_node;

                    if (Node_edges.getChildNode("edge", i).getChildNode("edge_spec").nChildNode("pool") > 0)
                        akt_node = Node_edges.getChildNode("edge", i).getChildNode("edge_spec").getChildNode("pool");

                    akt_node.getChildNode(pID.c_str()).deleteText();
                    ostringstream os;
                    os << scientific << setprecision(5) << agelemek.at(j)->Get_dprop(pID);
                    akt_node.getChildNode(pID.c_str()).addText(os.str().c_str());
                }

                // Jelleggorbes fojtas
                if (agelemek.at(j)->Get_Tipus() == "Jelleggorbes fojtas") {
                    XMLNode akt_node;

                    if (Node_edges.getChildNode("edge", i).getChildNode("edge_spec").nChildNode("valve") > 0)
                        akt_node = Node_edges.getChildNode("edge", i).getChildNode("edge_spec").getChildNode("valve");

                    akt_node.getChildNode(pID.c_str()).deleteText();
                    ostringstream os;
                    os << scientific << setprecision(5) << agelemek.at(j)->Get_dprop(pID);
                    akt_node.getChildNode(pID.c_str()).addText(os.str().c_str());

                    if (pID == "position") {
                        akt_node.getChildNode("position").deleteText();
                        ostringstream os;
                        os << scientific << setprecision(5) << agelemek.at(j)->Get_dprop("position");
                        akt_node.getChildNode("position").addText(os.str().c_str());

                        akt_node.getChildNode("adzeta").deleteText();
                        os.str("");
                        os << scientific << setprecision(5) << agelemek.at(j)->Get_dprop("adzeta");
                        akt_node.getChildNode("adzeta").addText(os.str().c_str());
                    }
                }
            }

        } //for (unsigned int j=0; j<agelemek.size(); j++)

        if (debug && !megvan)
            cout << "nincs meg, atugrom...";
    } // for (int i=0; i<ag_db; i++)


    for (int i = 0; i < csp_db; i++) {
        id = Node_nodes.getChildNode("node", i).getChildNode("id").getText();
        megvan = false;
        if (debug)
            cout << endl << "\tcsp id : " << id << "  ";
        for (unsigned int j = 0; j < cspok.size(); j++) {
            if ((id == (cspok.at(j)->Get_nev())) && (id == eID)) {
                XMLNode akt_node = Node_nodes.getChildNode("node", i);
                akt_node.getChildNode(pID.c_str()).deleteText();
                ostringstream os;
                os << scientific << setprecision(5) << cspok.at(j)->Get_dprop(pID);
                akt_node.getChildNode(pID.c_str()).addText(os.str().c_str());
            }

        } //for (unsigned int j=0; j<agelemek.size(); j++)

        if (debug && !megvan)
            cout << "nincs meg, atugrom...";
    } // for (int i=0; i<ag_db; i++)


    xMainNode.writeToFile(xml_fnev);
    //cout<<"  kesz\n\n";
}

//--------------------------------------------------------------------------------
void data_io::load_ini_values(vector<Csomopont *> &cspok, vector<Agelem *> &agelemek) {

    XMLNode xMainNode = XMLNode::openFileHelper(xml_fnev, "staci");

    XMLNode Node_settings = xMainNode.getChildNode("settings");
    XMLNode Node_nodes = xMainNode.getChildNode("nodes");
    XMLNode Node_edges = xMainNode.getChildNode("edges");

    if (debug)
        cout << endl
             << "cpp_xml_debug= true => inicializacios fajl olvasasa kozbeni debug bekapcsolasa..."
             << endl;

    // Csomopontok es agak szamanak es nevenek kiolvasasa...
    csp_db = Node_nodes.nChildNode("node");
    ag_db = Node_edges.nChildNode("edge");
    if (debug)
        cout << endl << "Ini fajl: " << xml_fnev << endl << "\tcsomopontok szama: "
             << csp_db << endl << "\tagak szama:        " << ag_db << endl;

    string id;
    double pressure;
    for (int i = 0; i < csp_db; i++) {
        id = Node_nodes.getChildNode("node", i).getChildNode("id").getText();
        for (unsigned int j = 0; j < cspok.size(); j++) {
            if (cspok.at(j)->Get_nev() == id) {
                pressure = atof(Node_nodes.getChildNode("node", i).getChildNode("pressure").getText()) / 1000 / 9.81;
                cspok.at(j)->Ini(1, pressure);
                if (debug)
                    cout << endl << "\t id: " << id << " =>  p=" << (pressure * 1000 * 9.81)
                         << "Pa = " << pressure << "vom";
            }
        }
    }

    double mass_flow_rate;
    for (int i = 0; i < ag_db; i++) {
        id = Node_edges.getChildNode("edge", i).getChildNode("id").getText();
        for (unsigned int j = 0; j < agelemek.size(); j++) {
            if (agelemek.at(j)->Get_nev() == id) {
                mass_flow_rate = atof(Node_edges.getChildNode("edge", i).getChildNode("mass_flow_rate").getText());
                if (debug)
                    cout << endl << "\t id: " << id << " => mp=" << mass_flow_rate
                         << " kg/s";
                agelemek.at(j)->Ini(1, mass_flow_rate);
            }
        }
    }
}

//--------------------------------------------------------------------------------
void data_io::save_transport(int mode, vector<Csomopont *> cspok,
                             vector<Agelem *> agelemek) {
    XMLNode xMainNode = XMLNode::openFileHelper(xml_fnev, "staci");

    XMLNode Node_settings = xMainNode.getChildNode("settings");
    XMLNode Node_nodes = xMainNode.getChildNode("nodes");
    XMLNode Node_edges = xMainNode.getChildNode("edges");

    if (debug)
        cout << endl << endl << "cpp_xml_debug=" << Node_settings.getChildNode("cpp_xml_debug").getText()
             << " => adatfajl irasa kozbeni debug bekapcsolasa..." << endl;

    // Csomopontok es agak szamanak es nevenek kiolvasasa...
    csp_db = Node_nodes.nChildNode("node");
    ag_db = Node_edges.nChildNode("edge");

    string tagname;
    if (mode == 1) {
        tagname = "travel_time";
        Node_settings.getChildNode("tsolution_exists").deleteText();
        Node_settings.getChildNode("tsolution_exists").addText("true");
    }
    if (mode == 2) {
        tagname = "concentration";
        Node_settings.getChildNode("csolution_exists").deleteText();
        Node_settings.getChildNode("csolution_exists").addText("true");

    }

    string id, is_endnode;
    bool megvan;
    for (int i = 0; i < csp_db; i++) {
        id = Node_nodes.getChildNode("node", i).getChildNode("id").getText();
        megvan = false;
        if (debug)
            cout << endl << "\tcsp id: " << id << "  ";
        is_endnode = Node_nodes.getChildNode("node", i).getChildNode("endnode").getText();
        if (is_endnode == "false") {
            for (unsigned int j = 0; j < cspok.size(); j++) {
                if (id == (cspok.at(j)->Get_nev())) {
                    Node_nodes.getChildNode("node", i).getChildNode(tagname.c_str()).deleteText();
                    double c = cspok.at(j)->Get_dprop("konc_atlag");
                    if (mode == 1)
                        c = c / 3600;

                    ostringstream os;
                    os << scientific << setprecision(5) << c;
                    Node_nodes.getChildNode("node", i).getChildNode(tagname.c_str()).addText(os.str().c_str());
                    megvan = true;
                    if (debug)
                        cout << " => " << tagname << ": " << c;
                }
            }
            if (debug && !megvan)
                cout << "nincs meg, atugrom...";
        } else {
            if (debug)
                cout << " ez endnode, nem erdekes...";
        }
    }

    double mp;
    for (int i = 0; i < ag_db; i++) {
        id = Node_edges.getChildNode("edge", i).getChildNode("id").getText();
        megvan = false;
        if (debug)
            cout << endl << "\tag id : " << id << "  ";
        for (unsigned int j = 0; j < agelemek.size(); j++) {
            if (id == (agelemek.at(j)->Get_nev())) {
                mp = agelemek.at(i)->Get_mp();
                Node_edges.getChildNode("edge", i).getChildNode(tagname.c_str()).deleteText();
                ostringstream os;
                double konc = agelemek.at(j)->mean(agelemek.at(j)->konc);
                if (mode == 1)
                    konc = konc / 3600;
                os << scientific << setprecision(5) << konc;
                Node_edges.getChildNode("edge", i).getChildNode(tagname.c_str()).addText(os.str().c_str());

                megvan = true;
            }
        }
        if (debug && !megvan)
            cout << "nincs meg, atugrom...";
    }
    xMainNode.writeToFile(xml_fnev);
}

//--------------------------------------------------------------------------------
string data_io::read_setting(string which) {
    XMLNode xMainNode = XMLNode::openFileHelper(xml_fnev, "staci");

    XMLNode Node_settings = xMainNode.getChildNode("settings");

    string out;
    if (Node_settings.nChildNode(which.c_str()) == 1) {
        /*XMLNode thisnode=xMainNode.getChildNode("settings").getChildNode(which.c_str());
        if (!thisnode.isEmpty())
            out=thisnode.getText();
        else
            out="";*/
        out = xMainNode.getChildNode("settings").getChildNode(which.c_str()).getText();
        //if (debug) cout<<endl<<"\t"<<which<<" => "<<out;
        return out;
    } else
        return "nincs ilyen node!";
}

//--------------------------------------------------------------------------------
double
data_io::string_to_double(const string &s, const string &elem_name, const string &tag_name, const double &def_value) {
    std::istringstream i(s);
    double x;
    if (s.empty()) {
        cout << endl << endl << "Error! Element: " << elem_name << ", tag: " << tag_name << " is empty." << endl;
        cout << endl << "  Returning default value of " << def_value << " but you should check the input file." << endl
             << endl;
        cin.get();
        return def_value;
    } else {
        if (!(i >> x)) {
            cout << endl << endl << "Error! Element: " << elem_name << ", tag: " << tag_name << " - non-numeric value: "
                 << s << endl;
            cout << endl << "  Returning default value of " << def_value << " but you should check the input file."
                 << endl << endl;
            cin.get();
            return def_value;
        } else {
            return x;
        }
    }
}
