using namespace std;
#include <ctime>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include "Staci.h"
#include "time.h"
#include <numeric>

void solve_hydraulics(Staci&  feladat, double& time2, double& time3);

int main(int argc, char* argv[]) {

	cout << endl<< "STACI v2.0";
	cout << endl<< " (c) BME Dept. of Hydrodynamic Systems";
	cout << endl<< " (c) C. Hos (csaba.hos@hds.bme.hu)";
	cout << endl<< " info: www.hds.bme.hu\\staci" << endl;

	clock_t startTime, endTime;
    vector<double> timing (4,0.0);  /*!< Vector containing timing data */

	startTime = clock();
	Staci feladat(argc, argv);
	endTime = clock();
	timing.at(0)=double( endTime - startTime ) / CLOCKS_PER_SEC;

	time_t start;
	time(&start);

	// cout<<endl<<endl<<"feladat.get_mode()="<<feladat.get_mode()<<endl;
	// cin.get();    

	// Nincs argumentum...
	//----------------------------------------
	if (feladat.get_mode()!=-1)
		startTime = clock();				
	feladat.build_system();
	endTime = clock();
	timing.at(1)=double( endTime - startTime ) / CLOCKS_PER_SEC;


	// Stacioner halozatszamitas
	//----------------------------------------
	if (feladat.get_mode()==0) {
		solve_hydraulics(feladat,timing.at(2),timing.at(3));
		//feladat.solve_residence_time();
	}

	// Tartozkodasi ido szamitasa
	//----------------------------------------
	if (feladat.get_mode()==1) {
		feladat.set_out_file(feladat.get_def_file()+".rot");
		feladat.set_progress_file(feladat.get_def_file()+".rpt");
		feladat.ini();
		feladat.set_up_transport();
		feladat.solve_transport(1);
		feladat.save_transport(1);

		// Sikeres volt a futas? (itt biztos igen...)
		ofstream outfile((feladat.get_def_file()+".rrt").c_str(), ios::trunc);
		outfile<<"ok\n";
		outfile.close();

	}

	// Koncentracioeloszlas szamitasa
	//----------------------------------------
	if (feladat.get_mode()==2) {
		feladat.set_out_file(feladat.get_def_file()+".roc");
		feladat.set_progress_file(feladat.get_def_file()+".rpc");
		feladat.ini();
		feladat.set_up_transport();
		feladat.solve_transport(2);
		feladat.save_transport(2);

		// Sikeres volt a futas? (itt biztos igen...)
		ofstream outfile((feladat.get_def_file()+".rrc").c_str(), ios::trunc);
		outfile<<"ok\n";
		outfile.close();
	}

	// Valamilyen adat megvaltoztatasa
	//----------------------------------------
	if (feladat.get_mode()==3) {
		// az adatok a parancssorbol jonnek
		// az elem neve: feladat.eID
		// az adat neve: feladat.pID
		// uj adat:      feladat.newValue
		// uj adatfajl:  feladat.new_def_file

		double old_data = feladat.m_get_dprop();

		feladat.m_set_dprop();
		double new_data = feladat.m_get_dprop();

		cout << endl<< endl << "elem                : " << feladat.element_ID;
		cout << endl << "property            : "<< feladat.property_ID;
		cout << endl << "jelenlegi ertek     : " << old_data;
		cout << endl << "uj ertek            : " << new_data << endl;
		cout << endl << "modositott adatfajl : " << feladat.new_def_file<<endl;

		feladat.copy_file(feladat.get_def_file(),feladat.new_def_file);
		feladat.set_res_file(feladat.new_def_file);
		feladat.save_mod_prop();
	}

	// Minden elem listazasa
	//----------------------------------------
	if (feladat.get_mode()==4) {
		feladat.list_all_elements();
	}

	// Valamilyen adat kiolvasasa
	//----------------------------------------
	if (feladat.get_mode()==5) {
		// az adatok a parancssorbol jonnek
		// az elem neve: feladat.eID
		// az adat neve: feladat.pID		

		//double data = feladat.m_get_dprop();
		cout << endl << "element_ID  : " << feladat.element_ID;
		cout << endl << "property_ID : " << feladat.property_ID;
		cout << endl << "value       : " << feladat.m_get_dprop() << endl<<endl;		
		
	}

	// Erzekenysegvizsgalat
	//----------------------------------------
	if (feladat.get_mode()==6) {
		// az adatok a parancssorbol jonnek
		// az elem neve: feladat.eID
		// az adat neve: feladat.pID		
		cout << endl << "element_ID  : " << feladat.element_ID;
		cout << endl << "property_ID : " << feladat.property_ID;

		feladat.set_do_save_file(false);
		solve_hydraulics(feladat,timing.at(2),timing.at(3));

		feladat.Compute_dxdmu();
		feladat.Print_dxdmu();

	}

	time_t stop;
	time(&stop);

	double sum_time =timing.at(0)+timing.at(1)+timing.at(2)+timing.at(3);

	printf("\nTiming:\n============================");
	printf("\n\t Reading data files   : %5.3f sec  (%4.1f %%)",timing.at(0),timing.at(0)/sum_time*100.);
	printf("\n\t Building the system  : %5.3f sec  (%4.1f %%)",timing.at(1),timing.at(1)/sum_time*100.);
	printf("\n\t Number cranching     : %5.3f sec  (%4.1f %%)",timing.at(2),timing.at(2)/sum_time*100.);
	printf("\n\t Writing result files : %5.3f sec  (%4.1f %%)",timing.at(3),timing.at(3)/sum_time*100.);
	printf("\n\t---------------------------------------------");
	printf("\n\t Overall              : %5.3f sec\n\n",sum_time);

	return 0;
}

void solve_hydraulics(Staci&  feladat, double& time2, double& time3) {

	clock_t startTime, endTime;

	feladat.ini();
//	feladat.list_system();

	if (feladat.Get_debug_level()==-1)
	{
		cout<<"\n============================================================";
		cout<<"\ndebug_level = -1   -->";
		cout<<"\n\nExporting node names to <nodelist.txt>...";
		cout<<"\nExporting connected nodes to <connected_nodes.txt>...";
		feladat.export_connected_nodes();
		cout<<" done.\n\nUse script connected_nodes.py to find isolated islands.";
		cout<<"\n============================================================\n\n";
	}

	startTime = clock();
	bool konv_ok = feladat.solve_system();
	endTime = clock();
	time2=double( endTime - startTime ) / CLOCKS_PER_SEC;

	ostringstream msg;
	string rrs_file = feladat.get_def_file()+".rrs";
	msg << "\n\nTrying to delete previous *.rrs file ... ";
	if ( remove( rrs_file.c_str()) != 0 )
		msg << " file not found, cannot delete it.";
	else
		msg << "file successfully deleted." ;
	cout << msg.str();

	if (!konv_ok) {
		// Put "ERROR" to <problem.spr>.rrs
		ofstream outfile(rrs_file.c_str(), ios::trunc);
		outfile<<"ERROR!\n";
		outfile.close();
		feladat.logfile_write("\n\nERROR!", 0);

		cout<<endl<<endl <<"***********************************************************************************";
		cout<<endl <<      "* !!! THE SOLVER DID NOT CONVERGE, LAST ITERATION RESULTS ARE IN PROBLEM FILE !!! *";
		cout<<endl <<      "***********************************************************************************" <<endl<<endl;

	} else {
		// Put "OK" to <problem.spr>.rrs
		ofstream outfile(rrs_file.c_str(), ios::trunc);
		outfile<<"OK\n";
		outfile.close();
		feladat.logfile_write("\n\nOK", 0);

		feladat.solve_residence_time();
	}
	startTime = clock();
	feladat.set_res_file(feladat.get_def_file());
	feladat.save_results(konv_ok);
	endTime = clock();
	time3=double( endTime - startTime ) / CLOCKS_PER_SEC;

}