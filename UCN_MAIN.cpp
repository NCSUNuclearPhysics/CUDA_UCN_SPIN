/* Include files ----------------------------------------------------------------- */


#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <string.h>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>
#include <map>
#include <set>
#include <stdio.h>
#include <stdlib.h>
#include <time.h> 

using namespace std;

#include "UCN_CUDA_WRAPPER.h"
#include "UCN_RUN.h"
#include "UCN_MAGNETIC.h"

int main()
{
	map<int, double> param_m_vi_CONST;
	map<int, int> param_m_vi_CONST_INT;
	param_m_vi_CONST[e_d_CONST_def_MAXERR] = 1e-10;
	param_m_vi_CONST[e_d_CONST_def_MAXERR1] = 1e0;
	c_UCN_RUN_DATA obj_MOAR_c_UCN_RUN_DATA_TEST_0 = fxn_PARAMETERIZED_CONSTRUCTOR_c_UCN_RUN_DATA(
		/* const int vi_INPUT_e_d_CONST_INT_numThreadsPerBlock_EXP_MAX, */ 5, 
		/* const int vi_INPUT_e_d_CONST_INT_numRecordsPerThread_EXP, */ 12, 
		/* const int vi_INPUT_e_d_CONST_INT_numThreads_EXP, */ 2,
		/* const int vi_INPUT_e_d_CONST_INT_numRecordsPerUpdate, */ 256,
		/* const int vi_INPUT_e_d_CONST_INT_numCyclesPerRecord, */ 1,
		param_m_vi_CONST_INT, 
		param_m_vi_CONST);
	obj_MOAR_c_UCN_RUN_DATA_TEST_0.fxn_FORMATTED_WRITE_MINIMUM("UCN", true);
	return 0;
	param_m_vi_CONST.erase(e_d_CONST_def_MAXERR);
	param_m_vi_CONST.erase(e_d_CONST_def_MAXERR1);
	param_m_vi_CONST[e_d_CONST_def_MAXERR] = 1e-8;
	param_m_vi_CONST[e_d_CONST_def_MAXERR1] = 1e-9;
	c_UCN_RUN_DATA obj_MOAR_c_UCN_RUN_DATA_TEST_1 = fxn_PARAMETERIZED_CONSTRUCTOR_c_UCN_RUN_DATA(
		/* const int vi_INPUT_e_d_CONST_INT_numThreadsPerBlock_EXP_MAX, */ 5, 
		/* const int vi_INPUT_e_d_CONST_INT_numRecordsPerThread_EXP, */ 16, 
		/* const int vi_INPUT_e_d_CONST_INT_numThreads_EXP, */ 2,
		/* const int vi_INPUT_e_d_CONST_INT_numRecordsPerUpdate, */ 1,
		/* const int vi_INPUT_e_d_CONST_INT_numCyclesPerRecord, */ 100,
		param_m_vi_CONST_INT, 
		param_m_vi_CONST);
	obj_MOAR_c_UCN_RUN_DATA_TEST_1.fxn_FORMATTED_WRITE("UCN", true);

	param_m_vi_CONST.erase(e_d_CONST_def_MAXERR);
	param_m_vi_CONST.erase(e_d_CONST_def_MAXERR1);
	param_m_vi_CONST[e_d_CONST_def_MAXERR] = 1e-7;
	param_m_vi_CONST[e_d_CONST_def_MAXERR1] = 1e-9;
	c_UCN_RUN_DATA obj_MOAR_c_UCN_RUN_DATA_TEST_2 = fxn_PARAMETERIZED_CONSTRUCTOR_c_UCN_RUN_DATA(
		/* const int vi_INPUT_e_d_CONST_INT_numThreadsPerBlock_EXP_MAX, */ 5, 
		/* const int vi_INPUT_e_d_CONST_INT_numRecordsPerThread_EXP, */ 16, 
		/* const int vi_INPUT_e_d_CONST_INT_numThreads_EXP, */ 2,
		/* const int vi_INPUT_e_d_CONST_INT_numRecordsPerUpdate, */ 1,
		/* const int vi_INPUT_e_d_CONST_INT_numCyclesPerRecord, */ 1000,
		param_m_vi_CONST_INT, 
		param_m_vi_CONST);
	obj_MOAR_c_UCN_RUN_DATA_TEST_2.fxn_FORMATTED_WRITE("UCN", true);

	param_m_vi_CONST.erase(e_d_CONST_def_MAXERR);
	param_m_vi_CONST.erase(e_d_CONST_def_MAXERR1);
	param_m_vi_CONST[e_d_CONST_def_MAXERR] = 1e-8;
	param_m_vi_CONST[e_d_CONST_def_MAXERR1] = 1e-10;
	c_UCN_RUN_DATA obj_MOAR_c_UCN_RUN_DATA_TEST_3 = fxn_PARAMETERIZED_CONSTRUCTOR_c_UCN_RUN_DATA(
		/* const int vi_INPUT_e_d_CONST_INT_numThreadsPerBlock_EXP_MAX, */ 5, 
		/* const int vi_INPUT_e_d_CONST_INT_numRecordsPerThread_EXP, */ 16, 
		/* const int vi_INPUT_e_d_CONST_INT_numThreads_EXP, */ 2,
		/* const int vi_INPUT_e_d_CONST_INT_numRecordsPerUpdate, */ 1,
		/* const int vi_INPUT_e_d_CONST_INT_numCyclesPerRecord, */ 1000,
		param_m_vi_CONST_INT, 
		param_m_vi_CONST);
	obj_MOAR_c_UCN_RUN_DATA_TEST_3.fxn_FORMATTED_WRITE("UCN", true);
	param_m_vi_CONST.erase(e_d_CONST_def_MAXERR);
	param_m_vi_CONST.erase(e_d_CONST_def_MAXERR1);
	return 0;

	// param_m_vi_CONST[e_d_CONST_def_MAXERR] = 1e-9;
	// param_m_vi_CONST[e_d_CONST_def_MAXERR1] = 1e-11;
	// c_UCN_RUN_DATA obj_MOAR_c_UCN_RUN_DATA_TEST = fxn_PARAMETERIZED_CONSTRUCTOR_c_UCN_RUN_DATA(
		// /* const int vi_INPUT_e_d_CONST_INT_numThreadsPerBlock_EXP_MAX, */ 5, 
		// /* const int vi_INPUT_e_d_CONST_INT_numRecordsPerThread_EXP, */ 6, 
		// /* const int vi_INPUT_e_d_CONST_INT_numThreads_EXP, */ 2,
		// /* const int vi_INPUT_e_d_CONST_INT_numRecordsPerUpdate, */ 1,
		// /* const int vi_INPUT_e_d_CONST_INT_numCyclesPerRecord, */ 100,
		// param_m_vi_CONST_INT, 
		// param_m_vi_CONST);
	// obj_MOAR_c_UCN_RUN_DATA_TEST.fxn_FORMATTED_WRITE("UCN", true);
	// return 0;
	// param_m_vi_CONST.erase(e_d_CONST_def_MAXERR);
	// param_m_vi_CONST.erase(e_d_CONST_def_MAXERR1);
	param_m_vi_CONST[e_d_CONST_def_MAXERR] = 1e-6;
	param_m_vi_CONST[e_d_CONST_def_MAXERR1] = 1e-8;
	c_UCN_RUN_DATA obj_MOAR_c_UCN_RUN_DATA_0 = fxn_PARAMETERIZED_CONSTRUCTOR_c_UCN_RUN_DATA(
		/* const int vi_INPUT_e_d_CONST_INT_numThreadsPerBlock_EXP_MAX, */ 5, 
		/* const int vi_INPUT_e_d_CONST_INT_numRecordsPerThread_EXP, */ 12, 
		/* const int vi_INPUT_e_d_CONST_INT_numThreads_EXP, */ 2,
		/* const int vi_INPUT_e_d_CONST_INT_numRecordsPerUpdate, */ 16,
		/* const int vi_INPUT_e_d_CONST_INT_numCyclesPerRecord, */ 100,
		param_m_vi_CONST_INT, 
		param_m_vi_CONST);
	obj_MOAR_c_UCN_RUN_DATA_0.fxn_FORMATTED_WRITE("UCN_N6", true);
	param_m_vi_CONST.erase(e_d_CONST_def_MAXERR);
	param_m_vi_CONST[e_d_CONST_def_MAXERR] = 1e-7;
	param_m_vi_CONST.erase(e_d_CONST_def_MAXERR1);
	param_m_vi_CONST[e_d_CONST_def_MAXERR1] = 1e-9;
	c_UCN_RUN_DATA obj_MOAR_c_UCN_RUN_DATA_1 = fxn_PARAMETERIZED_CONSTRUCTOR_c_UCN_RUN_DATA(
		/* const int vi_INPUT_e_d_CONST_INT_numThreadsPerBlock_EXP_MAX, */ 5, 
		/* const int vi_INPUT_e_d_CONST_INT_numRecordsPerThread_EXP, */ 12, 
		/* const int vi_INPUT_e_d_CONST_INT_numThreads_EXP, */ 2,
		/* const int vi_INPUT_e_d_CONST_INT_numRecordsPerUpdate, */ 16,
		/* const int vi_INPUT_e_d_CONST_INT_numCyclesPerRecord, */ 100,
		param_m_vi_CONST_INT, 
		param_m_vi_CONST);
	obj_MOAR_c_UCN_RUN_DATA_1.fxn_FORMATTED_WRITE("UCN_N7", true);
	param_m_vi_CONST.erase(e_d_CONST_def_MAXERR);
	param_m_vi_CONST[e_d_CONST_def_MAXERR] = 1e-8;
	param_m_vi_CONST.erase(e_d_CONST_def_MAXERR1);
	param_m_vi_CONST[e_d_CONST_def_MAXERR1] = 1e-10;
	c_UCN_RUN_DATA obj_MOAR_c_UCN_RUN_DATA_2 = fxn_PARAMETERIZED_CONSTRUCTOR_c_UCN_RUN_DATA(
		/* const int vi_INPUT_e_d_CONST_INT_numThreadsPerBlock_EXP_MAX, */ 5, 
		/* const int vi_INPUT_e_d_CONST_INT_numRecordsPerThread_EXP, */ 12, 
		/* const int vi_INPUT_e_d_CONST_INT_numThreads_EXP, */ 2,
		/* const int vi_INPUT_e_d_CONST_INT_numRecordsPerUpdate, */ 16,
		/* const int vi_INPUT_e_d_CONST_INT_numCyclesPerRecord, */ 100,
		param_m_vi_CONST_INT, 
		param_m_vi_CONST);
	obj_MOAR_c_UCN_RUN_DATA_2.fxn_FORMATTED_WRITE("UCN_N8", true);
	param_m_vi_CONST.erase(e_d_CONST_def_MAXERR);
	param_m_vi_CONST[e_d_CONST_def_MAXERR] = 1e-9;
	param_m_vi_CONST.erase(e_d_CONST_def_MAXERR1);
	param_m_vi_CONST[e_d_CONST_def_MAXERR1] = 1e-11;
	c_UCN_RUN_DATA obj_MOAR_c_UCN_RUN_DATA_3 = fxn_PARAMETERIZED_CONSTRUCTOR_c_UCN_RUN_DATA(
		/* const int vi_INPUT_e_d_CONST_INT_numThreadsPerBlock_EXP_MAX, */ 5, 
		/* const int vi_INPUT_e_d_CONST_INT_numRecordsPerThread_EXP, */ 12, 
		/* const int vi_INPUT_e_d_CONST_INT_numThreads_EXP, */ 2,
		/* const int vi_INPUT_e_d_CONST_INT_numRecordsPerUpdate, */ 16,
		/* const int vi_INPUT_e_d_CONST_INT_numCyclesPerRecord, */ 100,
		param_m_vi_CONST_INT, 
		param_m_vi_CONST);
	obj_MOAR_c_UCN_RUN_DATA_3.fxn_FORMATTED_WRITE("UCN_N9", true);
	param_m_vi_CONST.erase(e_d_CONST_def_MAXERR);
	param_m_vi_CONST[e_d_CONST_def_MAXERR] = 1e-10;
	param_m_vi_CONST.erase(e_d_CONST_def_MAXERR1);
	param_m_vi_CONST[e_d_CONST_def_MAXERR1] = 1e-12;
	c_UCN_RUN_DATA obj_MOAR_c_UCN_RUN_DATA_4 = fxn_PARAMETERIZED_CONSTRUCTOR_c_UCN_RUN_DATA(
		/* const int vi_INPUT_e_d_CONST_INT_numThreadsPerBlock_EXP_MAX, */ 5, 
		/* const int vi_INPUT_e_d_CONST_INT_numRecordsPerThread_EXP, */ 12, 
		/* const int vi_INPUT_e_d_CONST_INT_numThreads_EXP, */ 2,
		/* const int vi_INPUT_e_d_CONST_INT_numRecordsPerUpdate, */ 16,
		/* const int vi_INPUT_e_d_CONST_INT_numCyclesPerRecord, */ 100,
		param_m_vi_CONST_INT, 
		param_m_vi_CONST);
	obj_MOAR_c_UCN_RUN_DATA_4.fxn_FORMATTED_WRITE("UCN_N10", true);
	param_m_vi_CONST.erase(e_d_CONST_def_MAXERR);
	param_m_vi_CONST[e_d_CONST_def_MAXERR] = 1e-11;
	param_m_vi_CONST.erase(e_d_CONST_def_MAXERR1);
	param_m_vi_CONST[e_d_CONST_def_MAXERR1] = 1e-13;
	c_UCN_RUN_DATA obj_MOAR_c_UCN_RUN_DATA_5 = fxn_PARAMETERIZED_CONSTRUCTOR_c_UCN_RUN_DATA(
		/* const int vi_INPUT_e_d_CONST_INT_numThreadsPerBlock_EXP_MAX, */ 5, 
		/* const int vi_INPUT_e_d_CONST_INT_numRecordsPerThread_EXP, */ 12, 
		/* const int vi_INPUT_e_d_CONST_INT_numThreads_EXP, */ 2,
		/* const int vi_INPUT_e_d_CONST_INT_numRecordsPerUpdate, */ 16,
		/* const int vi_INPUT_e_d_CONST_INT_numCyclesPerRecord, */ 100,
		param_m_vi_CONST_INT, 
		param_m_vi_CONST);
	obj_MOAR_c_UCN_RUN_DATA_5.fxn_FORMATTED_WRITE("UCN_N11", true);
	param_m_vi_CONST.erase(e_d_CONST_def_MAXERR);
	param_m_vi_CONST[e_d_CONST_def_MAXERR] = 1e-12;
	param_m_vi_CONST.erase(e_d_CONST_def_MAXERR1);
	param_m_vi_CONST[e_d_CONST_def_MAXERR1] = 1e-14;
	c_UCN_RUN_DATA obj_MOAR_c_UCN_RUN_DATA_6 = fxn_PARAMETERIZED_CONSTRUCTOR_c_UCN_RUN_DATA(
		/* const int vi_INPUT_e_d_CONST_INT_numThreadsPerBlock_EXP_MAX, */ 5, 
		/* const int vi_INPUT_e_d_CONST_INT_numRecordsPerThread_EXP, */ 12, 
		/* const int vi_INPUT_e_d_CONST_INT_numThreads_EXP, */ 2,
		/* const int vi_INPUT_e_d_CONST_INT_numRecordsPerUpdate, */ 16,
		/* const int vi_INPUT_e_d_CONST_INT_numCyclesPerRecord, */ 100,
		param_m_vi_CONST_INT, 
		param_m_vi_CONST);
	obj_MOAR_c_UCN_RUN_DATA_6.fxn_FORMATTED_WRITE("UCN_N12", true);
	

	// int vi_FINAL_RECORD = obj_MOAR_c_UCN_RUN_DATA_0.fvci_h_CONST(e_d_CONST_INT_numRecordsPerThread) - 1;
	// for(int i = 0; i<obj_MOAR_c_UCN_RUN_DATA_0.fvci_h_CONST(e_d_CONST_INT_numThreads); i++)
	// {
		// double vd_START_T = obj_MOAR_c_UCN_RUN_DATA_0.fxn_vd_IO(i, 0, e_d_IO_T);
		// double vd_END_T = obj_MOAR_c_UCN_RUN_DATA_0.fxn_vd_IO(i, vi_FINAL_RECORD, e_d_IO_T);
		// double vd_START_SPINNOR_1 = obj_MOAR_c_UCN_RUN_DATA_0.fxn_vd_IO(i, 0, e_d_IO_SPINNOR_1);
		// double vd_END_SPINNOR_1 = obj_MOAR_c_UCN_RUN_DATA_0.fxn_vd_IO(i, vi_FINAL_RECORD, e_d_IO_SPINNOR_1);
		// double vd_START_SPINNOR_2 = obj_MOAR_c_UCN_RUN_DATA_0.fxn_vd_IO(i, 0, e_d_IO_SPINNOR_2);
		// double vd_END_SPINNOR_2 = obj_MOAR_c_UCN_RUN_DATA_0.fxn_vd_IO(i, vi_FINAL_RECORD, e_d_IO_SPINNOR_2);
		
		// double vd_DELTA_T = vd_END_T - vd_START_T;
		// double vd_START_PHI = atan(vd_START_SPINNOR_2/vd_START_SPINNOR_1);
		// double vd_END_PHI = atan(vd_END_SPINNOR_2/vd_END_SPINNOR_1);
		// double vd_DELTA_PHI = vd_END_PHI - vd_START_PHI;
		// double vd_DELTA_PHI_DELTA_T = vd_DELTA_PHI/vd_DELTA_T;
		
		// cout << "\nNEUTRON NUMBER " << i << " has a ds/dt of " << vd_DELTA_PHI_DELTA_T << " radians per second";
		// cout << "\n\tEqual to (" << vd_END_PHI << "-" << vd_START_PHI << ")/(" << vd_END_T << "-" << vd_START_T << ") or " << vd_DELTA_PHI << "/" <<  vd_DELTA_T;
		// cout << "\n\tGotten via a starting and ending (X,Y) of (" << vd_START_SPINNOR_1 << "," << vd_START_SPINNOR_2 << ") and (" << vd_END_SPINNOR_1 << "," << vd_END_SPINNOR_2 << ")";
	// }
	// param_m_vi_CONST[e_d_CONST_def_MAXERR1] = 1e-11;
	// c_UCN_RUN_DATA obj_MOAR_c_UCN_RUN_DATA_1 = fxn_PARAMETERIZED_CONSTRUCTOR_c_UCN_RUN_DATA(
	// 5, 12, 0, 256, 1, param_m_vi_CONST_INT, param_m_vi_CONST);
	// param_m_vi_CONST.clear();
	// param_m_vi_CONST[e_d_CONST_def_MAXERR1] = 1e-12;
	// c_UCN_RUN_DATA obj_MOAR_c_UCN_RUN_DATA_2 = fxn_PARAMETERIZED_CONSTRUCTOR_c_UCN_RUN_DATA(
	// 5, 12, 0, 256, 1, param_m_vi_CONST_INT, param_m_vi_CONST);
	// param_m_vi_CONST.clear();
	// param_m_vi_CONST.clear();
	// param_m_vi_CONST[e_d_CONST_def_MAXERR1] = 1e-13;
	// c_UCN_RUN_DATA obj_MOAR_c_UCN_RUN_DATA_3 = fxn_PARAMETERIZED_CONSTRUCTOR_c_UCN_RUN_DATA(
	// 5, 12, 0, 256, 1, param_m_vi_CONST_INT, param_m_vi_CONST);
	// param_m_vi_CONST.clear();
	// param_m_vi_CONST.clear();
	// param_m_vi_CONST[e_d_CONST_def_MAXERR1] = 1e-14;
	// c_UCN_RUN_DATA obj_MOAR_c_UCN_RUN_DATA_4 = fxn_PARAMETERIZED_CONSTRUCTOR_c_UCN_RUN_DATA(
	// 5, 12, 0, 256, 1, param_m_vi_CONST_INT, param_m_vi_CONST);
	// param_m_vi_CONST.clear();
	return 0;
}

