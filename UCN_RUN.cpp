#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <ctime>
#include <iomanip>
#include <cmath>
#include <limits>
#include <string>
#include <vector>
#include <map>
#include <time.h>
using namespace std;

#include "UCN_CUDA_WRAPPER.h"
#include "UCN_RUN.h"
// #include "UCN_STRING.h"

int ipow(int base, int exp)
{
	int result = -1;
	if(exp<0) result = 0;
	else if(exp==0) result = 1;
	else if(exp>0)
	{
		result = 1;
		for(int i = exp; i>0; i--)
		{
				result *= base;
		}
	}
  return result;
}
// c_UCN_RUN_DATA fxn_PROXY_CONSTRUCTOR_c_UCN_RUN_DATA(
	// map<int, string> & r_m_vi_cstr_IO,
	// map<int, string> & r_m_vi_cstr_IO_INT,
	// map<int, string> & r_m_vi_cstr_CONST,
	// map<int, string> & r_m_vi_cstr_CONST_INT,
	// map<int, string> & r_m_vi_cstr_RKQS_ERROR,
	// map<int, double> & r_m_vi_CONST,
	// map<int, int> & r_m_vi_CONST_INT)
// {
	// const int numThreads_TEMP = r_m_vi_CONST_INT[e_d_CONST_INT_numThreads];
	// const int numRecordsPerThread_TEMP = r_m_vi_CONST_INT[e_d_CONST_INT_numRecordsPerThread];
	// return c_UCN_RUN_DATA(
		// r_m_vi_cstr_IO,
		// r_m_vi_cstr_IO_INT,
		// r_m_vi_cstr_CONST,
		// r_m_vi_cstr_CONST_INT,
		// r_m_vi_cstr_RKQS_ERROR,
		// r_m_vi_CONST_INT,
		// r_m_vi_CONST,
		// numThreads_TEMP,
		// numRecordsPerThread_TEMP);
// }

c_UCN_RUN_DATA::c_UCN_RUN_DATA(
	map<int, string> & r_m_vi_cstr_IO,
	map<int, string> & r_m_vi_cstr_IO_INT,
	map<int, string> & r_m_vi_cstr_CONST,
	map<int, string> & r_m_vi_cstr_CONST_INT,
	map<int, string> & r_m_vi_cstr_RKQS_ERROR,
	map<int, int> & r_m_vi_CONST_INT,
	map<int, double> & r_m_vi_CONST,
	const int param_numThreads,
	const int param_numRecordsPerThread) :
	m_vi_cstr_IO(r_m_vi_cstr_IO),
	m_vi_cstr_IO_INT(r_m_vi_cstr_IO_INT),
	m_vi_cstr_CONST(r_m_vi_cstr_CONST),
	m_vi_cstr_CONST_INT(r_m_vi_cstr_CONST_INT),
	m_vi_cstr_RKQS_ERROR(r_m_vi_cstr_RKQS_ERROR),
	m_vi_CONST_INT(r_m_vi_CONST_INT),
	m_vi_CONST(r_m_vi_CONST),
	vci_numThreads(param_numThreads),
	vci_numRecordsPerThread(param_numRecordsPerThread),
	vci_numInts(param_numThreads*param_numRecordsPerThread*e_d_IO_INT_LAST),
	vci_numDoubles(param_numThreads*param_numRecordsPerThread*e_d_IO_LAST)
{
	
	// map<int, string> m_vi_cstr_IO = r_m_vi_cstr_IO;
	// map<int, string> m_vi_cstr_IO_INT = r_m_vi_cstr_IO_INT;
	// map<int, string> m_vi_cstr_CONST = r_m_vi_cstr_CONST;
	// map<int, string> m_vi_cstr_CONST_INT = r_m_vi_cstr_CONST_INT;
	// map<int, double> m_vi_CONST = r_m_vi_CONST;
	// map<int, int> m_vi_CONST_INT = r_m_vi_CONST_INT;
	
	time_t rawtime;
	struct tm * s_tm_timeinfo_INPUT;
	struct tm * s_tm_timeinfo_START;
	struct tm * s_tm_timeinfo_FINISH;
	
	// double param_CONST[e_d_CONST_LAST];
	// for(int i = 0; i<e_d_CONST_LAST; i++) param_CONST[i] = m_vi_CONST_INT[i];
	
	// int param_CONST_INT[e_d_CONST_INT_LAST];
	// for(int i = 0; i<e_d_CONST_LAST; i++) param_CONST[i] = m_vi_CONST_INT[i];
	
  // const int vci_numDoubles = m_vi_CONST_INT[e_d_CONST_INT_numRecords]*e_d_IO_LAST;
	// const int vci_numInts = m_vi_CONST_INT[e_d_CONST_INT_numRecords]*e_d_IO_INT_LAST;
	
	const int vci_numDoubleBytes = vci_numDoubles*sizeof(double);
	const int vci_numIntBytes = vci_numInts*sizeof(int);
	
	double vd_IO[vci_numDoubles];
	int vi_IO[vci_numInts];
	
	double h_CONST[e_d_CONST_LAST];
	int h_CONST_INT[e_d_CONST_INT_LAST];
	
	for(int k = 0; k<e_d_CONST_INT_LAST; k++)
	{
		h_CONST_INT[k] = m_vi_CONST_INT[k];
	}
	for(int k = 0; k<e_d_CONST_LAST; k++)
	{
		h_CONST[k] = m_vi_CONST[k];
	}
	printf("\n%d doubles and %d ints\n\n,",vci_numDoubles,vci_numInts);
	printf("\n%d double Bytes and %d int Bytes\n\n,",vci_numDoubleBytes,vci_numIntBytes);
	enum e_DATA {e_DATA_FAKE, e_DATA_CUDA, e_DATA_PREVIOUS, e_DATA_LAST};
	int e_v_DATA_OPTION = e_DATA_CUDA;
	if(e_v_DATA_OPTION==e_DATA_FAKE)
	{
		for(int i = 0; i<m_vi_CONST_INT[e_d_CONST_INT_numThreads]; i++)
		{
			for(int j = 0; j<m_vi_CONST_INT[e_d_CONST_INT_numRecordsPerThread]; j++)
			{
				for(int k = 0; k<e_d_IO_INT_LAST; k++)
				{
					int vi_PLACEHOLDER_LOCAL = (i*m_vi_CONST[e_d_CONST_INT_numRecordsPerThread]*e_d_IO_INT_LAST + j);
					vi_IO[vi_PLACEHOLDER_LOCAL] = 10000*i + 100*j + k;
				}
				for(int k = 0; k<e_d_IO_LAST; k++)
				{
					int vi_PLACEHOLDER_LOCAL = (i*m_vi_CONST[e_d_CONST_INT_numRecordsPerThread]*e_d_IO_LAST + j);
					vd_IO[vi_PLACEHOLDER_LOCAL] = 10000*i + 100*j + k;
				}
			}
		}
	}
	else if(e_v_DATA_OPTION==e_DATA_CUDA)
	{
		srand (time(NULL));
		double vd_zoutafp =  1.7;
		double vd_zstart  = -.06;
		cout << "\nSTART OF CUDA DATA OPTION";
		for(int i = 0; i<vci_numThreads; i++)
		{
			
			int vi_THREAD_INT_OFFSET = i*e_d_IO_INT_LAST;
			int vi_THREAD_DOUBLE_OFFSET = i*e_d_IO_LAST;
			
			int vi_IO_THREAD_INT_OFFSET = i*vci_numRecordsPerThread*e_d_IO_INT_LAST;
			int vi_IO_THREAD_DOUBLE_OFFSET = i*vci_numRecordsPerThread*e_d_IO_LAST;
			
			double vd_RANDOM_0 = ( (rand() % 32767) / 32767.0);
			double vd_RANDOM_1 = ( (rand() % 32767) / 32767.0);
			double vd_RANDOM_2 = ( (rand() % 32767) / 32767.0);
			double vd_RANDOM_3 = ( (rand() % 32767) / 32767.0);
			double vd_RANDOM_4 = ( (rand() % 32767) / 32767.0);
			double vd_RANDOM_5 = ( (rand() % 32767) / 32767.0);
			double vd_R_TEMP = m_vi_CONST[e_d_CONST_guide_1_radius] * vd_RANDOM_0;
			double vd_PHI_TEMP = m_vi_CONST[e_d_CONST_def_PI] * vd_RANDOM_1;
			double vd_THETA_TEMP = 2.0 * m_vi_CONST[e_d_CONST_def_PI] * vd_RANDOM_2;
			double vd_PHI_R_TEMP = m_vi_CONST[e_d_CONST_def_PI] * vd_RANDOM_3;
			double vd_THETA_R_TEMP = 2.0 * m_vi_CONST[e_d_CONST_def_PI] * vd_RANDOM_4;
			double vd_SX_TEMP = cos(vd_THETA_TEMP)*cos(vd_PHI_TEMP);
			double vd_SY_TEMP = sin(vd_THETA_TEMP)*cos(vd_PHI_TEMP);
			double vd_SZ_TEMP = sin(vd_PHI_TEMP);
			double vd_xin = vd_R_TEMP * cos(vd_THETA_R_TEMP);
			double vd_yin = vd_R_TEMP * sin(vd_THETA_R_TEMP);
			double vd_xout = 0;
			double vd_yout = 0;
			double vd_range = sqrt( (vd_xout-vd_xin)*(vd_xout-vd_xin)
				+(vd_yout-vd_yin)*(vd_yout-vd_yin)
				+(vd_zoutafp-vd_zstart)*(vd_zoutafp-vd_zstart));
			double vd_V_TEMP = 7.0 * pow( vd_RANDOM_3, (1.0/3.0) ); 
			double vd_VX_TEMP = vd_V_TEMP * (vd_xout-vd_xin)/vd_range;
			double vd_VY_TEMP = vd_V_TEMP * (vd_yout-vd_yin)/vd_range;
			double vd_VZ_TEMP = vd_V_TEMP * (vd_zoutafp-vd_zstart)/vd_range;
			
			vi_IO[vi_IO_THREAD_INT_OFFSET + e_d_IO_INT_RETURN_VALUE_XV] = 0;
			vi_IO[vi_IO_THREAD_INT_OFFSET + e_d_IO_INT_RETURN_VALUE_S] = 0; 
			vi_IO[vi_IO_THREAD_INT_OFFSET + e_d_IO_INT_ERROR] = 0;
			vi_IO[vi_IO_THREAD_INT_OFFSET + e_d_IO_INT_THREAD] = i;
			vi_IO[vi_IO_THREAD_INT_OFFSET + e_d_IO_INT_RECORD] = 0;
			vi_IO[vi_IO_THREAD_INT_OFFSET + e_d_IO_INT_RKQS_STEPS] = 0;
			vi_IO[vi_IO_THREAD_INT_OFFSET + e_d_IO_INT_RKQS_ERROR] = 0;
			
			vd_IO[vi_IO_THREAD_DOUBLE_OFFSET + e_d_IO_T] = 0.0;
			vd_IO[vi_IO_THREAD_DOUBLE_OFFSET + e_d_IO_X] = vd_xin;
			vd_IO[vi_IO_THREAD_DOUBLE_OFFSET + e_d_IO_Y] = vd_yin;
			vd_IO[vi_IO_THREAD_DOUBLE_OFFSET + e_d_IO_Z] = vd_zstart;
			vd_IO[vi_IO_THREAD_DOUBLE_OFFSET + e_d_IO_VX] = vd_VX_TEMP;
			vd_IO[vi_IO_THREAD_DOUBLE_OFFSET + e_d_IO_VY] = vd_VY_TEMP;
			vd_IO[vi_IO_THREAD_DOUBLE_OFFSET + e_d_IO_VZ] = vd_VZ_TEMP;
			vd_IO[vi_IO_THREAD_DOUBLE_OFFSET + e_d_IO_SPINNOR_0] = sqrt(
				pow(vd_SX_TEMP, 2.0) + 
				pow(vd_SY_TEMP, 2.0) + 
				pow(vd_SZ_TEMP, 2.0)
			);
			vd_IO[vi_IO_THREAD_DOUBLE_OFFSET + e_d_IO_SPINNOR_1] = vd_SX_TEMP;
			vd_IO[vi_IO_THREAD_DOUBLE_OFFSET + e_d_IO_SPINNOR_2] = vd_SY_TEMP;
			vd_IO[vi_IO_THREAD_DOUBLE_OFFSET + e_d_IO_SPINNOR_3] = vd_SZ_TEMP;
			
			vd_IO[vi_IO_THREAD_DOUBLE_OFFSET + e_d_IO_ERR_X] = 0.0; 
			vd_IO[vi_IO_THREAD_DOUBLE_OFFSET + e_d_IO_ERR_Y] = 0.0; 
			vd_IO[vi_IO_THREAD_DOUBLE_OFFSET + e_d_IO_ERR_Z] = 0.0; 
			vd_IO[vi_IO_THREAD_DOUBLE_OFFSET + e_d_IO_ERR_VX] = 0.0; 
			vd_IO[vi_IO_THREAD_DOUBLE_OFFSET + e_d_IO_ERR_VY] = 0.0; 
			vd_IO[vi_IO_THREAD_DOUBLE_OFFSET + e_d_IO_ERR_VZ] = 0.0; 
			vd_IO[vi_IO_THREAD_DOUBLE_OFFSET + e_d_IO_ERR_SPINNOR_0] = 0.0; 
			vd_IO[vi_IO_THREAD_DOUBLE_OFFSET + e_d_IO_ERR_SPINNOR_1] = 0.0; 
			vd_IO[vi_IO_THREAD_DOUBLE_OFFSET + e_d_IO_ERR_SPINNOR_2] = 0.0; 
			vd_IO[vi_IO_THREAD_DOUBLE_OFFSET + e_d_IO_ERR_SPINNOR_3] = 0.0;
			
			vd_IO[vi_IO_THREAD_DOUBLE_OFFSET + e_d_IO_HNEXT] = fmin(m_vi_CONST[e_d_CONST_h1], m_vi_CONST[e_d_CONST_h1_SPIN]);
			
			
			vd_IO[vi_IO_THREAD_DOUBLE_OFFSET + e_d_IO_SX] = vd_SX_TEMP; 
			vd_IO[vi_IO_THREAD_DOUBLE_OFFSET + e_d_IO_SY] = vd_SY_TEMP; 
 			vd_IO[vi_IO_THREAD_DOUBLE_OFFSET + e_d_IO_SZ] = 0.0; 
			
			printf("\nIO >>> OFFSET_INT:%d|OFFSET_DOUBLE:%d|ID:%d|RECORD:%d|TIME:%f|XYZ=(%f,%f,%f)|VEL=(%f,%f,%f)\nSPIN:(%f,%f,%f,%f)|HNEXT:%f",
				vi_IO_THREAD_INT_OFFSET,
				vi_IO_THREAD_DOUBLE_OFFSET,
				vi_IO[vi_IO_THREAD_INT_OFFSET + e_d_IO_INT_THREAD], 
				vi_IO[vi_IO_THREAD_INT_OFFSET + e_d_IO_INT_RECORD], 
				vd_IO[vi_IO_THREAD_DOUBLE_OFFSET + e_d_IO_T], 
				vd_IO[vi_IO_THREAD_DOUBLE_OFFSET + e_d_IO_X], vd_IO[vi_IO_THREAD_DOUBLE_OFFSET + e_d_IO_Y], vd_IO[vi_IO_THREAD_DOUBLE_OFFSET + e_d_IO_Z], 
				vd_IO[vi_IO_THREAD_DOUBLE_OFFSET + e_d_IO_VX], vd_IO[vi_IO_THREAD_DOUBLE_OFFSET + e_d_IO_VY], vd_IO[vi_IO_THREAD_DOUBLE_OFFSET + e_d_IO_VZ], 
				vd_IO[vi_IO_THREAD_DOUBLE_OFFSET + e_d_IO_SPINNOR_0], vd_IO[vi_IO_THREAD_DOUBLE_OFFSET + e_d_IO_SPINNOR_1], vd_IO[vi_IO_THREAD_DOUBLE_OFFSET + e_d_IO_SPINNOR_2], vd_IO[vi_IO_THREAD_DOUBLE_OFFSET + e_d_IO_SPINNOR_3], 
				vd_IO[vi_IO_THREAD_DOUBLE_OFFSET + e_d_IO_HNEXT]);
		}
		for(int i = 0; i<vci_numDoubles; i++) v_IO.push_back(vd_IO[i]);
		for(int i = 0; i<vci_numInts; i++) v_IO_INT.push_back(vi_IO[i]);
		// fxn_FORMATTED_WRITE("tempo", true);
		// fxn_FORMATTED_WRITE_ARRAY(vd_IO, vi_IO, "tempo_a1", false);
		v_IO.clear();
		v_IO_INT.clear();
		time(&rawtime);
		s_tm_timeinfo_START = gmtime ( &rawtime );
		GENERIC_WRAPPER_MULTI(
			h_CONST,
			h_CONST_INT,
			vd_IO,
			vi_IO);
		time(&rawtime);
		s_tm_timeinfo_FINISH = gmtime ( &rawtime );
		for(int i = 0; i<vci_numDoubles; i++) v_IO.push_back(vd_IO[i]);
		for(int i = 0; i<vci_numInts; i++) v_IO_INT.push_back(vi_IO[i]);
		cout << "\n final = " << e_d_IO_LAST << "\nmain vector is of size = (" << v_IO.size() << "\\" << v_IO_INT.size() << ")";
		
		// fxn_FORMATTED_WRITE_MINIMUM("UCN", true);
		
		fxn_FORMATTED_WRITE("UCN", true);
		// fxn_FORMATTED_WRITE_ARRAY(vd_IO, vi_IO, "tantrum_a1", false);
		
	}
	// else if(e_v_DATA_OPTION==e_DATA_PREVIOUS) {
		// ifstream ifp_INPUT;
		// char *cstring_INPUT_FILENAME = "UCN_OUTPUT_v_1_0_1_GMT_6_25_2014_2_29_43";
		// ifp_INPUT.open(cstring_INPUT_FILENAME, ofstream::in | ofstream::binary);
		// int *p_vi_vc_PH;
		// ifp_INPUT.read((char*) &p_vi_vc_PH,sizeof(int));
		// ifp_INPUT.read((char*) &p_vi_vc_PH,sizeof(int));
		// ifp_INPUT.read((char*) &p_vi_vc_PH,sizeof(int));
		
		// int vd_CONST_streamsize = sizeof(double)*e_d_CONST_LAST;
		// char h_CONST_vc_PH[vd_CONST_streamsize];
		// ifp_INPUT.read(h_CONST_vc_PH,vd_CONST_streamsize);
		// for(int vi_INDEX = 0; (vi_INDEX*sizeof(double))<vd_CONST_streamsize; vi_INDEX++){
			// double *p_vd_TEMP;
			// p_vd_TEMP = reinterpret_cast<double*>(&h_CONST_vc_PH[vi_INDEX*sizeof(double)]);
			// m_vi_CONST[vi_INDEX] = *p_vd_TEMP;
		// }
		// int vi_CONST_streamsize = sizeof(int)*e_d_CONST_INT_LAST;
		// char h_CONST_INT_vc_PH[vi_CONST_streamsize];
		// ifp_INPUT.read(h_CONST_INT_vc_PH,vi_CONST_streamsize);
		// for(int vi_INDEX = 0; (vi_INDEX*sizeof(int))<vi_CONST_streamsize; vi_INDEX++){
			// int *p_vi_TEMP;
			// p_vi_TEMP = reinterpret_cast<int*>(&h_CONST_INT_vc_PH[vi_INDEX*sizeof(int)]);
			// m_vi_CONST[vi_INDEX] = *p_vi_TEMP;
		// }
		// int vd_IO_streamsize = sizeof(double)*m_vi_CONST_INT[e_d_CONST_INT_numThreads]*m_vi_CONST_INT[e_d_CONST_INT_numRecordsPerThread]*e_d_IO_LAST;
		// char h_IO_vc_PH[vd_IO_streamsize];
		// ifp_INPUT.read(h_IO_vc_PH,vd_IO_streamsize);
		// for(int vi_INDEX = 0; (vi_INDEX*sizeof(double))<vd_IO_streamsize; vi_INDEX++){
			// double *p_vd_TEMP;
			// p_vd_TEMP = reinterpret_cast<double*>(&h_IO_vc_PH[vi_INDEX*sizeof(double)]);
			// vd_IO[vi_INDEX] = *p_vd_TEMP;
		// }
		// int vi_IO_streamsize = sizeof(int)*m_vi_CONST_INT[e_d_CONST_INT_numThreads]*m_vi_CONST_INT[e_d_CONST_INT_numRecordsPerThread]*e_d_IO_INT_LAST;
		// char h_IO_INT_vc_PH[vi_IO_streamsize];
		// ifp_INPUT.read(h_IO_INT_vc_PH,vi_IO_streamsize);
		// for(int vi_INDEX = 0; (vi_INDEX*sizeof(int))<vi_IO_streamsize; vi_INDEX++){
			// int *p_vi_TEMP;
			// p_vi_TEMP = reinterpret_cast<int*>(&h_IO_INT_vc_PH[vi_INDEX*sizeof(int)]);
			// vi_IO[vi_INDEX] = *p_vi_TEMP;
		// }
		// ifp_INPUT.close();
		
		// ifp_INPUT.open(cstring_INPUT_FILENAME, ifstream::in | ifstream::binary);
		
		// int vd_streamsize = sizeof(double)*m_vi_CONST_INT[e_d_CONST_INT_numThreads]*m_vi_CONST_INT[e_d_CONST_INT_numRecordsPerThread]*e_d_IO_LAST;
		// char vd_IO_vc_PH[vd_streamsize];h_IO
		// ifp_INPUT.read(vd_IO_vc_PH,vd_streamsize);
		// for(int vi_INDEX = 0; (vi_INDEX*sizeof(double))<vd_streamsize; vi_INDEX++){
			// double *p_vd_TEMP;
			// p_vd_TEMP = reinterpret_cast<double*>(&vd_IO_vc_PH[vi_INDEX*sizeof(double)]);
			// vd_IO[vi_INDEX] = *p_vd_TEMP;
		// }
		// ifp_INPUT.close();
	// }
	// char cstring_DOUBLE_INPUT_FILENAME[64];
	// sprintf (cstring_DOUBLE_INPUT_FILENAME,"INPUT_DOUBLE_v_1_0_0_GMT_%d_%d_%d_%d_%d_%d.bin",s_tm_timeinfo_FINISH->tm_mon+1,(s_tm_timeinfo_FINISH->tm_mday),(s_tm_timeinfo_FINISH->tm_year+1900),s_tm_timeinfo_FINISH->tm_hour,s_tm_timeinfo_FINISH->tm_min,s_tm_timeinfo_FINISH->tm_sec);
	// ofstream ofp_OUTPUT;
	// ofp_OUTPUT.open(cstring_DOUBLE_INPUT_FILENAME, ofstream::trunc | ofstream::out | ofstream::binary);
	// int vd_streamsize = sizeof(double)*m_vi_CONST_INT[e_d_CONST_INT_numThreads]*m_vi_CONST_INT[e_d_CONST_INT_numRecordsPerThread]*e_d_IO_LAST;
  // ofp_OUTPUT.write(reinterpret_cast<char*>(&vd_IO),vd_streamsize);
	// ofp_OUTPUT.close();
	// char cstring_INT_OUTPUT_FILENAME[64];
	// sprintf (cstring_INT_OUTPUT_FILENAME,"OUTPUT_INT_v_1_0_0_GMT_%d_%d_%d_%d_%d_%d.bin",s_tm_timeinfo_FINISH->tm_mon+1,(s_tm_timeinfo_FINISH->tm_mday),(s_tm_timeinfo_FINISH->tm_year+1900),s_tm_timeinfo_FINISH->tm_hour,s_tm_timeinfo_FINISH->tm_min,s_tm_timeinfo_FINISH->tm_sec);
	// ofp_OUTPUT.open(cstring_INT_OUTPUT_FILENAME, ofstream::trunc | ofstream::out | ofstream::binary);
	// int vi_streamsize = sizeof(int)*m_vi_CONST_INT[e_d_CONST_INT_numThreads]*m_vi_CONST_INT[e_d_CONST_INT_numRecordsPerThread]*e_d_IO_INT_LAST;
	// ofp_OUTPUT.write(reinterpret_cast<char*>(&vi_IO),vi_streamsize);
	// ofp_OUTPUT.close();

}
int c_UCN_RUN_DATA::fxn_V3_OUTPUT_CSV()
{
	// map<int, string>*> m_vi_cstr_IO = r_m_vi_cstr_IO;
	// map<int, string>*> m_vi_cstr_IO_INT = r_m_vi_cstr_IO_INT;
	// map<int, string>*> m_vi_cstr_CONST = r_m_vi_cstr_CONST;
	// map<int, string>*> m_vi_cstr_CONST_INT = r_m_vi_cstr_CONST_INT;
	// map<int, double> m_vi_CONST = r_m_vi_CONST;
	// map<int, int> m_vi_CONST_INT = r_m_vi_CONST_INT;
	static int param_PLACEHOLDER = 0;
	struct tm * s_tm_timeinfo_TEMP;
	int i, j, k;
	time_t rawtime;
	time(&rawtime);
	s_tm_timeinfo_TEMP = gmtime ( &rawtime );
	char cstring_DOUBLE_OUTPUT_FILENAME[64]; // = "tempo_tantrum.txt";
	sprintf (cstring_DOUBLE_OUTPUT_FILENAME,
		"UCN_OUTPUT_v_1_1_0_GMT_%d_%d_%d_%d_%d_%d_%d.csv",
		s_tm_timeinfo_TEMP->tm_mon+1,
		(s_tm_timeinfo_TEMP->tm_mday),
		(s_tm_timeinfo_TEMP->tm_year+1900),
		s_tm_timeinfo_TEMP->tm_hour,
		s_tm_timeinfo_TEMP->tm_min,
		s_tm_timeinfo_TEMP->tm_sec,
		param_PLACEHOLDER);
	cout << "\nATTENTION THIS IS TESTING OF VECTOR SIZE " << v_IO.size() << " AND " << v_IO_INT.size();
	cout << "\nATTENTION THIS IS i<" << m_vi_CONST_INT[e_d_CONST_INT_numThreads] << ", j<" << m_vi_CONST_INT[e_d_CONST_INT_numRecordsPerThread];
	cout << "\nATTENTION THIS IS CONST (" << m_vi_cstr_CONST.size() << ")/(" << e_d_CONST_LAST << ")";
	cout << "\nATTENTION THIS IS CONST (" << m_vi_cstr_CONST_INT.size() << ")/(" << e_d_CONST_INT_LAST << ")";
	cout << "\nATTENTION THIS IS CONST (" << m_vi_cstr_IO_INT.size() << ")/(" << e_d_IO_INT_LAST << ")";
	cout << "\nATTENTION THIS IS CONST (" << m_vi_cstr_IO.size() << ")/(" << e_d_IO_LAST << ")";
	cout << "\nATTENTION THIS IS CONST k<" << e_d_CONST_LAST << "/k<" << e_d_CONST_INT_LAST << " AND IO k<" << e_d_IO_LAST << "/k<" << e_d_IO_INT_LAST;
	// sprintf (cstring_DOUBLE_OUTPUT_FILENAME,
	// "OUT_%d_%d_%d_%d_%d_%d_%d.csv",
	// (s_tm_timeinfo_FINISH->tm_mon+1),
	// (s_tm_timeinfo_FINISH->tm_mday),
	// (s_tm_timeinfo_FINISH->tm_year+1900),
	// s_tm_timeinfo_FINISH->tm_hour,
	// s_tm_timeinfo_FINISH->tm_min, 
	// s_tm_timeinfo_FINISH->tm_sec,
	// param_PLACEHOLDER);
	// sprintf (cstring_DOUBLE_OUTPUT_FILENAME,
	// "OUT_v_1_0_0_%d.csv",
	// param_PLACEHOLDER);
	ofstream ofp_OUTPUT;
	ofp_OUTPUT.open(cstring_DOUBLE_OUTPUT_FILENAME, ofstream::trunc | ofstream::out);
	bool vb_START_FLAG = true;
	// for(k = 0; k<e_d_CONST_LAST; k++)
	// {
	 // cout << "\n" << m_vi_cstr_CONST[k] << ",";
	// }
	// cout << "WHAT THE FUCK IS GOING ON?????";
	for(i = 0; i<m_vi_CONST_INT[e_d_CONST_INT_numThreads]; i++)
	{
		
		for(j = 0; j<m_vi_CONST_INT[e_d_CONST_INT_numRecordsPerThread]; j++)
		{
			ofp_OUTPUT << "*";
			if(vb_START_FLAG)
			{
				for(k = 0; k<e_d_CONST_INT_LAST; k++)
				{
				 ofp_OUTPUT << m_vi_cstr_CONST_INT[k] << ",";
				}
				for(k = 0; k<e_d_CONST_LAST; k++)
				{
				 ofp_OUTPUT << m_vi_cstr_CONST[k] << ",";
				}
				ofp_OUTPUT << "\n";
				for(int k = 0; k<e_d_CONST_INT_LAST; k++)
				{
				 ofp_OUTPUT << m_vi_CONST_INT[k] << ",";
				}
				for(k = 0; k<e_d_CONST_LAST; k++)
				{
				 ofp_OUTPUT << m_vi_CONST[k] << ",";
				}
				ofp_OUTPUT << "\n";
				for(k = 0; k<e_d_IO_INT_LAST; k++)
				{
				 ofp_OUTPUT << m_vi_cstr_IO_INT[k] << ",";
				}
				for(k = 0; k<e_d_IO_LAST; k++)
				{
				 ofp_OUTPUT << m_vi_cstr_IO[k] << ",";
				}
				ofp_OUTPUT << "\n";
				vb_START_FLAG = false;
			}
			for(k = 0; k<e_d_IO_INT_LAST; k++) 
			{
				// ofp_OUTPUT << "h";
				ofp_OUTPUT << v_IO_INT[i*m_vi_CONST_INT[e_d_CONST_INT_numRecordsPerThread]*e_d_IO_INT_LAST+j*e_d_IO_INT_LAST+k] << ",";
			}
			for(k = 0; k<e_d_IO_LAST; k++)
			{
				// ofp_OUTPUT << "H";
				ofp_OUTPUT << v_IO[i*m_vi_CONST_INT[e_d_CONST_INT_numRecordsPerThread]*e_d_IO_LAST+j*e_d_IO_LAST+k] << ",";
			}
			ofp_OUTPUT << "\n";
		}
	}
	ofp_OUTPUT.close();
	param_PLACEHOLDER++;
	return 0;
}
double c_UCN_RUN_DATA::fxn_vd_IO(int param_THREAD, int param_RECORD, int param_INDEX)
{
	int vi_THREAD_OFFSET = param_THREAD*fvci_h_CONST(e_d_CONST_INT_numRecordsPerThread)*e_d_IO_LAST;
	int vi_RECORD_OFFSET = param_RECORD*e_d_IO_LAST;
	int vi_TOTAL_OFFSET = vi_THREAD_OFFSET + vi_RECORD_OFFSET + param_INDEX;
	return v_IO[vi_TOTAL_OFFSET];
}
int c_UCN_RUN_DATA::fxn_vi_IO(int param_THREAD, int param_RECORD, int param_INDEX)
{
	int vi_THREAD_OFFSET = param_THREAD*m_vi_CONST[e_d_CONST_INT_numRecordsPerThread]*e_d_IO_INT_LAST;
	int vi_RECORD_OFFSET = param_RECORD*e_d_IO_INT_LAST;
	int vi_TOTAL_OFFSET = vi_THREAD_OFFSET + vi_RECORD_OFFSET + param_INDEX;
	return v_IO_INT[vi_TOTAL_OFFSET];
}
const int c_UCN_RUN_DATA::fvci_h_CONST(int vi_INDEX)
{
	return m_vi_CONST_INT[vi_INDEX];
}
const double c_UCN_RUN_DATA::fvcd_h_CONST(int vi_INDEX) 
{
	return m_vi_CONST[vi_INDEX];
}

int c_UCN_RUN_DATA::fxn_FORMATTED_WRITE(const char *vc_a1_FILENAME, bool vb_TIMESTAMP_FLAG)
{
	static int vsi_PLACEHOLDER = 0;
	char cstring_OUTPUT_FILENAME[64]; // = "tempo_tantrum.txt";
	if(vb_TIMESTAMP_FLAG)
	{
		time_t l_rawtime;
		time(&l_rawtime);
		struct tm * s_tm_timeinfo_TEMP = gmtime(&l_rawtime);
		sprintf (cstring_OUTPUT_FILENAME,
		"%s_%d_%d_%d_%d_%d_%d_%d.csv",
		vc_a1_FILENAME,
		s_tm_timeinfo_TEMP->tm_mon+1,
		(s_tm_timeinfo_TEMP->tm_mday),
		(s_tm_timeinfo_TEMP->tm_year+1900),
		s_tm_timeinfo_TEMP->tm_hour,
		s_tm_timeinfo_TEMP->tm_min,
		s_tm_timeinfo_TEMP->tm_sec,
		vsi_PLACEHOLDER);
		vsi_PLACEHOLDER++;
	}
	else
	{
		sprintf (cstring_OUTPUT_FILENAME,
			"%s.csv",
			vc_a1_FILENAME);
	}
	
	ofstream ofp_INPUT;
	ofp_INPUT.open(cstring_OUTPUT_FILENAME, ofstream::trunc | ofstream::out);
	ofp_INPUT << std::fixed;
	ofp_INPUT.precision(20);
	for(int i = 0; i<m_vi_CONST_INT[e_d_CONST_INT_numThreads]; i++)
	{ 
		for(int j = 0; j<m_vi_CONST_INT[e_d_CONST_INT_numRecordsPerThread]; j++)
		{ 
			if(i==0 && j==0)
			{
				for(int k = 0; k<e_d_CONST_INT_LAST; k++)
				{
					ofp_INPUT << m_vi_cstr_CONST_INT[k] << ",";
				}
				for(int k = 0; k<e_d_CONST_LAST; k++)
				{
					ofp_INPUT << m_vi_cstr_CONST[k] << ",";
				}
				ofp_INPUT << "\n";
				for(int k = 0; k<e_d_CONST_INT_LAST; k++)
				{
					ofp_INPUT << m_vi_CONST_INT[k] << ",";
				}
				for(int k = 0; k<e_d_CONST_LAST; k++)
				{
					ofp_INPUT << m_vi_CONST[k] << ",";
				}
				ofp_INPUT << "\n";
				for(int k = 0; k<e_d_IO_INT_LAST; k++)
				{
					ofp_INPUT << m_vi_cstr_IO_INT[k] << ",";
				}
				for(int k = 0; k<e_d_IO_LAST; k++)
				{
					ofp_INPUT << m_vi_cstr_IO[k] << ",";
				}
				ofp_INPUT << "\n";
			}
			for(int k = 0; k<e_d_IO_INT_LAST; k++)
			{
				int vi_INDEX_TEMP = (i*m_vi_CONST_INT[e_d_CONST_INT_numRecordsPerThread]+j)*e_d_IO_INT_LAST+k;
				if(k==e_d_IO_INT_RKQS_ERROR)
				{
					ofp_INPUT << m_vi_cstr_RKQS_ERROR[(v_IO_INT[vi_INDEX_TEMP])] << ",";
				}
				else 
				{
					ofp_INPUT << v_IO_INT[vi_INDEX_TEMP] << ",";
				}
			}
			for(int k = 0; k<e_d_IO_LAST; k++)
			{
				int vi_INDEX_TEMP = (i*m_vi_CONST_INT[e_d_CONST_INT_numRecordsPerThread]+j)*e_d_IO_LAST+k;
				ofp_INPUT << v_IO[vi_INDEX_TEMP] << ",";
			}
			ofp_INPUT << "\n";
		}
	}
	ofp_INPUT.close();
	return 0;
}

int c_UCN_RUN_DATA::fxn_FORMATTED_WRITE_ARRAY(
	double *vd_IO, 
	int *vi_IO, 
	const char *vc_a1_FILENAME, 
	bool vb_TIMESTAMP_FLAG)
{
	static int vsi_PLACEHOLDER = 0;
	char cstring_OUTPUT_FILENAME[64]; // = "tempo_tantrum.txt";
	if(vb_TIMESTAMP_FLAG)
	{
		time_t l_rawtime;
		time(&l_rawtime);
		struct tm * s_tm_timeinfo_TEMP = gmtime(&l_rawtime);
		sprintf (cstring_OUTPUT_FILENAME,
		"UCN_v_1_2_0_GMT_%s_%d_%d_%d_%d_%d_%d_%d.csv",
		vc_a1_FILENAME,
		s_tm_timeinfo_TEMP->tm_mon+1,
		(s_tm_timeinfo_TEMP->tm_mday),
		(s_tm_timeinfo_TEMP->tm_year+1900),
		s_tm_timeinfo_TEMP->tm_hour,
		s_tm_timeinfo_TEMP->tm_min,
		s_tm_timeinfo_TEMP->tm_sec,
		vsi_PLACEHOLDER);
		vsi_PLACEHOLDER++;
	}
	else
	{
		sprintf (cstring_OUTPUT_FILENAME,
			"%s.csv",
			vc_a1_FILENAME);
	}
	ofstream ofp_INPUT;
	ofp_INPUT.open(cstring_OUTPUT_FILENAME, ofstream::trunc | ofstream::out);
	for(int i = 0; i<m_vi_CONST_INT[e_d_CONST_INT_numThreads]; i++)
	{ 
		for(int j = 0; j<m_vi_CONST_INT[e_d_CONST_INT_numRecordsPerThread]; j++)
		{ 
			if(i==0 && j==0)
			{
				for(int k = 0; k<e_d_CONST_INT_LAST; k++)
				{
					ofp_INPUT << m_vi_cstr_CONST_INT[k] << ",";
				}
				for(int k = 0; k<e_d_CONST_LAST; k++)
				{
					ofp_INPUT << m_vi_cstr_CONST[k] << ",";
				}
				ofp_INPUT << "\n";
				for(int k = 0; k<e_d_CONST_INT_LAST; k++)
				{
					ofp_INPUT << m_vi_CONST_INT[k] << ",";
				}
				for(int k = 0; k<e_d_CONST_LAST; k++)
				{
					ofp_INPUT << m_vi_CONST[k] << ",";
				}
				ofp_INPUT << "\n";
				for(int k = 0; k<e_d_IO_INT_LAST; k++)
				{
					ofp_INPUT << m_vi_cstr_IO_INT[k] << ",";
				}
				for(int k = 0; k<e_d_IO_LAST; k++)
				{
					ofp_INPUT << m_vi_cstr_IO[k] << ",";
				}
				ofp_INPUT << "\n";
			}
			for(int k = 0; k<e_d_IO_INT_LAST; k++)
			{
				int vi_INDEX_TEMP = (i*m_vi_CONST_INT[e_d_CONST_INT_numRecordsPerThread]+j)*e_d_IO_INT_LAST+k;
				ofp_INPUT << vi_IO[vi_INDEX_TEMP] << ",";
			}
			for(int k = 0; k<e_d_IO_LAST; k++)
			{
				int vi_INDEX_TEMP = (i*m_vi_CONST_INT[e_d_CONST_INT_numRecordsPerThread]+j)*e_d_IO_LAST+k;
				ofp_INPUT << vd_IO[vi_INDEX_TEMP] << ",";
			}
			ofp_INPUT << "\n";
		}
	}
	ofp_INPUT.close();
	return 0;
}

int c_UCN_RUN_DATA::fxn_FORMATTED_WRITE_MINIMUM(const char *vc_a1_FILENAME, bool vb_TIMESTAMP_FLAG)
{
	static int vsi_PLACEHOLDER = 0;
	char cstring_OUTPUT_FILENAME[64]; // = "tempo_tantrum.txt";
	if(vb_TIMESTAMP_FLAG)
	{
		time_t l_rawtime;
		time(&l_rawtime);
		struct tm * s_tm_timeinfo_TEMP = gmtime(&l_rawtime);
		sprintf (cstring_OUTPUT_FILENAME,
		"%s_%d_%d_%d_%d_%d_%d_%d.csv",
		vc_a1_FILENAME,
		s_tm_timeinfo_TEMP->tm_mon+1,
		(s_tm_timeinfo_TEMP->tm_mday),
		(s_tm_timeinfo_TEMP->tm_year+1900),
		s_tm_timeinfo_TEMP->tm_hour,
		s_tm_timeinfo_TEMP->tm_min,
		s_tm_timeinfo_TEMP->tm_sec,
		vsi_PLACEHOLDER);
		vsi_PLACEHOLDER++;
	}
	else
	{
		sprintf (cstring_OUTPUT_FILENAME,
			"%s.csv",
			vc_a1_FILENAME);
	}
	ofstream ofp_INPUT;
	ofp_INPUT.open(cstring_OUTPUT_FILENAME, ofstream::trunc | ofstream::out);
	ofp_INPUT << std::fixed;
	ofp_INPUT.precision(20);
	for(int i = 0; i<m_vi_CONST_INT[e_d_CONST_INT_numThreads]; i++)
	{ 
		for(int j = 0; j<m_vi_CONST_INT[e_d_CONST_INT_numRecordsPerThread]; j++)
		{ 
			if(i==0 && j==0)
			{
				ofp_INPUT << m_vi_cstr_IO_INT[e_d_IO_INT_THREAD] << ","; 
				ofp_INPUT << m_vi_cstr_IO_INT[e_d_IO_INT_RECORD] << ","; 
				ofp_INPUT << m_vi_cstr_IO[e_d_IO_T] << ","; 
				ofp_INPUT << m_vi_cstr_IO[e_d_IO_X] << ","; 
				ofp_INPUT << m_vi_cstr_IO[e_d_IO_Y] << ","; 
				ofp_INPUT << m_vi_cstr_IO[e_d_IO_Z] << ","; 
				ofp_INPUT << m_vi_cstr_IO[e_d_IO_VX] << ","; 
				ofp_INPUT << m_vi_cstr_IO[e_d_IO_VY] << ","; 
				ofp_INPUT << m_vi_cstr_IO[e_d_IO_VZ] << ","; 
				ofp_INPUT << m_vi_cstr_IO[e_d_IO_SPINNOR_1] << ","; 
				ofp_INPUT << m_vi_cstr_IO[e_d_IO_SPINNOR_2] << ","; 
				ofp_INPUT << m_vi_cstr_IO[e_d_IO_SPINNOR_3] << ","; 
				ofp_INPUT << m_vi_cstr_IO[e_d_IO_BX] << ","; 
				ofp_INPUT << m_vi_cstr_IO[e_d_IO_BY] << ","; 
				ofp_INPUT << m_vi_cstr_IO[e_d_IO_BZ] << ","; 
			}
			int vi_INDEX_TEMP = (i*m_vi_CONST_INT[e_d_CONST_INT_numRecordsPerThread]+j)*e_d_IO_INT_LAST+ e_d_IO_INT_THREAD;
			ofp_INPUT << v_IO_INT[vi_INDEX_TEMP] << ",";
			vi_INDEX_TEMP = (i*m_vi_CONST_INT[e_d_CONST_INT_numRecordsPerThread]+j)*e_d_IO_INT_LAST+ e_d_IO_INT_RECORD;
			ofp_INPUT << v_IO_INT[vi_INDEX_TEMP] << ",";
			vi_INDEX_TEMP = (i*m_vi_CONST_INT[e_d_CONST_INT_numRecordsPerThread]+j)*e_d_IO_LAST + e_d_IO_T; 
			ofp_INPUT << v_IO[vi_INDEX_TEMP] << ",";
			vi_INDEX_TEMP = (i*m_vi_CONST_INT[e_d_CONST_INT_numRecordsPerThread]+j)*e_d_IO_LAST + e_d_IO_X; 
			ofp_INPUT << v_IO[vi_INDEX_TEMP] << ",";
			vi_INDEX_TEMP = (i*m_vi_CONST_INT[e_d_CONST_INT_numRecordsPerThread]+j)*e_d_IO_LAST + e_d_IO_Y; 
			ofp_INPUT << v_IO[vi_INDEX_TEMP] << ",";
			vi_INDEX_TEMP = (i*m_vi_CONST_INT[e_d_CONST_INT_numRecordsPerThread]+j)*e_d_IO_LAST + e_d_IO_Z; 
			ofp_INPUT << v_IO[vi_INDEX_TEMP] << ",";
			ofp_INPUT << v_IO[vi_INDEX_TEMP] << ",";
			vi_INDEX_TEMP = (i*m_vi_CONST_INT[e_d_CONST_INT_numRecordsPerThread]+j)*e_d_IO_LAST + e_d_IO_VX; 
			ofp_INPUT << v_IO[vi_INDEX_TEMP] << ",";
			vi_INDEX_TEMP = (i*m_vi_CONST_INT[e_d_CONST_INT_numRecordsPerThread]+j)*e_d_IO_LAST + e_d_IO_VY; 
			ofp_INPUT << v_IO[vi_INDEX_TEMP] << ",";
			vi_INDEX_TEMP = (i*m_vi_CONST_INT[e_d_CONST_INT_numRecordsPerThread]+j)*e_d_IO_LAST + e_d_IO_VZ; 
			ofp_INPUT << v_IO[vi_INDEX_TEMP] << ",";
			vi_INDEX_TEMP = (i*m_vi_CONST_INT[e_d_CONST_INT_numRecordsPerThread]+j)*e_d_IO_LAST + e_d_IO_SPINNOR_1; 
			ofp_INPUT << v_IO[vi_INDEX_TEMP] << ",";
			vi_INDEX_TEMP = (i*m_vi_CONST_INT[e_d_CONST_INT_numRecordsPerThread]+j)*e_d_IO_LAST + e_d_IO_SPINNOR_2; 
			ofp_INPUT << v_IO[vi_INDEX_TEMP] << ",";
			vi_INDEX_TEMP = (i*m_vi_CONST_INT[e_d_CONST_INT_numRecordsPerThread]+j)*e_d_IO_LAST + e_d_IO_SPINNOR_3; 
			ofp_INPUT << v_IO[vi_INDEX_TEMP] << ",";
			vi_INDEX_TEMP = (i*m_vi_CONST_INT[e_d_CONST_INT_numRecordsPerThread]+j)*e_d_IO_LAST + e_d_IO_BX; 
			ofp_INPUT << v_IO[vi_INDEX_TEMP] << ",";
			vi_INDEX_TEMP = (i*m_vi_CONST_INT[e_d_CONST_INT_numRecordsPerThread]+j)*e_d_IO_LAST + e_d_IO_BY; 
			ofp_INPUT << v_IO[vi_INDEX_TEMP] << ",";
			vi_INDEX_TEMP = (i*m_vi_CONST_INT[e_d_CONST_INT_numRecordsPerThread]+j)*e_d_IO_LAST + e_d_IO_BZ; 
			ofp_INPUT << v_IO[vi_INDEX_TEMP] << ",";
			ofp_INPUT << "\n";
		}
	}
	ofp_INPUT.close();
	return 0;
}


int fxn_STANDARD_INPUT_DATA(
	map<int, string> & r_m_vi_cstr_IO_INPUT, 
	map<int, string> & r_m_vi_cstr_IO_INT_INPUT, 
	map<int, string> & r_m_vi_cstr_CONST_INPUT, 
	map<int, string> & r_m_vi_cstr_CONST_INT_INPUT, 
	map<int, string> & r_m_vi_cstr_RKQS_ERROR, 
	map<int, int> & r_m_vi_CONST_INT, 
	map<int, double> & r_m_vi_CONST, 
	const int param_numThreads, 
	int *param_a1_INPUT_INT, //[param_numThreads*e_d_INPUT_INT_LAST], 
	double *param_a1_INPUT) //[param_numThreads*e_d_INPUT_LAST])
{
	srand (time(NULL));
	double vd_zoutafp =  1.7;
	double vd_zstart  = -.06;
	// unsigned vi_seed1 = std::chrono::system_clock::now().time_since_epoch().count();
	// std::mt19937 mt19937_LOCAL (vi_seed1);
	// std::uniform_real_distribution<double> distribution_vd_0_1(0.0,1.0);
	for(int i = 0; i<r_m_vi_CONST_INT[e_d_CONST_INT_numThreads]; i++)
	{
		double vd_RANDOM_0 = ( (rand() % 32767) / 32767.0);
		double vd_RANDOM_1 = ( (rand() % 32767) / 32767.0);
		double vd_RANDOM_2 = ( (rand() % 32767) / 32767.0);
		double vd_RANDOM_3 = ( (rand() % 32767) / 32767.0);
		double vd_R_TEMP = r_m_vi_CONST[e_d_CONST_guide_1_radius] * vd_RANDOM_0;
		double vd_PHI_TEMP = 2.0 * r_m_vi_CONST[e_d_CONST_def_PI] * vd_RANDOM_1;
		double vd_THETA_TEMP = 2.0 * r_m_vi_CONST[e_d_CONST_def_PI] * vd_RANDOM_2;
		double vd_SX_TEMP = cos(vd_THETA_TEMP);
		double vd_SY_TEMP = sin(vd_THETA_TEMP);
		double vd_xin = vd_R_TEMP * cos(vd_PHI_TEMP);
		double vd_yin = vd_R_TEMP * sin(vd_PHI_TEMP);
		double vd_xout = 0;
		double vd_yout = 0;
		double vd_range = sqrt( (vd_xout-vd_xin)*(vd_xout-vd_xin)
			+(vd_yout-vd_yin)*(vd_yout-vd_yin)
			+(vd_zoutafp-vd_zstart)*(vd_zoutafp-vd_zstart));
		double vd_V_TEMP = 7.0 * pow( vd_RANDOM_3, (1.0/3.0) ); 
		double vd_VX_TEMP = vd_V_TEMP * (vd_xout-vd_xin)/vd_range;
		double vd_VY_TEMP = vd_V_TEMP * (vd_yout-vd_yin)/vd_range;
		double vd_VZ_TEMP = vd_V_TEMP * (vd_zoutafp-vd_zstart)/vd_range;
		
		// cout << "\n";
		// cout << 0.0;
		// cout << vd_xin; 
		// cout << vd_yin; 
		// cout << vd_zstart; 
		// cout << vd_VX_TEMP; 
		// cout << vd_VY_TEMP; 
		// cout << vd_VZ_TEMP; 
		// cout << 0.0; 
		// cout << 0.0; 
		// cout << 1.0; 
		int vi_THREAD_INT_OFFSET = i*e_d_IO_INT_LAST;
		int vi_THREAD_DOUBLE_OFFSET = i*e_d_IO_LAST;
		
		param_a1_INPUT_INT[vi_THREAD_INT_OFFSET + e_d_INPUT_INT_THREAD] = i;
		param_a1_INPUT_INT[vi_THREAD_INT_OFFSET + e_d_INPUT_INT_RECORD] = 0;
		
		param_a1_INPUT[vi_THREAD_DOUBLE_OFFSET + e_d_INPUT_T] = 0.0;
		param_a1_INPUT[vi_THREAD_DOUBLE_OFFSET + e_d_INPUT_X] = vd_xin; 
		param_a1_INPUT[vi_THREAD_DOUBLE_OFFSET + e_d_INPUT_Y] = vd_yin; 
		param_a1_INPUT[vi_THREAD_DOUBLE_OFFSET + e_d_INPUT_Z] = vd_zstart; 
		param_a1_INPUT[vi_THREAD_DOUBLE_OFFSET + e_d_INPUT_VX] = vd_VX_TEMP; 
		param_a1_INPUT[vi_THREAD_DOUBLE_OFFSET + e_d_INPUT_VY] = vd_VY_TEMP; 
		param_a1_INPUT[vi_THREAD_DOUBLE_OFFSET + e_d_INPUT_VZ] = vd_VZ_TEMP; 
		param_a1_INPUT[vi_THREAD_DOUBLE_OFFSET + e_d_INPUT_SPINNOR_0] = 1.0; 
		param_a1_INPUT[vi_THREAD_DOUBLE_OFFSET + e_d_INPUT_SPINNOR_1] = (1.0/sqrt(2.0)); 
		param_a1_INPUT[vi_THREAD_DOUBLE_OFFSET + e_d_INPUT_SPINNOR_2] = (1.0/sqrt(2.0)); 
		param_a1_INPUT[vi_THREAD_DOUBLE_OFFSET + e_d_INPUT_SPINNOR_3] =  0.0;
		param_a1_INPUT[vi_THREAD_DOUBLE_OFFSET + e_d_INPUT_HNEXT] = fmin(r_m_vi_CONST[e_d_CONST_h1], r_m_vi_CONST[e_d_CONST_h1_SPIN]);
		
		printf("\nID:%d|RECORD:%d|TIME:%f|XYZ=(%f,%f,%f)|VEL=(%f,%f,%f)|SPIN:(%f,%f,%f,%f)|HNEXT:%f",
			param_a1_INPUT_INT[vi_THREAD_INT_OFFSET + e_d_INPUT_INT_THREAD], 
			param_a1_INPUT_INT[vi_THREAD_INT_OFFSET + e_d_INPUT_INT_RECORD], 
			param_a1_INPUT[vi_THREAD_DOUBLE_OFFSET + e_d_INPUT_T], 
			param_a1_INPUT[vi_THREAD_DOUBLE_OFFSET + e_d_INPUT_X], param_a1_INPUT[vi_THREAD_DOUBLE_OFFSET + e_d_INPUT_Y], param_a1_INPUT[vi_THREAD_DOUBLE_OFFSET + e_d_INPUT_Z], 
			param_a1_INPUT[vi_THREAD_DOUBLE_OFFSET + e_d_INPUT_VX], param_a1_INPUT[vi_THREAD_DOUBLE_OFFSET + e_d_INPUT_VY], param_a1_INPUT[vi_THREAD_DOUBLE_OFFSET + e_d_INPUT_VZ], 
			param_a1_INPUT[vi_THREAD_DOUBLE_OFFSET + e_d_INPUT_SPINNOR_0], param_a1_INPUT[vi_THREAD_DOUBLE_OFFSET + e_d_INPUT_SPINNOR_1], param_a1_INPUT[vi_THREAD_DOUBLE_OFFSET + e_d_INPUT_SPINNOR_2], param_a1_INPUT[vi_THREAD_DOUBLE_OFFSET + e_d_INPUT_SPINNOR_3], 
			param_a1_INPUT[vi_THREAD_DOUBLE_OFFSET + e_d_INPUT_HNEXT]);
	}
	return 0;
}
c_UCN_RUN_DATA fxn_PARAMETERIZED_CONSTRUCTOR_c_UCN_RUN_DATA(
	const int param_e_d_CONST_INT_numThreadsPerBlock_EXP_MAX, 
	const int param_e_d_CONST_INT_numRecordsPerThread_EXP, 
	const int param_e_d_CONST_INT_numThreads_EXP, 
	const int param_e_d_CONST_INT_numRecordsPerUpdate, 
	const int param_e_d_CONST_INT_numCyclesPerRecord,
	map<int, int> param_m_vi_CONST_INT, 
	map<int, double> param_m_vi_CONST)
{ 
	const int vi_INPUT_e_d_CONST_INT_numThreadsPerBlock_EXP_MAX = 
		( param_e_d_CONST_INT_numThreadsPerBlock_EXP_MAX>=0 ? param_e_d_CONST_INT_numThreadsPerBlock_EXP_MAX : 0 );
	const int vi_INPUT_e_d_CONST_INT_numRecordsPerThread_EXP = 
		( param_e_d_CONST_INT_numRecordsPerThread_EXP>=0 ? param_e_d_CONST_INT_numRecordsPerThread_EXP : 0 );
	const int vi_INPUT_e_d_CONST_INT_numThreads_EXP = 
		( param_e_d_CONST_INT_numThreads_EXP>=0 ? param_e_d_CONST_INT_numThreads_EXP : 0 );
	const int vi_INPUT_e_d_CONST_INT_numRecordsPerUpdate = 
		( param_e_d_CONST_INT_numRecordsPerUpdate>=1 ? param_e_d_CONST_INT_numRecordsPerUpdate : 1 );
	const int vi_INPUT_e_d_CONST_INT_numCyclesPerRecord = 
		( param_e_d_CONST_INT_numCyclesPerRecord>=1 ? param_e_d_CONST_INT_numCyclesPerRecord : 1 );
	
	cout << vi_INPUT_e_d_CONST_INT_numThreads_EXP << "THIS CORRECT????";
	map<int, string> m_vi_cstr_IO;
	map<int, string> m_vi_cstr_IO_INT;
	map<int, string> m_vi_cstr_CONST;
	map<int, string> m_vi_cstr_CONST_INT;
	map<int, string> m_vi_cstr_RKQS_ERROR;
	map<int, int> m_vi_CONST_INT;
	map<int, int> m_vi_CONST_INT_DEFAULT;
	map<int, double> m_vi_CONST;
	map<int, double> m_vi_CONST_DEFAULT;
	
	m_vi_cstr_RKQS_ERROR.insert(std::pair<int,string>(e_RKQS_ERROR_NONE,"NONE"));
	m_vi_cstr_RKQS_ERROR.insert(std::pair<int,string>(e_RKQS_ERROR_STEPSIZE_UNDERFLOW,"STEPSIZE_UNDERFLOW"));
	m_vi_cstr_RKQS_ERROR.insert(std::pair<int,string>(e_RKQS_ERROR_XV_BOUNDS,"XV_BOUNDS"));
	m_vi_cstr_RKQS_ERROR.insert(std::pair<int,string>(e_RKQS_ERROR_SPIN_BOUNDS,"SPIN_BOUNDS"));
	m_vi_cstr_RKQS_ERROR.insert(std::pair<int,string>(e_RKQS_ERROR_COMBINED_BOUNDS,"COMBINED_BOUNDS"));
	m_vi_cstr_RKQS_ERROR.insert(std::pair<int,string>(e_RKQS_ERROR_REVERSED_INTERVAL,"REVERSED_INTERVAL"));
	m_vi_cstr_RKQS_ERROR.insert(std::pair<int,string>(e_RKQS_ERROR_UNKNOWN,"UNKNOWN"));
	
	m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_T,"T"));
	m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_X,"X"));
	m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_Y,"Y"));
	m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_Z,"Z"));
	m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_VX,"VX"));
	m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_VY,"VY"));
	m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_VZ,"VZ"));
	m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_SPINNOR_0,"SPINNOR_0"));
	m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_SPINNOR_1,"SPINNOR_1"));
	m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_SPINNOR_2,"SPINNOR_2"));
	m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_SPINNOR_3,"SPINNOR_3"));
	m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_RED_VX,"RED_VX"));
	m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_RED_VY,"RED_VY"));
	m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_RED_VZ,"RED_VZ"));  
	m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_AX,"AX"));
	m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_AY,"AY"));
	m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_AZ,"AZ"));
	m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_DDT_SPINNOR_0,"DDT_SPINNOR_0")); 
  m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_DDT_SPINNOR_1,"DDT_SPINNOR_1")); 
  m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_DDT_SPINNOR_2,"DDT_SPINNOR_2")); 
  m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_DDT_SPINNOR_3,"DDT_SPINNOR_3")); 
	m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_SCAL_X,"SCAL_X"));
	m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_SCAL_Y,"SCAL_Y"));
	m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_SCAL_Z,"SCAL_Z"));
	m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_SCAL_VX,"SCAL_VX"));
	m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_SCAL_VY,"SCAL_VY"));
	m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_SCAL_VZ,"SCAL_VZ"));
	m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_SCAL_SPINNOR_0,"SCAL_SPINNOR_0")); 
  m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_SCAL_SPINNOR_1,"SCAL_SPINNOR_1")); 
  m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_SCAL_SPINNOR_2,"SCAL_SPINNOR_2")); 
  m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_SCAL_SPINNOR_3,"SCAL_SPINNOR_3"));
	m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_ERR_X,"ERR_X"));
	m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_ERR_Y,"ERR_Y"));
	m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_ERR_Z,"ERR_Z"));
	m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_ERR_VX,"ERR_VX"));
	m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_ERR_VY,"ERR_VY"));
	m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_ERR_VZ,"ERR_VZ"));
	m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_ERR_SPINNOR_0,"ERR_SPINNOR_0")); 
  m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_ERR_SPINNOR_1,"ERR_SPINNOR_1")); 
  m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_ERR_SPINNOR_2,"ERR_SPINNOR_2")); 
  m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_ERR_SPINNOR_3,"ERR_SPINNOR_3")); 
	m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_BX,"BX"));
	m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_BY,"BY"));
	m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_BZ,"BZ"));
	m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_dB_XDX,"dB_XDX"));
	m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_dB_XDY,"dB_XDY"));
	m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_dB_XDZ,"dB_XDZ"));
	m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_dB_YDX,"dB_YDX"));
	m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_dB_YDY,"dB_YDY"));
	m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_dB_YDZ,"dB_YDZ"));
	m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_dB_ZDX,"dB_ZDX"));
	m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_dB_ZDY,"dB_ZDY"));
	m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_dB_ZDZ,"dB_ZDZ"));
	m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_POLARIZATION,"POLARIZATION"));
  m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_HNEXT,"HNEXT")); 
	m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_SX,"SX"));
	m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_SY,"SY"));
	m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_SZ,"SZ"));
	// m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_EPSILON_MAX_XV,"EPSILON_MAX_XV")); 
  // m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_EPSILON_MAX_S,"EPSILON_MAX_S")); 
  // m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_ERROR_XV,"ERROR_XV")); 
  // m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_ERROR_S,"ERROR_S")); 
	// m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_HNEXT_XV,"HNEXT_XV")); 
	// m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_HNEXT_S,"HNEXT_S")); 
	// m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_HNEXT_XV_PGROW_BOUNDED,"HNEXT_XV_PGROW_BOUNDED")); 
	// m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_HNEXT_S_PGROW_BOUNDED,"HNEXT_S_PGROW_BOUNDED")); 
	// m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_HNEXT_XV_PGROW_UNBOUNDED,"HNEXT_XV_PGROW_UNBOUNDED")); 
	// m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_HNEXT_S_PGROW_UNBOUNDED,"HNEXT_S_PGROW_UNBOUNDED")); 
	// m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_HNEXT_XV_PSHRNK_BOUNDED,"HNEXT_XV_PSHRNK_BOUNDED")); 
	// m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_HNEXT_S_PSHRNK_BOUNDED,"HNEXT_S_PSHRNK_BOUNDED")); 
	// m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_HNEXT_XV_PSHRNK_UNBOUNDED,"HNEXT_XV_PSHRNK_UNBOUNDED")); 
	// m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_HNEXT_S_PSHRNK_UNBOUNDED,"HNEXT_S_PSHRNK_UNBOUNDED")); 
	// m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_HCURRENT,"HCURRENT")); 
 	// m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_EXTRA_0,"EXTRA_0")); 
	// m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_EXTRA_1,"EXTRA_1")); 
	// m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_EXTRA_2,"EXTRA_2")); 
	// m_vi_cstr_IO.insert(std::pair<int,string>(e_d_IO_EXTRA_3,"EXTRA_3")); 
	
	m_vi_cstr_IO_INT.insert(std::pair<int,string>(e_d_IO_INT_RETURN_VALUE_XV,"RETURN_VALUE_XV")); 
	m_vi_cstr_IO_INT.insert(std::pair<int,string>(e_d_IO_INT_RETURN_VALUE_S,"RETURN_VALUE_S")); 
	m_vi_cstr_IO_INT.insert(std::pair<int,string>(e_d_IO_INT_ERROR,"ERROR")); 
	m_vi_cstr_IO_INT.insert(std::pair<int,string>(e_d_IO_INT_THREAD,"THREAD")); 
	m_vi_cstr_IO_INT.insert(std::pair<int,string>(e_d_IO_INT_RECORD,"RECORD")); 
	m_vi_cstr_IO_INT.insert(std::pair<int,string>(e_d_IO_INT_RKQS_STEPS,"RKQS_STEPS"));
	m_vi_cstr_IO_INT.insert(std::pair<int,string>(e_d_IO_INT_RKQS_ERROR,"RKQS_ERROR"));
	
	m_vi_cstr_CONST_INT.insert(std::pair<int,string>(e_d_CONST_INT_numRecordsPerUpdate,"numRecordsPerUpdate"));
	m_vi_CONST_INT_DEFAULT.insert(std::pair<int,int>(e_d_CONST_INT_numRecordsPerUpdate, vi_INPUT_e_d_CONST_INT_numRecordsPerUpdate));//64[
	m_vi_cstr_CONST_INT.insert(std::pair<int,string>(e_d_CONST_INT_numCyclesPerRecord,"numCyclesPerRecord"));
	m_vi_CONST_INT_DEFAULT.insert(std::pair<int,int>(e_d_CONST_INT_numCyclesPerRecord, vi_INPUT_e_d_CONST_INT_numCyclesPerRecord)); //8096[
	m_vi_cstr_CONST_INT.insert(std::pair<int,string>(e_d_CONST_INT_def_CUDA_Mag,"def_CUDA_Mag"));
	m_vi_CONST_INT_DEFAULT.insert(std::pair<int,int>(e_d_CONST_INT_numCyclesPerRecord,1)); //8096[
	m_vi_cstr_CONST_INT.insert(std::pair<int,string>(e_d_CONST_INT_def_CUDA_dB,"def_CUDA_dB"));
	m_vi_CONST_INT_DEFAULT.insert(std::pair<int,int>(e_d_CONST_INT_def_CUDA_dB,5));
	m_vi_cstr_CONST_INT.insert(std::pair<int,string>(e_d_CONST_INT_FLAG_GRAVITY,"FLAG_GRAVITY"));
	m_vi_CONST_INT_DEFAULT.insert(std::pair<int,int>(e_d_CONST_INT_FLAG_GRAVITY,1));
	m_vi_cstr_CONST_INT.insert(std::pair<int,string>(e_d_CONST_INT_FLAG_SPRING,"FLAG_SPRING"));
	m_vi_CONST_INT_DEFAULT.insert(std::pair<int,int>(e_d_CONST_INT_FLAG_SPRING,1));
	m_vi_cstr_CONST_INT.insert(std::pair<int,string>(e_d_CONST_INT_FLAG_MAGNETIC,"FLAG_MAGNETIC"));
	m_vi_CONST_INT_DEFAULT.insert(std::pair<int,int>(e_d_CONST_INT_FLAG_MAGNETIC,2));
	m_vi_cstr_CONST_INT.insert(std::pair<int,string>(e_d_CONST_INT_FLAG_RF,"e_d_CONST_INT_FLAG_RF"));
	m_vi_CONST_INT_DEFAULT.insert(std::pair<int,int>(e_d_CONST_INT_FLAG_RF,0));
	m_vi_cstr_CONST_INT.insert(std::pair<int,string>(e_d_CONST_INT_def_CUDA_derivs_6,"def_CUDA_derivs_6"));
	m_vi_CONST_INT_DEFAULT.insert(std::pair<int,int>(e_d_CONST_INT_def_CUDA_derivs_6,3));
	m_vi_cstr_CONST_INT.insert(std::pair<int,string>(e_d_CONST_INT_def_PERFECT_POLARIZATION,"def_PERFECT_POLARIZATION"));
	m_vi_CONST_INT_DEFAULT.insert(std::pair<int,int>(e_d_CONST_INT_def_PERFECT_POLARIZATION,0));
	m_vi_cstr_CONST_INT.insert(std::pair<int,string>(e_d_CONST_INT_def_SPIN_INTEGRATOR_FLAG,"def_SPIN_INTEGRATOR_FLAG"));
	m_vi_CONST_INT_DEFAULT.insert(std::pair<int,int>(e_d_CONST_INT_def_SPIN_INTEGRATOR_FLAG,1));
	m_vi_cstr_CONST_INT.insert(std::pair<int,string>(e_d_CONST_INT_def_MAXSTP,"def_MAXSTP"));
	m_vi_CONST_INT_DEFAULT.insert(std::pair<int,int>(e_d_CONST_INT_def_MAXSTP,1000000));
	m_vi_cstr_CONST_INT.insert(std::pair<int,string>(e_d_CONST_INT_def_MAXSTP1,"def_MAXSTP1"));
	m_vi_CONST_INT_DEFAULT.insert(std::pair<int,int>(e_d_CONST_INT_def_MAXSTP1,1000000));
	m_vi_cstr_CONST_INT.insert(std::pair<int,string>(e_d_CONST_INT_def_FLAG_TEST_NEUTRON,"def_FLAG_TEST_NEUTRON"));
	m_vi_CONST_INT_DEFAULT.insert(std::pair<int,int>(e_d_CONST_INT_def_FLAG_TEST_NEUTRON,0));
	m_vi_cstr_CONST_INT.insert(std::pair<int,string>(e_d_CONST_INT_def_CONST_ACCEL,"def_CONST_ACCEL"));
	m_vi_CONST_INT_DEFAULT.insert(std::pair<int,int>(e_d_CONST_INT_def_CONST_ACCEL,0));
	m_vi_cstr_CONST_INT.insert(std::pair<int,string>(e_d_CONST_INT_def_NMAX,"def_NMAX"));
	m_vi_CONST_INT_DEFAULT.insert(std::pair<int,int>(e_d_CONST_INT_def_NMAX,0));
	m_vi_cstr_CONST_INT.insert(std::pair<int,string>(e_d_CONST_INT_def_IA,"def_IA"));
	m_vi_CONST_INT_DEFAULT.insert(std::pair<int,int>(e_d_CONST_INT_def_IA,16807));
	m_vi_cstr_CONST_INT.insert(std::pair<int,string>(e_d_CONST_INT_def_IM,"def_IM"));
	m_vi_CONST_INT_DEFAULT.insert(std::pair<int,int>(e_d_CONST_INT_def_IM,2147483647));
	m_vi_cstr_CONST_INT.insert(std::pair<int,string>(e_d_CONST_INT_def_IQ,"def_IQ"));
	m_vi_CONST_INT_DEFAULT.insert(std::pair<int,int>(e_d_CONST_INT_def_IQ,127773));
	m_vi_cstr_CONST_INT.insert(std::pair<int,string>(e_d_CONST_INT_def_IR,"def_IR"));
	m_vi_CONST_INT_DEFAULT.insert(std::pair<int,int>(e_d_CONST_INT_def_IR,2836));
	m_vi_cstr_CONST_INT.insert(std::pair<int,string>(e_d_CONST_INT_def_NTAB,"def_NTAB"));
	m_vi_CONST_INT_DEFAULT.insert(std::pair<int,int>(e_d_CONST_INT_def_NTAB,32));
	m_vi_cstr_CONST_INT.insert(std::pair<int,string>(e_d_CONST_INT_DIM_0_EXP,"DIM_0_EXP"));
	m_vi_CONST_INT_DEFAULT.insert(std::pair<int,int>(e_d_CONST_INT_DIM_0_EXP,0));
	m_vi_cstr_CONST_INT.insert(std::pair<int,string>(e_d_CONST_INT_DIM_1_EXP,"DIM_1_EXP"));
	m_vi_CONST_INT_DEFAULT.insert(std::pair<int,int>(e_d_CONST_INT_DIM_1_EXP,0));
	m_vi_cstr_CONST_INT.insert(std::pair<int,string>(e_d_CONST_INT_DIM_2_EXP,"DIM_2_EXP"));
	m_vi_CONST_INT_DEFAULT.insert(std::pair<int,int>(e_d_CONST_INT_DIM_2_EXP,6));
	m_vi_cstr_CONST_INT.insert(std::pair<int,string>(e_d_CONST_INT_DIM_BFIELD_FLAG,"DIM_BFIELD_FLAG"));
	m_vi_CONST_INT_DEFAULT.insert(std::pair<int,int>(e_d_CONST_INT_DIM_BFIELD_FLAG,0));
	m_vi_cstr_CONST_INT.insert(std::pair<int,string>(e_d_CONST_INT_DIM_0,"DIM_0"));
	m_vi_CONST_INT_DEFAULT.insert(std::pair<int,int>(e_d_CONST_INT_DIM_0,ipow(2,m_vi_CONST_INT[e_d_CONST_INT_DIM_0_EXP])));
	m_vi_cstr_CONST_INT.insert(std::pair<int,string>(e_d_CONST_INT_DIM_1,"DIM_1"));
	m_vi_CONST_INT_DEFAULT.insert(std::pair<int,int>(e_d_CONST_INT_DIM_1,ipow(2,m_vi_CONST_INT[e_d_CONST_INT_DIM_1_EXP])));
	m_vi_cstr_CONST_INT.insert(std::pair<int,string>(e_d_CONST_INT_DIM_2,"DIM_2"));
	m_vi_CONST_INT_DEFAULT.insert(std::pair<int,int>(e_d_CONST_INT_DIM_2,ipow(2,m_vi_CONST_INT[e_d_CONST_INT_DIM_2_EXP])));
	int vi_numBlocks_EXP_TEMP, vi_numThreadsPerBlock_EXP_TEMP;
	
	cout << "Threads = " << vi_INPUT_e_d_CONST_INT_numThreads_EXP << "\n";
	cout << "minThreadsPerBlock_EXP = " << vi_INPUT_e_d_CONST_INT_numThreadsPerBlock_EXP_MAX << "\n";
	if(vi_INPUT_e_d_CONST_INT_numThreads_EXP<vi_INPUT_e_d_CONST_INT_numThreadsPerBlock_EXP_MAX)
	{
		// vi_numBlocks_EXP_TEMP = 0;
		// vi_numThreadsPerBlock_EXP_TEMP = vi_INPUT_e_d_CONST_INT_numThreads_EXP;
		vi_numThreadsPerBlock_EXP_TEMP = vi_INPUT_e_d_CONST_INT_numThreads_EXP;
		vi_numBlocks_EXP_TEMP = 0;
	}
	else
	{
		// vi_numBlocks_EXP_TEMP = vi_INPUT_e_d_CONST_INT_numThreads_EXP - vi_INPUT_e_d_CONST_INT_numThreadsPerBlock_EXP_MAX;
		// vi_numThreadsPerBlock_EXP_TEMP = vi_INPUT_e_d_CONST_INT_numThreadsPerBlock_EXP_MAX;
		vi_numThreadsPerBlock_EXP_TEMP = vi_INPUT_e_d_CONST_INT_numThreadsPerBlock_EXP_MAX;
		vi_numBlocks_EXP_TEMP = (vi_INPUT_e_d_CONST_INT_numThreads_EXP - vi_INPUT_e_d_CONST_INT_numThreadsPerBlock_EXP_MAX);
	}
	cout << "Blocks_EXP = " << vi_numBlocks_EXP_TEMP << "\n";
	cout << "ThreadsPerBlock_EXP = " << vi_numThreadsPerBlock_EXP_TEMP << "\n";
	m_vi_cstr_CONST_INT.insert(std::pair<int,string>(e_d_CONST_INT_numThreadsPerBlock_EXP,"numThreadsPerBlock_EXP"));
	m_vi_CONST_INT_DEFAULT.insert(std::pair<int,int>(e_d_CONST_INT_numThreadsPerBlock_EXP, vi_numThreadsPerBlock_EXP_TEMP));
	m_vi_cstr_CONST_INT.insert(std::pair<int,string>(e_d_CONST_INT_numBlocks_EXP,"numBlocks_EXP"));
	m_vi_CONST_INT_DEFAULT.insert(std::pair<int,int>(e_d_CONST_INT_numBlocks_EXP,vi_numBlocks_EXP_TEMP));
	m_vi_cstr_CONST_INT.insert(std::pair<int,string>(e_d_CONST_INT_numRecordsPerThread_EXP,"numRecordsPerThread_EXP"));
	m_vi_CONST_INT_DEFAULT.insert(std::pair<int,int>(e_d_CONST_INT_numRecordsPerThread_EXP, vi_INPUT_e_d_CONST_INT_numRecordsPerThread_EXP));
	m_vi_cstr_CONST_INT.insert(std::pair<int,string>(e_d_CONST_INT_numRecords_EXP,"numRecords_EXP"));
	m_vi_CONST_INT_DEFAULT.insert(std::pair<int,int>(e_d_CONST_INT_numRecords_EXP, vi_INPUT_e_d_CONST_INT_numThreads_EXP + vi_INPUT_e_d_CONST_INT_numRecordsPerThread_EXP));
	m_vi_cstr_CONST_INT.insert(std::pair<int,string>(e_d_CONST_INT_numThreads_EXP,"numThreads_EXP"));
	m_vi_CONST_INT_DEFAULT.insert(std::pair<int,int>(e_d_CONST_INT_numThreads_EXP, vi_INPUT_e_d_CONST_INT_numThreads_EXP));
	
	cout << "\nATTTTTTTTEEENTION:\n\nthis value " << vi_INPUT_e_d_CONST_INT_numThreads_EXP << " should be the same as this value " <<  m_vi_CONST_INT_DEFAULT[e_d_CONST_INT_numThreads_EXP] << "\n\n";
	
	m_vi_cstr_CONST_INT.insert(std::pair<int,string>(e_d_CONST_INT_numThreadsPerBlock,"numThreadsPerBlock"));
	m_vi_CONST_INT_DEFAULT.insert(std::pair<int,int>(e_d_CONST_INT_numThreadsPerBlock,ipow(2,m_vi_CONST_INT_DEFAULT[e_d_CONST_INT_numThreadsPerBlock_EXP])));
	m_vi_cstr_CONST_INT.insert(std::pair<int,string>(e_d_CONST_INT_numRecordsPerThread,"numRecordsPerThread"));
	m_vi_cstr_CONST_INT.insert(std::pair<int,string>(e_d_CONST_INT_numBlocks,"numBlocks"));
	m_vi_CONST_INT_DEFAULT.insert(std::pair<int,int>(e_d_CONST_INT_numRecordsPerThread,ipow(2,m_vi_CONST_INT_DEFAULT[e_d_CONST_INT_numRecordsPerThread_EXP])));
	m_vi_CONST_INT_DEFAULT.insert(std::pair<int,int>(e_d_CONST_INT_numBlocks,ipow(2,m_vi_CONST_INT_DEFAULT[e_d_CONST_INT_numBlocks_EXP])));
	m_vi_cstr_CONST_INT.insert(std::pair<int,string>(e_d_CONST_INT_numRecords,"numRecords"));
	m_vi_CONST_INT_DEFAULT.insert(std::pair<int,int>(e_d_CONST_INT_numRecords,ipow(2,m_vi_CONST_INT_DEFAULT[e_d_CONST_INT_numRecords_EXP])));
	m_vi_cstr_CONST_INT.insert(std::pair<int,string>(e_d_CONST_INT_numThreads,"numThreads"));
	m_vi_CONST_INT_DEFAULT.insert(std::pair<int,int>(e_d_CONST_INT_numThreads,ipow(2,m_vi_CONST_INT_DEFAULT[e_d_CONST_INT_numThreads_EXP])));
	
	cout << "Threads: EXP=" << m_vi_CONST_INT_DEFAULT[e_d_CONST_INT_numThreads_EXP] << "\tACTUAL=" << m_vi_CONST_INT_DEFAULT[e_d_CONST_INT_numThreads] << "\n";	
	
	m_vi_cstr_CONST_INT.insert(std::pair<int,string>(e_d_CONST_INT_numRecordsPerReverse,"numRecordsPerReverse"));
	m_vi_CONST_INT_DEFAULT.insert(std::pair<int,int>(e_d_CONST_INT_numRecordsPerReverse,1));
	m_vi_cstr_CONST_INT.insert(std::pair<int,string>(e_d_CONST_INT_numRecordsPerReverse_SPIN,"numRecordsPerReverse_SPIN"));
	m_vi_CONST_INT_DEFAULT.insert(std::pair<int,int>(e_d_CONST_INT_numRecordsPerReverse_SPIN,0)); 
	m_vi_cstr_CONST_INT.insert(std::pair<int,string>(e_d_CONST_INT_FLAG_CLASSICAL_SPIN,"FLAG_CLASSICAL_SPIN"));
	m_vi_CONST_INT_DEFAULT.insert(std::pair<int,int>(e_d_CONST_INT_FLAG_CLASSICAL_SPIN,1)); 
	
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_tframe,"tframe"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_tframe,1e-4));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_tframe_SPIN,"tframe_SPIN"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_tframe_SPIN,1e-6)); // 1e-8));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_eps,"eps"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_eps,1e-10));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_eps_SPIN,"eps_SPIN"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_eps_SPIN,1e-12));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_h1,"h1"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_h1,1e-5));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_h1_SPIN,"h1_SPIN"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_h1_SPIN,1e-9));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_hmin,"hmin"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_hmin,1e-9));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_hmin_SPIN,"hmin_SPIN"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_hmin_SPIN,1e-12));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_SCALE,"def_SCALE"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_SCALE,1.0));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_SCALE_SPIN,"def_SCALE_SPIN"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_SCALE_SPIN,1.0));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_TIME_STEP_H,"def_TIME_STEP_H"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_TIME_STEP_H,1e-6));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_TIME_STEP_MIN,"def_TIME_STEP_MIN"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_TIME_STEP_MIN,1e-12));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_TIME_STEP_MAX,"def_TIME_STEP_MAX"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_TIME_STEP_MAX,1e-5));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_POS_X_TEST_NEUTRON,"def_POS_X_TEST_NEUTRON"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_POS_X_TEST_NEUTRON,0));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_POS_Y_TEST_NEUTRON,"def_POS_Y_TEST_NEUTRON"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_POS_Y_TEST_NEUTRON,0));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_POS_Z_TEST_NEUTRON,"def_POS_Z_TEST_NEUTRON"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_POS_Z_TEST_NEUTRON,0));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_VEL_X_TEST_NEUTRON,"def_VEL_X_TEST_NEUTRON"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_VEL_X_TEST_NEUTRON,0));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_VEL_Y_TEST_NEUTRON,"def_VEL_Y_TEST_NEUTRON"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_VEL_Y_TEST_NEUTRON,9.8));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_VEL_Z_TEST_NEUTRON,"def_VEL_Z_TEST_NEUTRON"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_VEL_Z_TEST_NEUTRON,0));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_ACCEL_X_CONST_ACCEL,"def_ACCEL_X_CONST_ACCEL"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_ACCEL_X_CONST_ACCEL,0));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_ACCEL_Y_CONST_ACCEL,"def_ACCEL_Y_CONST_ACCEL"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_ACCEL_Y_CONST_ACCEL,-9.8));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_ACCEL_Z_CONST_ACCEL,"def_ACCEL_Z_CONST_ACCEL"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_ACCEL_Z_CONST_ACCEL,0));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_SAFETY,"def_SAFETY"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_SAFETY,0.9));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_PGROW,"def_PGROW"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_PGROW,-0.2));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_PSHRNK,"def_PSHRNK"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_PSHRNK,-0.25));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_ERRCON,"def_ERRCON"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_ERRCON,1.89e-4));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_MAXERR,"def_MAXERR"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_MAXERR,1e-10));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_TINY,"def_TINY"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_TINY,1e-12)); 
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_SAFETY1,"def_SAFETY1"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_SAFETY1,0.9));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_PGROW1,"def_PGROW1"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_PGROW1,-0.2));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_PSHRNK1,"def_PSHRNK1"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_PSHRNK1,-0.25));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_ERRCON1,"def_ERRCON1"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_ERRCON1,1.89e-4)); 
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_TINY1,"def_TINY1"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_TINY1,1e-13));  
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_MAXERR1,"def_MAXERR1"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_MAXERR1,1e-12));	
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_TMAX,"def_TMAX"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_TMAX,60));	
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_MASS,"def_MASS"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_MASS,1.67492729e-27));	
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_MOMENT,"def_MOMENT"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_MOMENT,-9.662364e-27));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_MOMENT_DIV_MASS,"def_MOMENT_DIV_MASS"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_MOMENT_DIV_MASS,(1.67492729/9.662364)));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_PI,"def_PI"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_PI,3.1415926536));	
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_MU,"def_MU"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_MU,1.256637061e-7));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_HBAR,"def_HBAR"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_HBAR,1.0545717e-34));   
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_pos_BASE,"def_pos_BASE"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_pos_BASE,1e0));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_vel_BASE,"def_vel_BASE"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_vel_BASE,3e0));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_CORRECTIVE_FACTOR_SPIN,"def_CORRECTIVE_FACTOR_SPIN"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_CORRECTIVE_FACTOR_SPIN,1));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_START_SPIN_X,"def_START_SPIN_X"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_START_SPIN_X,0));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_START_SPIN_Y,"def_START_SPIN_Y"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_START_SPIN_Y,0));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_START_SPIN_Z,"def_START_SPIN_Z"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_START_SPIN_Z,0));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_HALBACH_NODES,"def_HALBACH_NODES"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_HALBACH_NODES,20));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_HALBACH_RADIUS,"def_HALBACH_RADIUS"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_HALBACH_RADIUS,2));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_HALBACH_CIRCUMFERENCE,"def_HALBACH_CIRCUMFERENCE"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_HALBACH_CIRCUMFERENCE,2*m_vi_CONST_DEFAULT[e_d_CONST_def_PI]*m_vi_CONST_DEFAULT[e_d_CONST_def_HALBACH_RADIUS]));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_HALBACH_WIDTH,"def_HALBACH_WIDTH"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_HALBACH_WIDTH,(m_vi_CONST_DEFAULT[e_d_CONST_def_HALBACH_CIRCUMFERENCE]/m_vi_CONST_DEFAULT[e_d_CONST_def_HALBACH_NODES])));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_HALBACH_K,"def_HALBACH_K"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_HALBACH_K,((2*m_vi_CONST_DEFAULT[e_d_CONST_def_PI])/m_vi_CONST_DEFAULT[e_d_CONST_def_HALBACH_WIDTH])));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_HALBACH_MAX_TESLA,"def_HALBACH_MAX_TESLA"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_HALBACH_MAX_TESLA,100));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_GRAVITY,"def_GRAVITY"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_GRAVITY,-9.8));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_SPRING,"def_SPRING"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_SPRING,1e4));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_DISTANCE_BETWEEN_COILS,"def_DISTANCE_BETWEEN_COILS"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_DISTANCE_BETWEEN_COILS,1.0));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_CURRENT_ANTI_HH,"def_CURRENT_ANTI_HH"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_CURRENT_ANTI_HH,1e13));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_AM,"def_AM"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_AM,(1.0/m_vi_CONST_DEFAULT[e_d_CONST_INT_def_IM])));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_NDIV,"def_NDIV"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_NDIV,(1+((m_vi_CONST_DEFAULT[e_d_CONST_INT_def_IM]-1)/m_vi_CONST_DEFAULT[e_d_CONST_INT_def_NTAB]))));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_EPS,"def_EPS"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_EPS,1.2e-7));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_def_RNMX,"def_RNMX"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_def_RNMX,(1.0-m_vi_CONST_DEFAULT[e_d_CONST_def_EPS])));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_guide_1_radius,"guide_1_radius"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_guide_1_radius,0.03));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_DEFAULT_X,"DEFAULT_X"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_DEFAULT_X,0.));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_DEFAULT_Y,"DEFAULT_Y"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_DEFAULT_Y,0.));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_DEFAULT_Z,"DEFAULT_Z"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_DEFAULT_Z,0.));
	m_vi_cstr_CONST.insert(std::pair<int,string>(e_d_CONST_RF_BFIELD_MAG,"RF_BFIELD_MAG"));
	m_vi_CONST_DEFAULT.insert(std::pair<int,double>(e_d_CONST_RF_BFIELD_MAG,2e-8));
	cout << "\n<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>";
	cout << "\nH1:" <<  m_vi_CONST_DEFAULT[e_d_CONST_h1];
	cout << "\n<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>";
	cout << "\nH1_SPIN:" <<  m_vi_CONST_DEFAULT[e_d_CONST_h1_SPIN];
	cout << "\n<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>";
	cout << "\nTHREADS_INPUT: N/A\tEXP:" <<  vi_INPUT_e_d_CONST_INT_numThreads_EXP;
	cout << "\n<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>";
	cout << "\nTHREADSPERBLOCK_MAX_INPUT: N/A\tEXP:" <<  vi_numThreadsPerBlock_EXP_TEMP;
	cout << "\n<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>";
	cout << "\nTHREADS:" << m_vi_CONST_INT_DEFAULT[e_d_CONST_INT_numThreads]  << "\tEXP:" <<  m_vi_CONST_INT_DEFAULT[e_d_CONST_INT_numThreads_EXP];
	cout << "\n<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>";
	cout << "\nTHREADSPERBLOCK:" << m_vi_CONST_INT_DEFAULT[e_d_CONST_INT_numThreadsPerBlock]  << "\tEXP:" <<  m_vi_CONST_INT_DEFAULT[e_d_CONST_INT_numThreadsPerBlock_EXP];
	cout << "\n<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>";
	cout << "\nBLOCKS:" << m_vi_CONST_INT_DEFAULT[e_d_CONST_INT_numBlocks]  << "\tEXP:" << m_vi_CONST_INT_DEFAULT[e_d_CONST_INT_numBlocks_EXP];
	cout << "\n<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>";
	
	int index = 0;
	cout << "==================" << param_m_vi_CONST_INT.size() << "\n";
	for(map<int,int>::iterator m_vi_vi_IT = m_vi_CONST_INT_DEFAULT.begin(); m_vi_vi_IT!=m_vi_CONST_INT_DEFAULT.end(); ++m_vi_vi_IT, index++)
  {
		cout << m_vi_cstr_CONST_INT[index] << " current value = " << m_vi_CONST_INT_DEFAULT[index] << " XXXXXXXXXX ";
		cout << m_vi_cstr_CONST_INT[index] << " next value = ";
		if(param_m_vi_CONST_INT.count(m_vi_vi_IT->first)==0)
		{
			m_vi_CONST_INT.insert ( pair<int,int>(m_vi_vi_IT->first, m_vi_vi_IT->second) );
			cout << (m_vi_vi_IT->second);
			// m_vi_CONST_INT.erase(m_vi_vi_IT->first);
			// m_vi_CONST_INT.insert ( pair<int,int>(m_vi_vi_IT->first, m_vi_vi_IT->second) );
		}
		else
		{
			map<int,int>::iterator m_vi_vi_IT_TEMP = param_m_vi_CONST_INT.find(m_vi_vi_IT->first);
			m_vi_CONST_INT.insert ( pair<int,int>((m_vi_vi_IT_TEMP->first), (m_vi_vi_IT_TEMP->second)) );
		}
	}
	index = 0;
	cout << "==================" << param_m_vi_CONST.size() << "\n";
	for(map<int,double>::iterator m_vi_vd_IT = m_vi_CONST_DEFAULT.begin(); m_vi_vd_IT!=m_vi_CONST_DEFAULT.end(); ++m_vi_vd_IT, index++)
  {
		cout << m_vi_cstr_CONST[index] << " current value = " << m_vi_CONST_DEFAULT[index] << " XXXXXXXXXX ";
		cout << m_vi_cstr_CONST[index] << " next value = ";
		
		if(param_m_vi_CONST.count(m_vi_vd_IT->first)==0)
		{
			m_vi_CONST.insert ( pair<int,double>( (m_vi_vd_IT->first), (m_vi_vd_IT->second) ) );
			cout << (m_vi_vd_IT->second);
			// m_vi_CONST.erase(m_vi_vd_IT->first);
			// m_vi_CONST.insert ( pair<int,double>(m_vi_vd_IT->first, m_vi_vd_IT->second) );
		}
		else
		{
			map<int,double>::iterator m_vi_vd_IT_TEMP = param_m_vi_CONST.find(m_vi_vd_IT->first);
			m_vi_CONST.insert ( pair<int,double>(m_vi_vd_IT_TEMP->first, m_vi_vd_IT_TEMP->second) );
			cout << ((param_m_vi_CONST.find(m_vi_vd_IT->first))->second);
		
		}
		cout << "\n";
	}
	// printf("e_d_IO_LAST=%d,",e_d_IO_LAST);
	// printf(",e_d_IO_INT_LAST=%d",e_d_IO_INT_LAST);
	// printf(",mapsize of m_vi_CONST=%d",m_vi_CONST.size());
	// printf(",mapsize of m_vi_CONST_INT=%d",m_vi_CONST_INT.size());
	// printf(",mapsize of m_vi_cstr_CONST=%d",m_vi_cstr_CONST.size());
	// printf(",mapsize of m_vi_cstr_CONST_INT=%d",m_vi_cstr_CONST_INT.size());
	// printf("e_d_CONST_LAST=%d",e_d_CONST_LAST);
	// printf(",e_d_CONST_INT_LAST=%d",e_d_CONST_INT_LAST);
	
	// for(map<int,int>::iterator mvi_IT = param_m_vi_CONST_INT.begin(); mvi_IT!=param_m_vi_CONST_INT.end(); ++mvi_IT)
  // {
		// cout << (m_vi_CONST_INT.find(mvi_IT->first))->second << "ehlloelhel\n";
		// if(m_vi_CONST_INT.count(mvi_IT->first)==1)
		// {
			// if(
			// m_vi_CONST_INT.erase(it->first);
			// m_vi_CONST_INT.insert ( pair<int,int>(it->first, it->second) );
		// }
	// }
	// for(map<int,double>::iterator it = param_m_vi_CONST.begin(); it!=param_m_vi_CONST.end(); ++it)
  // {
		// if(m_vi_CONST_INT.count(it->first)==1)
		// {
			// m_vi_CONST_INT.erase(it->first);
			// m_vi_CONST_INT.insert ( pair<int,double>(it->first, it->second) );
		// }
	// }
	
	cout << "\n<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>";
	cout << "\nH1:" <<  m_vi_CONST[e_d_CONST_h1];
	cout << "\n<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>";
	cout << "\nH1_SPIN:" <<  m_vi_CONST[e_d_CONST_h1_SPIN];
	cout << "\n<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>";
	
	for(int index = 0; index<e_d_CONST_LAST; index++) 
	{
		cout << "\n" << m_vi_cstr_CONST[index] << "=" << m_vi_CONST[index];
	}
	for(int index = 0; index<e_d_CONST_INT_LAST; index++) 
	{
		cout << "\n" << m_vi_cstr_CONST_INT[index] << "=" << m_vi_CONST_INT[index];
	}
	cout << "\n" << m_vi_cstr_CONST_INT[e_d_CONST_INT_numThreadsPerBlock] << "=" << m_vi_CONST_INT[e_d_CONST_INT_numThreadsPerBlock];
	cout << "\n" << m_vi_cstr_CONST_INT[e_d_CONST_INT_numThreadsPerBlock] << "=" << m_vi_CONST_INT[e_d_CONST_INT_numThreadsPerBlock];
	cout << "\n" << m_vi_cstr_CONST_INT[e_d_CONST_INT_numThreadsPerBlock] << "=" << m_vi_CONST_INT[e_d_CONST_INT_numThreadsPerBlock];
	cout << "\n" << m_vi_cstr_CONST_INT[e_d_CONST_INT_numThreadsPerBlock] << "=" << m_vi_CONST_INT[e_d_CONST_INT_numThreadsPerBlock];
	cout << "\n" << m_vi_cstr_CONST_INT[e_d_CONST_INT_numThreadsPerBlock] << "=" << m_vi_CONST_INT[e_d_CONST_INT_numThreadsPerBlock];
	int FIIK, i, j, k, l, m;
	int vi_PLACEHOLDER = 0;
	
	map<int, string> & r_m_vi_cstr_IO = m_vi_cstr_IO;
	map<int, string> & r_m_vi_cstr_IO_INT = m_vi_cstr_IO_INT;
	map<int, string> & r_m_vi_cstr_CONST = m_vi_cstr_CONST;
	map<int, string> & r_m_vi_cstr_CONST_INT = m_vi_cstr_CONST_INT;
	map<int, string> & r_m_vi_cstr_RKQS_ERROR = m_vi_cstr_RKQS_ERROR;
	map<int, int> & r_m_vi_CONST_INT = m_vi_CONST_INT;
	map<int, double> & r_m_vi_CONST = m_vi_CONST;
	
	const int numThreads_TEMP = r_m_vi_CONST_INT[e_d_CONST_INT_numThreads];
	const int numRecordsPerThread_TEMP = r_m_vi_CONST_INT[e_d_CONST_INT_numRecordsPerThread];
	const int numRecords_TEMP = numThreads_TEMP*numRecordsPerThread_TEMP;
	
	// int vi_a1_INPUT[numThreads_TEMP*e_d_IO_INT_LAST];
	// double vd_a1_INPUT[numThreads_TEMP*e_d_IO_LAST];
	
	// fxn_STANDARD_INPUT_DATA(
		// r_m_vi_cstr_IO,
		// r_m_vi_cstr_IO_INT,
		// r_m_vi_cstr_CONST,
		// r_m_vi_cstr_CONST_INT,
		// r_m_vi_cstr_RKQS_ERROR,
		// r_m_vi_CONST_INT,
		// r_m_vi_CONST,
		// numThreads_TEMP,
		// vi_a1_INPUT, 
		// vd_a1_INPUT);
		// cout << "\nJUST BEFORE CONSTRUCTOR";
	// for(int i = 0; i<numThreads_TEMP; i++)
	// {
		// int vi_THREAD_INT_OFFSET = i*m_vi_CONST_INT[e_d_CONST_INT_numRecordsPerThread]*e_d_IO_INT_LAST;
		// int vi_THREAD_DOUBLE_OFFSET = i*m_vi_CONST_INT[e_d_CONST_INT_numRecordsPerThread]*e_d_IO_LAST;		
		// printf("\nID:%d|RECORD:%d|TIME:%f|XYZ=(%f,%f,%f)|VEL=(%f,%f,%f)|SPIN:(%f,%f,%f,%f)|HNEXT:%f",
			// vi_a1_INPUT[vi_THREAD_INT_OFFSET + e_d_IO_INT_THREAD], 
			// vi_a1_INPUT[vi_THREAD_INT_OFFSET + e_d_IO_INT_RECORD], 
			// vd_a1_INPUT[vi_THREAD_DOUBLE_OFFSET + e_d_IO_T], 
			// vd_a1_INPUT[vi_THREAD_DOUBLE_OFFSET + e_d_IO_X], vd_a1_INPUT[vi_THREAD_DOUBLE_OFFSET + e_d_IO_Y], vd_a1_INPUT[vi_THREAD_DOUBLE_OFFSET + e_d_IO_Z], 
			// vd_a1_INPUT[vi_THREAD_DOUBLE_OFFSET + e_d_IO_VX], vd_a1_INPUT[vi_THREAD_DOUBLE_OFFSET + e_d_IO_VY], vd_a1_INPUT[vi_THREAD_DOUBLE_OFFSET + e_d_IO_VZ], 
			// vd_a1_INPUT[vi_THREAD_DOUBLE_OFFSET + e_d_IO_SPINNOR_0], vd_a1_INPUT[vi_THREAD_DOUBLE_OFFSET + e_d_IO_SPINNOR_1], vd_a1_INPUT[vi_THREAD_DOUBLE_OFFSET + e_d_IO_SPINNOR_2], vd_a1_INPUT[vi_THREAD_DOUBLE_OFFSET + e_d_IO_SPINNOR_3], 
			// vd_a1_INPUT[vi_THREAD_DOUBLE_OFFSET + e_d_IO_HNEXT]);
	// }
	return c_UCN_RUN_DATA(
		r_m_vi_cstr_IO,
		r_m_vi_cstr_IO_INT,
		r_m_vi_cstr_CONST,
		r_m_vi_cstr_CONST_INT,
		r_m_vi_cstr_RKQS_ERROR,
		r_m_vi_CONST_INT,
		r_m_vi_CONST,
		numThreads_TEMP,
		numRecordsPerThread_TEMP); //, 
		// vi_a1_INPUT, 
		// vd_a1_INPUT);
		
	// return fxn_PROXY_CONSTRUCTOR_c_UCN_RUN_DATA(
			// r_m_vi_cstr_IO,
			// r_m_vi_cstr_IO_INT,
			// r_m_vi_cstr_CONST,
			// r_m_vi_cstr_CONST_INT,
			// r_m_vi_cstr_RKQS_ERROR,
			// r_m_vi_CONST_INT,
			// r_m_vi_CONST);
}

