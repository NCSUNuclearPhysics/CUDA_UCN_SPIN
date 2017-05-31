
 
#include "UCN_ENUM.h"
/* User Variables ---------------------------------------------------------------- */

class c_UCN_RUN_DATA {
    const double *h_CONST; // [e_d_CONST_LAST];
    const int *h_CONST_INT; // [e_d_CONST_INT_LAST];
    const int vci_numThreads;
    const int vci_numRecordsPerThread;
    const int vci_numDoubles;
    const int vci_numInts;
    struct tm * s_tm_timeinfo_START;
    struct tm * s_tm_timeinfo_FINISH;
    vector<double> v_IO;
    vector<int> v_IO_INT;
  public:
    map<int, string> m_vi_cstr_IO;
    map<int, string> m_vi_cstr_IO_INT;
    map<int, string> m_vi_cstr_CONST;
    map<int, string> m_vi_cstr_CONST_INT;
    map<int, string> m_vi_cstr_RKQS_ERROR;
    map<int, double> m_vi_CONST;
    map<int, int> m_vi_CONST_INT;
    int fxn_FORMATTED_WRITE(
      const char *vc_a1_FILENAME, 
      bool vb_TIMESTAMP_FLAG);
    int fxn_FORMATTED_WRITE_ARRAY(
      double *vd_IO, 
      int *vi_IO, 
      const char *vc_a1_FILENAME, 
      bool vb_TIMESTAMP_FLAG);
    int fxn_FORMATTED_WRITE_MINIMUM(const char *vc_a1_FILENAME, bool vb_TIMESTAMP_FLAG);
    const int fvci_h_CONST(int vi_INDEX);
    const double fvcd_h_CONST(int vi_INDEX);
    int fxn_V3_OUTPUT_CSV();
    double fxn_vd_IO(int param_THREAD, int param_RECORD, int param_INDEX);
    int fxn_vi_IO(int param_THREAD, int param_RECORD, int param_INDEX);
    // c_UCN_RUN_DATA(
      // map<int, string> & r_m_vi_cstr_IO_INPUT,
      // map<int, string> & r_m_vi_cstr_IO_INT_INPUT,
      // map<int, string> & r_m_vi_cstr_CONST_INPUT,
      // map<int, string> & r_m_vi_cstr_CONST_INT_INPUT,
      // map<int, string> & r_m_vi_cstr_RKQS_ERROR,
      // map<int, int> & r_m_vi_CONST_INT_INPUT,
      // map<int, double> & r_m_vi_CONST_INPUT,
      // const int param_numThreads,
      // const int param_numRecordsPerThread);
    c_UCN_RUN_DATA(
      map<int, string> & r_m_vi_cstr_IO_INPUT,
      map<int, string> & r_m_vi_cstr_IO_INT_INPUT,
      map<int, string> & r_m_vi_cstr_CONST_INPUT,
      map<int, string> & r_m_vi_cstr_CONST_INT_INPUT,
      map<int, string> & r_m_vi_cstr_RKQS_ERROR,
      map<int, int> & r_m_vi_CONST_INT_INPUT,
      map<int, double> & r_m_vi_CONST_INPUT,
      const int param_numThreads,
      const int param_numRecordsPerThread);
};
// c_UCN_RUN_DATA::vsi_obj_ID = 0;
// c_UCN_RUN_DATA fxn_PROXY_CONSTRUCTOR_c_UCN_RUN_DATA(
  // map<int, string> & r_m_vi_cstr_IO,
  // map<int, string> & r_m_vi_cstr_IO_INT,
  // map<int, string> & r_m_vi_cstr_CONST,
  // map<int, string> & r_m_vi_cstr_CONST_INT,
  // map<int, string> & r_m_vi_cstr_RKQS_ERROR,
  // map<int, int> & r_m_vi_CONST_INT, 
  // map<int, double> & r_m_vi_CONST);
int fxn_STANDARD_INPUT_DATA(
  map<int, string> & r_m_vi_cstr_IO_INPUT, 
  map<int, string> & r_m_vi_cstr_IO_INT_INPUT, 
  map<int, string> & r_m_vi_cstr_CONST_INPUT, 
  map<int, string> & r_m_vi_cstr_CONST_INT_INPUT, 
  map<int, string> & r_m_vi_cstr_RKQS_ERROR, 
  map<int, int> & r_m_vi_CONST_INT_INPUT, 
  map<int, double> & r_m_vi_CONST_INPUT, 
  // const int param_numThreads, 
  // const int param_numRecordsPerThread, 
  int *param_a1_INPUT_INT, 
  double *param_a1_INPUT);
c_UCN_RUN_DATA fxn_PARAMETERIZED_CONSTRUCTOR_c_UCN_RUN_DATA(
  const int vi_INPUT_e_d_CONST_INT_numThreadsPerBlock_EXP_MAX, 
  const int vi_INPUT_e_d_CONST_INT_numRecordsPerThread_EXP, 
  const int vi_INPUT_e_d_CONST_INT_numThreads_EXP, 
  const int vi_INPUT_e_d_CONST_INT_numRecordsPerUpdate, 
  const int vi_INPUT_e_d_CONST_INT_numCyclesPerRecord,
  map<int, int> m_vi_CONST_INT,
  map<int, double> m_vi_CONST);