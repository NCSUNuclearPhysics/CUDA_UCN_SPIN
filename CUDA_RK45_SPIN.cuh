#include "UCN_ENUM.h"

/* define a debugging variable which will turn on and off printf's
  in each subroutine.  Different levels, 0 is nothing, 1 is simple
  messages, 2 is info about variables, .... */
// __constant__ const int def_DEBUGGING = 0;
// __constant__ const int int_DEBUG_SDM = 0;

// FUNCTION DECLARATIONS

__device__ /*__host__*/ int CUDA_ipow( int base, int exp);
__device__ /*__host__*/ int getGlobal_blockId_3D();
__device__ /*__host__*/ int getGlobalIdx_3D_3D();
__device__ /*__host__*/ int fvi_ISOLATE_INT_RANGE( int vi_INPUT, int vi_MSB, int  vi_LSB);
__device__ /*__host__*/ double CUDA_phi( double x, double y);
__device__ /*__host__*/ double CUDA_polcalc_XVS( double bb1x, double bb1y, double bb1z, double spinor[]);
__device__ /*__host__*/ double CUDA_theta( double x, double y, double z);
__device__ /*__host__*/ double CUDA_ro( double x, double y);
__device__ /*__host__*/ double CUDA_r( double x, double y, double z);
__device__ /*__host__*/ int CUDA_OUTPUT_OFFSET( int vi_RECORD);
__device__ /*__host__*/ int CUDA_OUTPUT_OFFSET_INT( int vi_RECORD);
__device__ /*__host__*/ int CUDA_INPUT_OFFSET();
__device__ /*__host__*/ int CUDA_INPUT_OFFSET_INT();
__device__ /*__host__*/ int CUDA_derivs_XVS( double t, double xvs[], double dxvsdt[], double BField[], double t_START, double inter[]);
__device__ /*__host__*/ int CUDA_rkck_XVS( double *d_IO, int *d_IO_INT, int vi_RECORD, double t, double xvs[], double dxvsdt[], double h, double xvserr[], double xvsout[], double BField[], double inter[]); 
__device__ /*__host__*/ int CUDA_rkqs_XVS( double *d_IO, int *d_IO_INT, int vi_RECORD, double xvs[], double dxvsdt[], double *t, double htry, double *hdid, double *hnext, double xvs_scal[], double BField[], double epsilon[], int *rkqs_steps, double inter[]); //);
__device__ /*__host__*/ int CUDA_RECORD_INT( int *d_IO_INT, int vi_RECORD, int e_d_IO_INT_PARAM, int vi_PARAM);
__device__ /*__host__*/ int CUDA_RECORD_DOUBLE( double *d_IO,  int vi_RECORD, int e_d_IO_PARAM, double vd_PARAM);
__device__ /*__host__*/ int CUDA_RECORD_XVS( double *d_IO, int *d_IO_INT, int *p_vi_RECORD, 
  double l_time_CURRENT, double l_xvs[], double l_epsilon[], double l_xvs_scal[], double l_dxvsdt[], 
  double l_BField[], double l_pol, int l_rkqs_TRIED);
__device__ /*__host__*/ int CUDA_rkqs_SINGLE_ATTEMPT_XVS(
  double *d_IO, int *d_IO_INT, int vi_RECORD, 
  double xvs[],   double dxvsdt[],   double *t,   double htry, 
  double *hdid,   double *hnext,   double xvs_scal[], double BField[], 
  double epsilon[], int *rkqs_TRIED);
__global__ void GENERIC_PIECEWISE_KERNEL_MULTI_XVS_RKQS_LOOP(
  double *d_IO, 
  int *d_IO_INT, 
  int numRecordsStart,
  int numRecordsEnd);
void GENERIC_RECORD_FRAME(
  int param_numBlocks,
  int param_numThreadsPerBlock,
  double *d_IO,
  int *d_IO_INT,
  int vi_RecordsStartCurrent,
  int vi_RecordsEndCurrent);
void GENERIC_MIDDLEMAN_MULTI(
  const double *h_CONST,
  const int *h_CONST_INT,
  double *h_IO,
  int *h_IO_INT);
