
/* Include files ----------------------------------------------------------------- */
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <cuda.h>
#include <cuda_runtime.h>
// #include "cuPrintf.cuh"
// #include "cuPrintf.cu"

/* 

    ###########################################################
    # This is version 2.0
    # CUDA Runge-Kutta 4th/5th Order
    #
    # Author : Spencer McBride Day Moore
    # Date   : May 15, 2015
    ###########################################################

    This is a standalone nvcc compiled program that implements a function in order
    to be implemented with double precision floating point numbers compute capability
    of CUDA device must be 2.0 or higher (-arch=sm_20 is the flag that goes with nvcc 
    for 2.0, -arch=sm_21 for 2.1 and so on). the #define variables are listed below 
    the header files.  Also the function fxn_dydx is the example given is only allowed
    single variable dependence in addition to time f_x(t,x)=dx/dt, f_y(t,y)=dy/dt etc. 
    and not f_x(t,x,y,z)=dx/dt etc..  This will need to be changed for usefulness in 
    the neutron transport code.  Currently the user can pass variables to the 
    individual threads via the constant double arrays cu_A and cu_B.  They are currently
    initialized ot random values between positive and negative one for testing purposes. 
    There is no adaptive step implemented although the stepwise error is outputted along
    with (t,x,y,z).
*/

#include "UCN_CUDA_ALL_KERNEL.cuh"

__constant__ int d_CONST_INT[e_d_CONST_INT_LAST];
__constant__ double d_CONST[e_d_CONST_LAST];

__device__ int CUDA_ipow( int base, int exp)
{
  int result = -1;
  if (exp<0) result = 0;
  else if (exp==0) result = 1;
  else if (exp>1)
  {
    result = 1;
    for (int i = exp; i>0; i--)
    {
        result *= base;
    }
  }
  return result;
}
__device__ int getGlobal_blockId_3D()
{
  int blockId = blockIdx.x + 
    (blockIdx.y * gridDim.x)  + (
    gridDim.x * gridDim.y * blockIdx.z);
  return blockId;
}
__device__ int getGlobalIdx_3D_3D()
{
  int blockId = getGlobal_blockId_3D();
  int threadId = 
    (blockId * blockDim.x * blockDim.y * blockDim.z) + 
    threadIdx.z * (blockDim.x * blockDim.y) + 
    threadIdx.y * blockDim.x + 
    threadIdx.x;
  return threadId;
  // int vi_RETURN_VALUE = ( blockIdx.x * blockDim.x + threadIdx.x );
  // return vi_RETURN_VALUE
}
__device__ int fvi_ISOLATE_INT_RANGE( int vi_INPUT, int vi_MSB, int vi_LSB)
{
  int vi_ABS_INPUT = vi_INPUT < 0 ? -vi_INPUT : vi_INPUT;
  int vi_INDEX, vi_MOD;
  for (vi_INDEX = 1, vi_MOD = 1 ; vi_INDEX<vi_MSB; vi_INDEX++)
  {
    vi_MOD = vi_MOD * 2;
  }
  int vi_INTERMEDIATE = ( vi_ABS_INPUT % vi_MOD );
  int vi_FINAL = ( vi_INTERMEDIATE >> ( vi_LSB - 1 ) );
  if ( vi_FINAL != ( vi_INTERMEDIATE * CUDA_ipow(2, vi_LSB - 1))) vi_FINAL = -1;
  return vi_FINAL;
}
__device__ /*__host__*/ double CUDA_theta( double x, double y, double z)
{
  double temp;
  if ( (z!=0) && ((x*x+y*y)!=0) ) temp = atan2(sqrt(x*x+y*y),z);
  else if (z==0) temp = d_CONST[e_d_CONST_def_PI]/2.0;
  else temp = d_CONST[e_d_CONST_def_TINY];  
  if (temp >= d_CONST[e_d_CONST_def_PI]) temp=2*d_CONST[e_d_CONST_def_PI]-temp;  
  return(temp);
}
__device__ /*__host__*/ double CUDA_ro( double x, double y)
{
  double solution;
  if (x*x+y*y != 0) solution = sqrt(x*x + y*y);
  else solution = ( d_CONST[e_d_CONST_def_TINY] );
  return solution;
}
__device__ /*__host__*/ double CUDA_r( double x, double y, double z)
{
  double solution;
  if ((x*x+y*y+z*z)!= 0) solution = (sqrt(x*x+y*y+z*z));
  else solution = ( d_CONST[e_d_CONST_def_TINY] );
  return solution;
}
__device__ /*__host__*/ double CUDA_phi( double x, double y)
{
  double solution; 
  if (x!=0)
    solution = atan2(y,x);
  else
    solution = ((d_CONST[e_d_CONST_def_PI]/2.0));
  return solution;
}
__device__ /*__host__*/ double CUDA_polcalc_XVS( 
  double bb1x, 
  double bb1y, 
  double bb1z, 
  double spinor[])
{
  double xhat1,yhat1,zhat1,snorm;
  /* const double xhat2,yhat2,zhat2; */
  double bxhat,byhat,bzhat;
  double bnorm,polarization;
  /* const double spinnorm,rnorm,rtest; */

  bnorm = sqrt(bb1x*bb1x + bb1y*bb1y + bb1z*bb1z);
  bxhat = bb1x/bnorm;
  byhat = bb1y/bnorm;
  bzhat = bb1z/bnorm;
  
  xhat1 = 0.;
  yhat1 = 0.;
  zhat1 = 0.;
  // CUDA_SPIN_XVS_XYZ( spin, spinor );
  snorm = (spinor[0]*spinor[0]) + (spinor[1]*spinor[1]) +  
          (spinor[2]*spinor[2]) + (spinor[3]*spinor[3]);
  xhat1 =( 2.*((spinor[0]*spinor[2]) + (spinor[1]*spinor[3])))/snorm;
  yhat1 = ( 2.*((spinor[0]*spinor[3]) - (spinor[1]*spinor[2])))/snorm;
  zhat1 = ((spinor[0]*spinor[0]) + (spinor[1]*spinor[1]) - 
     (spinor[2]*spinor[2]) - (spinor[3]*spinor[3]))/snorm;
  polarization = xhat1*bxhat + yhat1*byhat + zhat1*bzhat;
  return(polarization);
}
__device__ /*__host__*/ int CUDA_derivs_XVS( 
  double t,
  double spinor[],
  double dspinordt[],
  double BField[])
{  
  double rfstr = d_CONST[e_d_CONST_RF_BFIELD_MAG];
  double  brf[4],omega,btoomega,zrf,zdist,zdelta;
  btoomega = d_CONST[e_d_CONST_def_MOMENT]/d_CONST[e_d_CONST_def_HBAR];  // coeff. for psi-dot 
  omega    = 2.0*d_CONST[e_d_CONST_def_MOMENT]*(1.0)/d_CONST[e_d_CONST_def_HBAR]; // reson. frequ. for 1.0 T field 
  zrf      = 1.14; // chosen to match crossing poconst int 
  // some comment explaining this. v 
	double t_INTERVAL = inter[e_inter_t_FINAL] - inter[e_inter_t_INITIAL];
	double t_STEP = t - inter[e_inter_t_INITIAL];
	double BField_STEP[3];
	spin_INTERVAL_3 = inter[e_inter_spin_FINAL_3] - inter[e_inter_spin_INITIAL_3];
  double BField_INTERVAL_X = inter[e_inter_BField_FINAL_X] - inter[e_inter_BField_INITIAL_X];
  double BField_INTERVAL_Y = inter[e_inter_BField_FINAL_Y] - inter[e_inter_BField_INITIAL_Y];
  double BField_INTERVAL_Z = inter[e_inter_BField_FINAL_Z] - inter[e_inter_BField_INITIAL_Z];
  BField[0] = (BField_INTERVAL_X/t_INTERVAL)*t_STEP + inter[e_inter_BField_INITIAL_X];
  BField[1] = (BField_INTERVAL_Y/t_INTERVAL)*t_STEP + inter[e_inter_BField_INITIAL_Y];
  BField[2] = (BField_INTERVAL_Z/t_INTERVAL)*t_STEP + inter[e_inter_BField_INITIAL_Z];
  
	zdist = (spin_INTERVAL_3/t_INTERVAL)*t_STEP+inter[e_inter_spin_INITIAL_3]; //spin[2];
  zdelta = (zdist - zrf)*(zdist -zrf)/(.05*.05);
  brf[0] = rfstr*cos(omega*(t-0.0)) * (1.0/((1.0+zdelta)*(1.0+zdelta))); 
  brf[1] = 0.;
    // exp(-(zdist-zrf)*(zdist-zrf)/(.05*.05)); 
  brf[2] = 0.;
  int i;
  double extra_brf[3];
  // CUDA_RF_BFIELD(t, spinor, extra_brf);
  if (d_CONST_INT[e_d_CONST_INT_FLAG_RF]==1)
  {
    for ( i = 0; i<3; i++) BField[i] += brf[i]; // extra_brf[i];
  }
  
  // Classical Description of Spin 
  // NB moment has "I=1/2" already... 
  // dspinordt[0] = 0.0;
  // dspinordt[1] = btoomega*(spinor[2]*BField[2] - spinor[3]*BField[1]);
  // dspinordt[2] = btoomega*(spinor[3]*BField[0] - spinor[1]*BField[2]);
  // dspinordt[3] = btoomega*(spinor[1]*BField[1] - spinor[2]*BField[0]);
  // ihd/dt = (moment)*(sigma)*B 
  // NB moment has "I=1/2" already... 
  dspinordt[0] = btoomega * ( (BField[2] * spinor[1]) + (BField[0] * spinor[3]) - (BField[1] * spinor[2]));
  dspinordt[1] = btoomega * (-(BField[2] * spinor[0]) - (BField[0] * spinor[2]) - (BField[1] * spinor[3]));
  dspinordt[2] = btoomega * (-(BField[2] * spinor[3]) + (BField[0] * spinor[1]) + (BField[1] * spinor[0]));
  dspinordt[3] = btoomega * ( (BField[2] * spinor[2]) - (BField[0] * spinor[0]) + (BField[1] * spinor[1]));
  return 0;
}
__device__ /*__host__*/ int CUDA_rkck_XVS( 
  double *d_IO, 
  int *d_IO_INT, 
  int vi_RECORD,
  double t,
  double spinor[],
  double dspinordt[],
  double h,
  double spinorerr[],
  double spinorout[],
  double BField[])
{
  // double BField[3][3];
  CUDA_derivs_XVS(t, spinor, dspinordt, BField);
  int i;
  /*   static const a2=(0.2, a3=(0.3, a4=(0.6, a5=(1.0, a6=(0.875; */
  const double b21=(0.2), b31=(3.0/40.0), b32=(9.0/40.0), b41=(0.3), b42=( -0.9), b43=(1.2);
  const double b51=( -11.0/54.0), b52=(2.5), b53=( -70.0/27.0), b54=(35.0/27.0); 
  const double b61=(1631.0/55296.0), b62=(175.0/512.0), b63=(575.0/13824.0), b64=(44275.0/110592.0), b65=(253.0/4096.0); 
  const double c1=(37.0/378.0), c3=(250.0/621.0), c4=(125.0/594.0), c6=(512.0/1771.0);
  const double dc5=( -277.0/14336.0), dc1=(c1-2825.0/27648.0),  dc3=(c3-18575.0/48384.0), dc4=(c4-13525.0/55296.0), dc6=(c6-0.25);
  double ak2[4], ak3[4], ak4[4], ak5[4], ak6[4], spinortemp[4];
  for (i = 0; i<4; i++)   /* First step */ spinortemp[i]=spinor[i]+b21*h*dspinordt[i];
  CUDA_derivs_XVS(t, spinor, ak2, BField);    /* Second step */
  for (i = 0; i<4; i++) spinortemp[i]=spinor[i]+h*(b31*dspinordt[i]+b32*ak2[i]);
  CUDA_derivs_XVS(t, spinor, ak3, BField);    /* Third step */
  for (i = 0; i<4; i++) spinortemp[i]=spinor[i]+h*(b41*dspinordt[i]+b42*ak2[i]+b43*ak3[i]);
  CUDA_derivs_XVS(t, spinor, ak4, BField);    /* Fourth step */
  for (i = 0; i<4; i++) spinortemp[i]=spinor[i]+h*(b51*dspinordt[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
  CUDA_derivs_XVS(t, spinortemp, ak5, BField);    /* Fifth step */
  for (i = 0; i<4; i++) spinortemp[i]=spinor[i]+h*(b61*dspinordt[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
  CUDA_derivs_XVS(t, spinor, ak6, BField);    /* Sixth step */
  for (i = 0; i<4; i++)  spinorout[i]=spinor[i]+h*(c1*dspinordt[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]); /* Accumulate increments with proper weights */
  for (i = 0; i<4; i++) spinorerr[i]=h*(dc1*dspinordt[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
  return 0;
}
__device__ /*__host__*/ int CUDA_rkqs_SINGLE_ATTEMPT_XVS(
  double *d_IO, 
  int *d_IO_INT,
  int vi_RECORD, 
  double spinor[],
   double dspinordt[],
   double *t,
   double htry,
  double *hdid,
   double *hnext,
   double spinor_scal[],
   double BField[],
   double epsilon[],
   int *rkqs_TRIED,
	 double inter[])
{
  int i, j, k;
  double hnext_xv;
  double hnext_s;
  double hcurrent = htry;
  double spinor_temp[4];
  double epsilon_temp[4];
  double error, error_xv, error_s;
  int return_value = -1;
  int return_value_xv = -1;
  int return_value_s = -1;
  double abs_max_xv, abs_max_s;
  double hnext_s_GB, hnext_s_GU, hnext_s_SB, hnext_s_SU;
  double hnext_xv_GB, hnext_xv_GU, hnext_xv_SB, hnext_xv_SU;
  CUDA_rkck_XVS( 
    d_IO, 
    d_IO_INT, 
    vi_RECORD, 
    *t, 
    spinor, 
    dspinordt, 
    hcurrent, 
    epsilon_temp, 
    spinor_temp, 
    BField);
  double epsilon_max_xv = 0.0;
  double epsilon_max_s = 0.0;
  for (i = 0; i<4; i++)
  {
		abs_max_s = fabs(epsilon_temp[i]/spinor_scal[i]);
		// if (abs_max_s<0) abs_max_s *= -1.0; 
		if (abs_max_s>epsilon_max_s) epsilon_max_s = abs_max_s;
  }
  error_s = epsilon_max_s/d_CONST[e_d_CONST_def_MAXERR1];
  
  /////////////////////////////////////////////////////////////////////////
  /////////////////// SPIN STEP COMPUTATION ////////////////////////
  /////////////////////////////////////////////////////////////////////////
  
  /* S: Step succeeded. Compute size of next step */
  hnext_s_GU = d_CONST[e_d_CONST_def_SAFETY1]*hcurrent*pow(error_s,d_CONST[e_d_CONST_def_PGROW1]);
  /* S: No more than a factor of 5 increase */
  hnext_s_GB = 5.0*hcurrent; // ((hnext_s_GU >= 0.0) ? fmin(hnext_s_GU, 5.0*hcurrent) : fmax(hnext_s_GU, 5.0*hcurrent));
  /* S: Truncation error too large. Reduce stepsize */
  hnext_s_SU = d_CONST[e_d_CONST_def_SAFETY1]*hcurrent*pow(error_s, d_CONST[e_d_CONST_def_PSHRNK1]);
  /* S: No more than a factor of 10 */
  hnext_s_SB = 0.1*hcurrent; // ((hcurrent >= 0.0) ? fmax(hnext_s_SU, 0.1*hcurrent) : fmin(hnext_s_SU, 0.1*hcurrent));
  
 /////////////////////////////////////////////////////////////////////////
  ////////////////// SPIN ERROR AND STEP DECISION /////////////////////////
  /////////////////////////////////////////////////////////////////////////
  
  if (error_s <= 1.0)
  {
    hnext_s = ((hcurrent >= 0.0) ? fmin(hnext_s_GU, hnext_s_GB) : fmax(hnext_s_GU, hnext_s_GB));
    return_value_s = 0;
  }
  else if (error_s > 1.0)
  {
    hnext_s = ((hcurrent >= 0.0) ? fmax(hnext_s_SU, hnext_s_SB) : fmin(hnext_s_SU, hnext_s_SB));
    return_value_s = 1;
  }
  
  /////////////////////////////////////////////////////////////////////////
  /////////////////// FINAL TIME STEP SIGN CHECK //////////////////////////
  /////////////////////////////////////////////////////////////////////////
  
  
  
  /////////////////////////////////////////////////////////////////////////
  /////////////////// FINAL TIME STEP DECISION ////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  
  /* nrerror("stepsize underflow in rkqs"); */ /* diag */
  if ((*t + hcurrent) == *t) return_value = e_RKQS_ERROR_STEPSIZE_UNDERFLOW;
  else if (return_value_xv==0 && return_value_s==0)
  {
    for (i=0; i<4; i++)
    {
      spinor[i] = spinor_temp[i];
      epsilon[i] = epsilon_temp[i];
    }
    *t += hcurrent;
    *hdid = hcurrent;
    *rkqs_TRIED = 0;
    if (hcurrent<=0.0)
    {
      *hnext = fmax ( hnext_s , fmax ( d_CONST[e_d_CONST_tframe], d_CONST[e_d_CONST_tframe_SPIN]));
    }
    else if (hcurrent>0.0) 
    {
      *hnext = fmin ( hnext_s , fmin ( d_CONST[e_d_CONST_tframe], d_CONST[e_d_CONST_tframe_SPIN]));
    }
    return_value = e_RKQS_ERROR_NONE;
  }
  else
  {
    *hdid = 0.0;
    (*rkqs_TRIED)++;
    if (return_value_s==1)
    {
      *hnext = hnext_s;
      return_value = e_RKQS_ERROR_SPIN_BOUNDS;
    }
    else return_value = e_RKQS_ERROR_UNKNOWN;
  }
  
  /////////////////////////////////////////////////////////////////////////
  /////////////////// FINAL TIME STEP DECISION ////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  
  // CUDA_RECORD_DOUBLE(d_IO, vi_RECORD, e_d_IO_HCURRENT, hcurrent);
  // CUDA_RECORD_DOUBLE(d_IO, vi_RECORD, e_d_IO_EPSILON_MAX_XV, epsilon_max_xv);
  // CUDA_RECORD_DOUBLE(d_IO, vi_RECORD, e_d_IO_EPSILON_MAX_S, epsilon_max_s);
  // CUDA_RECORD_DOUBLE(d_IO, vi_RECORD, e_d_IO_ERROR_XV, error_xv);
  // CUDA_RECORD_DOUBLE(d_IO, vi_RECORD, e_d_IO_ERROR_S, error_s);
  // CUDA_RECORD_DOUBLE(d_IO, vi_RECORD, e_d_IO_HNEXT, *hnext);
  // CUDA_RECORD_DOUBLE(d_IO, vi_RECORD, e_d_IO_HNEXT_XV, hnext_xv);
  // CUDA_RECORD_DOUBLE(d_IO, vi_RECORD, e_d_IO_HNEXT_S, hnext_s);
  // CUDA_RECORD_DOUBLE(d_IO, vi_RECORD, e_d_IO_HNEXT_XV_PGROW_BOUNDED, hnext_xv_GB);
  // CUDA_RECORD_DOUBLE(d_IO, vi_RECORD, e_d_IO_HNEXT_S_PGROW_BOUNDED, hnext_s_GB);
  // CUDA_RECORD_DOUBLE(d_IO, vi_RECORD, e_d_IO_HNEXT_XV_PGROW_UNBOUNDED, hnext_xv_GU);
  // CUDA_RECORD_DOUBLE(d_IO, vi_RECORD, e_d_IO_HNEXT_S_PGROW_UNBOUNDED, hnext_s_GU);
  // CUDA_RECORD_DOUBLE(d_IO, vi_RECORD, e_d_IO_HNEXT_XV_PSHRNK_BOUNDED, hnext_xv_SB);
  // CUDA_RECORD_DOUBLE(d_IO, vi_RECORD, e_d_IO_HNEXT_S_PSHRNK_BOUNDED, hnext_s_SB);
  // CUDA_RECORD_DOUBLE(d_IO, vi_RECORD, e_d_IO_HNEXT_XV_PSHRNK_UNBOUNDED, hnext_xv_SU);
  // CUDA_RECORD_DOUBLE(d_IO, vi_RECORD, e_d_IO_HNEXT_S_PSHRNK_UNBOUNDED, hnext_s_SU);
  // CUDA_RECORD_INT(d_IO_INT, vi_RECORD, e_d_IO_INT_RETURN_VALUE_XV, return_value_xv);
  // CUDA_RECORD_INT(d_IO_INT, vi_RECORD, e_d_IO_INT_RETURN_VALUE_S, return_value_s);
  // CUDA_RECORD_INT(d_IO_INT, vi_RECORD, e_d_IO_INT_RKQS_ERROR, return_value);
  return return_value;
}
__device__ /*__host__*/ int CUDA_RECORD_INT( 
  int *d_IO_INT, int vi_RECORD, int e_d_IO_INT_PARAM, int vi_PARAM)
{
  int vi_INT_IO_OFFSET = (getGlobalIdx_3D_3D()*d_CONST_INT[e_d_CONST_INT_numRecordsPerThread] + vi_RECORD)*e_d_IO_INT_LAST;
  d_IO_INT[vi_INT_IO_OFFSET + e_d_IO_INT_PARAM] = vi_PARAM;
  return 0;
}
__device__ /*__host__*/ int CUDA_RECORD_DOUBLE( 
  double *d_IO,  int vi_RECORD, int e_d_IO_PARAM, double vd_PARAM)
{
  int vi_DOUBLE_IO_OFFSET = (getGlobalIdx_3D_3D()*d_CONST_INT[e_d_CONST_INT_numRecordsPerThread] + vi_RECORD)*e_d_IO_LAST;
  d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_PARAM] = vd_PARAM;
  return 0;
}
__device__ /*__host__*/ int CUDA_RECORD_XVS( 
  double *d_IO, int *d_IO_INT,  int *p_vi_RECORD, 
  double l_time_CURRENT, double l_spinor[], double l_epsilon[], double l_spinor_scal[], double l_dspinordt[], 
  double l_BField[], double l_pol, int l_rkqs_TRIED, double l_hnext)
{
  int vi_DOUBLE_IO_OFFSET = (getGlobalIdx_3D_3D()*d_CONST_INT[e_d_CONST_INT_numRecordsPerThread] + (*p_vi_RECORD))*e_d_IO_LAST;
  int vi_INT_IO_OFFSET = (getGlobalIdx_3D_3D()*d_CONST_INT[e_d_CONST_INT_numRecordsPerThread] + (*p_vi_RECORD))*e_d_IO_INT_LAST;
  
  d_IO_INT[vi_INT_IO_OFFSET + e_d_IO_INT_ERROR] = 0;
  d_IO_INT[vi_INT_IO_OFFSET + e_d_IO_INT_THREAD] = getGlobalIdx_3D_3D();
  d_IO_INT[vi_INT_IO_OFFSET + e_d_IO_INT_RECORD] = (*p_vi_RECORD);
  d_IO_INT[vi_INT_IO_OFFSET + e_d_IO_INT_RKQS_STEPS] = l_rkqs_TRIED; 
  
  d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_T] = l_time_CURRENT;
  d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_BX] = l_BField[0];
  d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_BY] = l_BField[1];
  d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_BZ] = l_BField[2];
  d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_SPINNOR_0] = l_spinor[0];
  d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_SPINNOR_1] = l_spinor[1];
  d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_SPINNOR_2] = l_spinor[2];
  d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_SPINNOR_3] = l_spinor[3];
  d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_DDT_SPINNOR_0] = l_dspinordt[0]; 
  d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_DDT_SPINNOR_1] = l_dspinordt[1]; 
  d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_DDT_SPINNOR_2] = l_dspinordt[2]; 
  d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_DDT_SPINNOR_3] = l_dspinordt[3]; 
  d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_ERR_SPINNOR_0] = l_epsilon[0]; 
  d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_ERR_SPINNOR_1] = l_epsilon[1]; 
  d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_ERR_SPINNOR_2] = l_epsilon[2]; 
  d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_ERR_SPINNOR_3] = l_epsilon[3];   
  d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_SCAL_SPINNOR_0] = l_spinor_scal[0]; 
  d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_SCAL_SPINNOR_1] = l_spinor_scal[1]; 
  d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_SCAL_SPINNOR_2] = l_spinor_scal[2]; 
  d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_SCAL_SPINNOR_3] = l_spinor_scal[3];
  d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_POLARIZATION] = l_pol;
  d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_HNEXT] = l_hnext;
  
  (*p_vi_RECORD)++;
  
  return 0;
}


__global__ void GENERIC_PIECEWISE_KERNEL_MULTI_XVS_RKQS_LOOP(
  double *d_IO, 
  int *d_IO_INT, 
  int numRecordsStart,
  int numRecordsEnd)
{
  int vi_RECORD_TEST_0, vi_RECORD_TEST_1, vi_RECORD_TEST_2, vi_RECORD_TEST_3,  vi_RECORD_TEST_4;
  int vi_RECORD_IO_OFFSET_END = 0, vi_IO_OFFSET_END, vi_IO_INT_OFFSET_END, vi_IO_DOUBLE_OFFSET_END;
  int vi_TESTOTESTO;
  int i, j, k, vi_RECORD, vi_RKQS_STEP, return_value_RKQS, l_rkqs_TRIED;
  int l_odeint_steps, vi_BREAK_FLAG, vi_REVERSE_FLAG, vi_INDEX;
  double l_time_CURRENT;
  double l_spinor[4];
  double l_dspinordt[4];
  double l_spinor_scal[4];
  double l_epsilon[4];
  double l_BField[3];
  double l_pol;
  double l_hnext, l_hdid, l_htry;
  int vi_RECORD_IO_OFFSET_START = getGlobalIdx_3D_3D()*d_CONST_INT[e_d_CONST_INT_numRecordsPerThread] + numRecordsStart;
  int vi_IO_OFFSET_START = vi_RECORD_IO_OFFSET_START*e_d_IO_LAST;
  int vi_IO_INT_OFFSET_START = vi_RECORD_IO_OFFSET_START*e_d_IO_INT_LAST;
  int vi_IO_DOUBLE_OFFSET_START = vi_RECORD_IO_OFFSET_START*e_d_IO_LAST;
  
  return_value_RKQS = d_IO[vi_IO_INT_OFFSET_START + e_d_IO_INT_ERROR];  
  l_rkqs_TRIED = 0;
  l_time_CURRENT = d_IO[vi_IO_DOUBLE_OFFSET_START + e_d_IO_T]; 
  for (i = 0 ; i<4; i++)
  {
    l_spinor[i] = d_IO[vi_IO_DOUBLE_OFFSET_START + e_d_IO_X + i];
    l_dspinordt[i] = d_IO[vi_IO_DOUBLE_OFFSET_START + e_d_IO_RED_VX + i];
    l_spinor_scal[i] = d_IO[vi_IO_DOUBLE_OFFSET_START + e_d_IO_SCAL_X + i];
    l_epsilon[i] = d_IO[vi_IO_DOUBLE_OFFSET_START + e_d_IO_ERR_X + i];
  }
  for (i = 0 ; i<3; i++)
  {
    l_BField[i] = d_IO[vi_IO_DOUBLE_OFFSET_START + e_d_IO_BX + i]; 
  }
  l_pol = d_IO[vi_IO_DOUBLE_OFFSET_START + e_d_IO_POLARIZATION]; 
  l_hnext = d_IO[vi_IO_DOUBLE_OFFSET_START + e_d_IO_HNEXT]; 
  
  double vd_TIME_START = l_time_CURRENT;
  
  for (i = 0 ; i<3; i++)
  {
    l_BField[i] = d_IO[vi_IO_DOUBLE_OFFSET_START + e_d_IO_BX + i]; 
    
  }
  l_pol = d_IO[vi_IO_DOUBLE_OFFSET_START + e_d_IO_POLARIZATION]; 
  l_hnext = d_IO[vi_IO_DOUBLE_OFFSET_START + e_d_IO_HNEXT]; 
  
  // CUDA_setspin( l_spinor, l_spin, l_spinnor, l_BField);
  int vi_CYCLE = 0;
  CUDA_derivs_XVS(l_time_CURRENT, l_spinor, l_dspinordt, l_BField);
  l_pol = CUDA_polcalc_XVS(l_BField[0], l_BField[1], l_BField[2], l_spinor);
  vi_REVERSE_FLAG = 0;
  // if (d_CONST[e_d_CONST_h1]<d_CONST[e_d_CONST_h1_SPIN]) l_htry = d_CONST[e_d_CONST_h1];
  // else l_htry = fmin(d_CONST[e_d_CONST_h1],d_CONST[e_d_CONST_h1_SPIN]);
  vi_RKQS_STEP = 0;
  return_value_RKQS = 0;
  while(vi_RECORD<numRecordsEnd)
  {
    l_htry = l_hnext;
    l_spinor_scal[0]=fabs(l_spinor[0])+fabs(l_dspinordt[0]*l_htry)+d_CONST[e_d_CONST_def_TINY1];
    l_spinor_scal[1]=fabs(l_spinor[1])+fabs(l_dspinordt[1]*l_htry)+d_CONST[e_d_CONST_def_TINY1];
    l_spinor_scal[2]=fabs(l_spinor[2])+fabs(l_dspinordt[2]*l_htry)+d_CONST[e_d_CONST_def_TINY1];
    l_spinor_scal[3]=fabs(l_spinor[3])+fabs(l_dspinordt[3]*l_htry)+d_CONST[e_d_CONST_def_TINY1];
    
    // if (0==0)/*(vi_RKQS_STEP>0 && (vi_RKQS_STEP%d_CONST_INT[e_d_CONST_INT_numCyclesPerRecord]))*/
      // {
        // if (d_CONST_INT[e_d_CONST_INT_numRecordsPerReverse]!=0)
        // {
          // if (d_CONST_INT[e_d_CONST_INT_numRecordsPerReverse_SPIN]==1 || 
            // (vi_RECORD%d_CONST_INT[e_d_CONST_INT_numRecordsPerReverse_SPIN])==0)
          // {
            // if (vi_REVERSE_FLAG==0) vi_REVERSE_FLAG = 1;
            // else vi_REVERSE_FLAG = 0;
          // }
        // }
        
      // }
    if (vi_RECORD==0)
    {
      CUDA_derivs_XVS(l_time_CURRENT, l_spinor, l_dspinordt, l_BField);
      l_pol = CUDA_polcalc_XVS(l_BField[0], l_BField[1], l_BField[2], l_spinor);
      CUDA_RECORD_XVS(d_IO, d_IO_INT,  &vi_RECORD, 
        l_time_CURRENT, l_spinor, l_epsilon, l_spinor_scal, l_dspinordt, 
        l_BField, l_pol, l_rkqs_TRIED, l_hnext);
    }
    else
    {
      return_value_RKQS = CUDA_rkqs_SINGLE_ATTEMPT_XVS(
        d_IO, 
        d_IO_INT, 
        vi_RECORD, 
        l_spinor, 
        l_dspinordt, 
        &l_time_CURRENT, 
        l_htry, 
        &l_hdid, 
        &l_hnext, 
        l_spinor_scal, 
        l_BField, 
        l_epsilon, 
        &l_rkqs_TRIED);
      if (return_value_RKQS==e_RKQS_ERROR_NONE)
      {
        l_time_CURRENT = l_time_CURRENT + l_hdid;
        CUDA_derivs_XVS(l_time_CURRENT, l_spinor, l_dspinordt, l_BField);
        l_pol = CUDA_polcalc_XVS(l_BField[0], l_BField[1], l_BField[2], l_spinor);
        // CUDA_RECORD_DOUBLE(d_IO, vi_RECORD, e_d_IO_EXTRA_1, l_hnext);
        vi_RKQS_STEP++;
        if (vi_RKQS_STEP>=d_CONST_INT[e_d_CONST_INT_numCyclesPerRecord])
        {
          // CUDA_RECORD_DOUBLE(d_IO, vi_RECORD, e_d_IO_EXTRA_0, l_time_CURRENT);
          CUDA_RECORD_XVS(d_IO, d_IO_INT,  &vi_RECORD, 
            l_time_CURRENT, l_spinor, l_epsilon, l_spinor_scal, l_dspinordt, 
            l_BField, l_pol, l_rkqs_TRIED, l_hnext);
          vi_RKQS_STEP = 0;
        }
        l_rkqs_TRIED = 0;
      }
      else l_rkqs_TRIED++;
      
      // CUDA_RECORD_DOUBLE(d_IO, vi_RECORD, e_d_IO_EXTRA_1, l_hnext);
      // CUDA_RECORD_DOUBLE(d_IO, vi_RECORD, e_d_IO_EXTRA_2, l_htry);
      // CUDA_RECORD_DOUBLE(d_IO, vi_RECORD, e_d_IO_EXTRA_3, l_hnext);
    }
  }
  
  // vi_RECORD_TEST_0 = vi_THREAD_IO_OFFSET;
  // vi_RECORD_TEST_1 = d_CONST_INT[e_d_CONST_INT_numRecordsPerThread];
  // vi_RECORD_TEST_2 = vi_THREAD_IO_OFFSET*d_CONST_INT[e_d_CONST_INT_numRecordsPerThread];
  // vi_RECORD_TEST_3 = numRecordsEnd;
  // vi_RECORD_TEST_4 = vi_THREAD_IO_OFFSET*d_CONST_INT[e_d_CONST_INT_numRecordsPerThread] + numRecordsEnd;
  
  int testsssss = getGlobalIdx_3D_3D();
  // vi_RECORD_IO_OFFSET_END = getGlobalIdx_3D_3D()*d_CONST_INT[e_d_CONST_INT_numRecordsPerThread] + numRecordsEnd;
  // vi_IO_INT_OFFSET_END = vi_RECORD_IO_OFFSET_END*e_d_IO_INT_LAST;
  vi_IO_INT_OFFSET_END = (testsssss * d_CONST_INT[e_d_CONST_INT_numRecordsPerThread] + numRecordsEnd) * e_d_IO_INT_LAST;
  // vi_IO_DOUBLE_OFFSET_END = vi_RECORD_IO_OFFSET_END*e_d_IO_LAST;
  vi_IO_DOUBLE_OFFSET_END = (testsssss * d_CONST_INT[e_d_CONST_INT_numRecordsPerThread] + numRecordsEnd)*e_d_IO_LAST;

  // d_IO[vi_IO_INT_OFFSET_END + e_d_IO_INT_RKQS_ERROR] = vi_TESTOTESTO;  
  // d_IO[vi_IO_INT_OFFSET_END + e_d_IO_INT_RKQS_STEPS] = l_rkqs_TRIED;
  
  d_IO[vi_IO_DOUBLE_OFFSET_END + e_d_IO_T] = l_time_CURRENT; 
  for (i = 0 ; i<4; i++)
  {
    d_IO[vi_IO_DOUBLE_OFFSET_END + e_d_IO_X + i] = l_spinor[i];
    d_IO[vi_IO_DOUBLE_OFFSET_END + e_d_IO_RED_VX + i] =   l_dspinordt[i];
    d_IO[vi_IO_DOUBLE_OFFSET_END + e_d_IO_SCAL_X + i] =   l_spinor_scal[i];
    d_IO[vi_IO_DOUBLE_OFFSET_END + e_d_IO_ERR_X + i] =   l_epsilon[i];
  }
  for (i = 0 ; i<3; i++)
  {
    d_IO[vi_IO_DOUBLE_OFFSET_END + e_d_IO_BX + i] = l_BField[i]; 
  }
  d_IO[vi_IO_DOUBLE_OFFSET_END + e_d_IO_POLARIZATION] = l_pol; 
  d_IO[vi_IO_DOUBLE_OFFSET_END + e_d_IO_HNEXT] = l_hnext; 
  
  return;
}
void GENERIC_RECORD_FRAME(
  int param_numBlocks,
  int param_numThreadsPerBlock,
  double *d_IO,
  int *d_IO_INT,
  int numRecordsStart,
  int numRecordsEnd)
{
  // cudaPrintfInit ();
    // Run kernel
  cudaThreadSynchronize();
  GENERIC_PIECEWISE_KERNEL_MULTI_XVS_RKQS_LOOP<<< param_numBlocks, param_numThreadsPerBlock >>>(
    d_IO, 
    d_IO_INT, 
    numRecordsStart,
    numRecordsEnd);
  // GENERIC_PIECEWISE_KERNEL_MULTI_XVS<<<h_CONST_INT[e_d_CONST_INT_numBlocks],h_CONST_INT[e_d_CONST_INT_numThreadsPerBlock]>>>(
    // d_IO, 
    // d_IO_INT, 
    // vi_RecordsStartCurrent,
    // vi_RecordsEndCurrent);
  // GENERIC_PIECEWISE_KERNEL_MULTI<<<h_CONST_INT[e_d_CONST_INT_numBlocks],h_CONST_INT[e_d_CONST_INT_numThreadsPerBlock]>>>(
    // d_IO, 
    // d_IO_INT, 
    // vi_RecordsStartCurrent,
    // vi_RecordsEndCurrent);
  // cudaPrintfDisplay(stdout,true);
  // cudaPrintEnd();
  cudaThreadSynchronize();
  return;
}

void GENERIC_MIDDLEMAN_MULTI(
  const double *h_CONST,
  const int *h_CONST_INT,
  double *h_IO,
  int *h_IO_INT)
{
  printf("\n\n\nelloel,lleoeooo%d shoudl be zero\n\n",e_RKQS_ERROR_NONE);
  // Establish Scope parameters for simulation: 
  // number of neutrons, for how long, how many records to keep, etc.
  
  // Copy passsed constant values from Host Memory (CPU front side bus RAM) to Device Constant Memory.
  // Device Constant Memory is limited in size but accessible with close to register level latency at the the thread level due to mandatory caching in every CUDA Multi-processor.
  
  int vi_INDEX;
  
  double h_UNOFFICIAL_CONST[e_d_CONST_LAST];
  for (vi_INDEX = 0; vi_INDEX<e_d_CONST_LAST; vi_INDEX++) h_UNOFFICIAL_CONST[vi_INDEX] = h_CONST[vi_INDEX];
  int h_UNOFFICIAL_CONST_INT[e_d_CONST_INT_LAST];
  for (vi_INDEX = 0; vi_INDEX<e_d_CONST_INT_LAST; vi_INDEX++) h_UNOFFICIAL_CONST_INT[vi_INDEX] = h_CONST_INT[vi_INDEX];
  
  const int numBytesCONST = e_d_CONST_LAST*sizeof(double);
  const int numBytesCONST_INT = e_d_CONST_INT_LAST*sizeof(int);
  
  int vi_ERROR = cudaMemcpyToSymbol(d_CONST,h_UNOFFICIAL_CONST,numBytesCONST);
  int vi_ERROR_INT = cudaMemcpyToSymbol(d_CONST_INT,h_UNOFFICIAL_CONST_INT,numBytesCONST_INT);
  
  // Move passed input parameters specific to each thread from Host Memory to Device Memory (on-card RAM)
  double *d_IO = NULL;
  int *d_IO_INT = NULL;
  const int numBytesIO = h_CONST_INT[e_d_CONST_INT_numRecords]*e_d_IO_LAST*sizeof(double);
  const int numBytesIO_INT = h_CONST_INT[e_d_CONST_INT_numRecords]*e_d_IO_INT_LAST*sizeof(int);
  cudaMalloc((void**)&d_IO, numBytesIO);
  cudaMalloc((void**)&d_IO_INT, numBytesIO_INT);
  cudaMemcpy(d_IO, h_IO, numBytesIO, cudaMemcpyHostToDevice);
  cudaMemcpy(d_IO_INT, h_IO_INT, numBytesIO_INT, cudaMemcpyHostToDevice);
  
  printf("\ncheck eeeee check check");
  cudaThreadSynchronize();
  // Allocate where output data goes with room for all threads
  // double *d_OUT = NULL;
  // int *d_OUT_INT = NULL;
  // const int numBytesOUT = h_CONST_INT[e_d_CONST_INT_numRecords]*e_d_IO_LAST*sizeof(double);
  // const int numBytesOUT_INT = h_CONST_INT[e_d_CONST_INT_numRecords]*e_d_IO_INT_LAST*sizeof(int);
  // cudaMalloc((void**)&d_OUT, numBytesOUT);
  // cudaMalloc((void**)&d_OUT_INT, numBytesOUT_INT);
  int vi_RecordsStartCurrent = 0;
  int vi_RecordsEndCurrent = 0;
  int vi_RecordsEndFinal = h_CONST_INT[e_d_CONST_INT_numRecordsPerThread];
  int vi_UPDATES;
  for (vi_UPDATES = 0; vi_RecordsEndCurrent<vi_RecordsEndFinal; vi_UPDATES++)
  {
    // Run kernel
    vi_RecordsStartCurrent = vi_RecordsEndCurrent;
    int vi_RecordsEndCurrent_TEMP =  vi_RecordsStartCurrent + h_CONST_INT[e_d_CONST_INT_numRecordsPerUpdate];
    if (vi_RecordsEndCurrent_TEMP>vi_RecordsEndFinal) vi_RecordsEndCurrent = vi_RecordsEndFinal;
    else vi_RecordsEndCurrent = vi_RecordsStartCurrent + h_CONST_INT[e_d_CONST_INT_numRecordsPerUpdate];
    // vi_RecordsStartCurrent = 0;
    // vi_RecordsEndCurrent = vi_RecordsEndFinal;
    printf("\nStarting Records %d-%d of %d... ",vi_RecordsStartCurrent,vi_RecordsEndCurrent,h_CONST_INT[e_d_CONST_INT_numRecordsPerThread]);
    GENERIC_RECORD_FRAME(
      h_CONST_INT[e_d_CONST_INT_numBlocks],
      h_CONST_INT[e_d_CONST_INT_numThreadsPerBlock],
      d_IO,
      d_IO_INT,
      vi_RecordsStartCurrent,
      vi_RecordsEndCurrent);
    // GENERIC_PIECEWISE_KERNEL_MULTI<<<h_CONST_INT[e_d_CONST_INT_numBlocks],h_CONST_INT[e_d_CONST_INT_numThreadsPerBlock]>>>(
      // d_IO, 
      // d_IO_INT, 
      // vi_RecordsStartCurrent,
      // vi_RecordsEndCurrent);
    printf("Completed");
  }
  printf("\nDay o day o daylight come and me wanna go %d threads",h_CONST_INT[e_d_CONST_INT_numThreads]);
  // Move results of output to Host Memory from Device Memory
  cudaMemcpy(h_IO, d_IO, numBytesIO, cudaMemcpyDeviceToHost);
  cudaMemcpy(h_IO_INT, d_IO_INT, numBytesIO_INT, cudaMemcpyDeviceToHost);
  
  cudaFree(d_IO);
  cudaFree(d_IO_INT);
  
  return;
}
