
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
		(blockIdx.y * gridDim.x)	+ (
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
__device__ /*__host__*/ double CUDA_SPIN_X( double spinor_0, double spinor_1, double spinor_2, double spinor_3 )
{
	double xhat1;
  if (d_CONST[e_d_CONST_INT_FLAG_CLASSICAL_SPIN]==0) xhat1 = spinor_1;
	else if (d_CONST[e_d_CONST_INT_FLAG_CLASSICAL_SPIN]==1)
	{
		double snorm = (spinor_0*spinor_0) + 	(spinor_1*spinor_1) + (spinor_2*spinor_2) + (spinor_3*spinor_3);
		xhat1 =( 2.*((spinor_0*spinor_2) + (spinor_1*spinor_3)))/snorm;
	}
	return xhat1;
}
__device__ /*__host__*/ double CUDA_SPIN_Y( double spinor_0, double spinor_1, double spinor_2, double spinor_3 )
{
	double yhat1;
  if (d_CONST[e_d_CONST_INT_FLAG_CLASSICAL_SPIN]==0) yhat1 = spinor_2;
	else if (d_CONST[e_d_CONST_INT_FLAG_CLASSICAL_SPIN]==1)\
	{
		double snorm = (spinor_0*spinor_0) + (spinor_1*spinor_1) +  
						(spinor_2*spinor_2) + (spinor_3*spinor_3);
		yhat1 = ( 2.*((spinor_0*spinor_3) - (spinor_1*spinor_2)))/snorm;
	}
  return yhat1;
}
__device__ /*__host__*/ double CUDA_SPIN_Z( double spinor_0, double spinor_1, double spinor_2, double spinor_3 )
{
  double zhat1;
  if (d_CONST[e_d_CONST_INT_FLAG_CLASSICAL_SPIN]==0) zhat1 = spinor_3;
	else if (d_CONST[e_d_CONST_INT_FLAG_CLASSICAL_SPIN]==1)
	{
		double snorm = (spinor_0*spinor_0) + (spinor_1*spinor_1) +  
          (spinor_2*spinor_2) + (spinor_3*spinor_3);
		zhat1 = ((spinor_0*spinor_0) + (spinor_1*spinor_1) - 
	   (spinor_2*spinor_2) - (spinor_3*spinor_3))/snorm;
	}
	return zhat1;
}
__device__ /*__host__*/ int CUDA_SPIN_XVS_XYZ( 
	double spin[],
	double xvs[])
{	
	spin[0] = CUDA_SPIN_X(xvs[6], xvs[7], xvs[8], xvs[9]);
	spin[1] = CUDA_SPIN_Y(xvs[6], xvs[7], xvs[8], xvs[9]);
	spin[2] = CUDA_SPIN_Z(xvs[6], xvs[7], xvs[8], xvs[9]);
}
__device__ /*__host__*/ double CUDA_SPIN_XVS_SINGLE( 
	double xvs[],
	int vi_X0_Y1_Z2)
{
	double vd_RETURN_VALUE;
	if (vi_X0_Y1_Z2==0) vd_RETURN_VALUE = CUDA_SPIN_X(xvs[6], xvs[7], xvs[8], xvs[9]);
	if (vi_X0_Y1_Z2==1) vd_RETURN_VALUE = CUDA_SPIN_Y(xvs[6], xvs[7], xvs[8], xvs[9]);
	if (vi_X0_Y1_Z2==2) vd_RETURN_VALUE = CUDA_SPIN_Z(xvs[6], xvs[7], xvs[8], xvs[9]);
	return vd_RETURN_VALUE;
}
__device__ /*__host__*/ double CUDA_polcalc_XVS( 
	double bb1x, 
	double bb1y, 
	double bb1z, 
	double xvs[])
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
	// CUDA_SPIN_XVS_XYZ( spin, xvs );
  snorm = (xvs[6]*xvs[6]) + (xvs[7]*xvs[7]) +  
          (xvs[8]*xvs[8]) + (xvs[9]*xvs[9]);
  xhat1 =( 2.*((xvs[6]*xvs[8]) + (xvs[7]*xvs[9])))/snorm;
  yhat1 = ( 2.*((xvs[6]*xvs[9]) - (xvs[7]*xvs[8])))/snorm;
  zhat1 = ((xvs[6]*xvs[6]) + (xvs[7]*xvs[7]) - 
	   (xvs[8]*xvs[8]) - (xvs[9]*xvs[9]))/snorm;
  polarization = xhat1*bxhat + yhat1*byhat + zhat1*bzhat;
  return(polarization);
}
__device__ /*__host__*/ int CUDA_Mag_ROHM_HOLLEY_XVS( double xvs[], double BField[] )
{
	double xi, yi, zi, l, zoff, zoff2, zcent, a0, rcoil, lcoil, rsol, bmax, b0, bpr, bpp, bppp, zrf, lsol,sinth0,costh0, bc0, sinthc0, costhc0, bcc0;
	zoff = 0.295;
  zoff2 = 1.595;
  rsol = .08;
  b0   = .975372;/* .93562; */

  a0 = .21;

  zcent = 0.1641;
  bmax = 7.0327;
  rcoil = .144;
  lcoil = .097;

  bpr  = 0.10749; 
  bpp = -0.102;  /*T/m^2*/
  bppp = 3.45;  /* T/m^3 */

  zrf = 1.14;
  lsol = (zoff2-zoff)/2.;
  sinth0 = rsol/sqrt(rsol*rsol+lsol*lsol);
  costh0 = lsol*sinth0/rsol;
  bc0 = b0/costh0;
 
  sinthc0 = rcoil/sqrt(rcoil*rcoil + lcoil*lcoil);
  costhc0 = lcoil*sinthc0/rcoil;
  bcc0  = bmax/costhc0;
  double rho;


  double zd1, anm,anmsq,rhnm1;
  double costh1,sinth1,costh2,sinth2;
  double sinth1sq,sinth2sq;
  double rhnm2;
  int result;

  double sinthc1,costhc1,sinthc2,costhc2;
  double sinthc1sq,sinthc2sq;
  double rhnmc;

  /*
  double sinthn1,costhn1,sinthn2,costhn2;
  double sinthn1sq,sinthn2sq;
  double rhnmn;
  */

  /* compute common expressions 
   (this is for legibility as well as optimization) */
  sinth1 = rsol/sqrt(rsol*rsol+(zoff2-xvs[2])*(zoff2-xvs[2]));
  costh1 = (zoff2 - xvs[2])*sinth1/rsol;

  sinth2 = rsol/sqrt(rsol*rsol+(xvs[2]-zoff)*(xvs[2]-zoff));
  costh2 = (xvs[2] - zoff)*sinth2/rsol;

  sinth1sq = sinth1*sinth1;
  sinth2sq = sinth2*sinth2;

  rho = sqrt(xvs[0]*xvs[0]+xvs[1]*xvs[1]);
  rhnm2 = rho/rsol;
  rhnmc = rho/rcoil;
  /*  rhnmn = rho/rnew; */

  sinthc1 = rcoil/sqrt(rcoil*rcoil+(zcent+lcoil-xvs[2])*(zcent+lcoil-xvs[2]));
  costhc1 = (zcent+lcoil - xvs[2])*sinthc1/rcoil;

  sinthc2 = rcoil/sqrt(rcoil*rcoil + (xvs[2] - (zcent - lcoil))*(xvs[2] - (zcent - lcoil)));
  costhc2 = (xvs[2]-(zcent - lcoil))*sinthc2/rcoil;

  sinthc1sq = sinthc1*sinthc1;
  sinthc2sq = sinthc2*sinthc2;
  /*
  sinthn1 = rnew/sqrt(rnew*rnew+(zcnew+lnew-xvs[2])*(zcnew+lnew-xvs[2]));
  costhn1 = (zcnew+lnew - xvs[2])*sinthn1/rnew;

  sinthn2 = rnew/sqrt(rnew*rnew + (xvs[2] - (zcnew - lnew))*(xvs[2] - (zcnew - lnew)));
  costhn2 = (xvs[2]-(zcnew - lnew))*sinthn2/rnew;

  sinthn1sq = sinthn1*sinthn1;
  sinthn2sq = sinthn2*sinthn2;
  */
  zd1 = sqrt(a0*a0+(xvs[2]-zcent)*(xvs[2]-zcent));
  anm = a0/zd1;
  anmsq = anm*anm;
  rhnm1 = rho/zd1;
  result = 0;
  BField[0] = bcc0*xvs[0]*(sinthc1*sinthc1sq-sinthc2*sinthc2sq)/(4.*rcoil)
    /*    + bcn0*xvs[0]*(sinthn1*sinthn1sq-sinthn2*sinthn2sq)/(4.*rnew)  */
    + bc0*xvs[0]*(sinth1*sinth1sq-sinth2*sinth2sq)/(4.*rsol)
    - xvs[0]*bpr/2. - xvs[0]*(xvs[2]-zrf)*bpp/2. + xvs[0]*(rho*rho/16. -(xvs[2]-zrf)*(xvs[2]-zrf)/4.)*bppp;
  BField[1] =  bcc0*xvs[1]*(sinthc1*sinthc1sq-sinthc2*sinthc2sq)/(4.*rcoil)
    /*   + bcn0*xvs[1]*(sinthn1*sinthn1sq-sinthn2*sinthn2sq)/(4.*rnew) */
    + bc0*xvs[1]*(sinth1*sinth1sq-sinth2*sinth2sq)/(4.*rsol)
    - xvs[1]*bpr/2. - xvs[1]*(xvs[2]-zrf)*bpp/2. + xvs[1]*(rho*rho/16. -(xvs[2]-zrf)*(xvs[2]-zrf)/4.)*bppp;
  BField[2] =  bcc0*((costhc1+costhc2)/2.+
		     .375*rhnmc*rhnmc*(sinthc1sq*sinthc1sq*costhc1+sinthc2sq*sinthc2sq*costhc2))
    /*    + bcn0*((costhn1+costhn2)/2.+
	  .375*rhnmn*rhnmn*(sinthn1sq*sinthn1sq*costhn1+sinthn2sq*sinthn2sq*costhn2)) */
    + bc0*((costh1+costh2)/2.
	   +.375*rhnm2*rhnm2*(sinth1sq*sinth1sq*costh1+sinth2sq*sinth2sq*costh2))
    + bpr*(xvs[2]-zrf) + ((xvs[2]-zrf)*(xvs[2]-zrf)/2. - rho*rho/4.)*bpp + ((xvs[2]-zrf)*((xvs[2]-zrf)*(xvs[2]-zrf)/6.-rho*rho/4.))*bppp;
 
 return( result );
}


__device__ /*__host__*/ int CUDA_dB_ROHM_HOLLEY_XVS( double xvs[], double BField[], double dField_1D_FLAT[])
{
	double xi, yi, zi, l, zoff, zoff2, zcent, a0, rcoil, lcoil, rsol, bmax, b0, bpr, bpp, bppp, zrf, lsol,sinth0,costh0, bc0, sinthc0, costhc0, bcc0;
	zoff = 0.295;
  zoff2 = 1.595;
  rsol = .08;
  b0   = .975372;/* .93562; */

  a0 = .21;

  zcent = 0.1641;
  bmax = 7.0327;
  rcoil = .144;
  lcoil = .097;

  bpr  = 0.10749; 
  bpp = -0.102;  /*T/m^2*/
  bppp = 3.45;  /* T/m^3 */

  zrf = 1.14;
  lsol = (zoff2-zoff)/2.;
  sinth0 = rsol/sqrt(rsol*rsol+lsol*lsol);
  costh0 = lsol*sinth0/rsol;
  bc0 = b0/costh0;
 
  sinthc0 = rcoil/sqrt(rcoil*rcoil + lcoil*lcoil);
  costhc0 = lcoil*sinthc0/rcoil;
  bcc0  = bmax/costhc0;	
  double rho;   
  double zd1sq,anmsq,rhnm1;
  double costh1,sinth1,costh2,sinth2;
  double sinth1sq,sinth2sq;
  double  rhnm2;

  double costhc1,sinthc1,costhc2,sinthc2;
  double sinthc1sq,sinthc2sq;
  double  rhnmc2;
  double BZ;
  double Bslope, Bdelta, Bsigma;
  
  int result;


  /*  !!!!!!! NOT UP TO DATE   !!!!!   */


  rho = sqrt(xvs[0]*xvs[0]+xvs[1]*xvs[1]);

  sinth1 = rsol/sqrt(rsol*rsol+(zoff2-xvs[2])*(zoff2-xvs[2]));
  costh1 = (zoff2 - xvs[2])*sinth1/rsol;
  sinth2 = rsol/sqrt(rsol*rsol+(xvs[2]-(zoff2-2.*lsol))*(xvs[2]-(zoff2-2.*lsol)));
  costh2 = (xvs[2]-(zoff2-2.*lsol))*sinth2/rsol;
  sinth1sq = sinth1*sinth1;
  sinth2sq = sinth2*sinth2;
 
  rhnm2 = rho/rsol;

  sinthc1 = rcoil/sqrt(rcoil*rcoil+(zcent+lcoil-xvs[2])*(zcent+lcoil-xvs[2]));
  costhc1 = (zcent + lcoil - xvs[2])*sinthc1/rcoil;
  sinthc2 = rcoil/sqrt(rcoil*rcoil+(xvs[2]-(zcent-lcoil))*(xvs[2]-(zcent-lcoil)));
  costhc2 = (xvs[2]-(zcent-lcoil))*sinthc2/rcoil;
  sinthc1sq = sinthc1*sinthc1;
  sinthc2sq = sinthc2*sinthc2;
  
  rhnmc2 = rho/rcoil;  

  zd1sq = (a0*a0+(xvs[2]-zcent)*(xvs[2]-zcent));   /* ! */
  rhnm1 = rho/sqrt(zd1sq);          /* ! */
  anmsq = a0*a0/zd1sq;             /* ! */
    /*	anmsq = anm*anm;*/
  
  BZ = bmax*anmsq*sqrt(anmsq);        /* ! */
  result = 0; /* if direct1,2 <0 or >3, we are at sea. */  
  /*   */
  Bslope = (xvs[2]-zcent)*BZ/zd1sq;      /* ! */
  Bdelta =  bc0*(sinth1*sinth1sq-sinth2*sinth2sq)/(4.*rsol) +
          bcc0*(sinthc1*sinthc1sq-sinthc2*sinthc2sq)/(4.*rcoil);
  /*  fixed sinth4 bug here */
  Bsigma = bc0*(sinth1sq*sinth1sq*costh1+sinth2sq*sinth2sq*costh2) +
    bcc0*(sinthc1sq*sinthc1sq*costhc1+sinthc2sq*sinthc2sq*costhc2);

	/* */
  dField_1D_FLAT[0] = Bdelta - bpr/2. - (xvs[2]-zrf)*bpp/2.+(3.*xvs[0]*xvs[0]+xvs[1]*xvs[1])*bppp/16.;
  dField_1D_FLAT[1] = xvs[0]*xvs[1]*bppp/8.;
  dField_1D_FLAT[2] =  .75 * xvs[0] * Bsigma - xvs[0]*(bpp+(xvs[2]-zrf)*bppp)/2.;
  dField_1D_FLAT[3] =  xvs[0]*xvs[1]*bppp/8.;
  dField_1D_FLAT[4] =  Bdelta - bpr/2. - (xvs[2]-zrf)*bpp/2.+(xvs[0]*xvs[0]+3.*xvs[1]*xvs[1])*bppp/16.;
  dField_1D_FLAT[5] =  .75 * xvs[1] * Bsigma - xvs[1]*(bpp+(xvs[2]-zrf)*bppp)/2.;
  dField_1D_FLAT[6] = - xvs[0]*(bpp+(xvs[2]-zrf)*bppp)/2.;
  dField_1D_FLAT[7] = - xvs[1]*(bpp+(xvs[2]-zrf)*bppp)/2.;
  dField_1D_FLAT[8] =  .5*(bc0/rsol)*( -sinth1*sinth1sq + sinth2*sinth2sq)+ 0.375*(bc0/rsol)
    *(rhnm2*rhnm2)*(sinth1sq*sinth1*(4.*sinth1sq-5.*sinth1sq*sinth1sq))
    - 0.375*(bc0/rsol)*(rhnm2*rhnm2)*(sinth2sq*sinth2*(4.*sinth2sq-5.*sinth2sq*sinth2sq))
    + .5*(bcc0/rcoil)*( -sinthc1*sinthc1sq + sinthc2*sinthc2sq)+ 0.375*(bcc0/rcoil)
    *(rhnmc2*rhnmc2)*(sinthc1sq*sinthc1*(4.*sinthc1sq-5.*sinthc1sq*sinthc1sq))
    - 0.375*(bcc0/rcoil)*(rhnmc2*rhnmc2)*(sinthc2sq*sinthc2*(4.*sinthc2sq-5.*sinthc2sq*sinthc2sq))
    + bpr + (xvs[2]-zrf)*bpp + (2.*(xvs[2]-zrf)*(xvs[2]-zrf)- (xvs[0]*xvs[0]+xvs[1]*xvs[1]))*bppp/4.;
  /* result = dField[comp][part]; */
	/* return(result); */
	return 0;
}
__device__ /*__host__*/ int CUDA_Mag_CONSTANT_XVS( double xvs[], double BField[])
{
	BField[0] = 0.0;
	BField[1] = 0.0;
	BField[2] = 1.0;
	return 0;
}
__device__ /*__host__*/ int CUDA_dB_CONSTANT_XVS( double xvs[], double BField[], double dField_1D_FLAT[])
{
	dField_1D_FLAT[0] = 0.0;
	dField_1D_FLAT[1] = 0.0;
	dField_1D_FLAT[2] = 0.0;
	dField_1D_FLAT[3] = 0.0;
	dField_1D_FLAT[4] = 0.0;
	dField_1D_FLAT[5] = 0.0;
	dField_1D_FLAT[6] = 0.0;
	dField_1D_FLAT[7] = 0.0;
	dField_1D_FLAT[8] = 0.0;
	return 0;
}
__device__ /*__host__*/ int CUDA_Mag_XVS( double xvs[], double BField[])
{
	if (d_CONST_INT[e_d_CONST_INT_FLAG_MAGNETIC]==1) return CUDA_Mag_ROHM_HOLLEY_XVS( xvs, BField );
	else if (d_CONST_INT[e_d_CONST_INT_FLAG_MAGNETIC]==2) return CUDA_Mag_CONSTANT_XVS( xvs, BField );
	else return -1;
}
__device__ /*__host__*/ int CUDA_dB_XVS( double xvs[], double BField[], double dField_1D_FLAT[])
{
	if (d_CONST_INT[e_d_CONST_INT_FLAG_MAGNETIC]==1) return CUDA_dB_ROHM_HOLLEY_XVS( xvs, BField, dField_1D_FLAT);
	else if (d_CONST_INT[e_d_CONST_INT_FLAG_MAGNETIC]==2) return CUDA_dB_CONSTANT_XVS( xvs, BField, dField_1D_FLAT);
	else return -1;
}
__device__ /*__host__*/ int CUDA_RF_BFIELD( double t, double xvs[],	double extra_brf[])
{
		double rfstr = d_CONST[e_d_CONST_RF_BFIELD_MAG];
		double  omega,btoomega,zrf,zdist,zdelta;
		btoomega = d_CONST[e_d_CONST_def_MOMENT]/d_CONST[e_d_CONST_def_HBAR];  // coeff. for psi-dot 
		omega    = 2.0*d_CONST[e_d_CONST_def_MOMENT]*(1.0)/d_CONST[e_d_CONST_def_HBAR]; // reson. frequ. for 1.0 T field 
		zrf      = 1.14; // chosen to match crossing poconst int 
		// some comment explaining this. v 
		zdist = xvs[8];
		zdelta = (zdist - zrf)*(zdist -zrf)/(.05*.05);
		extra_brf[0] = rfstr*cos(omega*(t-0.0)) * (1.0/((1.0+zdelta)*(1.0+zdelta))); 
		extra_brf[1] = 0.;
			// exp(-(zdist-zrf)*(zdist-zrf)/(.05*.05)); 
		extra_brf[2] = 0.;
		
		return 0;
}
__device__ /*__host__*/ int CUDA_derivs_XVS( 
	double t,
	double xvs[],
	double dxvsdt[],
	double BField[],
	double dField_1D_FLAT[])
{	
	double x,y,z;	/* Cartesian coordinates */
  double rad, th;	/* Polar coordinates */
  double B; 	/* Magnitude of B */
  int Bgood;
  double spin[3];
	CUDA_SPIN_XVS_XYZ(spin, xvs);
	x=xvs[0];
  y=xvs[1];
  z=xvs[2];
  rad=CUDA_r(x,y,z);
  th=CUDA_theta(x,y,z);
  dxvsdt[0]=xvs[3]; /* derivative of x = v in x direction */
  dxvsdt[1]=xvs[4]; /* ditto for y */
  dxvsdt[2]=xvs[5]; /* ditto for z */
  dxvsdt[3]=0; /* initializing acceleration in x direction */
  dxvsdt[4]=0; /* ditto for y */
  dxvsdt[5]=0; /* ditto for z */
  Bgood = CUDA_Mag_XVS(xvs, BField);
  B = sqrt( BField[0]*BField[0]+BField[1]*BField[1]+BField[2]*BField[2]);
  Bgood = CUDA_dB_XVS(xvs, BField, dField_1D_FLAT);
	// if (d_CONST_INT[e_d_CONST_INT_def_PERFECT_POLARIZATION] && 
		// isfinite(BField[0]) && 
		// isfinite(BField[1]) && 
		// isfinite(BField[2]) && 
		// (B>0 && isfinite(B)))
	// {
		// spin[0] = BField[0]/B;
		// spin[1] = BField[1]/B;
		// spin[2] = BField[2]/B;
	// }
  if (d_CONST_INT[e_d_CONST_INT_FLAG_MAGNETIC]==1 && B>0 && isfinite(B))
	{
		dxvsdt[3] += - d_CONST[e_d_CONST_def_CORRECTIVE_FACTOR_SPIN]*(d_CONST[e_d_CONST_def_MOMENT_DIV_MASS])*spin[0] / B * 
			(BField[0] * dField_1D_FLAT[0] +     BField[1] * dField_1D_FLAT[1] + BField[2] * dField_1D_FLAT[2]);
		dxvsdt[4] += - d_CONST[e_d_CONST_def_CORRECTIVE_FACTOR_SPIN]*(d_CONST[e_d_CONST_def_MOMENT_DIV_MASS])*spin[1] / B * 
			(BField[0] * dField_1D_FLAT[3] + BField[1] * dField_1D_FLAT[4] +  BField[2] * dField_1D_FLAT[5]);
		dxvsdt[5] += - d_CONST[e_d_CONST_def_CORRECTIVE_FACTOR_SPIN]*(d_CONST[e_d_CONST_def_MOMENT_DIV_MASS])*spin[2] / B * 
			(BField[0] * dField_1D_FLAT[6] + BField[1] * dField_1D_FLAT[7] + BField[2] * dField_1D_FLAT[8]);
  }
	// dxvsdt[3] += - (d_CONST[e_d_CONST_def_MOMENT_DIV_MASS])*spin[0] * 
		// (cos(d_CONST[e_d_CONST_def_HALBACH_K]*xvs[0]) * dField_1D_FLAT[0] + sin(d_CONST[e_d_CONST_def_HALBACH_K]*xvs[0]) * dField_1D_FLAT[1]);
	// dxvsdt[4] += - (d_CONST[e_d_CONST_def_MOMENT_DIV_MASS])*spin[1] * 
		// (cos(d_CONST[e_d_CONST_def_HALBACH_K]*xvs[0]) * dField_1D_FLAT[3] + sin(d_CONST[e_d_CONST_def_HALBACH_K]*xvs[0]) * dField_1D_FLAT[4] );
	// dxvsdt[5] += - (d_CONST[e_d_CONST_def_MOMENT_DIV_MASS])*spin[2] * 
		// (cos(d_CONST[e_d_CONST_def_HALBACH_K]*xvs[0]) * dField_1D_FLAT[6] + sin(d_CONST[e_d_CONST_def_HALBACH_K]*xvs[0]) * dField_1D_FLAT[7] );
	// if (B>0 && isfinite(B)) dxvsdt[4] += exp(-2*d_CONST[e_d_CONST_def_HALBACH_K]*xvs[1])*spin[1]*d_CONST[e_d_CONST_def_HALBACH_K]*d_CONST[e_d_CONST_def_HALBACH_MAX_TESLA]*d_CONST[e_d_CONST_def_HALBACH_MAX_TESLA]*d_CONST[e_d_CONST_def_CORRECTIVE_FACTOR_SPIN]*(d_CONST[e_d_CONST_def_MOMENT_DIV_MASS]);
	/* Gravity is in +/- y direction depending on sign */ 
  if (d_CONST_INT[e_d_CONST_INT_FLAG_GRAVITY]==1) dxvsdt[4]  += d_CONST[e_d_CONST_def_GRAVITY];
	/* Ground spring begins at y=0 and below */
	if (d_CONST_INT[e_d_CONST_INT_FLAG_SPRING]==1)
	{
		if (xvs[1]<0) dxvsdt[4] += -xvs[1]*d_CONST[e_d_CONST_def_SPRING];
	}
	double rfstr = d_CONST[e_d_CONST_RF_BFIELD_MAG];
	double  brf[4],omega,btoomega,zrf,zdist,zdelta;
	btoomega = d_CONST[e_d_CONST_def_MOMENT]/d_CONST[e_d_CONST_def_HBAR];  // coeff. for psi-dot 
	omega    = 2.0*d_CONST[e_d_CONST_def_MOMENT]*(1.0)/d_CONST[e_d_CONST_def_HBAR]; // reson. frequ. for 1.0 T field 
	zrf      = 1.14; // chosen to match crossing poconst int 
	// some comment explaining this. v 
	zdist = spin[2];
	zdelta = (zdist - zrf)*(zdist -zrf)/(.05*.05);
	brf[0] = rfstr*cos(omega*(t-0.0)) * (1.0/((1.0+zdelta)*(1.0+zdelta))); 
	brf[1] = 0.;
		// exp(-(zdist-zrf)*(zdist-zrf)/(.05*.05)); 
	brf[2] = 0.;
	int i;
	double extra_brf[3];
	CUDA_RF_BFIELD(t, xvs, extra_brf);
	if (d_CONST_INT[e_d_CONST_INT_FLAG_RF]==1)
	{
		for ( i = 0; i<3; i++) BField[i] += extra_brf[i];
	}
	
  // Classical Description of Spin 
  // NB moment has "I=1/2" already... 
	dxvsdt[6] = 0.0;
	dxvsdt[7] = btoomega*(xvs[8]*BField[2] - xvs[9]*BField[1]);
	dxvsdt[8] = btoomega*(xvs[9]*BField[0] - xvs[7]*BField[2]);
	dxvsdt[9] = btoomega*(xvs[7]*BField[1] - xvs[8]*BField[0]);
	// ihd/dt = (moment)*(sigma)*B 
  // NB moment has "I=1/2" already... 
  // dxvsdt[6] = btoomega * ( (BField[2] * xvs[7]) + (BField[0] * xvs[9]) - (BField[1] * xvs[8]));
  // dxvsdt[7] = btoomega * (-(BField[2] * xvs[6]) - (BField[0] * xvs[8]) - (BField[1] * xvs[9]));
  // dxvsdt[8] = btoomega * (-(BField[2] * xvs[9]) + (BField[0] * xvs[7]) + (BField[1] * xvs[6]));
  // dxvsdt[9] = btoomega * ( (BField[2] * xvs[8]) - (BField[0] * xvs[6]) + (BField[1] * xvs[7]));
  return 0;
}
__device__ /*__host__*/ int CUDA_rkck_XVS( 
	double *d_IO, 
	int *d_IO_INT, 
	int vi_RECORD,
	double t,
	double xvs[],
	double dxvsdt[],
	double h,
	double xvserr[],
	double xvsout[],
	double BField[], 
	double dField_1D_FLAT[])
{
	// double BField[3], dField_1D_FLAT[9];
  CUDA_derivs_XVS(t, xvs, dxvsdt, BField, dField_1D_FLAT);
  int i;
  /* 	static const a2=(0.2, a3=(0.3, a4=(0.6, a5=(1.0, a6=(0.875; */
  const double b21=(0.2), b31=(3.0/40.0), b32=(9.0/40.0), b41=(0.3), b42=( -0.9), b43=(1.2);
  const double b51=( -11.0/54.0), b52=(2.5), b53=( -70.0/27.0), b54=(35.0/27.0); 
  const double b61=(1631.0/55296.0), b62=(175.0/512.0), b63=(575.0/13824.0), b64=(44275.0/110592.0), b65=(253.0/4096.0); 
  const double c1=(37.0/378.0), c3=(250.0/621.0), c4=(125.0/594.0), c6=(512.0/1771.0);
  const double dc5=( -277.0/14336.0), dc1=(c1-2825.0/27648.0),  dc3=(c3-18575.0/48384.0), dc4=(c4-13525.0/55296.0), dc6=(c6-0.25);
  double ak2[10], ak3[10], ak4[10], ak5[10], ak6[10], xvstemp[10];
  for (i = 0; i<10; i++)	 /* First step */ xvstemp[i]=xvs[i]+b21*h*dxvsdt[i];
  CUDA_derivs_XVS(t, xvs, ak2, BField, dField_1D_FLAT);		/* Second step */
  for (i = 0; i<10; i++) xvstemp[i]=xvs[i]+h*(b31*dxvsdt[i]+b32*ak2[i]);
  CUDA_derivs_XVS(t, xvs, ak3, BField, dField_1D_FLAT);		/* Third step */
  for (i = 0; i<10; i++) xvstemp[i]=xvs[i]+h*(b41*dxvsdt[i]+b42*ak2[i]+b43*ak3[i]);
  CUDA_derivs_XVS(t, xvs, ak4, BField, dField_1D_FLAT);		/* Fourth step */
  for (i = 0; i<10; i++) xvstemp[i]=xvs[i]+h*(b51*dxvsdt[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
  CUDA_derivs_XVS(t, xvstemp, ak5, BField, dField_1D_FLAT);		/* Fifth step */
  for (i = 0; i<10; i++) xvstemp[i]=xvs[i]+h*(b61*dxvsdt[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
  CUDA_derivs_XVS(t, xvs, ak6, BField, dField_1D_FLAT);		/* Sixth step */
  for (i = 0; i<10; i++)	xvsout[i]=xvs[i]+h*(c1*dxvsdt[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]); /* Accumulate increments with proper weights */
  for (i = 0; i<10; i++) xvserr[i]=h*(dc1*dxvsdt[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
  return 0;
}
__device__ /*__host__*/ int CUDA_rkqs_SINGLE_ATTEMPT_XVS(
	double *d_IO, 
	int *d_IO_INT,
	int vi_RECORD, 
	double xvs[],
 	double dxvsdt[],
 	double *t,
 	double htry,
	double *hdid,
 	double *hnext,
 	double xvs_scal[],
 	double BField[],
	double dField_1D_FLAT[],
 	double epsilon[],
 	int *rkqs_TRIED)
{
	int i, j, k;
	double hnext_xv;
	double hnext_s;
  double hcurrent = htry;
	double xvs_temp[10];
  double epsilon_temp[10];
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
		xvs, 
		dxvsdt, 
		hcurrent, 
		epsilon_temp, 
		xvs_temp, 
		BField, 
		dField_1D_FLAT);
	double epsilon_max_xv = 0.0;
	double epsilon_max_s = 0.0;
	for (i = 0; i<10; i++)
	{
		if (i<6)
		{
			abs_max_xv = fabs(epsilon_temp[i]/xvs_scal[i]);
			// if (abs_max_xv<0) abs_max_xv *= -1.0; 
			if (abs_max_xv>epsilon_max_xv) epsilon_max_xv = abs_max_xv;
			
		}
		else
		{
			abs_max_s = fabs(epsilon_temp[i]/xvs_scal[i]);
			// if (abs_max_s<0) abs_max_s *= -1.0; 
			if (abs_max_s>epsilon_max_s) epsilon_max_s = abs_max_s;
		}
	}
	error_xv = epsilon_max_xv/d_CONST[e_d_CONST_def_MAXERR];
	error_s = epsilon_max_s/d_CONST[e_d_CONST_def_MAXERR1];
	
	/////////////////////////////////////////////////////////////////////////
	/////////////////// XV AND SPIN STEP COMPUTATION ////////////////////////
	/////////////////////////////////////////////////////////////////////////
	
	/* XV: Step succeeded. Compute size of next step */
	hnext_xv_GU = d_CONST[e_d_CONST_def_SAFETY]*hcurrent*pow(error_xv,d_CONST[e_d_CONST_def_PGROW]);
	/* XV: No more than a factor of 5 increase */
	hnext_xv_GB = 5.0*hcurrent; // ((hnext_xv_GU >= 0.0) ? fmin(hnext_xv_GU, 5.0*hcurrent) : fmax(hnext_xv_GU, 5.0*hcurrent));
	/* XV: Truncation error too large. Reduce stepsize */
	hnext_xv_SU = d_CONST[e_d_CONST_def_SAFETY]*hcurrent*pow(error_xv, d_CONST[e_d_CONST_def_PSHRNK]);
	/* XV: No more than a factor of 10 */
	hnext_xv_SB = 0.1*hcurrent; // ((hnext_xv_SU >= 0.0) ? fmax(hnext_xv_SU, 0.1*hcurrent) : fmin(hnext_xv_SU, 0.1*hcurrent));
	/* S: Step succeeded. Compute size of next step */
	hnext_s_GU = d_CONST[e_d_CONST_def_SAFETY1]*hcurrent*pow(error_s,d_CONST[e_d_CONST_def_PGROW1]);
	/* S: No more than a factor of 5 increase */
	hnext_s_GB = 5.0*hcurrent; // ((hnext_s_GU >= 0.0) ? fmin(hnext_s_GU, 5.0*hcurrent) : fmax(hnext_s_GU, 5.0*hcurrent));
	/* S: Truncation error too large. Reduce stepsize */
	hnext_s_SU = d_CONST[e_d_CONST_def_SAFETY1]*hcurrent*pow(error_s, d_CONST[e_d_CONST_def_PSHRNK1]);
	/* S: No more than a factor of 10 */
	hnext_s_SB = 0.1*hcurrent; // ((hcurrent >= 0.0) ? fmax(hnext_s_SU, 0.1*hcurrent) : fmin(hnext_s_SU, 0.1*hcurrent));
	
	/////////////////////////////////////////////////////////////////////////
	/////////////////// XV ERROR AND STEP DECISION //////////////////////////
	/////////////////////////////////////////////////////////////////////////
	
	if (error_xv <= 1.0)
	{
		hnext_xv = ((hcurrent >= 0.0) ? fmin(hnext_xv_GU, hnext_xv_GB) : fmax(hnext_xv_GU, hnext_xv_GB));
    return_value_xv = 0;
  }
	else if (error_xv > 1.0)
	{
		hnext_xv = ((hcurrent >= 0.0) ? fmax(hnext_xv_SU, hnext_xv_SB) : fmin(hnext_xv_SU, hnext_xv_SB));
		return_value_xv = 1;
	}
	
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
		for (i=0; i<10; i++)
		{
			xvs[i] = xvs_temp[i];
			epsilon[i] = epsilon_temp[i];
		}
		*t += hcurrent;
		*hdid = hcurrent;
		*rkqs_TRIED = 0;
		if (hcurrent<=0.0)
		{
			*hnext = fmax ( fmax (hnext_xv, hnext_s) , fmax ( d_CONST[e_d_CONST_tframe], d_CONST[e_d_CONST_tframe_SPIN]));
		}
		else if (hcurrent>0.0) 
		{
			*hnext = fmin ( fmin (hnext_xv, hnext_s) , fmin ( d_CONST[e_d_CONST_tframe], d_CONST[e_d_CONST_tframe_SPIN]));
		}
		return_value = e_RKQS_ERROR_NONE;
	}
	else
	{
		*hdid = 0.0;
		(*rkqs_TRIED)++;
		if (return_value_xv==1 && return_value_s==0) 
		{
			*hnext = hnext_xv;
			return_value = e_RKQS_ERROR_XV_BOUNDS;
		}
		else if (return_value_xv==0 && return_value_s==1)
		{
			*hnext = hnext_s;
			return_value = e_RKQS_ERROR_SPIN_BOUNDS;
		}
		else if (return_value_xv==1 && return_value_s==1)
		{
			if ((hnext_xv<=0.0 && hnext_s<=0.0) || (hnext_xv>0.0 && hnext_s>0.0>=0.0))
			{
				return_value = e_RKQS_ERROR_COMBINED_BOUNDS;
				if (hnext_xv>0.0 && hnext_s>0.0) *hnext = fmin (hnext_xv, hnext_s);
				else if (hnext_xv<=0.0 && hnext_s<=0.0) *hnext = fmin (hnext_xv, hnext_s);
			}
			else return_value = e_RKQS_ERROR_REVERSED_INTERVAL;
			
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
	double l_time_CURRENT, double l_xvs[], double l_epsilon[], double l_xvs_scal[], double l_dxvsdt[], 
	double l_BField[], double l_dField_1D_FLAT[], double l_pol, int l_rkqs_TRIED, double l_hnext)
{
	int vi_DOUBLE_IO_OFFSET = (getGlobalIdx_3D_3D()*d_CONST_INT[e_d_CONST_INT_numRecordsPerThread] + (*p_vi_RECORD))*e_d_IO_LAST;
	int vi_INT_IO_OFFSET = (getGlobalIdx_3D_3D()*d_CONST_INT[e_d_CONST_INT_numRecordsPerThread] + (*p_vi_RECORD))*e_d_IO_INT_LAST;
	
	d_IO_INT[vi_INT_IO_OFFSET + e_d_IO_INT_ERROR] = 0;
	d_IO_INT[vi_INT_IO_OFFSET + e_d_IO_INT_THREAD] = getGlobalIdx_3D_3D();
	d_IO_INT[vi_INT_IO_OFFSET + e_d_IO_INT_RECORD] = (*p_vi_RECORD);
  d_IO_INT[vi_INT_IO_OFFSET + e_d_IO_INT_RKQS_STEPS] = l_rkqs_TRIED; 
	
	d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_T] = l_time_CURRENT;
	d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_X] = l_xvs[0];
	d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_Y] = l_xvs[1];
	d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_Z] = l_xvs[2];
	d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_VX] = l_xvs[3];
	d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_VY] = l_xvs[4];
	d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_VZ] = l_xvs[5];
	d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_SX] = CUDA_SPIN_XVS_SINGLE(l_xvs, 0);
	d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_SY] = CUDA_SPIN_XVS_SINGLE(l_xvs, 1);
	d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_SZ] = CUDA_SPIN_XVS_SINGLE(l_xvs, 2);
	d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_ERR_X] = l_epsilon[0];
	d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_ERR_Y] = l_epsilon[1];
	d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_ERR_Z] = l_epsilon[2];
	d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_ERR_VX] = l_epsilon[3];
	d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_ERR_VY] = l_epsilon[4];
	d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_ERR_VZ] = l_epsilon[5];
	d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_SCAL_X] = l_xvs_scal[0];
	d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_SCAL_Y] = l_xvs_scal[1];
	d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_SCAL_Z] = l_xvs_scal[2];
	d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_SCAL_VX] = l_xvs_scal[3];
	d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_SCAL_VY] = l_xvs_scal[4];
	d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_SCAL_VZ] = l_xvs_scal[5];
	d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_RED_VX] = l_dxvsdt[0];
	d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_RED_VY] = l_dxvsdt[1];
	d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_RED_VZ] = l_dxvsdt[2];  
	d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_AX] = l_dxvsdt[3];
	d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_AY] = l_dxvsdt[4];
	d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_AZ] = l_dxvsdt[5];
	d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_BX] = l_BField[0];
	d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_BY] = l_BField[1];
	d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_BZ] = l_BField[2];
	d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_dB_XDX] = l_dField_1D_FLAT[0];
	d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_dB_XDY] = l_dField_1D_FLAT[1];
	d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_dB_XDZ] = l_dField_1D_FLAT[2];
	d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_dB_YDX] = l_dField_1D_FLAT[3];
	d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_dB_YDY] = l_dField_1D_FLAT[4];
	d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_dB_YDZ] = l_dField_1D_FLAT[5];
	d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_dB_ZDX] = l_dField_1D_FLAT[6];
	d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_dB_ZDY] = l_dField_1D_FLAT[7];
	d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_dB_ZDZ] = l_dField_1D_FLAT[8];
	d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_SPINNOR_0] = l_xvs[6];
	d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_SPINNOR_1] = l_xvs[7];
	d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_SPINNOR_2] = l_xvs[8];
	d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_SPINNOR_3] = l_xvs[9];
	d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_DDT_SPINNOR_0] = l_dxvsdt[6]; 
	d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_DDT_SPINNOR_1] = l_dxvsdt[7]; 
  d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_DDT_SPINNOR_2] = l_dxvsdt[8]; 
  d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_DDT_SPINNOR_3] = l_dxvsdt[9]; 
	d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_ERR_SPINNOR_0] = l_epsilon[6]; 
  d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_ERR_SPINNOR_1] = l_epsilon[7]; 
  d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_ERR_SPINNOR_2] = l_epsilon[8]; 
  d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_ERR_SPINNOR_3] = l_epsilon[9]; 	
	d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_SCAL_SPINNOR_0] = l_xvs_scal[6]; 
  d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_SCAL_SPINNOR_1] = l_xvs_scal[7]; 
  d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_SCAL_SPINNOR_2] = l_xvs_scal[8]; 
  d_IO[vi_DOUBLE_IO_OFFSET + e_d_IO_SCAL_SPINNOR_3] = l_xvs_scal[9];
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
	int vi_RECORD_TEST_0, vi_RECORD_TEST_1, vi_RECORD_TEST_2, vi_RECORD_TEST_3,	vi_RECORD_TEST_4;
	int vi_RECORD_IO_OFFSET_END = 0, vi_IO_OFFSET_END, vi_IO_INT_OFFSET_END, vi_IO_DOUBLE_OFFSET_END;
	int vi_TESTOTESTO;
	int i, j, k, vi_RECORD, vi_RKQS_STEP, return_value_RKQS, l_rkqs_TRIED;
	int l_odeint_steps, vi_BREAK_FLAG, vi_REVERSE_FLAG, vi_INDEX;
	double l_time_CURRENT;
	double l_xvs[10];
	double l_dxvsdt[10];
	double l_xvs_scal[10];
	double l_epsilon[10];
	double l_BField[3];
	double l_dField_1D_FLAT[9];
	double l_pol;
	double l_hnext, l_hdid, l_htry;
	int vi_RECORD_IO_OFFSET_START = getGlobalIdx_3D_3D()*d_CONST_INT[e_d_CONST_INT_numRecordsPerThread] + numRecordsStart;
	int vi_IO_OFFSET_START = vi_RECORD_IO_OFFSET_START*e_d_IO_LAST;
	int vi_IO_INT_OFFSET_START = vi_RECORD_IO_OFFSET_START*e_d_IO_INT_LAST;
	int vi_IO_DOUBLE_OFFSET_START = vi_RECORD_IO_OFFSET_START*e_d_IO_LAST;
	
	return_value_RKQS = d_IO[vi_IO_INT_OFFSET_START + e_d_IO_INT_ERROR];  
	l_rkqs_TRIED = 0;
	l_time_CURRENT = d_IO[vi_IO_DOUBLE_OFFSET_START + e_d_IO_T]; 
	for (i = 0 ; i<10; i++)
	{
		l_xvs[i] = d_IO[vi_IO_DOUBLE_OFFSET_START + e_d_IO_X + i];
		l_dxvsdt[i] = d_IO[vi_IO_DOUBLE_OFFSET_START + e_d_IO_RED_VX + i];
		l_xvs_scal[i] = d_IO[vi_IO_DOUBLE_OFFSET_START + e_d_IO_SCAL_X + i];
		l_epsilon[i] = d_IO[vi_IO_DOUBLE_OFFSET_START + e_d_IO_ERR_X + i];
	}
	for (i = 0 ; i<3; i++)
	{
		l_BField[i] = d_IO[vi_IO_DOUBLE_OFFSET_START + e_d_IO_BX + i]; 
	}
	for (i = 0 ; i<9; i++)
	{
		l_dField_1D_FLAT[i] = d_IO[vi_IO_DOUBLE_OFFSET_START + e_d_IO_dB_XDX + i]; 
	}
	l_pol = d_IO[vi_IO_DOUBLE_OFFSET_START + e_d_IO_POLARIZATION]; 
	l_hnext = d_IO[vi_IO_DOUBLE_OFFSET_START + e_d_IO_HNEXT]; 
	
	double vd_TIME_START = l_time_CURRENT;
	
	for (i = 0 ; i<3; i++)
	{
		l_BField[i] = d_IO[vi_IO_DOUBLE_OFFSET_START + e_d_IO_BX + i]; 
		
	}
	for (i = 0 ; i<9; i++)
	{
		l_dField_1D_FLAT[i] = d_IO[vi_IO_DOUBLE_OFFSET_START + e_d_IO_dB_XDX + i]; 
		
	}
	l_pol = d_IO[vi_IO_DOUBLE_OFFSET_START + e_d_IO_POLARIZATION]; 
	l_hnext = d_IO[vi_IO_DOUBLE_OFFSET_START + e_d_IO_HNEXT]; 
	
	// CUDA_setspin( l_xvs, l_spin, l_spinnor, l_BField);
	int vi_CYCLE = 0;
	CUDA_derivs_XVS(l_time_CURRENT, l_xvs, l_dxvsdt, l_BField, l_dField_1D_FLAT);
	l_pol = CUDA_polcalc_XVS(l_BField[0], l_BField[1], l_BField[2], l_xvs);
	vi_REVERSE_FLAG = 0;
	// if (d_CONST[e_d_CONST_h1]<d_CONST[e_d_CONST_h1_SPIN]) l_htry = d_CONST[e_d_CONST_h1];
	// else l_htry = fmin(d_CONST[e_d_CONST_h1],d_CONST[e_d_CONST_h1_SPIN]);
	vi_RKQS_STEP = 0;
	return_value_RKQS = 0;
	while(vi_RECORD<numRecordsEnd)
	{
		l_htry = l_hnext;
		l_xvs_scal[0]=fabs(l_xvs[0])+fabs(l_dxvsdt[0]*l_htry)+d_CONST[e_d_CONST_def_TINY];
		l_xvs_scal[1]=fabs(l_xvs[1])+fabs(l_dxvsdt[1]*l_htry)+d_CONST[e_d_CONST_def_TINY];
		l_xvs_scal[2]=fabs(l_xvs[2])+fabs(l_dxvsdt[2]*l_htry)+d_CONST[e_d_CONST_def_TINY];
		l_xvs_scal[3]=fabs(l_xvs[3])+fabs(l_dxvsdt[3]*l_htry)+d_CONST[e_d_CONST_def_TINY];
		l_xvs_scal[4]=fabs(l_xvs[4])+fabs(l_dxvsdt[4]*l_htry)+d_CONST[e_d_CONST_def_TINY];
		l_xvs_scal[5]=fabs(l_xvs[5])+fabs(l_dxvsdt[5]*l_htry)+d_CONST[e_d_CONST_def_TINY];
		l_xvs_scal[6]=fabs(l_xvs[6])+fabs(l_dxvsdt[6]*l_htry)+d_CONST[e_d_CONST_def_TINY1];
		l_xvs_scal[7]=fabs(l_xvs[7])+fabs(l_dxvsdt[7]*l_htry)+d_CONST[e_d_CONST_def_TINY1];
		l_xvs_scal[8]=fabs(l_xvs[8])+fabs(l_dxvsdt[8]*l_htry)+d_CONST[e_d_CONST_def_TINY1];
		l_xvs_scal[9]=fabs(l_xvs[9])+fabs(l_dxvsdt[9]*l_htry)+d_CONST[e_d_CONST_def_TINY1];
		
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
			CUDA_derivs_XVS(l_time_CURRENT, l_xvs, l_dxvsdt, l_BField, l_dField_1D_FLAT);
			l_pol = CUDA_polcalc_XVS(l_BField[0], l_BField[1], l_BField[2], l_xvs);
			CUDA_RECORD_XVS(d_IO, d_IO_INT,  &vi_RECORD, 
				l_time_CURRENT, l_xvs, l_epsilon, l_xvs_scal, l_dxvsdt, 
				l_BField, l_dField_1D_FLAT, l_pol, l_rkqs_TRIED, l_hnext);
		}
		else
		{
			return_value_RKQS = CUDA_rkqs_SINGLE_ATTEMPT_XVS(
				d_IO, 
				d_IO_INT, 
				vi_RECORD, 
				l_xvs, 
				l_dxvsdt, 
				&l_time_CURRENT, 
				l_htry, 
				&l_hdid, 
				&l_hnext, 
				l_xvs_scal, 
				l_BField, 
				l_dField_1D_FLAT, 
				l_epsilon, 
				&l_rkqs_TRIED);
			if (return_value_RKQS==e_RKQS_ERROR_NONE)
			{
				l_time_CURRENT = l_time_CURRENT + l_hdid;
				CUDA_derivs_XVS(l_time_CURRENT, l_xvs, l_dxvsdt, l_BField, l_dField_1D_FLAT);
				l_pol = CUDA_polcalc_XVS(l_BField[0], l_BField[1], l_BField[2], l_xvs);
				// CUDA_RECORD_DOUBLE(d_IO, vi_RECORD, e_d_IO_EXTRA_1, l_hnext);
				vi_RKQS_STEP++;
				if (vi_RKQS_STEP>=d_CONST_INT[e_d_CONST_INT_numCyclesPerRecord])
				{
					// CUDA_RECORD_DOUBLE(d_IO, vi_RECORD, e_d_IO_EXTRA_0, l_time_CURRENT);
					CUDA_RECORD_XVS(d_IO, d_IO_INT,  &vi_RECORD, 
						l_time_CURRENT, l_xvs, l_epsilon, l_xvs_scal, l_dxvsdt, 
						l_BField, l_dField_1D_FLAT, l_pol, l_rkqs_TRIED, l_hnext);
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
	for (i = 0 ; i<10; i++)
	{
		d_IO[vi_IO_DOUBLE_OFFSET_END + e_d_IO_X + i] = l_xvs[i];
		d_IO[vi_IO_DOUBLE_OFFSET_END + e_d_IO_RED_VX + i] = 	l_dxvsdt[i];
		d_IO[vi_IO_DOUBLE_OFFSET_END + e_d_IO_SCAL_X + i] = 	l_xvs_scal[i];
		d_IO[vi_IO_DOUBLE_OFFSET_END + e_d_IO_ERR_X + i] = 	l_epsilon[i];
	}
	for (i = 0 ; i<3; i++)
	{
		d_IO[vi_IO_DOUBLE_OFFSET_END + e_d_IO_BX + i] = l_BField[i]; 
	}
	for (i = 0 ; i<9; i++)
	{
		d_IO[vi_IO_DOUBLE_OFFSET_END + e_d_IO_dB_XDX + i] = l_dField_1D_FLAT[i]; 
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
