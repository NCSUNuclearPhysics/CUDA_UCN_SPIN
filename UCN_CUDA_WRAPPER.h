#include "UCN_ENUM.h"

// int MAIN_CONST_INTERFACE(int param_d_CONST_INT[e_d_CONST_INT_LAST], double param_d_CONST[e_d_CONST_LAST]);
// 
// double MAIN_CUDA_theta(int param_d_CONST_INT[e_d_CONST_INT_LAST], double param_d_CONST[e_d_CONST_LAST], double x, double y, double z);
// 
// double MAIN_CUDA_ro(int param_d_CONST_INT[e_d_CONST_INT_LAST], double param_d_CONST[e_d_CONST_LAST], double x, double y);
// 
// double MAIN_CUDA_r(int param_d_CONST_INT[e_d_CONST_INT_LAST], double param_d_CONST[e_d_CONST_LAST], double x, double y, double z);
// 
// double MAIN_CUDA_phi(int param_d_CONST_INT[e_d_CONST_INT_LAST], double param_d_CONST[e_d_CONST_LAST], double x, double y);
// 
// int MAIN_CUDA_Mag_HALBACH_SURFACE(int param_d_CONST_INT[e_d_CONST_INT_LAST], double param_d_CONST[e_d_CONST_LAST], double xv[], double BField[]);
//  
// int MAIN_CUDA_Mag(int param_d_CONST_INT[e_d_CONST_INT_LAST], double param_d_CONST[e_d_CONST_LAST], double xv[], double BField[]);
// 
// int MAIN_CUDA_dB_HALBACH_SURFACE(int param_d_CONST_INT[e_d_CONST_INT_LAST], double param_d_CONST[e_d_CONST_LAST], double xv[], double BField[], double dField_1D_FLAT[]);
// 
// int MAIN_CUDA_dB(int param_d_CONST_INT[e_d_CONST_INT_LAST], double param_d_CONST[e_d_CONST_LAST], double xv[], double BField[], double dField_1D_FLAT[]);
// 
// int MAIN_CUDA_derivs_6_BACKUP(int param_d_CONST_INT[e_d_CONST_INT_LAST], double param_d_CONST[e_d_CONST_LAST], double xv[], double dxvdt[], 
//   double spin[], double BField[], double dField_1D_FLAT[]);
// 
// int MAIN_CUDA_derivs_6(int param_d_CONST_INT[e_d_CONST_INT_LAST], double param_d_CONST[e_d_CONST_LAST], double y[], double dydx[], 
//   double spin[], double BField[], double dField_1D_FLAT[]);
// 
// int MAIN_CUDA_rkck_6(int param_d_CONST_INT[e_d_CONST_INT_LAST], double param_d_CONST[e_d_CONST_LAST], double y[], double dydx[], double *h, double yerr[], double yout[], 
//   double spin[], double BField[], double dField_1D_FLAT[]);
// 
// int MAIN_CUDA_rkqs_6(int param_d_CONST_INT[e_d_CONST_INT_LAST], double param_d_CONST[e_d_CONST_LAST], double y[], double dydt[], 
//   double *x, double htry, double eps, double *hdid, double *hnext, 
//   double yscal[], double spin[], double BField[], double dField_1D_FLAT[], double yerr[]);
// 
// int MAIN_CUDA_odeint_6(int param_d_CONST_INT[e_d_CONST_INT_LAST], double param_d_CONST[e_d_CONST_LAST], double ystart[], double spin[], 
//   double x1, double x2, double eps, double epsilon[], double h1, double hmin, 
//   double BField[], double dField_1D_FLAT[]);

int GENERIC_WRAPPER_MULTI(
  const double *h_CONST,
  const int *h_CONST_INT,
  double *h_IO,
  int *h_IO_INT);
// void GENERIC_WRAPPER_BFIELD(
  // const double *h_CONST,
  // const int *h_CONST_INT,
  // double *h_IO);