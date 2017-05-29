
/* Include files ----------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <ctime>
#include <iomanip>
#include <cmath>
#include <limits>
#include <string>

#include <cuda.h> 
#include <cuda_runtime.h>
#include "UCN_CUDA_ALL_KERNEL.cuh"
#include "UCN_CUDA_WRAPPER.h"

int GENERIC_WRAPPER_MULTI(
	const double *h_CONST,
	const int *h_CONST_INT,
  double *h_IO,
  int *h_IO_INT)
{
	GENERIC_MIDDLEMAN_MULTI(
		h_CONST,
		h_CONST_INT,
		h_IO,
		h_IO_INT);
	return 0;
}