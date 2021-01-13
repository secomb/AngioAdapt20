/**************************************************************************
bicgstabBLASDinit - double precision
initialize bicgstabBLASD
TWS, March 2011
Cuda 10.1 Version, August 2019
**************************************************************************/
//#include <shrUtils.h>
//#include <cutil_inline.h>
#include <cuda_runtime.h>
#include <cublas.h>
#include "nrutil.h"

void bicgstabBLASDinit(int n)
{
	extern int useGPU;
	extern double *h_x, *h_b, *h_a, *h_rs;
	extern double *d_a, *d_res, *d_x, *d_b;
	extern double *d_r, *d_rs, *d_v, *d_s, *d_t, *d_p, *d_er;
 	const int nmem0 = sizeof(double);	//needed for malloc
 	const int nmem1 = n*sizeof(double);	//needed for malloc
	const int nmem2 = n*n*sizeof(double);

	h_a = dvector(0,n*n-1);
	h_x = dvector(0,n-1);
	h_b = dvector(0,n-1);
	h_rs = dvector(0,n-1);
	
	cudaSetDevice( useGPU-1 );
	cublasInit(); 

	cudaMalloc((void **)&d_a, nmem2);
	cudaMalloc((void **)&d_r, nmem1);
	cudaMalloc((void **)&d_rs, nmem1);
	cudaMalloc((void **)&d_v, nmem1);
	cudaMalloc((void **)&d_s, nmem1);
	cudaMalloc((void **)&d_t, nmem1);
	cudaMalloc((void **)&d_p, nmem1);
	cudaMalloc((void **)&d_er, nmem1);
	cudaMalloc((void **)&d_x, nmem1);
	cudaMalloc((void **)&d_b, nmem1);
	cudaMalloc((void **)&d_res, nmem0);
}
