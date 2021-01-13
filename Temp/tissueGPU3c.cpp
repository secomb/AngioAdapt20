/**************************************************************************
tissueGPU3.cpp
program to call tissueGPU3.cu on GPU
TWS, January 2012
Cuda 10.1 Version, August 2019
**************************************************************************/
#include <shrUtils.h>
#include <cutil_inline.h>
//#include <cuda_runtime.h>

extern "C" void tissueGPU3(float *d_tissxyz, float *d_vessxyz, float *d_pt000, float *d_qv000, float *d_hex_norm, float aphexp,
		int nnt, int nnv, int is2d, float req, float r2d, int stepGPU, int hexperiodic);

void tissueGPU3c(int stepGPU)
{
	extern int nnt,nnv,is2d,hexperiodic;
	extern float *d_tissxyz,*d_vessxyz,*d_pt000,*d_qv000,*qv000,*pt000,*d_hex_norm,aphexp,req,r2d;
	cudaError_t error;
	error = cudaMemcpy(d_qv000, qv000, nnv*sizeof(float), cudaMemcpyHostToDevice);
	tissueGPU3(d_tissxyz,d_vessxyz,d_pt000,d_qv000,d_hex_norm,aphexp,nnt,nnv,is2d,req,r2d,stepGPU,hexperiodic);
	error = cudaMemcpy(pt000, d_pt000, nnt*sizeof(float), cudaMemcpyDeviceToHost);
}
