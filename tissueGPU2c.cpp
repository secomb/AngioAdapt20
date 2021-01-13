/**************************************************************************
tissueGPU2.cpp
program to call tissueGPU2.cu on GPU
TWS, January 2012
Cuda 10.1 Version, August 2019
**************************************************************************/
//#include <shrUtils.h>
//#include <cutil_inline.h>
#include <cuda_runtime.h>

extern "C" void tissueGPU2(float* d_tissxyz, float* d_vessxyz, float* d_pv000, float* d_qt000, float* d_hex_norm, float aphexp,
	int nnt, int nnv, int is2d, float req, float r2d, int stepGPU, int hexperiodic);

void tissueGPU2c(int stepGPU)
{
	extern int nnt, nnv, is2d, hexperiodic;
	extern float* d_tissxyz, * d_vessxyz, * d_pv000, * d_qt000, * d_hex_norm, * qt000, * pv000, req, r2d, aphexp;
	cudaMemcpy(d_qt000, qt000, nnt * sizeof(float), cudaMemcpyHostToDevice);
	tissueGPU2(d_tissxyz, d_vessxyz, d_pv000, d_qt000, d_hex_norm, aphexp, nnt, nnv, is2d, req, r2d, stepGPU, hexperiodic);
	cudaMemcpy(pv000, d_pv000, nnv * sizeof(float), cudaMemcpyDeviceToHost);
}
