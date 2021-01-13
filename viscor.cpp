/************************************************************************
viscor - for AngioAdapt07
Computation of segment viscosity based on Pries et al. Circ Res. 75,
904-915, 1994. Diameter corrected for differences in MCV
TWS October 07
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

float viscor(float d, float h)
{
	extern float mcvcorr, diamthresh, optw, vplas;
	extern float* cpar, * viscpar;

	float dcorr, c, eta45, etarel, hdfac, hdref = 0.45F;
	dcorr = d * mcvcorr;
	if (dcorr < diamthresh) dcorr = diamthresh;
	c = (cpar[1] + exp(cpar[2] * dcorr)) * (-1. + 1. / (1. + pow(10., cpar[3]) * pow(dcorr, cpar[4])))
		+ 1. / (1. + pow(10., cpar[3]) * pow(dcorr, cpar[4]));
	eta45 = viscpar[1] * exp(viscpar[2] * dcorr) + viscpar[3] + viscpar[4] * exp(viscpar[5] * pow(dcorr, viscpar[6]));
	hdfac = (pow(1.F - h, c) - 1.) / (pow(1. - hdref, c) - 1.);
	etarel = (1. + (eta45 - 1.) * hdfac * SQR(dcorr / (dcorr - optw))) * SQR(dcorr / (dcorr - optw));
	return etarel * vplas;
}
