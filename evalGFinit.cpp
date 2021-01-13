/*****************************************************
evalGFinit - Initialize parameters needed for GF evaluation
This speeds up evalGF
TWS, March 2018
Two growth factors - December 2019
******************************************************/
#include <math.h>
#include "nrutil.h"
#include <stdio.h>

float bessi0(float x);
float bessk0(float x);
float bessi1(float x);
float bessk1(float x);

void evalGFinit() {
	extern int is2d;
	extern float pi1, w2d, req;
	extern float lambdaGF2, cfacGF2, ssGF2, bbGF2, ccGF2;
	extern float lambdaGF3, cfacGF3, ssGF3, bbGF3, ccGF3;
	extern float** tissparam, * diff;

	float DG, KG;

	DG = diff[2];
	KG = tissparam[2][2];
	lambdaGF2 = sqrt(KG / DG);
	if (is2d) {
		ssGF2 = 1. / (pi1 * req * req * w2d * KG);
		cfacGF2 = bessi0(lambdaGF2 * req) * bessk1(lambdaGF2 * req) + bessi1(lambdaGF2 * req) * bessk0(lambdaGF2 * req);
		bbGF2 = ssGF2 * bessi1(lambdaGF2 * req) / cfacGF2;
		ccGF2 = ssGF2 * bessk1(lambdaGF2 * req) / cfacGF2;
	}
	else {
		ssGF2 = 3. / (4. * pi1 * req * req * req * KG);
		bbGF2 = ssGF2 * req * (cosh(lambdaGF2 * req) - sinh(lambdaGF2 * req) / (lambdaGF2 * req));
		ccGF2 = ssGF2 * req * (1. + 1. / (lambdaGF2 * req)) * exp(-lambdaGF2 * req);
	}
	DG = diff[3];
	KG = tissparam[2][3];
	lambdaGF3 = sqrt(KG / DG);
	if (is2d) {
		ssGF3 = 1. / (pi1 * req * req * w2d * KG);
		cfacGF3 = bessi0(lambdaGF3 * req) * bessk1(lambdaGF3 * req) + bessi1(lambdaGF3 * req) * bessk0(lambdaGF3 * req);
		bbGF3 = ssGF3 * bessi1(lambdaGF3 * req) / cfacGF3;
		ccGF3 = ssGF3 * bessk1(lambdaGF3 * req) / cfacGF3;
	}
	else {
		ssGF3 = 3. / (4. * pi1 * req * req * req * KG);
		bbGF3 = ssGF3 * req * (cosh(lambdaGF3 * req) - sinh(lambdaGF3 * req) / (lambdaGF3 * req));
		ccGF3 = ssGF3 * req * (1. + 1. / (lambdaGF3 * req)) * exp(-lambdaGF3 * req);
	}
}
