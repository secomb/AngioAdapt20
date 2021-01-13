/*****************************************************
Tissue uptake rates of solutes as a function of levels
Must be provided for each application.  TWS November 07.
Revised TWS 2010
Version with two growth factors - December 2019
******************************************************/
#include <math.h>
#include "nrutil.h"
#include <stdio.h>

void tissrate(int nsp, float* p, float* mtiss, float* mptiss)
{
	extern float** tissparam;

	int isp;
	float m0 = tissparam[1][1];
	float pcr = tissparam[2][1];
	float pgf = tissparam[3][1];
	float gfrate2 = tissparam[1][2];
	float gfrate3 = tissparam[1][3];
	float ngf2 = tissparam[3][2];
	float ngf3 = tissparam[3][3];

	for (isp = 1; isp <= nsp; isp++) {
		switch (isp) {
		case 1:
			if (p[isp] >= 0.) {
				mtiss[1] = -m0 * p[1] / (p[1] + pcr);
				mptiss[1] = -m0 * pcr / SQR(p[1] + pcr);
			}
			else {
				mtiss[1] = -m0 * p[1] / pcr;
				mptiss[1] = -m0 / pcr;
			}
			break;
		case 2:	//note - linear uptake not included - used in Green's function calculation
			if (p[1] >= 0.) mtiss[2] = gfrate2 / (1. + powf(p[1] / pgf, ngf2));	//Hill-type model, March 2018
			else mtiss[2] = gfrate2;
			mptiss[2] = 0.;
			break;
		case 3:	//note - linear uptake used in Green's function calculation
			mtiss[3] = 0.;
			mptiss[3] = 0.;
			break;
		default:
			printf("*** Error: isp is too large in tissrate\n");
		}
	}
}
