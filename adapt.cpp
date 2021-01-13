/************************************************************************
adapt.cpp
Diameter adaptation, drop out of segments.  TWS October 07
Modified to use 2001 adaptation model, December 2010
Edited TWS Feb. 2017, June 2020
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void adapt()
{
	extern int nseg, nnod, eliminate, inftrans, nsminus, nspruned, * segtyp, adaptd;
	extern float diamthresh, km1, ks1, tauref1, timestep1, eqtime;
	extern float* segpress, * tau, * adaptpar, * condup, * conddown, * diam, * stot, * ksseg;
	extern float* stau, * spress, * uptrans, * downtrans, ** signalvar;

	int iseg;
	float tauexp;

	for (iseg = 1; iseg <= nseg; iseg++) {
		if (segtyp[iseg] == 3 || segtyp[iseg] == 5) {	//don't adapt type 4 segments
			stau[iseg] = log10(tau[iseg] + tauref1);
			tauexp = adaptpar[1] + adaptpar[2] / (1. + pow(segpress[iseg] / adaptpar[3], adaptpar[4]));	//TWS 2020
			spress[iseg] = -log10(tauexp);
			if (inftrans == 1) {	//downstream convected and upstream conducted responses, inverse diameter dependence, June 2020
				downtrans[iseg] = km1 * conddown[iseg] / diam[iseg];
				uptrans[iseg] = km1 * condup[iseg] / diam[iseg];
			}
			stot[iseg] = stau[iseg] + spress[iseg] + uptrans[iseg] + downtrans[iseg] - ksseg[iseg];
			if (adaptd == 1 || stot[iseg] > 0.) diam[iseg] += diam[iseg] * stot[iseg] * timestep1 / eqtime;	//Jan. 2020

			//test if diameters too small
			if (diam[iseg] < diamthresh) {
				if (eliminate) {
					printf(" %i", iseg);
					diam[iseg] = 0.;
					segtyp[iseg] = 10;
					nspruned++;
					nsminus++;
				}
				else diam[iseg] = diamthresh;	//use this option to inhibit pruning
			}
		}
		else {	//nonflowing and type 4 segments
			stau[iseg] = 0.;
			spress[iseg] = 0.;
			uptrans[iseg] = 0.;
			downtrans[iseg] = 0;
			stot[iseg] = 0.;
		}
		//save values for plots
		signalvar[iseg][1] = stau[iseg];	//shear
		signalvar[iseg][2] = spress[iseg];	//pressure
		signalvar[iseg][3] = -ksseg[iseg];	//shrinking tendency
		signalvar[iseg][4] = uptrans[iseg];	//upstream
		signalvar[iseg][5] = downtrans[iseg];//downstream
		signalvar[iseg][6] = stot[iseg];	//total
	}
}
