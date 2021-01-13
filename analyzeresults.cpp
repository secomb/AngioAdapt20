/************************************************************************
analyzeresults for AngioAdapt
TWS, June 2020
*************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void histogram(float* var, float* weight, int n, const char filename[]);

void analyzeresults()
{
	extern int nnod, nseg, nnt, nnv, phaseseparation, * ista, * iend, * mainseg;
	extern int* nodtyp, * segtyp, ** nodseg;
	extern float pi1;
	extern float* q, * qq, * diam, * lseg, * tau, * segpress, * hd, * dtmin, * sourcefac, * ds;
	extern float** pt, ** pvseg;
	extern double* nodpress, * cond;

	int i, iseg, inod, maxdim;
	float* histogramdisplay, * histogramweight, * length_weight;

	maxdim = IMAX(nnod, nseg);
	maxdim = IMAX(maxdim, nnv);
	maxdim = IMAX(maxdim, nnt);
	histogramdisplay = vector(1, maxdim);
	histogramweight = vector(1, maxdim);
	length_weight = vector(1, nnod);

	//node weighting factors: 1/2 the sum of the lengths of the segments connected to the node
	for (inod = 1; inod <= nnod; inod++) {
		length_weight[inod] = 0;
		for (i = 1; i <= nodtyp[inod]; i++) {
			iseg = nodseg[i][inod];
			if (segtyp[iseg] >= 3 && segtyp[iseg] <= 5) length_weight[inod] += 0.5 * lseg[iseg];
		}
	}

	//Node pressures
	i = 0;
	for (inod = 1; inod <= nnod; inod++) if (length_weight[inod] > 0.) {
		i++;
		histogramdisplay[i] = nodpress[inod];
		histogramweight[i] = length_weight[inod];
	}
	histogram(histogramdisplay, histogramweight, i, "Histogram-pressures.out");

	//Segment flows
	i = 0;
	for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] >= 3 && segtyp[iseg] <= 5 && qq[iseg] > 0.) {
		i++;
		histogramdisplay[i] = log10(qq[iseg]);
		histogramweight[i] = lseg[iseg];
	}
	histogram(histogramdisplay, histogramweight, i, "Histogram-logflows.out");

	//Sement velocities in mm/s
	i = 0;
	for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] >= 3 && segtyp[iseg] <= 5 && qq[iseg] > 0.) {
		i++;
		histogramdisplay[i] = log10(4000. / 60. * qq[iseg] / pi1 / SQR(diam[iseg]));
		histogramweight[i] = lseg[iseg];
	}
	histogram(histogramdisplay, histogramweight, i, "Histogram-logvelocities.out");

	//Shear stress
	i = 0;
	for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] >= 3 && segtyp[iseg] <= 5 && fabs(tau[iseg]) > 0.) {
		i++;
		histogramdisplay[i] = log10(fabs(tau[iseg]));
		histogramweight[i] = lseg[iseg];
	}
	histogram(histogramdisplay, histogramweight, i, "Histogram-logstress.out");

	//vessel PO2
	i = 0;
	for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] >= 3 && segtyp[iseg] <= 5) {
		i++;
		histogramdisplay[i] = pvseg[iseg][1];
		histogramweight[i] = lseg[iseg];
	}
	histogram(histogramdisplay, histogramweight, i, "HistogramPO2v.out");

	//tissue PO2
	for (i = 1; i <= nnt; i++) {
		histogramdisplay[i] = FMIN(pt[i][1], 99.9);
		histogramweight[i] = sourcefac[i];
	}
	histogram(histogramdisplay, histogramweight, nnt, "HistogramPO2t.out");

	//distance to nearest vessel
	for (i = 1; i <= nnt; i++) {
		histogramdisplay[i] = dtmin[i];
		histogramweight[i] = sourcefac[i];
	}
	histogram(histogramdisplay, histogramweight, nnt, "HistogramDistance.out");

	//hematocrit
	if (phaseseparation) {
		i = 0;
		for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] >= 3 && segtyp[iseg] <= 5) {
			i++;
			histogramdisplay[i] = hd[iseg];
			histogramweight[i] = lseg[iseg];
		}
		histogram(histogramdisplay, histogramweight, i, "Histogram-hematocrit.out");
	}

	free_vector(histogramdisplay, 1, maxdim);
	free_vector(histogramweight, 1, maxdim);
	free_vector(length_weight, 1, nnod);
}