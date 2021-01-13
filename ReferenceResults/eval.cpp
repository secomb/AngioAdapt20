/*****************************************************
Evaluate solute field from source strengths.  TWS November 07.
Version 1.0, May 1, 2008.
With 9/11/08 update.
Updated June 2016 for faster execution.
Note: result for isp=2 is not used, kept for consistency
******************************************************/
#include <math.h>
#include "nrutil.h"
#include <stdio.h>

float length(float* x);
float length0(float* x);

float* eval(int slsegdiv, float req, float* x)
{
	extern int nnt, nnv, nseg, nnod, nsp, hexperiodic, is2d;
	extern int* mainseg, * segtyp, * permsolute, ** tisspoints;
	extern float fac, w2d, r2d, aphexp;
	extern float* axt, * ayt, * azt, * ds, * diff, * g0, * y, * p, * P1, * midpt;
	extern float** qt, ** scos, ** qv, ** ax, ** hex_norm;

	float dist, gtt, gtv;
	int i, j, k, iseg, itp, isp, idomain, ndomain;

	if (hexperiodic) {
		for (i = 1; i <= 3; i++) y[i] = x[i] - midpt[i];
		length(y);		//remap x into base domain - gives periodic solute fields
		for (i = 1; i <= 3; i++) x[i] = y[i] + midpt[i];
	}

	for (isp = 1; isp <= nsp; isp++) p[isp] = g0[isp];		//initialize to g0
	for (itp = 1; itp <= nnt; itp++) {
		P1[1] = x[1] - axt[tisspoints[1][itp]];
		P1[2] = x[2] - ayt[tisspoints[2][itp]];
		P1[3] = x[3] - azt[tisspoints[3][itp]];
		if (hexperiodic) ndomain = 7;
		else ndomain = 1;
		gtt = 0.;
		for (idomain = 0; idomain < ndomain; idomain++) {
			for (k = 1; k <= 3; k++) {
				y[k] = P1[k];
				if (idomain > 0) y[k] += 2. * aphexp * hex_norm[idomain][k];
			}
			dist = length0(y);
			if (dist <= req) {
				if (is2d) gtt += fac / w2d * (2. * log(r2d / req) + 1. - SQR(dist / req));
				else gtt += fac * (1.5 - 0.5 * SQR(dist / req)) / req;
			}
			else {
				if (is2d) gtt += 2. * fac / w2d * log(r2d / dist);
				else gtt += fac / dist;
			}
		}
		for (isp = 1; isp <= nsp; isp++)	p[isp] += gtt / diff[isp] * qt[itp][isp];//add contributions from tissue sources
	}
	//add contributions from vessel sources. Subdivide subsegments. Vessel point is at midpoint of subsegment
	for (i = 1; i <= nnv; i++) {
		iseg = mainseg[i];
		for (k = 1; k <= slsegdiv; k++) {
			for (j = 1; j <= 3; j++)	P1[j] = ax[j][i] + scos[j][iseg] * ds[iseg] * (-0.5 + (k - 0.5) / slsegdiv) - x[j];
			if (hexperiodic) ndomain = 7;
			else ndomain = 1;
			gtv = 0.;
			for (idomain = 0; idomain < ndomain; idomain++) {
				for (k = 1; k <= 3; k++) {
					y[k] = P1[k];
					if (idomain > 0) y[k] += 2. * aphexp * hex_norm[idomain][k];
				}
				dist = length0(y);
				if (dist <= req) {
					if (is2d) gtv += fac / w2d * (2. * log(r2d / req) + 1. - SQR(dist / req));
					else gtv += fac * (1.5 - 0.5 * SQR(dist / req)) / req;
				}
				else {
					if (is2d) gtv += 2. * fac / w2d * log(r2d / dist);
					else gtv += fac / dist;
				}
			}
			for (isp = 1; isp <= nsp; isp++) if (permsolute[isp]) p[isp] += gtv / diff[isp] * qv[i][isp] / slsegdiv;
		}
	}
	return p;
}
