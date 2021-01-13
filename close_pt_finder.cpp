/*****************************************************************************************
close_point_finder.cpp
This portion takes a point and finds the closest point on all other segments. If the point is
within tolerance_1, the segment number and closest point are returned.
tolerance_1: if closest point on another seg within this distance=>intercept at point
JA, conceived May 31 2004
Updated TWS, February 2018
**********************************************************************************************/

#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

float distance(float* x, float* y);
float length(float* x);

int close_pt_finder(float* new_pt, float* int_pt, int inod_from, int iseg_from)
{
	extern int nseg, nnod, hexperiodic;
	extern int* segtyp, * ista, * iend, * nodtyp, * boundseg;
	extern int** nodseg, ** nodnod;
	extern float tolerance_1, aphexp;
	extern float* lseg, * P1, * P2, * P3, ** cnode, ** hex_norm;

	int i, k, iseg, intercept_flag;
	float tprime, dist, dist_prime;

	intercept_flag = 0;
	dist = 1.e6;

	for (iseg = 1; iseg <= nseg; iseg++) {
		if (segtyp[iseg] == 3 || segtyp[iseg] == 4 || segtyp[iseg] == 10) goto skipiseg;
		if (inod_from > 0) for (k = 1; k <= nodtyp[inod_from]; k++) {	//exclude connected segments and segments connected to those
			if (nodseg[k][inod_from] == iseg) goto skipiseg;	//connected
			for (i = 1; i <= nodtyp[nodnod[k][inod_from]]; i++)
				if (nodseg[i][nodnod[k][inod_from]] == iseg) goto skipiseg;	//one segment away
		}
		if (iseg_from > 0) {		//exclude segment and connected segments
			if (iseg_from == iseg) goto skipiseg;
			for (i = 1; i <= nodtyp[ista[iseg_from]]; i++) if (nodseg[i][ista[iseg_from]] == iseg) goto skipiseg;
			for (i = 1; i <= nodtyp[iend[iseg_from]]; i++) if (nodseg[i][iend[iseg_from]] == iseg) goto skipiseg;
		}
		for (i = 1; i <= 3; i++) {
			P1[i] = cnode[i][iend[iseg]] - cnode[i][ista[iseg]];
			P2[i] = new_pt[i] - cnode[i][ista[iseg]];
		}
		if (hexperiodic) {
			length(P1);		//remap into base domain
			length(P2);
		}
		tprime = 0.;
		for (i = 1; i <= 3; i++) tprime += P1[i] * P2[i] / SQR(lseg[iseg]);
		if (tprime < 0.) tprime = 0.;
		if (tprime > 1.) tprime = 1.;
		for (i = 1; i <= 3; i++) P3[i] = cnode[i][ista[iseg]] + tprime * P1[i];	//closest point on segment to new_pt
		dist_prime = distance(P3, new_pt);
		if (dist_prime < dist && dist_prime <= tolerance_1) {
			intercept_flag = iseg;
			dist = dist_prime;
			for (i = 1; i <= 3; i++) int_pt[i] = P3[i];
		}
	skipiseg:;
	}
	return intercept_flag;
}
