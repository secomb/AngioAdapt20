/************************************************************************
output - for AngioAdapt07
TWS December 07
Revised December 2010
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void picturenetwork(float* nodvar, float* segvar, const char fname[]);

void output()
{
	extern int nseg, nnod, * nodname, * segname;
	extern float* diam, * condup, * conddown, * stot, * metsig, * q, * hd, * nodvar, * segvar, * ksseg;
	int iseg, inod;

	for (iseg = 1; iseg <= nseg; iseg++) segvar[iseg] = iseg;
	for (inod = 1; inod <= nnod; inod++) nodvar[inod] = inod;
	picturenetwork(nodvar, segvar, "Network.ps");

	for (iseg = 1; iseg <= nseg; iseg++) segvar[iseg] = segname[iseg];
	for (inod = 1; inod <= nnod; inod++) nodvar[inod] = nodname[inod];
	picturenetwork(nodvar, segvar, "NetworkNames.ps");
}
