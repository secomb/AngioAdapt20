/************************************************************************
conductconvect - for AngioAdapt10
Calculate downstream convected signal and upstream conducted response, with exponential decay.
condup and conddown represent total strength of signal.
Signal is generated in proportion to metpara*diam, and integrated along vessel.
Based on Pries et al 2001 version of adaptation model
Note: no dependence on diameters in upstream conduction: "summation or equal partition"
However, dependence on segment length is included to correct conceptual problem in original version
TWS November 2010
Version to allow option of diameter-weighted upstream conduction - see conductup1
**********************************************************************
Version of June 2020:
Signal = extravascular GF concentration
This is assumed to equal the influx rate of GF to the lumen per unit length.
(Both wall thickness and circumference are assumed proportional to diameter.)
conddown represents convective flux of GF, converted to wall response via Hill equation
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void putrank();
void evalGF(int isp, float req, float* x, float lambdaGF, float cfacGF, float ssGF, float bbGF, float ccGF,
	float* GF, float* GFx, float* GFy, float* GFz);
void conductup1();
void convectdown1();
float length(float* x);

void conductconvect()
{
	extern int nseg, hexperiodic;
	extern int* segtyp, * ista, * iend;
	extern float req, lambdaGF2, cfacGF2, ssGF2, bbGF2, ccGF2;
	extern float* metsig, * x, * midpt, ** cnode;

	int iseg, j;
	float GF, GFx, GFy, GFz;

	for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] >= 3 && segtyp[iseg] <= 5) {
		for (j = 1; j <= 3; j++) x[j] = (cnode[j][ista[iseg]] + cnode[j][iend[iseg]]) * 0.5 - midpt[j];
		if (hexperiodic) length(x);	//map into base domain
		for (j = 1; j <= 3; j++) x[j] += midpt[j];
		evalGF(2, req, x, lambdaGF2, cfacGF2, ssGF2, bbGF2, ccGF2, &GF, &GFx, &GFy, &GFz);
		metsig[iseg] = GF;	//simple version, June 2020
	}
	putrank();
	convectdown1();
	conductup1();
}

void convectdown1()
{
	extern int nnodfl, nseg;
	extern int* nodrank, * nodtyp, * nodout, * segtyp, ** nodseg;
	extern float Qref1, extlength1in, Sd50;
	extern float* conddown, * lseg, * qq, * metsig, ** pevseg;

	int inod, in, iseg, i;
	float signal, qsum;

	//all signals initially evaluated at upstream nodes
	for (in = 1; in <= nnodfl; in++) {
		inod = nodrank[in];
		//boundary inflow node: assume a further segment with same metsig and length extlength1in
		if (abs(nodtyp[inod]) == 1 && nodout[inod] == 1) {
			iseg = nodseg[1][inod];
			conddown[iseg] = metsig[iseg] * extlength1in;
		}
		if (nodtyp[inod] >= 2) {
			signal = 0.;	//input segments sum to give input of conddown
			for (i = nodout[inod] + 1; i <= nodtyp[inod]; i++) {	//inflow segments
				iseg = nodseg[i][inod];
				signal += conddown[iseg] + metsig[iseg] * lseg[iseg];
			}
			qsum = 0.;	//sum output flows and distribute convected signal
			for (i = 1; i <= nodout[inod]; i++) {
				iseg = nodseg[i][inod];
				qsum += qq[iseg];
			}
			for (i = 1; i <= nodout[inod]; i++) {
				iseg = nodseg[i][inod];
				conddown[iseg] = signal * qq[iseg] / qsum;
			}
		}
	}
	//calculate average value for segment, based on value at upstream node
	for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] >= 3 && segtyp[iseg] <= 5) {
		conddown[iseg] += metsig[iseg] * lseg[iseg] / 2.;
		signal = conddown[iseg] / (qq[iseg] + Qref1);		//convert to signal acting on vessel wall
		conddown[iseg] = signal / (signal + Sd50);			//saturate, max value 1
	}
}

void conductup1()
{
	extern int nnodfl, nseg;
	extern int* nodrank, * nodtyp, * nodout, * segtyp, ** nodseg;
	extern float* condup, * lseg, * conddown, * diam;
	extern float Su50, Sumax, L1, extlength1out;

	int inod, in, iseg, i;
	float signal, dsum, sigfac;

	//all signals initially evaluated at downstream nodes
	for (in = nnodfl; in >= 1; in--) {
		inod = nodrank[in];
		//boundary outflow nodes - assume a further segment with same conddown and length extlength1out
		if (abs(nodtyp[inod]) == 1 && nodout[inod] == 0) {
			iseg = nodseg[1][inod];
			sigfac = exp(-extlength1out / L1);
			condup[iseg] = conddown[iseg] * L1 * (1. - sigfac);
		}
		if (nodtyp[inod] >= 2) {
			//draining segments sum to give input of condup, weighted by diameters if diam
			signal = 0.;
			for (i = 1; i <= nodout[inod]; i++) {
				iseg = nodseg[i][inod];
				sigfac = exp(-lseg[iseg] / L1);
				signal += condup[iseg] * sigfac + conddown[iseg] * L1 * (1. - sigfac);
			}
			dsum = 0.;	//distribute conducted signal, with diameter weighting if diamweight = 1
			for (i = nodout[inod] + 1; i <= nodtyp[inod]; i++) {
				iseg = nodseg[i][inod];
				dsum += diam[iseg];	//diameter weighted
			}
			for (i = nodout[inod] + 1; i <= nodtyp[inod]; i++) {
				iseg = nodseg[i][inod];
				condup[iseg] = signal * diam[iseg] / dsum;	//diameter weighted
			}
		}
	}
	//calculate average value for segment, based on value at downstream node
	for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] >= 3 && segtyp[iseg] <= 5) {
		sigfac = exp(-lseg[iseg] / L1);
		signal = condup[iseg] * L1 / lseg[iseg] * (1. - sigfac)
			+ conddown[iseg] * L1 * (1. - L1 / lseg[iseg] * (1. - sigfac));
		condup[iseg] = Sumax * signal / (signal + Su50);		//saturate, max value Sumax
	}
	// For last equation, consider dCu/dx = Cd - Cu/L1, Cu(0) = Cu0 on [0,L]
	// Solution is Cu = L1*Cd + (Cu0 - L1 Cd)*exp(-x/L1)
	// Average of this over [0, L] is Cu0*L1/L*(1 - exp(-L/L1)) + Cd*L1*(1 - L1/L*(1 - exp(-L/L1)) 
}
