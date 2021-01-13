/************************************************************************
adapt.cpp
Diameter adaptation, drop out of segments.  TWS October 07
Modified to use 2001 adaptation model, December 2010
Diameter adaptation for type 5 but not for type 4
Diammax added May 2012--diameters "capped" as in parallel
vessel network networks growing from venous to arteriolar
side (higher threshold in angioparams) produced vessels too
large for Green's function to handle. JPA May 2012
Edited TWS Feb. 2017
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void adapt()
{
	extern int nseg,nnod,eliminate,inftrans,nsminus,nspruned,*segtyp,*ista,*iend,adaptd;
	extern float diamthresh,eqtime,kp1,km1,kc1,ks1,ranks1,J01, tauref1, tauref2, timestep1,thresh,threshn,thresh50,req;
	extern float kmetdiam;
	extern float *segpress,*tau,*adaptpar,*condup,*conddown,*diam,*stot,*ksseg;
	extern float *stau,*spress,*uptrans,*downtrans,*x,*nodvar,**signalvar;
	extern float **cnode;
	extern float *metsig;	//*****************not currently used

	int iseg,ii,diammaxcount;
	float tauexp,pressl,diammax;
	float klocal = 0.;		//*********************
	float metdiamfac; //diamref = 6., metdiamfac;
//	float nmetdiam = 0.40;
//	float kmetdiam = 12.; // 0.5; // 0.3;	// 0.7375; //Created by TS 10/2019, replaces diamref, nmetdiam in new metdiamfac. JA

	for(iseg=1; iseg<=nseg; iseg++) for(ii=1; ii<=6; ii++) signalvar[iseg][ii] = 0.;

	diammax = 44.0;//50.0 too large for Green's function
	diammaxcount = 0;

	for(iseg=1; iseg<=nseg; iseg++){
		if(segtyp[iseg] == 3 || segtyp[iseg] == 5){	//don't adapt type 4 segments
			stau[iseg] = log10(tau[iseg]/(1. + tau[iseg]/tauref2) + tauref1);	//TWS 2020
			pressl = segpress[iseg];
			if(pressl > 100.) pressl = 100.;
			if(pressl < 10.) pressl = 10.;
			tauexp = adaptpar[1] + adaptpar[2] * exp(adaptpar[3] * pow(log10(log10(pressl)), adaptpar[4]));
			spress[iseg] = -kp1*log10(tauexp);
			downtrans[iseg] = 0.;
			uptrans[iseg] = 0.;
			metdiamfac = 1. / (1. + kmetdiam * diam[iseg]);		//TWS 2020
			if(inftrans == 1){	//downstream convected and upstream conducted responses
				downtrans[iseg] = metdiamfac * km1 * conddown[iseg];
				uptrans[iseg] = metdiamfac * km1 * kc1 * condup[iseg] / (condup[iseg] + J01);
			}
			stot[iseg] = stau[iseg] + spress[iseg] + uptrans[iseg] + downtrans[iseg] + klocal * metsig[iseg] - ksseg[iseg];
			if (adaptd == 1 || stot[iseg] > 0.) diam[iseg] += diam[iseg] * stot[iseg] * timestep1 / eqtime;	//Jan. 2020

			//test if diameters too large, too small
			if(diam[iseg] > diammax){
				diam[iseg] = diammax;
				diammaxcount++;
				printf(" diam%i",iseg);
			}
			if(diam[iseg] < diamthresh){
				if(eliminate == 1){
					printf(" %i", iseg);
					diam[iseg] = 0.;
					segtyp[iseg] = 10;
					nspruned++;
					nsminus++;
				}
				else diam[iseg] = diamthresh;	//use this option to inhibit pruning
			}
		}
		else{	//nonflowing and type 4 segments
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
