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
float length(float *x);

void conductconvect()
{
	extern int nnodfl,nseg,num_runs_actual,slsegdiv,hexperiodic;
	extern int *nodtyp,*nodtyp,*nodfltyp,*segtyp,*lowflow,*ista,*iend;
	extern float kp1,km1,kc1,ks1,L1,Qref1,pO2ref1,nmeta1,thresh,threshn,thresh50,req;//TWS2010,JPA2013
	extern float lambdaGF2, cfacGF2, ssGF2, bbGF2, ccGF2;
	extern float *lseg,*metsig,*p,*x,*midpt;
	extern float **pvseg,**pevseg,**cnode;

	int iseg,j;
	float pv;
	//float k1tissue = 0.5, thresh1 = 0.25;
	float k1tissue = 1., thresh1 = 0.25;
	float GF,GFx,GFy,GFz;

	putrank();

	//local signal.  For low-flow segments, use extravascular PO2
	for(iseg=1; iseg<=nseg; iseg++)	if(segtyp[iseg] >= 3 && segtyp[iseg] <= 5){
		if(lowflow[iseg]) pv = pevseg[iseg][1];
		else pv = pvseg[iseg][1];
		//Include growth factor effect in metsig.  TWS 2010
		for(j=1; j<=3; j++) x[j] = (cnode[j][ista[iseg]] + cnode[j][iend[iseg]])*0.5 - midpt[j];
		if(hexperiodic) length(x);	//map into base domain
		for(j=1; j<=3; j++) x[j] += midpt[j];
		evalGF(2, req, x, lambdaGF2, cfacGF2, ssGF2, bbGF2, ccGF2, &GF, &GFx, &GFy, &GFz);
		if(pv < 0.) metsig[iseg] = 1.;
		else metsig[iseg] =1. - 1./(1. + (1. - k1tissue)*powf(pO2ref1/pv,nmeta1) + k1tissue*powf(GF/thresh1,nmeta1));
		//Hill-type model, March 2019. Combine GF effect with O2 effect: keeps metsig in range 0 to 1 
	}
	convectdown1();
	conductup1();
}

void convectdown1()
{
	extern int nnodfl,nseg;
	extern int *nodrank,*nodtyp,*nodout,*segtyp,**nodseg;
	extern float Qref1,extlength1in;
	extern float *conddown,*lseg,*qq,*metsig,**pevseg;

	int inod,in,iseg,i;
	float signal,qsum;

	//all signals initially evaluated at upstream nodes
	for(in=1; in<=nnodfl; in++){
		inod = nodrank[in];
		//boundary inflow node: assume a further segment with same metsig and length extlength1in
		if(abs(nodtyp[inod]) == 1 && nodout[inod] == 1){ 
			iseg = nodseg[1][inod];
			conddown[iseg] = metsig[iseg]*extlength1in;
		}
		if(nodtyp[inod] >= 2){ 
			//input segments sum to give input of conddown
			signal = 0.;
			for(i=nodout[inod]+1; i<=nodtyp[inod]; i++){
				iseg = nodseg[i][inod];
				signal += conddown[iseg] + metsig[iseg]*lseg[iseg];
			}
			//sum output flows and distribute convected signal
			qsum = 0.;
			for(i=1; i<=nodout[inod]; i++){
				iseg = nodseg[i][inod];
				qsum += qq[iseg];
			}
			for(i=1; i<=nodout[inod]; i++){
				iseg = nodseg[i][inod];
				conddown[iseg] = signal*qq[iseg]/qsum;
			}
		}
	}
	//calculate average value for segment, based on value at upstream node.  Apply log function
	for(iseg=1; iseg<=nseg; iseg++) if(segtyp[iseg] >= 3 && segtyp[iseg] <= 5){
		conddown[iseg] += metsig[iseg]*lseg[iseg]/2.;
		conddown[iseg] = log10(1. + conddown[iseg]/(qq[iseg] + Qref1));
	}
}

void conductup1()
{
	extern int nnodfl,nseg;
	extern int *nodrank,*nodtyp,*nodout,*segtyp,**nodseg;
	extern float *condup,*lseg,*conddown,*diam;
	extern float L1,extlength1out;//TWS2010

	int diamweight = 1;	//set diamweight = 1 for diameter weighted upstream conduction - TWS2014
	int inod,in,iseg,i;
	float signal,dsum,sigfac;

	//all signals initially evaluated at downstream nodes
	for(in=nnodfl; in>=1; in--){
		inod = nodrank[in];
		//boundary outflow nodes - assume a further segment with same conddown and length extlength1out
		if(abs(nodtyp[inod]) == 1 && nodout[inod] == 0){ 
			iseg = nodseg[1][inod];
			sigfac = exp(-extlength1out/L1);
			condup[iseg] = conddown[iseg]*L1*(1.-sigfac);
		}
		if(nodtyp[inod] >= 2){ 
			//draining segments sum to give input of condup, weighted by diameters if diam
			signal = 0.;
			for(i=1; i<=nodout[inod]; i++){
				iseg = nodseg[i][inod];
				sigfac = exp(-lseg[iseg]/L1);
				signal += condup[iseg]*sigfac + conddown[iseg]*L1*(1.-sigfac);
			}
			//distribute conducted signal, with diameter weighting if diamweight = 1 - TWS 2014
			dsum = 0.;
			for(i=nodout[inod]+1; i<=nodtyp[inod]; i++){
				iseg = nodseg[i][inod];
				if(diamweight == 1) dsum += diam[iseg];	//diameter weighted
				else dsum += 1.;	//diameter independent
			}
			for(i=nodout[inod]+1; i<=nodtyp[inod]; i++){
				iseg = nodseg[i][inod];
				if(diamweight == 1) condup[iseg] = signal*diam[iseg]/dsum;	//diameter weighted
				else condup[iseg] = signal/dsum; //diameter independent
			}
		}
	}
	//calculate average value for segment, based on value at downstream node
	for(iseg=1; iseg<=nseg; iseg++) if(segtyp[iseg] >= 3 && segtyp[iseg] <= 5){
		sigfac = exp(-lseg[iseg]/L1);
		condup[iseg] = condup[iseg]*L1/lseg[iseg]*(1.-sigfac)
			+ conddown[iseg]*L1*(1.-L1/lseg[iseg]*(1.-sigfac));
	}
}
