/*****************************************************
Evaluate histograms of solute levels, distance to nearest vessels and angles.
TWS December 07. Comments updated J. Alberding May 09. Revised TWS January 2011
******************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <math.h>
#include "nrutil.h"
#include <stdio.h>

void histogram(void)
{
	extern int nnod,nnt,nsp,*nodtyp,**nodnod;
	extern float pi1,*pmin,*pmax,*pmean,*dtmin;
	extern float **pt,**cnode;

	int i,itp,isp,nstat,inod,j,m,testgt180,trinode_count;
	float step,dev,dtmin_max;
	float *samp,*stat,*cumu,*theta,*alpha;
	FILE *ofp;

	ofp = fopen("Histograms.out", "w");
	for(isp=1; isp<=nsp; isp++){
		if(isp == 1) step = 2.;
		if(isp == 2) step = 0.05;
		nstat = IMIN(IMAX(pmax[isp]/step + 2,41),101);// 41 <= nstat <= 101
		samp = vector(1,nstat);
		stat = vector(1,nstat);
		cumu = vector(1,nstat);
		for(i=1; i<=nstat; i++){
			samp[i] = step*(i - 1);
			stat[i] = 0;
		}
		dev = 0.;
		for(itp=1; itp<=nnt; itp++){
			dev += SQR(pmean[isp] - pt[itp][isp]);
			for(i=1; i<=nstat; i++) if(pt[itp][isp] <= samp[i]){
				stat[i]++;
				goto binned;
			}
			binned:;
		}
		dev = sqrt(dev/nnt);
		for(i=1; i<=nstat; i++) stat[i] = stat[i]*100./nnt;
		cumu[1] = stat[1];
		for(i=2; i<=nstat; i++) cumu[i] = cumu[i-1] + stat[i];
		fprintf(ofp,"Histogram data for solute %i\n", isp);
		fprintf(ofp,"value  percent cumul. percent\n");
		for(i=1; i<=nstat; i++) fprintf(ofp,"%g %7.2f %7.2f\n", samp[i],stat[i],cumu[i]);
		fprintf(ofp,"Mean solute = %f deviation = %f min = %f max  = %f\n", pmean[isp],dev,pmin[isp],pmax[isp]);
		fprintf(ofp,"-----------------------------------------------------------\n");
		free_vector(samp,1,nstat);
		free_vector(stat,1,nstat);
		free_vector(cumu,1,nstat);
	}

	step = 2.;
	dtmin_max = 0.;
	for(itp=1; itp<=nnt; itp++)	if(dtmin[itp] > dtmin_max) dtmin_max = dtmin[itp];
	nstat = IMIN(IMAX(dtmin_max/step + 2,41),101);// 41 <= nstat <= 101
	samp = vector(1,nstat);
	stat = vector(1,nstat); 
	cumu = vector(1,nstat);
	for(i=1; i<=nstat; i++){
		samp[i] = step*(i - 1);
		stat[i] = 0;
	}
	for(itp=1; itp<=nnt; itp++){
		for(i=1; i<=nstat; i++)	if(dtmin[itp] < samp[i]){
			stat[i]++;
			goto binned1;
		}
		binned1:;
	}
	for(i=1; i<=nstat; i++) stat[i] = stat[i]*100./nnt;
	cumu[1] = 0.;
	for(i=2; i<=nstat; i++) cumu[i] = cumu[i-1] + stat[i];
	fprintf(ofp,"Histogram data for distance to nearest vessel\n");
	fprintf(ofp,"value  percent cumul. percent\n");
	for(i=1; i<=nstat; i++) fprintf(ofp,"%g %7.2f %7.2f\n", samp[i],stat[i],cumu[i]);
	fprintf(ofp,"distmax  = %f\n",dtmin_max);
	fprintf(ofp,"-----------------------------------------------------------\n");
	free_vector(samp,1,nstat);
	free_vector(stat,1,nstat); 
	free_vector(cumu,1,nstat);

//histogram of angles
	theta = vector(1,3);
	alpha = vector(1,3);
	step = 10.;	//every 10 degrees
	nstat = 360./step + 1.001;
	stat = vector(1,nstat);
	samp = vector(1,nstat);
	cumu = vector(1,nstat);
	for(i=1; i<=nstat; i++){
		samp[i] = step*(i - 1);
		stat[i] = 0;
	}
	trinode_count = 0;
	//This only works for 2D networks, does not work for hexperiodic *************************************
	for(inod=1; inod<=nnod; inod++) if(nodtyp[inod] == 3){	//angles about nodtyp == 3 only
		trinode_count++;
		for(j=1; j<=3; j++)	theta[j] = 180./pi1*
			atan2(cnode[2][nodnod[j][inod]] - cnode[2][inod],cnode[1][nodnod[j][inod]] - cnode[1][inod]);
//sort segment angles in ascending order
		if(theta[2] > theta[3]){
			alpha[1] = theta[2];
			theta[2] = theta[3];
			theta[3] = alpha[1];
		}
		if(theta[1] > theta[2]){
			alpha[1] = theta[1];
			theta[1] = theta[2];
			theta[2] = alpha[1];
		}
		if(theta[2] > theta[3]){
			alpha[1] = theta[2];
			theta[2] = theta[3];
			theta[3] = alpha[1];
		}
//angles between segments, all positive
		alpha[1] = theta[2] - theta[1];
		alpha[2] = theta[3] - theta[2];
		alpha[3] = 360. - alpha[1] - alpha[2];
		testgt180 = 0;
		for(m=1; m<=3; m++){
			if(alpha[m] > 180.) testgt180++;
			for(j=1; j<=nstat; j++) if(alpha[m] <= samp[j]){
				stat[j]++;
				goto binned2;
			}
binned2:;
		}
	}
	for(i=1; i<=nstat; i++) stat[i] = stat[i]*100./(3.*trinode_count);
	cumu[1] = 0.;
	for(i=2; i<=nstat; i++) cumu[i] = cumu[i-1] + stat[i];
	fprintf(ofp,"Histogram data for angles at nodtyp = 3\n");
	fprintf(ofp,"value  percent cumul. percent\n");
	for(i=1; i<=nstat; i++) fprintf(ofp,"%g %7.2f %7.2f\n", samp[i],stat[i],cumu[i]);
	fprintf(ofp,"-----------------------------------------------------------\n");
	fclose(ofp);
	free_vector(theta,1,3);
	free_vector(alpha,1,3);
	free_vector(samp,1,nstat);
	free_vector(stat,1,nstat); 
	free_vector(cumu,1,nstat);
}
