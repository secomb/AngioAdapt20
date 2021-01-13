/**********************************************************************************
netstats.cpp
Program to create statistics on existing segments
Node classification (from putrank)
1: diverging
2: converging
0: other
Segment classification
0: O = other     (any other at start or end of the segment)
1: 2->1, 2->0, 0->1, 0->0 M = mesh (c,d)
2: 2->2 V = venule   (c=converging,c)
3: 1->2, 0->2, 1->0 C = capillary (d,c)
4: 1->1 A = arteriole (d=diverging node,d)
JP Alberding, 02-09. Comments updated J. Alberding May 2009
Edited by TWS, November 2010, Comment updated J. Alberding April 2011
Reworked TWS Feb. 2017
*********************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "nrutil.h"

float length(float *x);

void netstats(float po2cutoff)
{
	extern int nseg,nnod,nnt,nnv,imain,nsplus,nsminus,nsspr,nspruned,imainmax;
	extern int *ista,*iend,*segname,*segtyp,*nodtyp,**nodnod,*nodclasstyp,*segclasstyp,**nodseg;
	extern float pi1,qinsum;
	extern float *P1,*P2,*diam,*q, *sourcefac;
	extern float **cnode,**pt;
	extern float cs,flowfac,*qtsum,*bchd,**tissparam,vol;	//added Feb. 2017
	extern float totalexist_length, hypox_fract, time1, extraction;

	int totalexist_count,existflow_count,iseg,inod,inod1,inod2,inod_test,previnod_test,inod1_orig,i,max = 200,itp,nodtypflow,ncount;
	int *nodtypdist;
	float totalexist_volume,existflowing_length,totalexist_flow_volume,volume_increase;
	float Aflowing_length,Vflowing_length,Cflowing_length,Mflowing_length,otherflowing_length;
	float rad,mesh_length,seg_length,totalpts;
	float demand,consumption;
	FILE *ofp;

	totalexist_length = 0.;
	totalexist_count = 0;
	totalexist_volume = 0.;
	existflow_count = 0;
	existflowing_length = 0.;
	Aflowing_length = 0.;
	Vflowing_length = 0.;
	Cflowing_length = 0.;
	Mflowing_length = 0.;
	otherflowing_length = 0.;
	totalexist_flow_volume = 0.;
	nodtypdist = ivector(1,6);
	mesh_length = 0.;

	//Hypoxic fraction < po2cutoff 
	hypox_fract = 0.;
	totalpts = 0.;
	for (itp = 1; itp <= nnt; itp++){
		totalpts += sourcefac[itp];
		if (pt[itp][1] <= po2cutoff) hypox_fract += sourcefac[itp];
	}
	hypox_fract = hypox_fract*100./totalpts;

	//Classify segments
	for(iseg=1; iseg<=nseg; iseg++) if(segtyp[iseg] >= 3 && segtyp[iseg] <= 5){
		segclasstyp[iseg] = 0;
		if(q[iseg] > 0.){
			inod1 = ista[iseg];
			inod2 = iend[iseg];
		}
		else{
			inod1 = iend[iseg];
			inod2 = ista[iseg];
		}
		inod1_orig = inod1;	//flow is from this node
		if(nodtyp[inod1] == 2){	//search upstream for a node with 1, 3 or more segments
			ncount = 0;
			previnod_test = inod2;
			while(nodtyp[inod1] == 2 && ncount < nseg){
				for(i=1; i<=2; i++) if(nodnod[i][inod1] != previnod_test) inod_test = nodnod[i][inod1];
				previnod_test = inod1;
				inod1 = inod_test;
				ncount++;
			}
		}
		if(nodtyp[inod2] == 2){//search downstream for a node with 1, 3 or more segments
			ncount = 0;
			previnod_test = inod1_orig;
			while(nodtyp[inod2] == 2 && ncount < nseg){
				for(i=1; i<=2; i++) if(nodnod[i][inod2] != previnod_test) inod_test = nodnod[i][inod2];
				previnod_test = inod2;
				inod2 = inod_test;
				ncount++;
			}
		}
		segclasstyp[iseg] = 0;
		if(nodclasstyp[inod1] == 0){
			if(nodclasstyp[inod2] == 0) segclasstyp[iseg] = 1;
			else if(nodclasstyp[inod2] == 1) segclasstyp[iseg] = 1;
			else if(nodclasstyp[inod2] == 2) segclasstyp[iseg] = 3;
		}
		else if(nodclasstyp[inod1] == 1){
			if(nodclasstyp[inod2] == 0) segclasstyp[iseg] = 3;
			else if(nodclasstyp[inod2] == 1) segclasstyp[iseg] = 4;
			else if(nodclasstyp[inod2] == 2) segclasstyp[iseg] = 3;
		}
		else if(nodclasstyp[inod1] == 2){
			if(nodclasstyp[inod2] == 0) segclasstyp[iseg] = 1;
			else if(nodclasstyp[inod2] == 1) segclasstyp[iseg] = 1;
			else if(nodclasstyp[inod2] == 2) segclasstyp[iseg] = 2;
		}
	}
	//Total vessel length and volume
	for(iseg=1; iseg<=nseg; iseg++) if(segtyp[iseg] != 10){
		for(i=1; i<=3; i++)	P1[i] = cnode[i][ista[iseg]] - cnode[i][iend[iseg]];
		seg_length = length(P1);
		rad = diam[iseg]/2.;
		volume_increase = pi1*rad*rad*seg_length;
		totalexist_count++;
		totalexist_length += seg_length;
		totalexist_volume += volume_increase;
		if(segtyp[iseg] >= 3 && segtyp[iseg] <= 5){
			existflow_count++;
			existflowing_length += seg_length;
			totalexist_flow_volume += volume_increase;
			if(segclasstyp[iseg] == 4) Aflowing_length += seg_length;
			if(segclasstyp[iseg] == 2) Vflowing_length += seg_length;
			if(segclasstyp[iseg] == 3) Cflowing_length += seg_length;
			if(segclasstyp[iseg] == 1) Mflowing_length += seg_length;
			if(segclasstyp[iseg] == 0) otherflowing_length += seg_length;
		}
	}

	ofp = fopen("LengthVolumeStats1.out","a");
	fprintf(ofp,"%5i %7.2f %8.1f %8.1f  %10i  %10i  %10.1f %10.1f  %10i\n",
		imain, time1, totalexist_length, existflowing_length, totalexist_count, existflow_count, totalexist_volume,
		totalexist_flow_volume, nnv);
	fclose(ofp);

	//Vessel types - length of each
	ofp = fopen("LengthVolumeStats2.out","a");
	fprintf(ofp,"%5i %7.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %7i %7i %7i %7i %7i\n",
		imain, time1, totalexist_length, Aflowing_length, Vflowing_length, Cflowing_length,
		Mflowing_length, otherflowing_length, nnv,nsplus,nsminus,nsspr,nspruned);
	fclose(ofp);

	ofp = fopen("VesselTypes.ps", "a");
	fprintf(ofp, "0 0 0 setrgbcolor\n");
	fprintf(ofp, "%i mx %g my cf\n",imain,totalexist_length/1000.);//total, black
	fprintf(ofp, "0 1 1 setrgbcolor\n");
	fprintf(ofp, "%i mx %g my cf\n",imain,existflowing_length/1000.);//total flowing, light blue
	fprintf(ofp, "0 1 0 setrgbcolor\n");
	fprintf(ofp, "%i mx %g my cf\n",imain,Mflowing_length/1000.);//mesh, green
	fprintf(ofp, "1 0 0 setrgbcolor\n");
	fprintf(ofp, "%i mx %g my cf\n",imain,Aflowing_length/1000.);//A, red
	fprintf(ofp, "0 0 1 setrgbcolor\n");
	fprintf(ofp, "%i mx %g my cf\n",imain,Vflowing_length/1000.);//V, blue
	fprintf(ofp, "1 0 1 setrgbcolor\n");
	fprintf(ofp, "%i mx %g my cf\n",imain,Cflowing_length/1000.);//C, light purple
	fclose(ofp);

	//Node stats start here - number of flowing segments attached to each node
	for(i=1; i<=6; i++) nodtypdist[i] = 0;
	for(inod=1; inod<=nnod; inod++){
		nodtypflow = 0;
		for(i=1; i<=nodtyp[inod]; i++) if(segtyp[nodseg[i][inod]] >= 3 && segtyp[nodseg[i][inod]] <= 5) nodtypflow++;
		nodtypflow = IMIN(nodtypflow,6);
		if(nodtypflow != 0) nodtypdist[nodtypflow]++;
	}
	ofp = fopen("NodtypStats.out","a");
	fprintf(ofp,"%5i %7.2f      ",imain, time1);
	for(i=1; i<=6; i++) fprintf(ofp,"%4i ",nodtypdist[i]);
	fprintf(ofp,"\n");
	fclose(ofp);

	//oxygen transport statistics
	demand = tissparam[1][1]*vol*nnt/flowfac;		//nl/min
	consumption = -qtsum[1]/flowfac;	//nl/min
	extraction = 100.*consumption/qinsum/cs/bchd[1]; //in percent
	ofp = fopen("OxygenStats.out","a");
	fprintf(ofp,"%5i %7.2f %f %f %f %f %f\n",imain, time1, demand,consumption,qinsum,extraction,hypox_fract);
	fclose(ofp);

	ofp = fopen("OxygenStats.ps","a");
	fprintf(ofp, "0 0 0 setrgbcolor\n");
	fprintf(ofp, "%i mx %g my cf\n",imain,qinsum);	//black - flow in nl/min
	fprintf(ofp, "1 0 0 setrgbcolor\n");
	fprintf(ofp, "%i mx %g my cf\n",imain,consumption);	//red - total oxygen consumption in nl/min
	fprintf(ofp, "0 0 1 setrgbcolor\n");
	fprintf(ofp, "%i mx %g my cf\n",imain,extraction);	//blue - oxygen extraction %
	fprintf(ofp, "0 1 0 setrgbcolor\n");
	fprintf(ofp, "%i mx %g my cf\n",imain,hypox_fract);	//green - hypoxic fraction %
	fprintf(ofp, "0 0 0 setrgbcolor\n");
	fclose(ofp);

	free_ivector(nodtypdist,1,6);
}