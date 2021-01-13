/************************************************************************
setup_hexperiodic - for AngioAdapt17GPU.  TWS 2018
Move startpoint of each segment into base domain
Make vector from startpoint to endpoint of each segment lie within base domain
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"


float length(float *x);

void remap_hexperiodic()
{
	extern int nnod,nseg,*ista,*iend,idomain,*boundseg;
	extern float *P1,*P2,*midpt,**cnode;

	int i,iseg,inod,inod1,inod2;

	//map nodes into base domain
	for(inod=1; inod<=nnod; inod++){
		for(i=1; i<=3; i++) P1[i] = cnode[i][inod] - midpt[i];
		length(P1);	//remaps P1 to base domain
		for(i=1; i<=3; i++) cnode[i][inod] = P1[i] + midpt[i];
	}		
	//test for domain crossings by segments
	for(iseg=1; iseg<=nseg; iseg++){
		inod1 = ista[iseg];
		inod2 = iend[iseg];
		for(i=1; i<=3; i++) P1[i] = cnode[i][inod2] - cnode[i][inod1];
		length(P1);	//remaps P1 to base domain
		for(i=1; i<=3; i++) P2[i] = cnode[i][inod1] + P1[i] - midpt[i];
		length(P2);	//value of idomain tests whether it is outside base domain
		boundseg[iseg] = idomain;
	}
}