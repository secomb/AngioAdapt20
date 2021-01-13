/************************************************************************
analyzenet - for AngioAdapt07.  TWS October 07
Set up nodtyp, nodseg, nodnod arrays based on flowing segments if allsegs = 0,
all segments if allsegs = 1
These values are subsequently overridden by putrank
Revised TWS2010.
Modified for hexagonal periodic domain, February 2018
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

float length(float *x);

void analyzenet(int allsegs)
{
	extern int nseg,nnod,nnodbc,nodsegm,actnum,nodnamemax,segnamemax,hexperiodic;
	extern int *bcnod,*bcnodname,*nodtyp,*segtyp,*nodname,*segname,*ista,*iend,*boundseg;
	extern int **segnodname,**nodnod,**nodseg;
	extern float aphex,aphexp,**hex_norm,**hex_tang;
	extern float tlength;
	extern float *P1,*diam,*lseg,*rseg,*q,*qq,*axt,*ayt,*azt;
	extern float **cnode,**scos;

	int k,iseg,inod,inod1,inod2,inodbc;
	int *include;
	include = ivector(1,nseg);

	for(inod=1; inod<=nnod; inod++) nodtyp[inod] = 0;
	for(iseg=1; iseg<=nseg; iseg++){
		include[iseg] = 0;
		if(allsegs == 1 && segtyp[iseg] != 10) include[iseg] = 1;
		if(allsegs == 0 && (segtyp[iseg] >= 3 && segtyp[iseg] <= 5)) include[iseg] = 1;
	}
	for(iseg=1; iseg<=nseg; iseg++)	if(include[iseg] == 1){	//Search for nodes corresponding to this segment
		for(inod=1; inod<=nnod; inod++) if(nodname[inod] == segnodname[1][iseg]){
			ista[iseg] = inod;
			goto foundit1;
		}
		printf("*** Error: No matching node found for nodname %i\n", segnodname[1][iseg]);
		foundit1:;
		for(inod=1; inod<=nnod; inod++) if(nodname[inod] == segnodname[2][iseg]){
			iend[iseg] = inod;
			goto foundit2;
		}
		printf("*** Error: No matching node found for nodname %i\n", segnodname[2][iseg]);
		foundit2:;
//Setup nodtyp, nodseg and nodnod
		inod1 = ista[iseg];
		inod2 = iend[iseg];
		nodtyp[inod1]++;
		nodtyp[inod2]++;
		if(nodtyp[inod1] > nodsegm) printf("*** Error: Too many segments connected to node %i\n", inod1);
		if(nodtyp[inod2] > nodsegm) printf("*** Error: Too many segments connected to node %i\n", inod2);
		nodseg[nodtyp[inod1]][inod1] = iseg;
		nodseg[nodtyp[inod2]][inod2] = iseg;
		nodnod[nodtyp[inod1]][inod1] = inod2;
		nodnod[nodtyp[inod2]][inod2] = inod1;
	}
	for (inodbc=1; inodbc<=nnodbc; inodbc++){
		bcnod[inodbc] = 0;
//Search for node corresponding to this node name
		for(inod=1; inod<=nnod; inod++) if(nodname[inod] == bcnodname[inodbc]){
			bcnod[inodbc] = inod;
			if(nodtyp[inod] != 1) printf("*** Error: Boundary node %i is not a 1-segment node\n", inod);
			goto foundit;
		}
		printf("*** Error: No matching node found for nodname %i, %i\n", nodname[inod],bcnodname[inodbc]);
		foundit:;
	}
	tlength = 0.;
	for(iseg=1; iseg<=nseg; iseg++) if(include[iseg] == 1){
		rseg[iseg] = diam[iseg]/2.0;
		qq[iseg] = fabs(q[iseg]);
		for(k=1; k<=3; k++)	P1[k] = cnode[k][iend[iseg]] - cnode[k][ista[iseg]];
		lseg[iseg] = length(P1);	//P1 may be modified if hexperiodic is active
		for(k=1; k<=3; k++)	scos[k][iseg] = P1[k]/lseg[iseg];
		tlength += lseg[iseg];
	}
	free_ivector(include,1,nseg);

//actnum = total number of active segments
	actnum = 0;
	for(iseg=1; iseg<=nseg; iseg++)	if(segtyp[iseg] == 0) actnum++;
//needed for naming new nodes and segments
	nodnamemax = 0;
	segnamemax = 0;
	for(iseg=1; iseg<=nseg; iseg++) if(segname[iseg] > segnamemax) segnamemax = segname[iseg];
	for(inod=1; inod<=nnod; inod++) if(nodname[inod] > nodnamemax) nodnamemax = nodname[inod];
}