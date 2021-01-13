/************************************************************************
segreclass - for AngioAdapt10
reclassify type 0 and 1 segments after running flowtest
Type 0 and type 1 segments are removed if they are disconnected from other segments.
Note that this will not detect disconnected type 1 loops - needs more work!
TWS January 08, comments updated J. Alberding May 2009
Revised TWS2011
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void segreclass()
{
	extern int nseg,nsminus;
	extern int *segtyp,*ista,*iend,*bcnod,*nodtyp;
	extern int **nodseg;
	extern float *diam;

	int iseg,ii,nodtyp_ista,nodtyp_iend,flag,loop_control;

	flag = 1;
	loop_control = 0;
	while(flag == 1){
		flag = 0;
		for(iseg=1; iseg<=nseg; iseg++)	if(segtyp[iseg] <= 1){
			nodtyp_ista = nodtyp[ista[iseg]];
			nodtyp_iend = nodtyp[iend[iseg]];
			for(ii=1; ii<=nodtyp[ista[iseg]]; ii++) if(segtyp[nodseg[ii][ista[iseg]]] == 10) nodtyp_ista--;
			for(ii=1; ii<=nodtyp[iend[iseg]]; ii++) if(segtyp[nodseg[ii][iend[iseg]]] == 10) nodtyp_iend--;
			if(nodtyp_ista == 1 || (segtyp[iseg] == 1 && nodtyp_iend == 1)){	//type 0 segments have nodtyp_iend = 1
				printf(" %i", iseg);
				segtyp[iseg] = 10;
				diam[iseg] = 0.;
				flag = 1;
				nsminus++;
			}
		}
		loop_control++;
		if(loop_control > 1000){
			printf("*** Warning: segreclass loop exceeded loop control value (>1000)\n");
			flag = 0;
		}
	}
	printf("\n");
}
