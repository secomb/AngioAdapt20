/************************************************************************
setupgreens - for AngioAdapt07.  TWS 2010
subdivide segments into small elements, nnv = total number of elements
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void setupgreens()
{
	extern int nseg,nnv,*segtyp,*nspoint,*istart;
	extern float maxl,*lseg,*ds;
	int iseg,m;
	nnv = 0;
	for(iseg=1; iseg<=nseg; iseg++) if(segtyp[iseg] >= 3 && segtyp[iseg] <= 5){
		m = int(lseg[iseg]/maxl + 1);
		nspoint[iseg] = m;
		istart[iseg] = nnv + 1; 
		ds[iseg] = lseg[iseg]/m;
		nnv += m;
	}
	printf("total vessel points = %i\n", nnv);
}