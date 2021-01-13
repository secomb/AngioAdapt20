/************************************************************************
flowtest_afterAngio
Use linear solver to identify flowing segments
Compare to flowtest_afterAdapt which has the ability to remove segments.
Use local pressure array to avoid disturbing actual values between iterations
TWS January 08, comment updated J. Alberding May 2009
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void solve(double *press);

void flowtest_afterAngio()
{
	extern int nseg,nnod,segfltot;
	extern int *segtyp,*ista,*iend;
	extern int **nodseg;
	extern float constvisc,facfp;
	extern float *lseg,*cond,*q;

	int iseg,inod;
	double *press;
	press = dvector(1,nnod);

	for(inod=1; inod<=nnod; inod++) press[inod] = 50.;
//Test connections only, assume same diameter for all segments
	for(iseg=1; iseg<=nseg; iseg++) if(segtyp[iseg] != 10) cond[iseg] = facfp*1.e4/lseg[iseg]/constvisc;
	solve(press);
	segfltot = 0;
	for(iseg=1; iseg<=nseg; iseg++) if(segtyp[iseg] != 10){
		q[iseg] = (press[ista[iseg]] - press[iend[iseg]])*cond[iseg];
		if(fabs(q[iseg]) > 0.0002){
			segfltot++;
			if(segtyp[iseg] != 3 && segtyp[iseg] != 4) segtyp[iseg] = 5;	//label as flowing, updated Feb. 2017
		}
		else if(segtyp[iseg] > 1)
			segtyp[iseg] = 1; //if flowing, label as non-flowing
	}
	free_dvector(press,1,nnod);
}
