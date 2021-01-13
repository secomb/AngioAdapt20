/************************************************************************
flowtest - for AngioAdapt10
Use linear solver to identify segments without flow and set segtyp to 1
Set viscosities to constant values for this calculation
Use actual diameters. January 2011.
Use local pressure and flow arrays to avoid disturbing actual values between iterations
TWS January 08, comments updated J. Alberding May 2009
Revised TWS2011.  
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void solve(double *press);

void flowtest(int allsegs)
{
	extern int nseg,nnod,nnodbc,segfltot,flag;
	extern int *segtyp,*nodtyp,*nodout,*ista,*iend,*bcnod;
	extern int **nodseg;
	extern float constvisc,facfp;
	extern float *lseg,*cond,*diam,*hd;

	int iseg,inod,inodbc,*include;
	double *press,*q;	//make q and press local to avoid interference with actual flows
	press = dvector(1,nnod);
	q = dvector(1,nseg);
	include = ivector(1,nseg);

	for(iseg=1; iseg<=nseg; iseg++){
		include[iseg] = 0;
		if(allsegs == 1 && segtyp[iseg] != 10) include[iseg] = 1;
		if(allsegs == 0 && (segtyp[iseg] >= 3 && segtyp[iseg] <= 5)) include[iseg] = 1;
	}

	for(inod=1; inod<=nnod; inod++) press[inod] = 50.;
	for(iseg=1; iseg<=nseg; iseg++) if(include[iseg] == 1) cond[iseg] = facfp*pow(diam[iseg],4)/lseg[iseg]/constvisc;
	solve(press);
	segfltot = 0;
	for(iseg=1; iseg<=nseg; iseg++) if(include[iseg] == 1){
		q[iseg] = (press[ista[iseg]] - press[iend[iseg]])*cond[iseg];
		if(segtyp[iseg] == 3 || segtyp[iseg] == 4) segfltot++;
		else if(fabs(q[iseg]) > 0.0002){
			segtyp[iseg] = 5;	//label as type 5
			for(inodbc=1;inodbc<=nnodbc; inodbc++)	//label boundary segments as type 4
				if(ista[iseg] == bcnod[inodbc] || iend[iseg] == bcnod[inodbc]) segtyp[iseg] = 4;
			segfltot++;
		}
		else if(segtyp[iseg] > 1){
			flag = 1;
			if(segtyp[iseg] == 5) segfltot--;	//modified Feb. 2017
			segtyp[iseg] = 1;
		}
	}

	free_dvector(press,1,nnod);
	free_dvector(q,1,nseg);
	free_ivector(include,1,nseg);
}
