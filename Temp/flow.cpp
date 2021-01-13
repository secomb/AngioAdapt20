/************************************************************************
flow - for AngioAdapt07
TWS October 07
Flows in nl/min, pressures in mmHg, diameters in um, viscosity in cP
For nonlinear iterations, if convergence is slow, hematocrit is increasingly underrelaxed.
This eventually forces convergence even when system is unstable due to rheological effects
(see work of R.L. Carr et al.).
Revised TWS2010
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

float viscor(float d, float h);
void analyzenet(int allsegs);
void solve(double *nodpress);
void picturenetwork(float *nodvar, float *segvar, const char fname[]);
void dishem_generalized(void);
void putrank(void);

void flow()
{
	extern int nseg,nnod,nnodbc,varyviscosity,phaseseparation,nnodfl,nsminus;
	extern int segfltot,nitmax1;
	extern int *segtyp,*segname,*nodtyp,*nodrank,*nodout,*nodname,*ista,*iend,*bcnod;
	extern int **nodseg;
	extern float constvisc,consthd,facfp,qtol,hdtol;
	extern float *diam,*q,*qq,*hd,*cond,*lseg,*diam,*segpress,*tau,*qold,*hdold,*nodvar,*segvar;
	extern float **cnode;
	extern double *nodpress;

	int iseg,inod,niter,errsegq,errseghd;
	int ii,in,nodt,nout;
	float visc,maxqerr,maxhderr,qchange,hdchange,relax,sumin,sumout;

	for(iseg=1; iseg<=nseg; iseg++) q[iseg] = 0.;
	for(inod=1; inod<=nnod; inod++) nodpress[inod] = 50.;
	relax = 1.;
	for(niter=1; niter<=nitmax1; niter++){
		if(niter%5 == 0) relax = 0.8*relax;
		for(iseg=1; iseg<=nseg; iseg++) if(segtyp[iseg] >= 3 && segtyp[iseg] <= 5){
			qold[iseg] = q[iseg];
			hdold[iseg] = hd[iseg];
			if(varyviscosity == 1){
				hd[iseg] = FMIN(hd[iseg],0.9);
				visc = viscor(diam[iseg],hd[iseg]);
			}
			else visc = constvisc;
			cond[iseg] = facfp*pow(diam[iseg],4)/lseg[iseg]/visc;
		}
		solve(nodpress);
		for(iseg=1; iseg<=nseg; iseg++) if(segtyp[iseg] >= 3 && segtyp[iseg] <= 5)
			q[iseg] = (nodpress[ista[iseg]] - nodpress[iend[iseg]])*cond[iseg];
		//test for flow conservation
		if(phaseseparation == 1) putrank();
		for(in=1; in<=nnodfl; in++){	//scan nodes in downstream order
			inod = nodrank[in];
			nodt = nodtyp[inod];
			nout = nodout[inod];
			if(nodt > 1){	//don't do for network boundary nodes
				sumin = 0.;
				for(ii=nout+1; ii<=nodt; ii++){ //inflows
					iseg = nodseg[ii][inod];
					sumin += fabs(q[iseg]);
				}
				sumout = 0.;	
				for(ii=1; ii<=nout; ii++){		//outflows	
					iseg = nodseg[ii][inod];
					sumout += fabs(q[iseg]);	//check conservation of flow and hematocrit
				}
				if(sumin + sumout != 0.){
					if(fabs(sumin-sumout) > 0.01){
						printf("\n");
						printf("*** Error: Flow conservation violation at node %i\n", nodname[inod]);
					}
				}
			}
		}
		//calculate segment hematocrits
		for(iseg=1; iseg<=nseg; iseg++) segvar[iseg] = q[iseg];
		for(inod=1; inod<=nnod; inod++) nodvar[inod] = inod;
		picturenetwork(nodvar,segvar,"Network_bfofflows.ps");

		if(phaseseparation == 1){
			putrank();
			dishem_generalized();
		}
		else for(iseg=1; iseg<=nseg; iseg++) if(segtyp[iseg] >= 3 && segtyp[iseg] <= 5) hd[iseg] = consthd;
		//compare hd and q with previous values
		maxqerr = 0.;
		maxhderr = 0.;
		errsegq = 0;
		errseghd = 0;
		for(iseg=1; iseg<=nseg; iseg++) if(segtyp[iseg] >= 3 && segtyp[iseg] <= 5){
			qchange = q[iseg] - qold[iseg];
			hdchange = hd[iseg] - hdold[iseg];
			hd[iseg] = hdold[iseg] + relax*hdchange;
			if(hd[iseg] > 1.F){
				printf("*** Warning: hd[iseg] > 1\n");
			}
			if(fabs(qchange) >= maxqerr){
				maxqerr = fabs(qchange);
				errsegq = iseg;
			}
			if(fabs(hdchange) >= maxhderr){
				maxhderr = fabs(hdchange);
				errseghd = iseg;
			}
		}
		if(maxqerr < qtol && maxhderr < hdtol) goto converged;
	}
	printf("*** Warning: Nonlinear iteration not converged\n");
	printf("*** Flow error = %f at segment %i, h'crit error = %f at segment %i\n",maxqerr,errsegq,maxhderr,errseghd);
converged:;

	//recalculate segment hematocrits - otherwise errors can occur if relax < 1.  TWS  April 2011
	if(phaseseparation == 1  && relax < 0.9999){
		putrank();
		dishem_generalized();
	}

	for (iseg=1; iseg<=nseg; iseg++) if(segtyp[iseg] >= 3 && segtyp[iseg] <= 5){
		segpress[iseg] = (nodpress[ista[iseg]] + nodpress[iend[iseg]])/2.;
		tau[iseg] = fabs(nodpress[ista[iseg]] - nodpress[iend[iseg]])*1333.*diam[iseg]/lseg[iseg]/4.;
		qq[iseg] = fabs(q[iseg]);
	}
}
