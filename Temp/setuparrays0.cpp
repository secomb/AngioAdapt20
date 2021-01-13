/************************************************************************
setuparrays0 - for AngioAdapt07.  TWS January 08
Set up arrays with fixed dimensions
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void setuparrays0()
{
	extern int nsp,mxx,myy,mzz,nsprmax;
	extern int *errvesselcount,*errtissuecount,*imaxerrvessel,*imaxerrtissue;
	extern int ***nbou,***nbou_old;
	extern float *pmin,*pmax,*pmean,*mtiss,*mptiss,*g0old,*ptt,*ptpt,*qtsum,*qvsum,*errvessel,*qvfac;
	extern float *errtissue,*dqvsumdg0,*dqtsumdg0,*pinit,*p,*epsvessel,*epstissue,*eps,*g0facnew;
	extern float *x,*y,*P1,*P2,*P3,*midpt;
	extern float **gradp;
	extern float ***dtt;

	errvesselcount = ivector(1,nsp);
	errtissuecount = ivector(1,nsp);
	nbou = i3tensor(1,mxx,1,myy,1,mzz);
	nbou_old = i3tensor(1,mxx,1,myy,1,mzz);

	pmin = vector(1,nsp);
	pmax = vector(1,nsp);
	pmean = vector(1,nsp);
	mtiss = vector(1,(nsp+1));
	mptiss = vector(1,(nsp+1));
	g0old = vector(1,nsp);
	g0facnew = vector(1,nsp);

	ptt = vector(1,nsp);
	ptpt = vector(1,nsp);
	qtsum = vector(1,nsp);
	qvsum = vector(1,nsp);
	qvfac = vector(1,nsp);
	errvessel = vector(1,nsp);
	errtissue = vector(1,nsp);
	dqvsumdg0 = vector(1,nsp);
	dqtsumdg0 = vector(1,nsp);
	pinit = vector(1,nsp);
	p = vector(1,nsp);
	epsvessel = vector(1,nsp);
	epstissue = vector(1,nsp);
	eps = vector(1,nsp);
	imaxerrvessel = ivector(1,nsp);
	imaxerrtissue = ivector(1,nsp);

	x = vector(1,3);
	y = vector(1,3);
	P1 = vector(1,3);
	P2 = vector(1,3);
	P3 = vector(1,3);
	midpt = vector(1,3);

	gradp = matrix(1,3,1,nsp);

	dtt = f3tensor(1,mxx,1,myy,1,mzz);
}
