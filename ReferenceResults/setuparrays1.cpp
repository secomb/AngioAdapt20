/************************************************************************
setuparrays1 - for AngioAdapt07.  TWS January 08
Set up arrays with dimensions nnod and nseg
Revised TWS2010
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void setuparrays1(int nseg, int nnod)
{
	extern int nsp, nodsegm, nseg_dim, nnod_dim;
	extern int* nspoint, * istart, * lowflow, * nk, * nodrank, * nodout, * nodtyp, * nodclasstyp, * segclasstyp;
	extern int* ista, * iend, * boundseg;
	extern int** nodnod, ** nodseg;
	extern float* segvar, * qq, * qqhd, * rseg, * ksseg, * cbar, * lseg, * ds, * nodvar, * qold, * hdold, * cond;
	extern float* segpress, * tau, * condup, * conddown, * metsig;
	extern float* stot, * stau, * spress, * uptrans, * downtrans, * segc;
	extern float** gamma1, ** qvseg, ** pvseg, ** pevseg, ** scos, ** signalvar;
	extern double* nodpress;

	FILE* ofp2;
	ofp2 = fopen("ArrayChanges.out", "a");

	nseg_dim = nseg;			//these are actual array dimensions.  TWS2010
	nnod_dim = nnod;

	ista = ivector(1, nseg); fprintf(ofp2, "ista,%i\n", nseg);
	iend = ivector(1, nseg); fprintf(ofp2, "iend,%i\n", nseg);
	lowflow = ivector(1, nseg); fprintf(ofp2, "lowflow,%i\n", nseg);
	nk = ivector(1, nnod); fprintf(ofp2, "nk,%i\n", nnod);//not nseg - error fixed 20 April 2010
	nodout = ivector(1, nnod); fprintf(ofp2, "nodout,%i\n", nnod);
	nodrank = ivector(1, nnod); fprintf(ofp2, "nodrank,%i\n", nnod);
	nodtyp = ivector(1, nnod); fprintf(ofp2, "nodtyp,%i\n", nnod);
	boundseg = ivector(1, nseg); fprintf(ofp2, "boundseg,%i\n", nseg);
	nspoint = ivector(1, nseg); fprintf(ofp2, "nspoint,%i\n", nseg);
	istart = ivector(1, nseg); fprintf(ofp2, "istart,%i\n", nseg);
	segclasstyp = ivector(1, nseg); fprintf(ofp2, "segclasstyp,%i\n", nseg);
	nodclasstyp = ivector(1, nnod); fprintf(ofp2, "nodclasstype,%i\n", nnod);

	nodnod = imatrix(1, nodsegm, 1, nnod); fprintf(ofp2, "nodnod,%i\n", nnod);
	nodseg = imatrix(1, nodsegm, 1, nnod); fprintf(ofp2, "nodseg,%i\n", nnod);

	cbar = vector(1, nseg); fprintf(ofp2, "cbar,%i\n", nseg);
	cond = vector(1, nseg); fprintf(ofp2, "cond,%i\n", nseg);
	conddown = vector(1, nseg); fprintf(ofp2, "conddown,%i\n", nseg);
	condup = vector(1, nseg); fprintf(ofp2, "condup,%i\n", nseg);
	downtrans = vector(1, nseg); fprintf(ofp2, "downtrans,%i\n", nseg);
	ds = vector(1, nseg); fprintf(ofp2, "ds,%i\n", nseg);
	hdold = vector(1, nseg); fprintf(ofp2, "hdold,%i\n", nseg);
	lseg = vector(1, nseg); fprintf(ofp2, "lseg,%i\n", nseg);
	metsig = vector(1, nseg); fprintf(ofp2, "metsig,%i\n", nseg);
	nodvar = vector(1, nnod); fprintf(ofp2, "nodvar,%i\n", nnod);
	qold = vector(1, nseg); fprintf(ofp2, "qold,%i\n", nseg);
	qq = vector(1, nseg); fprintf(ofp2, "qq,%i\n", nseg);
	qqhd = vector(1, nseg); fprintf(ofp2, "qqhd,%i\n", nseg);
	rseg = vector(1, nseg); fprintf(ofp2, "rseg,%i\n", nseg);
	ksseg = vector(1, nseg); fprintf(ofp2, "ksseg,%i\n", nseg);
	segpress = vector(1, nseg); fprintf(ofp2, "segpress,%i\n", nseg);
	segvar = vector(1, nseg); fprintf(ofp2, "segvar,%i\n", nseg);
	spress = vector(1, nseg); fprintf(ofp2, "spress,%i\n", nseg);
	stau = vector(1, nseg); fprintf(ofp2, "stau,%i\n", nseg);
	stot = vector(1, nseg); fprintf(ofp2, "stot,%i\n", nseg);
	tau = vector(1, nseg); fprintf(ofp2, "tau,%i\n", nseg);
	uptrans = vector(1, nseg); fprintf(ofp2, "uptrans,%i\n", nseg);
	segc = vector(1, nseg); fprintf(ofp2, "segc,%i\n", nseg);//TWS2010

	nodpress = dvector(1, nnod); fprintf(ofp2, "nodpress,%i\n", nnod);

	signalvar = matrix(1, nseg, 1, 6); fprintf(ofp2, "signalvar,%i\n", nseg);
	gamma1 = matrix(1, nseg, 1, nsp); fprintf(ofp2, "gamma1,%i\n", nseg);
	pevseg = matrix(1, nseg, 1, nsp); fprintf(ofp2, "pevseg,%i\n", nseg);
	pvseg = matrix(1, nseg, 1, nsp); fprintf(ofp2, "pvseg,%i\n", nseg);
	qvseg = matrix(1, nseg, 1, nsp); fprintf(ofp2, "qvseg,%i\n", nseg);
	scos = matrix(1, 3, 1, nseg); fprintf(ofp2, "scos,%i\n", nseg);

	fclose(ofp2);
}
