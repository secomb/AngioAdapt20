/************************************************************************
setuparrays2 - for AngioAdapt07.  TWS January 08
Set up arrays with dimensions nnv and nnt
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void setuparrays2(int nnv, int nnt, int actnum)
{
	extern int nsp, nnod, nseg, nnv_dim, nnt_dim, nnt_old_dim, actnum_dim, nsprmax;
	extern int* indx, * mainseg, * activenode, * new_segnum, * new_segtyp, * new_active_segs, ** tisspoints;
	extern float* rhs, * qvtemp, * dtmin, * sourcefac, * sourcefacGF;
	extern float** qv, ** pv, ** pev, ** pvt, ** pvprev, ** qvprev, ** cv, ** dcdp, ** cv0, ** conv0, ** gvv;
	extern float** qt, ** qt_old, ** pt, ** pt_old, ** ptprev, ** ptv, ** qcoeff1, ** ax, ** al, ** activegr;
	extern double** rhsg, ** mat;
	extern double* rhsl, * matx;

	FILE* ofp2;
	ofp2 = fopen("ArrayChanges.out", "a");

	nnv_dim = nnv;	//these are the actual dimensions of the arrays
	nnt_dim = nnt;
	nnt_old_dim = nnt;
	actnum_dim = actnum;
	int actnum_temp = IMAX(actnum, 1);

	indx = ivector(1, nnv); fprintf(ofp2, "indx,%i\n", nnv);
	mainseg = ivector(1, nnv); fprintf(ofp2, "mainseg,%i\n", nnv);
	tisspoints = imatrix(1, 3, 1, nnt); fprintf(ofp2, "tisspoints,%i\n", nnt);
	activenode = ivector(1, actnum_temp); fprintf(ofp2, "activenode,%i\n", actnum_temp);
	dtmin = vector(1, nnt); fprintf(ofp2, "dtmin,%i\n", nnt);
	sourcefac = vector(1, nnt); fprintf(ofp2, "sourcefac,%i\n", nnt);
	sourcefacGF = vector(1, nnt); fprintf(ofp2, "sourcefacGF,%i\n", nnt);
	rhs = vector(1, nnv); fprintf(ofp2, "rhs,%i\n", nnv);
	qvtemp = vector(1, nnv); fprintf(ofp2, "qvtemp,%i\n", nnv);

	ax = matrix(1, 3, 1, nnv); fprintf(ofp2, "ax,%i\n", nnv);
	cv = matrix(1, nnv, 1, nsp); fprintf(ofp2, "cv,%i\n", nnv);
	cv0 = matrix(1, nnv, 1, nsp); fprintf(ofp2, "cv0,%i\n", nnv);
	conv0 = matrix(1, nnv, 1, nsp); fprintf(ofp2, "conv0,%i\n", nnv);
	dcdp = matrix(1, nnv, 1, nsp); fprintf(ofp2, "dcdp,%i\n", nnv);
	gvv = matrix(1, nnv, 1, nnv); fprintf(ofp2, "gvv,%i\n", nnv);
	mat = dmatrix(1, nnv + 1, 1, nnv + 1); fprintf(ofp2, "mat,%i\n", nnv + 1);
	pev = matrix(1, nnv, 1, nsp); fprintf(ofp2, "pev,%i\n", nnv);
	pt = matrix(1, nnt, 1, nsp); fprintf(ofp2, "pt,%i\n", nnt);
	pt_old = matrix(1, nnt, 1, nsp); fprintf(ofp2, "pt_old,%i\n", nnt);
	ptprev = matrix(1, nnt, 1, nsp); fprintf(ofp2, "ptprev,%i\n", nnt);
	ptv = matrix(1, nnt, 1, nsp); fprintf(ofp2, "ptv,%i\n", nnt);
	pv = matrix(1, nnv, 1, nsp); fprintf(ofp2, "pv,%i\n", nnv);
	pvprev = matrix(1, nnv, 1, nsp); fprintf(ofp2, "pvprev,%i\n", nnv);
	pvt = matrix(1, nnv, 1, nsp); fprintf(ofp2, "pvt,%i\n", nnv);
	qcoeff1 = matrix(1, nnt, 1, nsp); fprintf(ofp2, "qcoeff1,%i\n", nnt);
	qt = matrix(1, nnt, 1, nsp); fprintf(ofp2, "qt,%i\n", nnt);
	qt_old = matrix(1, nnt, 1, nsp); fprintf(ofp2, "qt_old,%i\n", nnt);
	qv = matrix(1, nnv, 1, nsp); fprintf(ofp2, "qv,%i\n", nnv);
	qvprev = matrix(1, nnv, 1, nsp); fprintf(ofp2, "qvprev,%i\n", nnv);
	rhsg = dmatrix(1, nnv + 1, 1, 2); fprintf(ofp2, "rhsg,%i\n", nnv + 1);
	rhsl = dvector(1, nnv + 1); fprintf(ofp2, "rhsl,%i\n", nnv + 1);
	matx = dvector(1, nnv + 1); fprintf(ofp2, "matx,%i\n", nnv + 1);
	al = matrix(1, nnv, 1, nnv); fprintf(ofp2, "al,%i\n", nnv);
	activegr = matrix(1, actnum_temp, 1, 6); fprintf(ofp2, "activegr,%i\n", actnum_temp);

	fclose(ofp2);
}
