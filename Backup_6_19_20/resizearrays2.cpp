/************************************************************************
resizearrays2 - for AngioAdapt07.  TWS January 08
Adjust arrays with dimensions nnv and nnt
Transfer data for qt, pt to pt_old, qt_old
to initialize next step
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void resizearrays2(int nnv, int nnt, int actnum)
{
	extern int nsp,nnod,nnt_old_dim,nnv_dim,nnt_dim,actnum_dim;
	extern int *indx,*mainseg,*activenode;
	extern int **tisspoints;
	extern float *rhs,*qvtemp,*dtmin, * sourcefac, * sourcefacGF;
	extern float **qv,**pv,**pev,**pvt,**pvprev,**qvprev,**cv,**dcdp,**cv0,**conv0,**al,**gvv;
	extern float **qt,**pt,**qt_old,**pt_old,**ptprev,**ptv,**qcoeff1,**ax;
	extern double *rhsl,*matx;
	extern double **rhsg,**mat;

//used for angio
	extern int *new_segnum,*new_segtyp,*new_active_segs;
	extern float **activegr;

	int i,isp;
	FILE *ofp2;
	ofp2 = fopen("ArrayChanges.out","a");

//save old values for use as new initial values
	if(nnt_dim > nnt_old_dim){
		free_matrix(qt_old,1,nnt_old_dim,1,nsp); fprintf(ofp2,"free_qt_old,%i\n",nnt_old_dim);
		qt_old = matrix(1,nnt_dim,1,nsp); fprintf(ofp2,"qt_old,%i\n",nnt_dim);

		free_matrix(pt_old,1,nnt_old_dim,1,nsp); fprintf(ofp2,"free_pt_old,%i\n",nnt_old_dim);
		pt_old = matrix(1,nnt_dim,1,nsp); fprintf(ofp2,"pt_old,%i\n",nnt_dim);

		nnt_old_dim = nnt_dim;
	}
	for(isp=1; isp<=nsp; isp++) for(i=1; i<=nnt_dim; i++){
		qt_old[i][isp] = qt[i][isp];
		pt_old[i][isp] = pt[i][isp];
	}

//update nnv arrays
	if(nnv > nnv_dim){
		free_matrix(al,1,nnv_dim,1,nnv_dim); fprintf(ofp2,"free_al,%i\n",nnv_dim);	//TWS2010
		al = matrix(1,nnv,1,nnv); fprintf(ofp2,"al,%i\n",nnv);

		free_ivector(indx,1,nnv_dim); fprintf(ofp2,"free_indx,%i\n",nnv_dim);
		indx = ivector(1,nnv); fprintf(ofp2,"indx,%i\n",nnv);

		free_ivector(mainseg,1,nnv_dim); fprintf(ofp2,"free_mainseg,%i\n",nnv_dim);
		mainseg = ivector(1,nnv); fprintf(ofp2,"mainseg,%i\n",nnv);

		free_vector(rhs,1,nnv_dim); fprintf(ofp2,"free_rhs,%i\n",nnv_dim);
		rhs = vector(1,nnv); fprintf(ofp2,"rhs,%i\n",nnv);

		free_vector(qvtemp,1,nnv_dim); fprintf(ofp2,"free_qvtemp,%i\n",nnv_dim);
		qvtemp = vector(1,nnv); fprintf(ofp2,"qvtemp,%i\n",nnv);

		free_dvector(rhsl,1,nnv_dim+1); fprintf(ofp2,"free_rhsl,%i\n",nnv_dim+1);
		rhsl = dvector(1,nnv+1); fprintf(ofp2,"rhsl,%i\n",nnv+1);

		free_dvector(matx,1,nnv_dim+1); fprintf(ofp2,"free_matx,%i\n",nnv_dim+1);
		matx = dvector(1,nnv+1); fprintf(ofp2,"matx,%i\n",nnv+1);

		free_matrix(ax,1,3,1,nnv_dim); fprintf(ofp2,"free_ax,%i\n",nnv_dim);
		ax = matrix(1,3,1,nnv); fprintf(ofp2,"ax,%i\n",nnv);

		free_matrix(cv,1,nnv_dim,1,nsp); fprintf(ofp2,"free_cv,%i\n",nnv_dim);
		cv = matrix(1,nnv,1,nsp); fprintf(ofp2,"cv,%i\n",nnv);

		free_matrix(cv0,1,nnv_dim,1,nsp); fprintf(ofp2,"free_cv0,%i\n",nnv_dim);
		cv0 = matrix(1,nnv,1,nsp); fprintf(ofp2,"cv0,%i\n",nnv);

		free_matrix(conv0,1,nnv_dim,1,nsp); fprintf(ofp2,"free_conv0,%i\n",nnv_dim);
		conv0 = matrix(1,nnv,1,nsp); fprintf(ofp2,"conv0,%i\n",nnv);

		free_matrix(dcdp,1,nnv_dim,1,nsp); fprintf(ofp2,"free_dcdp,%i\n",nnv_dim);
		dcdp = matrix(1,nnv,1,nsp); fprintf(ofp2,"dcdp,%i\n",nnv);

		free_matrix(gvv,1,nnv_dim,1,nnv_dim); fprintf(ofp2,"free_gvv,%i\n",nnv_dim);
		gvv = matrix(1,nnv,1,nnv); fprintf(ofp2,"gvv,%i\n",nnv);

		free_dmatrix(mat,1,nnv_dim+1,1,nnv_dim+1); fprintf(ofp2,"free_mat,%i\n",nnv_dim+1);
		mat = dmatrix(1,nnv+1,1,nnv+1); fprintf(ofp2,"mat,%i\n",nnv+1);

		free_matrix(pev,1,nnv_dim,1,nsp); fprintf(ofp2,"free_pev,%i\n",nnv_dim);
		pev = matrix(1,nnv,1,nsp); fprintf(ofp2,"pev,%i\n",nnv);

		free_matrix(pv,1,nnv_dim,1,nsp); fprintf(ofp2,"free_pv,%i\n",nnv_dim);
		pv = matrix(1,nnv,1,nsp); fprintf(ofp2,"pv,%i\n",nnv);

		free_matrix(pvt,1,nnv_dim,1,nsp); fprintf(ofp2,"free_pvt,%i\n",nnv_dim);
		pvt = matrix(1,nnv,1,nsp); fprintf(ofp2,"pvt,%i\n",nnv);

		free_matrix(pvprev,1,nnv_dim,1,nsp); fprintf(ofp2,"free_pvprev,%i\n",nnv_dim);
		pvprev = matrix(1,nnv,1,nsp); fprintf(ofp2,"pvprev,%i\n",nnv);

		free_matrix(qv,1,nnv_dim,1,nsp); fprintf(ofp2,"free_qv,%i\n",nnv_dim);
		qv = matrix(1,nnv,1,nsp); fprintf(ofp2,"qv,%i\n",nnv);

		free_matrix(qvprev,1,nnv_dim,1,nsp); fprintf(ofp2,"free_qvprev,%i\n",nnv_dim);
		qvprev = matrix(1,nnv,1,nsp); fprintf(ofp2,"qvprev,%i\n",nnv);

		free_dmatrix(rhsg,1,nnv_dim+1,1,2); fprintf(ofp2,"free_rhsg,%i\n",nnv_dim+1);
		rhsg = dmatrix(1,nnv+1,1,2); fprintf(ofp2,"rhsg,%i\n",nnv+1);

		nnv_dim = nnv;
	}
//update nnt arrays
	if(nnt > nnt_dim){
		free_imatrix(tisspoints,1,3,1,nnt_dim); fprintf(ofp2,"free_tisspoints,%i\n",nnt_dim);
		tisspoints = imatrix(1,3,1,nnt); fprintf(ofp2,"tisspoints,%i\n",nnt);

		free_matrix(pt,1,nnt_dim,1,nsp); fprintf(ofp2,"free_pt,%i\n",nnt_dim);
		pt = matrix(1,nnt,1,nsp); fprintf(ofp2,"pt,%i\n",nnt);

		free_matrix(ptprev,1,nnt_dim,1,nsp); fprintf(ofp2,"free_ptprev,%i\n",nnt_dim);
		ptprev = matrix(1,nnt,1,nsp); fprintf(ofp2,"ptprev,%i\n",nnt);

		free_matrix(ptv,1,nnt_dim,1,nsp); fprintf(ofp2,"free_ptv,%i\n",nnt_dim);
		ptv = matrix(1,nnt,1,nsp); fprintf(ofp2,"ptv,%i\n",nnt);

		free_matrix(qcoeff1,1,nnt_dim,1,nsp); fprintf(ofp2,"free_qcoeff1,%i\n",nnt_dim);
		qcoeff1 = matrix(1,nnt,1,nsp); fprintf(ofp2,"qcoeff1,%i\n",nnt);

		free_matrix(qt,1,nnt_dim,1,nsp); fprintf(ofp2,"free_qt,%i\n",nnt_dim);
		qt = matrix(1,nnt,1,nsp); fprintf(ofp2,"qt,%i\n",nnt);

		free_vector(dtmin, 1, nnt_dim); fprintf(ofp2, "free_dtmin,%i\n", nnt_dim);
		dtmin = vector(1, nnt); fprintf(ofp2, "dtmin,%i\n", nnt);

		free_vector(sourcefac, 1, nnt_dim); fprintf(ofp2, "free_sourcefac,%i\n", nnt_dim);
		sourcefac = vector(1, nnt); fprintf(ofp2, "sourcefac,%i\n", nnt);

		free_vector(sourcefacGF, 1, nnt_dim); fprintf(ofp2, "free_sourcefacGF,%i\n", nnt_dim);
		sourcefacGF = vector(1, nnt); fprintf(ofp2, "sourcefacGF,%i\n", nnt);

		nnt_dim = nnt;

		 
	}
//update actnum arrays
	if(actnum > actnum_dim){
		free_ivector(activenode,1,actnum_dim); fprintf(ofp2,"free_activenode,%i\n",actnum_dim);
		activenode = ivector(1,actnum); fprintf(ofp2,"activenode,%i\n",actnum);

		free_matrix(activegr,1,actnum_dim,1,6); fprintf(ofp2,"free_activegr,%i\n",actnum_dim);
		activegr = matrix(1,actnum,1,6); fprintf(ofp2,"activegr,%i\n",actnum);

		actnum_dim = actnum;
	}
	fclose(ofp2);
}
