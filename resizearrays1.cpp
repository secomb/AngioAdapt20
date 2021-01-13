/************************************************************************
resizearrays1 - for AngioAdapt07.  TWS January 08
Adjust arrays to new dimensions nnod and nseg
Transfer old data to new arrays for following variables:
ista, iend, nodname, segname, segtyp, segnodname
diam, ksseg, q, hd, cnode, nodpress
Requires temporary arrays for these
Also save values of istart
Revised TWS2010
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void resizearrays1(int nseg, int nnod)
{
	extern int nsp, nodsegm, nspint, nseg_dim, nnod_dim;
	extern int* nspoint, * istart, * lowflow, * nk, * nodrank, * nodout, * nodtyp, * nodclasstyp, * segclasstyp;
	extern int* ista, * iend, * segname, * segtyp, * nodname, * boundseg;
	extern int** nodnod, ** nodseg, ** segnodname;
	extern float* diam, * ksseg, * q, * hd, * segvar, * qq, * qqhd, * rseg, * cbar, * lseg, * ds, * nodvar, * qold, * hdold, * cond;
	extern float* segpress, * tau, * condup, * conddown, * metsig, * segc;
	extern float* stot, * stau, * spress, * uptrans, * downtrans;
	extern float** cnode, ** gamma1, ** qvseg, ** pvseg, ** pevseg, ** scos, ** signalvar;
	extern double* nodpress;

	int iseg, inod, i, * ista_t, * iend_t, * nodname_t, * segname_t, * segtyp_t, ** segnodname_t, * boundseg_t;
	float* diam_t, * ds_t, * ksseg_t, * q_t, * hd_t, ** cnode_t;
	double* nodpress_t;

	FILE* ofp2;
	ofp2 = fopen("ArrayChanges.out", "a");

	//resize nseg arrays
	if (nseg > nseg_dim) {
		//set up temporary arrays
		boundseg_t = ivector(1, nseg_dim); fprintf(ofp2, "boundseg_t,%i\n", nseg_dim);
		ista_t = ivector(1, nseg_dim); fprintf(ofp2, "ista_t,%i\n", nseg_dim);
		iend_t = ivector(1, nseg_dim); fprintf(ofp2, "iend_t,%i\n", nseg_dim);
		segname_t = ivector(1, nseg_dim); fprintf(ofp2, "segname_t,%i\n", nseg_dim);
		segtyp_t = ivector(1, nseg_dim); fprintf(ofp2, "segtyp_t,%i\n", nseg_dim);
		segnodname_t = imatrix(1, 2, 1, nseg_dim); fprintf(ofp2, "segnodname_t,%i\n", nseg_dim);
		diam_t = vector(1, nseg_dim); fprintf(ofp2, "diam_t,%i\n", nseg_dim);
		ds_t = vector(1, nseg_dim); fprintf(ofp2, "ds_t,%i\n", nseg_dim);
		ksseg_t = vector(1, nseg_dim); fprintf(ofp2, "ksseg_t,%i\n", nseg_dim);
		q_t = vector(1, nseg_dim); fprintf(ofp2, "q_t,%i\n", nseg_dim);
		hd_t = vector(1, nseg_dim); fprintf(ofp2, "hd_t,%i\n", nseg_dim);
		//save variables to be transferred to next run
		for (iseg = 1; iseg <= nseg_dim; iseg++) {
			boundseg_t[iseg] = boundseg[iseg];
			ista_t[iseg] = ista[iseg];
			iend_t[iseg] = iend[iseg];
			segname_t[iseg] = segname[iseg];
			segtyp_t[iseg] = segtyp[iseg];
			for (i = 1; i <= 2; i++) segnodname_t[i][iseg] = segnodname[i][iseg];
			diam_t[iseg] = diam[iseg];
			ds_t[iseg] = ds[iseg];
			ksseg_t[iseg] = ksseg[iseg];
			q_t[iseg] = q[iseg];
			hd_t[iseg] = hd[iseg];
		}
		free_ivector(boundseg, 1, nseg_dim); fprintf(ofp2, "free_boundseg,%i\n", nseg_dim);
		boundseg = ivector(1, nseg); fprintf(ofp2, "boundseg,%i\n", nseg);

		free_ivector(iend, 1, nseg_dim); fprintf(ofp2, "free_iend,%i\n", nseg_dim);
		iend = ivector(1, nseg); fprintf(ofp2, "iend,%i\n", nseg);

		free_ivector(ista, 1, nseg_dim); fprintf(ofp2, "free_ista,%i\n", nseg_dim);
		ista = ivector(1, nseg); fprintf(ofp2, "ista,%i\n", nseg);

		free_ivector(istart, 1, nseg_dim); fprintf(ofp2, "free_istart,%i\n", nseg_dim);
		istart = ivector(1, nseg); fprintf(ofp2, "istart,%i\n", nseg);

		free_ivector(lowflow, 1, nseg_dim); fprintf(ofp2, "free_lowflow,%i\n", nseg_dim);
		lowflow = ivector(1, nseg); fprintf(ofp2, "lowflow,%i\n", nseg);

		free_ivector(nspoint, 1, nseg_dim); fprintf(ofp2, "free_nspoint,%i\n", nseg_dim);
		nspoint = ivector(1, nseg); fprintf(ofp2, "nspoint,%i\n", nseg);

		free_ivector(segname, 1, nseg_dim); fprintf(ofp2, "free_segname,%i\n", nseg_dim);
		segname = ivector(1, nseg); fprintf(ofp2, "segname,%i\n", nseg);

		free_ivector(segtyp, 1, nseg_dim); fprintf(ofp2, "free_segtyp,%i\n", nseg_dim);
		segtyp = ivector(1, nseg); fprintf(ofp2, "segtyp,%i\n", nseg);

		free_ivector(segclasstyp, 1, nseg_dim); fprintf(ofp2, "free_segclasstyp,%i\n", nseg_dim);
		segclasstyp = ivector(1, nseg); fprintf(ofp2, "segclasstyp,%i\n", nseg);

		free_imatrix(segnodname, 1, 2, 1, nseg_dim); fprintf(ofp2, "free_segnodname,%i\n", nseg_dim);
		segnodname = imatrix(1, 2, 1, nseg); fprintf(ofp2, "segnodname,%i\n", nseg);

		free_vector(cbar, 1, nseg_dim); fprintf(ofp2, "free_cbar,%i\n", nseg_dim);
		cbar = vector(1, nseg); fprintf(ofp2, "cbar,%i\n", nseg);

		free_vector(segc, 1, nseg_dim); fprintf(ofp2, "free_segc,%i\n", nseg_dim);
		segc = vector(1, nseg); fprintf(ofp2, "segc,%i\n", nseg);//TWS2010

		free_vector(cond, 1, nseg_dim); fprintf(ofp2, "free_cond,%i\n", nseg_dim);
		cond = vector(1, nseg); fprintf(ofp2, "cond,%i\n", nseg);

		free_vector(condup, 1, nseg_dim); fprintf(ofp2, "free_condup,%i\n", nseg_dim);
		condup = vector(1, nseg); fprintf(ofp2, "condup,%i\n", nseg);

		free_vector(conddown, 1, nseg_dim); fprintf(ofp2, "free_conddown,%i\n", nseg_dim);
		conddown = vector(1, nseg); fprintf(ofp2, "conddown,%i\n", nseg);

		free_vector(diam, 1, nseg_dim); fprintf(ofp2, "free_diam,%i\n", nseg_dim);
		diam = vector(1, nseg); fprintf(ofp2, "diam,%i\n", nseg);

		free_vector(ksseg, 1, nseg_dim); fprintf(ofp2, "free_ksseg,%i\n", nseg_dim);
		ksseg = vector(1, nseg); fprintf(ofp2, "ksseg,%i\n", nseg);

		free_vector(downtrans, 1, nseg_dim); fprintf(ofp2, "free_downtrans,%i\n", nseg_dim);
		downtrans = vector(1, nseg); fprintf(ofp2, "downtrans,%i\n", nseg);

		free_vector(ds, 1, nseg_dim); fprintf(ofp2, "free_ds,%i\n", nseg_dim);
		ds = vector(1, nseg); fprintf(ofp2, "ds,%i\n", nseg);

		free_vector(hd, 1, nseg_dim); fprintf(ofp2, "free_hd,%i\n", nseg_dim);
		hd = vector(1, nseg); fprintf(ofp2, "hd,%i\n", nseg);

		free_vector(hdold, 1, nseg_dim); fprintf(ofp2, "free_hdold,%i\n", nseg_dim);
		hdold = vector(1, nseg); fprintf(ofp2, "hdold,%i\n", nseg);

		free_vector(lseg, 1, nseg_dim); fprintf(ofp2, "free_lseg,%i\n", nseg_dim);
		lseg = vector(1, nseg); fprintf(ofp2, "lseg,%i\n", nseg);

		free_vector(metsig, 1, nseg_dim); fprintf(ofp2, "free_metsig,%i\n", nseg_dim);
		metsig = vector(1, nseg); fprintf(ofp2, "metsig,%i\n", nseg);

		free_vector(q, 1, nseg_dim); fprintf(ofp2, "free_q,%i\n", nseg_dim);
		q = vector(1, nseg); fprintf(ofp2, "q,%i\n", nseg);

		free_vector(qq, 1, nseg_dim); fprintf(ofp2, "free_qq,%i\n", nseg_dim);
		qq = vector(1, nseg); fprintf(ofp2, "qq,%i\n", nseg);

		free_vector(qqhd, 1, nseg_dim); fprintf(ofp2, "free_qqhd,%i\n", nseg_dim);
		qqhd = vector(1, nseg); fprintf(ofp2, "qqhd,%i\n", nseg);

		free_vector(qold, 1, nseg_dim); fprintf(ofp2, "free_qold,%i\n", nseg_dim);
		qold = vector(1, nseg); fprintf(ofp2, "qold,%i\n", nseg);

		free_vector(rseg, 1, nseg_dim); fprintf(ofp2, "free_rseg,%i\n", nseg_dim);
		rseg = vector(1, nseg); fprintf(ofp2, "rseg,%i\n", nseg);

		free_vector(segvar, 1, nseg_dim); fprintf(ofp2, "free_segvar,%i\n", nseg_dim);
		segvar = vector(1, nseg); fprintf(ofp2, "segvar,%i\n", nseg);

		free_vector(segpress, 1, nseg_dim); fprintf(ofp2, "free_segpress,%i\n", nseg_dim);
		segpress = vector(1, nseg); fprintf(ofp2, "segpress,%i\n", nseg);

		free_vector(spress, 1, nseg_dim); fprintf(ofp2, "free_spress,%i\n", nseg_dim);
		spress = vector(1, nseg); fprintf(ofp2, "spress,%i\n", nseg);

		free_vector(stot, 1, nseg_dim); fprintf(ofp2, "free_stot,%i\n", nseg_dim);
		stot = vector(1, nseg); fprintf(ofp2, "stot,%i\n", nseg);

		free_vector(stau, 1, nseg_dim); fprintf(ofp2, "free_stau,%i\n", nseg_dim);
		stau = vector(1, nseg); fprintf(ofp2, "stau,%i\n", nseg);

		free_vector(tau, 1, nseg_dim); fprintf(ofp2, "free_tau,%i\n", nseg_dim);
		tau = vector(1, nseg); fprintf(ofp2, "tau,%i\n", nseg);

		free_vector(uptrans, 1, nseg_dim); fprintf(ofp2, "free_uptrans,%i\n", nseg_dim);
		uptrans = vector(1, nseg); fprintf(ofp2, "uptrans,%i\n", nseg);

		free_matrix(signalvar, 1, nseg_dim, 1, 6); fprintf(ofp2, "free_signalvar,%i\n", nseg_dim);
		signalvar = matrix(1, nseg, 1, 6); fprintf(ofp2, "signalvar,%i\n", nseg);

		free_matrix(gamma1, 1, nseg_dim, 1, nsp); fprintf(ofp2, "free_gamma1,%i\n", nseg_dim);
		gamma1 = matrix(1, nseg, 1, nsp); fprintf(ofp2, "gamma1,%i\n", nseg);

		free_matrix(pvseg, 1, nseg_dim, 1, nsp); fprintf(ofp2, "free_pvseg,%i\n", nseg_dim);
		pvseg = matrix(1, nseg, 1, nsp); fprintf(ofp2, "pvseg,%i\n", nseg);

		free_matrix(pevseg, 1, nseg_dim, 1, nsp); fprintf(ofp2, "free_pevseg,%i\n", nseg_dim);
		pevseg = matrix(1, nseg, 1, nsp); fprintf(ofp2, "pevseg,%i\n", nseg);

		free_matrix(qvseg, 1, nseg_dim, 1, nsp); fprintf(ofp2, "free_qvseg,%i\n", nseg_dim);
		qvseg = matrix(1, nseg, 1, nsp); fprintf(ofp2, "qvseg,%i\n", nseg);

		free_matrix(scos, 1, 3, 1, nseg_dim); fprintf(ofp2, "free_scos,%i\n", nseg_dim);
		scos = matrix(1, 3, 1, nseg); fprintf(ofp2, "scos,%i\n", nseg);

		//transfer information into new arrays
		for (iseg = 1; iseg <= IMIN(nseg_dim, nseg); iseg++) {
			boundseg[iseg] = boundseg_t[iseg];
			ista[iseg] = ista_t[iseg];
			iend[iseg] = iend_t[iseg];
			segname[iseg] = segname_t[iseg];
			segtyp[iseg] = segtyp_t[iseg];
			diam[iseg] = diam_t[iseg];
			ds[iseg] = ds_t[iseg];
			ksseg[iseg] = ksseg_t[iseg];
			q[iseg] = q_t[iseg];
			hd[iseg] = hd_t[iseg];
			for (i = 1; i <= 2; i++) segnodname[i][iseg] = segnodname_t[i][iseg];
		}
		//free temporary arrays
		free_ivector(boundseg_t, 1, nseg_dim); fprintf(ofp2, "free_boundseg_t,%i\n", nseg_dim);
		free_ivector(ista_t, 1, nseg_dim); fprintf(ofp2, "free_ista_t,%i\n", nseg_dim);
		free_ivector(iend_t, 1, nseg_dim); fprintf(ofp2, "free_iend_t,%i\n", nseg_dim);
		free_ivector(segname_t, 1, nseg_dim); fprintf(ofp2, "free_segname_t,%i\n", nseg_dim);
		free_ivector(segtyp_t, 1, nseg_dim); fprintf(ofp2, "free_segtyp_t,%i\n", nseg_dim);
		free_imatrix(segnodname_t, 1, 2, 1, nseg_dim); fprintf(ofp2, "free_segnodname_t,%i\n", nseg_dim);
		free_vector(diam_t, 1, nseg_dim); fprintf(ofp2, "free_diam_t,%i\n", nseg_dim);
		free_vector(ds_t, 1, nseg_dim); fprintf(ofp2, "free_ds_t,%i\n", nseg_dim);
		free_vector(ksseg_t, 1, nseg_dim); fprintf(ofp2, "free_ksseg_t,%i\n", nseg_dim);
		free_vector(q_t, 1, nseg_dim); fprintf(ofp2, "free_q_t,%i\n", nseg_dim);
		free_vector(hd_t, 1, nseg_dim); fprintf(ofp2, "free_hd_t,%i\n", nseg_dim);
		nseg_dim = nseg;
	}

	//resize nnod arrays
	if (nnod > nnod_dim) {
		//set up temporary arrays
		nodname_t = ivector(1, nnod_dim); fprintf(ofp2, "nodname_t,%i\n", nnod_dim);
		cnode_t = matrix(1, 3, 1, nnod_dim); fprintf(ofp2, "cnode_t,%i\n", nnod_dim);
		nodpress_t = dvector(1, nnod_dim); fprintf(ofp2, "nodpress_t,%i\n", nnod_dim);
		for (inod = 1; inod <= nnod_dim; inod++) {
			nodname_t[inod] = nodname[inod];
			nodpress_t[inod] = nodpress[inod];
			for (i = 1; i <= 3; i++) cnode_t[i][inod] = cnode[i][inod];
		}
		free_ivector(nk, 1, nnod_dim); fprintf(ofp2, "free_nk,%i\n", nnod_dim);
		nk = ivector(1, nnod); fprintf(ofp2, "nk,%i\n", nnod);

		free_ivector(nodname, 1, nnod_dim); fprintf(ofp2, "free_nodname,%i\n", nnod_dim);
		nodname = ivector(1, nnod); fprintf(ofp2, "nodname,%i\n", nnod);

		free_ivector(nodout, 1, nnod_dim); fprintf(ofp2, "free_nodout,%i\n", nnod_dim);
		nodout = ivector(1, nnod); fprintf(ofp2, "nodout,%i\n", nnod);

		free_ivector(nodrank, 1, nnod_dim); fprintf(ofp2, "free_nodrank,%i\n", nnod_dim);
		nodrank = ivector(1, nnod); fprintf(ofp2, "nodrank,%i\n", nnod);

		free_ivector(nodtyp, 1, nnod_dim); fprintf(ofp2, "free_nodtyp,%i\n", nnod_dim);
		nodtyp = ivector(1, nnod); fprintf(ofp2, "nodtyp,%i\n", nnod);

		free_ivector(nodclasstyp, 1, nnod_dim); fprintf(ofp2, "free_nodclasstyp,%i\n", nnod_dim);
		nodclasstyp = ivector(1, nnod); fprintf(ofp2, "nodclasstyp,%i\n", nnod);

		free_imatrix(nodnod, 1, nodsegm, 1, nnod_dim); fprintf(ofp2, "free_nodnod,%i\n", nnod_dim);
		nodnod = imatrix(1, nodsegm, 1, nnod); fprintf(ofp2, "nodnod,%i\n", nnod);

		free_imatrix(nodseg, 1, nodsegm, 1, nnod_dim); fprintf(ofp2, "free_nodseg,%i\n", nnod_dim);
		nodseg = imatrix(1, nodsegm, 1, nnod); fprintf(ofp2, "nodseg,%i\n", nnod);

		free_vector(nodvar, 1, nnod_dim); fprintf(ofp2, "free_nodvar,%i\n", nnod_dim);
		nodvar = vector(1, nnod); fprintf(ofp2, "nodvar,%i\n", nnod);

		free_dvector(nodpress, 1, nnod_dim); fprintf(ofp2, "free_nodpress,%i\n", nnod_dim);
		nodpress = dvector(1, nnod); fprintf(ofp2, "nodpress,%i\n", nnod);

		free_matrix(cnode, 1, 3, 1, nnod_dim); fprintf(ofp2, "free_cnode,%i\n", nnod_dim);
		cnode = matrix(1, 3, 1, nnod); fprintf(ofp2, "cnode,%i\n", nnod);

		//transfer information into new arrays
		for (inod = 1; inod <= nnod_dim; inod++) {
			nodname[inod] = nodname_t[inod];
			nodpress[inod] = nodpress_t[inod];
			for (i = 1; i <= 3; i++) cnode[i][inod] = cnode_t[i][inod];
		}
		//free temporary arrays
		free_ivector(nodname_t, 1, nnod_dim); fprintf(ofp2, "free_nodname_t,%i\n", nnod_dim);
		free_dvector(nodpress_t, 1, nnod_dim); fprintf(ofp2, "free_nodpress_t,%i\n", nnod_dim);
		free_matrix(cnode_t, 1, 3, 1, nnod_dim); fprintf(ofp2, "free_cnode_t,%i\n", nnod_dim);
		nnod_dim = nnod;
	}
	fclose(ofp2);
}
