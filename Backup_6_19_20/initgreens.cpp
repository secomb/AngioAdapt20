/****************************************************************************
initgreens - initial tissue source strengths, given uniform solute field, initial g0
initial vessel source strengths based on uniform efflux rate from all vessels
if values are not available from previous call to greens
TWS Jan 08
Version 2.0, May 1, 2010.
Version 3.0, May 17, 2011.
*******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void tissrate(int nsp, float *p, float *mtiss, float *mptiss);

void initgreens()
{
	extern int nnt,nnt_dim,nnv,nsp,imain,mxx,myy,mzz,nseg;
	extern int *lowflow,*permsolute,*nspoint,*segtyp,*istart;
	extern int ***nbou,***nbou_old;
	extern int *oxygen,*mainseg; //added TWS2010

	extern float vol,errfac,tlength;

	extern float *mtiss,*mptiss,*epsvessel,*epstissue,*eps,*errvessel,*errtissue,*pinit,*p;
	extern float *qq,*g0,*qtsum,*pref,*ds,*lseg,*hd;
	extern float **qt,**qt_old,**qv,**pt,**pt_old,**pv,**tissparam;

	int isp,i,j,k,itp,itp_old,iseg;
	float tlengthq,tlengthqhd;

	for(isp=1; isp<=nsp; isp++) pinit[isp] = FMIN(g0[isp],pref[isp]);//changed 7/2008 as epstissue[2] was too low
	tissrate(nsp,pinit,mtiss,mptiss);
	for(isp=1; isp<=nsp; isp++){
		qtsum[isp] = 0.;
		for(i=1; i<=mxx; i++) for(j=1; j<=myy; j++) for(k=1; k<=mzz; k++){
			itp = nbou[i][j][k];
			itp_old = nbou_old[i][j][k];
//add  nnt_dim to extern, then use if to see itp,itpold > nnt,nnt_dim
			if(itp > 0){
				if(imain == 1 || itp_old == 0){
					qt[itp][isp] = mtiss[isp]*vol;
					pt[itp][isp] = pinit[isp];
				}
				else{
					qt[itp][isp] = qt_old[itp_old][isp];
					pt[itp][isp] = pt_old[itp_old][isp];
				}
				qtsum[isp] += qt[itp][isp];
			}
		}
	}
	tlength = 0.;
	tlengthq = 0.;
	tlengthqhd = 0.;
	for(iseg=1; iseg<=nseg; iseg++) if(segtyp[iseg] >= 3 && segtyp[iseg] <= 5){
		tlength += lseg[iseg];
		tlengthq += lseg[iseg]*qq[iseg];
		tlengthqhd += lseg[iseg]*qq[iseg]*hd[iseg];//added 8/10
	}
	for(isp=1; isp<=nsp; isp++){
		for(i=1; i<=nnv; i++){
			qv[i][isp] = 0.;
			pv[i][isp] = 0.;
			if(permsolute[isp] == 1){
				if(oxygen[isp] == 1){
					qv[i][isp] = -qtsum[isp]*ds[mainseg[i]]*qq[mainseg[i]]*hd[mainseg[i]]/tlengthqhd;//modified 8/10
					if(lowflow[mainseg[i]] == 1) qv[i][isp] = 0.;  //low q*hd
				}
				else qv[i][isp] = -qtsum[isp]*ds[mainseg[i]]*qq[mainseg[i]]/tlengthq;//modified 8/09
				pv[i][isp] = pinit[isp];
			}
		}
	}
//set error bounds, proportional to errfac
	for(isp=1; isp<=nsp; isp++){
		epsvessel[isp] = fabs(qtsum[isp])/nnv*errfac;
		eps[isp] = pref[isp]*errfac;
		epstissue[isp] = tissparam[1][isp]*vol*errfac;  //June 2009 - requires tissparam[1][isp] to be max rate
	}
}