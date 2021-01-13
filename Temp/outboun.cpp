/****************************************************
outboun.cpp
----------------------------------------------------
outbounmethod = 1: finds the smallest convex region inside the cuboid
which excludes tissue node points that have a distance to the nearest
vessel greater than a value specified by the user (lb).  Any point
outside a region between two planes containing the required points is
excluded.  This is repeated for multiple plane orientations, determined by am.
----------------------------------------------------
outbounmethod = 2: finds all tissue points within a distance lb of the vessels, but
does not make a convex region.  Fills in 'holes' in interior, whether 2D or 3D.
----------------------------------------------------
outbounmethod = 3: as in outbounmethod 1, but remove a strip at top and bottom of rectangular domain
----------------------------------------------------
outbounmethod = 4: retinal angiogenesis model. Consider a quadrant of an annulus
with inner radius fixed, but with outer radius increasing to a maximum based on most
distant vessel from origin plus outer bound distance . 
--------------------------------------------------------------------------------------------------------
outbounmethod = 5: retinal angiogenesis model, full disk centered on domain.
As in outbounmethod 2, but exclude regions r < IR and r > OR
----------------------------------------------------
Outbounmethod = 6: brain angiogenesis model, cylinder in z-direction
With rounded rims as sharp edges lead to hypoxia. October 2017
----------------------------------------------------
Outbounmethod = 7: brain angiogenesis model with periodic boundary conditions,
hexagonal cylinder in z-direction, February 2018
----------------------------------------------------
Output:	nnt, total tissue node points inside the region.
nbou > 1 if inside region, value gives tissue point index
TWS 2010
Version 3.0, May 17, 2011.
******************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

int outboun(int outbounmethod)
{
	extern int nseg,mxx,myy,mzz,dogreensflag,num_runs_actual,*segtyp,*ista,*iend,***nbou;
	extern float lb,*axt,*ayt,*azt,alx,aly,alz,**scos,*lseg,*rseg,*midpt;
	extern float IR,OR;
	extern float aphex,aphexp,**cnode,**hex_norm;
	
	int i,j,k,iseg,jseg,a,b,c,am,nnt;
	int ii,jj,kk,imin,imax,jmin,jmax,kmin,kmax;
	float aa,bb,cc,abc,ddmin,ddmax,dd,t;
	float *tisspt, *orig;
	float ds2,de2,d2max,dmax,disp2,d,x11,x22,x33;
	float width6 = 50., radius6 = 200.;	//width of extra growth space, radius of chamfer for outbounmethod = 6

	tisspt = vector(1,3);
	orig = vector(1,3);

	if(outbounmethod == 1 || outbounmethod == 3){
		for(i=1; i<=mxx; i++) for(j=1; j<=myy; j++) for(k=1; k<=mzz; k++) nbou[i][j][k] = 1;
		am = 6;
		for(a=0; a<=am; a++) for(b=-am; b<=am; b++) for(c=-am; c<=am; c++){
			if(a != 0 || b != 0 || c != 0){
				aa = 0.;
				bb = 0.;
				cc = 0.;
				if(a != 0) aa = 1.0/a;
				if(b != 0) bb = 1.0/b;
				if(c != 0) cc = 1.0/c;
				abc = sqrt(SQR(aa) + SQR(bb) + SQR(cc));
				ddmin =  1.e8;
				ddmax = -1.e8;
				for(iseg=1; iseg<=nseg; iseg++)	if(segtyp[iseg] >= 3 && segtyp[iseg] <= 5){
					dd = (aa*cnode[1][ista[iseg]] + bb*cnode[2][ista[iseg]] + cc*cnode[3][ista[iseg]])/abc;
					ddmin = FMIN(dd - rseg[iseg] - lb,ddmin);
					ddmax = FMAX(dd + rseg[iseg] + lb,ddmax);
					dd = (aa*cnode[1][iend[iseg]] + bb*cnode[2][iend[iseg]] + cc*cnode[3][iend[iseg]])/abc;
					ddmin = FMIN(dd - rseg[iseg] - lb,ddmin);
					ddmax = FMAX(dd + rseg[iseg] + lb,ddmax);
				}
				for(i=1; i<=mxx; i++) for(j=1; j<=myy; j++) for(k=1; k<=mzz; k++){
					t = (aa*axt[i] + bb*ayt[j] + cc*azt[k])/abc;
					if(t > ddmax || t < ddmin) nbou[i][j][k] = 0;
				}
			}
		}
	}
	if(outbounmethod == 2 || outbounmethod == 5){
		for(i=1; i<=mxx; i++) for(j=1; j<=myy; j++) for(k=1; k<=mzz; k++){
			nbou[i][j][k] = 0;
			for(jseg=1; jseg<=nseg; jseg++) if(segtyp[jseg] >= 3 && segtyp[jseg] <= 5){
				x11 = (axt[i]-cnode[1][ista[jseg]])*scos[2][jseg]-(ayt[j]-cnode[2][ista[jseg]])*scos[1][jseg];
				x22 = (ayt[j]-cnode[2][ista[jseg]])*scos[3][jseg]-(azt[k]-cnode[3][ista[jseg]])*scos[2][jseg];
				x33 = (azt[k]-cnode[3][ista[jseg]])*scos[1][jseg]-(axt[i]-cnode[1][ista[jseg]])*scos[3][jseg];
				disp2 = SQR(x11) + SQR(x22) + SQR(x33);
				ds2 = SQR(axt[i]-cnode[1][ista[jseg]]) + SQR(ayt[j]-cnode[2][ista[jseg]]) + SQR(azt[k]-cnode[3][ista[jseg]]);
				de2 = SQR(axt[i]-cnode[1][iend[jseg]]) + SQR(ayt[j]-cnode[2][iend[jseg]]) + SQR(azt[k]-cnode[3][iend[jseg]]);
				if(FMAX(ds2,de2)-disp2 > SQR(lseg[jseg])) d = sqrt(FMIN(ds2,de2)) - rseg[jseg];
				else d = sqrt(disp2) - rseg[jseg];
				if(d < lb) nbou[i][j][k] = 1;
			}
		}
		for(i=1; i<=mxx; i++) for(j=1; j<=myy; j++) for(k=1; k<=mzz; k++) if(nbou[i][j][k] == 0){
			imin = mxx;
			imax = 1;
			for(ii=1; ii<=mxx; ii++) if(nbou[ii][j][k] == 1){
				imin = IMIN(ii,imin);
				imax = IMAX(ii,imax);
			}
			jmin = myy;
			jmax = 1;
			for(jj=1; jj<=myy; jj++) if(nbou[i][jj][k] == 1){
				jmin = IMIN(jj,jmin);
				jmax = IMAX(jj,jmax);
			}
			kmin = mzz;
			kmax = 1;
			for(kk=1; kk<=mzz; kk++) if(nbou[i][j][kk] == 1){
				kmin = IMIN(kk,kmin);
				kmax = IMAX(kk,kmax);
			}
			if(i >= imin && i <= imax) if(j >= jmin && j <= jmax) if(k >= kmin && k <= kmax) nbou[i][j][k] = 1;
		}
		if(outbounmethod == 5){
			for(i=1; i<=mxx; i++) for(j=1; j<=myy; j++) for(k=1; k<=mzz; k++){
				d = SQR(axt[i] - midpt[1]) + SQR(ayt[j] - midpt[2]);
				if(d < SQR(IR)) nbou[i][j][k] = 0;  		//exclude region near origin
				if(d > SQR(OR)) nbou[i][j][k] = 0;  		//exclude region at outer edge
			}
		}	
	}
	if(outbounmethod == 4){	//this method for the retinal angiogenesis problem with one-quarter disk
		d2max = 0.;
		for(iseg=1; iseg<=nseg; iseg++)	if(segtyp[iseg] != 10){
			ds2 = SQR(cnode[1][ista[iseg]]) + SQR(cnode[2][ista[iseg]]) + SQR(cnode[3][ista[iseg]]);
			if(ds2 > d2max) d2max = ds2;
			de2 = SQR(cnode[1][iend[iseg]]) + SQR(cnode[2][iend[iseg]]) + SQR(cnode[3][iend[iseg]]);
			if(de2 > d2max) d2max = de2;
		}
		dmax = FMIN(FMAX(sqrt(d2max),2000.),OR);	//choose distance to the most distant vessel plus 250
		for(i=1; i<=mxx; i++) for(j=1; j<=myy; j++) for(k=1; k<=mzz; k++){
			nbou[i][j][k] = 0;
			d = sqrt(SQR(axt[i]) + SQR(ayt[j]) + SQR(azt[k]));
			if(d < dmax && d > IR) nbou[i][j][k] = 1;
		}
	}
	if(outbounmethod == 6){	//this method for brain angiogenesis problem with a cylinder of tissue
		//width6 = 50., radius6 = 50.;	//width of extra growth space, radius of chamfer for outbounmethod = 6
		dmax = FMIN(alx,aly)/2. - width6;	//allow a region where vessels can grow outside tissue region
		for(i=1; i<=mxx; i++) for(j=1; j<=myy; j++) for(k=1; k<=mzz; k++){
			nbou[i][j][k] = 0;
			d = sqrt(SQR(axt[i] - midpt[1]) + SQR(ayt[j] - midpt[2]));
			if(d < dmax && azt[k] > width6){
				nbou[i][j][k] = 1;
				//remove tissue points to round off the edges of the cylindrical region
				if(azt[k] > alz - radius6 && SQR(d - (dmax - radius6)) + SQR(azt[k] - (alz - radius6)) > SQR(radius6)) nbou[i][j][k] = 0;
				if(azt[k] < width6 + radius6 && SQR(d - (dmax - radius6)) + SQR(width6 + radius6 - azt[k]) > SQR(radius6)) nbou[i][j][k] = 0;
			}
		}
	}
	if(outbounmethod == 7){	//this method for brain angiogenesis problem with periodic network
		for(i=1; i<=mxx; i++) for(j=1; j<=myy; j++) for(k=1; k<=mzz; k++){
			nbou[i][j][k] = 1;
			for(ii=1; ii<=6; ii++){
				d = (axt[i] - midpt[1])*hex_norm[ii][1] + (ayt[j] - midpt[2])*hex_norm[ii][2];
				if(d > aphexp) nbou[i][j][k] = 0;
			}
		}
	}
	nnt = 0;
	for(i=1; i<=mxx; i++) for(j=1; j<=myy; j++) for(k=1; k<=mzz; k++) if(nbou[i][j][k] == 1){
		nnt++;
		nbou[i][j][k] = nnt;	//assign nbou values from 1 to nnt
	}
	free_vector(tisspt,1,3);
	free_vector(orig,1,3);

	printf("Tissue node points used in simulation  = %i\n", nnt);
	printf("V of simulation region/V of cubic region = %f\n", nnt/(1.0*mxx*myy*mzz));
	return nnt;
}