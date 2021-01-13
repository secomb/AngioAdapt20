/**********************************************************
contourO2GF.cpp - generate data for contour plot.
Produces a single image with both oxygen and GF over threshold.
October 2010.
Hexagonal mesh added, July 2016
Multiple slice capability added February 2017
Modified for better results in 3D. Feb. 2017 and August 2017
***********************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void contr_shade(FILE *ofp, int m, int n, float scalefac, int nl,
		   float xmin, float xmax, float ymin, float ymax, float *cl, float **zv,
		   int showscale, int lowcolor, int hatch, int plotcontour);
float *eval(int slsegdiv, float req, float *x);
void evalGF(int isp, float req, float* x, float lambdaGF, float cfacGF, float ssGF, float bbGF, float ccGF,
	float* GF, float* GFx, float* GFy, float* GFz);
void evalVAF(int isp, float req, float* x, float lambdaGF, float cfacGF, float ssGF, float bbGF, float ccGF,
	float* GF, float* GFx, float* GFy, float* GFz);

void contourO2GF(int imain)
{
	extern int max, nsp, nseg, * segtyp, * nl, * ista, * iend;
	extern int* slsegdiv, * nsl1, * nsl2, ncontourplanes;	//needed for multiple contour plots per page
	extern int hexperiodic, * boundseg;
	extern float pi1, req, thresh;
	extern float lambdaGF2, cfacGF2, ssGF2, bbGF2, ccGF2;
	extern float lambdaGF3, cfacGF3, ssGF3, bbGF3, ccGF3;
	extern float* x, * y, * p, * diam, ** cnode, ** pvseg, * midpt;
	extern float** xsl0, ** xsl1, ** xsl2, * clmin, * clint, * cl, ** zv, *** psl;
	extern float ahex;	//needed for hexagon guidance
	extern float offset;	//needed for multiple contour plots per page
	extern float aphex, aphexp, ** hex_norm, ** hex_tang;

	int i, iseg, isp, isl1, isl2, icp, idomain, ndomain, nsp0;
	int ilevel, nlevel = 100;
	int plotVAF = 0;		//only plot VAF if needed
	float xmin, ymin, xs1, ys1, xs2, ys2, xzmin, xzmax, red, green, blue, xz;
	float GF, GFx, GFy, GFz;
	float xmax, ymax, scalefac, ** cos, zcoord, zbottom, ztop, zmin, zmax;
	char fname[80];
	FILE* ofp, * ofp1;

	cos = matrix(1, 3, 1, 3);
	nsp0 = nsp;

	if (hexperiodic) ndomain = 7;
	else ndomain = 1;

	//create file name - need 3-digit frame number
	sprintf(fname, "Current\\ContourO2GF%03i.ps", imain);
	ofp = fopen(fname, "w");
	fprintf(ofp, "%%!PS-Adobe-2.0\n");
	if (nsp == 3 && plotVAF) {
		sprintf(fname, "Current\\ContourVAF%03i.ps", imain);
		ofp1 = fopen(fname, "w");
		fprintf(ofp1, "%%!PS-Adobe-2.0\n");
	}

	for(icp=1; icp<=ncontourplanes; icp++){
		xmax = sqrt(SQR(xsl1[1][icp]-xsl0[1][icp]) + SQR(xsl1[2][icp]-xsl0[2][icp]) + SQR(xsl1[3][icp]-xsl0[3][icp]));
		ymax = sqrt(SQR(xsl2[1][icp]-xsl0[1][icp]) + SQR(xsl2[2][icp]-xsl0[2][icp]) + SQR(xsl2[3][icp]-xsl0[3][icp]));
		if(icp == 1) scalefac = 350./FMAX(xmax,ymax);
		offset = (ncontourplanes - icp)*250.; //375 multiplier, changed to make both frames fit in Vdub, August 2017 JPA

		//Calculate P on a planar slice through the region, specified by three corners and number of points along each edge
		for(isl1=1; isl1<=nsl1[icp]; isl1++) for(isl2=1; isl2<=nsl2[icp]; isl2++){
			for(i=1; i<=3; i++) x[i] = xsl0[i][icp] + (isl1-1)*(xsl1[i][icp]-xsl0[i][icp])/(nsl1[icp]-1) + (isl2-1)*(xsl2[i][icp]-xsl0[i][icp])/(nsl2[icp]-1);
			nsp = 1;	//don't use eval for growth factor
			p = eval(slsegdiv[icp], req, x);	//if hexperiodic, x is remapped into base domain here
			nsp = nsp0;
			psl[isl1][isl2][1] = p[1];
			evalGF(2, req, x, lambdaGF2, cfacGF2, ssGF2, bbGF2, ccGF2, &GF, &GFx, &GFy, &GFz);
			psl[isl1][isl2][2] = GF;
			if (nsp == 3 && plotVAF) {
				evalVAF(3, req, x, lambdaGF3, cfacGF3, ssGF3, bbGF3, ccGF3, &GF, &GFx, &GFy, &GFz);
				psl[isl1][isl2][3] = GF;
			}
		}
		xmin = 0.;
		ymin = 0.;
		for(i=1; i<=3; i++){	//set up matrix of direction cosines
			cos[1][i] = (xsl1[i][icp]-xsl0[i][icp])/xmax;
			cos[2][i] = (xsl2[i][icp]-xsl0[i][icp])/ymax;
		}
		cos[3][1] = cos[1][2]*cos[2][3] - cos[1][3]*cos[2][2];
		cos[3][2] = cos[1][3]*cos[2][1] - cos[1][1]*cos[2][3];
		cos[3][3] = cos[1][1]*cos[2][2] - cos[1][2]*cos[2][1];

	//Determine range of positions along viewing axis
		zmin = 1.e6;
		zmax = -1.e6;
		for(iseg=1; iseg<=nseg; iseg++) if(segtyp[iseg] < 10){	//updated Feb. 2017
			for(idomain=0; idomain<ndomain; idomain++){
				zcoord = 0.;
				for(i=1; i<=3; i++){
					x[i] = cnode[i][ista[iseg]];
					y[i] = cnode[i][iend[iseg]];
					if(hexperiodic == 1 && boundseg[iseg] > 0) y[i] += 2.*aphexp*hex_norm[boundseg[iseg]][i];	//segments that cross boundary
					if(hexperiodic == 1 && idomain > 0){
						x[i] += 2.*aphexp*hex_norm[idomain][i];	//replicate network in each domain
						y[i] += 2.*aphexp*hex_norm[idomain][i];
					}
					zcoord += (x[i] + y[i])/2.*cos[3][i];
				}
				zmin = FMIN(zmin,zcoord-1.);
				zmax = FMAX(zmax,zcoord+1.);
			}
		}

		for(isp=1; isp<=nsp; isp++){
			for(isl1=1; isl1<=nsl1[icp]; isl1++) for(isl2=1; isl2<=nsl2[icp]; isl2++) zv[isl1][isl2] = psl[isl1][isl2][isp];
			for(i=1; i<=nl[isp]; i++) cl[i] = clmin[isp] + (i-1)*clint[isp];
			if(isp == 1){
				cl[1] = 1.;		//override contour levels to give contours at p = 1, 2 and 5
				cl[2] = 2.;
				cl[3] = 5.;
				contr_shade(ofp,nsl1[icp],nsl2[icp],scalefac,nl[isp],xmin,xmax,ymin,ymax,cl,zv,1,1,0,0);
			}
			if(isp == 2){
				cl[1] = thresh;		//contour at threshold for sprout formation, with diagonal hatching
				contr_shade(ofp,nsl1[icp],nsl2[icp],scalefac,1,xmin,xmax,ymin,ymax,cl,zv,0,0,1,1);
			}
			if (isp == 3 && plotVAF) {
				contr_shade(ofp1, nsl1[icp], nsl2[icp], scalefac, nl[isp], xmin, xmax, ymin, ymax, cl, zv, 1, 1, 0, 0);
			}
		}

	//Plot projection of network in contour plane, in order according to rotated z-coordinate
	//plot vessel colors according to pvseg[1] and segtyp
		fprintf(ofp, "/sl {setlinewidth} def\n");
		fprintf(ofp, "/sc {setrgbcolor} def\n");
		fprintf(ofp, "1 setlinecap\n");
		if (nsp == 3 && plotVAF) {
			fprintf(ofp1, "/sl {setlinewidth} def\n");
			fprintf(ofp1, "/sc {setrgbcolor} def\n");
			fprintf(ofp1, "1 setlinecap\n");
		}
		xzmin = clmin[1];
		xzmax = clmin[1] + (nl[1]-1)*clint[1];
		for(ilevel=1; ilevel<=nlevel; ilevel++){
			zbottom = zmin + (ilevel-1)*(zmax - zmin)/nlevel;
			ztop = zmin + ilevel*(zmax - zmin)/nlevel;
			for(iseg=1; iseg<=nseg; iseg++) if(segtyp[iseg] < 10){
				for(idomain=0; idomain<ndomain; idomain++){
					zcoord = 0.;
					for(i=1; i<=3; i++){
						x[i] = cnode[i][ista[iseg]];
						y[i] = cnode[i][iend[iseg]];
						if(hexperiodic == 1 && boundseg[iseg] > 0) y[i] += 2.*aphexp*hex_norm[boundseg[iseg]][i];	//segments that cross boundary
						if(hexperiodic == 1 && idomain > 0){
							x[i] += 2.*aphexp*hex_norm[idomain][i];	//replicate network in each domain
							y[i] += 2.*aphexp*hex_norm[idomain][i];
						}
						zcoord += (x[i] + y[i])/2.*cos[3][i];
					}
					if(zcoord >= zbottom && zcoord < ztop){		
						if(xzmin != xzmax) xz = (pvseg[iseg][1] - xzmin)/(xzmax - xzmin);
						else xz = 0.75;
						if(segtyp[iseg] >= 3 && segtyp[iseg] <= 5){	//Set up colors using Matlab 'jet' scheme based on z value
							red  = FMIN(FMAX(1.5-4.*fabs(xz-0.75), 0.), 1.);
							green= FMIN(FMAX(1.5-4.*fabs(xz-0.5), 0.), 1.);
							blue = FMIN(FMAX(1.5-4.*fabs(xz-0.25), 0.), 1.);
						}
						else if(segtyp[iseg] == 1){		//purple, sprout
							red  = 1.;
							green= 0.;
							blue = 1.;
						}
						else if(segtyp[iseg] == 0){		//redder purple, tip of sprout
							red  = 0.5;
							green= 0.;
							blue = 1.;
						}
						else{							//black, error
							red  = 0.;
							green= 0.;
							blue = 0.;
						}
						xs1 = 0.;
						ys1 = 0.;
						xs2 = 0.;
						ys2 = 0.;
						for(i=1; i<=3; i++){
							xs1 += (x[i]-xsl0[i][icp])*(xsl1[i][icp]-xsl0[i][icp])/xmax;
							ys1 += (x[i]-xsl0[i][icp])*(xsl2[i][icp]-xsl0[i][icp])/ymax;
							xs2 += (y[i]-xsl0[i][icp])*(xsl1[i][icp]-xsl0[i][icp])/xmax;
							ys2 += (y[i]-xsl0[i][icp])*(xsl2[i][icp]-xsl0[i][icp])/ymax;
						}
						//Plot vessels slightly larger in black to outline - moved here Feb. 2017
						fprintf(ofp, "%g sl\n", scalefac* diam[iseg] + 2.);
						fprintf(ofp, "0 0 0 sc\n");
						fprintf(ofp, "%g mx %g my m ", xs1, ys1);
						fprintf(ofp, "%g mx %g my l s \n", xs2, ys2);
						if (nsp == 3 && plotVAF) {
							fprintf(ofp1, "%g sl\n", scalefac * diam[iseg] + 2.);
							fprintf(ofp1, "0 0 0 sc\n");
							fprintf(ofp1, "%g mx %g my m ", xs1, ys1);
							fprintf(ofp1, "%g mx %g my l s \n", xs2, ys2);
						}
						//Plot vessels in color
						fprintf(ofp, "%g sl\n", scalefac* diam[iseg]);
						fprintf(ofp, "%f %f %f sc\n", red, green, blue);
						fprintf(ofp, "%g mx %g my m ", xs1, ys1);
						fprintf(ofp, "%g mx %g my l s \n", xs2, ys2);
						if (nsp == 3 && plotVAF) {
							fprintf(ofp1, "%g sl\n", scalefac * diam[iseg]);
							fprintf(ofp1, "%f %f %f sc\n", red, green, blue);
							fprintf(ofp1, "%g mx %g my m ", xs1, ys1);
							fprintf(ofp1, "%g mx %g my l s \n", xs2, ys2);
						}
					}
				}
			}
		}
		//plot hexagons in gray, requires that plot no. 1 is in x-y plane
		if(hexperiodic && icp == 1){
			fprintf(ofp,"%g sl\n",1.);
			fprintf(ofp,"0.5 0.5 0.5 sc\n");
			for(idomain=0; idomain<ndomain; idomain++){
				for(i=1; i<=6; i++){
					xs1 = aphexp*hex_norm[i][1] - 0.5*aphex*hex_tang[i][1] + midpt[1];
					xs2 = aphexp*hex_norm[i][1] + 0.5*aphex*hex_tang[i][1] + midpt[1];
					ys1 = aphexp*hex_norm[i][2] - 0.5*aphex*hex_tang[i][2] + midpt[2];
					ys2 = aphexp*hex_norm[i][2] + 0.5*aphex*hex_tang[i][2] + midpt[2];
					if(idomain > 0){
						xs1 += 2.*aphexp*hex_norm[idomain][1];
						xs2 += 2.*aphexp*hex_norm[idomain][1];
						ys1 += 2.*aphexp*hex_norm[idomain][2];
						ys2 += 2.*aphexp*hex_norm[idomain][2];
					}
					xs1 = (xs1-xsl0[1][icp])*(xsl1[1][icp]-xsl0[1][icp])/xmax;
					ys1 = (ys1-xsl0[2][icp])*(xsl2[2][icp]-xsl0[2][icp])/ymax;
					xs2 = (xs2-xsl0[1][icp])*(xsl1[1][icp]-xsl0[1][icp])/xmax;
					ys2 = (ys2-xsl0[2][icp])*(xsl2[2][icp]-xsl0[2][icp])/ymax;
					fprintf(ofp, "%g mx %g my m ", xs1,ys1);
					fprintf(ofp, "%g mx %g my l s \n", xs2,ys2);
				}
			}
		}
	}
	fprintf(ofp, "showpage\n");
	fclose(ofp);
	if (nsp == 3 && plotVAF) {
		fprintf(ofp1, "showpage\n");
		fclose(ofp1);
	}
	free_matrix(cos,1,3,1,3);
}