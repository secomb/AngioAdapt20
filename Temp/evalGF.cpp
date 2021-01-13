/*****************************************************
eval - Evaluate GF from source strengths, using Green's
function for diffusion with linear uptake.
TWS October 2010.
JPA February 2015
Updated to include gradients in GF, TWS May 2016
Updated to include 3D and 2D versions, TWS September 2016
Evaluation of constants moved to evalGFinit, March 2018
******************************************************/
#include <math.h>
#include "nrutil.h"
#include <stdio.h>

float bessi0(float x);
float bessk0(float x);
float bessi1(float x);
float bessk1(float x);
void tissrate(int nsp, float *p, float *mtiss, float *mptiss);
float length(float *x);

void GFGF2D(float req, float *y, float lambdaGF, float cfacGF, float ssGF, float bbGF, float ccGF,
	float *GFGF, float *GFGFx, float *GFGFy, float *GFGFz){
	/**********************************
	2D Greens function for growth factor
	with linear uptake in tissue
	***********************************/
	extern float pi1;
	float dist,gff,gfp;

	//lambdaGF = sqrt(KG/DG);
	//ssGF = 1./(pi1*req*req*w2d*KG); 
	//cfacGF = bessi0(lambdaGF*req)*bessk1(lambdaGF*req) + bessi1(lambdaGF*req)*bessk0(lambdaGF*req);
	//bbGF = ssGF*bessi1(lambdaGF*req)/cfacGF;
	//ccGF = ssGF*bessk1(lambdaGF*req)/cfacGF;
	dist = sqrt(SQR(y[1]) + SQR(y[2]));
	if(dist < req){
		gff = ssGF - ccGF*bessi0(dist*lambdaGF);
		gfp = -lambdaGF*ccGF*bessi1(dist*lambdaGF);
	}
	else{
		gff = bbGF*bessk0(dist*lambdaGF);
		gfp = -lambdaGF*bbGF*bessk1(dist*lambdaGF);
	}
	*GFGF = gff;
	if(dist == 0.){
		*GFGFx = 0.;
		*GFGFy = 0.;
	}
	else{
		*GFGFx = gfp*y[1]/dist;
		*GFGFy = gfp*y[2]/dist;
	}
	*GFGFz = 0.;
}

void GFGF3D(float req, float *y, float lambdaGF, float cfacGF, float ssGF, float bbGF, float ccGF, float *GFGF, float *GFGFx, float *GFGFy, float *GFGFz){
	/**********************************
	3D Greens function for growth factor
	with linear uptake in tissue
	***********************************/
	extern float pi1;
	float dist,gff,gfp;

	//lambdaGF = sqrt(KG/DG);
	//ssGF = 3./(4.*pi1*req*req*req*KG); 
	//bbGF = ssGF*req*(cosh(lambdaGF*req) - sinh(lambdaGF*req)/(lambdaGF*req));
	//ccGF = ssGF*req*(1.+1./(lambdaGF*req))*exp(-lambdaGF*req);
	dist = sqrt(SQR(y[1]) + SQR(y[2]) + SQR(y[3]));
	if(dist == 0.){
		gff = ssGF - ccGF*lambdaGF;
		gfp = 0.;
	}
	else if(dist < req){
		gff = ssGF - ccGF*sinh(dist*lambdaGF)/dist;
		gfp = -ccGF*(lambdaGF*cosh(dist*lambdaGF)/dist - sinh(dist*lambdaGF)/SQR(dist));
	}
	else{
		gff = bbGF*exp(-dist*lambdaGF)/dist;
		gfp = -bbGF*exp(-dist*lambdaGF)*(lambdaGF/dist + 1./SQR(dist));
	}
	*GFGF = gff;
	if(dist == 0.){
		*GFGFx = 0.;
		*GFGFy = 0.;
		*GFGFz = 0.;
	}
	else{
		*GFGFx = gfp*y[1]/dist;
		*GFGFy = gfp*y[2]/dist;
		*GFGFz = gfp*y[3]/dist;
	}
}
	
void evalGF(int isp, float req, float *x, float lambdaGF, float cfacGF, float ssGF, float bbGF, float ccGF, 
	float *GF, float *GFx, float *GFy, float *GFz)
{
	extern int nnt,is2d,nsp,**tisspoints,hexperiodic;
	extern int retinaflag,mxx,myy,mzz,***nbou;	//needed for retina version
	extern float pi1,w2d,alz,*axt,*ayt,*azt,*y,*P3,*midpt;
	extern float **qt,**tissparam,*diff;
	extern float IR,OR,*ptpt,**pt,*mtiss,*mptiss,vol;	//needed for retina version
	extern float aphex,aphexp,**hex_norm;	//needed for hexperiodic

	int i,itp,ndomain,idomain,idomain1,ndomain1;
	int ix,iy,iz,k;	//needed for retina version
	float GFGF,GFGFx,GFGFy,GFGFz;
	float dist1,qt1;	//needed for retina version

	ndomain1 = 3;	//3 to include top and bottom reflections, 1 to exclude

	if(hexperiodic){
		for(i=1; i<=3; i++) y[i] = x[i] - midpt[i];
		length(y);		//remap x into base domain - gives periodic solute fields
		for(i=1; i<=3; i++) x[i] = y[i] + midpt[i];
	}

	*GF = 0.;
	*GFx = 0.;
	*GFy = 0.;
	*GFz = 0.;	//fixed Feb. 2017
	if(retinaflag){		//This version gives a source of GF at all points in annulus, even if not tissue points
		for(ix=1; ix<=mxx; ix++) for(iy=1; iy<=myy; iy++) for(iz=1; iz<=mzz; iz++){
			dist1 = SQR(axt[ix] - midpt[1]) + SQR(ayt[iy] - midpt[2]);	//distance from midpoint of domain
			if(dist1 >= SQR(IR) && dist1 <= SQR(OR)){
				itp = nbou[ix][iy][iz];
				if(itp == 0) ptpt[1] = 0.;	//maximal GF production rate for points not in tissue domain
				else ptpt[1] = pt[itp][1];			
				tissrate(nsp,ptpt,mtiss,mptiss);
				qt1 = mtiss[isp]*vol;	//GF production rate
				y[1] = x[1] - axt[ix];
				y[2] = x[2] - ayt[iy];
				y[3] = x[3] - azt[iz];
				if(is2d) GFGF2D(req, y, lambdaGF, cfacGF, ssGF, bbGF, ccGF, &GFGF, &GFGFx, &GFGFy, &GFGFz);
				else GFGF3D(req, y, lambdaGF, cfacGF, ssGF, bbGF, ccGF, &GFGF, &GFGFx, &GFGFy, &GFGFz);
				*GF += qt1*GFGF;
				*GFx += qt1*GFGFx;
				*GFy += qt1*GFGFy;
				*GFz += qt1*GFGFz;
			}
		}
	}
	else{
		for(itp=1; itp<=nnt; itp++){
			//ptpt[1] = pt[itp][1];			
			//tissrate(nsp,ptpt,mtiss,mptiss);
			//qt1 = mtiss[isp]*vol;	//GF production rate
			qt1 = qt[itp][isp];	//see calculation of qt in main.cpp
			P3[1] = x[1] - axt[tisspoints[1][itp]];
			P3[2] = x[2] - ayt[tisspoints[2][itp]];
			P3[3] = x[3] - azt[tisspoints[3][itp]];
			if(hexperiodic)	ndomain = 7;
			else ndomain = 1;
			for(idomain1=1; idomain1<=ndomain1; idomain1++){
				for(idomain=0; idomain<ndomain; idomain++){
					for(k=1; k<=3; k++){
						y[k] = P3[k];
						if(idomain > 0 && idomain <=6) y[k] += 2.*aphexp*hex_norm[idomain][k];
					}
					if(idomain1 == 2) y[3] = x[3] + azt[tisspoints[3][itp]];	//reflection from bottom of domain
					if(idomain1 == 3) y[3] = x[3] - 2*alz + azt[tisspoints[3][itp]]; //reflection from top of domain
					if(is2d) GFGF2D(req, y, lambdaGF, cfacGF, ssGF, bbGF, ccGF, &GFGF, &GFGFx, &GFGFy, &GFGFz);
					else GFGF3D(req, y, lambdaGF, cfacGF, ssGF, bbGF, ccGF, &GFGF, &GFGFx, &GFGFy, &GFGFz);
					*GF += qt1*GFGF;
					*GFx += qt1*GFGFx;
					*GFy += qt1*GFGFy;
					*GFz += qt1*GFGFz;
				}
			}
		}
	}
}

void evalVAF(int isp, float req, float* x, float lambdaGF, float cfacGF, float ssGF, float bbGF, float ccGF,
	float* GF, float* GFx, float* GFy, float* GFz)
{
	extern int nnv, is2d, hexperiodic, *mainseg;
	extern float pi1, w2d, alz, **ax, * y, * P3, * midpt, *ds, *rseg;
	extern float** tissparam;
	extern float aphex, aphexp, ** hex_norm;	//needed for hexperiodic

	int i, j, iseg, ndomain, idomain, idomain1, ndomain1;
	float GFGF, GFGFx, GFGFy, GFGFz;
	float dist1, qv1;	//needed for retina version

	ndomain1 = 3;	//3 to include top and bottom reflections, 1 to exclude

	if (hexperiodic) {
		for (i = 1; i <= 3; i++) y[i] = x[i] - midpt[i];
		length(y);		//remap x into base domain - gives periodic solute fields
		for (i = 1; i <= 3; i++) x[i] = y[i] + midpt[i];
	}
	*GF = 0.;
	*GFx = 0.;
	*GFy = 0.;
	*GFz = 0.;	
	for (i = 1; i <= nnv; i++) {
		iseg = mainseg[i];
		qv1 = tissparam[1][isp] * 2. * pi1 * ds[iseg] * rseg[iseg];  //source strength * area of subsegment
		for (j = 1; j <= 3; j++) P3[j] = x[j] - ax[j][i];
		if (hexperiodic) ndomain = 7;
		else ndomain = 1;
		for (idomain1 = 1; idomain1 <= ndomain1; idomain1++) {
			for (idomain = 0; idomain < ndomain; idomain++) {
				for (j = 1; j <= 3; j++) {
					y[j] = P3[j];
					if (idomain > 0 && idomain <= 6) y[j] += 2. * aphexp * hex_norm[idomain][j];
				}
				if (idomain1 == 2) y[3] = x[3] + ax[3][i];	//reflection from bottom of domain
				if (idomain1 == 3) y[3] = x[3] - 2 * alz + ax[3][i]; //reflection from top of domain
				if (is2d) GFGF2D(req, y, lambdaGF, cfacGF, ssGF, bbGF, ccGF, &GFGF, &GFGFx, &GFGFy, &GFGFz);
				else GFGF3D(req, y, lambdaGF, cfacGF, ssGF, bbGF, ccGF, &GFGF, &GFGFx, &GFGFy, &GFGFz);
				*GF += qv1 * GFGF;
				*GFx += qv1 * GFGFx;
				*GFy += qv1 * GFGFy;
				*GFz += qv1 * GFGFz;
			}
		}
	}
}
