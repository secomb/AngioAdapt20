/***********************************************************************************************
actgrowth.cpp - identification of growing sprouts
Program to grow a new vessel from an active segment.
The factors that determine the direction are:
-The previous vessel direction, calculated in the main program.  A new active
segment tends to continue in the direction of the parent segment.  A sprout 
forms perpendicular to the parent segment's direction, randomly choosing 
between the angles -90 and 90 degrees from the parent segment's longitudinal axis in 2d
or randomly perpendicular to the axis in 3d.
-A random angular component, given by a gaussian distribution about the direction computed from the
previous/perpendicular, ICD, and cgfactor_grad component directions.
-The direction signalled from the integrated conical detector (ICeD).
-a directional component following the gfactor gradient, from low to high concentration - currently disabled.

The factors determining the magnitude of the new growth are:
-At the inception of this program, new growth is set to a constant length.
Another program then tests for intercepts while incrementing along this length in
discrete steps.  Segment length may then be truncated if intercepts occur.
-There is an option to calculate a variable new growth length depending on growth factor, 
based upon a formula in Tong and Yuan, Microvascular Research, 61, 14-27, 2001.
J. Alberding, August 28 2006
Updated Oct. 2016 to use angio3d to determine 2D vs 3D angiogenesis
***************************************************************************************************/

#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"
#include <time.h>

float randgauss();
float length(float *x);
float length0(float *x);
float ICeD(float *direct, float *xp, float *signal, float r0, float theta0);
void evalGF(int isp, float req, float* x, float lambdaGF, float cfacGF, float ssGF, float bbGF, float ccGF,
	float* GF, float* GFx, float* GFy, float* GFz);
void evalVAF(int isp, float req, float* x, float lambdaGF, float cfacGF, float ssGF, float bbGF, float ccGF,
	float* GF, float* GFx, float* GFy, float* GFz);
//void hexdomaintest(int iact);

void actgrowth()
{
	extern int dogreensflag,hexperiodic,nsp;	
	extern int actnum,nseg,nnod,mxx,myy,mzz,angio3d;
	extern int *activenode,*segname,*segtyp,*nodname,*ista,*iend,**nodseg,*nodtyp;
	extern float req,tolerance_1,tolerance_2,variance,Vmax,increment, kgGF, kgVAF, kv,ko,r0,theta0,timestep1;
	extern float lambdaGF2, cfacGF2, ssGF2, bbGF2, ccGF2;
	extern float lambdaGF3, cfacGF3, ssGF3, bbGF3, ccGF3;
	extern float *p,*pref,*P1,*P2;
	extern float **cnode,**gradp,**activegr;
	
	int i,iseg,iact;
	float length1,magp1,magp2,angle,angle1,angle2,angle3,angle4;
	float *direct1,*direct,*endpt,*signal;
	float GF,GFx,GFy,GFz;

	direct1 = vector(1,3);
	direct = vector(1,3);
	endpt = vector(1,3);
	signal = vector(1,3);

	//continue growth of an active segment
	iact = 0;
	for(iseg=1; iseg<=nseg; iseg++) if(segtyp[iseg] == 0){
		iact++;
		activenode[iact] = iend[iseg];
	}
	if(iact != actnum) printf("*** Inconsistent number of active segments %i != %i\n",iact, actnum); 
	for(iact=1; iact<=actnum; iact++){
		iseg = nodseg[1][activenode[iact]];		
		for(i=1; i<=3; i++){
			endpt[i] = cnode[i][iend[iseg]];		//endpt - endpoint of active segment
			direct[i] = cnode[i][iend[iseg]] - cnode[i][ista[iseg]];	//direct - unit vector in direction of active segment
		}
		length1 = length(direct);		//note that direct can be modified if hexperiodic is active
		for(i=1; i<=3; i++) direct[i] = direct[i]/length1;					
//************** generate random change of growth angle
		if(angio3d){			//3d version, updated Jan 2011
			angle1 = randgauss()*sqrt(variance);
			angle2 = randgauss()*sqrt(variance);
			angle3 = sqrt(SQR(angle1) + SQR(angle2));
			angle4 = atan2(angle1,angle2);
			//construct a unit vector perpendicular to the segment
			if(direct[1] != 0.0){
				P1[1] = -direct[3];
				P1[2] = 0.;
				P1[3] = direct[1];
			}
			else{
				P1[1] = 0.;
				P1[2] = direct[3];
				P1[3] = -direct[2];
			}
			//construct a second unit vector perpendicular to both
			P2[1] = direct[2]*P1[3] - direct[3]*P1[2];
			P2[2] = direct[3]*P1[1] - direct[1]*P1[3];
			P2[3] = direct[1]*P1[2] - direct[2]*P1[1];
			magp1 = length0(P1);
			magp2 = length0(P2);
			for(i=1; i<=3; i++){				//normalize
				P1[i] = P1[i]/magp1;
				P2[i] = P2[i]/magp2;
			}
			//calculate active segment direction
			for(i=1; i<=3; i++) direct1[i] = direct[i]*cos(angle3)
				+ P1[i]*sin(angle3)*cos(angle4) + P2[i]*sin(angle3)*sin(angle4);
		}
		else{			//2d version
			angle = randgauss()*sqrt(variance);
			direct1[1] = direct[1]*cos(angle) - direct[2]*sin(angle);
			direct1[2] = direct[1]*sin(angle) + direct[2]*cos(angle);
			direct1[3] = direct[3];
		}

		//obtain signal from cone detector, using all segments except type 10
		ICeD(direct,endpt,signal,r0,theta0);

		//GF gradient vector
		if (dogreensflag && nsp >= 2) {
			evalGF(2, req, endpt, lambdaGF2, cfacGF2, ssGF2, bbGF2, ccGF2, &GF, &GFx, &GFy, &GFz);
			gradp[1][2] = GFx;
			gradp[2][2] = GFy;
			gradp[3][2] = GFz;
		}
		
		//VAF gradient vector
		if (dogreensflag && nsp >= 3) {
			evalVAF(3, req, endpt, lambdaGF3, cfacGF3, ssGF3, bbGF3, ccGF3, &GF, &GFx, &GFy, &GFz);
			gradp[1][3] = GFx;
			gradp[2][3] = GFy;
			gradp[3][3] = GFz;
		}

		//combine with previous direction in 3D (even for 2D problem)
		//create unit vector direct1 in growth direction
		for(i=1; i<=3; i++){
			direct1[i] += kv*signal[i];						//contribution from cone detector
			if (dogreensflag && nsp >= 2) direct1[i] += kgGF * gradp[i][2];	//contribution from GF gradient - up gradient
			if (dogreensflag && nsp >= 3) direct1[i] -= kgVAF * gradp[i][3];	//contribution from VAF gradient - down gradient
		}

		length1 = length0(direct1);
		for(i=1; i<=3; i++){			//define endpoint of new segment
			activegr[iact][i] = endpt[i];
			activegr[iact][i+3] = endpt[i] + Vmax*timestep1*direct1[i]/length1;
		}
	}
	free_vector(direct1,1,3);
	free_vector(signal,1,3);
	free_vector(direct,1,3);
	free_vector(endpt,1,3);
}
