/******************************************************************************
sprgrowth.cpp - identification of new sprouts
New version for AngioAdapt17GPU
Combines old functions of prespr_estimat, sprfind, sprgrowth and sprupdate
For each existing segment, the following steps are performed:
1) determine whether it sprouts:
-find growth factor concentration at midpoint of the segment.
-if growth factor concentration > threshold, assigns a probability of sprouting = smax*f(thresh)*seg_length*dt
-find a random number between 0 and 1.  If that is < probability, a sprout can form.
2) determine where it sprouts from:
-find a point on the segment with uniform probability
-test if point is too close to an existing node - then move to that node
-nullify the sprout if it tries to sprout from an input or output
-nullify the sprout if it tries to sprout from a node that has three or more branches
-create a new start node and subdivide parent segment if needed
3) choose direction of sprout:
-in 3d, sprouts form perpendicularly with a random angle about the parent segment
-in 2d, sprouts form perpendicularly to the left or right according to GF gradient
J. Alberding, August 2006
Revised TWS February 2018
********************************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include "nrutil.h"

void evalGF(int isp, float req, float* x, float lambdaGF, float cfacGF, float ssGF, float bbGF, float ccGF,
	float* GF, float* GFx, float* GFy, float* GFz);
int close_pt_finder(float *new_pt, float *int_pt, int inod_from, int iseg_from);
int newnode(float *int_pt, int iseg_to);
float length(float *x);
float length0(float *x);
float distance(float *x, float *y);

void sprgrowth()
{
	extern int angio3d,dogreensflag;
	extern int nseg,nsprout,nsprmax,sprmaxflag,nnodbc,actnum;
	extern int *activenode,*segtyp,*ista,*iend,*nodtyp,*bcnod;
	extern int **nodseg;
	extern float thresh,threshn,thresh50,req,smax;
	extern float lambdaGF2, cfacGF2, ssGF2, bbGF2, ccGF2;
	extern float lambdaGF3, cfacGF3, ssGF3, bbGF3, ccGF3;
	extern float timestep1,tolerance_2,pi1,Vmax;
	extern float *x,*y,*P1,*P2,*lseg,**cnode,**activegr;

	int iseg,i,inod1,inod2,inod_from,inodbc,nseg0,iact;
	float GF,GFx,GFy,GFz;
	float *V1,*V2,*gdir,*pdir;
	float rnum,prob,lfac,thresh_func;
	float magv1,magv2,orthotest,magpdir,angle2,dist1,dist2,factor;

	V1 = vector(1,3);
	V2 = vector(1,3);	
	pdir = vector(1,3);
	gdir = vector(1,3);

	nsprout = 0;
	iact = actnum;
	nseg0 = nseg;	//do not consider segments created during sprout formation
	for(iseg=1; iseg<=nseg0; iseg++) if(segtyp[iseg] <= 1 || segtyp[iseg] == 5){
		inod1 = ista[iseg];
		inod2 = iend[iseg];
		for(i=1; i<=3; i++){
			P2[i] = cnode[i][inod2] - cnode[i][inod1];
			P1[i] = cnode[i][inod1];
		}
		length(P2);	//remaps to base domain
		for(i=1; i<=3; i++)	P2[i] += cnode[i][inod1];	//may be outside domain
		//********************** determine whether it sprouts *****************************
		for(i=1; i<=3; i++) x[i] = (P1[i] + P2[i]) * 0.5;
		if(dogreensflag == 0) GF = 1.25 * thresh;
		else evalGF(2, req, x, lambdaGF2, cfacGF2, ssGF2, bbGF2, ccGF2, &GF, &GFx, &GFy, &GFz);
		if(GF <= thresh) goto nosprout;
		thresh_func = 1. - 1. / (1. + pow((GF - thresh) / thresh50, threshn));
		prob = smax*lseg[iseg]*thresh_func*timestep1;		//Calculate probability of sprouting
		rnum = rand()*1./RAND_MAX;
		if(rnum >= prob) goto nosprout;
		//********************** determine where it sprouts from **************************
		lfac = rand()*1./RAND_MAX;	//find a point on each segment with uniform probability
		for(i=1; i<=3; i++) y[i] = (1. - lfac)*P1[i] + lfac*P2[i];
		dist1 = lseg[iseg]*lfac;
		dist2 = lseg[iseg]*(1. - lfac);
		inod_from = -1;	//start by assuming need to create a new node
		if(dist1 < tolerance_2 || dist2 < tolerance_2){	//modified June 2016 to avoid losing many segments
			if(dist1 <= dist2){
				if(nodtyp[inod1] <= 2){	//move to this node if it is type 2
					inod_from = inod1;
					for(i=1; i<=3; i++) y[i] = P1[i];
				}
				else if(lseg[iseg] > 2.*tolerance_2){	//move a distance tolerance_2 away from the node				
					factor = tolerance_2/lseg[iseg];
					for(i=1; i<=3; i++) y[i] = P1[i]*(1. - factor) + P2[i]*factor;
				}
				else goto nosprout;
			}
			else{
				if(nodtyp[inod2] <= 2){	//move to this node if it is type 2
					inod_from = inod2;
					for(i=1; i<=3; i++) y[i] = P2[i];
				}
				else if(lseg[iseg] > 2.*tolerance_2){	//move tolerance_2 away from the node
					factor = tolerance_2/lseg[iseg];
					for(i=1; i<=3; i++) y[i] = P2[i]*(1. - factor) + P1[i]*factor;
				}
				else goto nosprout;
			}
			for(inodbc=1; inodbc<=nnodbc; inodbc++) if(inod_from == bcnod[inodbc]) goto nosprout;	//no sprout from an I/O node
		}	
		if(close_pt_finder(y, x, 0, iseg)) goto nosprout;		//check whether sprout point is too close to an existing vessel
		if(inod_from < 0) inod_from = newnode(y, iseg);  //if needed, create a new start node and subdivide parent segment
		//********************** choose direction of sprout *****************************
		for(i=1; i<=3; i++)	pdir[i] = P1[i] - P2[i];
		magpdir = length(pdir);	//if hexperiodic, pdir may be modified
		for(i=1; i<=3; i++) pdir[i] = pdir[i]/magpdir;	//make a unit vector
		if(angio3d){			//3d sprouts. Find random angle between 0 and 2*pi to rotate the sprout direction
			if(pdir[1] != 0.){	//first perpendicular vector
				V1[1] = -pdir[3];
				V1[2] = 0.;
				V1[3] = pdir[1];
			}
			else{
				V1[1] = 0.;
				V1[2] = pdir[3];
				V1[3] = -pdir[2];
			}
			V2[1] = pdir[2]*V1[3] - pdir[3]*V1[2];	//second perpendicular vector
			V2[2] = pdir[3]*V1[1] - pdir[1]*V1[3];
			V2[3] = pdir[1]*V1[2] - pdir[2]*V1[1];
			magv1 = length0(V1);
			magv2 = length0(V2);
			rnum = rand()*1./RAND_MAX;
			angle2 = rnum*2.*pi1;
			for(i=1; i<=3; i++) gdir[i] = V1[i]/magv1*cos(angle2) + V2[i]/magv2*sin(angle2);
		}
		else{			//choose direction up gradient in growth factor
			gdir[1] = -pdir[2];
			gdir[2] = pdir[1];
			gdir[3] = 0.0;
			evalGF(2, req,y, lambdaGF2, cfacGF2, ssGF2, bbGF2, ccGF2, &GF,&GFx,&GFy,&GFz);
			if(GFx*gdir[1] + GFy*gdir[2] < 0.) for(i=1; i<=3; i++) gdir[i] = -gdir[i];
		}
		//test orthogonality of sprout
		orthotest = fabs(pdir[1]*gdir[1] + pdir[2]*gdir[2] + pdir[3]*gdir[3]);
		if(orthotest > 0.01) printf("*** Non-orthogonal sprout created at seg %i\n", iseg);
		//********************** determine projected endpoint of sprout *****************************
		if(nsprout >= nsprmax){
			sprmaxflag = 1;
			printf("*** Error: too many sprouts, nsprout = %i, max = %i\n",nsprout,nsprmax);
		}
		else{
			nsprout++;
			iact++;
			for(i=1; i<=3; i++){			//define endpoint of new segment
				activegr[iact][i] = y[i];
				activegr[iact][i+3] = y[i] + Vmax*timestep1*gdir[i];
			}
			activenode[iact] = inod_from;
		}
nosprout:;
	}
	actnum = iact;
	free_vector(pdir,1,3);
	free_vector(gdir,1,3);
	free_vector(V1,1,3);
	free_vector(V2,1,3);
}
