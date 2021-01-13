/****************************************************************
tension.cpp
This program assigns a set tension to all segments. Physically,
real vessels maintain tension so as to distribute their trinode
angles in a bell shaped curve about 2pi/3, and to avoid sharp angles
at vertices (< pi/4, opposite anterior angles).
Specifically, the program:
1. Scans for all nodes with two or more segments only.
2. Assigns a tension to each segment with direction along the segment.
3. If resultant magnitude > lambda (tbd), moves the node to a point
1 micron in the resultant direction.
4. Tests the new segments formed by the new node and the adjacent
nodes to see if they are < the minimum length.
5. If any segment is shorter than the minimum length, the line between
the adjacent node and the modified node becomes a new segment and
the original segments are deactivated (type 10)
J. Alberding, June 13, 2007
Revised TWS 2010
Revised TWS April 2020 to include segment length in criterion for moving nodes
********************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include "nrutil.h"

void sht_seg_fixer(int inod_to, int seg_no, int inod);
void prep_for_adapt(int inod_to);
float length(float *x);
float length0(float *x);
float distance(float *P1, float *P2);

void tension()
{
	extern int nseg,nnod,imain;
	extern int *segname,*segtyp,*ista,*iend,*nodtyp;
	extern int **nodseg,**nodnod;
	extern float lambdatension,tensionrate,timestep1,tolerance_2;
	extern float *x,*P1,*P2,*diam, * lseg;
	extern float **cnode;

	int seg_no,seg_num,i,ii,j,inod,iseg,inod_to,inod1,tensflag,vessmovecount;
	int notinod_to;
	float magpdir,magtenres,diamsum,fac,diamlengthsum;
	float *new_new_pt,*pdir,*tenres;
	float testlength, znod;
	float bottomzonedecay = 50.;	//used to keep vessels from being pulled out of the bottom zone
	float bottomzoneforce = 3.;
	FILE *ofp;

	new_new_pt = vector(1,3);
	pdir = vector(1,3);
	tenres = vector(1,3);
	vessmovecount = 0;

	for(inod=1; inod<=nnod; inod++) if(nodtyp[inod] >= 2){
		for(j=1; j<=3; j++) tenres[j] = 0.;
		tensflag = 1;
		for(j=1; j<=nodtyp[inod]; j++) if(segtyp[nodseg[j][inod]] == 4) tensflag = 0;//skip nodes connected to type 4
		if(tensflag == 1){
			diamsum = 0.;
			diamlengthsum = 0.;
			for(i=1; i<=nodtyp[inod]; i++){
				iseg = nodseg[i][inod];
				if (segtyp[iseg] != 10) {
					diamsum += diam[iseg];
					diamlengthsum += diam[iseg] * lseg[iseg];
				for(j=1; j<=3; j++) pdir[j] = cnode[j][nodnod[i][inod]] - cnode[j][inod];
				magpdir = length(pdir);	//if hexperiodic, pdir may be changed
					if (magpdir != 0.) for (j = 1; j <= 3; j++) tenres[j] += diam[iseg] * pdir[j] / magpdir;
				//drag down nodes near the bottom, to avoid shrinkage of network. Feb. 2020
				znod = cnode[3][inod];
				if (bottomzoneforce > 0. && znod > 0.)  tenres[3] -= bottomzoneforce * exp(-znod / bottomzonedecay);
			}
			}
			for(j=1; j<=3; j++) tenres[j] = tenres[j] * diamsum / diamlengthsum;
			magtenres = length0(tenres);
			fac = 0.;
			if(magtenres > lambdatension){		//Continuous dependence on magtenres.  TWS 2010
				fac = timestep1*tensionrate*(1. - lambdatension/magtenres);
				vessmovecount++;
			}
			for(j=1; j<=3; j++)	new_new_pt[j] = cnode[j][inod] + fac*tenres[j];
			//update node_table again.  Test each segment attached to this node to see if the new_new_point is within
			//tolerance 2 of the other nodes for these segments
			for(i=1; i<=nodtyp[inod]; i++) if(segtyp[nodseg[i][inod]] != 10){
				seg_no = nodseg[i][inod];
				inod1 = nodnod[i][inod];
				for(j=1; j<=3; j++) P1[j] = cnode[j][inod1];
				if(distance(new_new_pt,P1) < tolerance_2){
					inod_to = inod1;
					sht_seg_fixer(inod_to, seg_no, inod);
					for(ii=1; ii<=nodtyp[inod_to]; ii++){
						seg_num = nodseg[ii][inod_to];
						if(segtyp[seg_num] != 10){
							for(i=1; i<=3; i++) x[i] = cnode[i][ista[seg_num]] - cnode[i][iend[seg_num]];
							testlength = length(x);
							if(testlength <= tolerance_2){
								if(nodtyp[ista[seg_num]] < nodtyp[iend[seg_num]]){
									notinod_to = ista[seg_num];
									inod_to = iend[seg_num];
								}
								else{
									notinod_to = iend[seg_num];
									inod_to = ista[seg_num];
								}
								sht_seg_fixer(inod_to, seg_num, notinod_to);
							}
						}
					}
					prep_for_adapt(inod_to);
					goto move_done;
				}
			}
			for(j=1; j<=3; j++) cnode[j][inod] = new_new_pt[j];
move_done:;
		}
	}
//do preadapt again
	for(inod=1; inod<=nnod; inod++) prep_for_adapt(inod);
	if(imain == 1) ofp = fopen("vmovecount.out", "w");
	else ofp = fopen("vmovecount.out", "a");
	fprintf(ofp,"At imain = %i, vessmovecount = %i\n", imain, vessmovecount);
	fclose(ofp);

	free_vector(new_new_pt,1,3);
	free_vector(pdir,1,3);
	free_vector(tenres,1,3);
}
