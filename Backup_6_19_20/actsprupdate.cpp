/**********************************************************************
actsprupdate.cpp
Checks both sprouts and active segments
Combines old actupdate and sprupdate.

Checks growing segments(type 0):

1. To assure that they will produce a segment length greater than
the minimum length allowed.
2. To see if the segment would cross the boundary of the box (makes a 
type 10 segment if yes).
3. To see if the segment comes within tolerance to intercept and connect
with another segment.

If an intercept occurs, actupdate checks:
	a. Is the intercept at a node point?  If so, actupdate checks if
		the node is an input/output node (yes => doesn't make new
		segment), and checks if node is another active segment's "to"
		node (takes other seg off of active list).  Otherwise,
		forms new segment w/o appending any new nodes.
	b. If does not intercept within tolerance at a node, creates new node at the intercept point
		and two new segments.  updates the old segment.

Otherwise, actupdate creates a new active segment at the end of the previous
active segment.

Numbering scheme for case when vessel intercepts another

                          |
						  * inod1
   start_pt   new_pt      |
   inod_prev  inod_from   | iseg_to
------*----------*........* inod_to = inod
       iseg_from    iseg  |
                          |
						  * inod2
						  |

J. Alberding, Sep. 19, 2006.  TWS January 08. Comment updated J. Alberding May 09
Revised TWS2010
actsprupdate.cpp - Feburary 2018
******************************************************************************/

#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

int intercept_check(float *start_pt, float *new_pt, float *int_pt, int inod_from);
float distance(float *x, float *y);
float length0(float *x);
int newnode(float *int_pt, int iseg_to);
float randgauss();

void actsprupdate()
{
	extern int mxx,myy,mzz,nseg,nnod,nnodbc,actnum,nodnamemax,segnamemax,nodsegm,nsprout,nsprmax;
	extern int nsplus,nseg_new,nnod_new,hexperiodic;
	extern int retinaflag; //for retinal annulus
	extern int *activenode,*segname,*segtyp,*ista,*iend,*nodname,*bcnod,*nodtyp,*sprseg_extra;
	extern int **nodseg,**nodnod,**segnodname;
	extern int ***nbou;
	extern float alx,aly,alz,tolerance_1,tolerance_2, new_diam, min_length,ks1,ranks1,thresh,req;
	extern float IR,OR;	//for retinal annulus
	extern float *diam,*q,*hd,*segvar,*nodvar,*ksseg,*P1,*P2,*midpt;
	extern float **activegr,**cnode;
	extern FILE *ofp2;

	int i,iseg,iact,inod,inod1,inod2,interceptflag,inodbc,ipt,jpt,kpt,newsprout;
	int inod_from,inod_to,iseg_from,inod_prev,iseg_to;
	float seg_length,dist1,dist2;
	float GF,GFx,GFy,GFz;
	float *start_pt,*new_pt,*int_pt;

	start_pt = vector(1,3);
	new_pt = vector(1,3);
	int_pt = vector(1,3);

	inod = nnod;
	iseg = nseg;
	for(iact=1; iact<=actnum; iact++){	//same procedure for both growing sprouts and new sprouts
		iseg_from = nodseg[1][activenode[iact]];
		newsprout = 1;
		if(segtyp[iseg_from] == 0){
			segtyp[iseg_from] = 1;	//if it is from a growing sprout
			newsprout = 0;
		}
		for(i=1; i<=3; i++){
			start_pt[i] = activegr[iact][i];
			new_pt[i] = activegr[iact][i+3];
		}
		seg_length = distance(start_pt,new_pt);
//test new segment length
		if(seg_length < min_length)	goto nonewseg;
//test for growth beyond the region of included tissue points
		if(retinaflag){
			if(new_pt[3] <= 0. || new_pt[3] >= alz) goto nonewseg;
			P1[1] = new_pt[1] - midpt[1];		//for retinal annulus - test vs. entire annulus
			P1[2] = new_pt[2] - midpt[2];
			P1[3] = 0.;
			//allow sprouts to grow outside the tissue region but within overall domain
			if(length0(P1) < IR || length0(P1) > OR) goto nonewseg;
		}
		else if(hexperiodic){
			if(new_pt[3] <= 0. || new_pt[3] >= alz) goto nonewseg;
		}
		else{				//test against included tissue points
			if(new_pt[1] <= 0. || new_pt[1] >= alx) goto nonewseg;
			if(new_pt[2] <= 0. || new_pt[2] >= aly) goto nonewseg;
			if(new_pt[3] <= 0. || new_pt[3] >= alz) goto nonewseg;
			ipt = (new_pt[1]/alx)*mxx + 1;
			jpt = (new_pt[2]/aly)*myy + 1;
			kpt = (new_pt[3]/alz)*mzz + 1;
			i = nbou[ipt][jpt][kpt];
			if (i == 0) goto nonewseg;	//no new segments outside tissue region
		}
//find starting node of new segment
		inod_prev = ista[iseg_from];
		inod_from = iend[iseg_from];
		interceptflag = intercept_check(start_pt, new_pt,int_pt,inod_from);
		if(interceptflag > 0){
			printf(" %i", interceptflag);
//find or create ending node of new segment.  Tests ending node for viability and how close to new node		
			iseg_to = interceptflag;
			inod1 = ista[iseg_to];
			inod2 = iend[iseg_to];
			for(i=1; i<=3; i++){
				P1[i] = cnode[i][ista[iseg_to]];
				P2[i] = cnode[i][iend[iseg_to]];
			}
			dist1 = distance(int_pt,P1);
			dist2 = distance(int_pt,P2);
			if(dist1 < tolerance_2 || dist2 < tolerance_2){//test whether intercept is too close to end nodes
				if(dist1 < dist2){
					inod_to = inod1;
					for(i=1; i<=3; i++) int_pt[i] = P1[i];
				}
				else{
					inod_to = inod2;
					for(i=1; i<=3; i++) int_pt[i] = P2[i];
				}			
				for(inodbc=1; inodbc<=nnodbc; inodbc++) if(inod_to == bcnod[inodbc]) goto nonewseg;//Intercepted an I/O node
				if(inod_to == inod_prev || iseg_to == iseg_from) goto nonewseg;//New segment folded back
			}
			else newnode(int_pt, iseg_to);	//split existing segment iseg_to into two segments, iseg_to and iseg
			seg_length = distance(start_pt,int_pt);
			if(seg_length < tolerance_2){	//extend existing segment to int_pt (no new segments, one node removed)
				iseg = nseg;
				inod = nnod;
				if(dist1 >= tolerance_2 && dist2 >= tolerance_2) inod_to =nnod;
				iend[iseg_from] = inod_to;
				segnodname[2][iseg_from] = nodname[inod_to];
				nodtyp[inod_to]++;
				if(nodtyp[inod_to] > nodsegm) printf("*** Error: Too many segments connected to node %i\n", inod_to);
				nodseg[nodtyp[inod_to]][inod_to] = iseg_from;
				nodnod[nodtyp[inod_to]][inod_to] = inod_prev;
				nodtyp[inod_from] = 0;
				segtyp[iseg_from] = 1;
			}
			else{	//add new segment to intersection point (no new nodes)
				iseg = nseg;
				inod = nnod;
				iseg++;
				segnamemax++;
				segname[iseg] = segnamemax;
				if(dist1 >= tolerance_2 && dist2 >= tolerance_2) inod_to = nnod;
				nodtyp[inod_to]++;
				if(nodtyp[inod_to] > nodsegm) printf("*** Error: Too many segments connected to node %i\n", inod_to);
				nodseg[nodtyp[inod_to]][inod_to] = iseg;
				nodnod[nodtyp[inod_to]][inod_to] = inod_from;
				nodtyp[inod_from]++;
				if(nodtyp[inod_to] > nodsegm) printf("*** Error: Too many segments connected to node %i\n", inod_to);
				nodseg[nodtyp[inod_from]][inod_from] = iseg;
				nodnod[nodtyp[inod_from]][inod_from] = inod_to;
				ista[iseg] = inod_from;		
				segnodname[1][iseg] = nodname[inod_from];
				iend[iseg] = inod_to;
				segnodname[2][iseg] = nodname[inod_to];
				diam[iseg] = new_diam;
				if(newsprout == 1) ksseg[iseg] = ks1 + randgauss() * ranks1;	//for a new sprout, choose ks from a Gaussian distribution
				else ksseg[iseg] = ksseg[iseg_from];	//continue with same ks
				q[iseg] = 0.;
				hd[iseg] = 0.;
				segtyp[iseg] = 1;
				nseg = iseg;
				if(nseg > nseg_new) printf("*** Error: nseg exceeds limit\n");
				nsplus++;
			}
		}
		else{	//continues with new active segment - no intercept - one new node, one new segment
			inod++;
			iseg++;
			nodnamemax++;
			segnamemax++;
			nodname[inod] = nodnamemax;
			segname[iseg] = segnamemax;
			for(i=1; i<=3; i++) cnode[i][inod] = new_pt[i];
			segtyp[iseg] = 0;
			segtyp[iseg_from] = 1;
			iend[iseg] = inod;
			segnodname[2][iseg] = nodname[inod];
			ista[iseg] = inod_from;
			segnodname[1][iseg] = nodname[inod_from];
			nodtyp[inod_from]++;
			if(nodtyp[inod_from] > nodsegm) printf("*** Error: Too many segments connected to node %i\n", inod_from);
			nodseg[nodtyp[inod_from]][inod_from] = iseg;
			nodnod[nodtyp[inod_from]][inod_from] = inod;
			nodtyp[inod] = 1;
			nodseg[nodtyp[inod]][inod] = iseg;
			nodnod[nodtyp[inod]][inod] = inod_from;
			diam[iseg] = new_diam;
			if(newsprout) ksseg[iseg] = ks1 + randgauss()*ranks1;	//for a new sprout, choose ks from a Gaussian distribution
			else ksseg[iseg] = ksseg[iseg_from];	//continue with same ks
			q[iseg] = 0.;
			hd[iseg] = 0.;
			nnod = inod;
			if(nnod > nnod_new) printf("*** Error: nnod exceeds limit\n");
			nseg = iseg;
			if(nseg > nseg_new) printf("*** Error: nseg exceeds limit\n");
			nsplus++;
		}
		nonewseg:;
	}
	free_vector(start_pt,1,3);
	free_vector(new_pt,1,3);
	free_vector(int_pt,1,3);
}