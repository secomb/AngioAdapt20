/****************************************************************************************
sht_seg_fixer.cpp
An unwanted node (vertex node) is moved to another node (inod_to) to eliminate a short segment.
Short seg is made type 10. The segments attached to the vertex node are moved to the inod_to node.
1. change ista, iend node
2. update nodseg and nodnod at inod_to
3. update nodseg and nodnod other nodes on segments
4. avoid duplicates or type 10 in nodnod/nodseg
5. nodtyp[vertex_node] = 0 (zero out nodseg, nodnod)
Tension sends inod_to and vertex node.  
Prep_for_adapt sends inod_to and notinod_to.
JA, 9/2008. Reworked by TWS, November 2010
*******************************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "nrutil.h"

void sht_seg_fixer(int inod_to, int seg_no, int vertex_node)
{
	extern int nseg,nnod;
	extern int *nodname,*segtyp,*ista,*iend,*nodtyp,**nodnod,**nodseg,**segnodname;

	int i,iseg,inew; 

	printf(" %i", seg_no);
	segtyp[seg_no] = 10;
	//segments connnected to inod_to: update nodnod and nodseg to remove segment
	inew = 0;
	for(i=1; i<=nodtyp[inod_to]; i++){
		iseg = nodseg[i][inod_to];
		if(segtyp[iseg] != 10){
			inew++;			//inew <= i so values will not be overwritten
			nodseg[inew][inod_to] = nodseg[i][inod_to];
			nodnod[inew][inod_to] = nodnod[i][inod_to];
		}
	}
	//segments connnected to vertex node: connect to inod_to, update nodnod and nodseg
	for(i=1; i<=nodtyp[vertex_node]; i++){
		iseg = nodseg[i][vertex_node];
		if(segtyp[iseg] != 10){
			inew++;
			nodseg[inew][inod_to] = nodseg[i][vertex_node];
			nodnod[inew][inod_to] = nodnod[i][vertex_node];
			if(ista[iseg] == vertex_node){
				ista[iseg] = inod_to;
				segnodname[1][iseg] = nodname[inod_to];
			}
			else{
				iend[iseg] = inod_to;
				segnodname[2][iseg] = nodname[inod_to];
			}
		}
	}
	nodtyp[inod_to] = inew;
	nodtyp[vertex_node] = 0;
}
