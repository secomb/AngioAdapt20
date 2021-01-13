/************************************************************************
putrank - generate list of nodes in order of flow direction
Considers only type 4 and 5 segments
nodrank --- if nodrank[i] < nodrank[j], node j is not upstream of node i
Also classifies nodes as converging (2), diverging (1) or other (0)
*************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void putrank(void)
{
	extern int nseg, nnod, nnodfl, nsegfl;
	extern int* nodrank, * nodtyp, * nodout, * nodname, * segtyp, * nodclasstyp, * nk, * ista, * iend;
	extern int** nodseg, ** nodnod;
	extern float* q, * nodvar, * segvar;

	int inod, j, iseg, nod1, nod2, flag;

	//construct node table; count outputs from node; output segments precede input segments
	for (inod = 1; inod <= nnod; inod++) {
		nodtyp[inod] = 0;
		nodout[inod] = 0;
	}
	nsegfl = 0;	//added TWS 2010
	for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] >= 3 && segtyp[iseg] <= 5) {	//output segments
		if (q[iseg] >= 0) {
			nod1 = ista[iseg];
			nod2 = iend[iseg];
		}
		else {
			nod1 = iend[iseg];
			nod2 = ista[iseg];
		}
		nodtyp[nod1]++;
		nodseg[nodtyp[nod1]][nod1] = iseg;
		nodnod[nodtyp[nod1]][nod1] = nod2;
		nodout[nod1]++;
		nsegfl++;
	}
	for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] >= 3 && segtyp[iseg] <= 5) {	//input segments
		if (q[iseg] >= 0) {
			nod1 = ista[iseg];
			nod2 = iend[iseg];
		}
		else {
			nod1 = iend[iseg];
			nod2 = ista[iseg];
		}
		nodtyp[nod2]++;
		nodseg[nodtyp[nod2]][nod2] = iseg;
		nodnod[nodtyp[nod2]][nod2] = nod1;
	}

	//assign low ranks to inflow nodes
	nnodfl = 0;
	for (inod = 1; inod <= nnod; inod++) {
		nk[inod] = 0;
		if (nodtyp[inod] == 1 && nodout[inod] == 1) {
			nnodfl++;
			nk[inod] = 1;
			nodrank[nnodfl] = inod;
		}
	}
	//assign increasing ranks to downstream connected nodes
	flag = 1;
	while (flag) {
		flag = 0;
		for (inod = 1; inod <= nnod; inod++) if (nk[inod] == 0 && nodtyp[inod] > 0) {
			for (j = nodout[inod] + 1; j <= nodtyp[inod]; j++) {
				iseg = nodseg[j][inod];	//inflow segment to node
				//process only nodes with all upstream nodes already ranked
				if (inod == iend[iseg] && (nk[ista[iseg]] == 0 || q[iseg] <= 0.)) goto skipnode;
				if (inod == ista[iseg] && (nk[iend[iseg]] == 0 || q[iseg] >= 0.)) goto skipnode;
			}
			nnodfl++;
			nk[inod] = 1;
			nodrank[nnodfl] = inod;
			flag = 1;
		skipnode:;
		}
	}
	for (iseg = 1; iseg <= nseg; iseg++) segvar[iseg] = iseg;
	for (inod = 1; inod <= nnod; inod++) nodvar[inod] = nk[inod];
	for (inod = 1; inod <= nnod; inod++) if (nodtyp[inod] != 0 && nk[inod] == 0)
		printf("*** Error: unprocessed node %i in putrank\n", inod);

	//setup nodclasstype for finding "mesh" segments
	for (inod = 1; inod <= nnod; inod++) {
		nodclasstyp[inod] = 0;
		if (nodtyp[inod] > 2) {
			if (nodout[inod] == nodtyp[inod] - 1) nodclasstyp[inod] = 1;	//diverging
			if (nodout[inod] == 1) nodclasstyp[inod] = 2;	//converging			
		}
		else if (nodtyp[inod] == 1) {
			if (nodout[inod] == 1) nodclasstyp[inod] = 1;	//classify inflow boundary node as diverging
			if (nodout[inod] == 0) nodclasstyp[inod] = 2;	//classify outflow boundary node as converging
		}
	}
}
