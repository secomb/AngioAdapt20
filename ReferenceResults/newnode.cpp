/*********************************************************************
newnode.cpp
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

Numbering scheme for case when vessel sprouts from another

						  |
						  * inod1
		   new_pt start_pt|
   inod_to = inod         | iseg_to
------.....*--------------* inod_from
						  | iseg
						  |
						  * inod2
						  |

This program will handle the addition of a new node midsegment in sprouting or intersected segments.
It updates nodtype, nodseg, nodnod, iseg, inod, nseg, nnod, segname, cnode,
as it revises geometry for the old segment(iseg_to) and the newly formed segment (iseg) and the new node.
NOTE: Test of the new node and/or new segment for viability should be done before newnode is called.
J. Alberding, Sep. 19, 2006.  TWS January 08. Revised TWS2010
******************************************************************************/

#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

int newnode(float* int_pt, int iseg_to)
{
	extern int nseg, nnod, nodnamemax, segnamemax, nseg_new, nnod_new, nsplus;
	extern int* segname, * segtyp, * ista, * iend, * nodname, * nodtyp;
	extern int** nodseg, ** nodnod, ** segnodname;
	extern float* diam, * ksseg, * q, * hd;
	extern float** cnode;
	extern FILE* ofp2;

	int i, iseg, inod, inod1, inod2, inod_to;

	inod1 = ista[iseg_to];
	inod2 = iend[iseg_to];
	inod = nnod;
	iseg = nseg;
	//split existing segment iseg_to into two segments, iseg_to and iseg
	inod++;
	iseg++;
	nodnamemax++;
	segnamemax++;
	nodname[inod] = nodnamemax;
	segname[iseg] = segnamemax;
	nodtyp[inod] = 2;
	inod_to = inod;
	for (i = 1; i <= 3; i++) cnode[i][inod] = int_pt[i];
	nodseg[1][inod] = iseg_to;
	nodseg[2][inod] = iseg;
	nodnod[1][inod] = ista[iseg_to];
	nodnod[2][inod] = iend[iseg_to];
	for (i = 1; i <= nodtyp[inod1]; i++) if (nodnod[i][inod1] == iend[iseg_to]) nodnod[i][inod1] = inod;
	for (i = 1; i <= nodtyp[inod2]; i++) if (nodseg[i][inod2] == iseg_to) nodseg[i][inod2] = iseg;
	for (i = 1; i <= nodtyp[inod2]; i++) if (nodnod[i][inod2] == ista[iseg_to]) nodnod[i][inod2] = inod;
	ista[iseg] = inod;
	segnodname[1][iseg] = nodnamemax;
	iend[iseg] = iend[iseg_to];
	segnodname[2][iseg] = nodname[iend[iseg]];
	iend[iseg_to] = inod;
	segnodname[2][iseg_to] = nodnamemax;
	diam[iseg] = diam[iseg_to];
	ksseg[iseg] = ksseg[iseg_to];
	q[iseg] = q[iseg_to];
	hd[iseg] = hd[iseg_to];

	if (segtyp[iseg_to] == 0) {	//only end portion is growth tip
		segtyp[iseg] = 0;
		segtyp[iseg_to] = 1;
	}
	else segtyp[iseg] = segtyp[iseg_to];	//TWS 2010
	nnod = inod;
	if (nnod > nnod_new) printf("*** Error: nnod exceeds limit\n");
	nseg = iseg;
	if (nseg > nseg_new) printf("*** Error: nseg exceeds limit\n");
	nsplus++;
	return inod;
}
