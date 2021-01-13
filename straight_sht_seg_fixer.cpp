/*************************************************************************************************
straight_sht_seg_fixer.cpp
This program examines the flowing (type 4 or 5) segments with endnodes of a nodtyp = 2. If the two
segments from the node are parallel, then:
   1. The routine estimates whether combining the two segments produces a segment with
	  a total length < or > 100 microns.
	  if total length > 100 microns, no combining of segments occurs.
	  if total length < 100 microns:
	2. The two segments are combined and the "new" segment is assigned the lower order of
		the original two segments.  The start/end node of the new segment are assigned
		from the  start/end nodes of the original segments.  the extra node is assigned
		a nodtyp = 0, and the original segment with the higher order is assigned a
		segtyp = 10.
	3. The diameter of the new segment is assigned the average of the original segment's diameters.
		The flow assigned is the higher of the flow values of the original segments.  The new segment
		segtyp = 5;
	4. Node tables and segnodname are updated by sht_seg_fixer.
inod = node to be removed if conditions listed above are satisfied
J. Alberding, March 19, 2009. Edited by TWS2010, 2018
**************************************************************************************************************/
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include "nrutil.h"

void sht_seg_fixer(int inod_to, int seg_no, int vertex_node);
float length(float* x);

void straight_sht_seg_fixer()
{
	extern int nseg, nnod, nnodbc;
	extern int* segname, * segtyp, * nodtyp, * bcnod, * ista, * iend;
	extern int** nodseg, ** nodnod;
	extern float* P1, * P2, * diam;
	extern float** cnode;

	int max_order, min_order, j, inod, inodbc, inod_to, inod1, inod2;
	float orthotest, total_length, length1, length2, dotproduct;

	for (inod = 1; inod <= nnod; inod++) {
		if (nodtyp[inod] == 2 && segtyp[nodseg[1][inod]] == 5 && segtyp[nodseg[2][inod]] == 5) {
			inod1 = nodnod[1][inod];
			inod2 = nodnod[2][inod];
			for (inodbc = 1; inodbc <= nnodbc; inodbc++)
				if (inod1 == bcnod[inodbc] || inod2 == bcnod[inodbc]) goto norecombine;
			for (j = 1; j <= 3; j++) {
				P1[j] = cnode[j][inod1] - cnode[j][inod];
				P2[j] = cnode[j][inod2] - cnode[j][inod];
			}
			length1 = length(P1);	//P1, P2 may be altered if hexperiodic is active
			length2 = length(P2);
			dotproduct = 0.;
			for (j = 1; j <= 3; j++) dotproduct += P1[j] * P2[j];
			total_length = length1 + length2;
			orthotest = 1. - fabs(dotproduct / length1 / length2);
			if (fabs(orthotest) < 1.e-4 && total_length < 100.) {
				max_order = nodseg[1][inod];
				min_order = nodseg[2][inod];
				if (min_order > max_order) {
					max_order = nodseg[2][inod];
					min_order = nodseg[1][inod];
				}
				if (ista[max_order] == inod) inod_to = iend[max_order];
				else inod_to = ista[max_order];
				sht_seg_fixer(inod_to, max_order, inod);
				diam[min_order] = (diam[min_order] + diam[max_order]) / 2.;
			}
		}
	norecombine:;
	}
	printf("\n");
}