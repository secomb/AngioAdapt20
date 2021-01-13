/***********************************************************************************
prep_for_adapt.cpp
1. Test for duplicate segments (diff iseg, same ista and iend nodes).
2. Test for reversed segments (diff iseg, ista and iend nodes reversed from another segment).
Tests 1 & 2 result in one duplicate or reversed segment being made into a type 10 seg.
3. Test for segments with duplicate iseg values (stops all activity, must be fixed).
4. Test for short segments (under tolerance_2 length) (move the node with the least number
of branch segments to the other end of the "short" segment).
J Alberding, Feb 2008.  Revised TWS January 2011
**************************************************************************************/
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include "nrutil.h"

void sht_seg_fixer(int inod_to, int seg_no, int vertex_node);
float length(float* P1);

void prep_for_adapt(int inod_to)
{
	extern int nseg, nnod, nnodbc, nsminus;
	extern int* segname, * nodname, * bcnodname, * segtyp, * ista, * iend, * nodtyp;
	extern int** nodnod, ** nodseg;
	extern float tolerance_2;
	extern float* diam, * P1;
	extern float** cnode;

	int i, ii, kk, same_nodes, seg_num, tseg_num, inew, notinod_to, flag;

	if (nodtyp[inod_to] == 0) return;

	do {
		flag = 0;
		for (ii = 1; ii < nodtyp[inod_to]; ii++) for (kk = ii + 1; kk <= nodtyp[inod_to]; kk++) {
			same_nodes = 0;
			seg_num = nodseg[ii][inod_to];
			tseg_num = nodseg[kk][inod_to];
			if (seg_num == tseg_num) {
				printf("*** Error: Duplicate numbered segments\n");
				return;
			}
			if (segtyp[seg_num] != 10 && segtyp[tseg_num] != 10) {
				same_nodes = 0;
				if (ista[tseg_num] == ista[seg_num] && iend[tseg_num] == iend[seg_num]) same_nodes = 1;
				if (ista[tseg_num] == iend[seg_num] && iend[tseg_num] == ista[seg_num]) same_nodes = 1;
			}
			if (same_nodes == 1) {
				printf(" dup%i/%i", segname[seg_num], segname[tseg_num]);
				segtyp[seg_num] = 10;	//remove segment seg_num
				diam[tseg_num] = FMAX(diam[seg_num], diam[tseg_num]);	//choose larger diameter
				nsminus++;
				inew = 0;
				for (i = 1; i <= nodtyp[iend[seg_num]]; i++) if (nodseg[i][iend[seg_num]] != seg_num) {	//re-evaluate nodnod and nodseg
					inew++;
					nodnod[inew][iend[seg_num]] = nodnod[i][iend[seg_num]];
					nodseg[inew][iend[seg_num]] = nodseg[i][iend[seg_num]];
				}
				nodtyp[iend[seg_num]]--;
				inew = 0;
				for (i = 1; i <= nodtyp[ista[seg_num]]; i++) if (nodseg[i][ista[seg_num]] != seg_num) {
					inew++;
					nodnod[inew][ista[seg_num]] = nodnod[i][ista[seg_num]];
					nodseg[inew][ista[seg_num]] = nodseg[i][ista[seg_num]];
				}
				nodtyp[ista[seg_num]]--;
			}
		}
		//check for short segments connected to inod_to. If any found, repeat above procedure
		for (ii = 1; ii <= nodtyp[inod_to]; ii++) {
			seg_num = nodseg[ii][inod_to];
			if (segtyp[seg_num] != 10) {
				for (i = 1; i <= 3; i++) P1[i] = cnode[i][ista[seg_num]] - cnode[i][iend[seg_num]];
				if (length(P1) <= tolerance_2) {	//note P1 may be changed if hexperiodic is active
					for (i = 1; i <= nnodbc; i++) if (bcnodname[i] == ista[seg_num] || bcnodname[i] == iend[seg_num]) {
						printf("*** Error: short segment %i attached to I/O node %i\n", seg_num, bcnodname[i]);
						return;
					}
					if (nodtyp[ista[seg_num]] < nodtyp[iend[seg_num]]) {
						notinod_to = ista[seg_num];
						inod_to = iend[seg_num];
					}
					else {
						notinod_to = iend[seg_num];
						inod_to = ista[seg_num];
					}
					sht_seg_fixer(inod_to, seg_num, notinod_to);
					flag = 1;
				}
			}
		}
	} while (flag);
	return;
}
