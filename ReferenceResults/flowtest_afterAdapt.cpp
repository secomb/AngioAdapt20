/************************************************************************
flowtest_afterAdapt - for AngioAdapt07
Use linear solver to identify and remove segments without flow
Set all diameters and viscosities to constant values for this calculation
Type 0 and type 1 segments are not removed, unless they are disconnected segment(s) from other segments.
	-this is the version that excludes/eliminates segments.  Compare to flowtest_afterAngio where there
	is no elimination of segments.
Use local pressure array to avoid disturbing actual values between iterations
TWS January 08,comments updated J. Alberding May 2009
Revised TWS2010
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void solve(double* press);

void flowtest_afterAdapt()
{
	extern int nseg, nnod, nnodbc, nsminus;
	extern int* segtyp, * ista, * iend, * bcnod, * nodtyp;
	extern int** nodseg;
	extern float constvisc, facfp;
	extern float* diam, * lseg, * cond;

	int iseg, inod, inodbc, ii, nodtyp_ista, nodtyp_iend, flag, loop_control;
	double* press, * q;	//make q and press local to avoid interference with actual flows
	press = dvector(1, nnod);
	q = dvector(1, nseg);

	for (inod = 1; inod <= nnod; inod++) press[inod] = 50.;
	for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] != 10) cond[iseg] = facfp * 1.e4 / lseg[iseg] / constvisc;
	solve(press);
	for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] != 10) q[iseg] = (press[ista[iseg]] - press[iend[iseg]]) * cond[iseg];
	for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] > 1 && segtyp[iseg] != 10) if (fabs(q[iseg]) < 0.0002) segtyp[iseg] = 1;

	//This handles "orphan" type 5 or 4 segments with no possible attachment to flowing segments
	flag = 1;
	loop_control = 0;
	while (flag) {
		flag = 0;
		for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] >= 3 && segtyp[iseg] <= 5) {
			for (inodbc = 1; inodbc <= nnodbc; inodbc++) if (ista[iseg] == bcnod[inodbc] || iend[iseg] == bcnod[inodbc]) goto skipseg;
			nodtyp_ista = nodtyp[ista[iseg]];
			nodtyp_iend = nodtyp[iend[iseg]];
			for (ii = 1; ii <= nodtyp[ista[iseg]]; ii++)
				if (segtyp[nodseg[ii][ista[iseg]]] != 4 && segtyp[nodseg[ii][ista[iseg]]] != 5) nodtyp_ista--;
			for (ii = 1; ii <= nodtyp[iend[iseg]]; ii++)
				if (segtyp[nodseg[ii][iend[iseg]]] != 4 && segtyp[nodseg[ii][iend[iseg]]] != 5) nodtyp_iend--;
			if (nodtyp_ista == 1 || nodtyp_iend == 1) {
				segtyp[iseg] = 1;
				flag = 1;
			}
		skipseg:;
		}
		//remove disconnected active segments (type 0 and 1) - find nodes with one segment
		//note this will not detect disconnected type 1 loops - needs more work
		for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] == 1 || segtyp[iseg] == 0) {
			nodtyp_ista = nodtyp[ista[iseg]];
			nodtyp_iend = nodtyp[iend[iseg]];
			for (ii = 1; ii <= nodtyp[ista[iseg]]; ii++) if (segtyp[nodseg[ii][ista[iseg]]] == 10) nodtyp_ista--;
			for (ii = 1; ii <= nodtyp[iend[iseg]]; ii++) if (segtyp[nodseg[ii][iend[iseg]]] == 10) nodtyp_iend--;
			if (nodtyp_ista == 1 || (segtyp[iseg == 1] && nodtyp_iend == 1)) {
				segtyp[iseg] = 10;
				diam[iseg] = 0.;
				flag = 1;
				nsminus++;
			}
		}
		loop_control++;
		if (loop_control > 150) {
			printf("*** Warning: adapt remove disconnected segs while loop exceeded loop control value(>150)\n");
			flag = 0;
		}
	}
	free_dvector(press, 1, nnod);
	free_dvector(q, 1, nseg);
}
