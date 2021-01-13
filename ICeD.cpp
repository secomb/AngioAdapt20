/****************************************************************************
ICeD.cpp
Integrated conical detector to simulate filopodia detection by a growing vessels.
Program to integrate signal contribution from vessels falling within a conical
section about the forward direction/start node of an actively growing segment.
The integral used is:
S = int[r(s)/|r(s)|*f(r(s))*g(theta(s))*ds
r(s) = vector from start pt. of active segment (actseg) to a point on the line segment
f(r(s)) = 1 - r(s)/r0, r0 = radius of cone sector = 50 microns
g(theta(s)) = 1 - (1 - cos theta)/(1 - cos thet0), theta_0 = max cone angle = pi/3.
direct = current direction
xp = endpoint of current segment
signal = vector representing signal
Jon Alberding, May 08 2006. Comments updated J. Alberding, May 09.  Revised TWS 2010
*************************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

float length(float* x);
float length0(float* x);
float distance(float* x, float* y);

float ICeD(float* direct, float* xp, float* signal, float r0, float theta0)
{
	extern int nseg, * ista, * iend, * segtyp;
	extern float** cnode, * P1, * P2, * P3;
	int i, ii, iseg, n;
	float length1, dist0, dist1, magd, magr, magrd, costheta, thetafac, rfac, intfac;
	float* rs;

	rs = vector(1, 3);

	magd = length0(direct);  //magnitude of current direction vector, moved out of loop, Feb. 2017

	for (i = 1; i <= 3; i++) signal[i] = 0.;
	//exclude types 3, 4 and 10. Feb. 2017
	for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] != 3 && segtyp[iseg] != 4 && segtyp[iseg] != 10) {
		for (i = 1; i <= 3; i++) {
			P1[i] = cnode[i][ista[iseg]] - xp[i];	//start point of target segment
			P2[i] = cnode[i][iend[iseg]] - xp[i];	//end point of target segment
			P3[i] = cnode[i][iend[iseg]] - cnode[i][ista[iseg]];
		}
		dist0 = length(P1);	//distance from start point of target segment, may change P1 if hexperiodic
		dist1 = length(P2);	//distance from end point of target segment, may change P2 if hexperiodic
		length1 = length(P3);	//length of segment, may change P3 if hexperiodic
		//prevent cone detector from detecting distant segments; don't use if already at the segment
		if ((dist0 < 2. * r0 || dist1 < 2. * r0) && dist0 > 0.1 && dist1 > 0.1) { //modified TWS, Feb. 2018
			//discretize line segment into n segments (5 micron each) - modified TWS, Feb. 2018
			n = IMAX(length1 / 5. + 0.5, 1);	//modfied Feb. 2018, was IMIN
			for (ii = 0; ii <= n; ii++) {
				magrd = 0.;
				for (i = 1; i <= 3; i++) rs[i] = P1[i] + ii * P3[i] / n;	//vector from current point to point on vessel
				magr = length(rs);	//distance from point on segment, may change rs if hexperiodic
				for (i = 1; i <= 3; i++) magrd += direct[i] * rs[i];	//dot product of rs and direct
				if (magr < r0 && magrd > 0.) {	//modified TWS, Feb. 2018
					costheta = magrd / magr / magd;	//cos(angle between growing segment and point on target segment)
					rfac = 1. - magr / r0;	//intensity of attraction
					thetafac = FMAX(1. - (1. - costheta) / (1. - cos(theta0)), 0.);	//define edge of cone
					intfac = 1. / n;
					if (ii == 0 || ii == n) intfac = 0.5 / n;
					for (i = 1; i <= 3; i++)	signal[i] += intfac * rs[i] / magr * rfac * thetafac;
				}
			}
		}
	}
	free_vector(rs, 1, 3);
	return 0.;
}
