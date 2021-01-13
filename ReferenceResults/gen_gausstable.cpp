/*******************************************************
randgauss.cpp - program to generate a variable from a Gaussian
distribution with s.d. = 1 using the polar or Box-Muller method.
See: C++ for Mathematicians: An Introduction for Students and
Professionals by Edward R. Scheinerman.  TWS - January 08
gen_gausstable is code to test distribution produced by randgauss
*******************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

float randgauss()
{
	float x, y, r, mu;
	do {
		x = 2. * rand() / float(RAND_MAX) - 1.;
		y = 2. * rand() / float(RAND_MAX) - 1.;
		r = x * x + y * y;
	} while (r >= 1.);
	mu = sqrt(-2.0 * log(r) / r);
	return mu * x;
}
/*
//this can be used for testing randgauss
void gen_gausstable()
{
	int i,ixg,resg[100];
	float xg;
	for(i=0; i<100; i++) resg[i] = 0;
	for(i=1; i<=1000; i++){
		xg = randgauss();
		if(xg > 0.) ixg = 10.*xg + 0.5;
		else  ixg = 10.*xg - 0.5;
		if(ixg > 49) ixg = 49;
		if(ixg < -49) ixg = -49;
		resg[ixg+50]++;
	}
	return;
}
*/