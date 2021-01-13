/*****************************************************
distance.cpp
Evaluate distance between two points x,y in a periodic hexagonal structure
If the vector y - x lies outside the unit hexagon, find direction where it
is most far outside, then map back. If necessary, repeat procedure.
TWS February 2018
******************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

float length(float* x)
{
	extern int hexperiodic, flag, idomain;
	extern float aphex, aphexp, ** hex_norm;

	int i, idomainmax, idomainsave;
	float dotproduct, dotproductmax, length2, length3;

	if (hexperiodic) {		//test for point outside base domain
		flag = 1;
		idomainsave = 0;
		do {
			flag = 0;
			idomain = 0;
			idomainmax = 0;
			dotproductmax = 0.;
			for (idomain = 1; idomain <= 6; idomain++) {	//determine which boundary is crossed by the largest amount
				dotproduct = 0.;
				for (i = 1; i <= 2; i++) dotproduct += x[i] * hex_norm[idomain][i];
				if (dotproduct > aphexp && dotproduct > dotproductmax) {
					dotproductmax = dotproduct;
					idomainmax = idomain;
					flag = 1;
				}
				if (flag) idomainsave = idomainmax;	//save last domain crossing identifier found in loop
			}
			//note that returned vector x is altered from original
			if (idomainmax != 0) for (i = 1; i <= 2; i++) x[i] -= 2. * aphexp * hex_norm[idomainmax][i];
			//find corresponding point in domain closer to base domain
		} while (flag);
		idomain = idomainsave;	//this value is available after running length
	}
	length2 = SQR(x[1]) + SQR(x[2]);
	length3 = length2 + SQR(x[3]);
	length2 = sqrt(length2);
	if (hexperiodic && length2 > aphex)
		printf("*** error in length calculation, %f\n", length2);
	length3 = sqrt(length3);
	return length3;
}

float length0(float* x)
{
	int i;
	float length;
	length = 0.;
	for (i = 1; i <= 3; i++) length += SQR(x[i]);
	length = sqrt(length);
	return length;
}

float distance(float* x, float* y)
{
	int i;
	float* P3, length1;
	P3 = vector(1, 3);
	for (i = 1; i <= 3; i++) P3[i] = x[i] - y[i];
	length1 = length(P3);
	free_vector(P3, 1, 3);
	return length1;
}