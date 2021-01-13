/********************************************************************************
intercept_check.cpp - increment an active segment
Uses close_pt_finder to test for proximity to existing segments
If an intersection of a segment is found, it returns a flag, segment number,
coordinates of intersection.
Flag:	no intercept: -1
		intercept: intersected segment number > 0
		intercept within distance < tol1: -1000 (short segment)
increment: amount seg increments by then checks for intercept (5.0)
tolerance_1: if increments to within this distance of new_point, stops incrementing.
			If it encounters a close segment (within tol1) while incrementing, intercepts
			the segment at its closest point.

Comments updated J. Alberding May 09
*********************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

int close_pt_finder(float* new_pt, float* int_pt, int inod_from, int iseg_from);
float distance(float* x, float* y);

int intercept_check(float* start_pt, float* new_pt, float* int_pt, int inod_from)
{
	extern float tolerance_1, increment;
	int i, intercept_flag;
	float more_inc, t_inc, d_new_start;
	float* new_new_pt;

	new_new_pt = vector(1, 3);
	d_new_start = distance(new_pt, start_pt);
	t_inc = increment / d_new_start;
	//more_inc varies from 0 at start_pt to 1 at new_pt
	more_inc = t_inc;	//don't start at very beginning of vessel - TWS July 2016
	intercept_flag = 0;
	while (more_inc <= 1.) {
		if ((1. - more_inc) * d_new_start < tolerance_1) more_inc = 1.;
		for (i = 1; i <= 3; i++) new_new_pt[i] = (1. - more_inc) * start_pt[i] + more_inc * new_pt[i];
		intercept_flag = close_pt_finder(new_new_pt, int_pt, inod_from, 0);
		if (intercept_flag == 0) more_inc += t_inc;
		else goto intercept;
	}
intercept:;
	free_vector(new_new_pt, 1, 3);
	return intercept_flag;
}
