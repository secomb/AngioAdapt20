/************************************************************************
setup_hexperiodic - for AngioAdapt17GPU.  TWS 2018
Variables neede for periodic hexagonal structure
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void setup_hexperiodic()
{
	extern float aphex, aphexp, alx, aly, * midpt, ** hex_norm, ** hex_tang;
	int i;

	midpt[1] = alx / 2.;
	midpt[2] = aly / 2.;
	midpt[3] = 0.;

	aphex = aly / 2.;
	aphexp = aphex * 0.5 * sqrt(3.);
	hex_norm = matrix(1, 6, 1, 3);
	hex_tang = matrix(1, 6, 1, 3);

	//normals to each side
	hex_norm[1][1] = 1.;
	hex_norm[1][2] = 0.;
	hex_norm[2][1] = 0.5;
	hex_norm[2][2] = 0.5 * sqrt(3.);
	hex_norm[3][1] = -0.5;
	hex_norm[3][2] = 0.5 * sqrt(3.);
	for (i = 1; i <= 3; i++) {
		hex_norm[i + 3][1] = -hex_norm[i][1];
		hex_norm[i + 3][2] = -hex_norm[i][2];
	}
	//tangents to each side
	for (i = 1; i <= 6; i++) {
		hex_norm[i][3] = 0.;
		hex_tang[i][1] = -hex_norm[i][2];
		hex_tang[i][2] = hex_norm[i][1];
		hex_tang[i][3] = 0.;
	}
}