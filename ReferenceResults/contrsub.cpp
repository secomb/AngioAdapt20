/*******************************************************************
C version of subroutine contr, contour program.
Data file should have following structure (separated by spaces or line feeds):
m n nl:  dimensions of array, no. of contour levels
xmin xmax ymin ymax:  boundaries of box
cl(nl):  array of contour levels
zv(m,n,nsp):  array of heights zv(1,1) zv(1,2) ... zv(1,n) zv(2,1) ...
Output to a postscript file.
TWS, November 1989.
Converted to C September 2007
*********************************************************************/

#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void contr(FILE* ofp, int isp, int m, int n, float scalefac, int nl,
	float xmin, float xmax, float ymin, float ymax, float* cl, float** zv)
{
	extern float thresh;

	FILE* ofptest;
	int i, j, k, iwsp, in, in1, nsq, nsq1, isq, isq1;
	int** nwsp;
	float c, d, di1, dj1;
	float* xv, * yv;
	float** wsp;
	const int ncmax = 10000;
	xv = vector(1, m);
	yv = vector(1, n);
	wsp = matrix(1, ncmax, 1, 2);
	nwsp = imatrix(1, ncmax, 1, 2);
	for (i = 1; i <= m; i++) xv[i] = xmin + (i - 1) * (xmax - xmin) / (m - 1);
	for (j = 1; j <= n; j++) yv[j] = ymin + (j - 1) * (ymax - ymin) / (n - 1);
	fprintf(ofp, "/mx {%g mul 10 add} def\n", scalefac);
	fprintf(ofp, "/my {%g mul 10 add} def\n", scalefac);
	fprintf(ofp, "/m {moveto} def\n");
	fprintf(ofp, "/l {lineto} def\n");
	fprintf(ofp, "newpath\n");
	fprintf(ofp, "%g mx %g my m\n", xmin, ymin);
	fprintf(ofp, "%g mx %g my l\n", xmax, ymin);
	fprintf(ofp, "%g mx %g my l\n", xmax, ymax);
	fprintf(ofp, "%g mx %g my l\n", xmin, ymax);
	fprintf(ofp, "closepath\n");
	fprintf(ofp, "stroke\n");
	fprintf(ofp, "0.5 setlinewidth\n");
	fprintf(ofp, "/Times-Roman findfont\n");
	fprintf(ofp, "8 scalefont\n");
	fprintf(ofp, "setfont\n");
	ofptest = fopen("ltthresh.out", "a");
	//Begin contour generation
	for (k = 1; k <= nl; k++) {
		iwsp = 0;
		c = cl[k];
		if (c <= thresh) fprintf(ofptest, "c = %g, under thresh\n", c);
		for (i = 1; i <= m; i++)	for (j = 1; j <= n; j++) {
			//Look at vertical segments of mesh to detect crossing points.  Label with index numbers of adjacent two squares
			d = zv[i][j] - c;
			if (j != n) {
				dj1 = zv[i][j + 1] - c;
				if ((d >= 0 && dj1 < 0) || (d < 0 && dj1 >= 0)) {
					iwsp++;
					if (iwsp > ncmax) printf("*** Error: ncmax too small in contr\n");
					wsp[iwsp][1] = xv[i];
					wsp[iwsp][2] = (d * yv[j + 1] - dj1 * yv[j]) / (d - dj1);
					if (i == 1) nwsp[iwsp][1] = 0;
					else nwsp[iwsp][1] = i - 1 + (j - 1) * m;
					if (i == m) nwsp[iwsp][2] = 0;
					else nwsp[iwsp][2] = i + (j - 1) * m;
				}
			}
			//Look at horizontal segments
			if (i != m) {
				di1 = zv[i + 1][j] - c;
				if ((d >= 0 && di1 < 0) || (d < 0 && di1 >= 0)) {
					iwsp++;
					wsp[iwsp][1] = (d * xv[i + 1] - di1 * xv[i]) / (d - di1);
					wsp[iwsp][2] = yv[j];
					if (j == 1) nwsp[iwsp][1] = 0;
					else nwsp[iwsp][1] = i + (j - 2) * m;
					if (j == n) nwsp[iwsp][2] = 0;
					else nwsp[iwsp][2] = i + (j - 1) * m;
				}
			}
		}
		//Work through crossing points in sequence.  If nwsp(1 or 2) is not zero, find another point with matching value.
		//Join them and set nwsp values to zero at both.  Continue around loops.
		for (in = 1; in <= iwsp; in++) for (isq = 1; isq <= 2; isq++) {
			nsq = nwsp[in][isq];
			if (nsq != 0) {
				fprintf(ofp, "newpath\n");
				if (c <= thresh && isp == 2) fprintf(ofp, "0 5 1 setrgbcolor\n");
				else if (isp == 2) fprintf(ofp, "0 1 0 setrgbcolor\n");
				fprintf(ofp, "%g mx %g my m\n", wsp[in][1], wsp[in][2]);
				fprintf(ofp, "(%g) show\n", c);  //label contour
				fprintf(ofp, "%g mx %g my m\n", wsp[in][1], wsp[in][2]); //return to starting point
				nwsp[in][isq] = 0;
				do {
					for (in1 = 1; in1 <= iwsp; in1++) {
						for (isq1 = 1; isq1 <= 2; isq1++) {
							nsq1 = nwsp[in1][isq1];
							if (nsq1 == nsq && nsq != 0) {
								fprintf(ofp, "%g mx %g my l\n", wsp[in1][1], wsp[in1][2]);
								nwsp[in1][isq1] = 0;
								nsq = nwsp[in1][3 - isq1];
								nwsp[in1][3 - isq1] = 0;
							}
						}
					}
				} while (nsq != 0);
				fprintf(ofp, "stroke\n");
			}
		}
	}
	fclose(ofptest);
	free_vector(xv, 1, m);
	free_vector(yv, 1, n);
	free_matrix(wsp, 1, ncmax, 1, 2);
	free_imatrix(nwsp, 1, ncmax, 1, 2);
}