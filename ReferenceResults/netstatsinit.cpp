/**********************************************************
netstatsinit.cpp - initialize output files for network statistics.
Also initialize postscript graph files.
TWS, Feb. 2017
***********************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void netstatsinit()
{
	extern int imainmax;

	int i;
	float scalefacx, scalefacy, xmax, ymax;
	FILE* ofp;

	ofp = fopen("LengthVolumeStats1.out", "w");
	fprintf(ofp, "imain   time  Length tot Lengthfl Total segs Total flow segs Total vol Total flow vol nnv\n");
	fclose(ofp);
	ofp = fopen("LengthVolumeStats2.out", "w");
	fprintf(ofp, "imain   time  Length tot ALengthfl VLengthfl CLengthfl MLengthfl otherLengthfl nnv nsplus  nsminus  nsspr  nspruned\n");
	fclose(ofp);
	ofp = fopen("NodtypStats.out", "w");
	fprintf(ofp, "imain   time  nodtyp: 1    2    3    4    5    6 or more\n");
	fclose(ofp);
	ofp = fopen("OxygenStats.out", "w");
	fprintf(ofp, "imain   time  Demand Consumption Inflow Extraction%% Hypoxic%%\n");
	fclose(ofp);

	xmax = imainmax;
	scalefacx = 400. / xmax;
	ymax = 100.;
	scalefacy = 500. / ymax;
	ofp = fopen("OxygenStats.ps", "w");
	fprintf(ofp, "%%!PS\n");
	fprintf(ofp, "/mx {%g mul 50 add} def\n", scalefacx);
	fprintf(ofp, "/my {%g mul 50 add} def\n", scalefacy);
	fprintf(ofp, "/m {moveto} def\n");
	fprintf(ofp, "/l {lineto} def\n");
	fprintf(ofp, "/cf {1 0 360 arc closepath fill} def\n");
	fprintf(ofp, "newpath\n");
	fprintf(ofp, "%g mx %g my m\n", 0., ymax);
	fprintf(ofp, "%g mx %g my l\n", 0., 0.);
	fprintf(ofp, "%g mx %g my l\n", xmax, 0.);
	fprintf(ofp, "stroke\n");
	fprintf(ofp, "/Times-Roman findfont 10 scalefont setfont\n");
	for (i = 0; i <= 10; i++) fprintf(ofp, "30 %f my m (%4.1f) show\n", i / 10. * ymax, i / 10. * ymax);	//y-axis labels
	for (i = 0; i <= 10; i++)	fprintf(ofp, "%f mx 40 m (%3.0f) show\n", i / 10. * xmax, i / 10. * xmax);	//x-axis labels
	fprintf(ofp, "stroke\n");
	fclose(ofp);

	ymax = 50.;
	scalefacy = 500. / ymax;
	ofp = fopen("Vesseltypes.ps", "w");
	fprintf(ofp, "%%!PS\n");
	fprintf(ofp, "/mx {%g mul 50 add} def\n", scalefacx);
	fprintf(ofp, "/my {%g mul 50 add} def\n", scalefacy);
	fprintf(ofp, "/m {moveto} def\n");
	fprintf(ofp, "/l {lineto} def\n");
	fprintf(ofp, "/cf {1 0 360 arc closepath fill} def\n");
	fprintf(ofp, "newpath\n");
	fprintf(ofp, "%g mx %g my m\n", 0., ymax);
	fprintf(ofp, "%g mx %g my l\n", 0., 0.);
	fprintf(ofp, "%g mx %g my l\n", xmax, 0.);
	fprintf(ofp, "stroke\n");
	fprintf(ofp, "/Times-Roman findfont 10 scalefont setfont\n");
	for (i = 0; i <= 10; i++) fprintf(ofp, "30 %f my m (%4.1f) show\n", i / 10. * ymax, i / 10. * ymax);	//y-axis labels
	for (i = 0; i <= 10; i++)	fprintf(ofp, "%f mx 40 m (%3.0f) show\n", i / 10. * xmax, i / 10. * xmax);	//x-axis labels
	fprintf(ofp, "stroke\n");	//needed to avoid artifact in plot
	fclose(ofp);
}