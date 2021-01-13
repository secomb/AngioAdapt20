/**********************************************************************************
adaptsignals.cpp
Program to create a plot of all adaptation signals vs. segment pressure.
TWS Feb. 2017
*********************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "nrutil.h"

void adaptsignals(int runs)
{
	extern int nseg;
	extern int *segtyp;
	extern float *diam,*segpress,**signalvar;	//added Feb. 2017

	int i,iseg;
	float scalefacx,scalefacy,xmax,ymax,xmin,ymin,xvar;
	char fname[80];
	FILE *ofp;

	xmax = 100.;
	xmin = 0.;
	scalefacx = 400./xmax;
	ymax = 4.;
	ymin = -4.;
	scalefacy = 500./(ymax - ymin);

	sprintf(fname,"Current\\AdaptSignals%03i.ps",runs);
	ofp = fopen(fname, "w");
	fprintf(ofp, "%%!PS\n");
	fprintf(ofp, "/mx {%g sub %g mul 50 add} def\n",xmin,scalefacx);
	fprintf(ofp, "/my {%g sub %g mul 50 add} def\n",ymin,scalefacy);
	fprintf(ofp, "/m {moveto} def\n");
	fprintf(ofp, "/l {lineto} def\n");
	fprintf(ofp, "/cf {1.5 0 360 arc closepath fill} def\n");
	fprintf(ofp, "newpath\n");
	fprintf(ofp, "%g mx %g my m\n",0.,ymin);	//y-axis
	fprintf(ofp, "%g mx %g my l\n",0.,ymax);
	fprintf(ofp, "stroke\n");
	fprintf(ofp, "%g mx %g my m\n",0.,0.);	//x-axis
	fprintf(ofp, "%g mx %g my l\n",xmax,0.);
	fprintf(ofp, "stroke\n");
	fprintf(ofp, "/Times-Roman findfont 10 scalefont setfont\n");
	for(i=0; i<=10; i++) fprintf(ofp,"%g mx %g my m (%4.1f) show\n",-xmax*0.05, ymin+i/10.*(ymax-ymin),ymin+i/10.*(ymax-ymin));	//y-axis labels
	for(i=0; i<=10; i++) fprintf(ofp,"%g mx %g my m (%3.0f) show\n", i/10.*xmax,-(ymax-ymin)*0.02,i/10.*xmax);	//x-axis labels
	fprintf(ofp, "stroke\n");
	for(iseg=1; iseg<=nseg; iseg++) if(segtyp[iseg] >= 3 && segtyp[iseg] <= 5){
		xvar = segpress[iseg];
		//xvar = diam[iseg]*4.;
		fprintf(ofp, "1 0 0 setrgbcolor\n");
		fprintf(ofp, "%g mx %g my cf\n",xvar,signalvar[iseg][1]);	//red - shear
		fprintf(ofp, "0 1 0 setrgbcolor\n");
		fprintf(ofp, "%g mx %g my cf\n",xvar,signalvar[iseg][2]);	//green - pressure
		fprintf(ofp, "0 0 1 setrgbcolor\n");
		fprintf(ofp, "%g mx %g my cf\n",xvar,signalvar[iseg][3]);	//blue - shrinking
		fprintf(ofp, "0.5 0 0.5 setrgbcolor\n");
		fprintf(ofp, "%g mx %g my cf\n",xvar,signalvar[iseg][4]);	//purple - upstream conducted
		fprintf(ofp, "0.5 0.5 0 setrgbcolor\n");
		fprintf(ofp, "%g mx %g my cf\n",xvar,signalvar[iseg][5]);	//olive - downstream convected
		fprintf(ofp, "0 0 0 setrgbcolor\n");
		fprintf(ofp, "%g mx %g my cf\n",xvar,signalvar[iseg][6]);	//black - total
	}
	fprintf(ofp, "showpage\n");
	fclose(ofp);
}
