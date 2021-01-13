/**********************************************************
signalplot.cpp - generate plots of adaptation signals vs. pressure.
TWS Jan. 08
***********************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void signalplot(const char fname[])
 {
	extern int nseg;
	extern int *segtyp;
	extern float ks1;
	extern float *segpress,*stau,*spress,*uptrans,*downtrans,*stot,*condup,*conddown;

	int iseg;
	float xs,ys,pmin,pmax,stotmax,stotmin,condmin,condmax,scalefacx,scalefacy;
	FILE *ofp;

	pmin = 1.e6;
	pmax = -1.e6;
	for(iseg=1; iseg<=nseg; iseg++) if(segtyp[iseg] >= 3 && segtyp[iseg] <= 5){
		if(segpress[iseg] > pmax) pmax = segpress[iseg];
		if(segpress[iseg] < pmin) pmin = segpress[iseg];
	}
	stotmin = -2.;
	stotmax = 2.;
	scalefacx = 300;
	scalefacy = 300;

	ofp = fopen(fname, "w");
	fprintf(ofp, "%%!PS\n");
	fprintf(ofp, "/mx {%g mul 150 add} def\n",scalefacx);
	fprintf(ofp, "/my {%g mul 500 add} def\n",scalefacy);
	fprintf(ofp, "/m {moveto} def\n");
	fprintf(ofp, "/l {lineto} def\n");
	fprintf(ofp, "/c {setrgbcolor} def\n");
	fprintf(ofp, "newpath\n");
	fprintf(ofp, "%g mx %g my m\n",0.,0.);
	fprintf(ofp, "%g mx %g my l\n",1.,0.);
	fprintf(ofp, "%g mx %g my l\n",1.,1.);
	fprintf(ofp, "%g mx %g my l\n",0.,1.);
	fprintf(ofp, "closepath\n");
	fprintf(ofp, "stroke\n");
	fprintf(ofp, "newpath\n");
	ys = (0. - stotmin)/(stotmax - stotmin);
	fprintf(ofp, "%g mx %g my m\n",0.,ys);
	fprintf(ofp, "%g mx %g my l\n",1.,ys);
	fprintf(ofp, "stroke\n");

	fprintf(ofp, "/Times-Roman findfont\n");
	fprintf(ofp, "10 scalefont\n");
	fprintf(ofp, "setfont\n");
	fprintf(ofp, "%g mx %g my m (%g) show\n",-0.05,0.,stotmin);
	fprintf(ofp, "%g mx %g my m (%g) show\n",-0.05,1.,stotmax);
	fprintf(ofp, "%g mx %g my m (%g) show\n",0.,-0.03,pmin);
	fprintf(ofp, "%g mx %g my m (%g) show\n",1.,-0.03,pmax);

	fprintf(ofp,"0.5 0.5 0.5 c\n");//gray
	ys = (ks1 - stotmin)/(stotmax - stotmin);
	fprintf(ofp, "%g mx %g my m\n",0.,ys);
	fprintf(ofp, "%g mx %g my l\n",1.,ys);
	fprintf(ofp, "stroke\n");
	fprintf(ofp,"0 0 0 c ");//black
	fprintf(ofp, "%g mx %g my m (stot) show\n", 0.05,0.01);
	fprintf(ofp,"1 0 0 c ");//red
	fprintf(ofp, "%g mx %g my m (stau) show\n", 0.2,0.01);
	fprintf(ofp,"0 1 0 c ");//green
	fprintf(ofp, "%g mx %g my m (spress) show\n", 0.35,0.01);
	fprintf(ofp,"0 0 1 c ");//blue
	fprintf(ofp, "%g mx %g my m (uptrans) show\n", 0.5,0.01);
	fprintf(ofp,"1 0 1 c ");//light green
	fprintf(ofp, "%g mx %g my m (downtrans) show\n", 0.65,0.01);

	for(iseg=1; iseg<=nseg; iseg++) if(segtyp[iseg] >= 3 && segtyp[iseg] <= 5){
		xs = (segpress[iseg] - pmin)/(pmax - pmin);
		fprintf(ofp,"0 0 0 c ");//black
		ys = (stot[iseg] - stotmin)/(stotmax - stotmin);
		fprintf(ofp, "%g mx %g my m ", xs + 0.004,ys);
		fprintf(ofp, "%g mx %g my 0.004 %g mul 0 360 arc fill\n", xs,ys,scalefacx);
		fprintf(ofp,"1 0 0 c ");//red
		ys = (stau[iseg] - stotmin)/(stotmax - stotmin);
		fprintf(ofp, "%g mx %g my m ", xs + 0.004,ys);
		fprintf(ofp, "%g mx %g my 0.004 %g mul 0 360 arc fill\n", xs,ys,scalefacx);
		fprintf(ofp,"0 1 0 c ");//green
		ys = (spress[iseg] - stotmin)/(stotmax - stotmin);
		fprintf(ofp, "%g mx %g my m ", xs + 0.004,ys);
		fprintf(ofp, "%g mx %g my 0.004 %g mul 0 360 arc fill\n", xs,ys,scalefacx);
		fprintf(ofp,"0 0 1 c ");//blue
		ys = (uptrans[iseg] - stotmin)/(stotmax - stotmin);
		fprintf(ofp, "%g mx %g my m ", xs + 0.004,ys);
		fprintf(ofp, "%g mx %g my 0.004 %g mul 0 360 arc fill\n", xs,ys,scalefacx);
		fprintf(ofp,"1 0 1 c ");//light green
		ys = (downtrans[iseg] - stotmin)/(stotmax - stotmin);
		fprintf(ofp, "%g mx %g my m ", xs + 0.004,ys);
		fprintf(ofp, "%g mx %g my 0.004 %g mul 0 360 arc fill\n", xs,ys,scalefacx);
	}
//show summed conducted signals before saturation relationship is applied
	condmin = 0.;
	condmax = 0.;
	for(iseg=1; iseg<=nseg; iseg++) if(segtyp[iseg] >= 3 && segtyp[iseg] <= 5){
		if(condup[iseg] > condmax) condmax = condup[iseg];
		if(conddown[iseg] > condmax) condmax = conddown[iseg];
	}
	if(condmax != 0){
		fprintf(ofp,"0 0 0 c ");//black
		fprintf(ofp, "/my {%g mul 100 add} def\n",scalefacy);
		fprintf(ofp, "newpath\n");
		fprintf(ofp, "%g mx %g my m\n",0.,0.);
		fprintf(ofp, "%g mx %g my l\n",1.,0.);
		fprintf(ofp, "%g mx %g my l\n",1.,1.);
		fprintf(ofp, "%g mx %g my l\n",0.,1.);
		fprintf(ofp, "closepath\n");
		fprintf(ofp, "stroke\n");
		fprintf(ofp, "%g mx %g my m (%g) show\n",-0.05,0.,condmin);
		fprintf(ofp, "%g mx %g my m (%g) show\n",-0.05,1.,condmax);
		fprintf(ofp, "%g mx %g my m (%g) show\n",0.,-0.03,pmin);
		fprintf(ofp, "%g mx %g my m (%g) show\n",1.,-0.03,pmax);

		fprintf(ofp,"0 0 0 c ");//black
		fprintf(ofp, "%g mx %g my m (condup) show\n", 0.05,0.01);
		fprintf(ofp,"1 0 0 c ");//red
		fprintf(ofp, "%g mx %g my m (conddown) show\n", 0.2,0.01);

		for(iseg=1; iseg<=nseg; iseg++) if(segtyp[iseg] >= 3 && segtyp[iseg] <= 5){
			xs = (segpress[iseg] - pmin)/(pmax - pmin);
			fprintf(ofp,"0 0 0 c ");//black
			ys = (condup[iseg] - condmin)/(condmax - condmin);
			fprintf(ofp, "%g mx %g my m ", xs + 0.004,ys);
			fprintf(ofp, "%g mx %g my 0.004 %g mul 0 360 arc fill\n", xs,ys,scalefacx);
			fprintf(ofp,"1 0 0 c ");//red
			ys = (conddown[iseg] - condmin)/(condmax - condmin);
			fprintf(ofp, "%g mx %g my m ", xs + 0.004,ys);
			fprintf(ofp, "%g mx %g my 0.004 %g mul 0 360 arc fill\n", xs,ys,scalefacx);
		}
	}
	fprintf(ofp, "showpage\n");
	fclose(ofp);
}
