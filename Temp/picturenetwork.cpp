/**********************************************************
picturenetwork.cpp - project network on z = 0 plane.  TWS Dec. 07
Uses parameters from CountourParams.dat
Labels nodes with nodevar and segments with segvar (must be float).
Generates a postscript file.
Added visualization of tissue points, April 2010.
Updated for 3D version, October 2017.
***********************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void picturenetwork(float *nodvar, float *segvar, const char fname[])
{
	extern int max,mxx,myy,mzz,nseg,nnod,is2d;
	extern int *segtyp,*nodtyp,*ista,*iend;
	extern int ***nbou;//added April 2010
	extern float *axt,*ayt,*azt;//added April 2010
	extern float *diam,**cnode,**xsl0,**xsl1,**xsl2;
	int i,j,k,iseg,inod,showpoint,ncontplne=2;
	float xmin,xmax,ymin,ymax,xs,ys,picfac;
	
	FILE *ofp;

	xmin = 0.;
	xmax = sqrt(SQR(xsl1[1][ncontplne]-xsl0[1][ncontplne]) + SQR(xsl1[2][ncontplne]-xsl0[2][ncontplne]) + SQR(xsl1[3][ncontplne]-xsl0[3][ncontplne]));
	ymin = 0.;
	ymax = sqrt(SQR(xsl2[1][ncontplne]-xsl0[1][ncontplne]) + SQR(xsl2[2][ncontplne]-xsl0[2][ncontplne]) + SQR(xsl2[3][ncontplne]-xsl0[3][ncontplne]));

	picfac = FMIN(500./xmax,700./ymax);//updated April 2010
	ofp = fopen(fname, "w");
	fprintf(ofp, "%%!PS-Adobe-2.0\n");
	fprintf(ofp, "%%%%Pages: 1\n");
	fprintf(ofp, "%%%%EndComments\n");
	fprintf(ofp, "%%%%Page: 1 1\n");
	fprintf(ofp, "/mx {%g mul 50 add} def\n",picfac);
	fprintf(ofp, "/my {%g mul 50 add} def\n",picfac);//updated April 2010
	fprintf(ofp, "/m {moveto} def\n");
	fprintf(ofp, "/l {lineto} def\n");
	fprintf(ofp, "/s {show} def\n");
	fprintf(ofp, "/sl {setlinewidth} def\n");
	fprintf(ofp, "/sc {setrgbcolor} def\n");
	fprintf(ofp, "newpath\n");
	fprintf(ofp, "%g mx %g my m\n",xmin,ymin);
	fprintf(ofp, "%g mx %g my l\n",xmax,ymin);
	fprintf(ofp, "%g mx %g my l\n",xmax,ymax);
	fprintf(ofp, "%g mx %g my l\n",xmin,ymax);
	fprintf(ofp, "closepath\n");
	fprintf(ofp, "stroke\n");
//show tissue points
	fprintf(ofp, "/Times-Roman findfont\n");
	fprintf(ofp, "8 scalefont\n");
	fprintf(ofp, "setfont\n");
	fprintf(ofp,"0 0 0 sc\n");//black
	if(is2d){	//simple 2D version
		for(i=1; i<=mxx; i++) for(j=1; j<=myy; j++){
			showpoint = 0;
			for(k=1; k<=mzz; k++) if(nbou[i][j][k] > 0) showpoint = 1;
			if(showpoint == 1) fprintf(ofp, "%g mx %g my m (.) s\n", axt[i],ayt[j]);
		}
	}
	else{	//full 3D version
		for(i=1; i<=mxx; i++) for(j=1; j<=myy; j++) for(k=1; k<=mzz; k++) if(nbou[i][j][k] > 0){
			xs = 0.;
			ys = 0.;
			xs += (axt[i]-xsl0[1][ncontplne])*(xsl1[1][ncontplne]-xsl0[1][ncontplne])/xmax;
			ys += (axt[i]-xsl0[1][ncontplne])*(xsl2[1][ncontplne]-xsl0[1][ncontplne])/ymax;
			xs += (ayt[j]-xsl0[2][ncontplne])*(xsl1[2][ncontplne]-xsl0[2][ncontplne])/xmax;
			ys += (ayt[j]-xsl0[2][ncontplne])*(xsl2[2][ncontplne]-xsl0[2][ncontplne])/ymax;
			xs += (azt[k]-xsl0[3][ncontplne])*(xsl1[3][ncontplne]-xsl0[3][ncontplne])/xmax;
			ys += (azt[k]-xsl0[3][ncontplne])*(xsl2[3][ncontplne]-xsl0[3][ncontplne])/ymax;
			fprintf(ofp, "%g mx %g my m (.) s\n", xs,ys);
		}
	}
	fprintf(ofp, "/Times-Roman findfont\n");
	fprintf(ofp, "4 scalefont\n");
	fprintf(ofp, "setfont\n");
//plot vessels
	//fprintf(ofp,"1 0 0 sc\n");//red
	for(iseg=1; iseg<=nseg; iseg++) if(segtyp[iseg] < 10){
//		fprintf(ofp,"%g sc\n",picfac*diam[iseg]);
		if(segtyp[iseg] == 0) fprintf(ofp,"0.5 0 1 sc\n");//dark purple
		else if(segtyp[iseg] == 1) fprintf(ofp,"1 0 1 sc\n");//purple
		else if(segtyp[iseg] == 3) fprintf(ofp,"0 0 1 sc\n");//blue
		else if(segtyp[iseg] == 4) fprintf(ofp,"0 1 0 sc\n");//green
		else if(segtyp[iseg] == 5) fprintf(ofp,"1 0 0 sc\n");//red
		else fprintf(ofp,"0 0 0 sc\n");//black
		xs = 0.;
		ys = 0.;
		for(i=1; i<=3; i++){
			xs += (cnode[i][ista[iseg]]-xsl0[i][ncontplne])*(xsl1[i][ncontplne]-xsl0[i][ncontplne])/xmax;
			ys += (cnode[i][ista[iseg]]-xsl0[i][ncontplne])*(xsl2[i][ncontplne]-xsl0[i][ncontplne])/ymax;
		}
		fprintf(ofp, "%g mx %g my m ", xs,ys);
		xs = 0.;
		ys = 0.;
		for(i=1; i<=3; i++){
			xs += (cnode[i][iend[iseg]]-xsl0[i][ncontplne])*(xsl1[i][ncontplne]-xsl0[i][ncontplne])/xmax;
			ys += (cnode[i][iend[iseg]]-xsl0[i][ncontplne])*(xsl2[i][ncontplne]-xsl0[i][ncontplne])/ymax;
		}

		fprintf(ofp, "%g mx %g my l ", xs,ys);
		fprintf(ofp, "stroke\n");
	}
//label nodes in black
	fprintf(ofp,"0 0 0 sc\n");//black
	for(inod=1; inod<=nnod; inod++) if(nodtyp[inod] > 0){
		xs = 0.;
		ys = 0.;
		for(i=1; i<=3; i++){
			xs += (cnode[i][inod]-xsl0[i][ncontplne])*(xsl1[i][ncontplne]-xsl0[i][ncontplne])/xmax;
			ys += (cnode[i][inod]-xsl0[i][ncontplne])*(xsl2[i][ncontplne]-xsl0[i][ncontplne])/ymax;
		}
		fprintf(ofp, "%g mx %g my m ", xs + 0.5/picfac,ys);
		fprintf(ofp, "(%g) show\n",nodvar[inod]);
	}
//label segments in blue
	fprintf(ofp,"0 0 1 sc\n");//blue
	for(iseg=1; iseg<=nseg; iseg++) if(segtyp[iseg] < 10){
		xs = 0.;
		ys = 0.;
		for(i=1; i<=3; i++){
			xs += ((cnode[i][ista[iseg]]+cnode[i][iend[iseg]])/2.-xsl0[i][ncontplne])*(xsl1[i][ncontplne]-xsl0[i][ncontplne])/xmax;
			ys += ((cnode[i][ista[iseg]]+cnode[i][iend[iseg]])/2.-xsl0[i][ncontplne])*(xsl2[i][ncontplne]-xsl0[i][ncontplne])/ymax;
		}
		fprintf(ofp, "%g mx %g my m ", xs + 0.5*picfac,ys);
		fprintf(ofp, "(%g) show\n",segvar[iseg]);
	}
	fprintf(ofp, "showpage\n");
	fclose(ofp);
}
