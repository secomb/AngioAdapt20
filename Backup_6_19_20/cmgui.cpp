/************************************************************************
cmgui - for Greens.  TWS December 2011
Based on code provided by Gib Bogle
Produces greens.exelem and greens.exnode
To visualize, use CMGUI (e.g. cmgui_install.msi) at:
http://sourceforge.net/projects/cmiss/files/cmgui/cmgui-wx-2.8.0/
Start CMGUI, File > Open > Com File
Navigate to the directory containing the data files and select greens.com.txt (as listed below).
At the bottom of the window that pops up, click the All button.  Change the view with the Layout,
Up and Front selections on the left hand side of the Graphics window.  Rotate by holding the left
button down and moving the mouse, pan with the middle button and zoom with the right button.
Use Graphics > Scene editor to change the display.  Note also the Spectrum editor.
calling code:

Modified Feb. 2017 for use with AngioAdapt
Note that a single elements file has to be created.
To allow smooth network growth, the node list has to be updated one time step before each step appears.
++++++++++++++++++greens.com.txt+++++++++++++++++++++++++++++++++
# Create a material in addition to the default.
gfx cre mat gold ambient 1 0.7 0 diffuse 1 0.7 0 specular 0.5 0.5 0.5 shininess 0.8

gfx create spectrum jet
gfx modify spectrum jet clear overwrite_colour
gfx modify spectrum jet linear range 0 1 red   colour_range 0 1 ambient diffuse component 1 
gfx modify spectrum jet linear range 0 1 green colour_range 0 1 ambient diffuse component 2
gfx modify spectrum jet linear range 0 1 blue  colour_range 0 1 ambient diffuse component 3

# Read in the reticular mesh (group vessels) and hide the axes.
gfx read nodes network001.exnode
gfx read elements network001.exelem

# The radius of the vessel is stored in component 1 of field
# 'vessel_radius', defined over the elements in the vessels group.

# Destroy the default lines.
gfx modify g_element vessels lines delete

gfx destroy node all
gfx modify g_element vessels general clear;
gfx modify g_element vessels cylinders coordinate coordinates tessellation default local circle_discretization 12 radius_scalar vessel_radius scale_factor 1 native_discretization NONE data node_colour spectrum jet
gfx modify g_element vessels node_points coordinate coordinates local glyph sphere general size "0*0*0" centre 0,0,0 font default orientation vessel_radius scale_factors "2*2*2" data node_colour spectrum jet

# Open the graphics window and turn on perspective (if desired).
gfx cre win 1
gfx mod win 1 view perspective
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <cstdio>
#include <math.h>
#include <string.h>
#include "nrutil.h"


void WriteExnodeHeader(FILE *exnode) // Write initial section of .exnode file
{
//    fprintf(exnode, "Region: /vessels\n");
    fprintf(exnode, "Group name: vessels\n");
    fprintf(exnode, " #Fields=3\n");
    fprintf(exnode, " 1) coordinates, coordinate, rectangular cartesian, #Components=3\n");
    fprintf(exnode, "  x.  Value index=1, #Derivatives=0, #Versions=1\n");
    fprintf(exnode, "  y.  Value index=2, #Derivatives=0, #Versions=1\n");
    fprintf(exnode, "  z.  Value index=3, #Derivatives=0, #Versions=1\n");
    fprintf(exnode, " 2) vessel_radius, coordinate, rectangular cartesian, #Components=1\n");
    fprintf(exnode, "  1.  Value index=4, #Derivatives=0, #Versions=1\n");
	fprintf(exnode, " 3) node_colour, coordinate, rectangular cartesian, #Components=3\n");
	fprintf(exnode, "  1.  Value index=5, #Derivatives=0, #Versions=1\n");
	fprintf(exnode, "  2.  Value index=6, #Derivatives=0, #Versions=1\n");
	fprintf(exnode, "  3.  Value index=7, #Derivatives=0, #Versions=1\n");
}

void WriteExelemHeader(FILE *exelem)  // Write initial section of .exelem file
{
 //   fprintf(exelem, "Region: /vessels\n");
    fprintf(exelem, "Group name: vessels\n");
    fprintf(exelem, " Shape.  Dimension=1\n");
    fprintf(exelem, " #Scale factor sets= 1\n");
    fprintf(exelem, "  l.Lagrange, #Scale factors= 2\n");
    fprintf(exelem, " #Nodes= 2\n #Fields=3\n");
    fprintf(exelem, " 1) coordinates, coordinate, rectangular cartesian, #Components=3\n");
    fprintf(exelem, "   x.  l.Lagrange, no modify, standard node based.\n     #Nodes= 2\n");
    fprintf(exelem, "      1.  #Values=1\n");
    fprintf(exelem, "       Value indices:     1\n       Scale factor indices:   1\n");
    fprintf(exelem, "      2.  #Values=1\n");
    fprintf(exelem, "       Value indices:     1\n       Scale factor indices:   2\n");
    fprintf(exelem, "   y.  l.Lagrange, no modify, standard node based.\n     #Nodes= 2\n");
    fprintf(exelem, "      1.  #Values=1\n");
    fprintf(exelem, "       Value indices:     1\n       Scale factor indices:   1\n");
    fprintf(exelem, "      2.  #Values=1\n");
    fprintf(exelem, "       Value indices:     1\n       Scale factor indices:   2\n");
    fprintf(exelem, "   z.  l.Lagrange, no modify, standard node based.\n     #Nodes= 2\n");
    fprintf(exelem, "      1.  #Values=1\n");
    fprintf(exelem, "       Value indices:     1\n       Scale factor indices:   1\n");
    fprintf(exelem, "      2.  #Values=1\n");
    fprintf(exelem, "       Value indices:     1\n       Scale factor indices:   2\n");
    fprintf(exelem, " 2) vessel_radius, coordinate, rectangular cartesian, #Components=1\n");
    fprintf(exelem, "   1.  l.Lagrange, no modify, standard node based.\n     #Nodes= 2\n");
    fprintf(exelem, "      1.  #Values=1\n");
    fprintf(exelem, "       Value indices:     1\n       Scale factor indices:   1\n");
    fprintf(exelem, "      2.  #Values=1\n");
    fprintf(exelem, "       Value indices:     1\n       Scale factor indices:   2\n");
    fprintf(exelem, " 3) node_colour, coordinate, rectangular cartesian, #Components=3\n");
    fprintf(exelem, "   1.  l.Lagrange, no modify, standard node based.\n");
    fprintf(exelem, "     #Nodes= 2\n");
    fprintf(exelem, "      1.  #Values=1\n");
    fprintf(exelem, "       Value indices:     1\n");
    fprintf(exelem, "       Scale factor indices:   1\n");
    fprintf(exelem, "      2.  #Values=1\n");
    fprintf(exelem, "       Value indices:     1\n");
    fprintf(exelem, "       Scale factor indices:   2\n");
    fprintf(exelem, "   2.  l.Lagrange, no modify, standard node based.\n");
    fprintf(exelem, "     #Nodes= 2\n");
    fprintf(exelem, "      1.  #Values=1\n");
    fprintf(exelem, "       Value indices:     1\n");
    fprintf(exelem, "       Scale factor indices:   1\n");
    fprintf(exelem, "      2.  #Values=1\n");
    fprintf(exelem, "       Value indices:     1\n");
    fprintf(exelem, "       Scale factor indices:   2\n");   
    fprintf(exelem, "   3.  l.Lagrange, no modify, standard node based.\n");
    fprintf(exelem, "     #Nodes= 2\n");
    fprintf(exelem, "      1.  #Values=1\n");
    fprintf(exelem, "       Value indices:     1\n");
    fprintf(exelem, "       Scale factor indices:   1\n");
    fprintf(exelem, "      2.  #Values=1\n");
    fprintf(exelem, "       Value indices:     1\n");
    fprintf(exelem, "       Scale factor indices:   2\n");
}

void cmgui(int imain)
{
	extern int max,mxx,myy,mzz,nseg,nnod,nsegprev,imainmax,hexperiodic;
	extern int *segtyp,*ista,*iend,*boundseg;
	extern float aphexp,*diam,*x,*y,**cnode,**pvseg,**hex_norm;				

	int iseg,i;
	float red,green,blue,xz,xzmin,xzmax;
	char fname[80];
	FILE *exelem,*exnode,*exnodeprev;

	if(imain == imainmax){	//cmgui wants a single element file for all time points
		exelem = fopen("Current\\network.exelem","w");
		WriteExelemHeader(exelem);
		for(iseg=1; iseg<=nseg; iseg++){
			fprintf(exelem, "Element: %d 0 0\n", iseg);
			fprintf(exelem, "  Nodes: %d %d\n", 2*iseg-1,2*iseg);
			fprintf(exelem, "  Scale factors: 1 1\n");
		}
		fclose(exelem);
	}

	//write current node file
	sprintf(fname,"Current\\network%03i.exnode",imain);
	exnode = fopen(fname,"w");
	WriteExnodeHeader(exnode);
	xzmin = 0.;	//set up fixed range of oxygen levels
	xzmax = 100.;
	//xzmin = 1.e6;	//set up range of oxygen levels
	//xzmax = -1.e6;
	//for(iseg=1; iseg<=nseg; iseg++) if(segtyp[iseg] >= 3 && segtyp[iseg] <= 5){
	//	xzmin = FMIN(xzmin,pvseg[iseg][1]);
	//	xzmax = FMAX(xzmax,pvseg[iseg][1]);
	//}
	for(iseg=1; iseg<=nseg; iseg++){
		for(i=1; i<=3; i++){
			x[i] = cnode[i][ista[iseg]];
			y[i] = cnode[i][iend[iseg]];
			if(hexperiodic == 1 && boundseg[iseg] > 0) y[i] += 2.*aphexp*hex_norm[boundseg[iseg]][i];	//segments that cross boundary
		}
		if(segtyp[iseg] != 10){
			if(segtyp[iseg] >= 3 && segtyp[iseg] <= 5){	//Set up colors using Matlab 'jet' scheme
				if(xzmin != xzmax) xz = (pvseg[iseg][1] - xzmin)/(xzmax - xzmin);
				else xz = 0.75;
				red  = FMIN(FMAX(1.5-4.*fabs(xz-0.75), 0.), 1.);
				green= FMIN(FMAX(1.5-4.*fabs(xz-0.5), 0.), 1.);
				blue = FMIN(FMAX(1.5-4.*fabs(xz-0.25), 0.), 1.);
			}
			else if(segtyp[iseg] == 1){		//purple, sprout
				red  = 1.;
				green= 0.;
				blue = 1.;
			}
			else if(segtyp[iseg] == 0){		//redder purple, tip of sprout
				red  = 0.5;
				green= 0.;
				blue = 1.;
			}
			else{							//black, error
				red  = 0.;
				green= 0.;
				blue = 0.;
			}
	//  write to nodes file
			fprintf(exnode, "Node: %d\n", 2*iseg-1);
			fprintf(exnode, "%6.1f %6.1f %6.1f\n", x[1],x[2],x[3]);
			fprintf(exnode, "%6.2f\n", diam[iseg]/2.);
			fprintf(exnode, "%6.2f %6.2f %6.2f\n", red,green,blue);
			fprintf(exnode, "Node: %d\n", 2*iseg);
			fprintf(exnode, "%6.1f %6.1f %6.1f\n", y[1],y[2],y[3]);
			fprintf(exnode, "%6.2f\n", diam[iseg]/2.);
			fprintf(exnode, "%6.2f %6.2f %6.2f\n", red,green,blue);
		}
		else{	//these are dummy nodes and segments but nodes need to be placed to avoid flying fragments
	//  write to nodes file
			fprintf(exnode, "Node: %d\n", 2*iseg-1);
			fprintf(exnode, "%6.1f %6.1f %6.1f\n", x[1],x[2],x[3]);
			fprintf(exnode, "0.\n");
			fprintf(exnode, "0. 0. 0.\n");
			fprintf(exnode, "Node: %d\n", 2*iseg);
			fprintf(exnode, "%6.1f %6.1f %6.1f\n", y[1],y[2],y[3]);
			fprintf(exnode, "0.\n");
			fprintf(exnode, "0. 0. 0.\n");
		}
	}
	fclose(exnode);
	//append additional data to previous node file
	if(imain > 0){
		sprintf(fname,"Current\\network%03i.exnode",imain-1);
		exnodeprev = fopen(fname,"a");
		for(iseg=nsegprev+1; iseg<=nseg; iseg++){	//put the nodes at the correct starting points to avoid flying fragments
			for(i=1; i<=3; i++){
				x[i] = cnode[i][ista[iseg]];
				y[i] = cnode[i][iend[iseg]];
				if(hexperiodic == 1 && boundseg[iseg] > 0) y[i] += 2.*aphexp*hex_norm[boundseg[iseg]][i];	//segments that cross boundary
			}
			fprintf(exnodeprev, "Node: %d\n", 2*iseg-1);
			fprintf(exnodeprev, "%6.1f %6.1f %6.1f\n", x[1],x[2],x[3]);
			fprintf(exnodeprev, "0.\n");
			fprintf(exnodeprev, "0. 0. 0.\n");
			fprintf(exnodeprev, "Node: %d\n", 2*iseg);
			fprintf(exnodeprev, "%6.1f %6.1f %6.1f\n", y[1],y[2],y[3]);
			fprintf(exnodeprev, "0.\n");
			fprintf(exnodeprev, "0. 0. 0.\n");
		}
		fclose(exnodeprev);
	}
	nsegprev = nseg;
}
