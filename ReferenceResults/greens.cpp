/************************************************************************
Green's function approach for multiple reacting species.
T.W. Secomb July 2007 - based on Greens.f by R. Hsu.
See http://www.physiology.arizona.edu/people/secomb/greens.html

Variable array sizes and bounds, using Numerical Recipes utilities.
Tissue-vessel and tissue-tissue matrices computed on the fly to save memory.
No nondimensionalization.  Lengths, diameters in microns, times in s.
Flows in nanoliters/min
Oxygen concentrations in cm^3/cm^3
Consumption rates in cm^3/cm^3/s - changed in V2
Partial pressures in mmHg
Diffusivities in cm2/s, converted to micron^2/s

Special parameters for oxygen
p50, fn --- parameters in the Hill equation
cs --- red blood cell oxygen binding capacity in cm^3 O2/cm^3
alphab --- average solubility in blood in cm^3 O2/cm^3/mmHg
gamma1 --- intravascular resistance, varies with vessel diameter, in mmHg.cm.s/cm^3 O2

Main variables:
  gvv --- Green's function matrix for vessels
  mat --- matrix for vessel strengths
  al --- matrix giving dependence of vessel convective fluxes on source strengths
  lseg --- segment length
  ds --- subsegment length
  qv --- oxygen efflux from subsegment
  pv --- PO2 in the subsegment
  cv --- oxygen concentration in the subsegment
  qt --- tissue source strength
  pvt --- PO2 on vessels due to source terms in tissue
  ptv --- PO2 in tissue due to source terms on vessels
  q --- flow rate, qq = abs(q)

Version 2.0, May 2010.
With 9/08 updates.  New vessel-vesel interaction coefficient. January 2009
With alternate terms for 2D version.  May 2009
With g0 computed as part of linear system, for permeable solutes.  March 2010
  g0method = 1:  include g0 in linear system to be solved - fastest method *****
  g0method = 2:  theoretical estimates of dqsum/dg0 - was used in Version 1
For impermeable solutes, method 2 is always used
With choice of Gauss-Jordan, LU or biconjugate gradient linear solvers.  March 2010
  linmethod = 1:  Gaussian elimination - was used in Version 1
  linmethod = 2:  LU decomposition
  linmethod = 3:  biconjugate gradient (iterative) - fastest method *****
Does not require that species 1 is oxygen. Allows for non-diffusing solutes.  April 2010
Creates log file. April 2010.
During tissue loop, scales vessel sources so that qvsum = qtsum, for faster convergence.  April 2010.
Modified for compatibility with Mac XCode compiler.  April 2010.
Includes intravascular resistance for all solutes.  May 2010.
Includes non-diffusible solutes.  May 2010.
Version 3.0, May 17, 2011.  Revised for angioadapt, January 2012
Uses convect instead of genalpha and genalphahd.
This gives improved convergence if hematocrit is non-uniform in network
Version 3.1, July 22, 2011 - includes 'greensverbose' option
**************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "nrutil.h"
#include <cuda_runtime.h>

void putrank(void);
void initgreens();
void blood(float c, float hem, float* p, float* pp);
float bloodconc(float p, float h);
void tissrate(int nsp, float* p, float* mtiss, float* mptiss);
void convect(int isp);		//new subroutine, August 2010, replaces genalpha and genalphahd
//void testconvect(int isp);
float* eval(int slsegdiv, float req, float* x);
float length(float* x);
float length0(float* x);

void gaussj(double** a, int n, double** b, int m);
void ludcmp(double** a, int n, int* indx, double* d);
void lubksb(double** a, int n, int* indx, double* b);
double bicgstab(double** a, double* b, double* x, int n, double eps, int itmax);

double bicgstabBLASD(double** a, double* b, double* x, int n, double eps, int itmax);
void bicgstabBLASDinit(int nnnvGPU);
void bicgstabBLASDend(int nnnvGPU);
void tissueGPUinit(int nntGPU, int nnvGPU);
void tissueGPUend(int nntGPU, int nnvGPU);
void tissueGPU1c(int stepGPU, int blockstep, int iblockstep);
void tissueGPU2c(int stepGPU);
void tissueGPU3c(int stepGPU);
void tissueGPUcopy();

void greens(void)
{
	extern int nmaxvessel, nmaxtissue, nmax, g0method, linmethod, imain, hexperiodic;
	extern int mxx, myy, mzz, nnt, nnv, nseg, nsp, nnodbc;
	extern int is2d; //needed for 2d version
	extern int nntGPU, nnvGPU, useGPU;//needed for GPU version
	extern int* mainseg, ** tisspoints, * permsolute, * segtyp, * ista, * iend;
	extern int* segname, * nspoint, * istart, * nodout, * bcnod, * lowflow;//not bcnodname.  April 2010
	extern int* errvesselcount, * errtissuecount, * boundseg;
	extern int* imaxerrvessel, * imaxerrtissue, * nresis;  //added April 2010
	extern int* oxygen, * diffsolute, * indx; //added April 2010
	extern int** nodseg;
	extern int*** nbou, *** nbou_old;

	extern float* dtmin;//added July 2011
	extern float* sourcefac, * sourcefacGF;	//added Feb 2020
	extern float p50, cs, req, q0, fac, flowfac, lowflowcrit, errfac/*,pcr,m0*/;
	extern float v, vol, vdom, tlength, pi1, aphexp;
	extern float alx, aly, alz, w2d, r2d; //needed for 2d version
	extern float* axt, * ayt, * azt, * ds, * diff, * pmin, * pmax, * pmean, * g0, * g0fac, * g0facnew, * pref;
	extern float* diam, * rseg, * q, * qq, * hd, * bchd, * qvtemp, * qvfac;
	extern float* x, * y, * lseg, * mtiss, * mptiss, * dqvsumdg0, * dqtsumdg0;
	extern float* epsvessel, * epstissue, * eps, * errvessel, * errtissue, * p;
	extern float* rhs, * rhstest, * g0old, * ptt, * ptpt, * qtsum, * qvsum;

	extern float** tissparam, ** cnode, ** scos, ** ax, ** bcp;
	extern float** qv, ** qt, ** pv, ** pev, ** pt, ** resisdiam, ** resis;
	extern float** qvseg, ** pvseg, ** pevseg, ** pvt, ** pvprev, ** qvprev, ** cv, ** dcdp;
	extern float** ptprev, ** ptv, ** gamma1, ** cv0, ** conv0, ** gvv, ** al, ** hex_norm;
	extern double** mat, ** rhsg, * rhsl, * matx;
	extern float*** dtt;
	extern float* pt000, * qt000, * pv000, * qv000, * dtt000;	//needed for GPU version
	int stepGPU = 8;

	int i, j, k, ix, iy, iz, jx, jy, jz, iseg, nt, ineg, ihigh, isp, imaxerr, ii;
	int ixdiff, iydiff, izdiff, isn, jseg, kmain, ktissue, kvessel, itp, jtp, convflag, convflagt, convflagv;
	int greensverbose = 0;//Added July 22 2011
	int bicgstabit = 2000;//5000; //2000 parameter for biconjugate gradient method.  March 2010 - increase for larger networks
	int blockstep = 4, iblockstep;		//use for tissue updating by blocks
	int ndomain, idomain;

	float x11, x22, x33, duration, rsegmax, dsmax, gvarfac;
	float gtt, gtv, gvt, disp2, ds2, dist, d, de2, dtave, den, dqsumdg0;
	float dif, err, qhdcrit;
	float lam, lam3d, lam2d, req2, r2d2;//modified 9/09
	double dd, bicgstaberr = 0.0001; //parameter for biconjugate gradient method.  March 2010

	FILE* ofp, * ofp1, * ofp2;	//General file, MasterLog file, GreensLog file
	clock_t tstart, tfinish, tstart1, tfinish1;

	req2 = SQR(req);
	w2d = alz;  //needed for 2d
	r2d2 = (SQR(alx) + SQR(aly) + SQR(alz));
	r2d = sqrt(r2d2);
	lam3d = 0.09;	//Underrelax iteration of tissue levels
	lam2d = 0.07;//0.1;	
	if (useGPU) {
		lam3d = 0.005;//0.0065;
		lam2d = 0.005;// best old value 0.00175 JA 5_2014; 0.01;	0.003
	}
	if (useGPU && nnv + 1 > nnvGPU) {
		bicgstabBLASDend(nnvGPU);
		tissueGPUend(nntGPU, nnvGPU);	//these have to be done also since bigstabBLASend interferes with TissueGPU
		nnvGPU = nnv + 100;	//increment by 100 to reduce number of times that this has to be done
		printf("*** more memory allocated on GPU ***\n");
		tissueGPUinit(nntGPU, nnvGPU);
		bicgstabBLASDinit(nnvGPU);
	}

	//setup mainseg (must be done after setuparrays2)
	for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] >= 3 && segtyp[iseg] <= 5)
		for (i = 0; i < nspoint[iseg]; i++) mainseg[istart[iseg] + i] = iseg;
	//identify vessel points
	for (i = 1; i <= nnv; i++) {
		iseg = mainseg[i];
		isn = i - istart[iseg];
		for (j = 1; j <= 3; j++) ax[j][i] = cnode[j][ista[iseg]] + scos[j][iseg] * ds[iseg] * (isn + 0.5);
	}
	//index tissue points
	for (i = 1; i <= mxx; i++) for (j = 1; j <= myy; j++) for (k = 1; k <= mzz; k++) {
		nt = nbou[i][j][k];
		if (nt > 0) {
			tisspoints[1][nt] = i;
			tisspoints[2][nt] = j;
			tisspoints[3][nt] = k;
			sourcefac[nt] = 1.;
			sourcefacGF[nt] = 1.;
			//if ((k - 0.5) * 20 <= mzz) sourcefac[nt] = 0.;	//turn off oxygen consumption near bottom surface
			//if ((mzz - k + 0.5) * 20 <= mzz) sourcefac[nt] = 0.;	//turn off oxygen consumption near top surface
			//if ((k - 0.5) * 20 <= mzz) sourcefacGF[nt] = 0.;	//turn off GF production near bottom surface
		}
	}
	//calculate the distance of tissue points to the nearest vessel
	dtave = 0.;
	for (itp = 1; itp <= nnt; itp++) {
		i = tisspoints[1][itp];
		j = tisspoints[2][itp];
		k = tisspoints[3][itp];
		dtmin[itp] = 1.e6;
		for (jseg = 1; jseg <= nseg; jseg++) if (segtyp[jseg] >= 3 && segtyp[jseg] <= 5) {
			for (ii = 1; ii <= 3; ii++) {
				x[ii] = cnode[ii][ista[jseg]];
				y[ii] = cnode[ii][iend[jseg]];
				if (hexperiodic == 1 && boundseg[jseg] > 0) y[ii] += 2. * aphexp * hex_norm[boundseg[jseg]][ii];	//segments that cross boundary
			}
			x11 = (axt[i] - x[1]) * scos[2][jseg] - (ayt[j] - x[2]) * scos[1][jseg];
			x22 = (ayt[j] - x[2]) * scos[3][jseg] - (azt[k] - x[3]) * scos[2][jseg];
			x33 = (azt[k] - x[3]) * scos[1][jseg] - (axt[i] - x[1]) * scos[3][jseg];
			disp2 = SQR(x11) + SQR(x22) + SQR(x33);
			ds2 = SQR(axt[i] - x[1]) + SQR(ayt[j] - x[2]) + SQR(azt[k] - x[3]);
			de2 = SQR(axt[i] - y[1]) + SQR(ayt[j] - y[2]) + SQR(azt[k] - y[3]);
			if (FMAX(ds2, de2) - disp2 > SQR(lseg[jseg])) d = FMAX(sqrt(FMIN(ds2, de2)) - rseg[jseg], 0.);
			else d = FMAX(sqrt(disp2) - rseg[jseg], 0.);
			if (d < dtmin[itp]) dtmin[itp] = d;
		}
		dtave += dtmin[itp];
	}
	dtave = dtave / nnt;
	vdom = nnt * vol;
	den = sqrt(vdom / tlength);
	if (greensverbose) printf("Average distance from tissue node to the nearest vessel = %f\n", dtave);
	if (greensverbose) printf("Sqrt(Tissue Volume/vessel length) = %f\n", den);
	//Calculate intravascular or wall transport resistance.  Zero unless specified in IntravascRes.dat.
	//If not oxygen, assume value from data is 1/(wall permeability in um/s)
	for (isp = 1; isp <= nsp; isp++) for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] >= 3 && segtyp[iseg] <= 5) {
		gamma1[iseg][isp] = 0.;
		if (nresis[isp] != 0) {
			gamma1[iseg][isp] = resis[1][isp];
			for (j = 2; j <= nresis[isp]; j++) if (diam[iseg] <= resisdiam[j][isp] && diam[iseg] > resisdiam[j - 1][isp])
				gamma1[iseg][isp] = resis[j - 1][isp] + (resis[j][isp] - resis[j - 1][isp])
				* (diam[iseg] - resisdiam[j - 1][isp]) / (resisdiam[j][isp] - resisdiam[j - 1][isp]);
			if (diam[iseg] > resisdiam[nresis[isp]][isp]) gamma1[iseg][isp] = resis[nresis[isp]][isp];
			if (oxygen[isp] != 1) gamma1[iseg][isp] = gamma1[iseg][isp] / 2. / pi1 / diam[isp];
		}
	}
	//vessel ~ vessel matrix elements gvv
	//Uses empirical fit to results from elliptical integral form for diagonal elements, updated 2009
	//if center of one segment lies within the other segment, calculate Gvv as for self-interaction term
	//based on larger radius and larger length. (Jan. 08)
	for (i = 1; i <= nnv; i++) {
		iseg = mainseg[i];
		for (j = 1; j <= nnv; j++) {
			jseg = mainseg[j];
			dsmax = FMAX(ds[iseg], ds[jseg]);
			rsegmax = FMAX(rseg[iseg], rseg[jseg]);
			for (k = 1; k <= 3; k++)	x[k] = ax[k][j] - ax[k][i];
			if (hexperiodic) ndomain = 7;
			else ndomain = 1;
			gvv[i][j] = 0.;
			for (idomain = 0; idomain < ndomain; idomain++) {
				for (k = 1; k <= 3; k++) {
					y[k] = x[k];
					if (idomain > 0) y[k] += 2. * aphexp * hex_norm[idomain][k];
				}
				dist = length0(y);
				if (dist < rsegmax) {
					gvarfac = 0.6 * exp(-0.45 * dsmax / rsegmax);
					//for distinct vessels close together, make distance rsegmax in following calculation, to improve convergence. TWS2011
					if (iseg != jseg) dist = rsegmax;
					gvv[i][j] += (1.298 / (1. + 0.297 * powf(dsmax / rsegmax, 0.838)) - gvarfac * SQR(dist / rsegmax)) * fac / rsegmax;
					//for 2D version, additional terms give effect of boundaries (reflection method)
					if (is2d) {
						gvv[i][j] -= fac * 2. / w2d * 0.926 * SQR(1. - 1. / (1. + 0.36 * dsmax / w2d)) * powf(1. + dsmax / w2d, 0.27);
						gvv[i][j] += fac * 2. / w2d * (log(r2d / w2d + 0.5 + 0.27 / r2d * w2d) - 0.117);
					}
				}
				else {
					if (is2d) gvv[i][j] += 2. * fac / w2d * log(r2d / dist);
					else gvv[i][j] += fac / dist;
				}
			}
		}
	}

	// tissue ~ vessel, vessel ~ tissue: compute matrix elements gvt, gtv on fly as needed
	// tissue ~ tissue: construct matrix of distances from a corner node
	//diagonal elements of tissue ~ tissue matrix gtt

	if (is2d) gtt = fac / w2d * (2. * log(r2d / req) + 0.5);
	else gtt = 1.2 * fac / req;
	for (jx = 1; jx <= mxx; jx++) for (jy = 1; jy <= myy; jy++) for (jz = 1; jz <= mzz; jz++) {
		x[1] = axt[1] - axt[jx];
		x[2] = ayt[1] - ayt[jy];
		x[3] = azt[1] - azt[jz];
		if (hexperiodic) ndomain = 7;
		else ndomain = 1;
		dtt[jx][jy][jz] = 0.;
		for (idomain = 0; idomain < ndomain; idomain++) {
			for (k = 1; k <= 3; k++) {
				y[k] = x[k];
				if (idomain > 0) y[k] += 2. * aphexp * hex_norm[idomain][k];
			}
			dist = length0(y);
			if (jx * jy * jz != 1 || idomain != 0) {
				if (is2d) dtt[jx][jy][jz] += 2. * fac / w2d * log(r2d / dist);
				else dtt[jx][jy][jz] += fac / dist;
			}
			else dtt[jx][jy][jz] += gtt;
		}
	}
	if (useGPU) tissueGPUcopy();	//this has to be done after dtt is set up

//detect and label vessels with very low q*hd - updated April 2010 - test if oxygen is one of the solutes
	qhdcrit = 0.;
	for (isp = 1; isp <= nsp; isp++) if (oxygen[isp] == 1) qhdcrit = lowflowcrit * tissparam[1][isp];
	for (iseg = 1; iseg <= nseg; iseg++) {
		lowflow[iseg] = 0;
		if ((qq[iseg] * (hd[iseg] + 0.01)) < qhdcrit) lowflow[iseg] = 1;//Added 0.01 to allow for high flow, zero hematocrit channels
	}
	initgreens();
	putrank();
	//	for(isp=1; isp<=nsp; isp++){//for testing purposes only
	//		convect(isp);
	//		testconvect(isp);	
	//	}
	//create log file
	ofp2 = fopen("GreensLog.txt", "w");
	fprintf(ofp2, "GreensLog.txt\n");
	fclose(ofp2);
	tstart = clock();
	//********************** start of main loop *****************************
	for (kmain = 1; kmain <= nmax; kmain++) {
		tstart1 = clock();
		if (greensverbose) printf("\n----- kmain = %i -----\n", kmain);
		else printf(" %i", kmain);
		for (isp = 1; isp <= nsp; isp++) {
			if (diffsolute[isp] == 1) {
				for (itp = 1; itp <= nnt; itp++) ptprev[itp][isp] = pt[itp][isp];
				for (i = 1; i <= nnv; i++) pvprev[i][isp] = pv[i][isp];
			}
			g0old[isp] = g0[isp];
		}
		//********************** start of vessel loop *****************************
		//compute contribution pvt from tissue source strengths qt
		if (useGPU) {
			for (isp = 1; isp <= nsp; isp++) {
				if (permsolute[isp]) {
					for (itp = 1; itp <= nnt; itp++) qt000[itp - 1] = qt[itp][isp];
					tissueGPU2c(stepGPU);	//Compute contribution from tissue source strengths
					for (i = 1; i <= nnv; i++) {
						if (is2d) pvt[i][isp] = fac * pv000[i - 1] / diff[isp] / w2d;
						else pvt[i][isp] = fac * pv000[i - 1] / diff[isp];
						if (g0method == 2) pvt[i][isp] += g0[isp];	//?????????????????????????????
					}
				}
				else if (diffsolute[isp]) for (i = 1; i <= nnv; i++) pvt[i][isp] = g0[isp];
			}
		}
		else {
			for (i = 1; i <= nnv; i++) {
				for (isp = 1; isp <= nsp; isp++) if (diffsolute[isp]) {
					if (g0method == 1 && permsolute[isp] == 1) pvt[i][isp] = 0.;
					else pvt[i][isp] = g0[isp];
				}
				for (itp = 1; itp <= nnt; itp++) {
					x[1] = ax[1][i] - axt[tisspoints[1][itp]];
					x[2] = ax[2][i] - ayt[tisspoints[2][itp]];
					x[3] = ax[3][i] - azt[tisspoints[3][itp]];
					if (hexperiodic) ndomain = 7;
					else ndomain = 1;
					gvt = 0.;
					for (idomain = 0; idomain < ndomain; idomain++) {
						for (k = 1; k <= 3; k++) {
							y[k] = x[k];
							if (idomain > 0) y[k] += 2. * aphexp * hex_norm[idomain][k];
						}
						dist = length0(y);
						if (dist <= req) {
							if (is2d) gvt += fac / w2d * (2. * log(r2d / req) + 1. - SQR(dist / req));
							else gvt += fac * (1.5 - 0.5 * SQR(dist / req)) / req;
						}
						else {
							if (is2d) gvt += 2. * fac / w2d * log(r2d / dist);
							else gvt += fac / dist;
						}
					}
					for (isp = 1; isp <= nsp; isp++) if (permsolute[isp] == 1) pvt[i][isp] += gvt / diff[isp] * qt[itp][isp];
				}
			}
		}
		//*********************************
		//compute contribution from vessel source strengths qv
		//for(i=1; i<=nnv; i++) for(j=1; j<=nnv; j++) for(isp=1; isp<=nsp; isp++) if(permsolute[isp]) pvt[i][isp] += gvv[i][j]/diff[isp]*qv[j][isp];
		//*********************************
//compute blood solute levels and PO2
		for (kvessel = 1; kvessel <= nmaxvessel; kvessel++) {
			convflagv = 1;
			for (isp = 1; isp <= nsp; isp++) {
				qvsum[isp] = 0.;
				dqvsumdg0[isp] = 0.;
				if (permsolute[isp] == 1) {
					ineg = 0;
					ihigh = 0;
					convect(isp);		//new subroutine, August 2010
					for (i = 1; i <= nnv; i++) {
						iseg = mainseg[i];
						qvprev[i][isp] = qv[i][isp];
						if (oxygen[isp] == 1) {
							if (lowflow[iseg] != 1) {//only do this if not a lowflow segment.  June 2009.
								if (cv[i][isp] < 0.) {
									ineg++;
									if (ineg == 1 && greensverbose) printf("*** Warning: cblood is negative -%i", segname[iseg]);
									if (ineg > 1 && greensverbose) printf("-%i", segname[iseg]);
								}
								if (cv[i][isp] > bloodconc(150., hd[iseg])) {	//August 17, 2010
									ihigh++;
									if (ihigh == 1 && greensverbose) printf("*** Warning: cblood is high +%i", segname[iseg]);
									if (ihigh > 1 && greensverbose) printf("+%i", segname[iseg]);
								}
								blood(cv[i][isp], hd[iseg], &pv[i][isp], &dcdp[i][isp]);
							}
						}
						else {
							pv[i][isp] = cv[i][isp];
							dcdp[i][isp] = 1.;
						}
					}
					if (ineg > 0 || ihigh > 0) if (greensverbose) printf("\n");
					////generate linear system to be solved
					//					//*********************************
					//					float err2 = 0.;
					//					for(i=1; i<=nnv; i++){
					//						iseg = mainseg[i];
					//						rhs[i] = pv[i][isp] - gamma1[iseg][isp]/ds[iseg]*qv[i][isp] - pvt[i][isp] - g0[isp];
					//						err2 += SQR(rhs[i]);
					//					}
					//					err2 = sqrt(err2);
					//					printf(" err2 = %f\n",err2);
					//					//*********************************
					for (i = 1; i <= nnv; i++) {
						iseg = mainseg[i];
						rhs[i] = pv[i][isp] - pvt[i][isp];
						//if (g0method == 2) rhs[i] -= g0[isp];
						for (j = 1; j <= nnv; j++) {
							jseg = mainseg[j];
							mat[i][j] = gvv[i][j] / diff[isp] + al[i][j] / dcdp[i][isp] / qq[iseg] / flowfac;
							if (i == j) mat[i][j] += gamma1[iseg][isp] / ds[iseg];
							rhs[i] += al[i][j] * qv[j][isp] / dcdp[i][isp] / qq[iseg] / flowfac;
							if (oxygen[isp] && lowflow[mainseg[i]]) {  //low q*hd
								if (i == j) mat[i][j] = 1.;
								else mat[i][j] = 0.;
							}
						}
						if (oxygen[isp] == 1 && lowflow[iseg] == 1) rhs[i] = qvprev[i][isp];  //low q*hd
					}
					//solve system of linear algebraic equations: Sum mat[i][j]*qv[j]= rhs[i]
					if (g0method == 1) {			//linear system including G0 as unknown
						for (i = 1; i <= nnv; i++) {
							if (oxygen[isp] == 1 && lowflow[mainseg[i]] == 1) mat[i][nnv + 1] = 0.;  //low q*hd
							else mat[i][nnv + 1] = 1.;
							mat[nnv + 1][i] = 1.;
							mat[nnv + 1][nnv + 1] = 0.;
						}
						if (linmethod == 1) {
							for (i = 1; i <= nnv; i++) rhsg[i][1] = rhs[i];
							rhsg[nnv + 1][1] = -qtsum[isp];
							gaussj(mat, nnv + 1, rhsg, 1);
							for (i = 1; i <= nnv; i++) {
								qv[i][isp] = rhsg[i][1];
								qvsum[isp] += qv[i][isp];
							}
							g0[isp] = rhsg[nnv + 1][1];
						}
						if (linmethod == 2) {
							ludcmp(mat, nnv + 1, indx, &dd);
							for (i = 1; i <= nnv; i++) rhsl[i] = rhs[i];
							rhsl[nnv + 1] = -qtsum[isp];
							lubksb(mat, nnv + 1, indx, rhsl);
							for (i = 1; i <= nnv; i++) {
								qv[i][isp] = rhsl[i];
								qvsum[isp] += qv[i][isp];
							}
							g0[isp] = rhsl[nnv + 1];
						}
						if (linmethod == 3) {
							for (i = 1; i <= nnv; i++) {
								rhsl[i] = rhs[i];
								matx[i] = qv[i][isp];
							}
							rhsl[nnv + 1] = -qtsum[isp];
							matx[nnv + 1] = g0[isp];
							if (useGPU) bicgstabBLASD(mat, rhsl, matx, nnv + 1, bicgstaberr, bicgstabit);
							else bicgstab(mat, rhsl, matx, nnv + 1, bicgstaberr, bicgstabit);
							for (i = 1; i <= nnv; i++) {
								qv[i][isp] = matx[i];
								qvsum[isp] += qv[i][isp];
							}
							g0[isp] += g0fac[isp] * (matx[nnv + 1] - g0[isp]);//Underrelax for configurations w/short initial vessels
						}
					}
					if (g0method == 2) {			//linear system NOT including G0 as unknown
						if (linmethod == 1) {
							for (i = 1; i <= nnv; i++) {
								rhsg[i][1] = rhs[i];
								rhsg[i][2] = -1.;
							}
							gaussj(mat, nnv, rhsg, 2);
							for (i = 1; i <= nnv; i++) {
								qv[i][isp] = rhsg[i][1];
								qvsum[isp] += qv[i][isp];
								if (oxygen[isp] != 1 || lowflow[mainseg[i]] != 1) dqvsumdg0[isp] += rhsg[i][2];
							}
						}
						if (linmethod == 2) {
							ludcmp(mat, nnv, indx, &dd);
							for (i = 1; i <= nnv; i++) rhsl[i] = rhs[i];
							lubksb(mat, nnv, indx, rhsl);
							for (i = 1; i <= nnv; i++) {
								qv[i][isp] = rhsl[i];
								qvsum[isp] += qv[i][isp];
							}
							for (i = 1; i <= nnv; i++) rhsl[i] = -1.;
							lubksb(mat, nnv, indx, rhsl);
							for (i = 1; i <= nnv; i++)
								if (oxygen[isp] != 1 || lowflow[mainseg[i]] != 1) dqvsumdg0[isp] += rhsl[i];
						}
						if (linmethod == 3) {
							for (i = 1; i <= nnv; i++) {
								rhsl[i] = rhs[i];
								matx[i] = qv[i][isp];
							}
							if (useGPU) bicgstabBLASD(mat, rhsl, matx, nnv, bicgstaberr, bicgstabit);
							else bicgstab(mat, rhsl, matx, nnv, bicgstaberr, bicgstabit);
							for (i = 1; i <= nnv; i++) {
								qv[i][isp] = matx[i];
								qvsum[isp] += qv[i][isp];
							}
						}
					}
					//for low q*hd segments, calculate efflux based on change in extravascular oxygen level
					//save values in qvtemp to avoid influence on eval, update qv, underrelax - July 2008
					for (i = 1; i <= nnv; i++) {
						iseg = mainseg[i];
						if (oxygen[isp] == 1 && lowflow[iseg] == 1) {
							for (j = 1; j <= 3; j++) x[j] = ax[j][i] - 0.5 * scos[j][iseg] * ds[iseg];
							p = eval(1, req, x);
							p[isp] = FMAX(p[isp], 0.);
							pv[i][isp] = p[isp] / 2.;
							qvtemp[i] = q[iseg] * flowfac * bloodconc(p[isp], hd[iseg]); //q here (not qq) April 2008
							for (j = 1; j <= 3; j++) x[j] = ax[j][i] + 0.5 * scos[j][iseg] * ds[iseg];
							p = eval(1, req, x);
							p[isp] = FMAX(p[isp], 0.);
							pv[i][isp] += p[isp] / 2.;
							qvtemp[i] -= q[iseg] * flowfac * bloodconc(p[isp], hd[iseg]);
						}
					}
					for (i = 1; i <= nnv; i++) if (oxygen[isp] == 1 && lowflow[mainseg[i]] == 1)
						qv[i][isp] = 0.5 * qvtemp[i] + 0.5 * qvprev[i][isp];
					errvessel[isp] = 0.;
					imaxerr = 0;
					errvesselcount[isp] = 0; //added June 2009
					for (i = 1; i <= nnv; i++) {
						dif = qv[i][isp] - qvprev[i][isp];
						//If qv is large, use relative rather than absolute error  - May 2008
						if (qv[i][isp] != 0.) dif = dif * FMIN(1., epsvessel[isp] / errfac / fabs(qv[i][isp]));
						if (fabs(dif) >= errvessel[isp]) {	
							imaxerrvessel[isp] = mainseg[i];
							errvessel[isp] = fabs(dif);
						}
						if (fabs(dif) > epsvessel[isp]) errvesselcount[isp]++;
					}
					if (greensverbose) printf("Solute %i: qtsum = %f, qvsum = %f\n", isp, qtsum[isp], qvsum[isp]);
					if (greensverbose) printf("Solute %i: kvessel = %i, errvessel_q = %f, imaxerr = %i, g0 = %f\n",
						isp, kvessel, errvessel[isp], imaxerrvessel[isp], g0[isp]);
					if (errvesselcount[isp] > 0) convflagv = 0;
				}
			}
			if (convflagv) goto vesselconv;
		}
		for (isp = 1; isp <= nsp; isp++) if (errvesselcount[isp] > 0)
			if (greensverbose) printf("*** Warning: solute %i, %i vessel source strengths not converged\n",
				isp, errvesselcount[isp]);
	vesselconv:;
		//********************** end of vessel loop *****************************	
		//********************** start of tissue loop *****************************
		//Compute tissue source strengths iteratively by successive relaxation: updated qt values are immediately used.
		//Continually scales up qv values so that their sum equals updated sum of qt values.  Added April 2010.
		//contribution ptv from vessel source strengths qv
		if (is2d) lam = lam2d;
		else lam = lam3d;
		if (useGPU) {
			for (isp = 1; isp <= nsp; isp++) {
				for (i = 1; i <= nnv; i++) qv000[i - 1] = qv[i][isp];
				tissueGPU3c(stepGPU);	//Compute contribution from vessel source strengths
				for (itp = 1; itp <= nnt; itp++) {
					if (is2d) ptv[itp][isp] = fac * pt000[itp - 1] / diff[isp] / w2d;
					else ptv[itp][isp] = fac * pt000[itp - 1] / diff[isp];
				}
			}
		}
		else {
			for (itp = 1; itp <= nnt; itp++) {
				for (isp = 1; isp <= nsp; isp++) ptv[itp][isp] = 0.;
				for (i = 1; i <= nnv; i++) {
					x[1] = ax[1][i] - axt[tisspoints[1][itp]];
					x[2] = ax[2][i] - ayt[tisspoints[2][itp]];
					x[3] = ax[3][i] - azt[tisspoints[3][itp]];
					if (hexperiodic) ndomain = 7;
					else ndomain = 1;
					gtv = 0.;
					for (idomain = 0; idomain < ndomain; idomain++) {
						for (k = 1; k <= 3; k++) {
							y[k] = x[k];
							if (idomain > 0) y[k] += 2. * aphexp * hex_norm[idomain][k];
						}
						dist = length0(y);
						if (dist <= req) {
							if (is2d) gtv += fac / w2d * (2. * log(r2d / req) + 1. - SQR(dist / req));
							else gtv += fac * (1.5 - 0.5 * SQR(dist / req)) / req;
						}
						else {
							if (is2d) gtv += 2. * fac / w2d * log(r2d / dist);
							else gtv += fac / dist;
						}
					}
					for (isp = 1; isp <= nsp; isp++) if (permsolute[isp]) ptv[itp][isp] += gtv / diff[isp] * qv[i][isp];
				}
			}
		}
		for (isp = 1; isp <= nsp; isp++) qvfac[isp] = 1.;
		for (ktissue = 1; ktissue <= nmaxtissue; ktissue++) {
			//Scale all av, qvsum and ptv values so that qvsum = qtsum.  April 2010.
			for (isp = 1; isp <= nsp; isp++) if (permsolute[isp] == 1 && g0method == 1) {
				qvfac[isp] = -qtsum[isp] / qvsum[isp];
				if (fabs(qvfac[isp]) > 2.) qvfac[isp] = 1.;  //avoid extreme values
				if (fabs(qvfac[isp]) < 0.5) qvfac[isp] = 1.;  //avoid extreme values
			}
			convflagt = 1;
			for (isp = 1; isp <= nsp; isp++) {
				qtsum[isp] = 0;
				errtissue[isp] = 0.;
				dqtsumdg0[isp] = 0.;
				errtissuecount[isp] = 0; //added June 2009
			}
			//contribution ptt from tissue source strengths qt
			if (useGPU) {	//GPU method NOT using successive updating of tissue source strengths: uses block updating. May 2016.
				for (iblockstep = 1; iblockstep <= blockstep; iblockstep++) {
					for (isp = 1; isp <= nsp; isp++) if (diffsolute[isp]) {
						for (itp = 1; itp <= nnt; itp++) qt000[itp - 1] = qt[itp][isp];
						tissueGPU1c(stepGPU, blockstep, iblockstep);	//Compute contribution from tissue source strengths. Use blockstep - May 2016
						for (itp = iblockstep; itp <= nnt; itp += blockstep)
							pt[itp][isp] = (1. - lam) * pt[itp][isp] + lam * (ptv[itp][isp] * qvfac[isp] + g0[isp] + pt000[itp - 1] / diff[isp]);	//underrelaxation
					}
					for (itp = iblockstep; itp <= nnt; itp += blockstep) {
						for (isp = 1; isp <= nsp; isp++)	ptpt[isp] = pt[itp][isp];
						tissrate(nsp, ptpt, mtiss, mptiss);
						for (isp = 1; isp <= nsp; isp++) {  //replace qt with value based on updated pt
							dif = mtiss[isp] * vol * sourcefac[itp] - qt[itp][isp];	//sourcefac Feb. 2020
							qt[itp][isp] += dif;
							qtsum[isp] += qt[itp][isp];
							if (diffsolute[isp]) dqtsumdg0[isp] += mptiss[isp] * vol * sourcefac[itp];
							else {	//non-diffusible - use Newton method to solve for pt.  May 2010.
								if (mptiss[isp] == 0.) printf("*** Error: mptiss[%i] = 0 at tissue point %i\n", isp, itp);
								else pt[itp][isp] -= mtiss[isp] / mptiss[isp];
							}
							if (fabs(dif) > errtissue[isp]) {
								errtissue[isp] = fabs(dif);
								imaxerrtissue[isp] = itp;
							}
							if (fabs(dif) > epstissue[isp] * FMAX(lam / 0.05, 1.)) errtissuecount[isp]++;//factor added January 2012 to allow for very small lam values 
						}
					}
				}
			}			
			else {	//non-GPU method using successive updating of tissue source strengths
				for (itp = 1; itp <= nnt; itp++) {
					ix = tisspoints[1][itp];
					iy = tisspoints[2][itp];
					iz = tisspoints[3][itp];
					for (isp = 1; isp <= nsp; isp++) ptt[isp] = 0.;	//all solutes
					for (jtp = 1; jtp <= nnt; jtp++) {
						jx = tisspoints[1][jtp];
						jy = tisspoints[2][jtp];
						jz = tisspoints[3][jtp];
						ixdiff = abs(ix - jx) + 1;
						iydiff = abs(iy - jy) + 1;
						izdiff = abs(iz - jz) + 1;
						for (isp = 1; isp <= nsp; isp++) if (diffsolute[isp]) ptt[isp] += dtt[ixdiff][iydiff][izdiff] * qt[jtp][isp];
					}
					for (isp = 1; isp <= nsp; isp++) {
						if (diffsolute[isp]) pt[itp][isp] = (1. - lam) * pt[itp][isp]
							+ lam * (ptv[itp][isp] * qvfac[isp] + g0[isp] + ptt[isp] / diff[isp]);//underrelaxation
						ptpt[isp] = pt[itp][isp];
					}
					tissrate(nsp, ptpt, mtiss, mptiss);
					for (isp = 1; isp <= nsp; isp++) {  //replace qt with value based on updated pt
						dif = mtiss[isp] * vol * sourcefac[itp] - qt[itp][isp];
						qt[itp][isp] += dif;
						qtsum[isp] += qt[itp][isp];
						if (diffsolute[isp]) dqtsumdg0[isp] += mptiss[isp] * vol * sourcefac[itp];
						else {	//non-diffusible - use Newton method to solve for pt.  May 2010.
							if (mptiss[isp] == 0.) printf("*** Error: mptiss[%i] = 0 at tissue point %i\n", isp, itp);
							else pt[itp][isp] -= mtiss[isp] / mptiss[isp];
						}
						if (fabs(dif) > errtissue[isp]) {
							errtissue[isp] = fabs(dif);
							imaxerrtissue[isp] = itp;
						}
						if (fabs(dif) > epstissue[isp] * FMAX(lam / 0.05, 1.)) errtissuecount[isp]++;//factor added January 2012 to allow for very small lam values 
					}
				}
				/*
				//using block updating, for testing purposes - May 2016
								lam = 0.025;  //0.00175;
								for(iblockstep=1; iblockstep<=blockstep; iblockstep++){
									for(itp=iblockstep; itp<=nnt; itp+=blockstep){
										ix = tisspoints[1][itp];
										iy = tisspoints[2][itp];
										iz = tisspoints[3][itp];
										for(isp=1; isp<=nsp; isp++)	ptt[isp] = 0.;//all solutes
										for(jtp=1; jtp<=nnt; jtp++){
											jx = tisspoints[1][jtp];
											jy = tisspoints[2][jtp];
											jz = tisspoints[3][jtp];
											ixdiff = abs(ix - jx) + 1;
											iydiff = abs(iy - jy) + 1;
											izdiff = abs(iz - jz) + 1;
											for(isp=1; isp<=nsp; isp++) if(diffsolute[isp] == 1) ptt[isp] += dtt[ixdiff][iydiff][izdiff]*qt[jtp][isp];
										}
										for(isp=1; isp<=nsp; isp++){
											if(diffsolute[isp] == 1) pt[itp][isp] = (1.-lam)*pt[itp][isp]
												+ lam*(ptv[itp][isp]*qvfac[isp] + g0[isp] + ptt[isp]/diff[isp]);//underrelaxation
										}
									}
									for(itp=iblockstep; itp<=nnt; itp+=blockstep){
										for(isp=1; isp<=nsp; isp++) ptpt[isp] = pt[itp][isp];
										tissrate(nsp,ptpt,mtiss,mptiss);
										for(isp=1; isp<=nsp; isp++){  //replace qt with value based on updated pt - all solutes
											dif = mtiss[isp]*vol - qt[itp][isp];
											qt[itp][isp] += dif;
											qtsum[isp] += qt[itp][isp];
											if(diffsolute[isp]) dqtsumdg0[isp] += mptiss[isp]*vol;
											else{	//non-diffusible - use Newton method to solve for pt.  May 2010.
												if(mptiss[isp] == 0.) printf("*** Error: mptiss[%i] = 0 at tissue point %i\n",isp,itp);
												else pt[itp][isp] -= mtiss[isp]/mptiss[isp];
											}
											if(fabs(dif) > errtissue[isp]){
												errtissue[isp] = fabs(dif);
												imaxerrtissue[isp] = itp;
											}
											if(fabs(dif) > epstissue[isp]) errtissuecount[isp]++;
										}
									}
								}
				*/
			}
			for (isp = 1; isp <= nsp; isp++) if (diffsolute[isp] == 1) {
				if (greensverbose) printf("Solute %i: qtsum = %f, qvsum = %f\n", isp, qtsum[isp], qvsum[isp] * qvfac[isp]);//May 2010
				if (greensverbose) printf("Solute %i: ktissue = %i, errtissue_q = %f, imaxerr = %i, g0 = %f\n",
					isp, ktissue, errtissue[isp], imaxerrtissue[isp], g0[isp]);
				if (errtissuecount[isp] > 0) convflagt = 0;
			}
			if (kmain > 1 && convflagt) goto tissueconv;  //force full number of iterations when kmain = 1.  May 2010
		}
		for (isp = 1; isp <= nsp; isp++) if (errtissuecount[isp] > 0)
			if (greensverbose) printf("*** Warning: solute %i, %i tissue source strengths not converged\n", isp, errtissuecount[isp]);
	tissueconv:;
		//Print log file.  April 2010
		ofp2 = fopen("GreensLog.txt", "a");
		kvessel = IMIN(kvessel, nmaxvessel);
		ktissue = IMIN(ktissue, nmaxtissue);
		fprintf(ofp2, "\n----- kmain = %i, kvessel = %i, ktissue = %i -----\n", kmain, kvessel, ktissue);
		for (isp = 1; isp <= nsp; isp++) {
			if (diffsolute[isp] == 1) fprintf(ofp2, "Solute %i: qtsum = %f, qvsum = %f, g0 = %f\n",
				isp, qtsum[isp], qvsum[isp] * qvfac[isp], g0[isp]);
			if (permsolute[isp] == 1) fprintf(ofp2, "Solute %i: errvessel_q = %f, imaxerr = %i\n",
				isp, errvessel[isp], segname[imaxerrvessel[isp]]);
			if (diffsolute[isp] == 1) fprintf(ofp2, "Solute %i: errtissue_q = %f, imaxerr = %i\n",
				isp, errtissue[isp], imaxerrtissue[isp]);
		}
		fclose(ofp2);
		//********************** end of tissue loop *****************************
		//Update g0.  If permsolute[isp] != 1, always use method 2.
		//Method 2 is based on derivative wrt g0 - new version September 2009 - automatic estimation of g0fac
		for (isp = 1; isp <= nsp; isp++) g0facnew[isp] = 0.;
		if (useGPU) {
			for (isp = 1; isp <= nsp; isp++) if (diffsolute[isp] == 1 && (g0method == 2 || permsolute[isp] == 0)) {
				for (itp = 1; itp <= nnt; itp++) {
					for (i = 1; i <= nsp; i++) ptpt[i] = pt[itp][i];
					tissrate(nsp, ptpt, mtiss, mptiss);
					qt000[itp - 1] = mptiss[isp] * sourcefac[itp];
				}
				tissueGPU1c(stepGPU, 1, 1);	//Compute contribution from tissue source strength derivatives
				for (itp = 1; itp <= nnt; itp++) g0facnew[isp] += pt000[itp - 1] / diff[isp] * vol;
			}
		}
		else {
			for (itp = 1; itp <= nnt; itp++) {
				for (isp = 1; isp <= nsp; isp++) ptpt[isp] = pt[itp][isp];
				tissrate(nsp, ptpt, mtiss, mptiss);
				ix = tisspoints[1][itp];
				iy = tisspoints[2][itp];
				iz = tisspoints[3][itp];
				for (jtp = 1; jtp <= nnt; jtp++) {
					jx = tisspoints[1][jtp];
					jy = tisspoints[2][jtp];
					jz = tisspoints[3][jtp];
					ixdiff = abs(ix - jx) + 1;
					iydiff = abs(iy - jy) + 1;
					izdiff = abs(iz - jz) + 1;
					for (isp = 1; isp <= nsp; isp++) if (diffsolute[isp] == 1)
						g0facnew[isp] += dtt[ixdiff][iydiff][izdiff] / diff[isp] * mptiss[isp] * vol * sourcefac[itp];
				}
			}
		}
		for (isp = 1; isp <= nsp; isp++) if (diffsolute[isp] == 1 && (g0method == 2 || permsolute[isp] == 0)) {
			g0facnew[isp] = 1. / (1. - g0facnew[isp] / nnt);
			dqsumdg0 = FMIN(dqtsumdg0[isp], 0.) * g0facnew[isp];
			if (fabs(dqsumdg0) > 1.e-6) {
				dif = (qvsum[isp] + qtsum[isp]) / dqsumdg0 * g0fac[isp];//This g0fac should normally be 1.0.  September 2009
				g0[isp] -= dif;
			}
		}
		//Convergence based on changes in pv, pt and g0
		convflag = 1;
		if (greensverbose) printf("\n");
		for (isp = 1; isp <= nsp; isp++) {
			err = 0.;
			imaxerr = 0;
			if (permsolute[isp] == 1) for (i = 1; i <= nnv; i++) {
				dif = fabs(pv[i][isp] - pvprev[i][isp]);
				if (dif > err) {
					imaxerr = mainseg[i];
					err = dif;
				}
			}
			errvessel[isp] = err;
			imaxerrvessel[isp] = imaxerr;
			err = 0.;
			imaxerr = 0;
			if (diffsolute[isp] == 1) for (itp = 1; itp <= nnt; itp++) {
				dif = fabs(pt[itp][isp] - ptprev[itp][isp]);
				if (dif > err) {
					imaxerr = itp;
					err = dif;
				}
			}
			errtissue[isp] = err;
			imaxerrtissue[isp] = imaxerr;
			if (errvessel[isp] > err) {
				imaxerr = imaxerrvessel[isp];
				err = errvessel[isp];
			}
			else imaxerr = -imaxerr;
			dif = fabs(g0[isp] - g0old[isp]);
			if (dif > err) {
				imaxerr = 0;
				err = dif;
			}
			if (greensverbose) printf("Solute %i: err = %f, imaxerr = %i (- for tissue point)\n", isp, err, imaxerr);
			if (greensverbose && imaxerr > 0) if (lowflow[imaxerr]) printf("Solute %i: max error is at a low-flow segment\n", isp);
			if (err > eps[isp]) convflag = 0;
		}
		//Print log file - April 2010
		ofp2 = fopen("GreensLog.txt", "a");
		for (isp = 1; isp <= nsp; isp++) {
			if (permsolute[isp] == 1) fprintf(ofp2, "Solute %i: errvessel_p = %f, imaxerr = %i\n",
				isp, errvessel[isp], segname[imaxerrvessel[isp]]);
			fprintf(ofp2, "Solute %i: errtissue_p = %f, imaxerr = %i\n",
				isp, errtissue[isp], imaxerrtissue[isp]);
		}
		fclose(ofp2);
		tfinish1 = clock();
		duration = (float)(tfinish1 - tstart1) / CLOCKS_PER_SEC;
		if (greensverbose) printf("\nkmain = %i, %2.3f seconds for step\n", kmain, duration);
		if (convflag && convflagv && convflagt) goto mainconv;
	}
	printf("\n*** Warning: tissue or vessel solute levels not converged\n");
	if (greensverbose) for (isp = 1; isp <= nsp; isp++) if (diffsolute[isp] == 1) {
		printf("Solute %i: qtsum = %f, qvsum = %f, g0 = %f\n", isp, qtsum[isp], qvsum[isp], g0[isp]);
		printf("Solute %i: kvessel = %i, errvessel_q = %f, imaxerr = %i\n", isp, kvessel, errvessel[isp], imaxerrvessel[isp]);
		printf("Solute %i: ktissue = %i, errtissue_q = %f, imaxerr = %i\n", isp, ktissue, errtissue[isp], imaxerrtissue[isp]);
	}
mainconv:;
	//********************** end of main loop *****************************
	tfinish = clock();
	duration = (float)(tfinish - tstart) / CLOCKS_PER_SEC;
	printf("\nGreens: %i iterations, %2.1f seconds for main loop\n", kmain, duration);
	ofp2 = fopen("GreensLog.txt", "a");
	fprintf(ofp2, "\n%i iterations, %2.1f seconds for main loop\n", kmain, duration);
	fclose(ofp2);
	ofp1 = fopen("MasterLog.out", "a");
	fprintf(ofp1, "%6i%8.2f%6i%6i%6i%6i   %6i   %6i   %6i   ", imain, duration, nseg, nnv, nnt, kmain, convflag, convflagv, convflagt);
	for (isp = 1; isp <= nsp; isp++) if (diffsolute[isp]) fprintf(ofp1, "%2i %8.2e %8.2e", isp, errvessel[isp], errtissue[isp]);
	fprintf(ofp1, "\n");
	fclose(ofp1);

	//Scale all qv values so that qvsum = qtsum.  April 2010.
	for (isp = 1; isp <= nsp; isp++) if (permsolute[isp] == 1 && g0method == 1) {
		qvsum[isp] *= qvfac[isp];
		for (i = 1; i <= nnv; i++) qv[i][isp] *= qvfac[isp];
	}
	//general output file, not read by readsources
	ofp = fopen("GreensRes.out", "w");
	fprintf(ofp, "%i %i %i %i %i %i\n", nnv, nseg, mxx, myy, mzz, nnt);
	fprintf(ofp, "Total flow rate into region q0 = %f nl/min\n", q0);
	for (isp = 1; isp <= nsp; isp++) fprintf(ofp, "g0[%i] = %f\n", isp, g0[isp]);
	fprintf(ofp, "\n");
	//extravascular solute levels
	for (isp = 1; isp <= nsp; isp++) if (diffsolute[isp] == 1) {
		fprintf(ofp, "\nSolute %i\n", isp);
		fprintf(ofp, "Segment");
		if (permsolute[isp] == 1) fprintf(ofp, "Efflux Pvessel Ptissue Cvessel");
		fprintf(ofp, "\n");
		for (i = 1; i <= nnv; i++) {
			pev[i][isp] = pvt[i][isp];
			if (g0method == 1 && permsolute[isp] == 1) pev[i][isp] += g0[isp];
			if (permsolute[isp] == 1) for (j = 1; j <= nnv; j++) pev[i][isp] += gvv[i][j] * qv[j][isp] / diff[isp];
		}
		for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] >= 3 && segtyp[iseg] <= 5) {
			qvseg[iseg][isp] = 0.;
			pevseg[iseg][isp] = 0.;
			pvseg[iseg][isp] = 0.;
		}
		for (i = 1; i <= nnv; i++) {
			iseg = mainseg[i];
			if (permsolute[isp] == 1) fprintf(ofp, "%4i %4i %10.4f %10.4f %10.4f %10.4f\n",
				i, iseg, qv[i][isp], pv[i][isp], pev[i][isp], cv[i][isp]);
			qvseg[iseg][isp] += qv[i][isp];
			pevseg[iseg][isp] += pev[i][isp] / nspoint[iseg];
			pvseg[iseg][isp] += pv[i][isp] / nspoint[iseg];
		}
		fprintf(ofp, "Solute %i: qtsum = %f, qvsum = %f\n", isp, qtsum[isp], qvsum[isp]);
	}
	for (isp = 1; isp <= nsp; isp++) if (permsolute[isp] == 1) {
		fprintf(ofp, "Solute %i: segment length pvseg pevseg qvseg gamma\n", isp);
		for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] >= 3 && segtyp[iseg] <= 5)
			fprintf(ofp, "%4i %10.4f %10.4f %10.4f %10.4f %10.4f\n",
				segname[iseg], lseg[iseg], pvseg[iseg][isp], pevseg[iseg][isp], qvseg[iseg][isp], gamma1[iseg][isp]);
	}
	fclose(ofp);

	//Write output files that allow restart of program without running greens
	//These files are read by readsources
	ofp = fopen("TissueSources.out", "w");
	fprintf(ofp, "%i %i %i %i %f\n", mxx, myy, mzz, nnt, req);
	fprintf(ofp, "X, Y, Z coords of source points");
	for (i = 1; i <= mxx; i++) {
		if (i % 10 == 1) fprintf(ofp, "\n");
		fprintf(ofp, " %g", axt[i]);
	}
	for (i = 1; i <= myy; i++) {
		if (i % 10 == 1) fprintf(ofp, "\n");
		fprintf(ofp, " %g", ayt[i]);
	}
	for (i = 1; i <= mzz; i++) {
		if (i % 10 == 1) fprintf(ofp, "\n");
		fprintf(ofp, " %g", azt[i]);
	}
	fprintf(ofp, "\nTissue point xyz indices");
	for (itp = 1; itp <= nnt; itp++) {
		if (itp % 10 == 1) fprintf(ofp, "\n");
		fprintf(ofp, "%4i %4i %4i", tisspoints[1][itp], tisspoints[2][itp], tisspoints[3][itp]);
	}
	for (isp = 1; isp <= nsp; isp++) if (diffsolute[isp] == 1) {
		fprintf(ofp, "\n%g = g0, source strengths for solute %i", g0[isp], isp);
		for (i = 1; i <= nnt; i++) {
			if (i % 10 == 1) fprintf(ofp, "\n");
			fprintf(ofp, " %10g", qt[i][isp]);
		}
	}
	for (isp = 1; isp <= nsp; isp++) if (diffsolute[isp] == 1) {
		fprintf(ofp, "\n%g = g0, past 1 timestep source strengths for solute %i", g0[isp], isp);
		for (i = 1; i <= nnt; i++) {
			if (i % 10 == 1) fprintf(ofp, "\n");
			fprintf(ofp, " %10g", qt[i][isp]);
		}
	}
	fclose(ofp);

	ofp = fopen("VesselSources.out", "w");
	fprintf(ofp, "%i %i\n", nseg, nnv);
	fprintf(ofp, "Segment start coords, direction cosines, length\n");
	for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] >= 3 && segtyp[iseg] <= 5)
		fprintf(ofp, "%10g %10g %10g %10g %10g %10g %10g\n",
			cnode[1][ista[iseg]], cnode[2][ista[iseg]], cnode[3][ista[iseg]], scos[1][iseg], scos[2][iseg], scos[3][iseg], ds[iseg]);
	fprintf(ofp, "Source point coords");
	for (i = 1; i <= nnv; i++) {
		if (i % 3 == 1) fprintf(ofp, "\n");
		fprintf(ofp, "  %10g %10g %10g", ax[1][i], ax[2][i], ax[3][i]);
	}
	fprintf(ofp, "\nMain segment numbers of subsegments");
	for (i = 1; i <= nnv; i++) {
		if (i % 20 == 1) fprintf(ofp, "\n");
		fprintf(ofp, " %i", mainseg[i]);
	}
	for (isp = 1; isp <= nsp; isp++) if (permsolute[isp] == 1) {
		fprintf(ofp, "\nSource strengths for solute %i", isp);
		for (i = 1; i <= nnv; i++) {
			if (i % 10 == 1) fprintf(ofp, "\n");
			fprintf(ofp, " %10g", qv[i][isp]);
		}
	}
	fclose(ofp);

	ofp = fopen("TissueLevels.out", "w");
	for (isp = 1; isp <= nsp; isp++) {
		pmax[isp] = -1.e8;
		pmean[isp] = 0.;
		pmin[isp] = 1.e8;
		fprintf(ofp, "Solute %i", isp);
		for (itp = 1; itp <= nnt; itp++) {
			if (itp % 10 == 1) fprintf(ofp, "\n");
			fprintf(ofp, "%12f ", pt[itp][isp]);
			pmean[isp] += pt[itp][isp];
			pmax[isp] = FMAX(pt[itp][isp], pmax[isp]);
			pmin[isp] = FMIN(pt[itp][isp], pmin[isp]);
		}
		pmean[isp] = pmean[isp] / nnt;
		fprintf(ofp, "\n%f %f %f Solute %i: pmean, pmin, pmax\n", pmean[isp], pmin[isp], pmax[isp], isp);
	}
	fclose(ofp);

	ofp = fopen("VesselLevels.out", "w");
	fprintf(ofp, "lowflow");
	i = 0;
	for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] >= 3 && segtyp[iseg] <= 5) {
		i++;
		if (i % 50 == 1) fprintf(ofp, "\n");
		fprintf(ofp, "%2i", lowflow[iseg]);
	}
	fprintf(ofp, "\n");
	for (isp = 1; isp <= nsp; isp++) {
		i = 0;
		fprintf(ofp, "Solute %i pvseg", isp);
		for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] >= 3 && segtyp[iseg] <= 5) {
			i++;
			if (i % 10 == 1) fprintf(ofp, "\n");
			fprintf(ofp, "%10.4f ", pvseg[iseg][isp]);
		}
		fprintf(ofp, "\n");
		fprintf(ofp, "Solute %i pevseg", isp);
		i = 0;
		for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] >= 3 && segtyp[iseg] <= 5) {
			i++;
			if (i % 10 == 1) fprintf(ofp, "\n");
			fprintf(ofp, "%10.4f ", pevseg[iseg][isp]);
		}
		fprintf(ofp, "\n");
	}
	fclose(ofp);
}
