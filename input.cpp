/************************************************************************
input - for AngioAdapt
TWS October 07
Pressures in mmHg, flows in um^3/s, viscosity in cP
Comments updated J. Alberding, May 2009
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void input()
{
	extern int nseg, nnod, nnnodbc, nmaxvessel, nmaxtissue, nmax, g0method, linmethod, useGPU, imainmax;
	extern int mxx, myy, mzz, nseg, nnod, nsp, nnodbc, nre, nodsegm, nsprmax;
	extern int nitmax1, nitmax, varyviscosity, phaseseparation, eliminate, seed;
	extern int adaptd, inftrans, angio, tensionm, outbounmethod;
	extern int ncontourplanes, * slsegdiv, * nsl1, * nsl2;
	extern int is2d, angio3d, hexperiodic;
	extern int* segname, * segtyp, * bcnod, * bcnodname, * bctyp, * nodname, * permsolute, * nl;
	extern int* diffsolute, * oxygen, * nresis; //added TWS2010
	extern int** segnodname;

	extern float pi1, alphab, p50, fn, cs, dalpha, facfp, vplas, mcvcorr, q0, errfac, lowflowcrit, time1;
	extern float plow, phigh, clowfac, chighfac, pphighfac;//added January 2012
	extern float alx, aly, alz, vol, lb, maxl, req;
	extern float tol, qtol, hdtol, omega, optw, optlam;
	extern float eqtime, diamthresh, qinsum_min;
	extern float smax, Vmax, min_length, thresh, variance;
	extern float new_diam, increment, tolerance_1, tolerance_2, threshn, thresh50;
	extern float kgGF, kgVAF, kv, lambdatension, tensionrate, r0, theta0, constvisc, consthd;
	extern float km1, ks1, ranks1, L1, Sumax, Su50, Sd50, tauref1, Qref1, extlength1in, extlength1out, timestep1;
	extern float* axt, * ayt, * azt, * ds, * g0, * diff, * g0fac, * pref;
	extern float* diam, * lseg, * rseg, * q, * qq, * hd, * bcprfl, * bchd;
	extern float* bifpar, * cpar, * viscpar, * fahrpar, * adaptpar;
	extern float** xsl0, ** xsl1, ** xsl2, * clmin, * clint, * cl, ** zv, *** psl;
	extern float** cnode, ** bcp, ** resis, ** resisdiam, ** tissparam;

	float totalq, mcv;
	int i, inodbc, iseg, isp, max = 200, nlmax;
	int icp, nsl1max, nsl2max;
	char bb[200];
	FILE* ifp, * ofp2;
	ofp2 = fopen("ArrayChanges.out", "a");

	bifpar = vector(1, 3);
	cpar = vector(1, 4);
	viscpar = vector(1, 6);
	fahrpar = vector(1, 4);
	adaptpar = vector(1, 4);

	ifp = fopen("RheolParams.dat", "r");
	fscanf(ifp, "%f %f %f%*[^\n]", &bifpar[1], &bifpar[2], &bifpar[3]);
	fscanf(ifp, "%f %f %f %f%*[^\n]", &cpar[1], &cpar[2], &cpar[3], &cpar[4]);
	fscanf(ifp, "%f %f %f %f %f %f%*[^\n]", &viscpar[1], &viscpar[2], &viscpar[3], &viscpar[4], &viscpar[5], &viscpar[6]);
	fscanf(ifp, "%f %f %f %f%*[^\n]", &fahrpar[1], &fahrpar[2], &fahrpar[3], &fahrpar[4]);
	fscanf(ifp, "%i %f %f%*[^\n]", &nitmax, &tol, &omega);
	fscanf(ifp, "%i %f %f%*[^\n]", &nitmax1, &qtol, &hdtol);
	fscanf(ifp, "%f %f%*[^\n]", &optw, &optlam);
	fscanf(ifp, "%f %f %f%*[^\n]", &constvisc, &vplas, &mcv);
	fscanf(ifp, "%f%*[^\n]", &consthd);
	fscanf(ifp, "%i%*[^\n]", &varyviscosity);
	fscanf(ifp, "%i", &phaseseparation);
	fclose(ifp);

	ifp = fopen("AdaptParams.dat", "r");
	fgets(bb, max, ifp);
	printf("%s\n", bb);
	fscanf(ifp, "%i%*[^\n]", &imainmax);
	fscanf(ifp, "%f%*[^\n]", &timestep1);
	fscanf(ifp, "%i%*[^\n]", &adaptd);
	fscanf(ifp, "%i%*[^\n]", &inftrans);
	fscanf(ifp, "%i%*[^\n]", &eliminate);
	fscanf(ifp, "%i%*[^\n]", &angio);
	fscanf(ifp, "%i%*[^\n]", &tensionm);
	fscanf(ifp, "%f %f %f %f%*[^\n]", &adaptpar[1], &adaptpar[2], &adaptpar[3], &adaptpar[4]);
	fscanf(ifp, "%f%*[^\n]", &km1);
	fscanf(ifp, "%f%*[^\n]", &Sumax);
	fscanf(ifp, "%f %f%*[^\n]", &Sd50, &Su50);
	fscanf(ifp, "%f%*[^\n]", &ks1);
	fscanf(ifp, "%f%*[^\n]", &ranks1);
	fscanf(ifp, "%f%*[^\n]", &tauref1);
	fscanf(ifp, "%f%*[^\n]", &Qref1);
	fscanf(ifp, "%f%*[^\n]", &L1);
	fscanf(ifp, "%f%*[^\n]", &extlength1in);
	fscanf(ifp, "%f%*[^\n]", &extlength1out);
	fscanf(ifp, "%f%*[^\n]", &diamthresh);
	fscanf(ifp, "%f%*[^\n]", &eqtime);
	fscanf(ifp, "%i%*[^\n]", &outbounmethod);
	fscanf(ifp, "%f%*[^\n]", &qinsum_min);
	fclose(ifp);

	ifp = fopen("AngioParams.dat", "r");
	fscanf(ifp, "%f%*[^\n]", &smax);
	fscanf(ifp, "%i%*[^\n]", &nsprmax);
	fscanf(ifp, "%f%*[^\n]", &Vmax);
	fscanf(ifp, "%f%*[^\n]", &min_length);
	fscanf(ifp, "%f%*[^\n]", &thresh);
	fscanf(ifp, "%f%*[^\n]", &variance);
	fscanf(ifp, "%f%*[^\n]", &new_diam);
	fscanf(ifp, "%f%*[^\n]", &increment);
	fscanf(ifp, "%f%*[^\n]", &tolerance_1);
	fscanf(ifp, "%f%*[^\n]", &tolerance_2);
	fscanf(ifp, "%f %f%*[^\n]", &threshn, &thresh50);
	fscanf(ifp, "%f %f%*[^\n]", &kgGF, &kgVAF);
	fscanf(ifp, "%f%*[^\n]", &kv);
	fscanf(ifp, "%f%*[^\n]", &lambdatension);
	fscanf(ifp, "%f%*[^\n]", &tensionrate);
	fscanf(ifp, "%i%*[^\n]", &seed);
	fscanf(ifp, "%f%*[^\n]", &r0);
	fscanf(ifp, "%f%*[^\n]", &theta0);
	fscanf(ifp, "%i%*[^\n]", &angio3d);
	fscanf(ifp, "%i%*[^\n]", &hexperiodic);
	fclose(ifp);

	ifp = fopen("SoluteParams.dat", "r");
	fgets(bb, max, ifp);
	printf("%s\n", bb);
	fscanf(ifp, "%i %i %i%*[^\n]", &g0method, &linmethod, &useGPU);
	fscanf(ifp, "%i %i %i%*[^\n]", &nmaxvessel, &nmaxtissue, &nmax);
	fscanf(ifp, "%f%*[^\n]", &errfac);
	fscanf(ifp, "%f%*[^\n]", &lowflowcrit);
	fscanf(ifp, "%f%*[^\n]", &p50);
	fscanf(ifp, "%f%*[^\n]", &fn);
	fscanf(ifp, "%f%*[^\n]", &cs);
	fscanf(ifp, "%f%*[^\n]", &alphab);
	fscanf(ifp, "%f%*[^\n]", &q0);
	fscanf(ifp, "%i", &nsp);
	fgets(bb, max, ifp);

	permsolute = ivector(1, nsp);
	diffsolute = ivector(1, nsp);
	oxygen = ivector(1, nsp);
	pref = vector(1, nsp);
	diff = vector(1, nsp);
	g0 = vector(1, nsp);
	g0fac = vector(1, nsp);
	tissparam = matrix(1, 3, 1, nsp);

	for (isp = 1; isp <= nsp; isp++) {
		fgets(bb, max, ifp);
		fscanf(ifp, "%i %i %i%*[^\n]", &permsolute[isp], &diffsolute[isp], &oxygen[isp]);
		fscanf(ifp, "%f%*[^\n]", &pref[isp]);
		fscanf(ifp, "%f%*[^\n]", &diff[isp]);
		diff[isp] = diff[isp] * 1.e8;
		for (i = 1; i <= 3; i++) fscanf(ifp, "%f%*[^\n]", &tissparam[i][isp]);
		fscanf(ifp, "%f%*[^\n]", &g0[isp]);
		fscanf(ifp, "%f", &g0fac[isp]);
		fgets(bb, max, ifp);
	}
	fclose(ifp);
	tissparam[1][1] = tissparam[1][1] / 6000.;	//correct units of oxygen consumption rate

//conversion factor for flow from pressure based on Poiseuille flows measured in nl/min
	facfp = pi1 * 1333. / 128. / 0.01 * 60. / 1.e6;
	mcvcorr = pow(92. / mcv, 0.33333);
	//parameters for blood
	plow = 0.1 * p50;
	phigh = 5. * p50;
	clowfac = cs * (1.0 - 1.0 / (1.0 + pow((plow / p50), fn)));
	chighfac = cs * (1.0 - 1.0 / (1.0 + pow((phigh / p50), fn)));
	pphighfac = cs * fn / p50 * pow(phigh / p50, (fn - 1)) / SQR(1. + pow(phigh / p50, fn));
	//intravascular diffusion resistance
	resisdiam = matrix(1, 20, 1, nsp);
	resis = matrix(1, 20, 1, nsp);
	nresis = ivector(1, nsp);

	ifp = fopen("IntravascRes.dat", "r");
	for (isp = 1; isp <= nsp; isp++) {
		fscanf(ifp, "%i", &nresis[isp]);
		if (nresis[isp] > 20) printf("*** Error: too many points in IntravascRes.dat, nresis = %i > 20\n", nresis[isp]);
		fgets(bb, max, ifp);
		if (nresis[isp] > 0) {
			fgets(bb, max, ifp);
			for (i = 1; i <= nresis[isp]; i++) fscanf(ifp, "%f %f", &resisdiam[i][isp], &resis[i][isp]);
		}
	}
	fclose(ifp);

	//network data file
	ifp = fopen("Network.dat", "r");
	time1 = 0.;
	fgets(bb, max, ifp);
	printf("%s\n", bb);
	fscanf(ifp, "%f %f %f%*[^\n]", &alx, &aly, &alz);
	fscanf(ifp, "%i %i %i%*[^\n]", &mxx, &myy, &mzz);
	fscanf(ifp, "%f%*[^\n]", &lb);
	fscanf(ifp, "%f%*[^\n]", &maxl);
	fscanf(ifp, "%i%*[^\n]", &nodsegm);
	fscanf(ifp, "%i%*[^\n]", &nseg);
	fgets(bb, max, ifp);
	fgets(bb, max, ifp);
	segname = ivector(1, nseg); fprintf(ofp2, "segname,%i\n", nseg);
	segtyp = ivector(1, nseg); fprintf(ofp2, "segtyp,%i\n", nseg);
	segnodname = imatrix(1, 2, 1, nseg); fprintf(ofp2, "segnodname,%i\n", nseg);
	diam = vector(1, nseg); fprintf(ofp2, "diam,%i\n", nseg);
	q = vector(1, nseg); fprintf(ofp2, "q,%i\n", nseg);
	hd = vector(1, nseg); fprintf(ofp2, "hd,%i\n", nseg);
	for (iseg = 1; iseg <= nseg; iseg++) fscanf(ifp, "%i %i %i %i %f %f %f%*[^\n]",
		&segname[iseg], &segtyp[iseg], &segnodname[1][iseg], &segnodname[2][iseg], &diam[iseg], &q[iseg], &hd[iseg]);
	//number of nodes in vessel network
	fgets(bb, max, ifp);
	fscanf(ifp, "%i%*[^\n]", &nnod);
	fgets(bb, max, ifp);
	fgets(bb, max, ifp);
	//coordinates of nodes
	nodname = ivector(1, nnod); fprintf(ofp2, "nodname,%i\n", nnod);
	cnode = matrix(1, 3, 1, nnod); fprintf(ofp2, "cnode,%i\n", nnod);
	for (i = 1; i <= nnod; i++)	fscanf(ifp, "%i %f %f %f*[^\n]", &nodname[i], &cnode[1][i], &cnode[2][i], &cnode[3][i]);
	//boundary nodes
	fscanf(ifp, "%i%*[^\n]", &nnodbc);
	fgets(bb, max, ifp);
	fgets(bb, max, ifp);
	bcnod = ivector(1, nnodbc);
	bcnodname = ivector(1, nnodbc);
	bctyp = ivector(1, nnodbc);
	bcprfl = vector(1, nnodbc);
	bchd = vector(1, nnodbc);
	bcp = matrix(1, nnodbc, 1, nsp);
	totalq = 0.;
	for (inodbc = 1; inodbc <= nnodbc; inodbc++) {
		fscanf(ifp, "%i %i %f %f", &bcnodname[inodbc], &bctyp[inodbc], &bcprfl[inodbc], &bchd[inodbc]);
		for (isp = 1; isp <= nsp; isp++) if (permsolute[isp] == 1) fscanf(ifp, "%f", &bcp[inodbc][isp]);
		if (bctyp[inodbc] == 2 && bcprfl[inodbc] > 0.) totalq = totalq + bcprfl[inodbc];
		fscanf(ifp, "%*[^\n]");		//TWS November 2010
	}
	fclose(ifp);
	//vol = volume represented by each tissue point; req = radius of equivalent sphere
	vol = alx * aly * alz / (mxx * myy * mzz);
	if (mzz == 1) req = pow(vol * 1. / alz / pi1, 0.5);//2d version, corrected January 2012
	else req = pow(vol * 0.75 / pi1, 0.333333);
	//compute coordinates of tissue points
	axt = vector(1, mxx);
	ayt = vector(1, myy);
	azt = vector(1, mzz);
	for (i = 1; i <= mxx; i++) axt[i] = (i - 0.5) * alx / mxx;
	for (i = 1; i <= myy; i++) ayt[i] = (i - 0.5) * aly / myy;
	for (i = 1; i <= mzz; i++) azt[i] = (i - 0.5) * alz / mzz;
	//Read parameters for slice on which P is computed for contour plot
	nl = ivector(1, nsp);
	xsl0 = matrix(1, 3, 1, 3);
	xsl1 = matrix(1, 3, 1, 3);
	xsl2 = matrix(1, 3, 1, 3);
	slsegdiv = ivector(1, 3);
	nsl1 = ivector(1, 3);
	nsl2 = ivector(1, 3);
	clmin = vector(1, nsp);
	clint = vector(1, nsp);

	ifp = fopen("ContourParams.dat", "r");
	fscanf(ifp, "%i%*[^\n]", &ncontourplanes);
	if (ncontourplanes == 1 || ncontourplanes == 2 || ncontourplanes == 3) {
		for (icp = 1; icp <= ncontourplanes; icp++)
			fscanf(ifp, "%f %f %f %i", &xsl0[1][icp], &xsl0[2][icp], &xsl0[3][icp], &slsegdiv[icp]);
		fscanf(ifp, "%*[^\n]");
		nsl1max = 1;
		for (icp = 1; icp <= ncontourplanes; icp++) {
			fscanf(ifp, "%f %f %f %i", &xsl1[1][icp], &xsl1[2][icp], &xsl1[3][icp], &nsl1[icp]);
			if (nsl1[icp] > nsl1max) nsl1max = nsl1[icp];
		}
		fscanf(ifp, "%*[^\n]");
		nsl2max = 1;
		for (icp = 1; icp <= ncontourplanes; icp++) {
			fscanf(ifp, "%f %f %f %i", &xsl2[1][icp], &xsl2[2][icp], &xsl2[3][icp], &nsl2[icp]);
			if (nsl2[icp] > nsl2max) nsl2max = nsl2[icp];
		}
		fscanf(ifp, "%*[^\n]");
		nlmax = 1;
		for (isp = 1; isp <= nsp; isp++) {
			fscanf(ifp, "%f %f %i%*[^\n]", &clmin[isp], &clint[isp], &nl[isp]);
			if (nl[isp] > nlmax) nlmax = nl[isp];
		}
	}
	else printf("*** Error: invalid value of ncontourplanes, must be 1, 2 or 3\n");
	fclose(ifp);

	cl = vector(1, nlmax);
	zv = matrix(1, nsl1max, 1, nsl2max);
	psl = f3tensor(1, nsl1max, 1, nsl2max, 1, nsp);

	fclose(ofp2);
}