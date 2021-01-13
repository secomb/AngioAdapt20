/************************************************************************
AngioAdapt19_GPU
program to simulate adaptation and angiogenesis in microvascular networks
Version in C initiated by T.W. Secomb, October 07.
Based on work of A.R. Pries, R. Hsu and J.P. Alberding
Revised TWS 2010, 2011, 2016, 2017
Version for periodic network structure - work started February 2018
Version for CUDA 10.1, August 2019
_____________________________________

Segment types, updated July 2016, Feb. 2017
segtyp = 0: new active segment (nonflowing)
segtyp = 1: existing active segment (nonflowing)
segtyp = 3: flowing, can adapt, can be pruned, cannot sprout (new 2016), cannot make new connections
segtyp = 4: flowing, cannot adapt, cannot be pruned, cannot sprout, cannot make new connections, includes boundary segments
segtyp = 5: flowing, can adapt, can be pruned, can sprout
segtyp = 10: pruned, not included in network
______________________________________

Numbering scheme for periodic network structure. February 2018.
Each segment and node has a "base number" and an "extended number."
The base number represents only the position in the base domain and is used for Greens calculations.
The extended number represents the position in a specific domain (base domain or neighboring domain).
Example: iseg = number in base domain, isegx = number in extended domain
isegx = 10*iseg + idom
where idom = 0 in base domain, idom = 1,2,3,.. in neighboring domains.
Note that iseg = isegx/10 and idom = isegx%10 in integer arithmetic.
To identify segments crossing the domain boundary, use a variable boundseg.
For a segment within the domain, boundseg = 0.
For a segment crossing boundary 1, boundseg = 1, etc.
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//#include <cutil_inline.h>
#include <cuda_runtime.h>
#include "nrutil.h"
#include <time.h>

#if defined(__linux__)
	// Requires c++17 support, should be included in all current linux releases
	#include <experimental/filesystem> 
	namespace fs = std::experimental::filesystem::v1;
#elif defined(__APPLE__)
	// Requires removal of the -lstdc++fs flag from makefile
	#include <filesystem>
	namespace fs = std::filesystem;
#elif defined(_WIN32)    //Windows version
	#include <Windows.h>
#endif

void input();
void analyzenet(int allsegs);
void picturenetwork(float* nodvar, float* segvar, const char fname[]);
void greens();
void contourO2GF(int imain);
void histogram();
void setuparrays0();
void setuparrays1(int nseg, int nnod);
void setuparrays2(int nnv, int nnt, int actnum);
void resizearrays1(int nseg_new, int nnod_new);
void resizearrays2(int nnv, int nnt, int actnum);
int outboun(int method);
void flowtest(int allsegs);
void segreclass();
void analyzenet(int allsegs);
void flow();
void setupgreens();
void greens();
void conductconvect();
void adapt();
void output();
void signalplot(const char fname[]);
void actgrowth();
void sprgrowth();
void actsprupdate();
void straight_sht_seg_fixer();
void picturenetwork(float* nodvar, float* segvar, const char fname[]);
void write_network(int imain);
void netstats(float po2cutoff);
void netstatsinit();
void histogram();
void contourO2GF(int imain);
void prep_for_adapt(int inod_to);
void tension();
void cmgui(int imain);
void adaptsignals(int runs);
void setup_hexperiodic();
void remap_hexperiodic();
float length(float* x);
void evalGFinit();
void evalGF(int isp, float req, float* x, float lambdaGF, float cfacGF, float ssGF, float bbGF, float ccGF,
	float* GF, float* GFx, float* GFy, float* GFz);
void tissrate(int nsp, float* p, float* mtiss, float* mptiss);

void bicgstabBLASDinit(int nnvGPU);
void bicgstabBLASDend(int nnvGPU);
void tissueGPUinit(int nntGPU, int nnvGPU);
void tissueGPUend(int nntGPU, int nnvGPU);

int max = 100, nseg, nnod, nsubseg, nt, nnodbc, niter, nnodfl, nsegfl, segfltot, imain, num_runs_actual, num_runs_actual_start, imainmax;
int mxx, myy, mzz, nntGPU, nre, nmaxvessel, nmaxtissue, nmax, nsp, rungreens, nodsegm, g0method, linmethod, nsprmax, sprmaxflag;
int nodnamemax, segnamemax, nnv_dim, nnt_dim, nnt_old_dim, nseg_dim, nnod_dim, nseg_new, nnod_new;
int nitmax1, nitmax, varyviscosity, phaseseparation, eliminate, nnt, nnv, varylb, outbounmethod;
int adaptd, inftrans, angio, tensionm, adjustks, is2d, angio3d;
int ncontourplanes, * slsegdiv, * nsl1, * nsl2;
int nsplus, nsspr, nspruned, nsminus;
int seed, actnum, actnum_dim, sprnum, nsprout, actupdateflag, runtyp, flag;
int dogreensflag;
int* activenode, * ista, * iend;
int* mainseg, * segtyp, * segname, * bcnodname, * bcnod, * bctyp, * nodname, * nodtyp, * lowflow, * nl, * indx;
int* ranknod, * nodrank, * nodfltyp, * nodout, * permsolute, * nspoint, * istart, * nk, * boundseg;
int* diffsolute, * oxygen, * imaxerrvessel, * imaxerrtissue, * nresis, * errvesselcount, * errtissuecount;
int* nodclasstyp, * segclasstyp;
int** nodnod, ** nodseg, ** nodflseg, ** segnodname, ** tisspoints;
int*** nbou, *** nbou_old;
int nsegprev;	//needed for cmgui, Feb. 2017

float c, cext, hext;
float pi1 = atan(1.0) * 4.0, fac = 1. / 4. / pi1, flowfac = 1.e6 / 60.;
float totalflow, maxerr, maxerr1, lowflowcrit;
float mcvcorr, diamthresh, vplas, facfp, tlength;
float dalphap, dalpha, constvisc, consthd, vdom;
float alx, aly, alz, vol, req, dcrit, maxl, clmax, lb, q0, errfac, qinsum;
float w2d, r2d, tol, qtol, hdtol, omega, optw, optlam;
float time1, eqtime, dalpha0, p50, fn, cs, dalphap0, tumorrad, alphab, lamg0;
float plow, phigh, clowfac, chighfac, pphighfac;//added January 2012
float smax_flow, smax_noflow, new_diam, GFmax, nmeta1;
float kp1, km1, kc1, ks1, ksnew, ranks1, L1, J01, tauref1, Qref1, pO2ref1, extlength1in, extlength1out, timestep1;
float thresh, threshn, thresh50, ks, kgGF, kgVAF, kv, ko, lambdatension, tensionrate, kmetdiam;
float min_length, tolerance_1, tolerance_2, increment, variance, Vmax, r0, theta0;
float offset, totalexist_length, hypox_fract, extraction;
float lambdaGF2, cfacGF2, ssGF2, bbGF2, ccGF2;
float lambdaGF3, cfacGF3, ssGF3, bbGF3, ccGF3;
float* axt, * ayt, * azt, * ds, * g0, * diff, * pmin, * pmax, * pmean, * pref, * g0fac, * qvfac, * g0facnew;
float* diam, * rseg, * ksseg, * q, * qq, * hd, * qqhd, * bcprfl, * bchd, * lseg, * tau, * qold, * hdold;
float* segpress, * cond, * condup, * conddown, * metsig, * stot;
float* bifpar, * cpar, * viscpar, * fahrpar, * adaptpar;
float* x, * y, * P1, * P2, * P3, * midpt, * cbar, * mtiss, * mptiss, * dqvsumdg0, * dqtsumdg0;
float* epsvessel, * epstissue, * eps, * errvessel, * errtissue, * pinit, * p, * dtmin, * sourcefac;
float* rhs, * rhstest, * g0old, * ptt, * ptpt, * qtsum, * qvsum, * qvtemp, * segvar, * nodvar, * segc;
float* stau, * spress, * uptrans, * downtrans;
float** xsl0, ** xsl1, ** xsl2, * clmin, * clint, * cl;
float** pvt, ** pvprev, ** qvprev, ** cv, ** dcdp, ** zv;
float** ptprev, ** ptv, ** gamma1, ** qcoeff1, ** cv0, ** conv0;
float** gvv, ** al, ** activegr;
float** scos, ** ax, ** cnode, ** bcp, ** resis, ** resisdiam;
float** qv, ** qt, ** qt_old, ** pv, ** pev, ** pt, ** pt_old;
float** tissparam, ** qvseg, ** pvseg, ** pevseg, ** gradp, ** signalvar;
float** hex_norm, ** hex_tang;	//needed for hexagonal periodic network
float*** dtt, *** psl;
double* nodpress, * rhsl, * matx;
double** mat, ** rhsg;

//Needed for GPU version
int useGPU, nnvGPU;  //useGPU can be 0, 1 (2, 3 or 4 - not working) - see SoluteParams.dat
double* h_x, * h_b, * h_a, * h_rs;
double* d_a, * d_res, * d_x, * d_b;
double* d_r, * d_rs, * d_v, * d_s, * d_t, * d_p, * d_er;

float* h_ats, * h_xts, * h_bts, * h_diagc;
float* h_rts, * h_pts, * h_apts;
float* d_ats, * d_xts, * d_bts, * d_diagc;
float* d_rts, * d_pts, * d_apts;

double* h_at, * h_xt, * h_bt;
double* h_rt, * h_pt, * h_apt;
double* d_at, * d_xt, * d_bt;
double* d_rt, * d_pt, * d_apt;

int* h_tisspoints, * d_tisspoints;
float* pt000, * qt000, * pv000, * qv000, * dtt000, * h_tissxyz, * h_vessxyz;
float* d_qt000, * d_pt000, * d_qv000, * d_pv000, * d_dtt000;
float* d_tissxyz, * d_vessxyz, * d_hex_norm, * h_hex_norm;

int retinaflag; //for retinal annulus
float IR, OR;	//for retinal annulus, mesh guidance
int hexperiodic, idomain;	//switch for hexagonal periodic network structure
float aphex, aphexp;	//side of hexagon in periodic structure, distance of midpoint from center

int main(int argc, char* argv[])
{
	int inod, inodbc, i, j, k, iseg, loop_control, nsp0, ilast10 = 0, itp, isp;
	float duration, nodpress0, q00, kschange;
	float dtmean, dtsd, totalpts, nGFhigh;
	float HFlast10 = 0., VLlast10 = 0., TFlast10 = 0., GFlast10 = 0.;
	float GF, GFx, GFy, GFz, GFhigh;
	clock_t tstart = clock(), tfinish;

	FILE* ofp, * ofp1, * ofp2, * ofp3;	//ShrinkAdjust,MasterLog,ArrayChanges,endnums 

#if defined(__unix__) //Create a Current subdirectory if it does not already exist.
	if (!fs::exists("Current")) fs::create_directory("Current");
#elif defined(_WIN32)
	DWORD ftyp = GetFileAttributesA("Current\\");
	if (ftyp != FILE_ATTRIBUTE_DIRECTORY) system("mkdir Current");
#endif

	ofp1 = fopen("MasterLog.out", "w");
	fprintf(ofp1, " imain runtime nseg   nnv  nnt   kmain convflag convflagv convflagt isp errvessel errtissue ...\n");
	fclose(ofp1);
	ofp2 = fopen("Arraychanges.out", "w");
	fprintf(ofp2, "List of array changes\n");
	fclose(ofp2);
	input();	//read data files
	nsp0 = nsp;
	setuparrays0();
	if (hexperiodic) setup_hexperiodic();
	evalGFinit();	//this evaluates coefficients needed in evalGF
	//******************************************************************
	float GFupdatefac = 1.0;	//values < 1 give time-delayed GF production
	//******************************************************************
	if (mzz == 1) is2d = 1;					//set is2d = 1 for 2d version of greens, 0 for 3d version
	else is2d = 0;
	if (num_runs_actual > 0) num_runs_actual_start = num_runs_actual;
	else num_runs_actual_start = 0;
	setuparrays1(nseg, nnod);
	for (iseg = 1; iseg <= nseg; iseg++) boundseg[iseg] = 0;		//for periodic network, initially no boundary crossings
	for (iseg = 1; iseg <= nseg; iseg++) ksseg[iseg] = ks1;	//initialize - for restarting, ksseg values should be stored and read
	analyzenet(1);							//analyzenet(0): analyze only flowing segments (type 4 and 5)
	sprmaxflag = 0;
	//************ option for changing tissue domain *****************
		//varylb = 0;	//overrides value from adaptparams
	//********************************* start main loop ****************
	nntGPU = mxx * myy * mzz;	//this is the maximum possible number of tissue points - needed for GPU version
	nnvGPU = 3000;	//start by assigning a good amount of space on GPU - may increase nnvGPU later
	if (useGPU) {
		bicgstabBLASDinit(nnvGPU);
		tissueGPUinit(nntGPU, nnvGPU);
	}
	setuparrays2(1, 1, 1);
	netstatsinit();
	for (imain = 1; imain <= imainmax; imain++) {
		if (seed == 1) srand(time(0));		//this has to be inside loop to allow restarting of runs
		else srand(seed);
		printf("\n*** Start main loop, imain = %i, time = %5.1f ***\n", imain, time1);
		printf("Starting analyzenet\n");
		analyzenet(1);						//analyze all segments so that picturenetwork and contourO2GF work
		printf("Starting pixnetwork_aft1analyznetmainloop\n");
		for (iseg = 1; iseg <= nseg; iseg++) segvar[iseg] = iseg;
		for (inod = 1; inod <= nnod; inod++) nodvar[inod] = inod;
		picturenetwork(nodvar, segvar, "Network_aft1anetmainloop.ps");
		flowtest(1);
		if (dogreensflag == 1) {	//remove disconnected type 1 segments
			printf("Starting segreclass");
			segreclass();
		}
		flag = 1;
		loop_control = 0;
		while (flag == 1) {
			flag = 0;
			analyzenet(0);
			flowtest(0);
			loop_control++;
			if (loop_control > (nseg + 150)) {
				printf("*** Warning: flowtest while loop exceeded loop control value(>nseg+150)\n");
				flag = 0;
			}
		}
		analyzenet(0);						//analyze only flowing segments (type 4 and 5) for flow calculation
		printf("Starting pixnetwork_bforeflow\n");
		for (iseg = 1; iseg <= nseg; iseg++) segvar[iseg] = iseg;
		for (inod = 1; inod <= nnod; inod++) nodvar[inod] = inod;
		picturenetwork(nodvar, segvar, "Network_bfflow.ps");
		printf("Starting flow\n");
		flow();
		qinsum = 0.;	//calculate total inflow to network
		for (inodbc = 1; inodbc <= nnodbc; inodbc++) {
			inod = bcnodname[inodbc];
			iseg = nodseg[1][inod];
			if (inod == nodname[ista[iseg]] && q[iseg] > 0.) qinsum += q[iseg];	//corrected Feb. 2017
			else if (inod == nodname[iend[iseg]] && q[iseg] < 0.) qinsum -= q[iseg];	//corrected Feb. 2017
		}
		printf("qinsum = %f\n", qinsum);
		dogreensflag = 0;	//don't run greens until vessels start flowing
		for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] == 5) dogreensflag = 1;
		if (imain == 1) {
			nodpress0 = nodpress[1];
			q00 = q[1];
			ofp = fopen("ShrinkAdjust.out", "w");
			fprintf(ofp, "Ptarget = %f\n", nodpress0);
			fprintf(ofp, "Qtarget = %f\n", q00);
			fprintf(ofp, "imain    time      Pin       Qin       ks       thresh     tlength\n");
			fclose(ofp);
		}
		else if (dogreensflag == 1 && imain > 1) {
			kschange = 0.;
			if (adjustks == 1) kschange = -(nodpress[1] - nodpress0) * 0.002;	//automatic adjustment of shrinking tendency based on pressure at main inflow
			if (adjustks == 2) kschange = (q[1] - q00) * 0.002;	//automatic adjustment of shrinking tendency based on pressure at main inflow, 0.002 changed for 6.0 consumption 6_2015JPA
			ks1 += kschange;
			if (ks1 <= 0.5) printf("*** Warning: low ks1 = %4.3f\n", ks1);
			for (iseg = 1; iseg <= nseg; iseg++) ksseg[iseg] += kschange;
		}
		ofp = fopen("ShrinkAdjust.out", "a");
		fprintf(ofp, "%4i  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f  %8.0f\n", num_runs_actual, time1, nodpress[1], q[1], ks1, thresh, tlength);
		fclose(ofp);
		printf("Starting setupgreens: ");
		setupgreens();						//this gives new value of nnv
		if (varylb) {
			lb = sqrt(2. * diff[1] * 30. / tissparam[1][1]);	//2D diffusion distance, assuming blood PO2 of 30 mmHg - adjusts according to
			printf("Outer bound distance of tisssue domain recalculated: lb = %f\n", lb);
		}
		nnt = outboun(outbounmethod);		//create array of tissue points inside domain, using method 1, 2 or 3: see outboun.cpp
		resizearrays2(nnv, nnt, actnum + nsprmax);
		analyzenet(1);
		write_network(-1);
		analyzenet(0);

		if (hexperiodic) remap_hexperiodic();

		if (dogreensflag == 1) {
			printf("Starting greens\n");
			nsp = 1;							//oxygen calculation only in greens
			greens();
			nsp = nsp0;							//restore value
			//This is need to give time-delayed GF production. March 2020.
			for (itp = 1; itp <= nnt; itp++) {
				ptpt[1] = pt[itp][1];
				tissrate(nsp, ptpt, mtiss, mptiss);	//GF production rate
				for (isp = 2; isp <= nsp; isp++) {
				if(imain == 1) qt[itp][isp] = mtiss[isp] * vol;
				else qt[itp][isp] += GFupdatefac * (mtiss[isp] * vol - qt[itp][isp]);
				}
			}
			printf("Starting histograms\n");
			histogram();
			printf("Starting contourO2GF\n");	//graphics of oxygen and growth factor
			contourO2GF(num_runs_actual);
			cmgui(num_runs_actual);				//make files for cmgui. Feb. 2017
			for (iseg = 1; iseg <= nseg; iseg++) segvar[iseg] = iseg;
			for (inod = 1; inod <= nnod; inod++) nodvar[inod] = inod;
			picturenetwork(nodvar, segvar, "Network_tissptsafternewflow.ps");
		}
		analyzenet(0);
		printf("Starting output\n");
		output();
		printf("Starting netstats\n");
		netstats(20.);							//argument is po2cutoff in mmHg
		printf("Starting write_network\n");
		analyzenet(1);							//added July 2013 because need for write-network that eliminates nodtyp =0 nodes
		write_network(-1);						//writes Networklast at each timestep
		//write_network(num_runs_actual);		//writes numbered network file in Current
		analyzenet(0);
		if (adaptd > 0 && dogreensflag == 1) {						//adaptation of diameters (adaptd = 1 or 2)
			printf("Starting conductconvect\n");
			conductconvect();
			printf("Starting adapt: segments killed");
			adapt();
			printf("\n");
			adaptsignals(num_runs_actual);		//make plot of signals vs. pressure. Feb. 2017
			printf("Starting analyzenet\n");	//analyze including all segments (except type 10)
			analyzenet(1);
			printf("Starting flowtest\n");
			flowtest(1);						//include all segments (except type 10)
			printf("Starting segreclass: segments killed");
			segreclass();						//remove disconnected type 1 segments
			printf("Starting pixnetwork_bforeflowtest\n");
			for (iseg = 1; iseg <= nseg; iseg++) segvar[iseg] = iseg;
			for (inod = 1; inod <= nnod; inod++) nodvar[inod] = inod;
			picturenetwork(nodvar, segvar, "Network_bfflowtest.ps");
		}
		if (angio == 1) {						//angiogenesis
			analyzenet(1);
			nseg_new = nseg + 2 * actnum + 3 * nsprmax;	//prepare arrays for new segments
			nnod_new = nnod + actnum + 2 * nsprmax;	//allow for intersections on two segments
			resizearrays1(nseg_new, nnod_new);

			//New program structure, February 2018
			printf("Starting analyzenet\n");
			analyzenet(1);
			if (actnum > 0) {
				printf("Starting actgrowth\n");
				actgrowth();
			}
			printf("Starting sprgrowth\n");
			sprgrowth();
			printf("Starting actsprupdate: intercepted segments");
			actsprupdate();
			printf("\n");
			printf("Starting analyzenet\n");
			analyzenet(1);
			printf("Starting prep_for_adapt: segments killed");
			for (inod = 1; inod <= nnod; inod++) prep_for_adapt(inod);
			printf("\n");
			printf("Starting analyzenet\n");
			analyzenet(1);						//analyzenet(1): include all segments except type 10
			printf("Starting flowtest\n");
			flowtest(1);						//see analyzenet above. remove disconnected type 1 segments
			if (dogreensflag == 1) {
				printf("Starting segreclass");
				segreclass();
			}
			printf("Starting straight_sht_seg_fixer: segments killed");
			straight_sht_seg_fixer();
		}
		if (tensionm == 1) {					//nodes moved due to vessel tension
			printf("Starting analyzenet\n");
			analyzenet(1);						//analyze including all segments (except type 10)
			printf("Starting tension: segments killed\n");
			tension();
			printf("\n");
			printf("Starting analyzenet\n");
			analyzenet(1);						//analyzenet(1): include all segments except type 10
			printf("Starting flowtest\n");
			flowtest(1);						//see analyzenet above
			if (dogreensflag == 1) {			//remove disconnected type 1 segments
				printf("Starting segreclass");
				segreclass();
			}
			printf("Starting pixnetwork_bforeflowtest\n");
			for (iseg = 1; iseg <= nseg; iseg++) segvar[iseg] = iseg;
			for (inod = 1; inod <= nnod; inod++) nodvar[inod] = inod;
			picturenetwork(nodvar, segvar, "Network_bfflowtest.ps");
		}
		analyzenet(1);
		write_network(-1);
		write_network(num_runs_actual);			//writes numbered network file in Current
		num_runs_actual++;
		seed += 7;
		if (imain == 101) timestep1 = 0.2;		//******************** increase time step ********************
		time1 += timestep1;
		for (i = 1; i <= mxx; i++) for (j = 1; j <= myy; j++) for (k = 1; k <= mzz; k++) nbou_old[i][j][k] = nbou[i][j][k];
		for (iseg = 1; iseg <= nseg; iseg++) segvar[iseg] = iseg;
		for (inod = 1; inod <= nnod; inod++) nodvar[inod] = inod;
		picturenetwork(nodvar, segvar, "Network_atend.ps");
		printf("*** End main loop, imain = %i, time = %5.1f ***\n", imain, time1);

		if (imain > imainmax - 10) {
			printf("Calculating growth factor levels at tissue points ...");
			nGFhigh = 0.;
			totalpts = 0.;
			for (itp = 1; itp <= nnt; itp++) {
				x[1] = axt[tisspoints[1][itp]];
				x[2] = ayt[tisspoints[2][itp]];
				x[3] = azt[tisspoints[3][itp]];
				evalGF(2, req, x, lambdaGF2, cfacGF2, ssGF2, bbGF2, ccGF2, &GF, &GFx, &GFy, &GFz);
				totalpts += sourcefac[itp];
				if (GF >= thresh) nGFhigh += sourcefac[itp];
			}
			GFhigh = nGFhigh * 100. / totalpts;
			ilast10++;
			HFlast10 += hypox_fract;
			VLlast10 += totalexist_length;
			TFlast10 += qinsum;
			GFlast10 += GFhigh;
			printf("done\n");
		}
	}
	HFlast10 = HFlast10 / ilast10;
	VLlast10 = VLlast10 / ilast10;
	TFlast10 = TFlast10 / ilast10;
	GFlast10 = GFlast10 / ilast10;
	//********************************* end main loop ****************
	printf("*** Completed main loop\n");
	tfinish = clock();
	duration = (float)(tfinish - tstart) / CLOCKS_PER_SEC;
	printf("%2.1f seconds for adaptation\n", duration);
	if (sprmaxflag > 0) printf("There were greater than max sprouts\n");
	fclose(ofp2);

	//statistics of distance from tissue points to nearest vessel
	dtmean = 0.;
	dtsd = 0.;
	totalpts = 0.;
	for (itp = 1; itp <= nnt; itp++) {
		totalpts += sourcefac[itp];
		dtmean += dtmin[itp] * sourcefac[itp];
		dtsd += SQR(dtmin[itp]) * sourcefac[itp];
	}
	dtmean = dtmean / totalpts;
	dtsd = sqrt(dtsd / totalpts - SQR(dtmean));

	ofp3 = fopen("endnums.out", "w");
	fprintf(ofp3, "%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n",
		HFlast10, VLlast10 / 1000., TFlast10, GFlast10, dtmean, dtsd, extraction);
	fprintf(ofp3, "Averages over last 10 time steps\n");
	fprintf(ofp3, "Hypoxic fraction = %f\n", HFlast10);
	fprintf(ofp3, "Vessel length = %f\n", VLlast10 / 1000.);
	fprintf(ofp3, "Total flow = %f\n", TFlast10);
	fprintf(ofp3, "HighGF fraction = %f\n", GFlast10);
	fprintf(ofp3, "Final values\n");
	fprintf(ofp3, "Hypoxic fraction = %f\n", hypox_fract);
	fprintf(ofp3, "Vessel length = %f\n", totalexist_length / 1000.);
	fprintf(ofp3, "Total flow = %f\n", qinsum);
	fprintf(ofp3, "HighGF fraction = %f\n", GFhigh);
	fprintf(ofp3, "Mean dist to vessel = %f\n", dtmean);
	fprintf(ofp3, "SD distance to vessel = %f\n", dtsd);
	fprintf(ofp3, "Extraction = %f\n", extraction);
	fclose(ofp3);

	printf("Averages over last 10 time steps\n");
	printf("Hypoxic fraction = %f\n", HFlast10);
	printf("Vessel length = %f\n", VLlast10 / 1000.);
	printf("Total flow = %f\n", TFlast10);
	printf("HighGF fraction = %f\n", GFlast10);
	printf("Final values\n");
	printf("Hypoxic fraction = %f\n", hypox_fract);
	printf("Vessel length = %f\n", totalexist_length / 1000.);
	printf("Total flow = %f\n", qinsum);
	printf("HighGF fraction = %f\n", GFhigh);
	printf("Mean dist to vessel = %f\n", dtmean);
	printf("SD distance to vessel = %f\n", dtsd);

	//close out cumulative postscript files
	ofp = fopen("OxygenStats.ps", "a");
	fprintf(ofp, "showpage\n");
	fclose(ofp);
	ofp = fopen("VesselTypes.ps", "a");
	fprintf(ofp, "showpage\n");
	fclose(ofp);

	if (useGPU) {
		tissueGPUend(nntGPU, nnvGPU);
		bicgstabBLASDend(nnvGPU);
	}

	return 0;
}
