# AngioAdapt20: Simulation of angiogenesis, remodeling and pruning in three dimensions

See Secomb, T.W., Alberding, J.P., Hsu, R., Dewhirst, M.W. and Pries, A.R. Angiogenesis: an adaptive biological patterning problem. PLoS Computational Biology 9:e1002983. doi:10.1371/journal.pcbi.1002983 (2013)  

See also manuscript "Simulation of angiogenesis in three dimensions: application to cerebral cortex" (in submission)

See https://github.com/secomb/FlowEstimateV1 for more details on calculation of blood flow in the network.

See https://github.com/secomb/GreensV4 for more details on Green's function method to calculate spatial distributions of oxygen and growth factor in tissue.

**Input files**

**Network.dat** This file specifies the network structure, the flow rate and hematocrit of all segments and the PO2 and solute concentrations of boundary nodes. (Note that these solute concentrations are used only for inflow nodes.)  
**SoluteParams.dat** This file gives solute transport and reaction parameters. Oxygen is handled differently than other solutes to account for binding by hemoglobin. Other solutes are assumed to be carried in solution and not to bind to blood components. This file contains oxygen transport parameters, which are needed only if the “oxygen” parameter is set to 1.  
**ContourParams.dat** At each time point, the program generates a postscript file showing the vessel network projected onto a single plane and the solute contours on that plane. The position and orientation of the plane are specified by ContourParams.dat, giving the coordinates of three corners of the plane. Also, the contour levels are specified.  
**IntravascRes.dat** Specifies intravascular resistance to solute transport from blood to tissue. For oxygen, this depends on the vessel diameter. For other solutes, a wall permeability is specified by the file.  
**RheolParams.dat** This file gives parameters for equations to describe the apparent viscosity of blood in microvessels in vivo as a function of vessel diameter and discharge hematocrit, and the partition of hematocrit in diverging bifurcations, as a function of the flow rates in the branches, the diameters of the branches, and the hematocrit arriving at the parent vessel. Set varyviscosity = 1 to get diameter-dependent viscosity. Set phaseseparation = 1 to compute phase separation in diverging bifurcations. This can lead to noncTonvergence of the method in some cases. Also, the number of iterations for linear and nonlinear loops are specified. These may need to be increased for very large networks.  
**AngioParams.dat**  This file gives the values of parameters associated with angiogenesis and vessel migration, including:
Time step  
Diameter of new sprouts  
Threshold GF concentration for sprouting  
Constant in sprouting probability function  
Maximum sprout formation probability  
Sprout growth rate  
Directional response to GF gradient  
Attraction constant to nearby vessels  
Maximum vessel sensing distance  
Maximum vessel sensing angle  
Variance of growth direction randomization  
Vessel migration parameters  
Threshold for migration  
Maximum migration velocity  
**AdaptParams.dat** This file gives the values of parameters associated with structural adaptation and pruning, including:  
Structural adaptation time scale  
Reference wall shear stress  
Metabolic sensitivity  
Shrinking tendency  
Vessel permeability to GF (or GF product)  
Reference flow rate for metabolic signal  
Convected response saturation constant  
Conducted response saturation constant  
Conducted response length constant  
Relative strength of conducted response 

**Output files**

The program creates a new folder called "Current" if it does not already exist. As the simulation proceeds, files are placed in this folder as follows, where xxx is the current time step number:  
**Networkxxx.dat** The current network file.  
**network.exelem** and **etworkxxx.exnode** These files are used to obtain 3D visualizations of the network using CMGUI. When CMGUI is started, use File - Open - com file, select greens.com.txt and hit “All” for the visualization. For details of cmgui, see: http://sourceforge.net/projects/cmiss/files/cmgui/cmgui-wx-2.8.0/  
**AdaptSignalsxxx.ps** This file gives a graphical representation of the current signals for diameter change, as a function of intrvascular pressure. The color code is: red - shear stress; green - pressure; blue - shringking tendency; purple - upstream conducted response; olive - downstream convected response; black - total signal.  
**ContourO2GF000.ps** This file shows a map of the current network projected onto the specified contour plane. Oxygen levels are shown by colored contours on the plane, and growth factor levels above threshold are indicated by diagonal hatching.  

At the end of the simulation, multiple numerical and graphical files are produced to show the final state of the network, such as histograms of vessel flow rate, flow velocity and shear stress, and vessel and tissue oxygen levels.

This code is free to use at your own risk. Feedback and/or acknowledgement are appreciated. Email secomb@u.arizona.edu.  

Note: This code makes use of nrutil.h and nrutil.cpp as placed in the public domain by Numerical Recipes at http://apps.nrbook.com/c/index.html.

Updated 12 January 2021