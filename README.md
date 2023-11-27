# PCE
MATLAB Code for the paper "Comparing intrusive and non-intrusive polynomial chaos for a class of exponential time differencing schemes"
This repository contains all the code necessary to reproduce the plots in the paper.
- All plots can be produced by running runall.m, and .txt files titled "errorarray...txt" will be generated, where the points contain specifics on the parameters used.
- Please note that running the entire runall.m takes a considerable amount of time. Especially some of the 2D iPCE simulations may take hours of computation time. It is therefore advisable to pick plots to reproduce and simply execute those lines of code you are interested in. (Be sure to initialize all the variables at the top of runall.m if you do this).
- The file PCE_diagram.pdf contains a chart sketching how the programs are connected, i.e. which programs call which other programs.
- The main programs which generate the plots are:
	1.) PCE_time_errorplot.m for iPCE simulations and their time-dependent error
	2.) Nonintrusive_PCE.m for niPCE simulations and their time-dependent error
	3.) PCE_solver.m for a single iPCE simulation
	4.) MonteCarlo_Nonlinear_Solver.m to generate a single MC, QMC or GQ niPCE simulation for mean or variance
	5.) Performance_Plot.m for iPCE and niPCE simulations, comparing to a reference solution and plotting the errors depending on M.

