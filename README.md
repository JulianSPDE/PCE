# PCE
1.) Supplementary figures and 
2.) MATLAB Code for the paper "Comparing intrusive and non-intrusive polynomial chaos for a class of exponential time differencing schemes"

1.)
The supplementary Figures and Tables are all contained in the file PCE_IJCM_Figures.pdf. They use the numbering from the paper, i.e. Equation (27) refers to the equation with that number in the main paper.


2.)
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
- PCE_supplementary_tables.pdf contains two tables showing the numbers of time steps and quadrature points used for all the simulations performed in the paper as well as a table showing the runtimes needed to produce the plots on our machine. The main takeaway from the runtime table is that the runtimes for iPCE and niPCE are roughly equal except for the cases with the cubic term, where iPCE starts to struggle with instabilities across the board and the number of time steps needs to increased substantially, and the case of D=! and the quadratic term for explicit Euler, where also the time step number has to be increased substantially due to instabilities.
