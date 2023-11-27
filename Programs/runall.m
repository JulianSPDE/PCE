
matlab2tikz_directory = strcat(pwd,'\matlab2tikz');
addpath(matlab2tikz_directory);

removefiles = 1; % If set to 1, all generated .mat, .txt, .tex and .png files are deleted from the folder in the end
printtikz = 1; % If set to 1, Tikz plots for Gray-Scott simulations will be generated in the folder
snaps = 50;  % Number of error snapshots in time
interval = [-1,1];  % Spatial interval
p = 128;   % Spatial resolution
q = 200;   % Number of quadrature points
N = [1,2,3,4,5];   % Vector with tested values for number of PCE coefficients
savefig = 0;    % 1 = save plots as png, 0 = don't save plots
% The different Mvecs contain the number of time steps for explicit Euler,
% ETDRDP-IF and ETDRK4, respectively.
%% 1D Intrusive PCE
% Without diffusion
% Linear, mean
%PCE_time_errorplot(snaps,'mean',2,interval,p,[1000,400,200],1000,q,N,0,0,6,1,'Legendre',0,0,'F',1,2,savefig,1);
PCE_time_errorplot(snaps,'mean',2,interval,p,[1000,200,100],1000,q,N,0,0,6,1,'Legendre',0,0,'F',1,2,savefig,1);
% Linear, variance
%PCE_time_errorplot(snaps,'variance',2,interval,p,[2000,400,200],1000,q,N,0,0,6,1,'Legendre',0,0,'F',1,2,savefig,1);
PCE_time_errorplot(snaps,'variance',2,interval,p,[2000,200,100],1000,q,N,0,0,6,1,'Legendre',0,0,'F',1,2,savefig,1);
% Quadratic
%PCE_time_errorplot(snaps,'mean',1,interval,p,[1000,200,100],1000,q,N,0,0,7,1,'Legendre',0,0,'F',1,2,savefig,1);
PCE_time_errorplot(snaps,'mean',1,interval,p,[1000,200,100],1000,q,N,0,0,7,1,'Legendre',0,0,'F',1,2,savefig,1);
% Cubic
%PCE_time_errorplot(snaps,'mean',10,interval,p,[1000,200,100],1000,q,N,0,0,8,1,'Legendre',0,0,'F',1,2,savefig,1);
PCE_time_errorplot(snaps,'mean',2,interval,p,[1000,200,100],1000,q,N,0,0,8,1,'Legendre',0,0,'F',1,2,savefig,1);

% With diffusion
% Linear
%PCE_time_errorplot(snaps,'mean',2,interval,p,[20000,5000,2000],20000,q,N,1,0,6,1,'Legendre',0,0,'F',1,2,savefig,1);
PCE_time_errorplot(snaps,'mean',2,interval,p,[20000,400,200],1000,q,N,1,0,6,1,'Legendre',0,0,'F',1,2,savefig,1);
% Quadratic
PCE_time_errorplot(snaps,'mean',2,interval,p,[20000,400,200],1000,q,N,1,0,7,1,'Legendre',0,0,'F',1,2,savefig,1);
% Cubic
%PCE_time_errorplot(snaps,'mean',2,interval,p,[20000,5000,2000],20000,q,N,1,0,8,1,'Legendre',0,0,'F',1,2,savefig,1);
PCE_time_errorplot(snaps,'mean',2,interval,p,[20000,400,200],1000,q,N,1,0,8,1,'Legendre',0,0,'F',1,2,savefig,1);
%% 1D Non-intrusive PCE
runs = 50;
% Without diffusion
% Linear, mean
Nonintrusive_PCE(runs,'mean',p,interval,2,[2000,200,100],0,0,0,0,'F',1,2,6,'Legendre',1,savefig,snaps)
Nonintrusive_PCE(runs,'variance',p,interval,2,[2000,200,100],0,0,0,0,'F',1,2,6,'Legendre',1,savefig,snaps)
% Quadratic
Nonintrusive_PCE(runs,'mean',p,interval,2,[2000,200,100],0,0,0,0,'F',1,2,7,'Legendre',1,savefig,snaps)
% Cubic
Nonintrusive_PCE(runs,'mean',p,interval,2,[500,200,100],0,0,0,0,'F',1,2,8,'Legendre',1,savefig,snaps)
% With diffusion
% Linear
Nonintrusive_PCE(runs,'mean',p,interval,2,[20000,200,100],1,0,0,0,'F',1,2,6,'Legendre',1,savefig,snaps)
% Quadratic
Nonintrusive_PCE(runs,'mean',p,interval,2,[20000,200,100],1,0,0,0,'F',1,2,7,'Legendre',1,savefig,snaps)
% Cubic
Nonintrusive_PCE(runs,'mean',p,interval,2,[20000,200,100],1,0,0,0,'F',1,2,8,'Legendre',1,savefig,snaps)


%% 2D Intrusive PCE

% Use this to retrieve the ETDRDPIF and ETDRK4 results (resolution p=512)
PCE_time_errorplot(snaps,'mean',1,interval,512,[20,20,20],1000,q,N,1,0,6,2,'Legendre',1,0,'F',1,2,savefig,1);
% Use this simulation to retrieve the explicit Euler results
PCE_time_errorplot(snaps,'mean',1,interval,128,[10000,50,50],1000,q,N,1,0,6,2,'Legendre',1,0,'F',1,2,savefig,1);
%PCE_time_errorplot(snapshots,moment,T,xL,pvec,Mvec,M_MC,runs,Nvec,D_u,D_v,system,dimension,polynomialtype,F,K,varname,varmin,varmax,savefig,relative)


%% Gray-Scott examples
% Error plot: Demonstration of failure of schemes
PCE_time_errorplot(snaps,'mean',500,interval,p,[5000,1000,500],1000,q,N,0.00002,0.00001,0,1,'Legendre',0,0.06,'F',0.04,0.05,savefig,1);
% plots in the paper


% Plots showing sensitivity in long-time integration
% If printtikz is activated, .tex files with the plots are generated in the
% folder
plot_snaps = 1;
Nonlinear_Solver('FD','EE',256,[-1,1],5000,100000,0.00002,0.00001,0.04,0.06,0,2,1,1,plot_snaps)
if printtikz
    matlab2tikz('EE_2.tex');
end
Nonlinear_Solver('FD','EE',1024,[-1,1],5000,100000,0.00002,0.00001,0.04,0.06,0,2,1,1,plot_snaps)
if printtikz
    matlab2tikz('EE_3.tex');
end
Nonlinear_Solver('FD','ETDRDPIF',256,[-1,1],5000,10000,0.00002,0.00001,0.04,0.06,0,2,1,1,plot_snaps);
if printtikz
    matlab2tikz('ETDRDP_2.tex');
end
Nonlinear_Solver('FD','ETDRDPIF',1024,[-1,1],5000,10000,0.00002,0.00001,0.04,0.06,0,2,1,1,plot_snaps);
if printtikz
    matlab2tikz('ETDRDP_4.tex');
end
Nonlinear_Solver('Spectral','ETDRK4',256,[-1,1],5000,10000,0.00002,0.00001,0.04,0.06,0,2,1,1,plot_snaps);
if printtikz
    matlab2tikz('Spectral_1.tex');
end
Nonlinear_Solver('Spectral','ETDRK4',1024,[-1,1],5000,10000,0.00002,0.00001,0.04,0.06,0,2,1,1,plot_snaps);
if printtikz
    matlab2tikz('Spectral_2.tex');
end


MonteCarlo_Nonlinear_Solver('GaussQuad','mean','Spectral','ETDRK4',20,512,[-1,1],1500,1500,0.00002,0.00001,0.04,0,'K',0.058,0.062,0,'Legendre',2,1,1)
if printtikz
    matlab2tikz('Spectral_niPCE_1.tex');
end
N_GS = 2;
PCE_solver('Spectral','ETDRK4','mean',256,[-1,1],200,20000,N_GS,0.00002,0.00001,0.04,0,'K',0.058,0.062,0,'Legendre',2,1,1,50);
if printtikz
    matlab2tikz('Spectral_iPCE_1.1.tex');
end
PCE_solver('FD','EE','mean',512,[-1,1],400,50000,N_GS,0.00002,0.00001,0.04,0,'K',0.058,0.062,0,'Legendre',2,1,1,50);
if printtikz
    matlab2tikz('EE_iPCE_1.tex');
end


%% Performance Plots

Performance_Plot('mean',2,[-1,1],128,[5000,2000,1000,500,200,100,50],20000,100,20,5,0,0,6,1,'Legendre',0,0,'F',1,2,0,1)
Performance_Plot('mean',2,[-1,1],128,[5000,2000,1000,500,200,100,50],20000,100,20,5,1,0,6,1,'Legendre',1,0,'F',1,2,0,1)
Performance_Plot('mean',0.4,[-1,1],128,[5000,2000,1000,500,200,100,50],20000,100,20,5,0,0,7,1,'Legendre',0,0,'F',1,2,0,1)
Performance_Plot('mean',2,[-1,1],128,[5000,2000,1000,500,200,100,50],20000,100,20,5,1,0,7,1,'Legendre',1,0,'F',1,2,0,1)
Performance_Plot('mean',2,[-1,1],128,[5000,2000,1000,500,200,100,50],20000,100,20,5,0,0,8,1,'Legendre',0,0,'F',1,2,0,1)
Performance_Plot('mean',2,[-1,1],128,[5000,2000,1000,500,200,100,50],20000,100,20,5,1,0,8,1,'Legendre',1,0,'F',1,2,0,1)

%% Run time test, showing iPCE runtimes depending on N
PCE_runtimetest;

%% Runtimes 1D intrusive PCE

r_6_1 = importdata('runtimes_mean_system=6_D=0.00000.mat');
r_6_2 = importdata('runtimes_mean_system=6_D=1.00000.mat');
r_7_1 = importdata('runtimes_mean_system=7_D=0.00000.mat');
r_7_2 = importdata('runtimes_mean_system=7_D=1.00000.mat');
r_8_1 = importdata('runtimes_mean_system=8_D=0.00000.mat');
r_8_2 = importdata('runtimes_mean_system=8_D=1.00000.mat');
fprintf("Runtimes intrusive PCE for N=%d:\n",N(end))
fprintf("System 6, D=0 / D=1:\n")
fprintf("EE: %.4f seconds, ETDRDP: %.4f seconds, ETDRK4: %.4f seconds\n",r_6_1(1),r_6_1(2),r_6_1(3))
fprintf("EE: %.4f seconds, ETDRDP: %.4f seconds, ETDRK4: %.4f seconds\n",r_6_2(1),r_6_2(2),r_6_2(3))
fprintf("System 7, D=0 / D=1:\n")
fprintf("EE: %.4f seconds, ETDRDP: %.4f seconds, ETDRK4: %.4f seconds\n",r_7_1(1),r_7_1(2),r_7_1(3))
fprintf("EE: %.4f seconds, ETDRDP: %.4f seconds, ETDRK4: %.4f seconds\n",r_7_2(1),r_7_2(2),r_7_2(3))
fprintf("System 8, D=0 / D=1:\n")
fprintf("EE: %.4f seconds, ETDRDP: %.4f seconds, ETDRK4: %.4f seconds\n",r_8_1(1),r_8_1(2),r_8_1(3))
fprintf("EE: %.4f seconds, ETDRDP: %.4f seconds, ETDRK4: %.4f seconds\n",r_8_2(1),r_8_2(2),r_8_2(3))


%% Runtimes 1D non-intrusive PCE
rn_6_1 = importdata('Nonintrusive_mean_runtimes_system=6_D=0.00000.mat');
rn_6_2 = importdata('Nonintrusive_mean_runtimes_system=6_D=1.00000.mat');
rn_7_1 = importdata('Nonintrusive_mean_runtimes_system=7_D=0.00000.mat');
rn_7_2 = importdata('Nonintrusive_mean_runtimes_system=7_D=1.00000.mat');
rn_8_1 = importdata('Nonintrusive_mean_runtimes_system=8_D=0.00000.mat');
rn_8_2 = importdata('Nonintrusive_mean_runtimes_system=8_D=1.00000.mat');
fprintf("Runtimes non-intrusive PCE (for either of MC, QMC, GQ) with %d realizations:\n",runs)
fprintf("System 6, D=0 / D=1:\n")
fprintf("EE: %.4f seconds, ETDRDP: %.4f seconds, ETDRK4: %.4f seconds\n",rn_6_1(2),rn_6_1(3),rn_6_1(4))
fprintf("EE: %.4f seconds, ETDRDP: %.4f seconds, ETDRK4: %.4f seconds\n",rn_6_2(2),rn_6_2(3),rn_6_2(4))
fprintf("System 7, D=0 / D=1:\n")
fprintf("EE: %.4f seconds, ETDRDP: %.4f seconds, ETDRK4: %.4f seconds\n",rn_7_1(2),rn_7_1(3),rn_7_1(4))
fprintf("EE: %.4f seconds, ETDRDP: %.4f seconds, ETDRK4: %.4f seconds\n",rn_7_2(2),rn_7_2(3),rn_7_2(4))
fprintf("System 8, D=0 / D=1:\n")
fprintf("EE: %.4f seconds, ETDRDP: %.4f seconds, ETDRK4: %.4f seconds\n",rn_8_1(2),rn_8_1(3),rn_8_1(4))
fprintf("EE: %.4f seconds, ETDRDP: %.4f seconds, ETDRK4: %.4f seconds\n",rn_8_2(2),rn_8_2(3),rn_8_2(4))

%% Deletion of created files (commented out, can be activated)
deletion_matfiles = "rm *.mat";
deletion_pngfiles = "rm *.png";
deletion_txtfiles = "rm *.txt";

if removefiles
    system(deletion_matfiles);
    system(deletion_pngfiles);
    system(deletion_txtfiles);
end

