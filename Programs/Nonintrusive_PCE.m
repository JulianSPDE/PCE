% this program uses the program Nonlinear_Solver.m to do a Nonintrusive PCE siulation
% with certain input parameters. 
% Inputs:
% order: Number 
% p: Spatial resolution along one axis
% xL: Interval along one axis, e.g. [-1,1]
% T: Final time
% M: Number of time steps
% D_u, F, K: parameters in the nonlinear system
% varname: The name of the random variable as a string; 'D_u', 'F' or 'K'
% varmin, varmax: The range of the random variable (normal distribution will be
% outside with probability 1:1 million
% system: The number of the nonlinear system 
% polynomialtype: Either 'Legendre' or 'Hermite'.
% dimension: Spatial dimension of the domain, can be 1 or 2
% savefig: If 1, plots at 10 time points (10 hard-coded) are saved.

% The aim of this program is to compute the integral from Eldred, Burkardt
% (2009):
% We're interested in alpha_0, because we want to compute the mean.
% 1/<psi_0^2> * int_\Omega R(xi)*psi_0(xi)*rho(xi) de xi,
% Where xi follows a given random distribution.
% xi is sampled according to a quadrature rule according to the orthogonal
% polynomials belonging to the given random distribution.

% Input 'runs': Number of samples computed
% moment: Either 'mean' or 'variance'
% Mvec: Three entries, numbers of time steps for explicit Euler, ETDRDPIF
% and ETDRK4
function Nonintrusive_PCE(runs,moment,p,xL,T,Mvec,D_u,D_v,F,K,varname,varmin,varmax,system,polynomialtype,dimension,savefig,snapshots)
    
    xx = linspace(xL(1),xL(2),p);
    referenceruns = 200;  % Gauss Quadrature runs for reference
    
    spatialdiscrlist = ["FD","FD","Spectral"];
    timesteplist = ["EE","ETDRDPIF","ETDRK4"];
    
    runtimes = zeros(4,1);
    %% Computing "exact" or reference solution
    % (Note: For the ETDRK4 scheme here, we use a bigger number of time steps for the reference
    % solution, the number of time steps used for explicit Euler)
    if system~=6 || dimension~=1
        tic
        M_reference = 1000;
        MonteCarlo_Nonlinear_Solver('GaussQuad',moment,'Spectral','ETDRK4',referenceruns,p,xL,T,M_reference,D_u,D_v,F,K,varname,varmin,varmax,system,polynomialtype,dimension,savefig,snapshots);
        runtimes(1) = toc;
        MC_filename = sprintf(strcat('GaussQuad','_',moment,'_',num2str(dimension),'D_',polynomialtype,'_runs=%d_system=%d_p=%d_xL=%.2f..%.2f_T=%d_M=%d_Du=%.5f_Dv=%.5f_F=%.4f_K=%.4f_',varname,'=%.6f..%.6f.mat'),referenceruns,system,p,xL(1),xL(2),T,M_reference,D_u,D_v,F,K,varmin,varmax);
        fprintf("Opening:\n")
        disp(MC_filename)
        MC = importdata(MC_filename);
        u_total = MC.u_total;
        exactsolution_u = u_total; % exactsolution_v = v_total;
    end
    if system==6 && dimension==1
        tic
        xxtemp = linspace(xL(1),xL(2),p+1);
        xxtemp(end)=[];
        tt = linspace(0,T,snapshots+1);
        tt(1)=[];
        [TT,XX] = meshgrid(tt,xxtemp);
        if strcmp(moment,'mean')==1
            MC_u = 1/(varmax-varmin)*cos(pi*XX)./TT.*(exp(-(D_u*pi^2+varmin).*TT)-exp(-(D_u*pi^2+varmax).*TT));
        end
        if strcmp(moment,'variance')==1
            MC_u = cos(pi*XX).^2.*(1./(2*(varmax-varmin)*TT).*(exp(-2*(D_u*pi^2+varmin).*TT)-exp(-2*(D_u*pi^2+varmax).*TT)) -(1./((varmax-varmin)*TT).*(exp(-(D_u*pi^2+varmin).*TT)-exp(-(D_u*pi^2+varmax).*TT))).^2);
        end
        exactsolution_u = MC_u;
        u_total = MC_u;
        runtimes(1) = toc;
    end
    
    for scheme=1:3
    spatialdiscr = spatialdiscrlist(scheme);
    timestepping = timesteplist(scheme);

    %% Computing approximations of the mean
    
    % MC
    mcruns = 10;
    if dimension==1
        u_MC_total = zeros(p,snapshots,mcruns);
    end
    if dimension==2
        u_MC_total = zeros(p,p,snapshots,mcruns);
    end
    tic
    for j=1:mcruns
        MonteCarlo_Nonlinear_Solver('MC',moment,spatialdiscr,timestepping,runs,p,xL,T,Mvec(scheme),D_u,D_v,F,K,varname,varmin,varmax,system,polynomialtype,dimension,savefig,snapshots);
        MC_filename = sprintf(strcat('MC','_',moment,'_',num2str(dimension),'D_',polynomialtype,'_runs=%d_system=%d_p=%d_xL=%.2f..%.2f_T=%d_M=%d_Du=%.5f_Dv=%.5f_F=%.4f_K=%.4f_',varname,'=%.6f..%.6f.mat'),runs,system,p,xL(1),xL(2),T,Mvec(scheme),D_u,D_v,F,K,varmin,varmax);
        fprintf("Opening:\n")
        disp(MC_filename)
        MC = importdata(MC_filename);
        u_MC_current = MC.u_total;
        if dimension==1
            u_MC_total(:,:,j) = u_MC_current;
            figure(6)
            plot(xx,u_MC_total(:,1,j))
        end
        if dimension==2
            u_MC_total(:,:,:,j) = u_MC_current;
        end
    end
    total_runtime_MC = toc/mcruns;
    % Sobol
    tic
    MonteCarlo_Nonlinear_Solver('Sobol',moment,spatialdiscr,timestepping,runs,p,xL,T,Mvec(scheme),D_u,D_v,F,K,varname,varmin,varmax,system,polynomialtype,dimension,savefig,snapshots);
    MC_filename = sprintf(strcat('Sobol','_',moment,'_',num2str(dimension),'D_',polynomialtype,'_runs=%d_system=%d_p=%d_xL=%.2f..%.2f_T=%d_M=%d_Du=%.5f_Dv=%.5f_F=%.4f_K=%.4f_',varname,'=%.6f..%.6f.mat'),runs,system,p,xL(1),xL(2),T,Mvec(scheme),D_u,D_v,F,K,varmin,varmax);
    fprintf("Opening:\n")
    disp(MC_filename)
    MC = importdata(MC_filename);
    u_Sobol_total = MC.u_total;
    total_runtime_Sobol = toc;
    % Gauss Quadrature
    tic
    MonteCarlo_Nonlinear_Solver('GaussQuad',moment,spatialdiscr,timestepping,runs,p,xL,T,Mvec(scheme),D_u,D_v,F,K,varname,varmin,varmax,system,polynomialtype,dimension,savefig,snapshots);
    MC_filename = sprintf(strcat('GaussQuad','_',moment,'_',num2str(dimension),'D_',polynomialtype,'_runs=%d_system=%d_p=%d_xL=%.2f..%.2f_T=%d_M=%d_Du=%.5f_Dv=%.5f_F=%.4f_K=%.4f_',varname,'=%.6f..%.6f.mat'),runs,system,p,xL(1),xL(2),T,Mvec(scheme),D_u,D_v,F,K,varmin,varmax);
    fprintf("Opening:\n")
    disp(MC_filename)
    MC = importdata(MC_filename);
    u_GQ_total = MC.u_total;
    total_runtime_GQ = toc;
        %{
        if j==1
            exactsolution_u = u_total; % exactsolution_v = v_total;
        end
        %}
    errorarray_MC = zeros(1,snapshots);
    figure(1)
    plot(xx,u_MC_total(:,10))
    disp(size(exactsolution_u))
    figure(2)
    plot(xx,exactsolution_u(:,10))
    figure(3)
    for i=1:mcruns
        if dimension==1
            errorarray_MC = errorarray_MC + sqrt(trapz(xx',abs(u_MC_total(:,:,i)-exactsolution_u).^2,1));
        end
        if dimension==2
            errorarray_MC = errorarray_MC + sqrt(trapz(xx',abs(u_MC_total(:,:,:,i)-exactsolution_u).^2,1));
        end
    end
    errorarray_MC = errorarray_MC/mcruns;
    if dimension==1
        errorarray_Sobol = sqrt(trapz(xx',abs(u_Sobol_total-exactsolution_u).^2,1));
        errorarray_GQ = sqrt(trapz(xx',abs(u_GQ_total-exactsolution_u).^2,1));
    end
    if dimension==2
        errorarray_Sobol = sqrt(trapz(xx',trapz(xx',abs(u_Sobol_total-exactsolution_u).^2,2),1));
        errorarray_GQ = sqrt(trapz(xx',trapz(xx',abs(u_GQ_total-exactsolution_u).^2,2),1));
    end
    if dimension==1
        errorarray_MC = errorarray_MC./sqrt(trapz(xx',exactsolution_u.^2,1));
        errorarray_Sobol = errorarray_Sobol./sqrt(trapz(xx',exactsolution_u.^2,1));
        errorarray_GQ = errorarray_GQ./sqrt(trapz(xx',exactsolution_u.^2,1));
    end
    if dimension==2
        errorarray_MC = errorarray_MC./sqrt(trapz(xx',trapz(xx',exactsolution_u.^2,2),1));
        errorarray_Sobol = errorarray_Sobol./sqrt(trapz(xx',trapz(xx',exactsolution_u.^2,2),1));
        errorarray_GQ = errorarray_GQ./sqrt(trapz(xx',trapz(xx',exactsolution_u.^2,2),1));
    end
    


    %% Plotting solutions (optional)
    D_u_scalar = D_u;
    D_v_scalar = D_v;
    F_scalar = F;
    K_scalar = K;
    filename = sprintf(strcat('Nonintrusive_PCE_GaussQuad_',moment,'_',num2str(dimension),'D_',polynomialtype,'_order=%d_system=%d_p=%d_xL=%.2f..%.2f_T=%d_M=%d_Du=%.5f_Dv=%.5f_F=%.4f_K=%.4f_',varname,'=%.6f..%.6f.mat'),runs,system,p,xL(1),xL(2),T,Mvec(scheme),D_u_scalar,D_v_scalar,F_scalar,K_scalar,varmin,varmax);
    %filename = sprintf(strcat(simtype,'_2D_V_%03d_',polynomialtype,'_t=%.4f_p=%d_xL=%.2f_T=%d_M=%d_Du=%.5f_F=%.4f_K=%.4f_',varname,'=%.6f..%.6f','system=%d.png'),fignum,t,p,xL,T,M,D_u,F,K,varmin,varmax,system);
    fprintf("Saved as:\n")
    disp(filename)
    save(filename,'u_total')
    fignum = 0;
    if savefig==1
        if dimension==1
            for i=1:snapshots
                t = T/snapshots*i;
                plot_1D('GaussQuad',moment,spatialdiscr,timestepping,u_total(:,i),u_total(:,i),fignum,polynomialtype,varname,t,p,xL,T,Mvec(scheme),0,D_u_scalar,D_v_scalar,F_scalar,K_scalar,varmin,varmax,system)
                fignum = fignum + 1;
            end
        end
        if dimension==2
            for i=1:snapshots
                t = T/snapshots*i;
                plot_2D('GaussQuad',spatialdiscr,timestepping,u_total(:,:,i),u_total(:,:,i),fignum,polynomialtype,varname,t,p,xL,T,Mvec(scheme),0,D_u_scalar,D_v_scalar,F_scalar,K_scalar,varmin,varmax,system,1)
                fignum = fignum + 1;
            end
        end
    end
        
    %% Plotting time-dependent error and saving it
    figure(1)
    clf;
    hold on
    xx_time = linspace(0,T,snapshots);
    plot(xx_time,errorarray_MC) % We only want the biggest order, we compare MC, QMC and GQ here
    plot(xx_time,errorarray_Sobol)
    plot(xx_time,errorarray_GQ)
    set(gca,'Yscale','log')
    hold off
    if savefig==1
        plotname = sprintf(strcat('Errorplot_Nonintrusive_PCE_',moment,'_',spatialdiscr,'_',timestepping,'_','system=%d_p=%d_xL=%.2f..%.2f_T=%d_M=%d_Du=%.5f_Dv=%.5f_F=%.4f_K=%.4f_',varname,'=%.6f..%.6f.png'),system,p,xL(1),xL(2),T,Mvec(scheme),D_u,D_v,F,K,varmin,varmax);
        saveas(gca,plotname)
    end
    clf
    errorarray = [errorarray_MC;errorarray_Sobol;errorarray_GQ];
    if strcmp(timestepping,'EE')==1
        errorarray_FD_EE = [xx_time',errorarray'];
        writematrix(errorarray_FD_EE,sprintf(strcat('errorarray_',moment,'_','Nonintrusive_FD_EE_system=%d_D=%.5f.txt'),system,D_u));
    end
    if strcmp(timestepping,'ETDRDPIF')==1
        errorarray_FD_ETDRDP = [xx_time',errorarray'];
        writematrix(errorarray_FD_ETDRDP,sprintf(strcat('errorarray_',moment,'_','Nonintrusive_FD_ETDRDP_system=%d_D=%.5f.txt'),system,D_u));
    end
    if strcmp(timestepping,'ETDRK4')==1
        errorarray_Spectral = [xx_time',errorarray'];
        writematrix(errorarray_Spectral,sprintf(strcat('errorarray_',moment,'_','Nonintrusive_Spectral_system=%d_D=%.5f.txt'),system,D_u));
    end
    total_runtime_average = (total_runtime_MC + total_runtime_Sobol + total_runtime_GQ)/3;
    
    runtimes(scheme+1) = total_runtime_average;
    %writematrix(errorarray_FD_EE,sprintf(strcat('Nonintrusive_errorarray_',spatialdiscr,'_',timestepping,'_system=%d_D=%.5f.txt'),system,D_u));
    end

    save(sprintf(strcat('Nonintrusive_',moment,'_runtimes_system=%d_D=%.5f.mat'),system,D_u),'runtimes');

end