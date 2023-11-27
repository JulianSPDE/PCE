% This program takes sets of input parameters and generates PCE and MC
% simulations for those parameters. Then, time-dependent error plots are
% generated, comparing each PCE simulation to the corresponding MC
% simulation and howing the evolution of the error over time.
% For each input, the simulations are done both for FD+EE and
% Spectral+ETDRK4.
% Last input: relative. For 1, relative error is plotted. For 0, absolute
% error is plotted
function PCE_time_errorplot(snapshots,moment,T,xL,pvec,Mvec,M_MC,runs,Nvec,D_u,D_v,system,dimension,polynomialtype,F,K,varname,varmin,varmax,savefig,relative)
figure(1)
clf;
xvec = linspace(0,T,snapshots);
plen = length(pvec); Nlen = length(Nvec);
type = "GaussQuad";
errorarray_Spectral = zeros(snapshots,Nlen);
errorarray_FD_EE = zeros(snapshots,Nlen);
errorarray_FD_ETDRDP = zeros(snapshots,Nlen);

runtimes = zeros(3,1);
% contains the runtimes for intrusive PCE for the maximum entered N, for
% EE, ETDRDP and ETDRK4.

%% Doing simulations

%spectral_on = 1; % Can switch off the spectral method simulations here

for i1=1:plen
		for i3=1:Nlen
			p = pvec(i1); N=Nvec(i3);
            if i3 == Nlen
                tic
            end
            PCE_solver('FD','EE',moment,p,xL,T,Mvec(1),N,D_u,D_v,F,K,varname,varmin,varmax,system,polynomialtype,dimension,1,savefig,snapshots);
            if i3 == Nlen
                runtimes(1) = toc;
                tic
            end
            PCE_solver('FD','ETDRDPIF',moment,p,xL,T,Mvec(2),N,D_u,D_v,F,K,varname,varmin,varmax,system,polynomialtype,dimension,1,savefig,snapshots);
            if i3 == Nlen
                runtimes(2) = toc;
                tic
            end
            PCE_solver('Spectral','ETDRK4',moment,p,xL,T,Mvec(3),N,D_u,D_v,F,K,varname,varmin,varmax,system,polynomialtype,dimension,1,savefig,snapshots);
            if i3 == Nlen
                runtimes(3) = toc;
            end
		end
    
    if system~=6
        MonteCarlo_Nonlinear_Solver(type,moment,'Spectral','ETDRK4',runs,pvec(i1),xL,T,M_MC,D_u,D_v,F,K,varname,varmin,varmax,system,polynomialtype,dimension,savefig,snapshots);
        %MonteCarlo_Nonlinear_Solver(type,moment,'FD','ETDRDPIF',runs,pvec(i1),xL,T,M_MC,D_u,D_v,F,K,varname,varmin,varmax,system,polynomialtype,dimension,savefig,snapshots)
    end
    
end

%% Plotting errors

legendname = [];
for i=1:Nlen
    Nstring = sprintf("N=%d",Nvec(i));
    legendname = [legendname,Nstring];
end
for i1=1:plen
        p = pvec(i1);
        spatialdiscr = 'Spectral'; timestepping = 'ETDRK4';
        %spatialdiscr = 'FD'; timestepping = 'ETDRDPIF';
        if system~=6 % Reference solution
            MC_filename = sprintf(strcat(type,'_',moment,'_',num2str(dimension),'D_',polynomialtype,'_runs=%d_system=%d_p=%d_xL=%.2f..%.2f_T=%d_M=%d_Du=%.5f_Dv=%.5f_F=%.4f_K=%.4f_',varname,'=%.6f..%.6f.mat'),runs,system,p,xL(1),xL(2),T,M_MC,D_u,D_v,F,K,varmin,varmax);
            fprintf("Opening:\n")
            disp(MC_filename)
            MC = importdata(MC_filename);
            MC_u = MC.u_total;
        end
        %% ETDRK4
        figure(11)
        hold on
        for i3=1:Nlen
            N = Nvec(i3);
            PCE_filename = sprintf(strcat('PCE_',spatialdiscr,'_',timestepping,'_',num2str(dimension),'D_',polynomialtype,'_system=%d_p=%d_xL=%.2f..%.2f_T=%d_M=%d_N=%d_Du=%.5f_Dv=%.5f_F=%.4f_K=%.4f_',varname,'=%.6f..%.6f.mat'),system,p,xL(1),xL(2),T,Mvec(3),N,D_u,D_v,F,K,varmin,varmax);
            fprintf("Opening:\n")
            disp(PCE_filename)
            PCE = importdata(PCE_filename);  
            PCE_u = PCE.U_array;
            disp(snapshots)
            disp(size(PCE_u))
            if system==6
                if dimension==1
                    xx = linspace(xL(1),xL(2),p+1);
                    xx(end)=[];
                    tt = linspace(0,T,snapshots+1);
                    tt(1)=[];
                    [TT,XX] = meshgrid(tt,xx);
                    if strcmp(moment,'mean')==1
                        MC_u = 1/(varmax-varmin)*cos(pi*XX)./TT.*(exp(-(D_u*pi^2+varmin).*TT)-exp(-(D_u*pi^2+varmax).*TT));
                    end
                    if strcmp(moment,'variance')==1
                        MC_u = cos(pi*XX).^2.*(1./(2*(varmax-varmin)*TT).*(exp(-2*(D_u*pi^2+varmin).*TT)-exp(-2*(D_u*pi^2+varmax).*TT)) -(1./((varmax-varmin)*TT).*(exp(-(D_u*pi^2+varmin).*TT)-exp(-(D_u*pi^2+varmax).*TT))).^2);
                    end
                    save('MC_u.mat','MC_u')
                end
                if dimension==2
                    xx = linspace(xL(1),xL(2),p+1);
                    xx(end)=[];
                    tt = linspace(0,T,snapshots+1);
                    tt(1)=[];
                    [X,Y] = meshgrid(xx,xx);
                    MC_u = zeros(p,p,snapshots);
                    if strcmp(moment,'mean')==1
                        for i=1:snapshots
                            MC_u(:,:,i) = 1/(varmax-varmin)*cos(pi*X).*cos(pi*Y)/tt(i).*(exp(-(2*D_u*pi^2+varmin)*tt(i))-exp(-(2*D_u*pi^2+varmax)*tt(i)));
                        end
                    end
                end
            end
            if strcmp(moment,"variance")==1 % Sum all coefficient functions starting from 1
                if dimension==1
                    PCE_u_temp = zeros(p,snapshots);
                    for j=1:N
                        PCE_u_temp = PCE_u_temp + PCE_u(j*p+1:(j+1)*p,:).^2;
                    end
                end
                if dimension==2
                    PCE_u_temp = zeros(snapshots,p,p);
                    for j=1:N
                        PCE_u_temp = PCE_u_temp + PCE_u(:,j*N+1:(j+1)*N,:).^2;
                    end
                end
                PCE_u = PCE_u_temp;
            end
            errorvec = geterrorvec(MC_u,PCE_u,xL,p,snapshots,dimension,relative);
            plot(xvec,errorvec)
            errorarray_Spectral(:,i3)=errorvec';
        end
        set(gca,'Yscale','log')
        title(strcat('u: ',spatialdiscr,' Random: ',varname,' Polynomial:',polynomialtype))
        legend(legendname)
        xlabel('T');
        ylabel('L2 Error u')
        if relative==1
            ylabel('rel. L2 Error u')
        end
        hold off
        plotname = sprintf(strcat('Errorplot_MC_vs_PCE_',moment,'_',spatialdiscr,'_',timestepping,'_','system=%d_p=%d_xL=%.2f..%.2f_T=%d_M=%d_Du=%.5f_Dv=%.5f_F=%.4f_K=%.4f_',varname,'=%.6f..%.6f.png'),system,p,xL(1),xL(2),T,Mvec(3),D_u,D_v,F,K,varmin,varmax');
        saveas(gcf,plotname)
        fprintf(strcat("Plot for ETDRK4 saved as: ",plotname,"\n"))
        %clf;
        %end
        %if spectral_on~=1
        %    spatialdiscr = 'FD'; timestepping = 'EE';
        %end
        %MC_filename = sprintf(strcat(type,'_',num2str(dimension),'D_',polynomialtype,'_runs=%d_system=%d_p=%d_xL=%.2f..%.2f_T=%d_M=%d_Du=%.5f_Dv=%.5f_F=%.4f_K=%.4f_',varname,'=%.6f..%.6f.mat'),runs,system,p,xL(1),xL(2),T,M_MC,D_u,D_v,F,K,varmin,varmax);
        %fprintf("Opening:\n")
        %disp(MC_filename)
        %MC = importdata(MC_filename);
        %MC_u = MC.u_total;
        
        %% Explicit Euler
        spatialdiscr = 'FD'; timestepping = 'EE';
        figure(12)
        hold on
        for i3=1:Nlen
            N = Nvec(i3);
            PCE_filename = sprintf(strcat('PCE_',spatialdiscr,'_',timestepping,'_',num2str(dimension),'D_',polynomialtype,'_system=%d_p=%d_xL=%.2f..%.2f_T=%d_M=%d_N=%d_Du=%.5f_Dv=%.5f_F=%.4f_K=%.4f_',varname,'=%.6f..%.6f.mat'),system,p,xL(1),xL(2),T,Mvec(1),N,D_u,D_v,F,K,varmin,varmax);
            fprintf("Opening:\n")
            disp(PCE_filename)
            PCE = importdata(PCE_filename);
            PCE_u = PCE.U_array;
            if system==6
                if dimension==1
                    xx = linspace(xL(1),xL(2),p+1);
                    xx(end)=[];
                    tt = linspace(0,T,snapshots+1);
                    tt(1)=[];
                    [TT,XX] = meshgrid(tt,xx);
                    if strcmp(moment,'mean')==1
                        MC_u = 1/(varmax-varmin)*cos(pi*XX)./TT.*(exp(-(D_u*pi^2+varmin).*TT)-exp(-(D_u*pi^2+varmax).*TT));
                    end
                    if strcmp(moment,'variance')==1
                        MC_u = cos(pi*XX).^2.*(1./(2*(varmax-varmin)*TT).*(exp(-2*(D_u*pi^2+varmin).*TT)-exp(-2*(D_u*pi^2+varmax).*TT)) -(1./((varmax-varmin)*TT).*(exp(-(D_u*pi^2+varmin).*TT)-exp(-(D_u*pi^2+varmax).*TT))).^2);
                    end
                    save('MC_u.mat','MC_u')
                end
                if dimension==2
                    xx = linspace(xL(1),xL(2),p+1);
                    xx(end)=[];
                    tt = linspace(0,T,snapshots+1);
                    tt(1)=[];
                    [X,Y] = meshgrid(xx,xx);
                    MC_u = zeros(p,p,snapshots);
                    if strcmp(moment,'mean')==1
                        for i=1:snapshots
                            MC_u(:,:,i) = 1/(varmax-varmin)*cos(pi*X).*cos(pi*Y)/tt(i).*(exp(-(2*D_u*pi^2+varmin)*tt(i))-exp(-(2*D_u*pi^2+varmax)*tt(i)));
                        end
                    end
                end
            end
            if strcmp(moment,"variance")==1 % Sum all coefficient functions starting from 1
                if dimension==1
                    PCE_u_temp = zeros(p,snapshots);
                    for j=1:N
                        PCE_u_temp = PCE_u_temp + PCE_u(j*p+1:(j+1)*p,:).^2;
                    end
                end
                if dimension==2
                    PCE_u_temp = zeros(snapshots,p,p);
                    for j=1:N
                        PCE_u_temp = PCE_u_temp + PCE_u(:,j*N+1:(j+1)*N,:).^2;
                    end
                end
                PCE_u = PCE_u_temp;
            end
            errorvec = geterrorvec(MC_u,PCE_u,xL,p,snapshots,dimension,relative);
            plot(xvec,errorvec)
            errorarray_FD_EE(:,i3)=errorvec';
        end
        set(gca,'Yscale','log')
        title(strcat('u: ',spatialdiscr,' (EE) Random: ',varname,' Polynomial:',polynomialtype))
        legend(legendname)
        xlabel('T');
        ylabel('L2 Error u')
        if relative==1
            ylabel('rel. L2 Error u')
        end
        hold off
        plotname = sprintf(strcat('Errorplot_MC_vs_PCE_',moment,'_',spatialdiscr,'_',timestepping,'_','system=%d_p=%d_xL=%.2f..%.2f_T=%d_M=%d_Du=%.5f_Dv=%.5f_F=%.4f_K=%.4f_',varname,'=%.6f..%.6f.png'),system,p,xL(1),xL(2),T,Mvec(1),D_u,D_v,F,K,varmin,varmax);
        saveas(gcf,plotname)
        fprintf(strcat("Plot for EE saved as: ",plotname,"\n"))
        %clf;

        %% ETDRDP
        spatialdiscr = 'FD'; timestepping = 'ETDRDPIF';
        figure(13)
        hold on
        for i3=1:Nlen
            N = Nvec(i3);
            PCE_filename = sprintf(strcat('PCE_',spatialdiscr,'_',timestepping,'_',num2str(dimension),'D_',polynomialtype,'_system=%d_p=%d_xL=%.2f..%.2f_T=%d_M=%d_N=%d_Du=%.5f_Dv=%.5f_F=%.4f_K=%.4f_',varname,'=%.6f..%.6f.mat'),system,p,xL(1),xL(2),T,Mvec(2),N,D_u,D_v,F,K,varmin,varmax);
            fprintf("Opening:\n")
            disp(PCE_filename)
            PCE = importdata(PCE_filename);                                                              
            PCE_u = PCE.U_array;
            
            if system==6
                if dimension==1
                    xx = linspace(xL(1),xL(2),p+1);
                    xx(end)=[];
                    tt = linspace(0,T,snapshots+1);
                    tt(1)=[];
                    [TT,XX] = meshgrid(tt,xx);
                    if strcmp(moment,'mean')==1
                        MC_u = 1/(varmax-varmin)*cos(pi*XX)./TT.*(exp(-(D_u*pi^2+varmin).*TT)-exp(-(D_u*pi^2+varmax).*TT));
                    end
                    if strcmp(moment,'variance')==1
                        MC_u = cos(pi*XX).^2.*(1./(2*(varmax-varmin)*TT).*(exp(-2*(D_u*pi^2+varmin).*TT)-exp(-2*(D_u*pi^2+varmax).*TT)) -(1./((varmax-varmin)*TT).*(exp(-(D_u*pi^2+varmin).*TT)-exp(-(D_u*pi^2+varmax).*TT))).^2);
                    end
                    save('MC_u.mat','MC_u')
                end
                if dimension==2
                    xx = linspace(xL(1),xL(2),p+1);
                    xx(end)=[];
                    tt = linspace(0,T,snapshots+1);
                    tt(1)=[];
                    [X,Y] = meshgrid(xx,xx);
                    MC_u = zeros(p,p,snapshots);
                    if strcmp(moment,'mean')==1
                        for i=1:snapshots
                            MC_u(:,:,i) = 1/(varmax-varmin)*cos(pi*X).*cos(pi*Y)/tt(i).*(exp(-(2*D_u*pi^2+varmin)*tt(i))-exp(-(2*D_u*pi^2+varmax)*tt(i)));
                        end
                    end
                end
            end 
            if strcmp(moment,"variance")==1 % Sum all coefficient functions starting from 1
                if dimension==1
                    PCE_u_temp = zeros(p,snapshots);
                    for j=1:N
                        PCE_u_temp = PCE_u_temp + PCE_u(j*p+1:(j+1)*p,:).^2;
                    end
                end
                if dimension==2
                    PCE_u_temp = zeros(snapshots,p,p);
                    for j=1:N
                        PCE_u_temp = PCE_u_temp + PCE_u(:,j*N+1:(j+1)*N,:).^2;
                    end
                end
                PCE_u = PCE_u_temp;
            end
            errorvec = geterrorvec(MC_u,PCE_u,xL,p,snapshots,dimension,relative);
            plot(xvec,errorvec)
            errorarray_FD_ETDRDP(:,i3)=errorvec';
        end
        set(gca,'Yscale','log')
        title(strcat('u: ',spatialdiscr,' (ETDRDP) Random: ',varname,' Polynomial:',polynomialtype))
        legend(legendname)
        xlabel('T');
        ylabel('L2 Error u')
        if relative==1
            ylabel('rel. L2 Error u')
        end
        hold off
        plotname = sprintf(strcat('Errorplot_MC_vs_PCE_',moment,'_',spatialdiscr,'_',timestepping,'_','system=%d_p=%d_xL=%.2f..%.2f_T=%d_M=%d_Du=%.5f_Dv=%.5f_F=%.4f_K=%.4f_',varname,'=%.6f..%.6f.png'),system,p,xL(1),xL(2),T,Mvec(2),D_u,D_v,F,K,varmin,varmax);
        saveas(gcf,plotname)
        fprintf(strcat("Plot for ETDRDP saved as: ",plotname,"\n"))
        %clf;
end
errorarray_FD_EE = [xvec',errorarray_FD_EE];
errorarray_FD_ETDRDP = [xvec',errorarray_FD_ETDRDP];
errorarray_Spectral = [xvec',errorarray_Spectral];
writematrix(errorarray_FD_EE,sprintf(strcat('errorarray_FD_EE_',moment,'_system=%d_D=%.5f.txt'),system,D_u));
writematrix(errorarray_FD_ETDRDP,sprintf(strcat('errorarray_FD_ETDRDP_',moment,'_system=%d_D=%.5f.txt'),system,D_u));
writematrix(errorarray_Spectral,sprintf(strcat('errorarray_Spectral_',moment,'_system=%d_D=%.5f.txt'),system,D_u));

save(sprintf(strcat('runtimes_',moment,'_system=%d_D=%.5f.mat'),system,D_u),'runtimes');
%MC_filename = sprintf(strcat('MC_',num2str(dimension),'D_',polynomialtype,'_system=%d_p=%d_xL=%.2f..%.2f_T=%d_M=%d_Du=%.5f_Dv=%.5f_F=%.4f_K=%.4f_',varname,'=%.6f..%.6f.mat'),system,p,xL(1),xL(2),T,M,D_u_scalar,D_v_scalar,F_scalar,K_scalar,varmin,varmax);
%PCE_filename = sprintf(strcat('PCE_',spatialdiscr,'_',timestepping,'_',num2str(dimension),'D_',polynomialtype,'_system=%d_p=%d_xL=%.2f..%.2f_T=%d_M=%d_N=%d_Du=%.5f_Dv=%.5f_F=%.4f_K=%.4f_',varname,'=%.6f..%.6f.mat'),system,p,xL(1),xL(2),T,M,N,D_u_scalar,D_v_scalar,F_scalar,K_scalar,varmin,varmax);
%plotname = sprintf(strcat('Errorplot_MC_vs_PCE_',spatialdiscr,'_',timestepping,'_','system=%d_p=%d_xL=%.2f..%.2f_T=%d_M=%d_Du=%.5f_Dv=%.5f_F=%.4f_K=%.4f_',varname,'=%.6f..%.6f.mat'),system,p,xL(1),xL(2),T,M,D_u_scalar,D_v_scalar,F_scalar,K_scalar,varmin,varmax');
end

