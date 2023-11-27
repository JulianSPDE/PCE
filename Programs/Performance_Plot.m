% This function runs several iPCE as well as niPCE simulations for
% different numbers of time steps and plots the relative L2 error at given
% time T, depending on the number M of time steps.

% All inputs:
% moment can be 'mean' or 'variance'
% T is the end time (in the paper T=2)
% xL is the space interval (in the paper [-1,1])
% p is the spatial resolution (in the paper 128)
% Mvec is the vector containing the different time step numbers which are
% compared
% M_MC, runs: The number of time steps and quadrature points used to compute reference solution
% q: Number quadrature points for other MC simulations
% D_u, D_v diffusion constants
% system: The equation which is solved (in the paper, system=6,7,8 for
% equation with linear, quadratic and cubic term)
% dimension: Spatial dimension
% polynomialtype: In our case, 'Legendre' since the random input is
% uniformly distributed
% F,K: Equation parameters; varname, varmin, varmax: The random variable
% and its minimum and maximum
        
function Performance_Plot(moment,T,xL,p,Mvec,M_MC,runs,q,N,D_u,D_v,system,dimension,polynomialtype,F,K,varname,varmin,varmax,savefig,relative)
snapshots = 1;
Mlen = length(Mvec);
% For each M, need simulation array at T=2 and then compute error between
% them... in each plot, we want six lines: EE, ETDRDP and ETDRK4 for iPCE
% and niPCE.

% Step 1. Reference solution
fprintf("Computing reference solution.\n")
    xx = linspace(xL(1),xL(2),p+1);
    xx(end)=[];
    tt = linspace(0,T,snapshots+1);
    tt(1)=[];
    [TT,XX] = meshgrid(tt,xx);
    if system~=6 || dimension~=1
        MonteCarlo_Nonlinear_Solver('GaussQuad',moment,'Spectral','ETDRK4',runs,p,xL,T,M_MC,D_u,D_v,F,K,varname,varmin,varmax,system,polynomialtype,dimension,savefig,snapshots)
        MC_filename = sprintf(strcat('GaussQuad','_',moment,'_',num2str(dimension),'D_',polynomialtype,'_runs=%d_system=%d_p=%d_xL=%.2f..%.2f_T=%d_M=%d_Du=%.5f_Dv=%.5f_F=%.4f_K=%.4f_',varname,'=%.6f..%.6f.mat'),runs,system,p,xL(1),xL(2),T,M_MC,D_u,D_v,F,K,varmin,varmax);
        fprintf("Opening:\n")
        disp(MC_filename)
        MC = importdata(MC_filename);
        reference_solution = MC.u_total;
    end
    if dimension==1 && system==6
        if strcmp(moment,'mean')==1
            reference_solution = 1/(varmax-varmin)*cos(pi*XX)./TT.*(exp(-(D_u*pi^2+varmin).*TT)-exp(-(D_u*pi^2+varmax).*TT));
        end
        if strcmp(moment,'variance')==1
            reference_solution = cos(pi*XX).^2.*(1./(2*(varmax-varmin)*TT).*(exp(-2*(D_u*pi^2+varmin).*TT)-exp(-2*(D_u*pi^2+varmax).*TT)) -(1./((varmax-varmin)*TT).*(exp(-(D_u*pi^2+varmin).*TT)-exp(-(D_u*pi^2+varmax).*TT))).^2);
        end
    end
    if dimension==1
        reference_norm = sqrt(trapz(xx,reference_solution.^2,1));
    end
    if dimension==2
        reference_norm = sqrt(trapz(xx,trapz(xx,reference_solution.^2,1),2));
    end
% Step 2. Computing the comparison arrays and the error vectors
    if dimension==1
        array_EE_iPCE = zeros(p,Mlen);      array_EE_niPCE = zeros(p,Mlen);
        array_ETDRDP_iPCE = zeros(p,Mlen);  array_ETDRDP_niPCE = zeros(p,Mlen);
        array_ETDRK4_iPCE = zeros(p,Mlen);  array_ETDRK4_niPCE = zeros(p,Mlen);
    end
    if dimension==2
        array_EE_iPCE = zeros(p,p,Mlen);      array_EE_niPCE = zeros(p,p,Mlen);
        array_ETDRDP_iPCE = zeros(p,p,Mlen);  array_ETDRDP_niPCE = zeros(p,p,Mlen);
        array_ETDRK4_iPCE = zeros(p,p,Mlen);  array_ETDRK4_niPCE = zeros(p,p,Mlen);
    end
    fprintf("Computing PCE simulations.\n")
    for i=1:Mlen
        fprintf("M=%d in progress.\n",Mvec(i))
        % iPCE
        fprintf("Computing iPCE EE.\n")
        [U_final_EE,~,~,~] = PCE_solver('FD','EE','mean',p,xL,T,Mvec(i),N,D_u,D_v,F,K,varname,varmin,varmax,system,polynomialtype,dimension,1,savefig,snapshots);
        fprintf("Computing iPCE ETDRDP.\n")
        [U_final_ETDRDP,~,~,~] = PCE_solver('FD','ETDRDPIF','mean',p,xL,T,Mvec(i),N,D_u,D_v,F,K,varname,varmin,varmax,system,polynomialtype,dimension,1,savefig,snapshots);
        fprintf("Computing iPCE ETDRK4.\n")
        [U_final_ETDRK4,~,~,~] = PCE_solver('Spectral','ETDRK4','mean',p,xL,T,Mvec(i),N,D_u,D_v,F,K,varname,varmin,varmax,system,polynomialtype,dimension,1,savefig,snapshots);

        if dimension==1
            if strcmp(moment,'mean')
                array_EE_iPCE(:,i) = U_final_EE(1:p);
                array_ETDRDP_iPCE(:,i) = U_final_ETDRDP(1:p);
                array_ETDRK4_iPCE(:,i) = U_final_ETDRK4(1:p);
            end
            if strcmp(moment,'variance')
                for j=1:N
                    array_EE_iPCE(:,i) = array_EE_iPCE(:,i) + U_final_EE(j*p+1:(j+1)*p).^2;
                    array_ETDRDP_iPCE(:,i) = array_ETDRDP_iPCE(:,i) + U_final_ETDRDP(j*p+1:(j+1)*p).^2;
                    array_ETDRK4_iPCE(:,i) = array_ETDRK4_iPCE(:,i) + U_final_ETDRK4(j*p+1:(j+1)*p).^2;
                end
            end
        end
        if dimension==2
            if strcmp(moment,'mean')
                array_EE_iPCE(:,:,i) = U_final_EE(1:p,1:p);
                array_ETDRDP_iPCE(:,:,i) = U_final_ETDRDP(1:p,1:p);
                array_ETDRK4_iPCE(:,:,i) = U_final_ETDRK4(1:p,1:p);
            end
            if strcmp(moment,'variance')
                array_EE_iPCE(:,:,i) = array_EE_iPCE(:,:,i) + U_final_EE(j*p+1:(j+1)*p,1:p).^2;
                array_ETDRDP_iPCE(:,:,i) = array_ETDRDP_iPCE(:,:,i) + U_final_ETDRDP(j*p+1:(j+1)*p,1:p).^2;
                array_ETDRK4_iPCE(:,:,i) = array_ETDRK4_iPCE(:,:,i) + U_final_ETDRK4(j*p+1:(j+1)*p,1:p).^2;
            end
        end
        % niPCE
        fprintf("Computing niPCE EE.\n")
        [u_total_EE,~] = MonteCarlo_Nonlinear_Solver('GaussQuad','mean','FD','EE',q,p,xL,T,Mvec(i),D_u,D_v,F,K,varname,varmin,varmax,system,polynomialtype,dimension,savefig,snapshots);
        fprintf("Computing niPCE ETDRDPIF.\n")
        [u_total_ETDRDP,~] = MonteCarlo_Nonlinear_Solver('GaussQuad','mean','FD','ETDRDPIF',q,p,xL,T,Mvec(i),D_u,D_v,F,K,varname,varmin,varmax,system,polynomialtype,dimension,savefig,snapshots);
        fprintf("Computing niPCE ETDRK4.\n")
        [u_total_ETDRK4,~] = MonteCarlo_Nonlinear_Solver('GaussQuad','mean','Spectral','ETDRK4',q,p,xL,T,Mvec(i),D_u,D_v,F,K,varname,varmin,varmax,system,polynomialtype,dimension,savefig,snapshots);
        if dimension==1
            array_EE_niPCE(:,i) = u_total_EE;
            array_ETDRDP_niPCE(:,i) = u_total_ETDRDP;
            array_ETDRK4_niPCE(:,i) = u_total_ETDRK4;
        end
        if dimension==2
            array_EE_niPCE(:,:,i) = u_total_EE;
            array_ETDRDP_niPCE(:,:,i) = u_total_ETDRDP;
            array_ETDRK4_niPCE(:,:,i) = u_total_ETDRK4;
        end
    end
% Step 3. Creating the error arrays and plot and saving the error arrays for TikZ

% Need six error lines: Three for iPCE and three for niPCE (for EE, ETDRDP, ETDRK4)
errorvector_EE_iPCE = zeros(1,Mlen);        errorvector_EE_niPCE = zeros(1,Mlen);
errorvector_ETDRDP_iPCE = zeros(1,Mlen);    errorvector_ETDRDP_niPCE = zeros(1,Mlen);
errorvector_ETDRK4_iPCE = zeros(1,Mlen);    errorvector_ETDRK4_niPCE = zeros(1,Mlen);

for i=1:Mlen
    if i==1
        figure(2)
        plot(xx,array_ETDRK4_iPCE(:,1)-array_ETDRK4_niPCE(:,1))
        %figure(3)
        %plot(xx,array_ETDRK4_niPCE(:,1))
    end
    if dimension==1
        errorvector_EE_iPCE(i) = sqrt(trapz(xx,(array_EE_iPCE(:,i)-reference_solution).^2,1));
        errorvector_ETDRDP_iPCE(i) = sqrt(trapz(xx,(array_ETDRDP_iPCE(:,i)-reference_solution).^2,1));
        errorvector_ETDRK4_iPCE(i) = sqrt(trapz(xx,(array_ETDRK4_iPCE(:,i)-reference_solution).^2,1));
        errorvector_EE_niPCE(i) = sqrt(trapz(xx,(array_EE_niPCE(:,i)-reference_solution).^2,1));
        errorvector_ETDRDP_niPCE(i) = sqrt(trapz(xx,(array_ETDRDP_niPCE(:,i)-reference_solution).^2,1));
        errorvector_ETDRK4_niPCE(i) = sqrt(trapz(xx,(array_ETDRK4_niPCE(:,i)-reference_solution).^2,1));
    end
    if dimension==2
        errorvector_EE_iPCE(i) = sqrt(trapz(xx,trapz(xx,(array_EE_iPCE(:,:,i)-reference_solution).^2,2),1));
        errorvector_ETDRDP_iPCE(i) = sqrt(trapz(xx,trapz(xx,(array_ETDRDP_iPCE(:,:,i)-reference_solution).^2,2),1));
        errorvector_ETDRK4_iPCE(i) = sqrt(trapz(xx,trapz(xx,(array_ETDRK4_iPCE(:,:,i)-reference_solution).^2,2),1));
        errorvector_EE_niPCE(i) = sqrt(trapz(xx,trapz(xx,(array_EE_niPCE(:,:,i)-reference_solution).^2,2),1));
        errorvector_ETDRDP_niPCE(i) = sqrt(trapz(xx,trapz(xx,(array_ETDRDP_niPCE(:,:,i)-reference_solution).^2,2),1));
        errorvector_ETDRK4_niPCE(i) = sqrt(trapz(xx,trapz(xx,(array_ETDRK4_niPCE(:,:,i)-reference_solution).^2,2),1));
    end
end

if relative==1
    errorvector_EE_iPCE = errorvector_EE_iPCE/reference_norm;
    errorvector_ETDRDP_iPCE = errorvector_ETDRDP_iPCE/reference_norm;
    errorvector_ETDRK4_iPCE = errorvector_ETDRK4_iPCE/reference_norm;
    errorvector_EE_niPCE = errorvector_EE_niPCE/reference_norm;
    errorvector_ETDRDP_niPCE = errorvector_ETDRDP_niPCE/reference_norm;
    errorvector_ETDRK4_niPCE = errorvector_ETDRK4_niPCE/reference_norm;
end

%We don't want too big numbers in our plots
errorvector_EE_iPCE(errorvector_EE_iPCE>10)=NaN; errorvector_EE_niPCE(errorvector_EE_niPCE>10)=NaN;
errorvector_ETDRDP_iPCE(errorvector_ETDRDP_iPCE>10)=NaN; errorvector_ETDRDP_niPCE(errorvector_ETDRDP_niPCE>10)=NaN;
errorvector_ETDRK4_iPCE(errorvector_ETDRK4_iPCE>10)=NaN; errorvector_ETDRK4_niPCE(errorvector_ETDRK4_niPCE>10)=NaN;

figure(1)
hold on
plot(Mvec,errorvector_EE_iPCE)
plot(Mvec,errorvector_ETDRDP_iPCE)
plot(Mvec,errorvector_ETDRK4_iPCE)
plot(Mvec,errorvector_EE_niPCE)
plot(Mvec,errorvector_ETDRDP_niPCE)
plot(Mvec,errorvector_ETDRK4_niPCE)
set(gca,'Yscale','log')
set(gca,'Xscale','log')
disp(errorvector_EE_iPCE)
disp(errorvector_ETDRDP_iPCE)
disp(errorvector_ETDRK4_iPCE)
disp(errorvector_EE_niPCE)
disp(errorvector_ETDRDP_niPCE)
disp(errorvector_ETDRK4_niPCE)
legend("EE iPCE","ETDRDP iPCE","ETDRK4 iPCE","EE niPCE","ETDRDP niPCE","ETDRK4 niPCE")
%legend("EE iPCE","EE niPCE")
errorarray = [Mvec',errorvector_EE_iPCE',errorvector_ETDRDP_iPCE',errorvector_ETDRK4_iPCE',...
              errorvector_EE_niPCE',errorvector_ETDRDP_niPCE',errorvector_ETDRK4_niPCE'];
writematrix(errorarray,sprintf(strcat('Performanceplot_',moment,'_system=%d_D=%.5f.txt'),system,D_u));



end

