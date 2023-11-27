% A little script to investigate convergence behaviour of MC and
% nonintrusive PCE

% MC 

function MC_PCE_convergence(spatialdiscr,timestepping,polynomialtype,runvector,p,xL,T,M,D_u,D_v,F,K,varname,varmin,varmax,system,dimension)
runlen = length(runvector); % 'runvector' contains the numbers of runs for different MC simulations, and the number of quadrature points for different Gauss quadrature simulations

Sobol_runvector = runvector;
Halton_runvector = runvector;
GaussQuad_runvector = runvector;
SC_runvector = runvector;

Gausslen = 500;
for i=1:length(runvector)
    if runvector(i)>Gausslen
        GaussQuad_runvector(i) = 0;
        SC_runvector(i) = 0;
    end
end
GaussQuad_runvector(GaussQuad_runvector==0)=[];
GaussQuadlen = length(GaussQuad_runvector);

mc_error = zeros(1,runlen-1);
gaussquad_error = zeros(1,GaussQuadlen-1);
sobol_error = zeros(1,runlen-1);
halton_error = zeros(1,runlen-1);
xx = linspace(xL(1),xL(2),p);
D_u_scalar = D_u; D_v_scalar = D_v; F_scalar = F; K_scalar = K;

%% Running simulations
%% MC
for i=1:runlen
    MonteCarlo_Nonlinear_Solver('MC',spatialdiscr,timestepping,runvector(i),p,xL,T,M,D_u,D_v,F,K,varname,varmin,varmax,system,polynomialtype,dimension,0,10);
    if i==runlen
        filename = sprintf(strcat('MC_',num2str(dimension),'D_',polynomialtype,'_runs=%d_system=%d_p=%d_xL=%.2f..%.2f_T=%d_M=%d_Du=%.5f_Dv=%.5f_F=%.4f_K=%.4f_',varname,'=%.6f..%.6f.mat'),runvector(i),system,p,xL(1),xL(2),T,M,D_u_scalar,D_v_scalar,F_scalar,K_scalar,varmin,varmax);
        MC_ref = importdata(filename);
        z_MC_ref = MC_ref.u_total;
    end
end

%% Gaussian quadrature
for i=1:GaussQuadlen
    MonteCarlo_Nonlinear_Solver('GaussQuad',spatialdiscr,timestepping,runvector(i),p,xL,T,M,D_u,D_v,F,K,varname,varmin,varmax,system,polynomialtype,dimension,0,10);
    if i==GaussQuadlen
        filename = sprintf(strcat('GaussQuad_',num2str(dimension),'D_',polynomialtype,'_runs=%d_system=%d_p=%d_xL=%.2f..%.2f_T=%d_M=%d_Du=%.5f_Dv=%.5f_F=%.4f_K=%.4f_',varname,'=%.6f..%.6f.mat'),runvector(i),system,p,xL(1),xL(2),T,M,D_u_scalar,D_v_scalar,F_scalar,K_scalar,varmin,varmax);
        GaussQuad_ref = importdata(filename);
        z_GaussQuad_ref = GaussQuad_ref.u_total;
    end
end



%% Sobol sequence
for i=1:runlen
    MonteCarlo_Nonlinear_Solver('Sobol',spatialdiscr,timestepping,runvector(i),p,xL,T,M,D_u,D_v,F,K,varname,varmin,varmax,system,polynomialtype,dimension,0,10);
    if i==runlen
        filename = sprintf(strcat('Sobol_',num2str(dimension),'D_',polynomialtype,'_runs=%d_system=%d_p=%d_xL=%.2f..%.2f_T=%d_M=%d_Du=%.5f_Dv=%.5f_F=%.4f_K=%.4f_',varname,'=%.6f..%.6f.mat'),runvector(i),system,p,xL(1),xL(2),T,M,D_u_scalar,D_v_scalar,F_scalar,K_scalar,varmin,varmax);
        Sobol_ref = importdata(filename);
        z_Sobol_ref = Sobol_ref.u_total;
    end
end

%% Halton sequence
for i=1:runlen
    MonteCarlo_Nonlinear_Solver('Halton',spatialdiscr,timestepping,runvector(i),p,xL,T,M,D_u,D_v,F,K,varname,varmin,varmax,system,polynomialtype,dimension,0,10);
    if i==runlen
        filename = sprintf(strcat('Halton_',num2str(dimension),'D_',polynomialtype,'_runs=%d_system=%d_p=%d_xL=%.2f..%.2f_T=%d_M=%d_Du=%.5f_Dv=%.5f_F=%.4f_K=%.4f_',varname,'=%.6f..%.6f.mat'),runvector(i),system,p,xL(1),xL(2),T,M,D_u_scalar,D_v_scalar,F_scalar,K_scalar,varmin,varmax);
        Halton_ref = importdata(filename);
        z_Halton_ref = Halton_ref.u_total;
    end
end


%% Importing and plotting data

% PCE simulation
for i=1:GaussQuadlen-1
    filename_GaussQuad = sprintf(strcat('GaussQuad_',num2str(dimension),'D_',polynomialtype,'_runs=%d_system=%d_p=%d_xL=%.2f..%.2f_T=%d_M=%d_Du=%.5f_Dv=%.5f_F=%.4f_K=%.4f_',varname,'=%.6f..%.6f.mat'),runvector(i),system,p,xL(1),xL(2),T,M,D_u_scalar,D_v_scalar,F_scalar,K_scalar,varmin,varmax);
    GaussQuad = importdata(filename_GaussQuad);
    z_GaussQuad = GaussQuad.u_total;
    if dimension==1
        z_GaussQuad = z_GaussQuad(:,end);
        z_GaussQuad_ref = z_GaussQuad_ref(:,end);
        gaussquad_error(i) = sqrt(trapz(xx,abs(z_GaussQuad-z_GaussQuad_ref).^2,1));
    end
    if dimension==2
        z_GaussQuad = z_GaussQuad(:,:,end);
        z_GaussQuad_ref = z_GaussQuad_ref(:,:,end);
        fprintf("GaussQuad error:\n")
        gaussquad_error(i) = sqrt(trapz(xx,trapz(xx,abs(z_GaussQuad-z_GaussQuad_ref).^2,2),1));
    end
end



% MC simulation
for i=1:runlen-1
    filename_MC = sprintf(strcat('MC_',num2str(dimension),'D_',polynomialtype,'_runs=%d_system=%d_p=%d_xL=%.2f..%.2f_T=%d_M=%d_Du=%.5f_Dv=%.5f_F=%.4f_K=%.4f_',varname,'=%.6f..%.6f.mat'),runvector(i),system,p,xL(1),xL(2),T,M,D_u_scalar,D_v_scalar,F_scalar,K_scalar,varmin,varmax);
    MC = importdata(filename_MC);
    z_MC = MC.u_total;
    if dimension==1
        z_MC_ref = z_MC_ref(:,end);
        z_MC = z_MC(:,end);
        disp(size(z_MC))
        disp(size(z_MC_ref))
        mc_error(i) = sqrt(trapz(xx,abs(z_MC-z_MC_ref).^2,1));
    end
    if dimension==2
        z_MC_ref = z_MC_ref(:,:,end);
        z_MC = z_MC(:,:,end);
        fprintf("MC error:\n")
        mc_error(i) = sqrt(trapz(xx,trapz(xx,abs(z_MC-z_MC_ref).^2,2),1)); % Switched MC_ref to PCE_ref, as it is clearly more exact
    end
end

% Sobol simulation
for i=1:runlen-1
    filename_Sobol = sprintf(strcat('Sobol_',num2str(dimension),'D_',polynomialtype,'_runs=%d_system=%d_p=%d_xL=%.2f..%.2f_T=%d_M=%d_Du=%.5f_Dv=%.5f_F=%.4f_K=%.4f_',varname,'=%.6f..%.6f.mat'),runvector(i),system,p,xL(1),xL(2),T,M,D_u_scalar,D_v_scalar,F_scalar,K_scalar,varmin,varmax);
    Sobol = importdata(filename_Sobol);
    z_Sobol = Sobol.u_total;
    if dimension==1
        z_Sobol_ref = z_Sobol_ref(:,end);
        z_Sobol = z_Sobol(:,end);
        sobol_error(i) = sqrt(trapz(xx,abs(z_Sobol-z_Sobol_ref).^2,1));
    end
    if dimension==2
        z_Sobol_ref = z_Sobol_ref(:,:,end);
        z_Sobol = z_Sobol(:,:,end);
        fprintf("Sobol error:\n")
        sobol_error(i) = sqrt(trapz(xx,trapz(xx,abs(z_Sobol-z_Sobol_ref).^2,2),1)); % Switched MC_ref to PCE_ref, as it is clearly more exact
    end
end

% Halton simulation
for i=1:runlen-1
    filename_Halton = sprintf(strcat('Halton_',num2str(dimension),'D_',polynomialtype,'_runs=%d_system=%d_p=%d_xL=%.2f..%.2f_T=%d_M=%d_Du=%.5f_Dv=%.5f_F=%.4f_K=%.4f_',varname,'=%.6f..%.6f.mat'),runvector(i),system,p,xL(1),xL(2),T,M,D_u_scalar,D_v_scalar,F_scalar,K_scalar,varmin,varmax);
    Halton = importdata(filename_Halton);
    z_Halton = Halton.u_total;
    if dimension==1
        z_Halton_ref = z_Halton_ref(:,end);
        z_Halton = z_Halton(:,end);
        halton_error(i) = sqrt(trapz(xx,abs(z_Halton-z_Halton_ref).^2,1));
    end
    if dimension==2
        z_Halton_ref = z_Halton_ref(:,:,end);
        z_Halton = z_Halton(:,:,end);
        fprintf("Halton error:\n")
        halton_error(i) = sqrt(trapz(xx,trapz(xx,abs(z_Halton-z_Halton_ref).^2,2),1)); % Switched MC_ref to PCE_ref, as it is clearly more exact
    end
end


%{
    if dimension==1
        mc_error(i) = sqrt(trapz(xx,abs(z_MC_ref-z_PCE_ref).^2,1));
    end
    if dimension==2
        mc_error(i) = sqrt(trapz(xx,trapz(xx,abs(z_MC_ref-z_PCE_ref).^2,2),1)); % Switched MC_ref to PCE_ref, as it is clearly more exact
    end
%}
runvector(end) = [];
figure(1)
plot(runvector,mc_error)
xlabel('MC runs')
ylabel('error')
set(gca,'Yscale','log')
saveas(gca,strcat(sprintf('Error_MC_%dD_T=%d_',dimension,T),polynomialtype,'.png'))

GaussQuad_runvector(end)=[];
figure(2)
plot(GaussQuad_runvector,gaussquad_error)
xlabel('Quadrature points')
ylabel('error')
set(gca,'Yscale','log')
saveas(gca,strcat(sprintf('Error_GaussQuad_%dD_T=%d_',dimension,T),polynomialtype,'.png'))

Sobol_runvector(end)=[];
figure(3)
plot(Sobol_runvector,sobol_error)
xlabel('Sobol points')
ylabel('error')
set(gca,'Yscale','log')
saveas(gca,strcat(sprintf('Error_Sobol_%dD_T=%d_',dimension,T),polynomialtype,'.png'))

Halton_runvector(end)=[];
figure(4)
plot(Halton_runvector,halton_error)
xlabel('Halton points')
ylabel('error')
set(gca,'Yscale','log')
saveas(gca,strcat(sprintf('Error_Halton_%dD_T=%d_',dimension,T),polynomialtype,'.png'))

end




