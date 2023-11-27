% Function which takes vectors U and V of discretized functions
% and plots them (or their expected value, in the PCE case, which is equal to u_0).
% Inputs:
% simtype: 'PCE' or 'MC', depending on whether MC or PCE simulations are
% being plotted. Discretized functions u and v as the input vectors U and
% V.
function plot_1D(simtype,moment,spatialdiscr,timestepping,U,V,fignum,polynomialtype,varname,t,p,xL,T,M,N,D_u,D_v,F,K,varmin,varmax,system)  
    if strcmp(moment,'mean')==1 || strcmp(simtype,'PCE')~=1
        U_0 = real(U(1:p));
        V_0 = real(V(1:p));
    end
    if strcmp(moment,'variance')==1 && strcmp(simtype,'PCE')==1
        U_0 = zeros(1,p); V_0 = zeros(1,p);
        for i=1:N
            U_0 = U_0 + U(i*p+1:(i+1)*p).^2;
            V_0 = V_0 + V(i*p+1:(i+1)*p).^2;
        end
    end
    %U_1 = U(discr+1:2*discr);
    xx = linspace(xL(1),xL(2),p+1);
    xx(end) = [];
    figure(2)
    clf;
    hold on
    plot(xx,U_0);
    if system~=6 && system~=7 && system~=8
        plot(xx,V_0);
    end
    xlim([xL(1),xL(2)])
    hold off
    %legend("U_0","V_0")
    title(sprintf("t=%.5f",t))
    if strcmp(simtype,'PCE')==1
        filename = sprintf(strcat(simtype,'_',spatialdiscr,'_',timestepping,'_1D_%03d',polynomialtype,'_t=%.4f_p=%d_xL=%.2f..%.2f_T=%d_M=%d_N=%d_Du=%.5f_Dv=%.5f_F=%.4f_K=%.4f_',varname,'=%.6f..%.6f','system=%d.png'),fignum,t,p,xL(1),xL(2),T,M,N,D_u,D_v,F,K,varmin,varmax,system);
    end
    if strcmp(simtype,'MC')==1 % same filename, just without N
        filename = sprintf(strcat(simtype,'_',spatialdiscr,'_',timestepping,'_1D_%03d',polynomialtype,'_t=%.4f_p=%d_xL=%.2f..%.2f_T=%d_M=%d_Du=%.5f_Dv=%.5f_F=%.4f_K=%.4f_',varname,'=%.6f..%.6f','system=%d.png'),fignum,t,p,xL(1),xL(2),T,M,D_u,D_v,F,K,varmin,varmax,system);
    end
    if strcmp(simtype,'GaussQuad')==1 % same filename, just without N
        filename = sprintf(strcat(simtype,'_',spatialdiscr,'_',timestepping,'_1D_%03d',polynomialtype,'_t=%.4f_p=%d_xL=%.2f..%.2f_T=%d_M=%d_Du=%.5f_Dv=%.5f_F=%.4f_K=%.4f_',varname,'=%.6f..%.6f','system=%d.png'),fignum,t,p,xL(1),xL(2),T,M,D_u,D_v,F,K,varmin,varmax,system);
    end
    if strcmp(simtype,'Sobol')==1 % same filename, just without N
        filename = sprintf(strcat(simtype,'_',spatialdiscr,'_',timestepping,'_1D_%03d',polynomialtype,'_t=%.4f_p=%d_xL=%.2f..%.2f_T=%d_M=%d_Du=%.5f_Dv=%.5f_F=%.4f_K=%.4f_',varname,'=%.6f..%.6f','system=%d.png'),fignum,t,p,xL(1),xL(2),T,M,D_u,D_v,F,K,varmin,varmax,system);
    end
    if strcmp(simtype,'Halton')==1 % same filename, just without N
        filename = sprintf(strcat(simtype,'_',spatialdiscr,'_',timestepping,'_1D_%03d',polynomialtype,'_t=%.4f_p=%d_xL=%.2f..%.2f_T=%d_M=%d_Du=%.5f_Dv=%.5f_F=%.4f_K=%.4f_',varname,'=%.6f..%.6f','system=%d.png'),fignum,t,p,xL(1),xL(2),T,M,D_u,D_v,F,K,varmin,varmax,system);
    end
    if strcmp(simtype,'Single')==1 % same filename, just without N
        filename = sprintf(strcat(simtype,'_',spatialdiscr,'_',timestepping,'_1D_%03d','_t=%.4f_p=%d_xL=%.2f..%.2f_T=%d_M=%d_Du=%.5f_Dv=%.5f_F=%.4f_K=%.4f_system=%d.png'),fignum,t,p,xL(1),xL(2),T,M,D_u,D_v,F,K,system);
    end
    saveas(gcf,filename)
    clear filename
end