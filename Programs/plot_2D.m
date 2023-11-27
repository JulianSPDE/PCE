
% Function which takes vectors U and V of discretized functions
% and plots them (or their expected value, in the PCE case, which is equal to u_0).
% Inputs:
% simtype: 'PCE' or 'MC', depending on whether MC or PCE simulations are
% being plotted. Discretized functions u and v as the input vectors U and
% V.

function plot_2D(simtype,spatialdiscr,timestepping,U,V,fignum,polynomialtype,varname,t,p,xL,T,M,N,D_u,D_v,F,K,varmin,varmax,system,anti_aliasing)
    U_0 = U(1:p^2);
    V_0 = V(1:p^2);
    xx = linspace(xL(1),xL(2),p+1);
    xx(end) = [];
    %% Plot U
    %figure(1);
    %f1 = figure('Visible', 'off', 'Renderer', 'painters');
    %set(f1, 'Visible', 'off');
    clf;
    imagesc(xx,xx,reshape(U_0,p,p))
    set(gca,'Ydir','normal')
    axis([xL(1),xL(2),xL(1),xL(2)])
    colorbar
    title(sprintf("t=%.5f",t))
    if strcmp(simtype,'PCE')==1 || strcmp(simtype,'GaussQuad')==1
        filename = sprintf(strcat(simtype,'_',spatialdiscr,'_',timestepping,'_2D_U_%04d_',polynomialtype,'_t=%.4f_p=%d_xL=%.2f..%.2f_T=%d_M=%d_N=%d_Du=%.5f_Dv=%.5f_F=%.4f_K=%.4f_',varname,'=%.6f..%.6f_system=%d.png'),fignum,t,p,xL(1),xL(2),T,M,N,D_u,D_v,F,K,varmin,varmax,system);
    end
    if strcmp(simtype,'MC')==1 % Same filename, just without N
        filename = sprintf(strcat(simtype,'_',spatialdiscr,'_',timestepping,'_2D_U_%04d_',polynomialtype,'_t=%.4f_p=%d_xL=%.2f..%.2f_T=%d_M=%d_Du=%.5f_Dv=%.5f_F=%.4f_K=%.4f_',varname,'=%.6f..%.6f_system=%d.png'),fignum,t,p,xL(1),xL(2),T,M,D_u,D_v,F,K,varmin,varmax,system);
    end
    if strcmp(simtype,'Single')==1 % Same filename, just without N
        filename = sprintf(strcat(simtype,'_',spatialdiscr,'_',timestepping,'_','AA=',num2str(anti_aliasing),'_2D_U_%04d_',polynomialtype,'_t=%.4f_p=%d_xL=%.2f..%.2f_T=%d_M=%d_Du=%.5f_Dv=%.5f_F=%.4f_K=%.4f_',varname,'=%.6f..%.6f_system=%d.png'),fignum,t,p,xL(1),xL(2),T,M,D_u,D_v,F,K,varmin,varmax,system);
    end
    saveas(gcf,filename)
    clear filename
    
    %% Plot V
    %figure(2)
    %f2 = figure('Visible', 'off', 'Renderer', 'painters');
    if system~=6 && system~=7 && system~=8
        clf;
        imagesc(xx,xx,reshape(V_0,p,p))
        set(gca,'Ydir','normal')
        axis([xL(1),xL(2),xL(1),xL(2)])
        colorbar
        title(sprintf("t=%.5f",t))
        %title(sprintf("t=%.5f",t))
        if strcmp(simtype,'PCE')==1 || strcmp(simtype,'GaussQuad')==1
            filename = sprintf(strcat(simtype,'_',spatialdiscr,'_',timestepping,'_2D_V_%04d_',polynomialtype,'_t=%.4f_p=%d_xL=%.2f..%.2f_T=%d_M=%d_N=%d_Du=%.5f_Dv=%.5f_F=%.4f_K=%.4f_',varname,'=%.6f..%.6f','system=%d.png'),fignum,t,p,xL(1),xL(2),T,M,N,D_u,D_v,F,K,varmin,varmax,system);
        end
        if strcmp(simtype,'MC')==1 % Same filename, just without N
            filename = sprintf(strcat(simtype,'_',spatialdiscr,'_',timestepping,'_2D_V_%04d_',polynomialtype,'_t=%.4f_p=%d_xL=%.2f..%.2f_T=%d_M=%d_Du=%.5f_Dv=%.5f_F=%.4f_K=%.4f_',varname,'=%.6f..%.6f','system=%d.png'),fignum,t,p,xL(1),xL(2),T,M,D_u,D_v,F,K,varmin,varmax,system);
        end
        if strcmp(simtype,'Single')==1 % Same filename, just without N
            filename = sprintf(strcat(simtype,'_',spatialdiscr,'_',timestepping,'_','AA=',num2str(anti_aliasing),'_2D_V_%04d_',polynomialtype,'_t=%.4f_p=%d_xL=%.2f..%.2f_T=%d_M=%d_Du=%.5f_Dv=%.5f_F=%.4f_K=%.4f_',varname,'=%.6f..%.6f','system=%d.png'),fignum,t,p,xL(1),xL(2),T,M,D_u,D_v,F,K,varmin,varmax,system);
        end
        saveas(gcf,filename)
    end
end