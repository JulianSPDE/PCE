%This program solves a system of two (possibly coupled, nonlinear) equations in either one or 
% two spatial dimensions. 
%% Inputs:
% spatialdiscr, timestepping: Spatial discretization and time stepping scheme. 
% Possible combinations: ('FD','EE'), ('FD','ETDRDPIF'), ('Spectral','ETDRK4'),
% p: Resolution in space
% xL: Length of the domain along one axis
% T: Final time
% M: Number of time steps
% D_u: First diffusion constant
% F,K: Parameters in the equation
% System: One of the systems of nonlinear equations. For the definitions of
% the different systems, check the solution routine below or the
% documentation in the PDF.
% Dimension: Can either be 1 or 2.
% plotfig: If it is 1, functions u and v at ten time points are plotted and saved
% (the number ten is hard-coded)
% snapshots: The total number of time instances where potentially something
% is plotted and the array at that time is stored
%% To be added: Input parameter 'time stepping scheme' or something like that to choose between explicit Euler,
%% implicit Euler or Crank-Nicolson
function [u_array,v_array] = Nonlinear_Solver(spatialdiscr,timestepping,p,xL,T,M,D_u,D_v,F,K,system,dimension,plotfig,anti_aliasing,snapshots)
                             
    u_array = [];
    v_array = [];
    %% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dt=T/M;  % time step
    xx=linspace(xL(1),xL(2),p+1);
    h=(xL(2)-xL(1))/(p); % spatial step size
    pp=p^2; 
    fignum = 0;
    snapshots = M/snapshots; % We use the loop with mod

%% Matrices setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(spatialdiscr,'FD')
        if strcmp(timestepping,'EE')
            if dimension == 1
                A=blktridiag(-2,1,1,p); % construct Laplacian matrix
                A(1,p)=1; A(p,1)=1; % periodic BCs
                A = A/(h^2);
            end
            if dimension == 2
                B_p = 1/(h^2)*blktridiag(2,-1,-1,p); 
                B_p(1,end) = -1/(h^2);
                B_p(end,1) = -1/(h^2);
                I_p = speye(p);
                A_1 = kron(I_p,B_p); A_2 = kron(B_p,I_p);
                A = -(A_1 + A_2);
            end
        end
        if strcmp(timestepping,'ETDRDPIF')==1
            if dimension == 1
                e=ones(p,1);
                B=spdiags([-e,2*e,-e],-1:1,p,p);
                B(1,p)=-1; B(p,1)=-1;
                % Now replacing with the matrix from Trefethen:
                %B1 = differential_matrix(p,xL)*2*pi/(xL(2)-xL(1));
                %B = B1*B1;
                %I=speye(p,p);
                D=speye(2,2);
                D(1,1)=D_u*dt/h^2; D(2,2)=D_v*dt/h^2;
                % Here also matrix replacement
                %D(1,1) = dt*D_u; D(2,2)=dt*D_v;
                % A1=kron(D,kron(I,B)); A2=kron(D,kron(B,I));
                A = kron(D,B); 
                I = speye(2*p,2*p);
                M1 = sparse(I+A); M2 = sparse(I+(1/3)*A); 
                M3 = sparse(I+(1/4)*A);
                %M11 = sparse(I+A); M22 = sparse(I+(1/3)*A); 
                %M33 = sparse(I+(1/4)*A);
                [L1,U1] =lu(M1); [L2,U2] =lu(M2); [L3,U3] =lu(M3);
                %[L11,U11] =lu(M11); [L22,U22] =lu(M22); [L33,U33] =lu(M33);            
            end
            if dimension == 2
                e=ones(p,1);
                B=spdiags([-e,2*e,-e],-1:1,p,p);
                B(1,p)=-1; B(p,1)=-1;
                % Now replacing with the matrix from Trefethen:
                %B1 = differential_matrix(p,xL)*2*pi/(xL(2)-xL(1));
                %B = B1*B1;
                I=speye(p,p);
                D=speye(2,2);
                D(1,1)=D_u*dt/h^2; D(2,2)=D_v*dt/h^2;
                % Here also matrix replacement
                %D(1,1) = dt*D_u; D(2,2)=dt*D_v;
                A1=kron(D,kron(I,B)); A2=kron(D,kron(B,I));
                pp2 = 2*p*p; I = speye(pp2); 
                M1 = sparse(I+A1); M2 = sparse(I+(1/3)*A1); 
                M3 = sparse(I+(1/4)*A1);
                M11 = sparse(I+A2); M22 = sparse(I+(1/3)*A2); 
                M33 = sparse(I+(1/4)*A2);
                [L1,U1] =lu(M1); [L2,U2] =lu(M2); [L3,U3] =lu(M3);
                [L11,U11] =lu(M11); [L22,U22] =lu(M22); [L33,U33] =lu(M33);
            end
        end
    end


    if strcmp(spatialdiscr,'Spectral')==1
        if dimension==1
            len = xL(2)-xL(1);
            wavenums = [0:p/2-1,0,-p/2+1:-1]'/(len/(2*pi));
            xi = wavenums;
            L_u = -D_u*xi.^2;
            L_v = -D_v*xi.^2;
            Fr=logical(zeros(p,1));
            Fr([p/2+1-round(p/6):p/2+round(p/6)])=1;
            alxi = Fr; ind = alxi;
            E_u = exp(dt*L_u); E_v = exp(dt*L_v);
            E2_u = exp(dt*L_u/2); E2_v = exp(dt*L_v/2);
            rou = 16;
            r = exp(1i*pi*((1:rou)-0.5)/rou);
            LR_u=dt*L_u(:,ones(rou,1))+r(ones(p,1),:);
            LR_v=dt*L_v(:,ones(rou,1))+r(ones(p,1),:);
            % ETDRK4 coefficients
            Q_u = dt*real(mean((exp(LR_u/2)-1)./LR_u,2));
            Q_v = dt*real(mean((exp(LR_v/2)-1)./LR_v,2));
            f1_u=dt*real(mean( (-4-LR_u+exp(LR_u).*(4-3*LR_u+LR_u.^2))./(LR_u.^3) ,2));
            f1_v=dt*real(mean( (-4-LR_v+exp(LR_v).*(4-3*LR_v+LR_v.^2))./(LR_v.^3) ,2));
            f2_u=dt*real(mean( (4+2*LR_u+exp(LR_u).*(-4+2*LR_u))./(LR_u.^3) ,2));     
            f2_v=dt*real(mean( (4+2*LR_v+exp(LR_v).*(-4+2*LR_v))./(LR_v.^3) ,2));       
            f3_u=dt*real(mean( (-4-3*LR_u-LR_u.^2+exp(LR_u).*(4-LR_u))./(LR_u.^3) ,2));      
            f3_v=dt*real(mean( (-4-3*LR_v-LR_v.^2+exp(LR_v).*(4-LR_v))./(LR_v.^3) ,2));               
        end
        if dimension==2
            len = xL(2)-xL(1);
            wavenums = [0:p/2-1,0,-p/2+1:-1]'/(len/(2*pi));
            [xi,eta] = ndgrid(wavenums,wavenums);
            L_u = -D_u*(eta.^2+xi.^2);
            L_v = -D_v*(eta.^2+xi.^2);
            L_u_vec = L_u(:); L_v_vec = L_v(:);
            Fr=logical(zeros(p,1));
            Fr([p/2+1-round(p/6):p/2+round(p/6)])=1;
            [alxi,aleta]=ndgrid(Fr,Fr);
            ind = alxi|aleta;
            E_u = exp(dt*L_u); E_v = exp(dt*L_v); % Here change back; matrix exponential needed?
            E2_u = exp(dt*L_u/2); E2_v = exp(dt*L_v/2);
            rou = 16; % Number of roots of unity
            r = exp(1i*pi*((1:rou)-0.5)/rou);
            LR_u=dt*L_u_vec(:,ones(rou,1))+r(ones(p^2,1),:);
            LR_v=dt*L_v_vec(:,ones(rou,1))+r(ones(p^2,1),:);
            % ETDRK4 coefficients
            Q_u = dt*real(mean((exp(LR_u/2)-1)./LR_u,2));
            Q_v = dt*real(mean((exp(LR_v/2)-1)./LR_v,2));
            f1_u=dt*real(mean( (-4-LR_u+exp(LR_u).*(4-3*LR_u+LR_u.^2))./(LR_u.^3) ,2));
            f1_v=dt*real(mean( (-4-LR_v+exp(LR_v).*(4-3*LR_v+LR_v.^2))./(LR_v.^3) ,2));
            f2_u=dt*real(mean( (4+2*LR_u+exp(LR_u).*(-4+2*LR_u))./(LR_u.^3) ,2));     
            f2_v=dt*real(mean( (4+2*LR_v+exp(LR_v).*(-4+2*LR_v))./(LR_v.^3) ,2));       
            f3_u=dt*real(mean( (-4-3*LR_u-LR_u.^2+exp(LR_u).*(4-LR_u))./(LR_u.^3) ,2));      
            f3_v=dt*real(mean( (-4-3*LR_v-LR_v.^2+exp(LR_v).*(4-LR_v))./(LR_v.^3) ,2));        
            f1_u = reshape(f1_u,p,p); f2_u = reshape(f2_u,p,p); f3_u = reshape(f3_u,p,p);
            f1_v = reshape(f1_v,p,p); f2_v = reshape(f2_v,p,p); f3_v = reshape(f3_v,p,p);
            Q_u = reshape(Q_u,p,p); Q_v = reshape(Q_v,p,p);
        end
    end
    
    %% Initial condition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(spatialdiscr,'FD')==1
        if dimension == 1
            w0=FUN_IC_1Ddistribution(p);
            %w0(1:p)=2*cos(pi*xx(1:end-1));
            %w0(p+1:2*p)=-2*cos(pi*xx(1:end-1));
            if system==6 || system==7 || system==8 || system==9
                w0(1:p) = cos(pi*xx(1:end-1));
                w0(p+1:end) = zeros(p,1);
                if system==9
                    w0(p+1:end) = -2*cos(pi*xx(1:end-1));
                end
                u0_fft = fft(w0(1:p)); v0_fft = fft(w0(p+1:2*p));
            end
        end
        if dimension == 2
            %[u0,v0] = fixed_rectangle(p);
            len = xL(2)-xL(1);
            [X,Y] = meshgrid(xx(1:end-1),xx(1:end-1));
            if system<6
                bump1 = exp(-150*((X-2/7).^2+(Y-2/7).^2));
                bump2 = exp(-150*((X+2/7).^2+(Y-2/7).^2));
                bump3 = exp(-150*((X-2/7).^2+(Y+2/7).^2));
                bump4 = exp(-150*((X+2/7).^2+(Y+2/7).^2));
                u0 = 1-0.25*(bump1+bump2+bump3+bump4);
                v0 = 0.25*(bump1+bump2+bump3+bump4);
                %u0 = 1-exp(-150*((X-len/2).^2+(Y-len/2).^2));
                %v0 = exp(-150*((X-len/2).^2+2*(Y-len/2).^2));  
                u0 = reshape(u0,p^2,1);
                v0 = reshape(v0,p^2,1);
                w0 = [u0;v0];
                w = zeros(2*p^2,1);
            end
            if system>=6
                u0 = cos(pi*X).*cos(pi*Y);
                u0 = reshape(u0,p^2,1);
                v0 = 0*u0;
                w0 = [u0;v0];
                w = zeros(2*p^2,1);
            end
        end
    end

    
    if strcmp(spatialdiscr,'Spectral')==1
        if dimension==1
            w0=FUN_IC_1Ddistribution(p);
            if system==6 || system==7 || system==8 || system==9
                w0(1:p) = cos(pi*xx(1:end-1));
                w0(p+1:end) = zeros(p,1);
                if system==9
                    w0(p+1:end) = -cos(pi*xx(1:end-1));
                end
            end
            %w0(1:p)=2*cos(pi*xx(1:end-1));
            %w0(p+1:2*p)=-2*cos(pi*xx(1:end-1));
            u0_fft = fft(w0(1:p)); v0_fft = fft(w0(p+1:2*p));
        end
        if dimension==2
            [X,Y] = meshgrid(xx(1:end-1),xx(1:end-1));
            if system<6
                bump1 = exp(-150*((X-2/7).^2+(Y-2/7).^2));
                bump2 = exp(-150*((X+2/7).^2+(Y-2/7).^2));
                bump3 = exp(-150*((X-2/7).^2+(Y+2/7).^2));
                bump4 = exp(-150*((X+2/7).^2+(Y+2/7).^2));
                u0 = 1-0.25*(bump1+bump2+bump3+bump4);
                v0 = 0.25*(bump1+bump2+bump3+bump4);
                %figure(1)
                %imagesc(linspace(-1,1,p),linspace(-1,1,p),u0)
                %figure(2)
                %imagesc(linspace(-1,1,p),linspace(-1,1,p),v0)
                %figure(3)
            end
            if system>=6
                u0 = cos(pi*X).*cos(pi*Y);
                u0 = reshape(u0,p^2,1);
                v0 = 0*u0;
                w0 = [u0;v0];
                w = zeros(2*p^2,1);
            end
            %u0 = 1-exp(-150*((X-len/2).^2+(Y-len/2).^2));
            %v0 = exp(-150*((X-len/2).^2+2*(Y-len/2).^2));
            %data=importdata('IC_xL2_n256.mat');
            %[u0,v0]
            %u0 = data(1:256^2);
            %v0 = data(256^2+1:end);
            u0 = reshape(u0,p,p);
            v0 = reshape(v0,p,p);
            u0_fft = fftn(u0);
            v0_fft = fftn(v0);
            %w0_fft = fftn(w0);
            %w = zeros(2*p^2,1);
        end
    end
    
    %% Time stepping %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Explicit Euler
    if strcmp(timestepping,'EE')==1
        if dimension==1
            w(:,1)=w0;
        end
        if dimension==2
            w = w0;
        end

        for k=1:M 
            t = k*dt;
            if dimension == 1
                w(:,k+1)=w(:,k)+Increment_Explicit_Euler_1D(dt,w(:,k),p,A,D_u,D_v,F,K,system);
            end
            if dimension == 2
                if system<6
                    w = w+Increment_Explicit_Euler_2D(dt,w,p,A,D_u,D_v,F,K,system);
                end
                if system>=6
                    w(1:pp) = w(1:pp)+Increment_Explicit_Euler_2D_onefunction(dt,w(1:pp),p,A,D_u,D_v,F,K,system);
                    %{
                    MM = Increment_Explicit_Euler_2D_onefunction(dt,w(1:pp),p,A,D_u,D_v,F,K,system);
                    save('MM_1.mat','MM')
                    M_U1 = D_u*A;
                    save('M_U1_1.mat','M_U1')
                    figure(19)
                    imagesc(linspace(-1,1,p),linspace(-1,1,p),reshape(MM,p,p))
                    figure(3)
                    %}
                end
            end
            if mod(k,snapshots)==0
                if dimension == 1
                    U = w(1:p,k+1); V = w(p+1:end,k+1);   
                    u_array = [u_array,U];  v_array = [v_array,V];
                    if plotfig == 1
                        plot_1D('Single','mean',spatialdiscr,timestepping,U,V,fignum,'_','_',t,p,xL,T,M,0,D_u,D_v,F,K,0,0,system);
                    end
                end
                if dimension == 2
                    U = w(1:pp); V = w(pp+1:end);
                    u_array = cat(3,u_array,reshape(U,p,p));
                    v_array = cat(3,v_array,reshape(V,p,p));
                    if plotfig == 1
                        plot_2D('Single',spatialdiscr,timestepping,U,V,fignum,'_','_',t,p,xL,T,M,0,D_u,D_v,F,K,0,0,system,anti_aliasing);
                    end
                end
                fignum = fignum + 1;
            end
        end
    end
    
    % ETDRDPIF
    if strcmp(timestepping,'ETDRDPIF')==1
        w(:,1)=w0;
        counter = 0;
        for k=1:M 
            counter = counter + 1;
            %disp(counter)
            t = k*dt;
            if dimension == 1
                u = w(1:p^dimension); v = w(p^dimension+1:end);
                % Fw = [-u.*v.^2+F*(1-u); u.*v.^2-(F+K)*v];
                Fw = [g_fun(system,u,v,F); h_fun(system,u,v,F,K)];
                w_star = U1\(L1\(w + dt*Fw));
                %w_star = U1\(L1\w_star);
                u_star = w_star(1:p^dimension); v_star = w_star(p^dimension+1:end);
                F_star = [g_fun(system,u_star,v_star,F); h_fun(system,u_star,v_star,F,K)];
                %a_1 = U2\(L2\w);
                %b_1 = U3\(L3\w);
                d_1 = U2\(L2\(9*w+2*dt*Fw + dt*F_star));
                d_2 = U3\(L3\(-8*w-(3/2)*dt*Fw-(dt/2)*F_star));               
                w = d_1 + d_2;
            end
            if dimension == 2
                u = w(1:p^dimension); v = w(p^dimension+1:end);
                % Fw = [-u.*v.^2+F*(1-u); u.*v.^2-(F+K)*v];
                Fw = [g_fun(system,u,v,F); h_fun(system,u,v,F,K)];
                w_star = U11\(L11\(w + dt*Fw));
                w_star = U1\(L1\w_star);
                u_star = w_star(1:p^dimension); v_star = w_star(p^dimension+1:end);
                F_star = [g_fun(system,u_star,v_star,F); h_fun(system,u_star,v_star,F,K)];
                a_1 = U2\(L2\w);
                b_1 = U3\(L3\w);
                c_1 = 9*a_1 - 8*b_1;
                a_2 = U2\(L2\Fw);
                b_2 = U3\(L3\Fw);
                c_2 = 9*a_2-8*b_2;
                d_1 = U22\(L22\(9*c_1+2*dt*c_2 + dt*F_star));
                d_2 = U33\(L33\(-8*c_1-(3/2)*dt*c_2-(dt/2)*F_star));               
                w = d_1 + d_2;
            end
            if mod(k,snapshots)==0
                if dimension == 1
                    U = w(1:p); V = w(p+1:end);   
                    u_array = [u_array,U];
                    v_array = [v_array,V];
                    if plotfig == 1
                        plot_1D('Single','mean',spatialdiscr,timestepping,U,V,fignum,'_','_',t,p,xL,T,M,0,D_u,D_v,F,K,0,0,system);
                    end
                end
                if dimension == 2
                    U = w(1:pp); V = w(pp+1:end);
                    %u_array = [u_array,U];
                    u_array = cat(3,u_array,reshape(U,p,p));
                    %v_array = [v_array,V];
                    v_array = cat(3,v_array,reshape(V,p,p));
                    if plotfig == 1
                        plot_2D('Single',spatialdiscr,timestepping,U,V,fignum,'_','_',t,p,xL,T,M,0,D_u,D_v,F,K,0,0,system,anti_aliasing);
                    end
                end
                fignum = fignum + 1;
            end
        end
    end

    % ETDRK4
    if strcmp(timestepping,'ETDRK4')==1
        if dimension==1
            u=u0_fft; v=v0_fft;
            for k=1:M
                t = k*dt;
                Nv_u = fft(g_fun(system,ifftn(u),ifftn(v),F));
                Nv_v = fft(h_fun(system,ifftn(u),ifftn(v),F,K));
                a_u = E2_u.*u + Q_u.*Nv_u; %a_uv = E2_u.*v + Q_u.*Nv_u; % First, second components
                a_v = E2_v.*v + Q_v.*Nv_v; %a_vv = E2_v.*v + Q_v.*Nv_v;
                Na_u = fft(g_fun(system,ifft(a_u),ifft(a_v),F));
                Na_v = fft(h_fun(system,ifft(a_u),ifft(a_v),F,K));
                b_u = E2_u.*u + Q_u.*Na_u;
                b_v = E2_v.*v + Q_v.*Na_v;
                Nb_u = fft(g_fun(system,ifft(b_u),ifft(b_v),F));
                Nb_v = fft(h_fun(system,ifft(b_u),ifft(b_v),F,K));
                c_u = E2_u.*a_u + Q_u.*(2*Nb_u-Nv_u);
                c_v = E2_v.*a_v + Q_v.*(2*Nb_v-Nv_v);
                Nc_u = fft(g_fun(system,ifft(c_u),ifft(c_v),F));
                Nc_v = fft(h_fun(system,ifft(c_u),ifft(c_v),F,K));
                u = E_u.*u + Nv_u.*f1_u + (Na_u+Nb_u).*f2_u + Nc_u.*f3_u;
                v = E_v.*v + Nv_v.*f1_v + (Na_v+Nb_v).*f2_v + Nc_v.*f3_v;
                if anti_aliasing==1
                    u(ind)=0;
                    v(ind)=0;
                end
                
                if mod(k,snapshots)==0
                    U = real(ifft(u)); V = real(ifft(v));
                    u_array = [u_array,U];
                    v_array = [v_array,V];
                    if plotfig == 1
                        plot_1D('Single','mean',spatialdiscr,timestepping,U,V,fignum,'_','_',t,p,xL,T,M,0,D_u,D_v,F,K,0,0,system)
                    end  
                end
            end
        end
        if dimension==2
            u=u0_fft; v=v0_fft; 
            for k=1:M 
                t = k*dt;
                % u, v in matrix form
                %u = reshape(w(1:p^2),p,p); v = reshape(w(p^2+1:2*p^2),p,p);
                Nv_u = fftn(g_fun(system,ifftn(u),ifftn(v),F));
                Nv_v = fftn(h_fun(system,ifftn(u),ifftn(v),F,K));
                a_u = E2_u.*u + Q_u.*Nv_u; %a_uv = E2_u.*v + Q_u.*Nv_u; % First, second components
                a_v = E2_v.*v + Q_v.*Nv_v; %a_vv = E2_v.*v + Q_v.*Nv_v;
                Na_u = fftn(g_fun(system,ifftn(a_u),ifftn(a_v),F));
                Na_v = fftn(h_fun(system,ifftn(a_u),ifftn(a_v),F,K));
                b_u = E2_u.*u + Q_u.*Na_u;
                b_v = E2_v.*v + Q_v.*Na_v;
                Nb_u = fftn(g_fun(system,ifftn(b_u),ifftn(b_v),F));
                Nb_v = fftn(h_fun(system,ifftn(b_u),ifftn(b_v),F,K));
                c_u = E2_u.*a_u + Q_u.*(2*Nb_u-Nv_u);
                c_v = E2_v.*a_v + Q_v.*(2*Nb_v-Nv_v);
                Nc_u = fftn(g_fun(system,ifftn(c_u),ifftn(c_v),F));
                Nc_v = fftn(h_fun(system,ifftn(c_u),ifftn(c_v),F,K));
                u = E_u.*u + Nv_u.*f1_u + (Na_u+Nb_u).*f2_u + Nc_u.*f3_u;
                v = E_v.*v + Nv_v.*f1_v + (Na_v+Nb_v).*f2_v + Nc_v.*f3_v;
                if anti_aliasing==1
                    u(ind)=0;
                    v(ind)=0;
                end
                
                if mod(k,snapshots)==0
                    U = reshape(real(ifftn(u)),p^2,1); V = reshape(real(ifftn(v)),p^2,1);
                    u_array = cat(3,u_array,reshape(U,p,p));
                    v_array = cat(3,v_array,reshape(U,p,p));
                    if plotfig == 1
                        plot_2D('Single',spatialdiscr,timestepping,U,V,fignum,'_','_',t,p,xL,T,M,0,D_u,D_v,F,K,0,0,system,anti_aliasing);
                    end
                    fignum = fignum + 1;
                end
            end
        end
    end

end

