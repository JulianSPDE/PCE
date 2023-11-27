% We try in this function to implement the scheme to solve the Gray-Scott
% system with a random component using an explicit Euler scheme for the
% time being.
% Inputs: 
% spatialdiscr: 'FD' or 'Spectral'.
% timestepping: 'EE' or 'ETDRK4' (Explicit Euler/Exponential time stepping
% Runge Kutta 4th order). The combinations 'FD'+'EE' or 'Spectral'+'ETDRK4'
% are possible.
% p: Number of spatial steps 
% xL: Vector [x1,x2] with start and end point of domain
% T: Final time
% M: Number of time steps
% N: Number of PCE polynomials
% D_u, D_v, F, K: Parameters in the equations
% varname: The name of the variable which is random; 'D_u', 'F' or 'K'
% varmin, varmax: Minimum and maximum of the random variable (normal
% distribution is outside with probability 1e-6)
% polynomialtype: The type of orthogonal polynomials: either 'Legendre' or 'Hermite'
% dimension: The spatial dimension of the problem, can be 1 or 2
% savefig: 1 if plots should be saved, 0 if not

% (Name: PCE Gray Scott Explicit Euler)
% Nonlinear_Solver('Spectral','ETDRK4',1024,[0,1.25],5000,5000,0.00002,0.00001,0.04,0.06,0,2,1,0)
function [U_final,V_final,U_array,V_array] = PCE_solver(spatialdiscr,timestepping,moment,p,xL,T,M,N,D_u,D_v,F,K,varname,varmin,varmax,system,polynomialtype,dimension,anti_aliasing,savefig,snapshots)
                                             
h = (xL(2)-xL(1))/p;
dt = T/M;
xx = linspace(xL(1),xL(2),p+1);
elapsedtime = 0;
exploded = 0;
snapshots = M/snapshots;

% Variables for filenames
D_u_scalar = D_u;
F_scalar = F;
K_scalar = K;

U_array = []; V_array = [];

% Generate bell curve lying inside the specified interval.
% A N(0,1) random variable must lie within
% [-a*(varmax-varmin)/2,a*(varmax-varmin)/2] with 0.999999 probability.
if strcmp(polynomialtype,'Hermite')==1
    fun = @(x) erf(x)-0.9999995;
    percentile = fzero(fun,1); % This is the value we need a*(varmax-varmin)/2 to have
    a = 2*percentile/(varmax-varmin); % Function: x/a + (varmax+varmin)/2
    bellcurvestring = strcat('@(x) x/',num2str(a),'+(',num2str(varmax),'+',num2str(varmin),')/2');
    % Resulting normal distribution: Mean is (varmin+varmax)/2, variance is
    % 1/(a^2).
end

%% Generating PCE matrices
% The entered function is of a standard normal r.v., i.e. entering 'x.^2'
% means that D_u ~ xi^2, where xi ~ N(0,1). Entering a constant will make
% randomness disappear
% Matrix D_u
%D_u = LegendreMatrix(N+1,0.00002,0.001);
%D_u = HermiteMatrix(N+1,'@(x) 0.00015*x+(0.001+0.00001)/2');
D_u = D_u*speye(N+1);
if strcmp(varname,'D_u')==1
    if strcmp(polynomialtype,'Legendre')==1
        D_u = LegendreMatrix(N+1,varmin,varmax);
    end
    if strcmp(polynomialtype,'Hermite')==1
        D_u = HermiteMatrix(N+1,bellcurvestring);
    end
end

% Matrix D_v
%D_v = LegendreMatrix(N+1,'@(x) 0.00001');
D_v_scalar = D_v;
D_v = D_v*speye(N+1);


% Matrix F
F = F*speye(N+1);
%Fvec = zeros(N+1,1);
%Fvec(1) = F(1,1);
Fvec = F_scalar*ones(N+1,1);
if strcmp(varname,'F')==1
    if strcmp(polynomialtype,'Legendre')==1
        F = LegendreMatrix(N+1,varmin,varmax);
        Fvec = LegendreVector(N+1,varmin,varmax);
    end
    if strcmp(polynomialtype,'Hermite')==1
        F = HermiteMatrix(N+1,bellcurvestring);
        Fvec = HermiteVector(N+1,bellcurvestring);
    end
end

% Matrix K
K = K*speye(N+1);
if strcmp(varname,'K')==1
    if strcmp(polynomialtype,'Legendre')==1
        K = LegendreMatrix(N+1,varmin,varmax);
    end
    if strcmp(polynomialtype,'Hermite')==1
        K = HermiteMatrix(N+1,bellcurvestring);
    end
end


%% Finite difference matrices
if strcmp(spatialdiscr,'FD')==1
    if dimension == 1
        A_p = 1/(h^2)*blktridiag(2,-1,-1,p); A_p(1,end) = -1/(h^2); A_p(end,1) = -1/(h^2);
        I_p = speye(p);
        % Matrices needed for explicit Euler
        M_U1 = kron(D_u,-A_p); M_U2 = kron(F,I_p); % Careful: We consider the positive Laplacian -> main diagonal negative
        M_V1 = kron(D_v,-A_p); M_V2 = kron(F+K,I_p); % In deterministic case, diagonal block matrices
        % 1 function setting
        A = D_u_scalar*A_p*dt;
        I = speye(p,p);
        M1 = sparse(I+A); M2 = sparse(I+(1/3)*A); 
        M3 = sparse(I+(1/4)*A);
        [L1,U1] =lu(M1); [L2,U2] =lu(M2); [L3,U3] =lu(M3);   
        L11 = speye(p); L22 = speye(p); L33 = speye(p);
        U11 = speye(p); U22 = speye(p); U33 = speye(p);
        % 2 function setting
        %AA = kron(D,A_p); 
        %I = speye(2*p,2*p);
        %MM1 = sparse(I+AA); MM2 = sparse(I+(1/3)*AA); MM3 = sparse(I+(1/4)*AA);
        %II1 = inv(MM1); II2 = inv(MM2); II3 = inv(MM3);
        if system~=6 && system~=7 && system~=8
            D=speye(2,2);
            D(1,1)=D_u_scalar*dt; D(2,2)=D_v_scalar*dt;
            A = kron(D,A_p); I = speye(2*p);
            M1 = sparse(I+A); M2 = sparse(I+(1/3)*A); M3 = sparse(I+(1/4)*A);
            [L1,U1] =lu(M1); [L2,U2] =lu(M2); [L3,U3] =lu(M3);   
            L11 = speye(p); L22 = speye(p); L33 = speye(p);
            U11 = speye(p); U22 = speye(p); U33 = speye(p);
        end
    end
    if dimension == 2
        % 1 function setting
        
        D=D_u_scalar;
        A_p = 1/(h^2)*blktridiag(2,-1,-1,p); A_p(1,end) = -1/(h^2); A_p(end,1) = -1/(h^2);
         A_p = D*A_p;

        if system==6 || system==7 || system==8
            I_p = speye(p); I_p2 = speye(p^2);
            A_1 = kron(I_p,A_p); A_2 = kron(A_p,I_p);
            %AA_1 = kron(speye(N+1),A_1); AA_2 = kron(speye(N+1),A_2);
            M_U1 = kron(D_u,-(A_1+A_2)); M_U2 = kron(F,I_p2); % Careful: We consider the positive Laplacian -> main diagonal negative
            M_V1 = kron(D_v,-(A_1+A_2)); M_V2 = kron(F+K,I_p2); % In deterministic case, diagonal block matrices    
            M1 = sparse(I_p2+dt*A_1); M2 = sparse(I_p2+(dt/3)*A_1); 
            M3 = sparse(I_p2+(dt/4)*A_1);
            M11 = sparse(I_p2+dt*A_2); M22 = sparse(I_p2+(dt/3)*A_2); 
            M33 = sparse(I_p2+(dt/4)*A_2);
            [L1,U1] =lu(M1); [L2,U2] =lu(M2); [L3,U3] =lu(M3);
            [L11,U11] =lu(M11); [L22,U22] =lu(M22); [L33,U33] =lu(M33);
        end
        
        % 2 function setting
        if system~=6 && system~=7 && system~=8  
            I_p = speye(p); I_p2 = speye(p^2);
            D=speye(2,2);
            D(1,1)=D_u_scalar*dt; D(2,2)=D_v_scalar*dt;
            A_1 = kron(D,kron(I_p,A_p)); A_2 = kron(D,kron(A_p,I_p));
            A_1_U = kron(D_u_scalar,kron(I_p,A_p)); A_2_U = kron(D_u_scalar,kron(A_p,I_p));
            A_1_V = kron(D_v_scalar,kron(I_p,A_p)); A_2_V = kron(D_v_scalar,kron(A_p,I_p));
            M_U1 = kron(speye(N+1),-(A_1_U+A_2_U)); M_U2 = kron(F,I_p2); % Careful: We consider the positive Laplacian -> main diagonal negative
            M_V1 = kron(speye(N+1),-(A_1_V+A_2_V)); M_V2 = kron(F+K,I_p2); % In deterministic case, diagonal block matrices    
            I_2p2 = speye(2*p^2);
            M1 = sparse(I_2p2+A_1); M2 = sparse(I_2p2+(1/3)*A_1); 
            M3 = sparse(I_2p2+(1/4)*A_1);
            M11 = sparse(I_2p2+A_2); M22 = sparse(I_2p2+(1/3)*A_2); 
            M33 = sparse(I_2p2+(1/4)*A_2);
            [L1,U1] =lu(M1); [L2,U2] =lu(M2); [L3,U3] =lu(M3);
            [L11,U11] =lu(M11); [L22,U22] =lu(M22); [L33,U33] =lu(M33);
        end
    end
end

%% Spectral discretization setup
if strcmp(spatialdiscr,'Spectral')==1
    if dimension==1
        len = xL(2)-xL(1);
       % w0=FUN_IC_1Ddistribution(p);
       % u0=w0(1:p); v0=w0(p+1:2*p);
        u0 = cos(pi*xx(1:end-1)); v0 = -0*cos(pi*xx(1:end-1));
        if savefig==1
            plot_1D('Single',moment,spatialdiscr,timestepping,u0,v0,1,polynomialtype,varname,0,p,xL,T,M,N,D_u_scalar,D_v_scalar,F_scalar,K_scalar,varmin,varmax,system)
        end
        wavenums = [0:p/2-1,0,-p/2+1:-1]'/(len/(2*pi));
        xi = wavenums;
        L = -(xi.^2);
        % Anti-aliasing indices
        Fr = logical(zeros(p,1));
        Fr([p/2+1-round(p/4) : p/2+round(p/4)])=1;
        ind = Fr;
        ind_stacked = repmat(ind,N+1,1);
        % Arrays for later
        L_diag = spdiags(L(:),[0],p,p); % Now a p x p matrix
        L_u = kron(D_u,speye(p))*kron(speye(N+1),L_diag); % (N+1)p x (N+1)p matrix
        L_v = kron(D_v,speye(p))*kron(speye(N+1),L_diag);
        % Matrices from Kassam paper
        L_u_vec = diag(L_u);
        L_v_vec = diag(L_v);
        M_r = 16; % M_r: Number of complex roots
        r = exp(1i*pi*((1:M_r)-0.5)/M_r);
        r_stacked = repmat(r(ones(p,1),:),N+1,1);
        LR_u = dt*L_u_vec(:,ones(M_r,1))+r_stacked;
        LR_v = dt*L_v_vec(:,ones(M_r,1))+r_stacked;
        % Size of each: p x 1
        Q_u = dt*real(mean( (exp(LR_u/2)-1)./LR_u,2));
        %Q_u_stacked = reshapecomp(Q_u,N,p);
        f1_u = dt*real(mean((-4-LR_u+exp(LR_u).*(4-3*LR_u+LR_u.^2))./LR_u.^3 ,2));
        %f1_u_stacked = reshapecomp(f1_u,N,p);
        f2_u = dt*real(mean((4+2*LR_u+exp(LR_u).*(-4+2*LR_u))./(LR_u.^3) ,2));
        %f2_u_stacked = reshapecomp(f2_u,N,p);
        f3_u = dt*real(mean((-4-3*LR_u-LR_u.^2+exp(LR_u).*(4-LR_u))./(LR_u.^3) ,2));       
        %f3_u_stacked = reshapecomp(f3_u,N,p);
        Q_v = dt*real(mean( (exp(LR_v/2)-1)./LR_v,2));
        %Q_v_stacked = reshapecomp(Q_v,N,p);
        f1_v = dt*real(mean((-4-LR_v+exp(LR_v).*(4-3*LR_v+LR_v.^2))./(LR_v.^3) ,2));
        %f1_v_stacked = reshapecomp(f1_v,N,p);
        f2_v = dt*real(mean((4+2*LR_v+exp(LR_v).*(-4+2*LR_v))./LR_v.^3 ,2));       
        %f2_v_stacked = reshapecomp(f2_v,N,p);
        f3_v = dt*real(mean((-4-3*LR_v-LR_v.^2+exp(LR_v).*(4-LR_v))./(LR_v.^3) ,2));      
        %f3_v_stacked = reshapecomp(f3_v,N,p);
        E_u_diag = exp(dt*L_u_vec);             E_u_2_diag = exp(dt*L_u_vec/2);
        %E_u_stacked = reshapecomp(E_u_diag,N,p);  E_u_2_stacked = reshapecomp(E_u_2_diag,N,p);
        E_v_diag = exp(dt*L_v_vec);             E_v_2_diag = exp(dt*L_v_vec/2);
        %E_v_stacked = reshapecomp(E_v_diag,N,p);  E_v_2_stacked = reshapecomp(E_v_2_diag,N,p);
    end
    if dimension==2
        % Initial condition: Bump function
        len = xL(2)-xL(1);
        %[X,Y] = meshgrid(xx(1:end-1),xx(1:end-1));
        %u0 = 1-exp(-150*((X-len/2).^2+(Y-len/2).^2));
        %v0 = exp(-150*((X-len/2).^2+2*(Y-len/2).^2));
        %if savefig==1
        %    plot_2D('PCE',spatialdiscr,timestepping,u0,v0,1,polynomialtype,varname,T,p,xL,T,M,N,D_u_scalar,D_v_scalar,F_scalar,K_scalar,varmin,varmax,system,0)
        %end
        wavenums = [0:p/2-1,0,-p/2+1:-1]'/(len/(2*pi));
        [xi,eta] = ndgrid(wavenums,wavenums);
        L = -(xi.^2+eta.^2);
        % Anti-aliasing indices
        Fr = logical(zeros(p,1));
        Fr([p/2+1-round(p/4) : p/2+round(p/4)])=1;
        [alxi,aleta] = ndgrid(Fr,Fr);
        ind = alxi | aleta;
        ind_stacked = repmat(ind,N+1,1);
        % Arrays for later
        L_diag = spdiags(L(:),[0],p^2,p^2); % Now a p^2 x p^2 matrix
        L_u = kron(D_u,speye(p^2))*kron(speye(N+1),L_diag); % (N+1)p^2 x (N+1)p^2 matrix
        L_v = kron(D_v,speye(p^2))*kron(speye(N+1),L_diag);
        % Matrices from Kassam paper
        L_u_vec = diag(L_u);
        L_v_vec = diag(L_v);
        M_r = 16; % M_r: Number of complex roots
        r = exp(1i*pi*((1:M_r)-0.5)/M_r);
        r_stacked = repmat(r(ones(p^2,1),:),N+1,1);
        LR_u = dt*L_u_vec(:,ones(M_r,1))+r_stacked;
        LR_v = dt*L_v_vec(:,ones(M_r,1))+r_stacked;
        % Size of each: p^2 x 1
        Q_u = dt*real(mean( (exp(LR_u/2)-1)./LR_u,2));
        Q_u_stacked = reshapecomp(Q_u,N,p);
        f1_u = dt*real(mean((-4-LR_u+exp(LR_u).*(4-3*LR_u+LR_u.^2))./LR_u.^3 ,2));
        f1_u_stacked = reshapecomp(f1_u,N,p);
        f2_u = dt*real(mean((4+2*LR_u+exp(LR_u).*(-4+2*LR_u))./(LR_u.^3) ,2));
        f2_u_stacked = reshapecomp(f2_u,N,p);
        f3_u = dt*real(mean((-4-3*LR_u-LR_u.^2+exp(LR_u).*(4-LR_u))./(LR_u.^3) ,2));       
        f3_u_stacked = reshapecomp(f3_u,N,p);
        Q_v = dt*real(mean( (exp(LR_v/2)-1)./LR_v,2));
        Q_v_stacked = reshapecomp(Q_v,N,p);
        f1_v = dt*real(mean((-4-LR_v+exp(LR_v).*(4-3*LR_v+LR_v.^2))./(LR_v.^3) ,2));
        f1_v_stacked = reshapecomp(f1_v,N,p);
        f2_v = dt*real(mean((4+2*LR_v+exp(LR_v).*(-4+2*LR_v))./LR_v.^3 ,2));       
        f2_v_stacked = reshapecomp(f2_v,N,p);
        f3_v = dt*real(mean((-4-3*LR_v-LR_v.^2+exp(LR_v).*(4-LR_v))./(LR_v.^3) ,2));      
        f3_v_stacked = reshapecomp(f3_v,N,p);
        E_u_diag = exp(dt*L_u_vec);             E_u_2_diag = exp(dt*L_u_vec/2);
        E_u_stacked = reshapecomp(E_u_diag,N,p);  E_u_2_stacked = reshapecomp(E_u_2_diag,N,p);
        E_v_diag = exp(dt*L_v_vec);             E_v_2_diag = exp(dt*L_v_vec/2);
        E_v_stacked = reshapecomp(E_v_diag,N,p);  E_v_2_stacked = reshapecomp(E_v_2_diag,N,p);
    end
end


%% Initial condition
if dimension == 1
    U_start = zeros((N+1)*p,1); % Including all points except for the start point
    V_start = zeros((N+1)*p,1);
    rng(17);
    if system < 6
        startvec = FUN_IC_1Ddistribution(p);
    end
    if system >= 6
        startvec = zeros(1,2*p);
        fprintf("startvector set.\n")
        startvec(1:p) = cos(pi*xx(1:end-1));
        startvec(p+1:2*p) = -0*cos(pi*xx(1:end-1));
    end
    U_start(1:p) = startvec(1:p);
    V_start(1:p) = startvec(p+1:end);
    xx = linspace(xL(1),xL(2),p+1);
    %U_start(1:p) = 2*cos(pi*xx(1:end-1));
    %V_start(1:p) = -2*cos(pi*xx(1:end-1));
end
if dimension == 2
    U_start = zeros((N+1)*p^2,1);
    V_start = zeros((N+1)*p^2,1);
    %[U,V] = fixed_rectangle(p);
    %U_start(1:p^2) = U;
    %V_start(1:p^2) = V;
    len = xL(2)-xL(1);
    [X,Y] = meshgrid(xx(1:end-1),xx(1:end-1)); 
    if system<6
        bump1 = exp(-150*((X-2/7).^2+(Y-2/7).^2));
        bump2 = exp(-150*((X+2/7).^2+(Y-2/7).^2));
        bump3 = exp(-150*((X-2/7).^2+(Y+2/7).^2));
        bump4 = exp(-150*((X+2/7).^2+(Y+2/7).^2));
        u0 = 1-0.25*(bump1+bump2+bump3+bump4);
        v0 = 0.25*(bump1+bump2+bump3+bump4);
        U_start(1:p^2) = reshape(u0,p^2,1);
        V_start(1:p^2) = reshape(v0,p^2,1);
    end
    if system>=6
        u0 = cos(pi*X).*cos(pi*Y);
        figure(19)
        imagesc(xx(1:end-1),xx(1:end-1),u0)
        figure(2)
        U_start(1:p^2) = reshape(u0,p^2,1);
    end
    
end



%% Explicit Euler scheme
if strcmp(timestepping,'EE')==1
    U_temp = U_start;
    V_temp = V_start;
    %plot_expectation(U_start,N,0,0,system);
    fignum = 1;
    for i=1:M
        % Standard Gray-Scott 
        %[F_UV,~,~] = uv2(U_temp,V_temp,N);
        %U_new = U_temp + dt*(M_U1*U_temp + kron(Fvec,ones(p^2,1)) - M_U2*U_temp - F_UV);
        %V_new = V_temp + dt*(M_V1*V_temp                          - M_V2*V_temp + F_UV);
        
        % System 0 (Gray-Scott)
        [U_new,V_new] = PCE_Explicit_Euler_Term(U_temp,V_temp,dt,N,p,polynomialtype,varmin,varmax,M_U1,M_U2,M_V1,M_V2,Fvec,dimension,system);
        U_temp = U_new;
        V_temp = V_new;  
        if mod(i,snapshots)==0
            t = i*dt;
            if dimension == 1
                U_array = [U_array,U_new];
                V_array = [V_array,U_new];
                if savefig == 1
                    plot_1D('PCE',moment,spatialdiscr,timestepping,U_new,V_new,fignum,polynomialtype,varname,t,p,xL,T,M,N,D_u_scalar,D_v_scalar,F_scalar,K_scalar,varmin,varmax,system);
                end
            end
            if dimension == 2
                U_array = cat(3,U_array,reshape(U_temp(1:p^2),p,p));
                V_array = cat(3,V_array,reshape(V_temp(1:p^2),p,p));
                if savefig == 1
                    plot_2D('PCE',spatialdiscr,timestepping,U_new,V_new,fignum,polynomialtype,varname,t,p,xL,T,M,N,D_u_scalar,D_v_scalar,F_scalar,K_scalar,varmin,varmax,system);
                end
            end
            fignum = fignum + 1;
        end
        if mod(i,round(M/5))==0
            fprintf("Time step %d of %d done.\n",i,M)
        end
        if exploded == 0 && max(U_new)>1000 % If scheme has exploded, save the time in a file
            t = i*dt;
            filename = sprintf(strcat('PCE_Explosiontime_',num2str(dimension),'D_',polynomialtype,'_t=%.4f_p=%d_xL=%.2f..%.2f_T=%d_M=%d_N=%d_Du=%.5f_Dv=%.5f_F=%.4f_K=%.4f_',varname,'=%.6f..%.6f','system=%d.mat'),t,p,xL(1),xL(2),T,M,N,D_u_scalar,D_v_scalar,F_scalar,K_scalar,varmin,varmax,system);
            save(filename,'t')
            exploded = 1;
        end
    end
    U_final = U_new;
    V_final = V_new;
end

%% ETDRDP scheme
if strcmp(timestepping,'ETDRDPIF')==1
    U_temp = U_start;
    V_temp = V_start;
    %plot_expectation(U_start,N,0,0,system);
    fignum = 1;
    for i=1:M
        % Standard Gray-Scott 
        %[F_UV,~,~] = uv2(U_temp,V_temp,N);
        %U_new = U_temp + dt*(M_U1*U_temp + kron(Fvec,ones(p^2,1)) - M_U2*U_temp - F_UV);
        %V_new = V_temp + dt*(M_V1*V_temp                          - M_V2*V_temp + F_UV);
        
        % System 0 (Gray-Scott)
        %[U_new,V_new] = PCE_ETDRDPIF_Term(U_temp,V_temp,dt,N,p,polynomialtype,varmin,varmax,M_U1,M_U2,M_V1,M_V2,Fvec,dimension,system);
        if system==0 || system==1 || system==2 || system==3 || system==4 || system==5 || system==9
            [U_new,V_new] = PCE_ETDRDPIF_Term(U_temp,V_temp,dt,N,p,polynomialtype,L1,U1,L2,U2,L3,U3,L11,U11,L22,U22,L33,U33,F,K,Fvec,dimension,system);
        end
        if system==6 || system==7 || system==8
            U_new = PCE_ETDRDPIF_Term_onefunction(U_temp,dt,N,p,polynomialtype,L1,U1,L2,U2,L3,U3,L11,U11,L22,U22,L33,U33,F,K,Fvec,dimension,system);
            V_new = 0*U_new;
        end
        U_temp = U_new;
        V_temp = V_new;  
        if mod(i,snapshots)==0
            t = i*dt;
            if dimension == 1
                U_array = [U_array,U_new];
                V_array = [V_array,U_new];
                if savefig == 1
                    plot_1D('PCE',moment,spatialdiscr,timestepping,U_new,V_new,fignum,polynomialtype,varname,t,p,xL,T,M,N,D_u_scalar,D_v_scalar,F_scalar,K_scalar,varmin,varmax,system);
                end
            end
            if dimension == 2
                U_array = cat(3,U_array,reshape(U_temp(1:p^2),p,p));
                V_array = cat(3,V_array,reshape(V_temp(1:p^2),p,p));
                if savefig == 1
                    plot_2D('PCE',spatialdiscr,timestepping,U_new,V_new,fignum,polynomialtype,varname,t,p,xL,T,M,N,D_u_scalar,D_v_scalar,F_scalar,K_scalar,varmin,varmax,system);
                end
            end
            fignum = fignum + 1;
        end
        if mod(i,round(M/5))==0
            fprintf("Time step %d of %d done.\n",i,M)
        end
        if exploded == 0 && max(U_new)>1000 % If scheme has exploded, save the time in a file
            t = i*dt;
            filename = sprintf(strcat('PCE_Explosiontime_',num2str(dimension),'D_',polynomialtype,'_t=%.4f_p=%d_xL=%.2f..%.2f_T=%d_M=%d_N=%d_Du=%.5f_Dv=%.5f_F=%.4f_K=%.4f_',varname,'=%.6f..%.6f','system=%d.mat'),t,p,xL(1),xL(2),T,M,N,D_u_scalar,D_v_scalar,F_scalar,K_scalar,varmin,varmax,system);
            save(filename,'t')
            exploded = 1;
        end
    end
    U_final = U_new;
    V_final = V_new;
end



%% ETDRK4 scheme
% We have a subtlety to take care of here: MATLAB needs a 2D array in order
% to perform a 2D DFT using the command fft. As a result, while in the main
% document, we write everything as (stacked) vectors with only one column,
% we here have to componentwise reshape these vectors into 2D arrays
% anytime we want to perform a DFT. 
% CAREFUL: Using reshape on a stacked vector of functions will yield wrong
% results. We use the function reshapecomp for this.
if strcmp(timestepping,'ETDRK4')==1
    if dimension==1
        %U_start = zeros(1,(N+1)*p); V_start = zeros(1,(N+1)*p);
        %U_start(1:p) = u0; V_start(1:p) = v0;
        U_temp_hat_stacked = fftcomp(U_start,N,p); V_temp_hat_stacked = fftcomp(V_start,N,p);
        fignum = 1;
        counter = 0;
        for i=1:M
            counter = counter + 1;
            t = i*dt;
            % Assembling nonlinear term N_UV
            U_temp = ifftcomp(U_temp_hat_stacked,N,p); V_temp = ifftcomp(V_temp_hat_stacked,N,p);
            [N_UV_stacked_1,N_UV_stacked_2] = PCE_N_UV(U_temp,V_temp,Fvec,F,K,N,p,polynomialtype,system,1,varmin,varmax);
            %disp(size(N_UV_stacked_1))
            Nv_u = fftcomp(N_UV_stacked_1,N,p); Nv_v = fftcomp(N_UV_stacked_2,N,p);
            a_U = E_u_2_diag.*U_temp_hat_stacked + Q_u.*Nv_u;
            a_V = E_v_2_diag.*V_temp_hat_stacked + Q_v.*Nv_v;
            [N_a_1,N_a_2] = PCE_N_UV(ifftcomp(a_U,N,p),ifftcomp(a_V,N,p),Fvec,F,K,N,p,polynomialtype,system,1,varmin,varmax);
            Na_u = fftcomp(N_a_1,N,p); Na_v = fftcomp(N_a_2,N,p);
            b_U = E_u_2_diag.*U_temp_hat_stacked + Q_u.*Na_u;
            b_V = E_v_2_diag.*V_temp_hat_stacked + Q_v.*Na_v;
            [N_b_1,N_b_2] = PCE_N_UV(ifftcomp(b_U,N,p),ifftcomp(b_V,N,p),Fvec,F,K,N,p,polynomialtype,system,1,varmin,varmax);
            Nb_u = fftcomp(N_b_1,N,p); Nb_v = fftcomp(N_b_2,N,p);
            c_U = E_u_2_diag.*a_U + Q_u.*(2*Nb_u-Nv_u);
            c_V = E_v_2_diag.*a_V + Q_v.*(2*Nb_v-Nv_v);
            [N_c_1,N_c_2] = PCE_N_UV(ifftcomp(c_U,N,p),ifftcomp(c_V,N,p),Fvec,F,K,N,p,polynomialtype,system,1,varmin,varmax);
            Nc_u = fftcomp(N_c_1,N,p); Nc_v = fftcomp(N_c_2,N,p);
            U_new_hat_stacked = E_u_diag.*U_temp_hat_stacked + Nv_u.*f1_u + ...
            (Na_u+Nb_u).*f2_u + Nc_u.*f3_u;
            V_new_hat_stacked = E_v_diag.*V_temp_hat_stacked + Nv_v.*f1_v + ...
            (Na_v+Nb_v).*f2_v + Nc_v.*f3_v;
            if anti_aliasing==1
                U_new_hat_stacked(ind_stacked) = 0;
                V_new_hat_stacked(ind_stacked) = 0;
            end
            U_temp_hat_stacked = U_new_hat_stacked; V_temp_hat_stacked = V_new_hat_stacked;
            
            if mod(i,snapshots)==0
                if savefig==1
                    plot_1D('PCE',moment,spatialdiscr,timestepping,real(ifftcomp(U_temp_hat_stacked,N,p)),real(ifftcomp(V_temp_hat_stacked,N,p)),1,polynomialtype,varname,t,p,xL,T,M,N,D_u_scalar,D_v_scalar,F_scalar,K_scalar,varmin,varmax,system)
                end
                U_add = real(ifftcomp(U_temp_hat_stacked,N,p));
                V_add = real(ifftcomp(V_temp_hat_stacked,N,p));
                U_array = [U_array,U_add];
                V_array = [V_array,V_add];
            end
        if mod(i,round(M/5))==0
            fprintf("Time step %d of %d done.\n",i,M)
        end    
        end
        U_final = U_add; V_final = V_add;
    end
    if dimension==2
        U_start = reshapecomp(U_start,N,p); V_start = reshapecomp(V_start,N,p);
        U_temp_hat_stacked = fftncomp(U_start,N,p); V_temp_hat_stacked = fftncomp(V_start,N,p);
        fignum = 1;
        counter = 0;
        for i=1:M
            counter = counter + 1;
            t = i*dt;
            % Assembling nonlinear term N_UV
            %fprintf("1\n")
            U_temp = ifftncomp(U_temp_hat_stacked,N,p); V_temp = ifftncomp(V_temp_hat_stacked,N,p);
            %fprintf("2\n")
            [N_UV_stacked_1,N_UV_stacked_2] = PCE_N_UV(U_temp,V_temp,Fvec,F,K,N,p,polynomialtype,system,2,varmin,varmax);
            %fprintf("3\n")
            Nv_u = fftncomp(N_UV_stacked_1,N,p); Nv_v = fftncomp(N_UV_stacked_2,N,p);
            %fprintf("4\n")
            a_U = E_u_2_stacked.*U_temp_hat_stacked + Q_u_stacked.*Nv_u;
            %fprintf("5\n")
            a_V = E_v_2_stacked.*V_temp_hat_stacked + Q_v_stacked.*Nv_v;
            %fprintf("6\n")
            [N_a_1,N_a_2] = PCE_N_UV(ifftncomp(a_U,N,p),ifftncomp(a_V,N,p),Fvec,F,K,N,p,polynomialtype,system,2,varmin,varmax);
            %fprintf("7\n")
            Na_u = fftncomp(N_a_1,N,p); Na_v = fftncomp(N_a_2,N,p);
            %fprintf("8\n")
            b_U = E_u_2_stacked.*U_temp_hat_stacked + Q_u_stacked.*Na_u;
            %fprintf("9\n")
            b_V = E_v_2_stacked.*V_temp_hat_stacked + Q_v_stacked.*Na_v;
            %fprintf("10\n")
            [N_b_1,N_b_2] = PCE_N_UV(ifftncomp(b_U,N,p),ifftncomp(b_V,N,p),Fvec,F,K,N,p,polynomialtype,system,2,varmin,varmax);
            %fprintf("11\n")
            Nb_u = fftncomp(N_b_1,N,p); Nb_v = fftncomp(N_b_2,N,p);
            %fprintf("12\n")
            c_U = E_u_2_stacked.*a_U + Q_u_stacked.*(2*Nb_u-Nv_u);
            %fprintf("13\n")
            c_V = E_v_2_stacked.*a_V + Q_v_stacked.*(2*Nb_v-Nv_v);
            %fprintf("14\n")
            [N_c_1,N_c_2] = PCE_N_UV(ifftncomp(c_U,N,p),ifftncomp(c_V,N,p),Fvec,F,K,N,p,polynomialtype,system,2,varmin,varmax);
            %fprintf("15\n")
            Nc_u = fftncomp(N_c_1,N,p); Nc_v = fftncomp(N_c_2,N,p);
            U_new_hat_stacked = E_u_stacked.*U_temp_hat_stacked + Nv_u.*f1_u_stacked + ...
            (Na_u+Nb_u).*f2_u_stacked + Nc_u.*f3_u_stacked;
            V_new_hat_stacked = E_v_stacked.*V_temp_hat_stacked + Nv_v.*f1_v_stacked + ...
            (Na_v+Nb_v).*f2_v_stacked + Nc_v.*f3_v_stacked;
            if anti_aliasing==1
                U_new_hat_stacked(ind_stacked) = 0;
                V_new_hat_stacked(ind_stacked) = 0;
            end
            U_temp_hat_stacked = U_new_hat_stacked; V_temp_hat_stacked = V_new_hat_stacked;
            if mod(i,snapshots)==0
                if savefig==1
                    plot_2D('PCE',spatialdiscr,timestepping,reshapecompinverse(real(ifftncomp(U_temp_hat_stacked,N,p)),N,p),reshapecompinverse(real(ifftncomp(V_temp_hat_stacked,N,p)),N,p),1,polynomialtype,varname,t,p,xL,T,M,N,D_u_scalar,D_v_scalar,F_scalar,K_scalar,varmin,varmax,system,0)
                end
                U_add = real(ifftncomp(U_temp_hat_stacked,N,p));
                V_add = real(ifftncomp(V_temp_hat_stacked,N,p));
                U_array = cat(3,U_array,U_add(1:p,1:p));
                V_array = cat(3,V_array,V_add(1:p,1:p));
            end
            if mod(i,round(M/5))==0
                fprintf("Time step %d of %d done.\n",i,M)
            end
        end
        U_final = real(ifftncomp(U_temp_hat_stacked,N,p)); V_final = real(ifftncomp(V_temp_hat_stacked,N,p));
        U_final = reshapecomp(U_final,N,p); V_final = reshapecomp(V_final,N,p);
        U_final = U_final(1:p,1:p);
        %U_final = real(ifftncomp(U_final,N,p)); V_final = real(ifftncomp(V_final,N,p));
        %U_final = reshapecompinverse(U_final,N,p); V_final = reshapecompinverse(V_final,N,p);
    end
end
filename = sprintf(strcat('PCE_',spatialdiscr,'_',timestepping,'_',num2str(dimension),'D_',polynomialtype,'_system=%d_p=%d_xL=%.2f..%.2f_T=%d_M=%d_N=%d_Du=%.5f_Dv=%.5f_F=%.4f_K=%.4f_',varname,'=%.6f..%.6f.mat'),system,p,xL(1),xL(2),T,M,N,D_u_scalar,D_v_scalar,F_scalar,K_scalar,varmin,varmax);
save(filename,'U_array','V_array')
fprintf("Saved as:\n")
disp(filename)
end


    