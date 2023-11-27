% this program uses the program Nonlinear_Solver.m to do Monte Carlo or
% other methods.
% Simulations with certain input parameters. 
% Inputs:
% type: Can be 'MC', 'GaussQuad', 'Sobol'
% moment: Can either be 'mean' or 'variance', according to this input, the sample mean or
% variance will be computed
% runs: Number of Monte Carlo runs
% p: Spatial resolution along one axis
% xL: Interval along one axis, e.g. [-1,1]
% T: Final time
% M: Number of time steps
% D_u, F, K: parameters in the nonlinear system
% varname: The name of the random variable as a string; 'D_u', 'F' or 'K'
% varmin, varmax: The range of the random variable (normal distribution will be
% outside with probability 1:1 million
% system: The number of the nonlinear system 
% polynomialtype: Either 'Legendre' or 'Hermite'. Careful: PCE coefficient
% functinos not correct for Hermite at the moment, usage of Hermite not
% advisable
% dimension: Spatial dimension of the domain, can be 1 or 2
% savefig: If 1, plots at 10 time points (10 hard-coded) are saved.
function [u_total,v_total] = MonteCarlo_Nonlinear_Solver(simtype,moment,spatialdiscr,timestepping,runs,p,xL,T,M,D_u,D_v,F,K,varname,varmin,varmax,system,polynomialtype,dimension,savefig,snapshots)
    variance_functions = 8; % Number of coefficient functions used to compute the variance
    %% Initializing arrays
    if dimension==1
        u_complete = zeros(p,snapshots,runs); 
        v_complete = zeros(p,snapshots,runs);
    end
    if dimension==2
        u_complete = zeros(p,p,snapshots,runs); 
        v_complete = zeros(p,p,snapshots,runs);
    end
    D_u_scalar = D_u;
    D_v_scalar = D_v;
    F_scalar = F;
    K_scalar = K;
    %% MC
    if strcmp(simtype,'MC')==1
        randomnums = zeros(1,runs);
        if strcmp(polynomialtype,'Hermite')==1
            fun = @(x) erf(x)-0.9999995;
            percentile = fzero(fun,1); % This is the value we need a*(varmax-varmin)/2 to have
            a = 2*percentile/(varmax-varmin); % Function: x/a + (varmax+varmin)/2
            bellcurvestring = strcat('@(x) x/',num2str(a),'+(',num2str(varmax),'+',num2str(varmin),')/2');
            bellcurvefun = str2func(bellcurvestring);
        end
        if strcmp(varname,'D_u')==1
            D_u = [varmin,varmax];
            S = D_u;
        end
        if strcmp(varname,'F')==1
            F = [varmin,varmax];
            S = F;
        end
        if strcmp(varname,'K')==1
            K = [varmin,varmax];
            S = K;
        end    
        for i=1:runs
            
            if strcmp(polynomialtype,'Legendre')==1
                R = rand*(S(2)-S(1))+S(1); % Uniform random number in the interval
            end
            if strcmp(polynomialtype,'Hermite')==1
                R = bellcurvefun(randn); % Normally distributed r.v. in the interval (with 0.999999 prob.)
            end
            randomnums(i) = R;
            if strcmp(varname,'D_u')==1
                [u_array,v_array] = Nonlinear_Solver(spatialdiscr,timestepping,p,xL,T,M,R,D_v,F,K,system,dimension,0,1,snapshots);
            end
            if strcmp(varname,'F')==1
                [u_array,v_array] = Nonlinear_Solver(spatialdiscr,timestepping,p,xL,T,M,D_u,D_v,R,K,system,dimension,0,1,snapshots);
            end
            if strcmp(varname,'K')==1
                [u_array,v_array] = Nonlinear_Solver(spatialdiscr,timestepping,p,xL,T,M,D_u,D_v,F,R,system,dimension,0,1,snapshots);
            end
            if dimension==1
                u_complete(:,:,i) = u_array; v_complete(:,:,i) = v_array;
            end
            if dimension==2
                u_complete(:,:,:,i) = u_array; v_complete(:,:,:,i) = v_array;
            end
            %u_total = u_total + u_array;
            %v_total = v_total + v_array;
            if mod(i,10)==0
                fprintf("Run %d of %d finished.\n",i,runs)
            end
        end
        if strcmp(moment,'mean')==1 % PCE simplifies to sample mean
            if dimension==1
                u_total = mean(u_complete,3);  v_total = mean(v_complete,3);
            end
            if dimension==2
                u_total = mean(u_complete,4);  v_total = mean(v_complete,4);
            end
        end
        if strcmp(moment,'variance')==1 % Need to compute coefficient functions u_i up to certain point,
            % Then compute sum of |u_i|^2.
            if dimension==1
                u_total = zeros(p,snapshots); v_total = zeros(p,snapshots);
                for ii=1:variance_functions % Loop starting from 1 (u_0 not considered)
                    u_ii = zeros(p,snapshots); v_ii = zeros(p,snapshots);
                    for j=1:runs
                        u_ii = u_ii + 0.5/runs*polyval(LegendrePoly_normalized(ii),(randomnums(j)-(S(1)+S(2))/2)/(S(2)-S(1))*2)*u_complete(:,:,j);
                        v_ii = v_ii + 0.5/runs*polyval(LegendrePoly_normalized(ii),(randomnums(j)-(S(1)+S(2))/2)/(S(2)-S(1))*2)*v_complete(:,:,j);
                        disp((randomnums(j)-(S(1)+S(2))/2)/(S(2)-S(1))*2)
                    end
                    u_total = u_total + u_ii.^2; v_total = v_total + v_ii.^2;
                end
            end
            if dimension==2
                u_total = zeros(p,p,snapshots);
                for ii=1:variance_functions % Loop starting from 1 (u_0 not considered)
                    u_ii = zeros(p,p,snapshots); v_ii = zeros(p,p,snapshots);
                    for j=1:runs
                        u_ii = u_ii + 0.5/runs*polyval(LegendrePoly_normalized(ii),(randomnums(j)-S(1))/(S(2)-S(1)))*u_complete(:,:,:,j);
                        v_ii = v_ii + 0.5/runs*polyval(LegendrePoly_normalized(ii),(randomnums(j)-S(1))/(S(2)-S(1)))*v_complete(:,:,:,j);
                    end
                    u_total = u_total + u_ii.^2;  v_total = v_total + v_ii.^2;
                end
            end
        end
    end

    %% Gaussian quadrature
    if strcmp(simtype,'GaussQuad')==1
        if strcmp(polynomialtype,'Legendre')==1 % Gauss-Legendre quadrature
            [quad_nodes,quad_weights] = lgwt(runs,varmin,varmax);
            quad_weights = quad_weights/(varmax-varmin); % A bug in the lgwt code?
        end
        if strcmp(polynomialtype,'Hermite')==1 % Gauss-Hermite quadrature
            [quad_nodes,quad_weights] = Gauss_Hermite_Quadrature(runs);
            meann = (varmin + varmax)/2;
            % Now compute simulation(node)*weight
            fun = @(x) erf(x)-0.9999995;
            percentile = fzero(fun,1); % This is the value we need a*(varmax-varmin)/2 to have
            a = 2*percentile/(varmax-varmin); % Function: x/a + (varmax+varmin)/2
            stddev = 1/a; % standard deviation for variable transformation
            quad_nodes = quad_nodes*sqrt(2)*stddev + meann; % sqrt(2) because probabilists Hermite (see end of Wikipedia article on Gauss-Hermite quadrature)
            quad_weights = quad_weights/sqrt(pi); % Weights adjusted due to variable substitution
        end

        % Computing integral
        if dimension==1
            u_total = zeros(p,snapshots); v_total = zeros(p,snapshots);
        end
        if dimension==2
            u_total = zeros(p,p,snapshots); v_total = zeros(p,p,snapshots);
        end
        if strcmp(moment,'mean')==1
            for i=1:runs
                R = quad_nodes(i);
                if strcmp(varname,'F')==1
                    [u_array,v_array] = Nonlinear_Solver(spatialdiscr,timestepping,p,xL,T,M,D_u,D_v,R,K,system,dimension,0,1,snapshots);
                end
                if strcmp(varname,'K')==1
                    [u_array,v_array] = Nonlinear_Solver(spatialdiscr,timestepping,p,xL,T,M,D_u,D_v,F,R,system,dimension,0,1,snapshots);
                end
                u_total = u_total + quad_weights(i)*u_array;
                v_total = v_total + quad_weights(i)*v_array;
            if mod(i,round(runs/5))==0
                fprintf("Quadrature point %d of %d done.\n",i,runs)
            end
            end
        end
        if strcmp(moment,'variance')==1 % Need to compute coefficient functions u_i up to certain point,
            % Then compute sum of |u_i^2|.
            for i=1:runs
                R = quad_nodes(i);
                if strcmp(varname,'F')==1
                    [u_array,v_array] = Nonlinear_Solver(spatialdiscr,timestepping,p,xL,T,M,D_u,D_v,R,K,system,dimension,0,1,snapshots);
                end
                if strcmp(varname,'K')==1
                    [u_array,v_array] = Nonlinear_Solver(spatialdiscr,timestepping,p,xL,T,M,D_u,D_v,F,R,system,dimension,0,1,snapshots);
                end
                if dimension==1
                    u_complete(:,:,i) = u_array;  v_complete(:,:,i) = v_array;
                end
                if dimension==2
                    u_complete(:,:,:,i) = u_array;  v_complete(:,:,:,i) = v_array;
                end
                if mod(i,round(runs/5))==0
                    fprintf("Quadrature point %d of %d done.\n",i,runs)
                end
            end
            if dimension==1
                u_total = zeros(p,snapshots);  v_total = zeros(p,snapshots);
                for ii=1:variance_functions % Loop starting from 1 (u_0 not considered)
                    u_ii = zeros(p,snapshots);  v_ii = zeros(p,snapshots);
                    for j=1:runs
                       u_ii = u_ii + quad_weights(j)*polyval(LegendrePoly_normalized(ii),(quad_nodes(j)-(varmin+varmax)/2)/(varmax-varmin)*2)*u_complete(:,:,j);
                       v_ii = v_ii + quad_weights(j)*polyval(LegendrePoly_normalized(ii),(quad_nodes(j)-(varmin+varmax)/2)/(varmax-varmin)*2)*v_complete(:,:,j);
                       %u_ii = u_ii + quad_weights(j)*polyval(LegendrePoly_normalized(ii),quad_nodes(j))*ones(p,snapshots,runs);%u_complete(:,:,j);
                       %v_ii = v_ii + quad_weights(j)*polyval(LegendrePoly_normalized(ii),quad_nodes(j))*ones(p,snapshots,runs);
                    end
                    u_total = u_total + u_ii.^2;  v_total = v_total + v_ii.^2;
                end
                figure(1)
                clf;
                hold on
                for i=1:5
                    plot(linspace(-1,1,p),u_complete(:,8,i))
                end
            end
            if dimension==2
                u_total = zeros(p,p,snapshots);  v_total = zeros(p,p,snapshots);
                for ii=1:variance_functions % Loop starting from 1 (u_0 not considered)
                    u_ii = zeros(p,p,snapshots);  v_ii = zeros(p,p,snapshots);
                    for j=1:runs
                        u_ii = u_ii + quad_weights(j)*polyval(LegendrePoly_normalized(ii),quad_nodes(j))*u_complete(:,:,:,j);
                        v_ii = v_ii + quad_weights(j)*polyval(LegendrePoly_normalized(ii),quad_nodes(j))*v_complete(:,:,:,j);
                    end
                    u_total = u_total + u_ii.^2;  v_total = v_total + v_ii.^2;
                end
            end
        end
   

    end
    

    %% Sobol sequence
    if strcmp(simtype,'Sobol')==1
        sobolpoints = i4_sobol_generate(1,runs,1);
        if strcmp(polynomialtype,'Hermite')==1
            fun = @(x) erf(x)-0.9999995;
            percentile = fzero(fun,1); % This is the value we need a*(varmax-varmin)/2 to have
            a = 2*percentile/(varmax-varmin); % Function: x/a + (varmax+varmin)/2
            bellcurvestring = strcat('@(x) x/',num2str(a),'+(',num2str(varmax),'+',num2str(varmin),')/2');
            bellcurvefun = str2func(bellcurvestring);
        end
        if strcmp(varname,'D_u')==1
            D_u = [varmin,varmax];
            S = D_u;
        end
        if strcmp(varname,'F')==1
            F = [varmin,varmax];
            S = F;
        end
        if strcmp(varname,'K')==1
            K = [varmin,varmax];
            S = K;
        end            
        for i=1:runs
            if strcmp(polynomialtype,'Legendre')==1
                R = sobolpoints(i)*(S(2)-S(1))+S(1); % Uniform random number in the interval
            end
            if strcmp(polynomialtype,'Hermite')==1
                R = bellcurvefun(sobolpoints(i)); % Normally distributed r.v. in the interval (with 0.999999 prob.)
            end
            if strcmp(varname,'D_u')==1
                [u_array,v_array] = Nonlinear_Solver(spatialdiscr,timestepping,p,xL,T,M,R,D_v,F,K,system,dimension,0,1,snapshots);
            end
            if strcmp(varname,'F')==1
                [u_array,v_array] = Nonlinear_Solver(spatialdiscr,timestepping,p,xL,T,M,D_u,D_v,R,K,system,dimension,0,1,snapshots);
            end
            if strcmp(varname,'K')==1
                [u_array,v_array] = Nonlinear_Solver(spatialdiscr,timestepping,p,xL,T,M,D_u,D_v,F,R,system,dimension,0,1,snapshots);
            end
            if dimension==1
                u_complete(:,:,i) = u_array; v_complete(:,:,i) = v_array;
            end
            if dimension==2
                u_complete(:,:,:,i) = u_array; v_complete(:,:,:,i) = v_array;
            end
            %u_total = u_total + u_array;
            %v_total = v_total + v_array;
            if mod(i,10)==0
                fprintf("Run %d of %d finished.\n",i,runs)
            end
        end
        if strcmp(moment,'mean')==1 % PCE simplifies to sample mean
            if dimension==1
                u_total = mean(u_complete,3);  v_total = mean(v_complete,3);
            end
            if dimension==2
                u_total = mean(u_complete,4);  v_total = mean(v_complete,4);
            end
        end
        if strcmp(moment,'variance')==1 % Need to compute coefficient functions u_i up to certain point,
            % Then compute sum of |u_i^2|.
            if dimension==1
                u_total = zeros(p,snapshots); v_total = zeros(p,snapshots);
                for ii=1:variance_functions % Loop starting from 1 (u_0 not considered)
                    u_ii = zeros(p,snapshots);  v_ii = zeros(p,snapshots);
                    for j=1:runs
                        u_ii = u_ii + 1/runs*polyval(LegendrePoly_normalized(ii),((sobolpoints(j)-0.5)/(varmax-varmin)*2))*u_complete(:,:,j);
                        v_ii = v_ii + 1/runs*polyval(LegendrePoly_normalized(ii),((sobolpoints(j)-0.5)/(varmax-varmin)*2))*v_complete(:,:,j);
                    end
                    u_total = u_total + u_ii.^2;  v_total = v_total + v_ii.^2;
                end
            end
            if dimension==2
                u_total = zeros(p,p,snapshots);
                for ii=1:variance_functions % Loop starting from 1 (u_0 not considered)
                    u_ii = zeros(p,p,snapshots);  v_ii = zeros(p,p,snapshots);
                    for j=1:runs
                        u_ii = u_ii + 0.5/runs*polyval(LegendrePoly_normalized(ii),((sobolpoints(j)-0.5)/(varmax-varmin)*2))*u_complete(:,:,:,j);
                        v_ii = v_ii + 0.5/runs*polyval(LegendrePoly_normalized(ii),((sobolpoints(j)-0.5)/(varmax-varmin)*2))*v_complete(:,:,:,j);
                    end
                    u_total = u_total + u_ii.^2;  v_total = v_total + v_ii.^2;
                end
            end
        end
    end

   %% Halton sequence
    if strcmp(simtype,'Halton')==1
        haltonpoints = HaltonSequence(runs,2); % Base two for now. Do other bases work better/worse?
        if strcmp(polynomialtype,'Hermite')==1
            fun = @(x) erf(x)-0.9999995;
            percentile = fzero(fun,1); % This is the value we need a*(varmax-varmin)/2 to have
            a = 2*percentile/(varmax-varmin); % Function: x/a + (varmax+varmin)/2
            bellcurvestring = strcat('@(x) x/',num2str(a),'+(',num2str(varmax),'+',num2str(varmin),')/2');
            bellcurvefun = str2func(bellcurvestring);
        end
        if strcmp(varname,'D_u')==1
            D_u = [varmin,varmax];
            S = D_u;
        end
        if strcmp(varname,'F')==1
            F = [varmin,varmax];
            S = F;
        end
        if strcmp(varname,'K')==1
            K = [varmin,varmax];
            S = K;
        end    
        for i=1:runs
            
            if strcmp(polynomialtype,'Legendre')==1
                R = haltonpoints(i)*(S(2)-S(1))+S(1); % Uniform random number in the interval
            end
            if strcmp(polynomialtype,'Hermite')==1
                R = bellcurvefun(haltonpoints(i)); % Normally distributed r.v. in the interval (with 0.999999 prob.)
            end
            if strcmp(varname,'D_u')==1
                [u_array,v_array] = Nonlinear_Solver(spatialdiscr,timestepping,p,xL,T,M,R,D_v,F,K,system,dimension,0,1,snapshots);
            end
            if strcmp(varname,'F')==1
                [u_array,v_array] = Nonlinear_Solver(spatialdiscr,timestepping,p,xL,T,M,D_u,D_v,R,K,system,dimension,0,1,snapshots);
            end
            if strcmp(varname,'K')==1
                [u_array,v_array] = Nonlinear_Solver(spatialdiscr,timestepping,p,xL,T,M,D_u,D_v,F,R,system,dimension,0,1,snapshots);
            end
            u_total = u_total + u_array;
            v_total = v_total + v_array;
            if mod(i,10)==0
                fprintf("Run %d of %d finished.\n",i,runs)
            end
        end
        u_total = u_total/runs;
        v_total = v_total/runs;
    end


    filename = sprintf(strcat(simtype,'_',moment,'_',num2str(dimension),'D_',polynomialtype,'_runs=%d_system=%d_p=%d_xL=%.2f..%.2f_T=%d_M=%d_Du=%.5f_Dv=%.5f_F=%.4f_K=%.4f_',varname,'=%.6f..%.6f.mat'),runs,system,p,xL(1),xL(2),T,M,D_u_scalar,D_v_scalar,F_scalar,K_scalar,varmin,varmax);
    %filename = sprintf(strcat(simtype,'_2D_V_%03d_',polynomialtype,'_t=%.4f_p=%d_xL=%.2f_T=%d_M=%d_Du=%.5f_F=%.4f_K=%.4f_',varname,'=%.6f..%.6f','system=%d.png'),fignum,t,p,xL,T,M,D_u,F,K,varmin,varmax,system);
    fprintf("Saved as:\n")
    disp(filename)
    save(filename,'u_total','v_total')
    fignum = 0;
    if savefig==1
        if dimension==1
            for i=1:snapshots
                t = T/snapshots*i;
                plot_1D(simtype,moment,spatialdiscr,timestepping,u_total(:,i),v_total(:,i),fignum,polynomialtype,varname,t,p,xL,T,M,0,D_u_scalar,D_v_scalar,F_scalar,K_scalar,varmin,varmax,system)
                fignum = fignum + 1;
            end
        end
        if dimension==2
            for i=1:snapshots
                t = T/snapshots*i;
                plot_2D(simtype,spatialdiscr,timestepping,u_total(:,:,i),v_total(:,:,i),fignum,polynomialtype,varname,t,p,xL,T,M,0,D_u_scalar,D_v_scalar,F_scalar,K_scalar,varmin,varmax,system,1)
                fignum = fignum + 1;
            end
        end
    end
    
end