% This program takes a PCE parameter N and a polynomial type and produces
% the multiplication tensor of third order, as suggested in the papers by
% Debusschere for example.
% Since this is for the PCE of third order, the tensor has four indices
% (three for the polynomial orders, one to store the different functions.)

function D = PCE_tensor(N,polynomialtype,varmin,varmax)
    
    D = zeros(N,N,N,N);
    if strcmp(polynomialtype,'Legendre')==1
        for eta=1:N
            for i=1:N
                for j=1:N
                    for k=1:N
                        intfun = conv((LegendrePoly_normalized(i-1)),(LegendrePoly_normalized(j-1))); % convolution, i.e. multiplication
                        intfun = conv(intfun,(LegendrePoly_normalized(k-1)));
                        intfun = conv(intfun,(LegendrePoly_normalized(eta-1)));
                        intfun = conv(intfun,[(varmax-varmin)/2,-((varmax-varmin)/2-varmax)]);
                        q = polyint(intfun);
                        intfun2 = conv((LegendrePoly_normalized(k-1)),(LegendrePoly_normalized(k-1)));
                        intfun2 = conv(intfun2,[(varmax-varmin)/2,-((varmax-varmin)/2-varmax)]);
                        q2 = polyint(intfun2);
                        D(eta,i,j,k) = 0.5*diff(polyval(q,[-1 1]));%/(0.5*diff(polyval(q2,[-1 1])));
                    end
                end
            end
        end
    end

    if strcmp(polynomialtype,'Hermite')==1
        fun = @(x) erf(x)-0.9999995;
        percentile = fzero(fun,1); % This is the value we need a*(varmax-varmin)/2 to have
        varmin = 2*percentile/(varmax-varmin); % Function: x/a + (varmax+varmin)/2
        bellcurvestring = strcat('@(x) x/',num2str(varmin),'+(',num2str(varmax),'+',num2str(varmin),')/2');
        g = str2func(bellcurvestring);
        %fun = @(x) erf(x)-0.9999995;
        %percentile = fzero(fun,1); % This is the value we need a*(varmax-varmin)/2 to have
        %a = 2*percentile/(varmax-varmin); % Function: x/a + (varmax+varmin)/2
        %mean = (varmin+varmax)/2;
        %variance = 1/(a^2);
        for eta=1:N
            for i=1:N
                for j=1:N
                    for k=1:N
                        norm_i = sqrt(factorial(i-1));
                        norm_j = sqrt(factorial(j-1));
                        norm_k = sqrt(factorial(k-1));
                        norm_eta = sqrt(factorial(eta-1));
                        %norm2 = integral(@(x)exp(-x.^2).*hermiteh(j,x).*hermiteh(j,x),-Inf,Inf);
                        intfun = @(x) exp(-x.^2/2)/sqrt(2*pi).*hermiteh(i-1,x)/norm_i.*hermiteh(j-1,x)/norm_j.*hermiteh(k-1,x)/norm_k.*hermiteh(eta-1,x)/norm_eta;
                        intfun2 = @(x) exp(-x.^2/2)/sqrt(2*pi).*hermiteh(eta-1,x)/norm_eta.*hermiteh(eta-1,x)/norm_eta;
                        D(eta,i,j,k) = integral(intfun,-Inf,Inf)/integral(intfun2,-Inf,Inf);
                    end
                end
            end
        end
    end

end