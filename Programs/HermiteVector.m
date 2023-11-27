% Input: (column) Vector length n,
% Function g input as string, e.g. '@(x) x.^2'
% This function uses the normalized Hermite polynomials
function v = HermiteVector(n,g)
    v = zeros(n,1);
    g = str2func(g);
    for i=0:n-1
        % The three functions in the integral:
        normi = sqrt(2*pi)*sqrt(factorial(i));
        %norm2 = integral(@(x)exp(-x.^2).*hermiteh(j,x).*hermiteh(j,x),-Inf,Inf);
        intfun = @(x) exp(-x.^2/2).*hermiteh_normalized(i,x)/normi.*g(x);
        v(i+1) = integral(intfun,-Inf,Inf);
    end
    %A = A;
    %save('HermiteMatrix.mat','A')
end



