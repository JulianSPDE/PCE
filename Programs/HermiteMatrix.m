
% Input: Matrix size n,
% Function g input as string, e.g. '@(x) x.^2'

function A = HermiteMatrix(n,g)
    A = zeros(n);
    g = str2func(g);
    for i=0:n-1
        for j=0:n-1
            % The three functions in the integral:
            normi = sqrt(sqrt(2*pi)*factorial(i));
            normj = sqrt(sqrt(2*pi)*factorial(j));
            %norm2 = integral(@(x)exp(-x.^2).*hermiteh(j,x).*hermiteh(j,x),-Inf,Inf);
            intfun = @(x) exp(-x.^2/2).*hermiteh(i,x)/normi.*hermiteh(j,x)/normj.*g(x);
            A(i+1,j+1) = integral(intfun,-Inf,Inf);
        end
    end
    %A = A;
    %save('HermiteMatrix.mat','A')
end


