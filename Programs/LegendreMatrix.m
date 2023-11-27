% this is a script which computes the scalar products  <x*Pm(x),Pn(x)>
% for Legendre polynomials Pn on a specified interval [a,b] up to a
% specified order N. The results are stored in a NxN matrix.

% By the substitution rule, the integral over the functions on [a,b] are
% the same as over [-1,1] multiplied with (b-a)/2

% Note: As opposed to the Hermite matrices, the feature of modifying the
% random variable by some function, like g(xi), is not implemented here
% yet. We only have the case of a uniformly distributed random variable on
% [a,b], where a and b can be changed
function A = LegendreMatrix(N,a,b)
    A = zeros(N);
    for i= 1:N
        for j=1:N
            intfun = conv((LegendrePoly_normalized(i-1)),(LegendrePoly_normalized(j-1))); % convolution, i.e. multiplication
            intfun = conv(intfun,[(b-a)/2,-((b-a)/2-b)]);
            q = polyint(intfun);
            A(i,j) = 0.5*diff(polyval(q,[-1 1]));
        end
    end
   % save('Amatrix.mat','A')
end
