function a=hermiteh(n,x)
% Non-normalized Hermite polynomials

% Probabilists polynomials using explicit formula from the Wikipedia
a = 0;
for m=0:floor(n/2)
    a = a + (-1)^m/(factorial(m)*factorial(n-2*m)*2^m)*x.^(n-2*m);
end
a = a*factorial(n);
end