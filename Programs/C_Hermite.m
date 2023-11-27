% This function takes as input three integers and outputs the coefficients
% which appear in the addition theorem for normalized Hermite polynomials.
% p must be smaller than both a and b.
function C = C_Hermite(a,b,p)
C = sqrt(nchoosek(a,p)*nchoosek(b,p)*nchoosek(a+b-2*p,a-p));
end

