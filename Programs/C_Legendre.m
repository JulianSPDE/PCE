% This function takes as input three integers and outputs the coefficients
% which appear in the addition theorem for normalized Hermite polynomials.
% p must be smaller than both a and b.
function C = C_Legendre(a,b,p)
    C = (A(p)*A(a-p)*A(b-p)/A(a+b-p))*(2*a+2*b-4*p+1)/(2*a+2*b-2*p+1)*sqrt((2*a+1)*(2*b+1)/((2*(a+b-2*p)+1)));
end


function A_r = A(r) % Input: nonnegative integer r
    A_r = rising_factorial(0.5,r)/factorial(r);
end

function c = rising_factorial(a,r) % Input: real number a, nonnegative integer r
    c = 1;
    for i=0:r-1
        c = c*(a+i);
    end
end



