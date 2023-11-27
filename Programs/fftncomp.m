% This function accepts a stack of spatially discretized PCE functions
% (i.e. the array has dimension (N+1)*p x p) and 
% returns the componentwise FFT (i.e. each component is transformed and put
% back into the resulting vector)

% Input: The stacked array U, the number of nonconstant polynomials N
function Utrans = fftncomp(U,N,p)
    Utrans = zeros((N+1)*p,p);
    for i=0:N
        Utrans(i*p+1:(i+1)*p,1:p) = fftn(U(i*p+1:(i+1)*p,1:p));
    end
end


