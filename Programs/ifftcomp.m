
% Takes a vector of N+1 stacked vectors  in Fourier space and returns the vector of stacked
% IFFT'd vectors

function Utrans = ifftcomp(U,N,p)
    Utrans = zeros((N+1)*p,1);
    for i=0:N
        Utrans(i*p+1:(i+1)*p) = ifft(U(i*p+1:(i+1)*p));
    end
end

