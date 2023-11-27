% Takes a vector of N+1 stacked vectors and returns the vector of stacked
% FFT'd vectors
function Utrans = fftcomp(U,N,p)
    Utrans = zeros((N+1)*p,1);
    for i=0:N
        Utrans(i*p+1:(i+1)*p) = fft(U(i*p+1:(i+1)*p));
    end
end
