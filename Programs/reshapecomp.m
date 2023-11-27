% Input: Vector of size (N+1)*p^2 which is to be reshaped into a stacked
% array

function U_stackedarray = reshapecomp(U_vector,N,p)
    U_stackedarray = zeros((N+1)*p,p);
    for i=0:N
        U_stackedarray(i*p+1:(i+1)*p,1:p) = reshape(U_vector(i*p^2+1:(i+1)*p^2),p,p);
    end
end