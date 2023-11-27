% Input: Stacked array of size (p*(N+1) x p) which is to be reshaped into a
% vector

function U_stackedvector = reshapecompinverse(U_array,N,p)
    U_stackedvector = zeros((N+1)*p^2,1);
    for i=0:N
        U_stackedvector(i*p^2+1:(i+1)*p^2) = reshape(U_array(i*p+1:(i+1)*p,1:p),p,p);
    end
end