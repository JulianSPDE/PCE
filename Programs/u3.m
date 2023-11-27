function [F_UV,tracker,trackernotzero] = u3(U,N,p,polynomialtype) % Input: vectors U, V of size p^d x (N+1) (column vectors!)
% (Basically U_0, U_1, ... U_N stacked on each other 
% This function takes vectorized functions u_i and v_i, i=0,...,N as they
% appear in a PCE expansion of the functions u and v, and outputs
% vectorized functions stored in the array F_U according to the PCE
% expansion of the function uv^2 (which will be relevant for the Gray-Scott
% system)
tracker = 0;
trackernotzero = 0;
pp=round(length(U)/(N+1));

% N+1 = size(U,2);
%C = @(x,y,z) x+y+z;
if strcmp(polynomialtype,'Hermite')==1
    C = str2func('C_Hermite');
end
if strcmp(polynomialtype,'Legendre')==1
    C = str2func('C_Legendre');
end
%C = str2func('C_Legendre_normalized');
F_UV = zeros(pp*(N+1),1); % same length as U, V
for eta=0:N
    for l=0:N
        for m=0:l
            for j=0:m
                for p=0:min(j,m-j)
                    for n=0:min(l-m,m-2*p)
                        F_UV(u_part(pp,eta)) = F_UV(u_part(pp,eta)) + U(u_part(pp,l-m)).*U(u_part(pp,j)).*U(u_part(pp,m-j)) ...
                            *C(j,m-j,p)*C(l-m,m-2*p,n)*krone(eta,l-2*p-2*n);
                        tracker = tracker + 1;
                        if krone(eta,l-2*p-2*n)~=0
                            trackernotzero = trackernotzero+1;
                        end
                    end
                end
            end
        end
    end
end

end

% Input: Big U or V vector with N+1 stacked vectorized functions, integer n
% between 0 and N
% Output: Indices of the nth part of the vector
function v = u_part(discr,n) 
    v = (n*discr+1:(n+1)*discr)';
end
