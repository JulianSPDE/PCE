function U2 = u2(U,N,polynomialtype) % Input: vectors U, V of size p^d x (N+1) (column vectors!)
% (Basically U_0, U_1, ... U_N stacked on each other 
% This function takes vectorized functions u_i and v_i, i=0,...,N as they
% appear in a PCE expansion of the functions u and v, and outputs
% vectorized functions stored in the array F_U according to the PCE
% expansion of the function uv^2 (which will be relevant for the Gray-Scott
% system)
pp = length(U)/(N+1); % pp: total number of grid points
% N+1 = size(U,2);
%C = @(x,y,z) x+y+z;
if strcmp(polynomialtype,'Hermite')==1
    C = str2func('C_Hermite');
end
if strcmp(polynomialtype,'Legendre')==1
    C = str2func('C_Legendre');
end
%C = str2func('C_Legendre_normalized');
U2 = zeros(pp*(N+1),1); % same length as U, V
for eta=0:N
    for l=0:N
        for m=0:l
                for p=0:min(m,l-m)
                    if eta == l-2*p
                        U2(u_part(pp,eta)) = U2(u_part(pp,eta)) + U(u_part(pp,m)).*U(u_part(pp,l-m)) ...
                            *C(m,l-m,p);
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
