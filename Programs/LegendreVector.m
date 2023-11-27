% Vector containing the expected values of L_i(x)*x for x~U(a,b), L_i
% defined on (a,b), for N polynomials L_i (from 0 to N-1)
function v = LegendreVector(N,a,b)
    v = zeros(N,1);
    for i = 1:N
        intfun = LegendrePoly(i-1);
        intfun = conv(intfun,[(b-a)/2,-((b-a)/2-b)]);
        q = polyint(intfun);
        v(i) = (1/2)*diff(polyval(q,[-1 1]));
    end
    %save('vector.mat','A')
end

