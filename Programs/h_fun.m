% These are the nonlinear functions h for the second component v in the
% coupled systems. The functions can be used as a function handle in
% other programs.

function result = h_fun(system,u,v,F,K)
    if system==0
        result = u.*v.^2-(F+K)*v;
    end   
    if system==1
        result = -u.*v -(F+K)*u;
    end
    if system==2
        result = -u-(F+K)*u;
    end
    if system==2.2
        result = -K*v;
    end
    if system==3
        result = -v-(F+K)*u;
    end
    if system==4
        result = -v+(F+K)*u.^2;
    end
    if system==5 || system==6 || system==7 || system==8
        result = 0*u;
    end
    if system==9
        result = u.*v-(F+K)*v;
    end
end