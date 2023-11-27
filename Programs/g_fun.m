% These are the nonlinear functions g for the first component u in the
% coupled systems. The functions can be used as a function handle in
% other programs.

function result = g_fun(system,u,v,F)
    if system==0
        result = -u.*v.^2+F*(1-u);
    end   
    if system==1
        result = u.*v +F*(1-v);
    end
    if system==2
        result = v+F*(1-v);
    end
    if system==2.2
        result = -F*u+v;
    end
    if system==3
        result = 0.01*u-F*(1-v);
    end
    if system==4
        result = 0.01*u-F*(1-v).*v;
    end
    if system==5
        result = 0*v;
    end
    if system==6 % ex2
        result = -F*u;
    end
    if system==7 % ex3
        result = -F*u.^2;
    end
    if system==8 % ex4
        result = -F*u.^3;
    end
    if system==9
        result = -u.*v+F*(1-u);
    end
    
end