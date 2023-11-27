function F=Increment_Explicit_Euler_2D(dt,w,p,A,D_u,Dv,f,k,system)
    u=w(1:p^2); v=w(p^2+1:end);
    
   % System 0 (Gray-Scott)
   if system == 0
       F=dt*[((D_u)*A*u) + (-u.*v.^2 +f*(1-u) ) ;
       ((Dv)*A*v) + ( u.*v.^2 -(f+k)*v ) ];
   end
   
   % System 1
   if system == 1
        F=dt*[((D_u)*A*u) + (+v.^1.*u.^1 +f*(1-v) ) ;
         ((Dv)*A*v) + ( -u.^1.*v.^1 -(f+k)*u ) ];
   end
   % System 2
   if system == 2
        F=dt*[((D_u)*A*u) + (+v.^1.*u.^0 +f*(1-v) ) ;
        ((Dv)*A*v) + ( -u.^1.*v.^0 -(f+k)*u ) ];
   end
   % System 3 
   if system == 3
        F=dt*[((D_u)*A*u) + (+0.01*u.^1.*u.^0 -f*(1-v.^1).*v ) ;
        ((Dv)*A*v) + ( -v.^1.*v.^0 +(f+k)*u.^1) ];   
   end
   % System 4 
   if system == 4
        F=dt*[((D_u)*A*u) + (+0.01*u.^1.*u.^0 -f*(1-v.^1).*v ) ;
        ((Dv)*A*v) + ( -v.^1.*v.^0 +(f+k)*u.^2) ];
   end
   % System 5 
   if system == 5
        F=dt*[((D_u)*A*u);
        ((Dv)*A*v)];
   end
   % System 6
   if system == 6
        F=dt*[((D_u)*A*u - f*u);
        0*v];
   end
   % 7
   if system == 7
        F=dt*[((D_u)*A*u - f*u.*u);
        0*v];
   end
   % 8
   if system == 8
        F=dt*[((D_u)*A*u - f*u.*u.*u);
        0*v];
   end
end