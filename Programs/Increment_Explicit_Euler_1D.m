function F=Increment_Explicit_Euler_1D(dt,w,p,A,D_u,D_v,f,k,system)
    u=w(1:p); v=w(p+1:end);
    
   % System 0 (Gray-Scott)
   if system == 0
       F=dt*[((D_u)*A*u) + (-u.*v.^2 +f*(1-u) ) ;
       ((D_v)*A*v) + ( u.*v.^2 -(f+k)*v ) ];
   end
   
   % System 1
   if system == 1
        F=dt*[((D_u)*A*u) + (+v.^1.*u.^1 +f*(1-v) ) ;
         ((D_v)*A*v) + ( -u.^1.*v.^1 -(f+k)*u ) ];
   end
   % System 2
   if system == 2
        F=dt*[((D_u)*A*u) + (+v.^1.*u.^0 +f*(1-v) ) ;
        ((D_v)*A*v) + ( -u.^1.*v.^0 -(f+k)*u ) ];
   end
   % System 3 
   if system == 3
        F=dt*[((D_u)*A*u) + (+0.01*u.^1.*u.^0 -f*(1-v.^1).*v ) ;
        ((D_v)*A*v) + ( -v.^1.*v.^0 +(f+k)*u.^1) ];   
   end
   % System 4 
   if system == 4
        F=dt*[((D_u)*A*u) + (+0.01*u.^1.*u.^0 -f*(1-v.^1).*v ) ;
        ((D_v)*A*v) + ( -v.^1.*v.^0 +(f+k)*u.^2) ];
   end
   % System 5 
   if system == 5
        F=dt*[((D_u)*A*u);
        ((D_v)*A*v)];
   end
   
   if system == 6
       F=dt*[D_u*A*u - f*u;
           zeros(p,1)];
   end

   if system == 7
       F=dt*[D_u*A*u - f*u.^2;
           zeros(p,1)];
   end
   if system == 8
       F=dt*[D_u*A*u - f*u.^3;
           zeros(p,1)];
   end
   if system == 9
       F=dt*[((D_u)*A*u) + (-u.*v +f*(1-u) ) ;
       ((D_v)*A*v) + ( u.*v -(f+k)*v ) ];
   end
end