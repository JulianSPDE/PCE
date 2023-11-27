function F=Increment_Explicit_Euler_2D_onefunction(dt,w,p,A,D_u,Dv,f,k,system)
    u=w(1:p^2);
    
   % System 6
   if system == 6
        F=dt*((D_u)*A*u - f*u);
        M1 = D_u*A; %disp(M1(1:10,1:10))
        %{
        fprintf("F:\n")
        disp(f)
        disp(dt)
        disp(size(u))
        %}
        %{
        MM = dt*((D_u)*A*u - f*u);
        disp(size(MM))
        disp(MM(1:10))
        figure(17)
        imagesc(linspace(-1,1,128),linspace(-1,1,128),reshape(MM,128,128))
        colorbar
        figure(1)
        %}
   end
   % 7
   if system == 7
        F=dt*((D_u)*A*u - f*u.*u);
   end
   % 8
   if system == 8
        F=dt*((D_u)*A*u - f*u.*u.*u);
   end
end