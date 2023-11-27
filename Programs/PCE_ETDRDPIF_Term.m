
% Input matrices: M_1 = I+k/3A, M_2 = I+kA with A=Laplacian matrix 1D

function [U_new,V_new] = PCE_ETDRDPIF_Term(U_temp,V_temp,dt,N,p,polynomialtype,L1,U1,L2,U2,L3,U3,L11,U11,L22,U22,L33,U33,F,K,Fvec,dimension,system)

        %% One spatial dimension (ETDRDP, without dimensional splitting)
        if dimension==1  % Not sure how it would look for dimension 2, would have to work this out separately
            LL1 = kron(speye(N+1),L1);  UU1 = kron(speye(N+1),U1);
            LL2 = kron(speye(N+1),L2);  UU2 = kron(speye(N+1),U2);
            LL3 = kron(speye(N+1),L3);  UU3 = kron(speye(N+1),U3);
            LL1_U = LL1(1:(N+1)*p,1:(N+1)*p); LL1_V = LL1((N+1)*p+1:end,(N+1)*p+1:end); 
            UU1_U = UU1(1:(N+1)*p,1:(N+1)*p); UU1_V = UU1((N+1)*p+1:end,(N+1)*p+1:end); 
            LL2_U = LL2(1:(N+1)*p,1:(N+1)*p); LL2_V = LL2((N+1)*p+1:end,(N+1)*p+1:end); 
            UU2_U = UU2(1:(N+1)*p,1:(N+1)*p); UU2_V = UU2((N+1)*p+1:end,(N+1)*p+1:end); 
            LL3_U = LL3(1:(N+1)*p,1:(N+1)*p); LL3_V = LL3((N+1)*p+1:end,(N+1)*p+1:end); 
            UU3_U = UU3(1:(N+1)*p,1:(N+1)*p); UU3_V = UU3((N+1)*p+1:end,(N+1)*p+1:end); 
            W_temp = [U_temp;V_temp];
            if system == 0
                f_U = -uv2(U_temp,V_temp,N,polynomialtype)+ kron(Fvec,ones(p^dimension,1)) -kron(F,speye(p))*U_temp;
                f_V = uv2(U_temp,V_temp,N,polynomialtype)-kron(F+K,speye(p))*V_temp;
            end
            if system == 1
                f_U = uv(U_temp,V_temp,N,polynomialtype)+ kron(Fvec,ones(p^dimension,1))-kron(F,speye(p))*V_temp;
                f_V = -uv(U_temp,V_temp,N,polynomialtype)-kron(F+K,speye(p))*U_temp;
            end
            if system == 2
                f_U = V_temp+ kron(Fvec,ones(p^dimension,1))-kron(F,speye(p))*V_temp;
                f_V = -U_temp-kron(F+K,speye(p))*U_temp;
            end
            if system == 3
                f_U = 0.01*U_temp- kron(Fvec,ones(p^dimension,1))+kron(F,speye(p))*u2(V_temp,N,p,polynomialtype);
                f_V = -V_temp+kron(F+K,speye(p))*U_temp;
            end
            if system == 4
                f_U = 0.01*U_temp- kron(Fvec,ones(p^dimension,1))+kron(F,speye(p))*u2(V_temp,N,p,polynomialtype);
                f_V = -V_temp+kron(F+K,speye(p))*u2(U_temp,N,p,polynomialtype);
            end
            if system == 5
                f_U = 0*U_temp;
                f_V = 0*V_temp;
            end   
            if system == 9
                f_U = uv(U_temp,V_temp,N,polynomialtype)+Fvec-kron(F,speye(p))*V_temp;
                f_V = -uv(U_temp,V_temp,N,polynomialtype)-kron(F+K,speye(p))*U_temp;
            end
            f_W = [f_U;f_V];
            star_rhs = W_temp + dt*f_W;
            W_star = UU1\(LL1\star_rhs);
            U_star = W_star(1:p*(N+1)); V_star = W_star(p*(N+1)+1:end);
           % star_rhs_U = U_temp + dt*f_U; star_rhs_V = V_temp + dt*f_V;
           % U_star = UU1_U\(LL1_U\star_rhs_U); V_star = UU1_V\(LL1_V\star_rhs_V);
            if system == 0
                f_U_star = -uv2(U_star,V_star,N,polynomialtype)+ kron(Fvec,ones(p^dimension,1))-kron(F,speye(p))*U_star;
                f_V_star = uv2(U_star,V_star,N,polynomialtype)-kron(F+K,speye(p))*V_star;
            end
            if system == 1
                f_U_star = uv(U_star,V_star,N,polynomialtype)+ kron(Fvec,ones(p^dimension,1))-kron(F,speye(p))*V_star;
                f_V_star = -uv(U_star,V_star,N,polynomialtype)-kron(F+K,speye(p))*U_star;
            end
            if system == 2
                f_U_star = V_star+ kron(Fvec,ones(p^dimension,1))-kron(F,speye(p))*V_star;
                f_V_star = -U_star-kron(F+K,speye(p))*U_star;
            end
            if system == 3
                f_U_star = 0.01*U_star- kron(Fvec,ones(p^dimension,1))+kron(F,speye(p))*u2(V_star,N,p,polynomialtype);
                f_V_star = -V_star+kron(F+K,speye(p))*U_star;
            end
            if system == 4
                f_U_star = 0.01*U_star- kron(Fvec,ones(p^dimension,1))+kron(F,speye(p))*u2(V_star,N,p,polynomialtype);
                f_V_star = -V_star+kron(F+K,speye(p))*u2(U_star,N,p,polynomialtype);
            end
            if system == 5
                f_U_star = 0*U_star;
                f_V_star = 0*V_star;
            end     
            if system == 9
                f_U_star = uv(U_star,V_star,N,polynomialtype)+ kron(Fvec,ones(p^dimension,1))-kron(F,speye(p))*V_star;
                f_V_star = -uv(U_star,V_star,N,polynomialtype)-kron(F+K,speye(p))*U_star;
            end      
        f_W_star = [f_U_star;f_V_star];
        rhs_1 = 9*W_temp+2*dt*f_W+dt*f_W_star;
        W_new_1 = UU2\(LL2\rhs_1);
        rhs_2 = 8*W_temp+1.5*dt*f_W+0.5*dt*f_W_star;
        W_new_2 = UU3\(LL3\rhs_2);
        W_new = W_new_1 - W_new_2;
        U_new = W_new(1:p*(N+1)); V_new = W_new(p*(N+1)+1:end);
        %rhs_1_U = 9*U_temp+2*dt*f_U+dt*f_U_star; rhs_1_V = 9*V_temp+2*dt*f_V+dt*f_V_star;
        %U_new_1 = UU2_U\(LL2_U\rhs_1_U); V_new_1 = UU2_V\(LL2_V\rhs_1_V);
        %rhs_2_U = 8*U_temp+1.5*dt*f_U+0.5*dt*f_U_star; rhs_2_V = 8*V_temp+1.5*dt*f_V+0.5*dt*f_V_star;
        %U_new_2 = UU3_U\(LL3_U\rhs_2_U); V_new_2 = UU3_V\(LL3_V\rhs_2_V);
        %U_new = U_new_1 - U_new_2; V_new = V_new_1 - V_new_2;
        end
        


        %% Two spatial dimensions (ETD-RDP-IF, with dimensional splitting)
        if dimension==2  % Not sure how it would look for dimension 2, would have to work this out separately
            LL1 = kron(speye(N+1),L1);  UU1 = kron(speye(N+1),U1);
            LL2 = kron(speye(N+1),L2);  UU2 = kron(speye(N+1),U2);
            LL3 = kron(speye(N+1),L3);  UU3 = kron(speye(N+1),U3);
            LL11 = kron(speye(N+1),L11); UU11 = kron(speye(N+1),U11);
            LL22 = kron(speye(N+1),L22); UU22 = kron(speye(N+1),U22);
            LL33 = kron(speye(N+1),L33); UU33 = kron(speye(N+1),U33);
            W_temp = [U_temp;V_temp];
            if system == 0
                f_U = uv2(U_temp,V_temp,N,polynomialtype)+kron(Fvec,ones(p^dimension,1))-kron(F,speye(p^2))*V_temp;
                f_V = -uv2(U_temp,V_temp,N,polynomialtype)-kron(F+K,speye(p^2))*U_temp;
            end
            if system == 1
                f_U = uv(U_temp,V_temp,N,polynomialtype)+kron(Fvec,ones(p^dimension,1))-kron(F,speye(p^2))*V_temp;
                f_V = -uv(U_temp,V_temp,N,polynomialtype)-kron(F+K,speye(p^2))*U_temp;
            end
            if system == 2
                f_U = V_temp+kron(Fvec,ones(p^dimension,1))-kron(F,speye(p^2))*V_temp;
                f_V = -U_temp-kron(F+K,speye(p^2))*U_temp;
            end
            if system == 3
                f_U = 0.01*U_temp-kron(Fvec,ones(p^dimension,1))+kron(F,speye(p^2))*u2(V_temp,N,p,polynomialtype);
                f_V = -V_temp+kron(F+K,speye(p^2))*U_temp;
            end
            if system == 4
                f_U = 0.01*U_temp-kron(Fvec,ones(p^dimension,1))+kron(F,speye(p^2))*u2(V_temp,N,p,polynomialtype);
                f_V = -V_temp+kron(F+K,speye(p^2))*u2(U_temp,N,p,polynomialtype);
            end
            if system == 5
                f_U = 0*U_temp;
                f_V = 0*V_temp;
            end
            if system == 6
                f_U = -kron(F,speye(p^2))*U_temp;
                f_V = 0*V_temp;
            end
            if system == 7
                f_U = -kron(F,speye(p^2))*u2(U_temp,N,polynomialtype);
                f_V = 0*V_temp;
            end     
            if system == 8
                f_U = -kron(F,speye(p^2))*u3(U_temp,N,p,polynomialtype);
                f_V = 0*V_temp;
            end      
            if system == 9
                f_U = uv(U_temp,V_temp,N,polynomialtype)+kron(Fvec,ones(p^dimension,1))-kron(F,speye(p^2))*V_temp;
                f_V = -uv(U_temp,V_temp,N,polynomialtype)-kron(F+K,speye(p^2))*U_temp;
            end
            f_W = [f_U;f_V];
            star_rhs = W_temp + dt*f_W;
            W_star = UU1\(LL1\star_rhs); W_star = UU11\(LL11\W_star);
            U_star = W_star(1:p^2*(N+1)); V_star = W_star(p^2*(N+1)+1:end);
            %star_rhs_U = U_temp + dt*f_U; star_rhs_V = V_temp + dt*f_V;
            %U_star = UU1\(LL1\star_rhs_U); U_star = UU11\(LL11\U_star);
            %V_star = UU1\(LL1\star_rhs_V); V_star = UU11\(LL11\V_star);
            if system == 0
                f_U_star = uv2(U_star,V_star,N,polynomialtype)+kron(Fvec,ones(p^dimension,1))-kron(F,speye(p^2))*V_star;
                f_V_star = -uv2(U_star,V_star,N,polynomialtype)-kron(F+K,speye(p^2))*U_star;
            end
            if system == 1
                f_U_star = uv(U_star,V_star,N,polynomialtype)+kron(Fvec,ones(p^dimension,1))-kron(F,speye(p^2))*V_star;
                f_V_star = -uv(U_star,V_star,N,polynomialtype)-kron(F+K,speye(p^2))*U_star;
            end
            if system == 2
                f_U_star = V_star+kron(Fvec,ones(p^dimension,1))-kron(F,speye(p^2))*V_star;
                f_V_star = -U_star-kron(F+K,speye(p^2))*U_star;
            end
            if system == 3
                f_U_star = 0.01*U_star-kron(Fvec,ones(p^dimension,1))+kron(F,speye(p^2))*u2(V_star,N,p,polynomialtype);
                f_V_star = -V_star+kron(F+K,speye(p^2))*U_star;
            end
            if system == 4
                f_U_star = 0.01*U_star-kron(Fvec,ones(p^dimension,1))+kron(F,speye(p^2))*u2(V_star,N,p,polynomialtype);
                f_V_star = -V_star+kron(F+K,speye(p^2))*u2(U_star,N,p,polynomialtype);
            end
            if system == 5
                f_U_star = 0*U_star;
                f_V_star = 0*V_star;
            end
            if system == 6
                f_U_star = -kron(F,speye(p^2))*U_star;
                f_V_star = 0*V_star;
            end
            if system == 7
                f_U_star = -kron(F,speye(p^2))*u2(U_star,N,polynomialtype);
                f_V_star = 0*V_star;
            end     
            if system == 8
                f_U_star = -kron(F,speye(p^2))*u3(U_star,N,p,polynomialtype);
                f_V_star = 0*V_star;
            end      
            if system == 9
                f_U_star = uv(U_star,V_star,N,polynomialtype)+kron(Fvec,ones(p^dimension,1))-kron(F,speye(p))*V_star;
                f_V_star = -uv(U_star,V_star,N,polynomialtype)-kron(F+K,speye(p))*U_star;
            end      
        %{
        a_U1 = UU2\(LL2\U_temp);  a_V1 = UU2\(LL2\V_temp);
        b_U1 = UU3\(LL3\U_temp);  b_V1 = UU3\(LL3\V_temp);
        c_U1 = 9*a_U1 - 8*b_U1;   c_V1 = 9*a_V1 - 8*b_V1;
        a_U2 = UU2\(LL2\f_U);     a_V2 = UU2\(LL2\f_V);  
        b_U2 = UU3\(LL3\f_U);     b_V2 = UU3\(LL3\f_V);   
        c_U2 = 9*a_U2-8*b_U2;     c_V2 = 9*a_V2-8*b_V2;
        d_U1 = UU22\(LL22\(9*c_U1+2*dt*c_U2 + dt*f_U_star));
        d_V1 = UU22\(LL22\(9*c_V1+2*dt*c_V2 + dt*f_V_star));
        d_U2 = UU33\(LL33\(-8*c_U1-(3/2)*dt*c_U2-(dt/2)*f_U_star));  
        d_V2 = UU33\(LL33\(-8*c_V1-(3/2)*dt*c_V2-(dt/2)*f_V_star));  
        U_new = d_U1 + d_U2; V_new = d_V1 + d_V2;
        %}
        f_W_star = [f_U_star;f_V_star];
        a_1 = UU2\(LL2\W_temp);
        b_1 = UU3\(LL3\W_temp);
        c_1 = 9*a_1 - 8*b_1;
        a_2 = UU2\(LL2\f_W);
        b_2 = UU3\(LL3\f_W);
        c_2 = 9*a_2-8*b_2;
        d_1 = UU22\(LL22\(9*c_1+2*dt*c_2 + dt*f_W_star));
        d_2 = UU33\(LL33\(-8*c_1-(3/2)*dt*c_2-(dt/2)*f_W_star));   
        W_new = d_1 + d_2;
        U_new = W_new(1:p^2*(N+1)); V_new = W_new(p^2*(N+1)+1:end);
        end

end