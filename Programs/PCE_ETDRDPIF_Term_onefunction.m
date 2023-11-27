
% Input matrices: M_1 = I+k/3A, M_2 = I+kA with A=Laplacian matrix 1D
                 
function U_new = PCE_ETDRDPIF_Term_onefunction(U_temp,dt,N,p,polynomialtype,L1,U1,L2,U2,L3,U3,L11,U11,L22,U22,L33,U33,F,K,Fvec,dimension,system)
        LL1 = kron(speye(N+1),L1);   UU1 = kron(speye(N+1),U1);
        LL2 = kron(speye(N+1),L2);   UU2 = kron(speye(N+1),U2);
        LL3 = kron(speye(N+1),L3);   UU3 = kron(speye(N+1),U3);
        LL11 = kron(speye(N+1),L11); UU11 = kron(speye(N+1),U11);
        LL22 = kron(speye(N+1),L22); UU22 = kron(speye(N+1),U22);
        LL33 = kron(speye(N+1),L33); UU33 = kron(speye(N+1),U33);
        if dimension==1
            if system == 6
                f_U = -kron(F,speye(p))*U_temp;
            end
            if system == 7
                f_U = -kron(F,speye(p))*u2(U_temp,N,polynomialtype);
            end     
            if system == 8
                f_U = -kron(F,speye(p))*u3(U_temp,N,p,polynomialtype);
            end      
            star_rhs = U_temp + dt*f_U;
            U_star = UU1\(LL1\star_rhs);
            if system == 6
                f_U_star = -kron(F,speye(p))*U_star;
            end
            if system == 7
                f_U_star = -kron(F,speye(p))*u2(U_star,N,polynomialtype);
            end     
            if system == 8
                f_U_star = -kron(F,speye(p))*u3(U_star,N,p,polynomialtype);
            end      
            rhs_1 = 9*U_temp+2*dt*f_U+dt*f_U_star;
            U_new_1 = UU2\(LL2\rhs_1);
            rhs_2 = 8*U_temp+1.5*dt*f_U+0.5*dt*f_U_star;
            U_new_2 = UU3\(LL3\rhs_2);
            U_new = U_new_1 - U_new_2;
        end
        if dimension==2
            % Need other matrices... LL11, UU11 etc.
            if system == 6
                f_U = -kron(F,speye(p^2))*U_temp;
            end
            if system == 7
                f_U = -kron(F,speye(p^2))*u2(U_temp,N,polynomialtype);
            end     
            if system == 8
                f_U = -kron(F,speye(p^2))*u3(U_temp,N,p,polynomialtype);
            end
            star_rhs = U_temp + dt*f_U;
            U_star = UU1\(LL1\star_rhs);
            U_star = UU11\(LL11\U_star);
            if system == 6
                f_U_star = -kron(F,speye(p^2))*U_star;
            end
            if system == 7
                f_U_star = -kron(F,speye(p^2))*u2(U_star,N,polynomialtype);
            end     
            if system == 8
                f_U_star = -kron(F,speye(p^2))*u3(U_star,N,p,polynomialtype);
            end      
            a_1 = UU2\(LL2\U_temp);
            b_1 = UU3\(LL3\U_temp);
            c_1 = 9*a_1 - 8*b_1;
            a_2 = UU2\(LL2\f_U);
            b_2 = UU3\(LL3\f_U);
            c_2 = 9*a_2-8*b_2;
            d_1 = UU22\(LL22\(9*c_1+2*dt*c_2 + dt*f_U_star));
            d_2 = UU33\(LL33\(-8*c_1-(3/2)*dt*c_2-(dt/2)*f_U_star));               
            U_new = d_1 + d_2;
        end
end