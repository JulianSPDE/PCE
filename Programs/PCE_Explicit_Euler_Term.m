
function [U_new,V_new] = PCE_Explicit_Euler_Term(U_temp,V_temp,dt,N,p,polynomialtype,varmin,varmax,M_U1,M_U2,M_V1,M_V2,Fvec,dimension,system)
        

        if system == 0
            [N_UV,~,~] = uv2(U_temp,V_temp,N,polynomialtype);
            U_new = U_temp + dt*(M_U1*U_temp + kron(Fvec,ones(p^dimension,1)) - M_U2*U_temp - N_UV);
            V_new = V_temp + dt*(M_V1*V_temp                          - M_V2*V_temp + N_UV);        
        end
        % System 1
        if system == 1
            UV = uv(U_temp,V_temp,N,polynomialtype);
            U_new = U_temp + dt*(M_U1*U_temp + kron(Fvec,ones(p^dimension,1)) - M_U2*V_temp + UV);
            V_new = V_temp + dt*(M_V1*V_temp                          - M_V2*U_temp - UV);
        end
        % System 2
        if system == 2
            U_new = U_temp + dt*(M_U1*U_temp + kron(Fvec,ones(p^dimension,1)) - M_U2*V_temp + V_temp);
            V_new = V_temp + dt*(M_V1*V_temp                          - M_V2*U_temp - U_temp);
        end
        % System 3
        if system == 3
            V2 = u2(V_temp,N,polynomialtype);
            U_new = U_temp + dt*(M_U1*U_temp                          - M_U2*(V_temp-V2) + 0.01*U_temp);
            V_new = V_temp + dt*(M_V1*V_temp                          + M_V2*U_temp - V_temp);
        end
        % System 4
        if system == 4
            U2 = u2(U_temp,N,polynomialtype);
            V2 = u2(V_temp,N,polynomialtype);
            U_new = U_temp + dt*(M_U1*U_temp                          - M_U2*(V_temp-V2) + 0.01*U_temp);
            V_new = V_temp + dt*(M_V1*V_temp                          + M_V2*U2 - V_temp);
        end
        % System 5
        if system == 5
            U_new = U_temp + dt*(M_U1*U_temp);
            V_new = V_temp + dt*(M_V1*V_temp);
        end
        % System 6
        if system == 6
            U_new = U_temp + dt*(M_U1*U_temp - M_U2*U_temp);
            V_new = 0*V_temp;
            %{
            fprintf("1st:\n")
            disp(M_U1(1:10,1:10))
            fprintf("2nd:\n")
            disp(M_U2(1:10,1:10))
            disp(dt)
            disp(size(U_temp))
            %}
            %{
            MM = dt*(M_U1*U_temp - M_U2*U_temp);
            save('MM_2.mat','MM')
            save('M_U1_2.mat','M_U1')
            save('M_U2_2.mat','M_U2')
            figure(18)
            imagesc(linspace(-1,1,128),linspace(-1,1,128),reshape(MM,128,128))
            colorbar
            figure(2)
            %}
        end
        % System 7
        if system == 7
            U2 = u2(U_temp,N,polynomialtype);
            U_new = U_temp + dt*(M_U1*U_temp - M_U2*U2);
            V_new = 0*V_temp;
        end
        % System 8        
        if system == 8
            % U2 = u2(U_temp,N,polynomialtype);
            % U3 = uv(U_temp,U2,N,polynomialtype);
            U3 = u3(U_temp,N,p,polynomialtype);
            U_new = U_temp + dt*(M_U1*U_temp -  M_U2*U3);
            V_new = 0*V_temp;
        end
        % System 9        
        if system == 9
            UV = uv(U_temp,V_temp,N,polynomialtype);
            U_new = U_temp + dt*(M_U1*U_temp + kron(Fvec,ones(p^dimension,1)) - M_U2*U_temp - UV);
            V_new = V_temp + dt*(M_V1*V_temp                          - M_V2*V_temp + UV);        
        end


end