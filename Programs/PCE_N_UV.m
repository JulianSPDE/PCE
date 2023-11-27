


function [N_UV_stacked_1,N_UV_stacked_2] = PCE_N_UV(U_temp,V_temp,Fvec,F,K,N,p,polynomialtype,system,dimension,a,b)
        
        if dimension==1
            if system==0       
                [UV2,~,~] = uv2(U_temp,V_temp,N,polynomialtype);
                N_UV_stacked_1 = -UV2 + kron(Fvec,ones(p,1)) - kron(F,speye(p))*U_temp;
                N_UV_stacked_2 = UV2 - kron(F+K,speye(p))*V_temp;
            end
            if system==1  
                UV = uv(U_temp,V_temp,N,polynomialtype);
                N_UV_stacked_1 = -UV + kron(Fvec,ones(p,1)) - kron(F,speye(p))*V_temp;
                N_UV_stacked_2 = -UV - kron(F+K,speye(p))*U_temp;
            end
            if system==2        
                N_UV_stacked_1 = V_temp + kron(Fvec,ones(p,1)) - kron(F,speye(p))*V_temp;
                N_UV_stacked_2 = -U_temp - kron(F+K,speye(p))*U_temp;
            end
            if system==3        
                N_UV_stacked_1 = 0.01*U_temp + kron(Fvec,ones(p,1)) - kron(F,speye(p))*V_temp;
                N_UV_stacked_2 = -V_temp - kron(F+K,speye(p))*U_temp;
            end
            if system==4   
                V2 = u2(V_temp,N,polynomialtype);
                U2 = u2(U_temp,N,polynomialtype);
                N_UV_stacked_1 = 0.01*U_temp + (kron(Fvec,ones(p,1)).*V_temp - kron(F,speye(p))*V2);
                N_UV_stacked_2 = -V_temp + kron(F+K,speye(p))*U2;
            end
            if system==5        
                N_UV_stacked_1 = zeros(1,(N+1)*p);
                N_UV_stacked_2 = zeros(1,(N+1)*p);
            end
            if system==6        
                N_UV_stacked_1 = zeros(1,(N+1)*p) - kron(F,speye(p))*U_temp;
                N_UV_stacked_2 = zeros(1,(N+1)*p);
            end
            if system==7       
                U2 = u2(U_temp,N,polynomialtype);
                N_UV_stacked_1 = zeros(1,(N+1)*p) - kron(F,speye(p))*U2;
                N_UV_stacked_2 = zeros(1,(N+1)*p);
            end
            if system==8    
                %U2 = u2(U_temp,N,polynomialtype);
                %U3 = uv(U_temp,U2,N,polynomialtype);
                U3 = u3(U_temp,N,p,polynomialtype);
                N_UV_stacked_1 = zeros(1,(N+1)*p) - kron(F,speye(p))*U3;
                N_UV_stacked_2 = zeros(1,(N+1)*p);
            end
            if system==9 
                UV = uv(U_temp,V_temp,N,polynomialtype);
                N_UV_stacked_1 = -UV + kron(Fvec,ones(p,1)) - kron(F,speye(p))*U_temp;
                N_UV_stacked_2 = UV - kron(F+K,speye(p))*V_temp;
            end

        end
        if dimension==2
            U_temp = reshapecompinverse(U_temp,N,p); V_temp = reshapecompinverse(V_temp,N,p);
            if system==0        
                [UV2,~,~] = uv2(U_temp,V_temp,N,polynomialtype);
                N_UV_stacked_1 = reshapecomp(-UV2 + kron(Fvec,sparse(ones(p^2,1))) - kron(F,speye(p^2))*U_temp ,N,p);
                N_UV_stacked_2 = reshapecomp(UV2 - kron(F+K,speye(p^2))*V_temp,N,p);
            end
            if system==1      
                UV = uv(U_temp,V_temp,N,polynomialtype);
                N_UV_stacked_1 = reshapecomp(-UV + kron(Fvec,sparse(ones(p^2,1))) - kron(F,speye(p^2))*V_temp ,N,p);
                N_UV_stacked_2 = reshapecomp(-UV - kron(F+K,speye(p^2))*U_temp,N,p);
            end
            if system==2        
                N_UV_stacked_1 = reshapecomp(V_temp + kron(Fvec,sparse(ones(p^2,1))) - kron(F,speye(p^2))*V_temp ,N,p);
                N_UV_stacked_2 = reshapecomp(-U_temp - kron(F+K,speye(p^2))*U_temp,N,p);
            end
            if system==3        
                N_UV_stacked_1 = reshapecomp(0.01*U_temp + kron(Fvec,sparse(ones(p^2,1))) - kron(F,speye(p^2))*V_temp ,N,p);
                N_UV_stacked_2 = reshapecomp(-V_temp - kron(F+K,speye(p^2))*U_temp,N,p);
            end
            if system==4
                U2 = u2(U_temp,N,polynomialtype);
                V2 = u2(V_temp,N,polynomialtype);
                N_UV_stacked_1 = reshapecomp(0.01*U_temp + (kron(Fvec,sparse(ones(p^2,1)))*V_temp - kron(F,speye(p^2))*V2),N,p);
                N_UV_stacked_2 = reshapecomp(-V_temp + kron(F+K,speye(p^2))*U2,N,p);
            end
            if system==5        
                N_UV_stacked_1 = sparse(zeros(1,(N+1)*p^2));
                N_UV_stacked_2 = sparse(zeros(1,(N+1)*p^2));
            end
            if system==6       
                N_UV_stacked_1 = reshapecomp(-kron(F,speye(p^2))*U_temp,N,p);
                N_UV_stacked_2 = reshapecomp(sparse(zeros(1,(N+1)*p^2)),N,p);
            end
            if system==7       
                U2 = u2(U_temp,N,polynomialtype);
                N_UV_stacked_1 = reshapecomp(-kron(F,speye(p^2))*U2,N,p);
                N_UV_stacked_2 = reshapecomp(sparse(zeros(1,(N+1)*p^2)),N,p);
            end
            if system==8    
                U3 = u3(U_temp,N,p,polynomialtype);
                N_UV_stacked_1 = reshapecomp(-kron(F,speye(p^2))*U3,N,p);
                N_UV_stacked_2 = reshapecomp(sparse(zeros(1,(N+1)*p^2),N,p));
            end
            if system==9 
                UV = uv(U_temp,V_temp,N,polynomialtype);
                N_UV_stacked_1 = reshapecomp(-UV + kron(Fvec,ones(p^2,1)) - kron(F,speye(p^2))*U_temp,N,p);
                N_UV_stacked_2 = reshapecomp(UV - kron(F+K,speye(p^2))*V_temp,N,p);
            end            
        end
end





