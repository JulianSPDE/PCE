% this function takes two arrays like in a MC vs PCE simulation, 
% And returns the L2 error vector, one entry per 'slice' or snapshot.
function vec = geterrorvec(A,B,xL,p,snapshots,dimension,relative)
    vec = zeros(1,snapshots);
    xx = linspace(xL(1),xL(2),p);
    fprintf("Computing error vector.\n")
    if dimension==1
        A = A(1:p,:); % We need only the first function (in the case of a PCE import)
        B = B(1:p,:);
        for i=1:snapshots
            %disp(size(xx))
            %fprintf("A:\n")
            %disp(size(A(:,i)))
            %fprintf("B:\n")
            %disp(size(B(:,i)))
            %figure(3)
            %clf
            %hold on
            %plot(xx,A(:,i))
            %plot(xx,B(:,i))
            %hold off
            
            %disp(A(:,i)'-B(:,i)')
            %v1 = A(:,i)'-cos(pi*xx);
            %v2 = B(:,i)'-cos(pi*xx);
            %fprintf("initial condition errors:\n")
           %disp(norm(v1,2))
            %disp(norm(v2,2))
            %disp()
            %disp(norm(A(:,i)-B(:,i),2))
            vec(i) = sqrt(trapz(xx,(abs(A(:,i)-B(:,i))).^2,1));
            if vec(i)>10000
                vec(i)=10000;
            end
            if relative==1
                vec(i) = vec(i)/sqrt(trapz(xx,(abs(A(:,i))).^2,1));
            end
        end
    end
    if dimension==2
        A = A(1:p,1:p,:);
        B = B(1:p,1:p,:);
        for i=1:snapshots
            vec(i) = sqrt(trapz(xx,trapz(xx,abs(A(:,:,i)-B(:,:,i)).^2,2),1));
            if vec(i)>10000
                vec(i)=10000;
            end
            if relative==1
                vec(i) = vec(i)/sqrt(trapz(xx,trapz(xx,abs(A(:,:,i)).^2,2),1));
            end
        end
    end
end