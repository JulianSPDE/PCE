% Testing how PCE runtimes scale with increasing N

Nmax = 9;
runtimevec = zeros(6,Nmax+1,3); % Six systems, Nmax+1 cases for N, three schemes
avg = 20; % Do this number of runs and average runtimes over them
spatialvector = ["FD","FD","Spectral"];
timesteppingvector = ["EE","ETDRDPIF","ETDRK4"];
M = 10;
for scheme=1:3
    spatial = spatialvector(scheme);
    timestep = timesteppingvector(scheme);
    for i=1:Nmax+1
       runtimetemp = zeros(1,6);
       for j=1:avg
           tic; PCE_solver(spatial,timestep,'mean',128,[-1,1],0.1,M,i-1,0,0,0,0,'F',1,2,6,'Legendre',1,1,0,1); 
           runtimetemp(1) = runtimetemp(1) + toc;
           tic; PCE_solver(spatial,timestep,'mean',128,[-1,1],0.1,M,i-1,1,0,0,0,'F',1,2,6,'Legendre',1,1,0,1); 
           runtimetemp(2) = runtimetemp(2) + toc;
           tic; PCE_solver(spatial,timestep,'mean',128,[-1,1],0.1,M,i-1,0,0,0,0,'F',1,2,7,'Legendre',1,1,0,1); 
           runtimetemp(3) = runtimetemp(3) + toc;
           tic; PCE_solver(spatial,timestep,'mean',128,[-1,1],0.1,M,i-1,1,0,0,0,'F',1,2,7,'Legendre',1,1,0,1); 
           runtimetemp(4) = runtimetemp(4) + toc;
           tic; PCE_solver(spatial,timestep,'mean',128,[-1,1],0.1,M,i-1,0,0,0,0,'F',1,2,8,'Legendre',1,1,0,1); 
           runtimetemp(5) = runtimetemp(5) + toc;
           tic; PCE_solver(spatial,timestep,'mean',128,[-1,1],0.1,M,i-1,1,0,0,0,'F',1,2,8,'Legendre',1,1,0,1); 
           runtimetemp(6) = runtimetemp(6) + toc;
       end
       for j=1:6
        runtimevec(j,i,scheme) = runtimetemp(j)/avg;
       end
    end
end

% Norming
for i=1:6
    for j=1:3
        runtimevec(i,:,j) = runtimevec(i,:,j)/runtimevec(i,1,j);
    end
end

figure(1)
title('EE')
hold on
for i=1:6
    plot(0:Nmax,runtimevec(i,:,1))
end
hold off
legend("6 D=0","6 D=1","7 D=0","7 D=1","8 D=0","8 D=1")

figure(2)
title('ETD-RDP')
hold on
for i=1:6
    plot(0:Nmax,runtimevec(i,:,2))
end
legend("6 D=0","6 D=1","7 D=0","7 D=1","8 D=0","8 D=1")


figure(3)
title('ETDRK4')
hold on
for i=1:6
    plot(0:Nmax,runtimevec(i,:,3))
end
legend("6 D=0","6 D=1","7 D=0","7 D=1","8 D=0","8 D=1")

runtimearray_EE = zeros(Nmax+1,7);
runtimearray_ETDRDP = zeros(Nmax+1,7);
runtimearray_ETDRK4 = zeros(Nmax+1,7);

runtimearray_EE(:,1) = (0:Nmax)';
runtimearray_ETDRDP(:,1) = (0:Nmax)';
runtimearray_ETDRK4(:,1) = (0:Nmax)';

for i=2:7
    runtimearray_EE(:,i) = runtimevec(i-1,:,1)';
    runtimearray_ETDRDP(:,i) = runtimevec(i-1,:,2)';
    runtimearray_ETDRK4(:,i) = runtimevec(i-1,:,3)';
end

writematrix(runtimearray_EE,sprintf('runtimearray_EE.txt'));
writematrix(runtimearray_ETDRDP,sprintf('runtimearray_ETDRDP.txt'));
writematrix(runtimearray_ETDRK4,sprintf('runtimearray_ETDRK4.txt'));

% save('PCE_runtimes.mat','runtimevec')


