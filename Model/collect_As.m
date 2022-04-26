% collect a lot of As
clear all 
N = 1500;
mu = 15;
L=ceil(mu*N/2);
rho0 = 0.01;
% rho0 = 0.005;
rhosi0 = rho0 * mu;
% SAM = 10000;
SAM = 500;
ii = 0;
for s = 1:SAM
    disp(s);
    A = ERG(N,L);
    x = InitState(rho0,N);
    si0 = ((x==0)*A)*(x==1)'/N;
    i0 = sum(x)/N;
    disp([rho0,i0])
    disp([rhosi0,si0])
    if rhosi0 == si0
        disp('got one')
        disp(ii);
        ii = ii+1;
        if ii == 1
            As = zeros(N,N,1);
            xs = zeros(N,1);
            As(:,:,1) = A;
            xs(:,1) = x;
        else
            B = zeros(N,N,1);
            ys = zeros(N,1);
            B(:,:,1) = A;
            ys(:,1) = x;
            xs = cat(2,xs,ys);
            As = cat(3,As,B);
        end
        disp([rhosi0,si0]);
    end
end

sizeAs = size(As);
for jj = 1:sizeAs(3)
    A = As(:,:,jj);
    x = xs(:,jj)';
    si0 = ((x==0)*A)*(x==1)'/N;
    i0 = sum(x)/N;
    disp([i0,si0]);
end
% 
% save('As_example.mat','As','xs','rho0','mu','N','jj','SAM')
% save('As_example_1000.mat','As','xs','rho0','mu','N','jj','SAM')
save('As_example_1500.mat','As','xs','rho0','mu','N','jj','SAM')
