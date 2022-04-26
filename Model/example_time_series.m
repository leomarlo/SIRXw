

load('As_example_1000.mat')
sizaAs = size(As);

tmax = 1000;
Nt = tmax;
tspan = [0:1:tmax];
gamma = 1/40;
delta = 1/10;
wplk = 0.001;
w = wplk/2;
k = wplk - w;

mu = 15;
beta = 0.004;

rhosi0 = rho0 * mu;
rhoss0 = (mu/2)- rhosi0;
ini=[1-rho0,rho0,0,rhosi0,rhoss0];
SA = 15;

ODE = @(x) PA_SIRXi_w(x,beta,gamma,w,k,delta);
[ts,ys] = ode45(@(t,x) ODE(x),tspan, ini);

mfs = zeros(Nt+1,2, SA);
sa = 1;
for sam = randsample(1:sizaAs(3),SA)
    A = As(:,:,sa);
    x = xs(:,sa);
    
    si0 = (x==0)'*A*x/N;
    disp(si0)
    result = SIRXi_w_fromAandx(A,x',beta,gamma,w,k,delta,tmax,Nt); 
    mfs(:,1,sa) = result.NIs/result.N;
    mfs(:,2,sa) = result.NRs/result.N;
    sa = sa+1;
end



pl = 1;
clear('pl');
figure; 
sa= 1;
sas = zeros(1,5);
for  sam = [12,     3 ,    7 ,    4,    15]
    pl(4) = plot(tspan, mfs(:,2,sam));
    pl(4).Color = [0.25,0.25,0.75];
    pl(4).LineWidth = 0.5;
    pl(4).DisplayName = '\rho_R Simulation';
    hold on;
    
    pl(2) = plot(tspan, mfs(:,1,sam));
    pl(2).Color = [0.75,0.25,0.25];
    pl(2).LineWidth = 0.5;
    pl(2).DisplayName = '\rho_I Simulation';
    
    sas(sa) = sam;
    sa = sa+1;
end
pl(3) = plot(tspan, ys(:,3));
pl(3).Color = [0,0,1];
pl(3).LineWidth = 3;
pl(3).DisplayName = '\rho_R Pair Approximation';
pl(1) = plot(tspan, ys(:,2));
pl(1).Color = [1,0,0];
pl(1).LineWidth = 3;
pl(1).DisplayName = '\rho_I Pair Approximation';

lgd = legend(pl);
lgd.Location = 'northwest';

ax = gca;
ax.XLim = [0,500];

folder='figures/';
filename=strcat('examples_plots_10Sept');
direction=strcat(folder,filename,'.png');
saveas(gcf,direction)


