%%% compare the effect of combinations

tmax=80;
Nt=tmax*2;

N=100;
mu=25;
beta=1/100;
gamma=1/40;
w=1/100;
rho0=1/100;
kap=12/10000;
delta=2/100;
% rhos,rhoi,rhor,rhosi,rhoss
rhoi=0.01;
rhos=1-rhoi;
rhor=0;
rhosi= rhoi*mu;
rhoss= (mu/2)-rho0*mu;
ini=[rhos,rhoi,rhor,rhosi,rhoss];

ODE = @(x) PA_SIRXi_w(x,beta,gamma,w,kap,delta);
[ts,xs] = ode45(@(t,x) ODE(x),[0 tmax],ini);
x.s = xs(:,1);% x(1) = rhos
x.i = xs(:,2);% x(2) = rhoi
x.r = xs(:,3);% x(3) = rhor
x.si = xs(:,4);% x(4) = rhosi
x.ss = xs(:,5);% x(5) = rhoss
    
figure;
pl(1) = plot(ts,x.s);
pl(1).Color = [100,0,232]/256;
pl(1).DisplayName = 'susceptible';
hold on
pl(2) = plot(ts,x.i);
pl(2).Color = [180,12,22]/256;
pl(2).DisplayName = 'infected';
pl(3) = plot(ts,x.r);
pl(3).Color = [10,112,122]/256;
pl(3).DisplayName = 'recovered';
pl(4) = plot(ts,x.si);
pl(4).Color = [10,212,122]/256;
pl(4).DisplayName = 'si-link-density';
pl(5) = plot(ts,x.ss);
pl(5).Color = [150,212,52]/256;
pl(5).DisplayName = 'ss-link-density';
ax = gca;
% ax.YLim = [0,10];
% ax.YScale='log';
lgd = legend(pl);
lgd.Location = 'northeast';
