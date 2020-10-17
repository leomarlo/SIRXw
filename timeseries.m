
beta=0.05; gamma=1/40; w=1/100; 
kap=12/10000; delta=2/100;
mu = 15;
rho0 = 0.001;
rhosi0 = rho0 * mu * (1-rho0);
rhoss0 = (mu/2)- rhosi0;
ini=[1-rho0,rho0,0,rhosi0,rhoss0];


wrange= [0/10: 1/10 :5/10];
kapparange = 5/10 - wrange;

figure;
for i=1:length(wrange)
    w = wrange(i);
    kappa = kapparange(i);
    ODE = @(x) PA_SIRXi_w(x,beta,gamma,w,kap,delta);
    [ts,xs] = ode45(@(t,x) ODE(x),[0 tmax],ini);
    
    plot(ts,xs)
    hold on;
    
end

pl = plot(ts,xs(:,2));
pl.Color = 'black';
pl.LineWidth = 3;