tmax = 500;
gamma = 1/40;
delta = 1/10;
mu = 15;
w = 0.023;
k = 0.006;
bet = 0.0025;

rho0 = 0.001;
rhosi0 = rho0 * mu * (1-rho0);
rhoss0 = (mu/2)- rhosi0;
ini=[1-rho0,rho0,0,rhosi0,rhoss0];

betcrit = @(w,k,gamma,mu) (w+k+gamma)/(mu+1);
kcrits = @(bet,gamma,w,mu) (bet*(mu-1) - gamma - w);

ODE = @(x) PA_SIRXi_w(x,bet,gamma,w,k,delta);
[ts,xs] = ode45(@(t,x) ODE(x),[0 tmax],ini);

fsi = -(bet+gamma+w+k)*xs(:,4) + bet*xs(:,4).*(2*xs(:,5)-xs(:,4))./xs(:,1);
fss = -2*bet*xs(:,4).*xs(:,5)./xs(:,1) + w*(xs(:,1)./(xs(:,1)+xs(:,3))).*xs(:,4);
diff = fsi - 2*fss;

figure;
plot(ts,diff)