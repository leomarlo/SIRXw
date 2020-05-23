tmax=500;
Nt=1000;

N=1000;
mu=20; 
beta=1/150;
gamma=1/40;
w=1/100;
rho0=1/100;
kap0=3/10000;
kap=12/10000;
r=1; % one standard deviation

result=SIRXw2(N,mu,beta,gamma,w,kap0,kap,r,rho0,tmax,Nt);
figure;
plot(result.cases)
figure;
plot(result.times,result.NIs/result.N)

clear('pl')
figure;
pl(1)=plot(result.times, result.NIs/result.N);pl(1).Color='red';pl(1).DisplayName='infected';
hold on;
pl(2)=plot(result.times, result.NRs/result.N);pl(2).Color='blue'; pl(2).DisplayName='recovered';
pl(3)=plot(result.times, result.NSs/result.N);pl(3).Color='green';pl(3).DisplayName='susceptible';
pl(4)=plot(result.times, result.NXss/result.N);pl(4).Color='cyan';pl(4).DisplayName='quarantined susceptibles';
pl(5)=plot(result.times, result.NXis/result.N);pl(5).Color='black';pl(5).DisplayName='quarantined infected';
% pl(6)=plot(result.times, (result.NSs+result.NIs+result.NRs+result.NXis+result.NXss)/result.N);
ax=gca;
ax.YLim=[0,1];
legend(pl);
% lgd.Location='northeast';