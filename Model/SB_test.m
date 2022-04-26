%%% test script 8. Jul
tmax = 900;
beta=0.004; gamma=1/40; w=1/100; 
kap=12/10000; delta=2/100;
betarange = [0.002:0.0001:0.006];
% rho0=1/100;
SA = 100; 
Nt = 1000;
N = 1000;
mu = 15;
rho0 = 0.01;

for sa = 1:6
    L=ceil(mu*N/2);
   
    A=ERG(N,L);
    x=InitState(rho0,N);
    si0=((x==0)*A)*(x==1)'/N;
    
    i0 =sum(x)/N;
    rhos=1-i0; rhor=0;
    rhoss= (mu/2)-si0;
    ini=[rhos,i0,rhor,si0,rhoss];
    
    result1 = SIRXi_w_fromAandx(A,x,beta,gamma,w,kap,delta,tmax,Nt); 
    ODE = @(x) PA_SIRXi_w(x,beta,gamma,w,kap,delta);
    [ts,xs] = ode45(@(t,x) ODE(x),[0 tmax],ini);
    Tt= min(Nt,length(ts));
    
    rhoss = result.NSSs/result.N;
    rhos_times_rhos = result.NSs .* result.NSs * mu/(2*result.N*result.N);
    rhor_2 = ((result.N - result.NRs).^2) * mu/(2*result.N*result.N);
    figure;
    pl(1) = plot(result.times,rhoss); 
    pl(1).Color = 'black';
    pl(1).LineStyle = '-';
    pl(1).LineWidth = 1;
    pl(1).DisplayName = '\rho_{SS} (simulation)';
    hold on; 
    pl(2) = plot(ts,xs(:,5));
    pl(2).Color = 'black';
    pl(2).LineStyle = '--';    
    pl(2).LineWidth = 1;
    pl(2).DisplayName = '\rho_{SS} (mean field)';
    
    pl(3) = plot(result.times,rhor_2);
    pl(3).Color = 'red';
    pl(3).LineStyle = '-';
    pl(3).LineWidth = 1;
    pl(3).DisplayName = '(\mu / 2 )(1-\rho_{R})^2 (simulation)';
    
    pl(4) = plot(ts,(mu/2)*(1-xs(:,3)).^2);
    pl(4).Color = 'red';
    pl(4).LineStyle = '--';  
    pl(4).LineWidth = 1;
    pl(4).DisplayName = '(\mu / 2 )(1-\rho_{R})^2 (mean field)';
    
    pl(5) = plot(result.times,rhos_times_rhos);
    pl(5).Color = 'blue';
    pl(5).LineStyle = '-';
    pl(5).LineWidth = 1;
    pl(5).DisplayName = '(\mu / 2 )\rho_{S}^2 (simulation)';
    
    pl(6) = plot(ts,(mu/2)*(xs(:,1).^2));
    pl(6).Color = 'blue';
    pl(6).LineStyle = '--';
    pl(6).LineWidth = 1;
    pl(6).DisplayName = '(\mu / 2 )\rho_{S}^2 (mean field)';
    
    ax= gca;
    ax.XLim=[0,tmax];
    
    lgd = legend(pl);
    lgd.Location = 'northeast';
    
    resolution=600;
    folder='figures/';
    filename=strcat('Jul09_MF_vs_simulation_approx',num2str(sa));
    direction=strcat(folder,filename,'.png');
    saveas(gcf,direction)
    
    figure; 
    ql(1) = plot(ts,(mu/2)*(1-xs(:,3)).^2);
    ql(1).Color = 'red';
    ql(1).LineStyle = '-';  
    ql(1).LineWidth = 1;
    ql(1).DisplayName = '(\mu / 2 )(1-\rho_{R})^2 (mean field)';
    hold on;
    ql(2) = plot(ts,xs(:,5));
    ql(2).Color = 'black';
    ql(2).LineStyle = '-';    
    ql(2).LineWidth = 1;
    ql(2).DisplayName = '\rho_{SS} (mean field)';
    ax= gca;
    ax.XLim=[0,tmax];
    
    lgd = legend(ql);
    lgd.Location = 'northeast';
    
    resolution=600;
    folder='figures/';
    filename=strcat('Jul09_MF_K',num2str(sa));
    direction=strcat(folder,filename,'.png');
    saveas(gcf,direction)
end