%% plot imax sim vs num

%% consider three parameter regimes
N = 500;
SA = 24;
rob_k = 2;
rob_perc = (rob_k*2 / SA ) * 100;
mu = 15;
tmax = 1000;
tspan = 0:1:tmax;
gamma = 1/40;
delta = 1/10;
params = [[0.004, 0.0005, 0.25];
          [0.004, 0.0005, 0.75];
          [0.004, 0.025, 0.25];
          [0.004, 0.025, 0.75];
          [0.015, 0.15, 0.25];
          [0.015, 0.15, 0.75]];

params = [[0.004, 0.0005, 0.25];
          [0.004, 0.0005, 0.75];
          [0.0035, 0.025, 0.25];
          [0.0035, 0.025, 0.75];
          [0.012, 0.15, 0.25];
          [0.012, 0.15, 0.75]];

betacrit = @(w,k,gamma,mu) (w+k+gamma)/(mu+1);
kcrits = @(beta,gamma,w,mu) (beta*(mu-1) - gamma - w);

sp = size(params);
storerhos_for_all = zeros(length(tspan),sp(1),SA);
marhos_for_all = zeros(length(tspan),sp(1),SA);
maxes_for_all = zeros(sp(1),SA);
rinfs_for_all = zeros(sp(1),SA);
% for i = 1:length(params)
tic;
time = toc;
for i = 3:6
    beta = params(i,1);
    wplk = params(i,2);
    k = wplk * params(i,3);
    w = wplk - k;
    
    disp(strcat('beta=',num2str(beta),', w+k=',num2str(wplk),', w=',num2str(w)))
    
    L=ceil(mu*N/2);

    A=ERG(N,L);
    rho0 = 0.005; % change back to 0.001
    x=InitState(rho0,N);
    si0=((x==0)*A)*(x==1)'/N;

    i0 =sum(x)/N;
    rhos=1-i0; rhor=0;
    rhoss= (mu/2)-si0;
    ini=[rhos,i0,rhor,si0,rhoss];
    storerhos = zeros(length(tspan),SA);
    marhos = zeros(length(tspan),SA);
    maxes = zeros(1,SA);
    rinfs = zeros(1,SA);
    
    ODE = @(x) PA_SIRXi_w(x,beta,gamma,w,k,delta);
    [~,xs] = ode45(@(t,x) ODE(x),tspan,ini);
    mfmin= xs(1,2);
    mfmax= max(xs(:,2));
    rinf_mf = xs(end,3);
    parfor sa = 1:SA
        disp(sa);
        result = SIRXi_w_fromAandx(A,x,beta,gamma,w,k,delta,tmax,tmax); 
        rhoi = result.NIs/result.N;
        marhoi = movmean(rhoi,65);
        storerhos(:,sa) = rhoi;
        marhos(:,sa) = marhoi
        mxrho = max(marhoi);
        rinfs(sa) = result.NRs(end)/result.N;
%         if (mxrho-mfmin)/(mfmax-mfmin) > 0.2
        disp(strcat('maximal rhoi=',num2str(max(marhoi)),', 1.2 times initial rhoi=', num2str((1+0.2)*marhoi(1)) ))
        disp(strcat('rinf_mf=',num2str(rinf_mf),', rinf_sim=',num2str(rinfs(sa))));
        if max(marhoi) > (1+0.0)*marhoi(1)
            maxes(sa) = mxrho;
        else
            maxes(sa) = NaN;
        end
    end
    maxes_for_all(i,:) = maxes;
    rinfs_for_all(i,:) = rinfs;
    for sa=1:SA
        marhos_for_all(:,i,sa) = marhos(:,sa);
        storerhos_for_all(:,i,sa) = storerhos(:,sa);
    end 
    legendplots = [];
    figure; 
    pl = plot(tspan,xs(:,2));
    pl.Color = 'black';
    pl.LineWidth = 1;
    pl.LineStyle = '-';
    hold on;
    for sa = 1:SA
        pls2 = plot(tspan,storerhos(:,sa));
        pls2.Color = [0.7,0.7,0.7];
        pls2.LineWidth = 0.5;
    end
    for sa = 1:SA
        pls = plot(tspan,marhos(:,sa));
        pls.Color = [0.5,0.2,0.2];
        pls.LineWidth = 1;
        qls = plot([0,tmax],[max(marhos(:,sa)),max(marhos(:,sa))]);
        if isnan(maxes(sa))
            qls.LineStyle = ':';
        else
            qls.LineStyle = '--';
        end
        qls.Color = [0.5,0.2,0.2];
    end
    pl = plot(tspan,xs(:,2));
    pl.Color = 'black';
    pl.LineWidth = 3;
    pl.LineStyle = '-';
    pl.DisplayName = 'PA';
    legendplots = [legendplots, pl];
    mxs = plot([0,tmax],[nanmean(maxes),nanmean(maxes)]);
    mxs.LineWidth = 2;
    mxs.LineStyle = '-';
    mxs.Color = [0.1,0.2,0.7];
    mxs.DisplayName = '\langle I_{max}\rangle';
    legendplots = [legendplots, mxs];
    if ~all(isnan(maxes))
        
        mxs_rob = plot([0,tmax],[trimmean(maxes(~isnan(maxes)),rob_perc),trimmean(maxes(~isnan(maxes)),rob_perc)]);
        mxs_rob.LineWidth = 2;
        mxs_rob.LineStyle = '-';
        mxs_rob.Color = [0.7,0.2,0.1];
        mxs_rob.DisplayName = '\langle I_{max}\rangle (robust)';
        
        legendplots = [legendplots, mxs_rob];
    end
    title(strcat('beta=',num2str(beta),', w+k=',num2str(wplk),', w=',num2str(w)))
    
    lgd = legend(legendplots);
    lgd.Location = 'northeast';
    resolution=300;

    folder='figures/';
    filename=strcat('examples_lowerbeta_robust_',num2str(i),'_29Aug');
    direction=strcat(folder,filename,'.png');
    saveas(gcf,direction)
   
end

% disp(maxes_for_all)
% disp(rinfs_for_all)

sime = toc;
disp(sime-time)

% 
% for i = 1:6
%     beta = params(i,1);
%     wplk = params(i,2);
%     k = wplk * params(i,3);
%     w = wplk - k;
%     betac = betacrit(w,k,gamma,mu);
%    
%     strcat('beta=',num2str(beta),', w+k=',num2str(wplk),', w=',num2str(w),', beta_c=',num2str(betac))
%     
% end