%%%%%% DIfferent Plots along constant w+kappa
% 
% 
tmax = 2000;
tspan = [0:1:tmax];
gamma = 1/40;
delta = 1/10;
mu = 15;
N = 1000;
wrange = 0:0.002:0.04;

rho0 = 0.001;
rhosi0 = rho0 * mu * (1-rho0);
rhoss0 = (mu/2)- rhosi0;
ini=[1-rho0,rho0,0,rhosi0,rhoss0];

brange = [0.001:0.0005:0.015];
wpluskrange = [5/10000:5/10000:2/10];
kfracrange = [0.25:0.25:0.75];
SA = 5;
% 
% % storesm = zeros(tmax+1,3,length(brange),length(wpluskrange),length(kfracrange), SA);
% storemf = zeros(tmax+1,3,length(brange),length(wpluskrange),length(kfracrange)); 
% maxmf = zeros(3,length(brange),length(wpluskrange),length(kfracrange)); 
% storecr = zeros(length(brange),length(wpluskrange));
% for bi = 1:length(brange)
%     beta = brange(bi);
%     disp(beta);
%     for wki = 1:length(wpluskrange)
%         wplusk = wpluskrange(wki);
%         critbeta = (wplusk + gamma )/(mu + 1);
%         storecr(bi,wki) = (critbeta < beta) * 1;
%         for kfi = 1:length(kfracrange)
%             kfrac = kfracrange(kfi);
%             k = wplusk * kfrac;
%             w = wplusk - k;
% %             for sa = 1:SA
% %             end
%             ODE = @(x) PA_SIRXi_w(x,beta,gamma,w,k,delta);
%             [ts,xs] = ode45(@(t,x) ODE(x),tspan,ini);
%             storemf(:,1,bi,wki,kfi) = xs(:,2);
%             storemf(:,2,bi,wki,kfi) = xs(:,4);
%             storemf(:,3,bi,wki,kfi) = xs(:,3);
%             maxmf(1,bi,wki,kfi) = ts(find(xs(:,2)==max(xs(:,2)),1,'first'));
%             maxmf(2,bi,wki,kfi) = ts(find(xs(:,4)==max(xs(:,4)),1,'first'));
%             maxmf(3,bi,wki,kfi) = ts(find(xs(:,3)==max(xs(:,3)),1,'first'));
%             fff = 1;
%         end
%     end
% end
% 
% %%
% % beta = 0.0025
% % wplusk = 0.004
% 
% % beta = 0.005
% % wplusk = 0.04
% 
% % beta = 0.013
% % wplusk = 0.03
% 
% % bi = 7
% % wki = 1
% 
% % bi 15
% % wki 150
% %%
% 
% 
% % bi = 29;
% % b = brange(bi);
% % wki = 300;
% % wk = wpluskrange(wki);
% % % kfi = 2;
% % disp(storecr(bi,wki))
% % disp(b)
% % disp(wk)



load('examplepaths24Aug.mat')
% 
clear('pls');
twosets = [[7,1];[7,50];[29,300]];
cols = {'green','blue','red','cyan'};
linestyles = {'-','--',':'};
I_SI_R = {'I','SI','R'};
nr = length(wpluskrange);
TN = 800;
ts = tspan(1:TN);
% cols = [[0:1:nr]/nr ; [0:1:nr]/nr; [nr:-1:0]/nr]';
figure; 
for ti = 1:length(twosets)
    params = twosets(ti,:);
    bi = params(1);
    beta = brange(bi);
    wki = params(2);
    wplusk = wpluskrange(wki);
    
    for kfi=length(kfracrange):-1:1
        pls(ti) = plot(ts,storemf(1:TN,1,bi,wki,kfi));
        pls(ti).Color = cols{ti};
        pls(ti).LineStyle = linestyles{kfi};
        pls(ti).LineWidth = 3.5;
        pls(ti).DisplayName = strcat('\beta=',num2str(beta),' \omega+\kappa=',num2str(wplusk));
        hold on;
        pr = plot([maxmf(1,bi,wki,kfi),maxmf(1,bi,wki,kfi)],[0,1]);
        pr.Color = cols{ti};
        pr.LineStyle = linestyles{kfi};
        pr.LineWidth = 1.5;

    end

end

ax = gca;
ax.YLim = [0,0.2];
ax.XLabel.String = 'Time';
ax.YLabel.String = 'Prevalence (\rho_{I})';

lgd = legend(pls);
lgd.Location = 'northeast';
resolution=300;

folder='figures/';
filename=strcat('examples_10Sept');
direction=strcat(folder,filename,'.png');
saveas(gcf,direction)


% some more testing

bi = 19;
TN = 1000;
ts = tspan(1:TN);
for wki=1:1
%     if mod(wki,6)==1
        figure;
        for kfi=length(kfracrange):-1:1
            plot(ts,storemf(1:TN,1,bi,wki,kfi));
            hold on;
        end
        ax = gca;
        ax.YLim = [0,0.07];
        ax.XLim = [0,max(ts)];
%     end
%     if wki>100
%         break
%     end
end

critwplk = @(beta,gamma,mu) beta * (mu-1) - gamma;
beta = 0.015; wplusk = 0.005; % high infection rate, away from transition. at 0.6 max
beta = 0.015; wplusk = 0.1;

% beta = 0.01; wplusk = 0.0025;
% beta = 0.01; wplusk = 0.075;
% beta = 0.0025; wplusk = 0.00025;
% beta = 0.0025; wplusk = 0.0075;


% 

params = [[0.01,0.0025]; [0.01, 0.075]; [0.0025, 0.00025]; [0.0025, 0.0075]];
    
Lp = length(params);

colors = {'black','red','blue'};
styles = {'-','--',':'};
lwidths = [3,1,2];
cols = zeros(2,2,3);
lwds = zeros(2,2,3);
stys = zeros(2,2,3);
TN = 1500;
ts = tspan(1:TN);
for bi = 1:2 
    for wki = 1:2
        for kfr = 1:3
            stys(bi,wki,kfr) = kfr;  
            cols(bi,wki,kfr) = wki;
            lwds(bi,wki,kfr) = bi;
        end
    end
end

kfrange = [0.25,0.5,0.75];
figure;
for j = 1:Lp
    beta = params(j,1);
    
    wplk = params(j,2);
    wki = mod(j,2)+1;
    bi = ceil(j/2);
    for kfi = 1:length(kfrange)
        kfr = kfrange(kfi);
        k = wplk * kfr;
        w = wplk - k;
        ODE = @(x) PA_SIRXi_w(x,beta,gamma,w,k,delta);
        [ts,xs] = ode45(@(t,x) ODE(x),tspan,ini);
        pl(j) = plot(ts,xs(:,2));
        pl(j).LineWidth = lwidths(lwds(bi,wki,kfi));
        pl(j).LineStyle = styles{stys(bi,wki,kfi)};
        pl(j).Color = colors{cols(bi,wki,kfi)};
        pl(j).DisplayName = strcat('\beat=',num2str(beta),', w+\kappa=',num2str(wplk));
        hold on;
    end
end
ax = gca; 
ax.XLim = [0, TN];
ax.YLim = [0, 0.1];

lgd = legend(pl);
lgd.Location = 'northeast';




params1 = [0.0025, 0.075]; 
params2 = [0.00025, 0.0075];
params = [params2;params1];
betas =  [0.0025, 0.01];
kfrange = [0.25,0.5,0.75];

tspan = 1:1:3000;
limis =[[0.05,0.5];[0.5,1]];
xlims = [2900,650];

tmax = 2000;
tspan = [0:1:tmax];
gamma = 1/40;
delta = 1/10;
mu = 15;
% N = 1000;
wrange = 0:0.002:0.04;

rho0 = 0.001;
rhosi0 = rho0 * mu * (1-rho0);
rhoss0 = (mu/2)- rhosi0;
ini=[1-rho0,rho0,0,rhosi0,rhoss0];


for bi = 1:length(betas)
    beta = betas(bi);
    figure;
    TN = xlims(bi);
    ts = tspan(1:TN);
    for wki = 1:2
        wplk = params(bi,wki);
        for kfi = 1:length(kfrange)
            kfr = kfrange(kfi);
            k = wplk * kfr;
            w = wplk - k;
            ODE = @(x) PA_SIRXi_w(x,beta,gamma,w,k,delta);
            [ts,xs] = ode45(@(t,x) ODE(x),tspan,ini);
            pl(wki) = plot(ts,xs(:,2));
            pl(wki).LineWidth = 2;
            pl(wki).LineStyle = styles{kfi};
            pl(wki).Color = colors{wki};
            pl(wki).DisplayName = strcat('\beat=',num2str(beta),', w+\kappa=',num2str(wplk));
            hold on;
        end
    end
    ax1 = gca;
    ax1.XLim = [0,TN];
    ax1.YLim = [0,limis(bi,1)];
    ax1.XLabel.String = 'time';
    ax1.YLabel.String = '\rho_{I}';
    ax1.YTick=[0:limis(bi,1)/5:limis(bi,1)];
    ax1.FontSize=16;
    ax1.FontName='CMU Serif';
    
    axes('Position',[.45 .4 .45 .48]);
    box on;
    for wki = 1:2
        wplk = params(bi,wki);
        for kfi = 1:length(kfrange)
            kfr = kfrange(kfi);
            k = wplk * kfr;
            w = wplk - k;
            ODE = @(x) PA_SIRXi_w(x,beta,gamma,w,k,delta);
            [ts,xs] = ode45(@(t,x) ODE(x),tspan,ini);
            ql(wki) = plot(ts,xs(:,3));
            ql(wki).LineWidth = 2;
            ql(wki).LineStyle = styles{kfi};
            ql(wki).Color = colors{wki};
            ql(wki).DisplayName = strcat('\beat=',num2str(beta),', w+\kappa=',num2str(wplk));
            hold on;
        end
    end
    ax2 = gca;
    ax2.XLim = [0,TN];
    ax2.YLim = [0,limis(bi,2)];
    ax2.XLabel.String = 'time';
    ax2.YLabel.String = '\rho_{R}';
    ax2.YTick=[0:limis(bi,2)/5:limis(bi,2)];
    ax2.FontSize=12;
    ax2.FontName='CMU Serif';
    
    folder='figures/';
    filename=strcat('examples_',num2str(bi),date);
    direction=strcat(folder,filename,'.png');
%     saveas(gcf,direction)

end

    

