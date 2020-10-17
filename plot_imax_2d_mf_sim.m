% 
% N = 500;
% SA = 120;
% rob_k = 8;
% rob_perc = (rob_k*2 / SA ) * 100;
% mu = 15;
% tmax = 1000;
% tspan = 0:1:tmax;
% gamma = 1/40;
% delta = 1/10;
% beta = 0.0025; % up to 0.01;
% % beta = 0.015; % up to 0.2;
% wrange = 0: 0.01/50 : 0.01;
% krange = 0: 0.01/50 : 0.01;
% % wrange = 0: 0.2/50 : 0.2;
% % krange = 0: 0.2/50 : 0.2;
% 
% mf = zeros(length(wrange),length(krange));
% sm = zeros(length(wrange),length(krange));
% 
% kcrits = @(beta,gamma,w,mu) (beta*(mu-1) - gamma - w);
% kcrits(beta,gamma,0,mu)
% 
% tic
% for wi=1:length(wrange)
%     w=wrange(wi);
%     s = toc;
%     for ki=1:length(krange)
%         ss = toc;
%         k=krange(ki);
%         disp(strcat('beta=',num2str(beta),', w=',num2str(w),', k=',num2str(k)))
%         L=ceil(mu*N/2);
%         A=ERG(N,L);
%         rho0 = 0.005; % change back to 0.001
%         x=InitState(rho0,N);
%         si0=((x==0)*A)*(x==1)'/N;
%         i0 =sum(x)/N;
%         rhos=1-i0; rhor=0;
%         rhoss= (mu/2)-si0;
%         ini=[rhos,i0,rhor,si0,rhoss];
%     
%         ODE = @(x) PA_SIRXi_w(x,beta,gamma,w,k,delta);
%         [~,xs] = ode45(@(t,x) ODE(x),tspan,ini);
%         mf(wi,ki) = max(xs(:,2));
%         
%         maxes = zeros(1,SA);
%         parfor sa = 1:SA
%             result = SIRXi_w_fromAandx(A,x,beta,gamma,w,k,delta,tmax,tmax); 
%             rhoi = result.NIs/result.N;
%             marhoi = movmean(rhoi,65);
%             if max(marhoi) > (1+0.0)*marhoi(1)
%                 maxes(sa) = max(marhoi);
%             else
%                 maxes(sa) = NaN;
%             end
%         end
%         if ~all(isnan(maxes))
%             nonnans = ~isnan(maxes);
%             sm(wi,ki) = trimmean(maxes(nonnans),rob_perc);
%         else
%             sm(wi,ki) = 0;
%         end
%         tt = toc;
%     end
%    
% end

load('Imax_2D_Aug_withhigherres.mat')

% % load('MFvsSIMin2Datrate.mat') %,'atrate_mf','atrate_sm')
% % load('MFvsSIMin2DatrateForImax.mat') %,'atrate_mf','atrate_sm')
% atrate_mf(atrate_mf==0)=nan;
% atrate_sm(atrate_sm==0)=nan;

LegendFontsizes = 12;
% gamma=1/40;
% bc = 0.0025;
wfine=linspace(0.0,0.01,1000);
kfine=(mu-1)*beta - gamma - wfine;

minval = 0.000;
maxval = 0.06;

figure; 
plot(wrange,zeros(1,length(wrange)));
hold on;
dw=wrange(2)-wrange(1);
dk=krange(2)-krange(1);
for wi=1:length(wrange)
    w = wrange(wi);
    for ki=1:length(krange)
        val = mf(wi,ki);
        k = krange(ki);
        r = rectangle('Position',[w k dw dk]');
        if isnan(val)
            r.FaceColor = 'white';
            r.EdgeColor = 'white';
        else
            scaled = (val-minval)/(maxval-minval);
            r.FaceColor = [1,0,0]*scaled +[0,0.5,1]*(1-scaled);
            r.EdgeColor = [1,0,0]*scaled +[0,0.5,1]*(1-scaled);
        end
        r.LineWidth = 0.0001;
    end
end

ql=plot(wfine,kfine);
ql.LineWidth=3;
ql.Color='black';
ql.DisplayName = 'critical curve (MF)';
ax=gca;
ax.YLim=[0,0.01];
ax.XLim=[0,0.01];
ax.XLabel.String='w';
ax.YLabel.String='\kappa';
title('I_{max} at \beta=0.0025 (Mean Field)')

colorMap = [linspace(0,1,256)',linspace(0.5,0,256)',linspace(1,0,256)'];
colormap(colorMap);
cb=colorbar;
cb.LineWidth=0.5;
cb.Ticks = [0:0.25:1];
cb.TicksMode ='manual';
cb.TickLabels = minval + [0:0.25:1]*(maxval-minval);

lgd=legend(ql);
lgd.Location='northeast';
lgd.Box='off';
lgd.BoxFace.ColorType='truecoloralpha';
lgd.FontSize=LegendFontsizes;

resolution=300;
folder='figures/';
filename=strcat('Imax_MF_date_01Sep');
direction=strcat(folder,filename,'.png');
saveas(gcf,direction)

clear('ql');

figure; 
plot(wrange,zeros(1,length(wrange)));
hold on;
for wi=1:length(wrange)
    w = wrange(wi);
    for ki=1:length(krange)
        val = sm(wi,ki);
        k = krange(ki);
        r = rectangle('Position',[w k dw dk]');
        if isnan(val)
            r.FaceColor = 'white';
            r.EdgeColor = 'white';
        else
            scaled = (val-minval)/(maxval-minval);
            r.FaceColor = [1,0,0]*scaled +[0,0.5,1]*(1-scaled);
            r.EdgeColor = [1,0,0]*scaled +[0,0.5,1]*(1-scaled);
        end
        r.LineWidth = 0.0001;
    end
end
ql=plot(wfine,kfine);
ql.LineWidth=3;
ql.Color='black';
ql.DisplayName = 'critical curve (MF)';
ax=gca;
ax.YLim=[0,0.01];
ax.XLim=[0,0.01];
ax.XLabel.String='w';
ax.YLabel.String='\kappa';
title('I_{max} at \beta=0.0025 (Simulation)')
colorMap = [linspace(0,1,256)',linspace(0.5,0,256)',linspace(1,0,256)'];
colormap(colorMap);
cb=colorbar;
cb.LineWidth=0.5;
cb.Ticks = [0:0.25:1];
cb.TicksMode ='manual';
cb.TickLabels = minval + [0:0.25:1]*(maxval-minval);

lgd=legend(ql);
lgd.Location='northeast';
lgd.Box='off';
lgd.BoxFace.ColorType='truecoloralpha';
lgd.FontSize=LegendFontsizes;

resolution=300;
folder='figures/';
filename=strcat('Imax_SIM_date_01Sep');

direction=strcat(folder,filename,'.png');
saveas(gcf,direction)