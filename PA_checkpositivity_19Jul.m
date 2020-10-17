
tmax = 1000;
gamma = 1/40;
delta = 1/10;
mu = 15;
wrange = 0:0.002:0.04;
krange = 0:0.002:0.1;
brange = [0.0024,0.0025,0.0026];

rho0 = 0.001;
rhosi0 = rho0 * mu * (1-rho0);
rhoss0 = (mu/2)- rhosi0;
ini=[1-rho0,rho0,0,rhosi0,rhoss0];

betacrit = @(w,k,gamma,mu) (w+k+gamma)/(mu+1);
kcrits = @(beta,gamma,w,mu) (beta*(mu-1) - gamma - w);

k_positives = ones(length(wrange),length(krange),length(brange));
k_sum = ones(length(wrange),length(krange),length(brange));
k_criticals = ones(length(wrange),length(brange));
k_doof = ones(length(wrange),length(krange),length(brange));
k_deef = ones(length(wrange),length(krange),length(brange));
k_neg = ones(length(wrange),length(krange),length(brange));
for bi = 1:length(brange)
    beta = brange(bi);
    disp('beta')
    disp(beta)
    for wi = 1:length(wrange)
        w = wrange(wi);
        disp('rewiring is')
        disp(w)
%         k_positive = NaN;
        k_criticals_k = ones(1,length(krange));
        k_positives_k = ones(1,length(krange));
        k_sum_k = ones(1,length(krange));
        k_neg_k = ones(1,length(krange));
        k_deef_k = ones(1,length(krange));
        parfor ki = 1:length(krange)
            k = krange(ki);
            ODE = @(x) PA_SIRXi_w(x,beta,gamma,w,k,delta);
            [ts,xs] = ode45(@(t,x) ODE(x),[0 tmax],ini);
            
            fsi = -(beta+gamma+w+k)*xs(:,4) + beta*xs(:,4).*(2*xs(:,5)-xs(:,4))./xs(:,1);
            fss = -2*beta*xs(:,4).*xs(:,5)./xs(:,1) + w*(xs(:,1)./(xs(:,1)+xs(:,3))).*xs(:,4);
            diff = fsi - 2*fss;
            doof = sum(diff(diff<0));
            deef = sum(diff(diff>0));
            sdiff = sum(diff);
            if sdiff < 0
                valid_ks = - sum(diff(diff>0))/sum(abs(diff));
            else
                valid_ks = - sum(diff(diff<0))/sum(abs(diff));
            end
            k_positives_k(ki) = valid_ks;
            k_sum_k(ki) = sdiff;
            k_neg_k(ki) = doof;
            k_deef_k(ki) = deef;
%             if abs(valid_ks) > 0.00000001
%                 k_positive=k;
%                 break;
%             end
            
        end
        k_positives(wi,:,bi) = k_positives_k;
        k_criticals(wi,bi) = kcrits(beta,gamma,w,mu);
        k_sum(wi,:,bi) = k_sum_k;
        k_doof(wi,:,bi) = k_neg_k;
        k_deef(wi,:,bi) = k_deef_k;
        k_neg(wi,:,bi) = k_neg_k;
    end
end
% 

bi =1;

ql =1;
clear('ql');
minval=-0.03;
maxval=0.05;

figure; 
plot(wrange,zeros(1,length(wrange)));
hold on;
dw=wrange(2)-wrange(1);
dk=krange(2)-krange(1);
for wi=1:length(wrange)
    w = wrange(wi);
    for ki=1:length(krange)
        k = krange(ki);
        r = rectangle('Position',[w k dw dk]');
        
        val = k_sum(wi,ki,2);
        if isnan(val)
            r.FaceColor = 'white';
            r.EdgeColor = 'white';
        else
            scaled = (val-minval)/(maxval-minval);
            if scaled>=1
                scaled = 1;
            end
            if scaled<=0
                scaled = 0;
            end
            r.FaceColor = [1,0,0]*scaled +[0,0.5,1]*(1-scaled);
            r.EdgeColor = [1,0,0]*scaled +[0,0.5,1]*(1-scaled);
        end
        r.LineWidth = 0.0001;
    end
end
% ql=plot(wfine,kapfine);
% ql.LineWidth=3;
% ql.Color='black';
% ql.DisplayName = 'critical curve (MF)';
ax=gca;
ax.YLim=[0,0.1];
ax.XLim=[0,0.04];
ax.XLabel.String='w';
ax.YLabel.String='\kappa';
title('total contribution at \beta=0.0025 (Simulation)')
colorMap = [linspace(0,1,256)',linspace(0.5,0,256)',linspace(1,0,256)'];
colormap(colorMap);
cb=colorbar;
cb.LineWidth=0.5;
cb.Ticks = [0:0.25:1];
cb.TicksMode ='manual';
cb.TickLabels = minval + [0:0.25:1]*(maxval-minval);
% 
% lgd=legend(ql);
% lgd.Location='northeast';
% lgd.Box='off';
% lgd.BoxFace.ColorType='truecoloralpha';
% lgd.FontSize=LegendFontsizes;

resolution=300;
folder='figures/';
filename=strcat('Simulation_2d_delta2_atrate_date_18Jun');



ql =1;
clear('ql');
minval=-0.005;
maxval=0.005;

figure; 
plot(wrange,zeros(1,length(wrange)));
hold on;
dw=wrange(2)-wrange(1);
dk=krange(2)-krange(1);
for wi=1:length(wrange)
    w = wrange(wi);
    for ki=1:length(krange)
        k = krange(ki);
        r = rectangle('Position',[w k dw dk]');
        
        val = k_deef(wi,ki,2);
        if isnan(val)
            r.FaceColor = 'white';
            r.EdgeColor = 'white';
        else
            scaled = (val-minval)/(maxval-minval);
            if scaled>=1
                scaled = 1;
            end
            if scaled<=0
                scaled = 0;
            end
            r.FaceColor = [1,0,0]*scaled +[0,0.5,1]*(1-scaled);
            r.EdgeColor = [1,0,0]*scaled +[0,0.5,1]*(1-scaled);
        end
        r.LineWidth = 0.0001;
    end
end
% ql=plot(wfine,kapfine);
% ql.LineWidth=3;
% ql.Color='black';
% ql.DisplayName = 'critical curve (MF)';
ax=gca;
ax.YLim=[0,0.1];
ax.XLim=[0,0.04];
ax.XLabel.String='w';
ax.YLabel.String='\kappa';
title('positive contribution at \beta=0.0025 (Simulation)')
colorMap = [linspace(0,1,256)',linspace(0.5,0,256)',linspace(1,0,256)'];
colormap(colorMap);
cb=colorbar;
cb.LineWidth=0.5;
cb.Ticks = [0:0.25:1];
cb.TicksMode ='manual';
cb.TickLabels = minval + [0:0.25:1]*(maxval-minval);
% 
% lgd=legend(ql);
% lgd.Location='northeast';
% lgd.Box='off';
% lgd.BoxFace.ColorType='truecoloralpha';
% lgd.FontSize=LegendFontsizes;

resolution=300;
folder='figures/';
filename=strcat('Simulation_2d_delta2_atrate_date_18Jun');



ql =1;
clear('ql');
minval=-0.005;
maxval=0.005;

figure; 
plot(wrange,zeros(1,length(wrange)));
hold on;
dw=wrange(2)-wrange(1);
dk=krange(2)-krange(1);
for wi=1:length(wrange)
    w = wrange(wi);
    for ki=1:length(krange)
        k = krange(ki);
        r = rectangle('Position',[w k dw dk]');
        
        val = k_doof(wi,ki,2);
        if isnan(val)
            r.FaceColor = 'white';
            r.EdgeColor = 'white';
        else
            scaled = (val-minval)/(maxval-minval);
            if scaled>=1
                scaled = 1;
            end
            if scaled<=0
                scaled = 0;
            end
            r.FaceColor = [1,0,0]*scaled +[0,0.5,1]*(1-scaled);
            r.EdgeColor = [1,0,0]*scaled +[0,0.5,1]*(1-scaled);
        end
        r.LineWidth = 0.0001;
    end
end
% ql=plot(wfine,kapfine);
% ql.LineWidth=3;
% ql.Color='black';
% ql.DisplayName = 'critical curve (MF)';
ax=gca;
ax.YLim=[0,0.1];
ax.XLim=[0,0.04];
ax.XLabel.String='w';
ax.YLabel.String='\kappa';
title('negative contribution at \beta=0.0025 (Simulation)')
colorMap = [linspace(0,1,256)',linspace(0.5,0,256)',linspace(1,0,256)'];
colormap(colorMap);
cb=colorbar;
cb.LineWidth=0.5;
cb.Ticks = [0:0.25:1];
cb.TicksMode ='manual';
cb.TickLabels = minval + [0:0.25:1]*(maxval-minval);
% 
% lgd=legend(ql);
% lgd.Location='northeast';
% lgd.Box='off';
% lgd.BoxFace.ColorType='truecoloralpha';
% lgd.FontSize=LegendFontsizes;

resolution=300;
folder='figures/';
filename=strcat('Simulation_2d_delta2_atrate_date_18Jun');



% direction=strcat(folder,filename,'.png');
% saveas(gcf,direction)
% 
% figure; 
% plot(wrange,k_positives(:,2));
% hold on;
% plot(wrange,k_criticals(:,2));
% ax=gca;
% ax.XLim=[0,0.055];
% ax.YLim=[0,0.055];