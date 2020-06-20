% load('sweepMay24N500.mat')

load('sweepJun06N500delta002.mat')

trans_thr = 0.05;

crit_mf = zeros(length(wrange),length(kaprange));
crit_sm = zeros(length(wrange),length(kaprange));
betasteps = 1000;
for wi=1:length(wrange)
    w=wrange(wi)
    for ki=1:length(kaprange)
        kap=kaprange(ki);
        rho_est  = ones(1,length(brange));
        weights = ones(1,length(brange));
        rinf_mf = zeros(1,length(brange));
        for bi = 1:length(brange)
            X = simRinf{wi,ki,bi};
            [xmean,xsdev,weight] = corrected(X,0.05,0.05);
            weights(bi) = weight/(xsdev+0.1);
            rho_est(bi) = xmean;
            
            rinf_mf(bi)=mfaRinf{wi,ki,bi};
        end
        
        
        % SIMULATION 
        sigmas = 2; % spikeyness of gaussian convolution 
        % (relative importance of the closeby data points, closeness to linear
        % interpolation)
        [ma,bs] = movingaverage(rho_est,weights,brange,betasteps,sigmas);
        crit_sm(wi,ki) = bs(find(ma>trans_thr,1,'first'));
        % MF threshold interpolate
        mf_interp = interp1q(brange',rinf_mf',bs'); 
        crit_mf(wi,ki) = bs(find(mf_interp>trans_thr,1,'first'));
        
    end
end
% 
% load('MFvsSIMin2D.mat') %'crit_mf','crit_sm')



minval = 0.0008;
maxval = 0.0033;
figure; 
plot(wrange,zeros(1,length(wrange)));
hold on;
dw=wrange(2)-wrange(1);
dk=kaprange(2)-kaprange(1);
for wi=1:length(wrange)
    w = wrange(wi);
    for ki=1:length(kaprange)
        val = crit_mf(wi,ki);
        scaled = (val-minval)/(maxval-minval);
        kap = kaprange(ki);
        r = rectangle('Position',[w kap dw dk]');
        r.FaceColor = [1,0,0]*scaled +[0,0.5,1]*(1-scaled);
        r.EdgeColor = [1,0,0]*scaled +[0,0.5,1]*(1-scaled);
        r.LineWidth = 0.0001;
    end
end
ax=gca;
ax.YLim=[0,0.01];
ax.XLim=[0,0.01];
ax.XLabel.String='w';
ax.YLabel.String='\kappa';
title('critical infection rate \beta (Mean Field)')

colorMap = [linspace(0,1,256)',linspace(0.5,0,256)',linspace(1,0,256)'];
colormap(colorMap);
cb=colorbar;
cb.LineWidth=0.5;
cb.Ticks = [0:0.2:1];
cb.TicksMode ='manual';
cb.TickLabels = minval + [0:0.2:1]*(maxval-minval);

resolution=300;
folder='figures/';
filename=strcat('MeanField_2d_002delta_date_18Jun');

direction=strcat(folder,filename,'.png');
saveas(gcf,direction)

figure; 
plot(wrange,zeros(1,length(wrange)));
hold on;
dw=wrange(2)-wrange(1);
dk=kaprange(2)-kaprange(1);
for wi=1:length(wrange)
    w = wrange(wi);
    for ki=1:length(kaprange)
        val = crit_sm(wi,ki);
        scaled = (val-minval)/(maxval-minval);
        kap = kaprange(ki);
        r = rectangle('Position',[w kap dw dk]');
        r.FaceColor = [1,0,0]*scaled +[0,0.5,1]*(1-scaled);
        r.EdgeColor = [1,0,0]*scaled +[0,0.5,1]*(1-scaled);
        r.LineWidth = 0.0001;
    end
end
ax=gca;
ax.YLim=[0,0.01];
ax.XLim=[0,0.01];
ax.XLabel.String='w';
ax.YLabel.String='\kappa';
title('critical infection rate \beta (Simulation)')

colorMap = [linspace(0,1,256)',linspace(0.5,0,256)',linspace(1,0,256)'];
colormap(colorMap);
cb=colorbar;
cb.LineWidth=0.5;
cb.Ticks = [0:0.2:1];
cb.TicksMode ='manual';
cb.TickLabels = minval + [0:0.2:1]*(maxval-minval);

resolution=300;
folder='figures/';
filename=strcat('Simulation_2d_002delta_date_18Jun');

direction=strcat(folder,filename,'.png');
saveas(gcf,direction)

